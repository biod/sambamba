/*
    This file is part of Sambamba.
    Copyright (C) 2013    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
module sambamba.utils.common.filtering;

import sambamba.utils.common.queryparser;
import std.algorithm;
import std.stdio;

mixin template GlobalFilter() {
    static __gshared Filter read_filter; 

    bool passing(BamRead read) {
        return read_filter.accepts(read);
    }

    auto filtered(R)(R reads) {
        return std.algorithm.filter!passing(reads);
    }
}

Filter createFilterFromQuery(string query) {
    if (query is null)
        return new NullFilter();
    auto query_grammar = new QueryGrammar();
    auto node = query_grammar.parse(query);
    auto condition_node = cast(ConditionNode) node;
    if (condition_node is null) {
        stderr.writeln("filter string must represent a condition");
        return null;
    }
    return condition_node.condition;
}

/**
  Set of filters for alignments. 
  All share a common interface and can be easily combined.
*/

import std.regex;
import std.algorithm;
import std.conv;
import bio.bam.read;
import bio.bam.tagvalue;
import bio.bam.validation.alignment;

/// Common interface for all filters
interface Filter {
    bool accepts(ref BamRead a) const;
}

/// Filter which accepts all alignments
final class NullFilter : Filter {
    bool accepts(ref BamRead a) const {
        return true;
    }
}

/// Validating filter
final class ValidAlignmentFilter : Filter {
    
    bool accepts(ref BamRead a) const {
        return isValid(a);
    }
}

/// Intersection of two filters
final class AndFilter : Filter {
    private Filter _a, _b;

    this(Filter a, Filter b) { _a = a; _b = b; }

    bool accepts(ref BamRead a) const {
        return _a.accepts(a) && _b.accepts(a);
    }
}

/// Union of two filters
final class OrFilter : Filter {
    private Filter _a, _b;

    this(Filter a, Filter b) { _a = a, _b = b; }

    bool accepts(ref BamRead a) const {
        return _a.accepts(a) || _b.accepts(a);
    }
}

/// Negation of a filter
final class NotFilter : Filter {
    private Filter _a;

    this(Filter a) { _a = a; }
    bool accepts(ref BamRead a) const {
        return !_a.accepts(a);
    }
}

/// Filter alignments which has $(D flagname) flag set
final class FlagFilter(string flagname) : Filter {
    bool accepts(ref BamRead a) const {
        mixin("return a." ~ flagname ~ ";");
    }
}

/// Filtering integer fields
final class IntegerFieldFilter(string op) : Filter {
    private long _value;
    private string _fieldname;
    this(string fieldname, long value) {
        _fieldname = fieldname;
        _value = value;
    }
    bool accepts(ref BamRead a) const {
        switch(_fieldname) {
            case "ref_id": mixin("return a.ref_id " ~ op ~ "_value;");
            case "position": mixin("return a.position " ~ op ~ "_value;");
            case "mapping_quality": mixin("return a.mapping_quality " ~ op ~ "_value;");
            case "sequence_length": mixin("return a.sequence_length " ~ op ~ "_value;");
            case "mate_ref_id": mixin("return a.mate_ref_id " ~ op ~ "_value;");
            case "mate_position": mixin("return a.mate_position " ~ op ~ "_value;");
            case "template_length": mixin("return a.template_length " ~ op ~ "_value;");
            default: throw new Exception("unknown integer field '" ~ _fieldname ~ "'");
        }
    }
}

final class TagExistenceFilter(string op) : Filter {
    static assert(op == "==" || op == "!=");
    private string _tagname;
    private static bool _should_exist = op == "!=";
    this(string tagname, typeof(null) dummy) {
        _tagname = tagname;
    }
    bool accepts(ref BamRead a) const {
        auto v = a[_tagname];
        if (_should_exist) 
            return !v.is_nothing;
        else 
            return v.is_nothing;
    }
}

/// Filtering integer tags
final class IntegerTagFilter(string op) : Filter {
    private long _value;
    private string _tagname;

    this(string tagname, long value) {
        _tagname = tagname;
        _value = value;
    }

    bool accepts(ref BamRead a) const {
        auto v = a[_tagname];
        if (!v.is_integer && !v.is_float) 
            return false;
        if (v.is_float) {
            mixin(`return cast(float)v` ~ op ~ `_value;`);
        } else {
            mixin(`return cast(long)v` ~ op ~ `_value;`);
        }
    }
}

/// Filtering string fields
final class StringFieldFilter(string op) : Filter {
    private string _value;
    private string _fieldname;
    this(string fieldname, string value) {
        _fieldname = fieldname;
        _value = value;
    }
    bool accepts(ref BamRead a) const {
        switch(_fieldname) {
            case "read_name": mixin("return a.name " ~ op ~ " _value;");
            case "sequence": mixin("return cmp(a.sequence, _value) " ~ op ~ " 0;");
            case "cigar": mixin("return a.cigarString() " ~ op ~ " _value;");
            default: throw new Exception("unknown string field '" ~ _fieldname ~ "'");
        }
    }
}

/// Filtering string and character tags
final class StringTagFilter(string op) : Filter {
    private string _value;
    private string _tagname;

    this(string tagname, string value) {
        _tagname = tagname;
        _value = value;
    }

    bool accepts(ref BamRead a) const {
        auto v = a[_tagname];
        if (v.is_string) {
            mixin(`return cast(string)v` ~ op ~ `_value;`);
        } else if (v.is_character) {
            if (_value.length != 1)
                return false; // doesn't make sense to compare char with string
            mixin(`return cast(char)v` ~ op ~ `_value[0];`);
        } else {
            return false;
        }
    }
}

/// Filtering string fields with a regular expression
final class RegexpFieldFilter : Filter {
    private string _fieldname;
    private Regex!char _pattern;
    
    this(string fieldname, Regex!char pattern) {
        _fieldname = fieldname; 
        _pattern = pattern;
    }

    bool accepts(ref BamRead a) const {
        switch(_fieldname) {
            case "read_name": return !match(a.name, cast()_pattern).empty;
            case "sequence": return !match(to!string(a.sequence), cast()_pattern).empty;
            case "cigar": return !match(a.cigarString(), cast()_pattern).empty;
            default: throw new Exception("unknown string field '" ~ _fieldname ~ "'");
        }
    }
}

/// Filtering string tags with a regular expression
final class RegexpTagFilter : Filter { 
    private string _tagname;
    private Regex!char _pattern;
    
    this(string tagname, Regex!char pattern) {
        _tagname = tagname;
        _pattern = pattern;
    }

    bool accepts(ref BamRead a) const {
        auto v = a[_tagname];
        if (!v.is_string) {
            return false;
        }
        return !match(cast(string)v, cast()_pattern).empty;
    }
}
