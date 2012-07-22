/**
  Set of filters for alignments. 
  All share a common interface and can be easily combined.
*/
module filter;

import std.regex;
import std.algorithm;
import std.conv;
import alignment;
import tagvalue;
import validation.alignment;

/// Common interface for all filters
interface Filter {
    bool accepts(ref Alignment a) const;
}

/// Filter which accepts all alignments
final class NullFilter : Filter {
    bool accepts(ref Alignment a) const {
        return true;
    }
}

/// Validating filter
final class ValidAlignmentFilter : Filter {
    
    bool accepts(ref Alignment a) const {
        return isValid(a);
    }
}

/// Intersection of two filters
final class AndFilter : Filter {
    private Filter _a, _b;

    this(Filter a, Filter b) { _a = a; _b = b; }

    bool accepts(ref Alignment a) const {
        return _a.accepts(a) && _b.accepts(a);
    }
}

/// Union of two filters
final class OrFilter : Filter {
    private Filter _a, _b;

    this(Filter a, Filter b) { _a = a, _b = b; }

    bool accepts(ref Alignment a) const {
        return _a.accepts(a) || _b.accepts(a);
    }
}

/// Negation of a filter
final class NotFilter : Filter {
    private Filter _a;

    this(Filter a) { _a = a; }
    bool accepts(ref Alignment a) const {
        return !_a.accepts(a);
    }
}

/// Filter alignments which has $(D flagname) flag set
final class FlagFilter(string flagname) : Filter {
    bool accepts(ref Alignment a) const {
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
    bool accepts(ref Alignment a) const {
        switch(_fieldname) {
            case "ref_id": mixin("return a.ref_id " ~ op ~ "_value;");
            case "position": mixin("return a.position " ~ op ~ "_value;");
            case "mapping_quality": mixin("return a.mapping_quality " ~ op ~ "_value;");
            case "sequence_length": mixin("return a.sequence_length " ~ op ~ "_value;");
            case "mate_ref_id": mixin("return a.next_ref_id " ~ op ~ "_value;");
            case "mate_position": mixin("return a.next_pos " ~ op ~ "_value;");
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
    bool accepts(ref Alignment a) const {
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

    bool accepts(ref Alignment a) const {
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
    bool accepts(ref Alignment a) const {
        switch(_fieldname) {
            case "read_name": mixin("return a.read_name " ~ op ~ " _value;");
            case "sequence": mixin("return cmp(a.sequence, _value) " ~ op ~ " 0;");
            case "cigar": mixin("return a.cigarString() " ~ op ~ " _value;");
            default: throw new Exception("unknown string field '" ~ _fieldname ~ "'");
        }
    }
}

/// Filtering string tags
final class StringTagFilter(string op) : Filter {
    private string _value;
    private string _tagname;

    this(string tagname, string value) {
        _tagname = tagname;
        _value = value;
    }

    bool accepts(ref Alignment a) const {
        auto v = a[_tagname];
        if (!v.is_string) {
            return false;
        }
        mixin(`return cast(string)v` ~ op ~ `_value;`);
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

    bool accepts(ref Alignment a) const {
        switch(_fieldname) {
            case "read_name": return !match(a.read_name, cast()_pattern).empty;
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

    bool accepts(ref Alignment a) const {
        auto v = a[_tagname];
        if (!v.is_string) {
            return false;
        }
        return !match(cast(string)v, cast()_pattern).empty;
    }
}
