/*
    This file is part of Sambamba.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

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
module samheader;

import reference;

import std.algorithm;
import std.conv;
import std.exception;
import std.array;
import std.range;
import std.traits;
import utils.format;
import utils.msgpack;

import std.stdio;

private {

    struct Field(string _name, string _abbr, T = string) {
        enum name = _name;
        enum abbr = _abbr;
        alias T FieldType;
    }

    mixin template structFields(T...) {
        static if (!(T.length == 0)) {
            mixin(T[0].FieldType.stringof ~ " " ~  T[0].name ~ ";");
            mixin structFields!(T[1..$]);
        }
    }

    string makeSwitchStatements(F...)() {
        /* certain assumptions about variable names are being made here,
         * namely, 'record' and 'contents'
         */
        char[] result;
        foreach (t; F) {
            result ~= `case "`~t.abbr~`":`~
                        `record.`~t.name~`=to!(`~t.FieldType.stringof~`)(contents);`~
                      `break;`.dup;
        }
        result ~= `default: break;`.dup;
        return cast(string)result;
    }

    auto fields(string header_line) {
        return splitter(header_line[3..$], '\t');
    }

    /*
        generates 'parse' method which parses given string and
        fills corresponding struct fields
     */
    mixin template parseStaticMethod(string struct_name, Field...) {

        static auto parse(string line) {
            mixin(struct_name ~ " record;");
            foreach (field; fields(line)) {
                if (field.length < 3) {
                    continue;
                }
                if (field[2] != ':') {
                    continue;
                }
                string contents = field[3..$];
                switch (field[0..2]) {
                    mixin(makeSwitchStatements!Field());
                }
            }
            return record;
        }
    }

    string serializeFields(Field...)() {
        static if (Field.length > 0) {
            char[] str = `if (`~Field[0].name~` != `~Field[0].FieldType.stringof~`.init) {`.dup;
            str ~= `putstring(stream, "\t` ~ Field[0].abbr ~ `:");`.dup;
            if (Field[0].FieldType.stringof == "string") {
                str ~= `putstring(stream, `~Field[0].name~`);`.dup;
            } else {
                str ~= `putinteger(stream, `~Field[0].name~`);`.dup;
            }
            str ~= `}`.dup;
            return str.idup ~ serializeFields!(Field[1..$]);
        } else {
            return "";
        }
    }

    /*
        generates 'serialize' method which converts a struct
        to SAM header line
     */
    mixin template serializeMethod(string line_prefix, Field...) {
        void serialize(S)(ref S stream) const {
            putstring(stream, line_prefix);
            mixin(serializeFields!Field());    
        }
    }

    string generateHashExpression(Field...)() {
        char[] res;
        foreach (t; Field) {
            res ~= "result = 31 * result + " ~
                   "typeid(" ~ t.name ~ ").getHash(&" ~ t.name ~ ");".dup;
        }
        return res.idup;
    }

    mixin template toHashMethod(string struct_name, string id_field, Field...) {
        static if (id_field != null) {
            hash_t toHash() const {
                hash_t result = 1;
                mixin(generateHashExpression!Field());    
                return result;
            }

            mixin("int opCmp(const ref " ~ struct_name ~ " other) " ~
                  "    const pure nothrow @safe" ~
                  "{" ~
                  "    return " ~ id_field ~ " < other." ~ id_field ~ " ? -1 : " ~
                  "           " ~ id_field ~ " > other." ~ id_field ~ " ? 1 : 0;" ~
                  "}");
        }
    }

    string opEqualsExpression(Field...)() {
        char[] result = Field[0].name ~ " == other.".dup ~ Field[0].name;
        foreach (t; Field[1..$]) {
            result ~= " && " ~ t.name ~ " == other.".dup ~ t.name;
        }
        return result.idup;
    }

    mixin template opEqualsMethod(string struct_name, Field...) {
        mixin("bool opEquals(const ref " ~ struct_name ~ " other)" ~
              "    pure const @safe nothrow" ~
              "{" ~
              "    return " ~ opEqualsExpression!Field() ~ ";" ~
              "}");

        mixin("bool opEquals(" ~ struct_name ~ " other)" ~
              "    pure const @safe nothrow" ~
              "{" ~
              "    return " ~ opEqualsExpression!Field() ~ ";" ~
              "}");
    }

    mixin template getSetIDMethods(string id_field) {
        static if (id_field != null) {
            auto getID() const pure nothrow @safe {
                mixin("return " ~ id_field ~";");
            }

            mixin("void setID(typeof("~id_field~") id) pure nothrow @safe { " ~ id_field ~ " = id; }");
        }
    }

    string generateToMsgpackMethod(Field...)() {
        char[] method = "packer.beginMap(" ~ to!string(Field.length) ~ ");".dup;
        foreach (t; Field) {
            method ~= "packer.pack(`" ~ t.abbr ~ "`);".dup;
            method ~= "packer.pack(" ~ t.name ~ ");".dup;
        }
        return method.idup;
    }

    mixin template toMsgpackMethod(Field...) {

        void toMsgpack(Packer)(ref Packer packer) const {
            mixin(generateToMsgpackMethod!Field());
        }
    }

    mixin template HeaderLineStruct(string struct_name, 
                                    string line_prefix,
                                    string id_field,
                                    Field...) 
    {
         mixin(`struct `~struct_name~`{ 
                    mixin structFields!Field;
                    mixin parseStaticMethod!(struct_name, Field);
                    mixin serializeMethod!(line_prefix, Field);
                    mixin toHashMethod!(struct_name, id_field, Field);
                    mixin opEqualsMethod!(struct_name, Field);
                    mixin getSetIDMethods!id_field;
                    mixin toMsgpackMethod!Field;
                }`);
    }

}

mixin HeaderLineStruct!("HdLine", "@HD", null,
          Field!("format_version", "VN"),
          Field!("sorting_order", "SO"));

mixin HeaderLineStruct!("SqLine", "@SQ", "name",
          Field!("name", "SN"),
          Field!("length", "LN", uint),
          Field!("assembly", "AS"),
          Field!("md5", "M5"),
          Field!("species", "SP"),
          Field!("uri", "UR"));

mixin HeaderLineStruct!("RgLine", "@RG", "identifier",
          Field!("identifier", "ID"),
          Field!("sequencing_center", "CN"),
          Field!("description", "DS"),
          Field!("date", "DT"),
          Field!("flow_order", "FO"),
          Field!("key_sequence", "KS"),
          Field!("library", "LB"),
          Field!("programs", "PG"),
          Field!("predicted_insert_size", "PI", uint),
          Field!("platform", "PL"),
          Field!("platform_unit", "PU"),
          Field!("sample", "SM"));

mixin HeaderLineStruct!("PgLine", "@PG", "identifier",
          Field!("identifier", "ID"),
          Field!("name", "PN"),
          Field!("command_line", "CL"),
          Field!("previous_program", "PP"),
          Field!("program_version", "VN")); // version is a keyword in D

unittest {
    import std.algorithm;
    import std.stdio;

    writeln("Testing @HD line parsing...");
    auto hd_line = HdLine.parse("@HD\tVN:1.0\tSO:coordinate");
    assert(hd_line.format_version == "1.0");
    assert(hd_line.sorting_order == "coordinate");

    writeln("Testing @SQ line parsing...");
    auto sq_line = SqLine.parse("@SQ\tSN:NC_007605\tLN:171823\tM5:6743bd63b3ff2b5b8985d8933c53290a\tUR:ftp://.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz\tAS:NCBI37\tSP:HUMAN");
    assert(sq_line.name == "NC_007605");
    assert(sq_line.length == 171823);
    assert(sq_line.md5 == "6743bd63b3ff2b5b8985d8933c53290a");
    assert(sq_line.uri.endsWith("hs37d5.fa.gz"));
    assert(sq_line.assembly == "NCBI37");
    assert(sq_line.species == "HUMAN");

    writeln("Testing @RG line parsing...");
    auto rg_line = RgLine.parse("@RG\tID:ERR016155\tLB:HUMgdtRAGDIAAPE\tSM:HG00125\tPI:488\tCN:BGI\tPL:ILLUMINA\tDS:SRP001294");
    assert(rg_line.identifier == "ERR016155");
    assert(rg_line.library == "HUMgdtRAGDIAAPE");
    assert(rg_line.sample == "HG00125");
    assert(rg_line.predicted_insert_size == 488);
    assert(rg_line.sequencing_center == "BGI");
    assert(rg_line.platform == "ILLUMINA");
    assert(rg_line.description == "SRP001294");

    writeln("Testing @PG line parsing...");
    auto pg_line = PgLine.parse("@PG\tID:bam_calculate_bq\tPN:samtools\tPP:bam_recalibrate_quality_scores\tVN:0.1.17 (r973:277)\tCL:samtools calmd -Erb $bam_file $reference_fasta > $bq_bam_file");
    assert(pg_line.identifier == "bam_calculate_bq");
    assert(pg_line.name == "samtools");
    assert(pg_line.previous_program == "bam_recalibrate_quality_scores");
    assert(pg_line.program_version == "0.1.17 (r973:277)");
    assert(pg_line.command_line.endsWith("$bq_bam_file"));
}

/// Common class for storing header lines
class HeaderLineDictionary(T) {

    /* D doesn't currently support invariant in this class 
       because values() has type auto :-(

    invariant() {
        assert(_index_to_id.length == _dict.length);
        assert(_id_to_index.length == _dict.length);
        foreach(id, index; _id_to_index) {
            assert(_index_to_id[index] == id);
        }
    }
    */

    ///
    T opIndex(string id) const {
        return _dict[id];
    }

    ///
    void opIndexAssign(T line, string id) {
        _dict[id] = line;
    }

    ///
    const(T)* opIn_r(string id) const {
        return id in _dict;
    }

    /// Append a line
    bool add(T line) {
        auto id = line.getID();
        if (id !in _dict) {
            _dict[id] = line;
            _id_to_index[id] = _index_to_id.length;
            _index_to_id ~= id;
            return true;
        }
        return false;
    }

    /// Remove a line with identifier $(D id).
    bool remove(string id) {
        if (id in _dict) {
            auto old_len = _dict.length;

            for (auto j = _id_to_index[id]; j < old_len - 1; ++j) {
                _index_to_id[j] = _index_to_id[j + 1];
                _id_to_index[_index_to_id[j]] = j;
            }

            _index_to_id.length = _index_to_id.length - 1;

            _dict.remove(id);
            _id_to_index.remove(id); 

            return true;
        }

        return false;
    }

    ///
    int opApply(scope int delegate(ref T line) dg) {
        foreach (size_t i; 0 .. _dict.length) {
            auto res = dg(_dict[_index_to_id[i]]);
            if (res != 0) {
                return res;
            }
        }
        return 0;
    }

    ///
    int opApply(scope int delegate(ref size_t index, ref T line) dg) {
        foreach (size_t i; 0 .. _dict.length) {
            auto res = dg(i, _dict[_index_to_id[i]]);
            if (res != 0) {
                return res;
            }
        }
        return 0;
    }

    /// Clear the dictionary
    void clear() {
        _dict = null;
        _id_to_index = null;
        _index_to_id.length = 0;
    }

    /// Returns: range of lines
    auto values() @property {
        return map!((size_t i) {
                        return _dict[_index_to_id[i]];
                    })(iota(_dict.length));
    }

    /// Returns: number of stored lines
    size_t length() @property const {
        return _dict.length;
    }

    protected {
        T[string] _dict;
        string[] _index_to_id;
        size_t[string] _id_to_index;
    }
}

/// Dictionary of @SQ lines.
final class SqLineDictionary : HeaderLineDictionary!SqLine
{
    ///
    SqLine getSequence(size_t index) {
        return _dict[_index_to_id[index]];
    }

    ///
    int getSequenceIndex(string sequence_name) {
        size_t* ind = sequence_name in _id_to_index;
        return (ind is null) ? -1 : cast(int)(*ind);
    }
}

/// Dictionary of @RG lines
alias HeaderLineDictionary!RgLine RgLineDictionary;

/// Dictionary of @PG lines
alias HeaderLineDictionary!PgLine PgLineDictionary;

/// Represents SAM header
class SamHeader {

    ///
    immutable DEFAULT_FORMAT_VERSION = "1.3";

    /// Construct empty SAM header
    this() {
        sequences = new SqLineDictionary();
        read_groups = new RgLineDictionary();
        programs = new PgLineDictionary();

        format_version = DEFAULT_FORMAT_VERSION;
    }

    /// Parse SAM header given in plain text.
    this(string header_text) {
        this();
        bool parsed_first_line = false;

        foreach (line; splitter(header_text, '\n')) {
            if (line.length < 3) {
                continue;
            }
            if (!parsed_first_line && line[0..3] == "@HD") {
                auto header_line = HdLine.parse(line);
                if (header_line.sorting_order.length > 0) {
                    try {
                        sorting_order = to!SortingOrder(header_line.sorting_order);
                    } catch (ConvException e) {
                        sorting_order = SortingOrder.unknown; 
                        // FIXME: should we do that silently?
                    }
                } else {
                    sorting_order = SortingOrder.unknown;
                }
                format_version = header_line.format_version;
            }
            switch (line[0..3]) {
                case "@SQ":
                    auto sq_line = SqLine.parse(line);
                    if (!sequences.add(sq_line)) {
                        stderr.writeln("duplicating @SQ line ",  sq_line.name);
                    }
                    break;
                case "@RG":
                    auto rg_line = RgLine.parse(line);
                    if (!read_groups.add(rg_line)) {
                        stderr.writeln("duplicating @RG line ",  rg_line.identifier);
                    }
                    break;
                case "@PG":
                    auto pg_line = PgLine.parse(line);
                    if (!programs.add(pg_line)) {
                        stderr.writeln("duplicating @PG line ", pg_line.identifier);
                    }
                    break;
                case "@HD":
                    break;
                case "@CO":
                    comments ~= line[4..$];
                    break;
                default:
                    assert(0);
            }

            parsed_first_line = true;
        }

        if (!parsed_first_line) {
            format_version = DEFAULT_FORMAT_VERSION;
        }
    }
       
    /// Format version
    string format_version;

    /// Sorting order
    SortingOrder sorting_order = SortingOrder.unknown;

    /// Dictionary of @SQ lines. 
    /// Removal is not allowed, you can only replace the whole dictionary.
    SqLineDictionary sequences;

    /// Dictionary of @RG lines
    RgLineDictionary read_groups;

    /// Dictionary of @PG lines
    PgLineDictionary programs;

    /// Array of @CO lines
    string[] comments;

    /// Zero-based index of sequence.
    /// If such sequence does not exist in the header, returns -1.
    int getSequenceIndex(string sequence_name) {
        return sequences.getSequenceIndex(sequence_name);
    }

    ///
    SqLine getSequence(size_t index) {
        return sequences.getSequence(index);
    }

    /// Get header text representation
    string text() @property {
        return toSam(this);
    }

}

/// Sorting order
enum SortingOrder {
    unknown,    ///
    unsorted,   ///
    coordinate, ///
    queryname   ///
}

/// Get SAM representation of $(D header)
string toSam(SamHeader header) {
    char[] buf;
    buf.reserve(65536);
    serialize(header, buf);
    return cast(string)buf;
}

/// Serialize $(D header) to $(D stream)
void serialize(S)(SamHeader header, ref S stream) {
    putstring(stream, "@HD\tVN:");
    putstring(stream, header.format_version);
    if (header.sorting_order != SortingOrder.unknown) {
        putstring(stream, "\tSO:");
        putstring(stream, to!string(header.sorting_order));
    }
    putcharacter(stream, '\n');
   
    for (size_t i = 0; i < header.sequences.length; i++) {
        auto sq_line = header.getSequence(i);
        sq_line.serialize(stream);
        putcharacter(stream, '\n');
    }

    foreach (rg_line; header.read_groups) {
        rg_line.serialize(stream);
        putcharacter(stream, '\n');
    }

    foreach (pg_line; header.programs) {
        pg_line.serialize(stream);
        putcharacter(stream, '\n');
    }

    foreach (comment; header.comments) {
        putstring(stream, "@CO\t");
        putstring(stream, comment);
        putcharacter(stream, '\n');
    }
}

unittest {
    auto header = new SamHeader();
    import std.stdio;
    assert(toSam(header) == "@HD\tVN:1.3\n");

    auto sequence = SqLine("abc", 123123);
    header.sequences.add(sequence);
    assert(toSam(header) == "@HD\tVN:1.3\n@SQ\tSN:abc\tLN:123123\n");

    header.sorting_order = SortingOrder.coordinate;
    header.format_version = "1.2";
    assert(toSam(header) == "@HD\tVN:1.2\tSO:coordinate\n@SQ\tSN:abc\tLN:123123\n");
    assert(header.getSequenceIndex("abc") == 0);
    assert(header.getSequenceIndex("bcd") == -1);

    header.sequences.clear();
    sequence = SqLine("bcd", 678);
    sequence.uri = "http://lorem.ipsum";
    header.sequences.add(sequence);
    header.format_version = "1.4";
    assert(toSam(header) == "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:bcd\tLN:678\tUR:http://lorem.ipsum\n");

    header.sequences.add(SqLine("def", 321));
    assert(header.getSequenceIndex("abc") == -1);
    assert(header.getSequenceIndex("bcd") == 0);
    assert(header.getSequenceIndex("def") == 1);

    header.sequences.remove("bcd");
    assert(header.getSequenceIndex("abc") == -1);
    assert(header.getSequenceIndex("bcd") == -1);
    assert(header.getSequenceIndex("def") == 0);

    assert(toSam(header) == "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:def\tLN:321\n");

    auto dict = new SqLineDictionary();
    dict.add(SqLine("yay", 111));
    dict.add(SqLine("zzz", 222));

    auto zzz = dict["zzz"];     // TODO: make 'dict["zzz"].uri = ...' work
    zzz.uri = "ftp://nyan.cat";
    dict["zzz"] = zzz;
    header.sequences = dict;

    assert(toSam(header) == 
      "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:yay\tLN:111\n@SQ\tSN:zzz\tLN:222\tUR:ftp://nyan.cat\n");
    assert(header.sequences == dict);

    header.sequences.remove("yay");
    header.sequences.remove("zzz");
    header.comments ~= "this is a comment";

    assert(toSam(header) == "@HD\tVN:1.4\tSO:coordinate\n@CO\tthis is a comment\n");
}
