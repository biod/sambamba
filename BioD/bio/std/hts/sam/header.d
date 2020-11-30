/*
    This file is part of BioD.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.std.hts.sam.header;

import bio.std.hts.thirdparty.msgpack;
import bio.core.utils.format;

import std.algorithm;
import std.conv;
import std.format;
import std.json;
import std.exception;
import std.array;
import std.range;
import std.traits;
import std.stdio : stderr;

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
          static if (t.FieldType.stringof == "string")
            result ~= `case "`~t.abbr~`":`~
                        `record.`~t.name~`=cast(string)(contents);`~
                      `break;`.dup;
          else
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
            str ~= `sink.write("\t` ~ Field[0].abbr ~ `:");`.dup;
            str ~= `sink.write(`~Field[0].name~`);`.dup;
            str ~= `}`.dup;
            return str.idup ~ serializeFields!(Field[1..$])();
        } else {
            return "";
        }
    }

    /*
        generates 'toSam' method which converts a struct
        to SAM header line
     */
    mixin template toSamMethod(string line_prefix, Field...) {
        void toSam(Sink)(auto ref Sink sink) const if (isSomeSink!Sink) {
            sink.write(line_prefix);
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
            hash_t toHash() const nothrow @safe{
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
                    mixin toSamMethod!(line_prefix, Field);
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
          Field!("uri", "UR"),
          Field!("alternate_haplotype", "AH"));

mixin HeaderLineStruct!("RgLine", "@RG", "identifier",
          Field!("identifier", "ID"),
          Field!("sequencing_center", "CN"),
          Field!("description", "DS"),
          Field!("date", "DT"),
          Field!("flow_order", "FO"),
          Field!("key_sequence", "KS"),
          Field!("library", "LB"),
          Field!("programs", "PG"),
          Field!("predicted_insert_size", "PI", int),
          Field!("platform", "PL"),
          Field!("platform_unit", "PU"),
          Field!("sample", "SM"),
          Field!("platform_model", "PM"));

mixin HeaderLineStruct!("PgLine", "@PG", "identifier",
          Field!("identifier", "ID"),
          Field!("name", "PN"),
          Field!("command_line", "CL"),
          Field!("previous_program", "PP"),
          Field!("program_version", "VN")); // version is a keyword in D

unittest {
    import std.algorithm;
    import std.stdio;

    // stderr.writeln("Testing @HD line parsing...");
    auto hd_line = HdLine.parse("@HD\tVN:1.0\tSO:coordinate");
    assert(hd_line.format_version == "1.0");
    assert(hd_line.sorting_order == "coordinate");

    // stderr.writeln("Testing @SQ line parsing...");
    auto sq_line = SqLine.parse("@SQ\tSN:NC_007605\tLN:171823\tM5:6743bd63b3ff2b5b8985d8933c53290a\tUR:ftp://.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz\tAS:NCBI37\tSP:HUMAN");
    assert(sq_line.name == "NC_007605");
    assert(sq_line.length == 171823);
    assert(sq_line.md5 == "6743bd63b3ff2b5b8985d8933c53290a");
    assert(sq_line.uri.endsWith("hs37d5.fa.gz"));
    assert(sq_line.assembly == "NCBI37");
    assert(sq_line.species == "HUMAN");

    // stderr.writeln("Testing @RG line parsing...");
    auto rg_line = RgLine.parse("@RG\tID:ERR016155\tLB:HUMgdtRAGDIAAPE\tSM:HG00125\tPI:488\tCN:BGI\tPL:ILLUMINA\tDS:SRP001294");
    assert(rg_line.identifier == "ERR016155");
    assert(rg_line.library == "HUMgdtRAGDIAAPE");
    assert(rg_line.sample == "HG00125");
    assert(rg_line.predicted_insert_size == 488);
    assert(rg_line.sequencing_center == "BGI");
    assert(rg_line.platform == "ILLUMINA");
    assert(rg_line.description == "SRP001294");

    // stderr.writeln("Testing @PG line parsing...");
    auto pg_line = PgLine.parse("@PG\tID:bam_calculate_bq\tPN:samtools\tPP:bam_recalibrate_quality_scores\tVN:0.1.17 (r973:277)\tCL:samtools calmd -Erb $bam_file $reference_fasta > $bq_bam_file");
    assert(pg_line.identifier == "bam_calculate_bq");
    assert(pg_line.name == "samtools");
    assert(pg_line.previous_program == "bam_recalibrate_quality_scores");
    assert(pg_line.program_version == "0.1.17 (r973:277)");
    assert(pg_line.command_line.endsWith("$bq_bam_file"));
}

// workaround for LDC bug #217
struct ValueRange(T) {
    this(T[string] dict, string[] ids) {
        _dict = dict;
        _ids = ids;
    }

    private {
        T[string] _dict;
        string[] _ids;
    }

    ref T front() @property { return _dict[_ids[0]]; }
    ref T back() @property { return _dict[_ids[$-1]]; }
    bool empty() @property { return _ids.length == 0; }
    void popFront() { _ids = _ids[1 .. $]; }
    void popBack() { _ids = _ids[0 .. $ - 1]; }
    ref T opIndex(size_t i) { return _dict[_ids[i]]; }
    size_t length() @property { return _ids.length; }
    ValueRange save() @property { return ValueRange(_dict, _ids[]); }
}

/// Common class for storing header lines
class HeaderLineDictionary(T) {

    invariant() {
        assert(_index_to_id.length == _dict.length);
        assert(_id_to_index.length == _dict.length);
        /*
        foreach(id, index; _id_to_index) {
            assert(_index_to_id[index] == id);
        }
        */
    }

    ///
    ref inout(T) opIndex(string id) inout {
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
    int opApply(scope int delegate(T line) dg) const {
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

    static assert(isRandomAccessRange!(ValueRange!T));

    /// Returns: range of lines
    ValueRange!T values() @property {
        return ValueRange!T(_dict, _index_to_id);
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
    ref inout(SqLine) getSequence(size_t index) inout {
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
    enum DEFAULT_FORMAT_VERSION = "1.3";

    /// Construct empty SAM header
    this() {
        sequences = new SqLineDictionary();
        read_groups = new RgLineDictionary();
        programs = new PgLineDictionary();

        format_version = DEFAULT_FORMAT_VERSION;
    }

    /// Parse SAM header given in plain text.
    this(string header_text) {
        read_groups = new RgLineDictionary();
        programs = new PgLineDictionary();
        format_version = DEFAULT_FORMAT_VERSION;

        import core.memory;
        core.memory.GC.disable();

        bool parsed_first_line = false;
        size_t n_sq_lines = 0;

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
            enforce(line[0] == '@', "Header lines must start with @");
            switch (line[1]) {
                case 'S':
                    enforce(line[2] == 'Q');
                    ++n_sq_lines;
                    break;
                case 'R':
                    enforce(line[2] == 'G');
                    auto rg_line = RgLine.parse(line);
                    if (!read_groups.add(rg_line)) {
                        stderr.writeln("duplicating @RG line ",  rg_line.identifier);
                    }
                    break;
                case 'P':
                    enforce(line[2] == 'G');
                    auto pg_line = PgLine.parse(line);
                    if (!programs.add(pg_line)) {
                        stderr.writeln("duplicating @PG line ", pg_line.identifier);
                    }
                    break;
                case 'H':
                    enforce(line[2] == 'D');
                    break;
                case 'C':
                    enforce(line[2] == 'O');
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

        _header_text = header_text;
        if (n_sq_lines <= 1000000)
            _parseSqLines(); // parse immediately for typical files

        core.memory.GC.enable();
    }

    /// Format version
    string format_version;

    /// Sorting order
    SortingOrder sorting_order = SortingOrder.unknown;

    /// Dictionary of @SQ lines.
    /// Removal is not allowed, you can only replace the whole dictionary.
    SqLineDictionary sequences() @property {
        if (_sequences is null)
            _parseSqLines();
        return _sequences;
    }

    void sequences(SqLineDictionary dict) @property {
        _sequences = dict;
    }

    private SqLineDictionary _sequences;
    private string _header_text;
    private void _parseSqLines() {
        import core.memory;
        core.memory.GC.disable();

        _sequences = new SqLineDictionary();

        foreach (line; splitter(_header_text, '\n')) {
            if (line.length < 3)
                continue;
            if (line[0 .. 3] != "@SQ")
                continue;

            auto sq_line = SqLine.parse(line);
            if (!_sequences.add(sq_line)) {
                stderr.writeln("duplicating @SQ line ",  sq_line.name);
            }
        }

        _header_text = null;
        core.memory.GC.enable();
    }

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

    /// Get header text representation in SAM format (includes trailing '\n')
    string text() @property {
        return to!string(this);
    }

    /// Header text representation in SAM ("%s") or JSON format ("%j").
    /// $(BR)
    /// Includes trailing '\n'.
    void toString(scope void delegate(const(char)[]) sink, FormatSpec!char fmt) {
        if (fmt.spec == 's')
            toSam(sink);
        else if (fmt.spec == 'j')
            toJson(sink);
        else
            throw new FormatException("unknown format specifier");
    }

    void toSam(Sink)(auto ref Sink sink) if (isSomeSink!Sink) {
        sink.write("@HD\tVN:");
        sink.write(format_version);
        if (sorting_order != SortingOrder.unknown) {
            sink.write("\tSO:");
            sink.write(to!string(sorting_order));
        }
        sink.write('\n');

        for (size_t i = 0; i < sequences.length; i++) {
            auto sq_line = getSequence(i);
            sq_line.toSam(sink);
            sink.write('\n');
        }

        foreach (rg_line; read_groups) {
            rg_line.toSam(sink);
            sink.write('\n');
        }

        foreach (pg_line; programs) {
            pg_line.toSam(sink);
            sink.write('\n');
        }

        foreach (comment; comments) {
            sink.write("@CO\t");
            sink.write(comment);
            sink.write('\n');
        }
    }

    void toJson(Sink)(auto ref Sink sink) if (isSomeSink!Sink) {
        JSONValue[string] result;

        result["format_version"] = format_version;

        if (sorting_order != SortingOrder.unknown) {
          result["sorting_order"] = sorting_order.to!string;
        }

        auto tmp = new JSONValue[sequences.length];
        for (auto i = 0; i < sequences.length; i++) {
          auto line = getSequence(i);
          JSONValue[string] sq;
          sq["sequence_name"] = line.name;
          sq["sequence_length"] = line.length;
          sq["assembly"] = line.assembly;
          sq["md5"] = line.md5;
          sq["species"] = line.species;
          sq["uri"] = line.uri;
          tmp[i].object = sq;
        }
        result["sq_lines"] = tmp.dup;
        tmp = null;

        auto tmp2 = new JSONValue[read_groups.length];
        foreach (i, line; read_groups) {
          JSONValue[string] sq;
          sq["identifier"] = line.identifier;
          sq["sequencing_center"] = line.sequencing_center;
          sq["description"] = line.description;
          sq["date"] = line.date;
          sq["flow_order"] = line.flow_order;
          sq["key_sequence"] = line.key_sequence;
          sq["library"] = line.library;
          sq["programs"] = line.programs;
          sq["predicted_insert_size"] = line.predicted_insert_size;
          sq["platform"] = line.platform;
          sq["platform_unit"] = line.platform_unit;
          sq["sample"] = line.sample;
          tmp2[i].object = sq;
        }
        result["rg_lines"] = tmp2;
        tmp2 = null;

        auto tmp3 = new JSONValue[programs.length];
        foreach (i, line; programs) {
          JSONValue[string] sq;
          sq["identifier"] = line.identifier;
          sq["program_name"] = line.name;
          sq["command_line"] = line.command_line;
          sq["previous_program"] = line.previous_program;
          sq["program_version"] = line.program_version;
          tmp3[i].object = sq;
        }
        result["pg_lines"] = tmp3;

        JSONValue json;
        json.object = result;
        static if (__VERSION__ < 2072)
            sink.write(toJSON(&json));
        else
          sink.write(toJSON(json,true));
    }

    /// Packs message in the following format:
    /// $(BR)
    /// MsgPack array with elements
    ///   $(OL
    ///     $(LI format version - string)
    ///     $(LI sorting order - string)
    ///     $(LI @SQ lines - array of dictionaries)
    ///     $(LI @RG lines - array of dictionaries)
    ///     $(LI @PG lines - array of dictionaries))
    /// $(BR)
    /// Dictionary keys are the same as in SAM format.
    void toMsgpack(Packer)(ref Packer packer) const {
        enforce(_sequences !is null, "failed to call msgpack");
        packer.beginArray(5);
        packer.pack(format_version);
        packer.pack(to!string(sorting_order));
        packer.beginArray(_sequences.length);
        foreach (sq; _sequences)
            packer.pack(sq);
        packer.beginArray(read_groups.length);
        foreach (rg; read_groups)
            packer.pack(rg);
        packer.beginArray(programs.length);
        foreach (pg; programs)
            packer.pack(pg);
    }
}

/// Sorting order
enum SortingOrder {
    unknown,    ///
    unsorted,   ///
    coordinate, ///
    queryname   ///
}

unittest {
    auto header = new SamHeader();
    import std.stdio;
    assert(header.text == "@HD\tVN:1.3\n");

    auto sequence = SqLine("abc", 123123);
    header.sequences.add(sequence);
    assert(header.text == "@HD\tVN:1.3\n@SQ\tSN:abc\tLN:123123\n");

    header.sorting_order = SortingOrder.coordinate;
    header.format_version = "1.2";
    assert(header.text == "@HD\tVN:1.2\tSO:coordinate\n@SQ\tSN:abc\tLN:123123\n");
    assert(header.getSequenceIndex("abc") == 0);
    assert(header.getSequenceIndex("bcd") == -1);

    header.sequences.clear();
    sequence = SqLine("bcd", 678);
    sequence.uri = "http://lorem.ipsum";
    header.sequences.add(sequence);
    header.format_version = "1.4";
    assert(header.text == "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:bcd\tLN:678\tUR:http://lorem.ipsum\n");

    header.sequences.add(SqLine("def", 321));
    assert(header.getSequenceIndex("abc") == -1);
    assert(header.getSequenceIndex("bcd") == 0);
    assert(header.getSequenceIndex("def") == 1);

    header.sequences.remove("bcd");
    assert(header.getSequenceIndex("abc") == -1);
    assert(header.getSequenceIndex("bcd") == -1);
    assert(header.getSequenceIndex("def") == 0);

    assert(header.text == "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:def\tLN:321\n");

    auto dict = new SqLineDictionary();
    dict.add(SqLine("yay", 111));
    dict.add(SqLine("zzz", 222));

    auto zzz = dict["zzz"];     // TODO: make 'dict["zzz"].uri = ...' work
    zzz.uri = "ftp://nyan.cat";
    dict["zzz"] = zzz;
    header.sequences = dict;

    assert(header.text ==
      "@HD\tVN:1.4\tSO:coordinate\n@SQ\tSN:yay\tLN:111\n@SQ\tSN:zzz\tLN:222\tUR:ftp://nyan.cat\n");
    assert(header.sequences == dict);

    header.sequences.remove("yay");
    header.sequences.remove("zzz");
    header.comments ~= "this is a comment";

    assert(header.text == "@HD\tVN:1.4\tSO:coordinate\n@CO\tthis is a comment\n");

    JSONValue[string] result;
    result["format_version"] = "1.2";
    assert(to!string(result) == "[\"format_version\":\"1.2\"]", to!string(result));
}
