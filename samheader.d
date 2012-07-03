module samheader;

import reference;

import std.algorithm;
import std.conv;
import std.exception;
import std.array;
import utils.format;

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

    mixin template HeaderLineStruct(string struct_name, 
                                    string line_prefix,
                                    Field...) 
    {
         mixin(`struct `~struct_name~`{ 
                    mixin structFields!Field;
                    mixin parseStaticMethod!(struct_name, Field);
                    mixin serializeMethod!(line_prefix, Field);
                }`);
    }

}

mixin HeaderLineStruct!("HdLine", "@HD",
          Field!("format_version", "VN"),
          Field!("sorting_order", "SO"));

mixin HeaderLineStruct!("SqLine", "@SQ",
          Field!("name", "SN"),
          Field!("length", "LN", uint),
          Field!("assembly", "AS"),
          Field!("md5", "M5"),
          Field!("species", "SP"),
          Field!("uri", "UR"));

mixin HeaderLineStruct!("RgLine", "@RG",
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

mixin HeaderLineStruct!("PgLine", "@PG",
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

class HeaderLineDictionary(T, alias getID) {
    T opIndex(string id) const {
        return _dict[id];
    }

    void opIndexAssign(T line, string id) {
        _dict[id] = line;
    }

    bool add(T line) {
        auto id = getID(line);
        if (id !in _dict) {
            _dict[id] = line;
            return true;
        }
        return false;
    }

    bool remove(string id) {
        version(GNU) {
            bool result = (id in _dict) !is null;
            _dict.remove(id);
            return result;
        } else {
            return _dict.remove(id);
        }
    }

    int opApply(int delegate(ref T line) dg) {
        foreach (ref T line; _dict.byValue()) {
            auto res = dg(line);
            if (res != 0) {
                return res;
            }
        }
        return 0;
    }

    void clear() {
        _dict = null;
    }

    /// Returns: range of lines
    auto values() @property const {
        version(GNU) {
            return _dict.values;
        } else {
            return _dict.byValue();
        }
    }

    /// Returns: number of stored lines
    size_t length() @property const {
        return _dict.length;
    }

    protected T[string] _dict;
}

private string sqLineGetId(const ref SqLine line) { return line.name; }
private string rgLineGetId(const ref RgLine line) { return line.identifier; }
private string pgLineGetId(const ref PgLine line) { return line.identifier; }

/// Dictionary of @SQ lines.
class SqLineDictionary : HeaderLineDictionary!(SqLine, sqLineGetId)
{

    invariant() {
        import std.algorithm;
        assert(_indices.length == _dict.length);
        assert(_names.length == _dict.length);
    }

    override bool add(SqLine line) {
        if (super.add(line)) {
            _indices[line.name] = _names.length;
            _names ~= line.name;
            return true;
        }
        return false;
    }

    override bool remove(string sequence_name) {
        auto old_len = _dict.length;
        if (super.remove(sequence_name)) {

            auto index = _indices[sequence_name];
            _indices.remove(sequence_name); 

            for (size_t j = index + 1; j < old_len; ++j) {
                auto name = _names[j];
                _names[j - 1] = _names[j];
                _indices[_names[j - 1]] = j - 1;
            }

            _names.length = _names.length - 1;

            return true;
        }

        return false;
    }

    override void clear() {
        super.clear();
        _names.length = 0;
        _indices = null;
    }

    SqLine getSequence(size_t index) {
        return _dict[_names[index]];
    }

    int getSequenceIndex(string sequence_name) {
        size_t* ind = sequence_name in _indices;
        return (ind is null) ? -1 : cast(int)(*ind);
    }
        
    private:
        string[] _names;
        size_t[string] _indices;
}

/// Dictionary of @RG lines
alias HeaderLineDictionary!(RgLine, rgLineGetId) RgLineDictionary;

/// Dictionary of @PG lines
alias HeaderLineDictionary!(PgLine, pgLineGetId) PgLineDictionary;

class SamHeader {

    immutable DEFAULT_FORMAT_VERSION = "1.3";

    this() {
        sequences = new SqLineDictionary();
        read_groups = new RgLineDictionary();
        programs = new PgLineDictionary();

        format_version = DEFAULT_FORMAT_VERSION;
    }

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

    SqLine getSequence(size_t index) {
        return sequences.getSequence(index);
    }

}

/// Sorting order
enum SortingOrder {
    unknown,    ///
    unsorted,   ///
    coordinate, ///
    queryname   ///
}

string toSam(SamHeader header) {
    char[] buf;
    buf.reserve(65536);
    serialize(header, buf);
    return cast(string)buf;
}

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
