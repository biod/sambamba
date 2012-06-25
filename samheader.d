module samheader;

import std.algorithm;
import std.conv;
import std.exception;
import reference;

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
                    mixin(makeSwitchStatements!(Field)());
                }
            }
            return record;
        }
    }

    mixin template HeaderLineStruct(string struct_name, Field...) {
         mixin(`struct `~struct_name~`{ 
                    mixin structFields!Field;
                    mixin parseStaticMethod!(struct_name, Field);
                }`);
    }

}

mixin HeaderLineStruct!("HdLine",
          Field!("format_version", "VN"),
          Field!("sorting_order", "SO"));

mixin HeaderLineStruct!("SqLine", 
          Field!("sequence_name", "SN"),
          Field!("sequence_length", "LN", uint),
          Field!("assembly", "AS"),
          Field!("md5", "M5"),
          Field!("species", "SP"),
          Field!("uri", "UR"));

mixin HeaderLineStruct!("RgLine",
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

mixin HeaderLineStruct!("PgLine",
          Field!("identifier", "ID"),
          Field!("program_name", "PN"),
          Field!("command_line", "CL"),
          Field!("previous_program", "PP"),
          Field!("program_version", "VN"));

unittest {
    import std.algorithm;
    import std.stdio;

    writeln("Testing @HD line parsing...");
    auto hd_line = HdLine.parse("@HD\tVN:1.0\tSO:coordinate");
    assert(hd_line.format_version == "1.0");
    assert(hd_line.sorting_order == "coordinate");

    writeln("Testing @SQ line parsing...");
    auto sq_line = SqLine.parse("@SQ\tSN:NC_007605\tLN:171823\tM5:6743bd63b3ff2b5b8985d8933c53290a\tUR:ftp://.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz\tAS:NCBI37\tSP:HUMAN");
    assert(sq_line.sequence_name == "NC_007605");
    assert(sq_line.sequence_length == 171823);
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
    assert(pg_line.program_name == "samtools");
    assert(pg_line.previous_program == "bam_recalibrate_quality_scores");
    assert(pg_line.program_version == "0.1.17 (r973:277)");
    assert(pg_line.command_line.endsWith("$bq_bam_file"));
}

/*
   The structure representing SAM header.

   SqLine, RgLine, and PgLine structs represent
   types of header lines described in the
   SAM/BAM specification. They are pretty much
   self-documented in the code but since they are
   declared via mixins, DDoc won't recognize them.
 */
struct SamHeader {
public:

    /* 
       Constructor taking header contents as a string
     */
    this(string header) {
        _header = header;
        parse();
    }

    /// Returns: whether there is a @HD line in the header
    bool hasHeaderLine() const { return _header_line != HdLine.init; }

    /// Returns: index of $(D name) in $(D sq_lines) or -1 if not found.
    int getReferenceSequenceId(string name) const { 
        const(size_t)* p_id = name in _sequences;
        if (p_id is null) return -1;
        return cast(int)*p_id;
    }

    /// Returns: array of SqLine structs representing @SQ lines
    SqLine[] sq_lines() @property { return _sq_lines; }

    /// Returns: array of RgLine structs representing @RG lines
    RgLine[] rg_lines() @property { return _rg_lines; }

    /// Returns: array of PgLine structs representing @PG lines
    PgLine[] pg_lines() @property { return _pg_lines; }

    /// Returns: format version if present in header, or null
    string format_version() @property const { return _header_line.format_version; }

    /// Returns: sorting order ('unknown', 'unsorted',
    ///          'queryname', or 'coordinate')
    string sorting_order() @property const { 
        if (!hasHeaderLine()) {
            return "unknown";
        } else {
            return _header_line.sorting_order; 
        }
    }

    /// Returns: urls of all fasta files encountered in @SQ lines
    string[] fasta_urls() @property { return _fasta_urls; }

    /// Returns: raw text of the header
    string text() @property const {
        return _header;
    }

private:
    string _header;

    HdLine _header_line;

    SqLine[] _sq_lines;
    RgLine[] _rg_lines;
    PgLine[] _pg_lines;

    string[] _fasta_urls;

    size_t[string] _sequences;

    void parse() {
        bool parsed_first_line = false;

        uint[string] _fasta_urls_dict; /// acts like a set

        foreach (line; splitter(_header, '\n')) {
            if (line.length < 3) {
                continue;
            }
            if (!parsed_first_line && line[0..3] == "@HD") {
                /* parse header line */
                _header_line = HdLine.parse(line);
            }
            switch (line[0..3]) {
                case "@SQ":
                    auto sq_line = SqLine.parse(line);
                    // set unique integer identifier for this sequence
                    _sequences[sq_line.sequence_name] = sq_lines.length;
                    _sq_lines ~= sq_line;
                    if (sq_line.uri != null) {
                        _fasta_urls_dict[sq_line.uri] = 1;
                    }
                    break;
                case "@RG":
                    _rg_lines ~= RgLine.parse(line);
                    break;
                case "@PG":
                    _pg_lines ~= PgLine.parse(line);
                    break;
                case "@HD":
                case "@CO":
                    break;
                default:
                    assert(0);
            }

            parsed_first_line = true;
        }

        _fasta_urls = _fasta_urls_dict.keys;
    }
}
