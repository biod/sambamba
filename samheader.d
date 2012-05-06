module samheader;

import std.algorithm;
import std.conv;

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
        /* certain assumptions about variable names are made here,
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
        return filter!"a[2]==':'"(splitter(header_line[3..$], '\t'));
    }

    mixin template parseStaticMethod(string struct_name, Field...) {

        static auto parse(string line) {
            mixin(struct_name ~ " record;");
            foreach (field; fields(line)) {
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
          Field!("program_version", "VN")); // 'version' is reserved word in D 

struct SamHeader {
public:

    this(string header) {
        this.header = header;
        parse();
    }

    SqLine[] sqLines() @property { return _sq_lines; }
    RgLine[] rgLines() @property { return _rg_lines; }
    PgLine[] pgLines() @property { return _pg_lines; }

private:
    string header;

    SqLine[] _sq_lines;
    RgLine[] _rg_lines;
    PgLine[] _pg_lines;

    void parse() {
        foreach (line; splitter(header, '\n')) {
            switch (line[0..3]) {
                case "@SQ":
                    _sq_lines ~= SqLine.parse(line);
                    break;
                case "@RG":
                    _rg_lines ~= RgLine.parse(line);
                    break;
                case "@PG":
                    _pg_lines ~= PgLine.parse(line);
                    break;
                default:
                    /* NYI */
                    break;
            }
        }
    }
}
