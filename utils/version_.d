module utils.version_;

immutable string VERSION = "0.6.8-pre1";
immutable string HEADER_VERSION = "1.0"; // goes in header for reproducibility between sambamba versions

import bio.sam.header;
import std.array : join;

SamHeader addPG(string tool, string[] args, SamHeader header) {
    auto pg_line = PgLine();
    pg_line.identifier = "sambamba";
    pg_line.command_line = tool ~ " " ~ join(args[1 .. $], " ");
    pg_line.program_version = HEADER_VERSION;

    if (header.programs.length > 0) {
        auto prev_id = header.programs.values.back.identifier;
        pg_line.previous_program = prev_id;
    }

    header.programs.add(pg_line);
    return header;
}
