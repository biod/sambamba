module utils.version_;

immutable string VERSION = "0.6.7-pre1";

import bio.sam.header;
import std.array : join;

SamHeader addPG(string tool, string[] args, SamHeader header) {
    auto pg_line = PgLine();
    pg_line.identifier = "sambamba";
    pg_line.command_line = tool ~ " " ~ join(args[1 .. $], " ");
    pg_line.program_version = VERSION;

    if (header.programs.length > 0) {
        auto prev_id = header.programs.values.back.identifier;
        pg_line.previous_program = prev_id;
    }

    header.programs.add(pg_line);
    return header;
}
