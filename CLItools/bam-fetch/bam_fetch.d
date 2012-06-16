import bamfile;
import region;
import sam.serialize;
import utils.format;

import std.stdio;
import std.c.stdio : stdout;
import std.array;
import std.getopt;
import std.algorithm;

void printUsage(string program) {
    writeln("Usage: " ~ program ~ " [options] <input.bam> region1 [...]");
    writeln();
    writeln("Options: -q, --quality-threshold=THRESHOLD");
    writeln("                    skip reads with mapping quality < THRESHOLD");
    writeln("         -r, --read-group=READGROUP");
    writeln("                    output only reads from read group READGROUP");
}

int quality_threshold = -1;
string read_group = null;

int main(string[] args) {
    try {

        getopt(args,
               "quality-threshold|q", &quality_threshold,
               "read-group|r",        &read_group);
        
        if (args.length < 3) {
            printUsage(args[0]);
            return 0;
        }

        auto bam = BamFile(args[1]);
        auto regions = map!parseRegion(args[2 .. $]);

		foreach (ref r; regions) {
			foreach (read; bam[r.reference][r.beg .. r.end]) {
				if (quality_threshold != -1 && read.mapping_quality < quality_threshold) {
					continue;
				}
				if (read_group !is null) {
					auto rg = read["RG"];
					if (!rg.is_string || (to!string(rg) != read_group)) {
						continue;
					}
				}
				serialize(read, bam.reference_sequences, stdout);
				putcharacter(stdout, '\n');
			}
		}
    } catch (Exception e) {
        writeln("bam-fetch: ", e.msg);
        return 1;
    }

    return 0;
}
