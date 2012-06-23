import bamfile;
import region;

import filter;
import serializer;

import utils.format;

import std.stdio;
import std.c.stdio : stdout;
import std.array;
import std.getopt;
import std.algorithm;

void printUsage(string program) {
    writeln("Usage: " ~ program ~ " [options] <input.bam> [region1 [...]]");
    writeln();
    writeln("Options: -q, --quality-threshold=THRESHOLD");
    writeln("                    skip reads with mapping quality < THRESHOLD");
    writeln("         -r, --read-group=READGROUP");
    writeln("                    output only reads from read group READGROUP");
    writeln("         -f, --format=sam|json");
    writeln("                    specify which format to use for output (default is SAM)");
    writeln("         -h, --with-header");
    writeln("                    print header before reads");
    writeln("         -H, --header");
    writeln("                    output only header");
    writeln("         -I, --reference-info");
    writeln("                    output only reference names and lengths in JSON");
    writeln("         -c, --count");
    writeln("                    output only count of matching records, hHI are ignored");
    writeln("         -v, --valid");
    writeln("                    output only valid alignments");
}

void outputReferenceInfoJson(ref BamFile bam) {
    
    putcharacter(stdout, '[');

    bool first = true;
    foreach (refseq; bam.reference_sequences) {
        if (first) {
            first = false;
        } else {
            putcharacter(stdout, ',');
        }
        putstring(stdout, `{"name":"`);
        foreach (char c; refseq.name) {
            if (c == '\\' || c == '"')
                putcharacter(stdout, '\\');
            putcharacter(stdout, c);
        }
        putstring(stdout, `","length":`);
        putinteger(stdout, refseq.length);
        putcharacter(stdout, '}');
    }

    putcharacter(stdout, ']');
    putcharacter(stdout, '\n');
}

ubyte quality_threshold = 0;
string read_group = null;
string format = "sam";
bool with_header;
bool header_only;
bool reference_info_only;
bool count_only;
bool skip_invalid_alignments;

int main(string[] args) {
    try {

        getopt(args,
               std.getopt.config.caseSensitive,
               "quality-threshold|q", &quality_threshold,
               "read-group|r",        &read_group,
               "format|f",            &format,
               "with-header|h",       &with_header,
               "header|H",            &header_only,
               "reference-info|I",    &reference_info_only,
               "count|c",             &count_only,
               "valid|v",             &skip_invalid_alignments);
        
        if (args.length < 2) {
            printUsage(args[0]);
            return 0;
        }

        auto bam = BamFile(args[1]);

        auto serializer = new Serializer(format);

        if (reference_info_only && !count_only) {
            outputReferenceInfoJson(bam);
            return 0;
        }

        if ((with_header || header_only) && !count_only) {
            serializer.writeln(bam.header);
        }
        
        if (header_only) return 0;

        Filter filter = new NullFilter();

        if (quality_threshold != 0) {
            filter = new AndFilter(filter, 
                                   new MappingQualityFilter(quality_threshold));
        }

        if (read_group !is null) {
            filter = new AndFilter(filter, 
                                   new ReadGroupFilter(read_group));
        }

        if (skip_invalid_alignments) {
            filter = new AndFilter(filter,
                                   new ValidAlignmentFilter());
        }

        static string processAlignments(string s)() {
            return `
            if (args.length == 2) {
                foreach (read; bam.alignments) {
                    if (filter.accepts(read)) ` ~ s ~ `;
                }
            } else {
                auto regions = map!parseRegion(args[2 .. $]);

                foreach (ref r; regions) {
                    foreach (read; bam[r.reference][r.beg .. r.end]) {
                        if (filter.accepts(read)) ` ~ s ~ `;
                    }
                }
            }`;
        }

        if (count_only) {
            uint count;
            mixin(processAlignments!"count += 1"());
            writeln(count);
        } else {
            mixin(processAlignments!"serializer.writeln(read, bam.reference_sequences)"());
        }

        
    } catch (Exception e) {
        writeln("sambamba: ", e.msg);
        return 1;
    }

    return 0;
}
