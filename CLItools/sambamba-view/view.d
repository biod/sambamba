module sambamba.view;

import bamfile;
import samfile;
import region;

import filter;
import serializer;
import queryparser;

import utils.format;

import std.stdio;
import std.c.stdio : stdout;
import std.array;
import std.getopt;
import std.algorithm;

void printUsage() {
    writeln("Usage: sambamba-view [options] <input.bam | input.sam> [region1 [...]]");
    writeln();
    writeln("Options: -F, --filter=FILTER");
    writeln("                    set custom filter for alignments");
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
    writeln("         -S, --sam-input");
    writeln("                    specify that input is in SAM format");
}

void outputReferenceInfoJson(T)(T bam) {
    
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

string format = "sam";
string query;
bool with_header;
bool header_only;
bool reference_info_only;
bool count_only;
bool skip_invalid_alignments;
bool is_sam;

version(standalone) {
    int main(string[] args) {
        return view_main(args);
    }
}

int view_main(string[] args) {
    try {

        getopt(args,
               std.getopt.config.caseSensitive,
               "filter|F",            &query,
               "format|f",            &format,
               "with-header|h",       &with_header,
               "header|H",            &header_only,
               "reference-info|I",    &reference_info_only,
               "count|c",             &count_only,
               "valid|v",             &skip_invalid_alignments,
               "sam-input|S",         &is_sam);
        
        if (args.length < 2) {
            printUsage();
            return 0;
        }

        if (!is_sam) {
            auto bam = BamFile(args[1]); 
            return sambambaMain(bam, args);
        } else {
            auto sam = SamFile(args[1]);
            return sambambaMain(sam, args);
        }
    } catch (Exception e) {
        writeln("sambamba: ", e.msg);

        version(development) {
            throw e; // rethrow to see detailed message
        }

        return 1;
    }
}

// In fact, $(D bam) is either BAM or SAM file
int sambambaMain(T)(T bam, string[] args) 
    if (is(T == SamFile) || is(T == BamFile)) 
{

    immutable is_sam = is(T == SamFile);

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

    if (skip_invalid_alignments) {
        filter = new AndFilter(filter, new ValidAlignmentFilter());
    }

    if (query !is null) {
        auto query_grammar = new QueryGrammar();
        auto node = query_grammar.parse(query);
        auto condition_node = cast(ConditionNode) node;
        if (condition_node is null) {
            stderr.writeln("filter string must represent a condition");
            return 1;
        }
        filter = new AndFilter(filter, condition_node.condition);
    }

    static string processAlignments(string s)() {
        return `
static if (is(T == SamFile)) {
        if (args.length > 2) {
            stderr.writeln("sorry, accessing regions is unavailable for SAM input");
            return 1;
        }
}
        if (args.length == 2) {
            foreach (read; bam.alignments) {
                if (filter.accepts(read)) ` ~ s ~ `;
            }
        } 

// for BAM, random access is available
static if (is(T == BamFile)) {
        if (args.length > 2) {
            auto regions = map!parseRegion(args[2 .. $]);

            foreach (ref r; regions) {
                foreach (read; bam[r.reference][r.beg .. r.end]) {
                    if (filter.accepts(read)) ` ~ s ~ `;
                }
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

    return 0;
}
