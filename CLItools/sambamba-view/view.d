module sambamba.view;

import bamfile;
import samfile;
import region;

import filtering;
import alignmentrangeprocessor;
import headerserializer;
import queryparser;

import utils.format;

import common.progressbar;

import std.stdio;
import std.c.stdio : stdout;
import std.array;
import std.traits;
import std.getopt;
import std.algorithm;

void printUsage() {
    writeln("Usage: sambamba-view [options] <input.bam | input.sam> [region1 [...]]");
    writeln();
    writeln("Options: -F, --filter=FILTER");
    writeln("                    set custom filter for alignments");
    writeln("         -f, --format=sam|bam|json");
    writeln("                    specify which format to use for output (default is SAM)");
    writeln("         -h, --with-header");
    writeln("                    print header before reads (always done for BAM output)");
    writeln("         -H, --header");
    writeln("                    output only header to stdout (if format=bam, the header is printed as SAM)");
    writeln("         -I, --reference-info");
    writeln("                    output to stdout only reference names and lengths in JSON");
    writeln("         -c, --count");
    writeln("                    output to stdout only count of matching records, hHI are ignored");
    writeln("         -v, --valid");
    writeln("                    output only valid alignments");
    writeln("         -S, --sam-input");
    writeln("                    specify that input is in SAM format");
    writeln("         -p, --show-progress");
    writeln("                    show progressbar in STDERR (works only for BAM files with no regions specified)");
    writeln("         -l, --compression-level");
    writeln("                    specify compression level (from 0 to 9, works only for BAM output)");
    writeln("         -o, --output-filename");
    writeln("                    specify output filename");
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

bool show_progress;

int compression_level = -1;
string output_filename;

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
               "sam-input|S",         &is_sam,
               "show-progress|p",     &show_progress,
               "compression-level|l", &compression_level,
               "output-filename|o",   &output_filename);
        
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
        writeln("sambamba-view: ", e.msg);

        version(development) {
            throw e; // rethrow to see detailed message
        }

        return 1;
    }
}

// TODO: mark pure functions/methods with 'pure' attribute
//       so that it becomes visible that accepts() is pure.
static __gshared Filter filter; 

bool passing(Alignment read) {
    return filter.accepts(read);
}

auto filtered(R)(R reads) {
    return std.algorithm.filter!passing(reads);
}

// In fact, $(D bam) is either BAM or SAM file
int sambambaMain(T)(T _bam, string[] args) 
    if (is(T == SamFile) || is(T == BamFile)) 
{

    immutable is_sam = is(T == SamFile);

    auto bam = _bam; // FIXME: uhm, that was a workaround for some closure-related bug

    if (reference_info_only && !count_only) {
        outputReferenceInfoJson(bam);
        return 0;
    }

    if ((with_header || header_only) && !count_only) {
        (new HeaderSerializer(format)).writeln(bam.header);
    }
    
    if (header_only) return 0;

    filter = new NullFilter();

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

    int processAlignments(AlignmentRangeProcessor)(AlignmentRangeProcessor processor) {
        static if (is(T == SamFile)) {
            if (args.length > 2) {
                stderr.writeln("sorry, accessing regions is unavailable for SAM input");
                return 1;
            }
        }

        if (args.length == 2) {
            
            static if (is(T == BamFile)) {
                if (show_progress) {
                    auto bar = new shared(ProgressBar)();
                    auto reads = bam.alignmentsWithProgress((lazy float p) { bar.update(p); });
                    processor.process(filtered(reads), bam);
                    bar.finish();
                } else {
                    processor.process(filtered(bam.alignments!withoutOffsets), bam);
                }
            } else { // SamFile
                processor.process(filtered(bam.alignments), bam);
            }
        } 

        // for BAM, random access is available
        static if (is(T == BamFile)) {
            if (args.length > 2) {
                auto regions = map!parseRegion(args[2 .. $]);

                alias ReturnType!(ReferenceSequence.opSlice) AlignmentRange;
                auto alignment_ranges = new AlignmentRange[regions.length];

                size_t i = 0;
                foreach (ref r; regions) {
                    alignment_ranges[i++] = bam[r.reference][r.beg .. r.end];
                }

                auto reads = joiner(alignment_ranges);
                processor.process(filtered(reads), bam);
            }
        }

        return 0;
    }

    if (count_only) {
        auto counter = new ReadCounter();

        if (processAlignments(counter))
            return 1;
        writeln(counter.number_of_reads);
    } else {
        switch (format) {
            case "bam":
                return processAlignments(new BamSerializer(output_filename, compression_level));
            case "sam":
                return processAlignments(new SamSerializer(output_filename));
            case "json":
                return processAlignments(new JsonSerializer(output_filename));
            default:
                stderr.writeln("output format must be one of sam, bam, json");
                return 1;
        }
    }

    return 0;
}