import bamfile;
import bai.read;
import std.stdio;
import std.algorithm;
import sam.serialize;

void main(string[] args) {
    auto bai = BaiFile(args[1]);
    auto bam = BamFile(args[2]);
    foreach (size_t i, index; bai.indices) {
        writeln("Reference sequence ", bam.reference_sequences[i].name,
                " (length: ", bam.reference_sequences[i].length, "):");
        writeln("\tBins:");
        foreach (bin; index.bins) {
            writeln("\t\t", bin.id);
            writeln("\t\tChunks:");
            foreach (chunk; bin.chunks) {
                writeln("\t\t\t[", chunk.beg, ", ", chunk.end, ")");
                try {
                writeln("\t\t\t", toSam(bam.getAlignmentAt(chunk.beg),
                                        bam.reference_sequences)[0..min(40, $)], "...");
                } catch (Throwable e) {
                    writeln("\t\t\t[error]...");
                }
                writeln();
            }
        }

        writeln("\tLinear index:");
        foreach (entry; index.ioffsets) {
            writeln("\t\t", entry);
        }
        writeln();
    }
}
