import samfile;
import sam.serialize;
import std.stdio;

void main(string[] args) {

    auto sam = SamFile(args[1]);

    foreach (read; sam.alignments) {
        if (read.read_name != "") {
            writeln(to_sam(read, sam.reference_sequences));
        }
    }
}
