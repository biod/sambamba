import samfile;
import sam.serialize;
import utils.format;
import std.c.stdio : stdout;

void main(string[] args) {

    auto sam = SamFile(args[1]);

    foreach (read; sam.alignments) {
        if (read.read_name != "") {
            serialize(read, sam.reference_sequences, stdout);
            putcharacter(stdout, '\n');
        }
    }
}
