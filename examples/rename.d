import bamfile;
import bamoutput;

import std.stream;
import std.range;
import std.typecons;

Alignment renamedAlignment(Tuple!(Alignment, int) ai) {
    ai[0].read_name = to!string(ai[1]);
    return ai[0];
}

void main(string[] args) {
    auto bam = BamFile(args[1]);

    auto alignments_with_indices = zip(bam.alignments, iota(1, int.max));
    auto alignments = map!renamedAlignment(alignments_with_indices);

    auto output = new BufferedFile("tmp.bam", FileMode.Out);
    scope(exit) output.close();

    writeBAM(output, bam.header.text, bam.reference_sequences, alignments, 0);
   
}
