import bamfile;

import std.datetime;
import std.stdio;

void main(string[] args) {
    StopWatch sw;
    sw.start();
    auto bam = BamFile(args[1]);
    auto count = 0;
    foreach (alignment; bam.alignments) {
        count += 1;
    }
    sw.stop();
    writeln("total time: ", sw.peek().nsecs, "ns");
    writeln("alignments found: ", count);

	import std.stream;
	import bamoutput;
	auto stream = new BufferedFile("tmp.bam", FileMode.Out);
	bam.rewind();
	writeBAM(stream, bam.header.text, bam.reference_sequences, bam.alignments, 0);
	stream.close();
}
