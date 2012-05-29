import bamfile;
import bai.read;
import bgzfrange;
import virtualoffset;

import std.stdio;

void main(string[] args) {
	auto bf = BamFile(args[1]);
	auto index = BaiFile(args[1] ~ ".bai");

/*	write(bf[VirtualOffset(2514054033992)].read_name);
*/
	foreach (ind; index.indices) {
		foreach (bin; ind.bins) {
            foreach (chunk; bin.chunks) {
                writeln(bin.level, "\t", chunk.beg, "\t", chunk.end);
            }
		}
	}

    //foreach(al; bf.fetch(19, 1916900, 1917000)) {
    //    writeln(al.position);
    //    writeln(al.sequence_length);
   // }
}

