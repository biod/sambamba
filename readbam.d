import bamfile;
import std.stdio;
void main(string[] args) {
	auto bam = BamFile(args[1]);
	foreach (alignment; bam.alignments) {
		foreach (k, v; alignment.tags) {
		}
	}
}
