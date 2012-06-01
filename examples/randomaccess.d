import bamfile;
import std.conv;
import std.stdio;

void main(string[] args) {
    BamFile bam;
    string chr;
    int beg;
    int end;
    try {
        bam = BamFile(args[1]);
        chr = args[2];
        beg = to!int(args[3]);
        end = to!int(args[4]);
    } catch (Throwable e) {
        writeln("usage: " ~ args[0] ~ " <input.bam> <chromosome> <begin> <end>");
        return;
    }

    foreach (alignment; bam[chr][beg .. end]) {
        writeln(alignment.read_name, " ", 
                to!string(alignment.sequence), " ",
                alignment.cigar_string(),
                " (begin: ", alignment.position,
                "; end: ", alignment.position + alignment.bases_covered(),
                ")");
    }
}
