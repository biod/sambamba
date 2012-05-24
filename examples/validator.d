import bamfile;
import validation.alignment;

import std.stdio;

class CigarStringValidator : AbstractAlignmentValidator {
    void onError(ref Alignment alignment, AlignmentError e) {
    }

    void onError(ref Alignment alignment, CigarError e) {
        writeln("Validation error: ", e);
        writeln("\t   cigar string: ", alignment.cigar_string);
        writeln("\tsequence length: ", alignment.sequence_length);
    }

    void onError(string key, ref Value value, TagError e) {}
}

void main(string[] args) {
    auto bf = BamFile(args[1]);
    auto validator = new CigarStringValidator();
    foreach (alignment; bf.alignments) {
        validator.validate(alignment);
    }
}
