import bamfile;
import validation.alignment;

import std.stdio;

class Validator : AbstractAlignmentValidator {
    void onError(ref Alignment alignment, AlignmentError e) {
        writeln(e, " in read '", alignment.read_name, "'");
    }

    void onError(ref Alignment alignment, CigarError e) {
        writeln("\tCIGAR error: ", e);
        writeln("\t\t   cigar string: ", alignment.cigar_string);
        writeln("\t\tsequence length: ", alignment.sequence_length);
    }

    void onError(string key, ref Value value, TagError e) {
        writeln("\tTag error: ", e , " ('", key, "' -> '", value.to_sam(), "')");
    }
}

void main(string[] args) {
    auto bf = BamFile(args[1]);
    auto validator = new Validator();
    foreach (alignment; bf.alignments) {
        validator.validate(alignment);
    }
}
