import bamfile;
import validation.alignment;

import sam.serialize;
import std.stdio;

class Validator : AbstractAlignmentValidator {
    bool onError(ref Alignment alignment, AlignmentError e) {
        writeln(e, " in read '", alignment.read_name, "'");
        return true; // continue checks
    }

    bool onError(ref Alignment alignment, CigarError e) {
        writeln("\tCIGAR error: ", e);
        writeln("\t\t   cigar string: ", alignment.cigarString());
        writeln("\t\tsequence length: ", alignment.sequence_length);
        return true;
    }

    bool onError(string key, ref Value value, TagError e) {
        writeln("\tTag error: ", e , " ('", key, "' -> '", toSam(value), "')");
        return true;
    }
}

void main(string[] args) {
    auto bf = BamFile(args[1]);
    auto validator = new Validator();
    foreach (alignment; bf.alignments) {
        validator.validate(alignment);
    }
}
