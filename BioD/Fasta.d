module BioD.Fasta;

import std.file;
import std.exception;
import std.algorithm;
import std.string;

struct FastaRecord {
    string header;
    string sequence;
}

auto fastaRecords(string filename) {

    static auto toFastaRecord(S)(S str) {
        auto res = findSplit(str, "\n");
        auto header = res[0];
        auto seq = res[2];
        return FastaRecord(header, removechars(seq, "\n"));
    }

    string text = cast(string)std.file.read(filename);

    enforce(text.length > 0 && text[0] == '>');
    text = text[1 .. $];

    auto records = splitter(text, '>');
    return map!toFastaRecord(records);
}
