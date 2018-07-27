import bio.bam.reader;
import bio.bam.writer;
import std.random;

void main(string[] args) {
  auto bam = new BamReader(args[1]);
  auto w = new BamWriter(args[2]);
  scope(exit) w.finish();

  auto gen = Random(unpredictableSeed);

  w.writeSamHeader(bam.header);
  w.writeReferenceSequenceInfo(bam.reference_sequences);

  foreach (r; bam.reads) {
    auto new_seq = new ubyte[r.sequence.length];
    foreach (ref x; new_seq)
      x = "ACGT"[uniform(0, 4, gen)];
    r.sequence = cast(string)new_seq;
    w.writeRecord(r);
  }
}
