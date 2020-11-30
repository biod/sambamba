module bio.std.experimental.hts.constants;

alias ulong size_d;
alias int RefId;      // -1 is invalid FIXME
alias uint GenomePos; // 32-bits, do check when reading!
alias ubyte MappingQuality; // 255 is invalid FIXME

struct GenomeLocation {
  RefId ref_id;
  GenomePos pos;
}
