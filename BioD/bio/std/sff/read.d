module bio.std.sff.read;

/// SFF record
struct SffRead {
    /// Read identifier
    string name;

    /// Homopolymer stretch estimates for each flow of the read
    ushort[] flowgram_values;

    ubyte[] flow_index_per_base;

    /// Basecalled nucleotide sequence
    char[] bases;

    /// Phred-scaled quality scores
    ubyte[] quality_scores;

    ushort clip_qual_left;
    ushort clip_qual_right;
    ushort clip_adapter_left;
    ushort clip_adapter_right;

    /// Record start offset in the file
    ulong file_offset;
}
