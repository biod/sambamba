module bio.core.kmer;

import bio.core.base;
import std.range;

/// Represents k-mer of ACGT bases of length no more than 32.
struct KMer(uint K) 
    if (K <= 32)
{
    private ulong _id;

    static Base5 code2base(int code) {
        return Base5("ACGT"[code]);
    }

    static int char2code(char base) {
        switch (base) {
            case 'A': return 0;
            case 'C': return 1;
            case 'G': return 2;
            case 'T': return 3;
            default: return -1;
        }
    }

    /// Unique ID
    ulong id() @property const {
        return _id;
    }

    /// Construct by ID
    this(S)(S id) 
        if (is(S == ulong))
    {
        _id = id;
    }
       
    /// Construct from sequence. Takes bases from the provided sequence
    /// until K symbols 'A/C/G/T' are found. That is, 'N' and other ambiguous
    /// bases are skipped.
    ///
    /// If sequence does not contain at least K bases 'A/C/G/T', the result of
    /// operation is undefined.
    this(S)(S sequence) 
        if (isInputRange!S) 
    {
        size_t i = 0;
        foreach (nuc; sequence) {
            _id <<= 2;
            ++i;
            switch (cast(char)nuc) {
                case 'A':
                    break;
                case 'C':
                    _id += 1;
                    break;
                case 'G':
                    _id += 2;
                    break;
                case 'T':
                    _id += 3;
                    break;
                default:
                    _id >>= 2;
                    --i;
                    break;
            }

            if (i == K)
                break;
        }
    }

    struct KMerSequence {
        this(ulong number) {
            _n = number;
        }

        private ulong _n;
        private size_t _len = K;

        bool empty() @property const { return _len == 0; }
        void popFront() { --_len; }
        void popBack() { --_len; _n >>= 2; }

        Base5 opIndex(size_t i) const {
            return code2base((_n >> (2 * (_len - i - 1))) & 3);
        }

        size_t length() @property const { return _len; }
        Base5 front() @property const { return opIndex(0); }
        Base5 back() @property const { return opIndex(_len - 1); }
        KMerSequence save() @property const { 
            KMerSequence _seq = void;
            _seq._n = _n;
            _seq._len = _len;
            return _seq;
        }
    }

    /// Sequence corresponding to the k-mer
    KMerSequence sequence() @property const {
        return KMerSequence(_id);
    }
}

unittest {
    import std.algorithm;
    auto kmer = KMer!10("AACGTACGTG");
    assert(equal(kmer.sequence, "AACGTACGTG"));

    assert(KMer!5(KMer!5(0b1011001001UL).sequence).id == 0b1011001001UL);
}
