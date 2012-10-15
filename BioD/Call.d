module BioD.Call;

import BioD.Base;
import BioD.Genotype;

/// A genotype call
struct Call(alias Gt, B) 
{
    alias Gt!B G;

    private {
        string _sample = void;
        string _chr = void;
        ulong _pos = void;
        B _refbase = void;
        G _gt = void;
        float _qual = void;
    }

    /// Constructor
    this(string sample, string chr, ulong pos,
         B refbase, G genotype, float quality=float.nan)
    {
        _sample = sample;
        _chr = chr;
        _pos = pos;
        _refbase = refbase;
        _gt = genotype;
        _qual = quality;
    }

    /// Sample name
    string sample() @property const {
        return _sample;
    }

    /// Chromosome name
    string chromosome() @property const {
        return _chr;
    }

    /// 0-based position on the reference
    ulong position() @property const {
        return _pos;
    }

    /// Reference base at the site
    B reference_base() @property const {
        return _refbase;
    }

    /// Most probable genotype
    ref const(G) genotype() @property const {
        return _gt;
    }

    ///
    alias genotype this; 

    /// Phred-scaled quality score. If unknown, set to NaN.
    float quality() @property const {
        return _qual;
    }

    /// Returns true if this call is not a reference one.
    bool is_variant() @property const {
        return _gt != G(_refbase);
    }
}

alias Call!(DiploidGenotype, Base5) DiploidCall5;
alias Call!(DiploidGenotype, Base16) DiploidCall;
alias DiploidCall DiploidCall16;

unittest {
    auto call = DiploidCall("NA01234", "chr10", 543210,
                            Base('T'), diploidGenotype(Base('C'), Base('T')),
                            47.0);

    assert(call.is_variant);
    assert(call.is_heterozygous);
    assert(call.reference_base == 'T');
}
