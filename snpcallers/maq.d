module snpcallers.maq;

/*
 * The code below is based on errmod.c from Samtools.
 */

import core.stdc.math;
import std.math : LN2, LN10;
import std.traits;
import std.range;
import std.algorithm;
import std.random;

import BioD.Base;
import BioD.Genotype;
import BioD.TinyMap;

struct BaseWithStrand {
    immutable ValueSetSize = Base.ValueSetSize * 2;
    private ubyte _code;
    ubyte internal_code() @property const {
        return _code;
    }

    static BaseWithStrand fromInternalCode(ubyte code) {
        BaseWithStrand bws;
        bws._code = code;
        return bws;
    }

    this(Base b, bool is_reverse) {
        _code = cast(ubyte)(b.internal_code * 2 + (is_reverse ? 1 : 0));
    }

    Base base() @property const {
        return Base.fromInternalCode(_code / 2);
    }

    bool is_reverse_strand() @property const {
        return (_code & 1) == 1;
    }
}

struct ReadBase {
    BaseWithStrand base_with_strand;
    alias base_with_strand this;
    private ubyte _quality;

    this(Base b, ubyte quality, bool is_reverse) {
        base_with_strand = BaseWithStrand(b, is_reverse);
        _quality = quality;
    }

    ubyte quality() @property const {
        return _quality;
    }
}

struct ErrorModelCoefficients {
    private {

        // _fk[n] = (1 - depcorr)^n * (1 - eta) + eta
        double[] _fk;

        // _beta[q << 16 | n << 8 | k ] = phred-scaled probability of ???
        double[] _beta;

        // _lhet[n << 8 | k] = log(1/2^n * choose(n, k))
        double[] _lhet;
        
        immutable Base[4] nucleotides = [Base('A'), Base('C'), Base('G'), Base('T')];
    }

    this(double depcorr, double eta) {
        _fk.length = 256;
        _beta.length = 256 * 256 * 64;
        _lhet.length = 256 * 256;

        foreach (n, ref v; _fk) {
            v = (1. - depcorr) ^^ n * (1. - eta) + eta;
        }

        // lC[n][k] = log(choose(n, k))
        double[256][256] lC;

        // lG[n] = logGamma(n + 1)
        double[256] lG;

        for (size_t n = 0; n <= 255; ++n) {
            lG[n] = core.stdc.math.lgamma(cast(double)(n + 1));
            for (size_t k = 0; k <= n / 2; ++k) {
                lC[n][n-k] = lC[n][k] = lG[n] - lG[k] - lG[n-k];

                // fill _lhet simultaneously
                _lhet[n << 8 | (n-k)] = _lhet[n << 8 | k] = lC[n][k] - n * LN2;
            }
        }

        for (size_t q = 1; q < 64; ++q) {
            real e = 10.0 ^^ (-cast(real)q / 10.0);
            real le = core.stdc.math.logl(e);
            real le1 = core.stdc.math.logl(1.0 - e);

            for (int n = 1; n <= 255; ++n) {
                real sum, sum1;
                sum = sum1 = 0;
                for (int k = n; k >= 0; --k) {
                    sum = sum1 + expl(lC[n][k] + k * le + (n-k) * le1);
                    _beta[q << 16 | n << 8 | k] = -10.0 / LN10 * logl(sum1 / sum);
                    sum1 = sum;
                }
            }
        }
    }

    double fk(size_t n) const {
        return _fk[n];
    }

    double beta(uint quality, size_t n, size_t k) const {
        return _beta[quality << 16 | n << 8 | k];
    }

    double lhet(size_t n, size_t k) const {
        return _lhet[n << 8 | k];
    }

    alias TinyMap!(DiploidGenotype, float, useDefaultValue) Dict;
    Dict computeLikelihoods(R)(R read_bases) const
        if (is(ElementType!R == ReadBase) && hasLength!R) 
    {
        // if there're more than 255 reads, subsample them
        ReadBase[255] buf;
        if (read_bases.length > buf.length) {
            copy(randomSample(read_bases, buf.length), buf[]);
        } else {
            copy(read_bases, buf[]);
        }
        auto bases = buf[0 .. min(read_bases.length, $)];

        sort!"a.quality < b.quality"(bases);

        auto w = TinyMap!(BaseWithStrand, uint, fillNoRemove)(0);
        auto c = TinyMap!(Base, uint, fillNoRemove)(0);
        auto fsum = TinyMap!(Base, double, fillNoRemove)(0.0);
        auto bsum = TinyMap!(Base, double, fillNoRemove)(0.0);

        foreach_reverse (ref read_base; bases) {
            auto quality = read_base.quality;
            if (quality < 4) quality = 4;
            if (quality > 63) quality = 63;
           
            auto bws = read_base.base_with_strand;
            auto b = bws.base;

            fsum[b] += fk(w[bws]);
            bsum[b] += fk(w[bws]) * beta(quality, bases.length, c[b]);
            c[b] += 1;
            w[bws] += 1;
        }

        alias DiploidGenotype G;

        auto q = Dict(float.min);

        foreach (i, b1; nucleotides) {
            float tmp1 = 0.0;
            int tmp2;
            float tmp3 = 0.0;

            // homozygous
            foreach (k, b2; nucleotides) {
                if (k != i) {
                    tmp1 += bsum[b2];
                    tmp2 += c[b2];
                    tmp3 += fsum[b2];
                }
            }

            int bar_e;
           
            if (tmp2 > 0) {
                q[G(b1, b1)] = tmp1;
            } else {
                q[G(b1, b1)] = 0.0;
            }

            // heterozygous
            for (size_t j = i + 1; j < nucleotides.length; ++j) {
                auto b2 = nucleotides[j];
                int cij = c[b1] + c[b2];
                tmp1 = tmp3 = 0.0;
                tmp2 = 0;
                foreach (k, b3; nucleotides) {
                    if (k != i && k != j) {
                        tmp1 += bsum[b3];
                        tmp2 += c[b3];
                        tmp3 += fsum[b3];
                    }
                }

                if (tmp2 > 0) {
                    q[G(b1, b2)] = q[G(b2, b1)] = tmp1 - 4.343 * lhet(cij, c[b2]);
                } else {
                    q[G(b1, b2)] = q[G(b2, b1)] = -4.343 * lhet(cij, c[b2]);
                }
            }

            foreach (k, b2; nucleotides) {
                auto g = G(b1, b2);
                if (g in q) {
                    if (q[g] < 0.0) q[g] = 0.0;
                }
            }
        }

        return q;
    }
}

class ErrorModel {
    
    private {
        float _depcorr;
        immutable _eta = 0.03;
        ErrorModelCoefficients _coef;
    }

    this(float depcorr) {
        _depcorr = depcorr;
        _coef = ErrorModelCoefficients(_depcorr, _eta);
    }

    const(ErrorModelCoefficients) coefficients() @property const {
        return _coef;
    }

    alias coefficients this;
}

void main(string[] args) {
    import bamfile;
    import reconstruct;
    import pileuprange;

    import std.stdio;
    import std.algorithm;
    import std.array;

    auto fn = args[1];
    auto reads = filter!"!a.is_unmapped"(BamFile(fn).alignments);

    static ReadBase toReadBase(R)(R read) {
        return ReadBase(Base(read.current_base),
                             min(read.current_base_quality, read.mapping_quality),
                             read.is_reverse_strand);
    }

    auto errmod = new ErrorModel(0.17);
    auto pileup = pileupWithReferenceBases(reads);
   
    auto pos = 0;
    ReadBase[8192] bases;
    foreach (column; pileup) {
        if (column.position > pos) {
            writeln(column.position);
            pos += 10000;
        }
        auto valid = filter!"a.mapping_quality != 255 && a.current_base != '-'"(column.reads);
        auto rbs = map!toReadBase(valid);

        size_t i = 0;
        while (i < 8192) {
            if (rbs.empty) break;
            bases[i++] = rbs.front;
            rbs.popFront();
        }

        auto read_bases = bases[0 .. i];
        if (read_bases.length == 0) continue;

        auto genotype_likelihoods = errmod.computeLikelihoods(read_bases);
        /*
        auto ref_base = Base(column.reference_base);
        auto g = DiploidGenotype(ref_base, ref_base);
        if (genotype_likelihoods[g] == 0.0)
            continue;

        writeln("--------------------------------------------------------------------");
        writeln("Reference base: ", column.reference_base);
        writeln("1-based position: ", column.position + 1);
        writeln("Read bases: ", map!"a.current_base"(valid));

        foreach (g, likelihood; genotype_likelihoods) {
            if (g.base1 == 'N' || g.base2 == 'N') continue;
            writeln("<", g.base1, "|", g.base2, "> ", likelihood);
        }*/
    }
}
