module bio.std.hts.snpcallers.maq;

/*
 * The code below is based on errmod.c from Samtools.
 */

import core.stdc.math;
import std.math : LN2, LN10, isNaN;
import std.traits;
import std.range;
import std.algorithm;
import std.random;
import std.typecons;

import bio.std.hts.bam.md.reconstruct;
import bio.std.hts.bam.pileup;

import bio.core.base;
import bio.core.genotype;
import bio.core.call;
import bio.core.tinymap;

struct BaseWithStrand {
    immutable ValueSetSize = Base.ValueSetSize * 2;
    private ubyte _code;
    ubyte internal_code() @property const {
        return _code;
    }

    static BaseWithStrand fromInternalCode(ubyte code) {
        BaseWithStrand bws = void;
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

        // _beta[q << 16 | n << 8 | k ] = see MAQ paper for meaning of \beta
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
            v = core.stdc.math.pow(1.0 - depcorr, cast(double)n) * (1.0 - eta) + eta;
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
                _lhet[n << 8 | (n-k)] = _lhet[n << 8 | k] = lC[n][k] - n * cast(double)LN2;
            }
        }

        for (size_t q = 1; q < 64; ++q) {
            real e = 10.0 ^^ (-(cast(real)q) / 10.0);
            real le = core.stdc.math.logl(e);
            real le1 = core.stdc.math.logl(1.0 - e);

            for (int n = 1; n <= 255; ++n) {
                real sum, sum1;
                sum = sum1 = 0.0;
                for (int k = n; k >= 0; --k) {
                    sum = sum1 + core.stdc.math.expl(lC[n][k] + k * le + (n-k) * le1);
                    _beta[q << 16 | n << 8 | k] = -10.0 / LN10 * core.stdc.math.logl(sum1 / sum);
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

    alias TinyMap!(DiploidGenotype!Base5, float, useDefaultValue) Dict;

    private immutable C = 10.0 / LN10;

    Dict computeLikelihoods(R)(R read_bases, bool symmetric=false) const
        if (is(ElementType!R == ReadBase) && hasLength!R) 
    {
        // if there're more than 255 reads, subsample them
        ReadBase[255] buf = void;
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

        alias diploidGenotype dG;

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

            auto b1_5 = cast(Base5)b1;
            if (tmp2 > 0) {
                q[dG(b1_5)] = tmp1;
            } else {
                q[dG(b1_5)] = 0.0;
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

                auto b2_5 = cast(Base5)b2;
                if (tmp2 > 0) {
                    q[dG(b2_5, b1_5)] = tmp1 - C * lhet(cij, c[b2]);
                } else {
                    q[dG(b2_5, b1_5)] = -C * lhet(cij, c[b2]);
                }

                if (symmetric) {
                    q[dG(b1_5, b2_5)] = q[dG(b2_5, b1_5)];
                }
            }

            foreach (k, b2; nucleotides) {
                auto g = dG(b1_5, cast(Base5)b2);
                if (g in q) {
                    if (q[g] < 0.0) q[g] = 0.0;
                }
            }
        }

        return q;
    }
}

// Encapsulates information about genotype likelihoods at a site.
struct GenotypeLikelihoodInfo {

    alias ErrorModelCoefficients.Dict ScoreDict;

    alias DiploidGenotype!Base5 Gt;

    this(ScoreDict dict) {

        _dict = dict;
        size_t k = 0;

        // copy all data into a buffer, combining that with insertion sort
        foreach (gt, score; _dict) {
            if (k == 0) {
                gt_buf[k++] = gt;
            } else {
                size_t j = k;
                while (j > 0 && _dict[gt_buf[j-1]] > score) {
                    gt_buf[j] = gt_buf[j-1];
                    --j;
                }
                gt_buf[j] = gt;
                ++k;
            }
        }

        assert(k >= 2);

        _count = cast(ubyte)k;
    }

    size_t count() @property const {
        return _count;
    }

    static struct GtInfo {
        private {
            Gt _gt;
            float _prob;
        } 

        Gt genotype() @property const {
            return _gt;
        }

        float score() @property const {
            return _prob;
        }
    }

    GtInfo opIndex(size_t index) {
        assert(index < count);
        auto gt = gt_buf[index];
        return GtInfo(gt, _dict[gt]);
    }

    private Gt[25] gt_buf;
    private ubyte _count;
    private ScoreDict _dict;
}

class ErrorModel {
    
    private {
        float _depcorr;
        float _eta;
        ErrorModelCoefficients _coef;
    }

    this(float depcorr, float eta=0.03) {
        _depcorr = depcorr;
        _eta = eta;
        _coef = ErrorModelCoefficients(_depcorr, _eta);
    }

    const(ErrorModelCoefficients) coefficients() @property const {
        return _coef;
    }

    alias coefficients this;
}

/// Class for calling SNPs using MAQ model.
///
/// Typical usage:
///     auto caller = new MaqSnpCaller();
///     caller.minimum_call_quality = 20.0f;
///     caller.minimum_base_quality = 13;
///     foreach (snp; caller.findSNPs(reads)) { ... }
///
final class MaqSnpCaller {
    
    private float _depcorr = 0.17;
    private float _eta = 0.03;
    private float _minimum_call_quality = 6.0;
    private ubyte _minimum_base_quality = 13;
    private bool _need_to_recompute_errmod = true;

    ///
    float depcorr() @property const {
        return _depcorr;
    }

    /// ditto
    void depcorr(float f) @property {
        _depcorr = f;
        _need_to_recompute_errmod = true;
    }

    ///
    float eta() @property const {
        return _eta;
    }

    ///
    void eta(float f) @property {
        _eta = f;
        _need_to_recompute_errmod = true;
    }
    
    /// Minimum call quality
    float minimum_call_quality() @property const {
        return _minimum_call_quality;
    }

    /// ditto
    void minimum_call_quality(float f) @property {
        _minimum_call_quality = f;
    }

    /// Discard reads with base quality less than this at a site
    ubyte minimum_base_quality() @property const {
        return _minimum_base_quality;
    }

    void minimum_base_quality(ubyte q) @property {
        _minimum_base_quality = q;
    }

    ErrorModel errmod() @property {
        if (_need_to_recompute_errmod) {
            synchronized {
                if (_need_to_recompute_errmod) {
                    _errmod = new ErrorModel(_depcorr, _eta);
                    _need_to_recompute_errmod = false;
                }
            }
        }
        return _errmod;
    }

    private ErrorModel _errmod;

    /// Get genotype likelihoods
    final GenotypeLikelihoodInfo genotypeLikelihoodInfo(C)(C column) {

        version(MaqCaller8192) {
            ReadBase[8192] buf = void;
        }

        size_t num_of_valid_bases = 0;

        foreach (read; column.reads) {

            version(MaqCaller8192) {
                if (num_of_valid_bases == 8192) break;
            }

            if (read.current_base_quality < minimum_base_quality)
                continue;
            if (read.current_base == '-')
                continue;

            version(MaqCaller8192) {
                buf[num_of_valid_bases] = ReadBase(Base(read.current_base),
                                                   min(read.current_base_quality, read.mapping_quality),
                                                   read.is_reverse_strand);
            }

            num_of_valid_bases++;
        }

        static struct ReadBaseRange(R) {
            private R _reads = void;
            private ubyte minimum_base_quality = void;
                                              
            this(R reads, ubyte minbq) { 
                _reads = reads; minimum_base_quality = minbq; _findNextValid();
            }

            ReadBase front() @property { 
                auto read = _reads.front;
                return ReadBase(Base(read.current_base), 
                                min(read.current_base_quality, read.mapping_quality),
                                read.is_reverse_strand);
            }
            bool empty() @property { return _reads.empty; }
            void popFront() { _reads.popFront(); _findNextValid(); }
            ReadBaseRange save() @property { return ReadBaseRange!R(_reads, minimum_base_quality); }

            private void _findNextValid() {
                while (!_reads.empty && 
                        (_reads.front.current_base_quality < minimum_base_quality ||
                         _reads.front.current_base == '-')) 
                {
                    _reads.popFront();
                }
            }
        }

        if (num_of_valid_bases == 0) {
            GenotypeLikelihoodInfo result;
            return result;
        }

        version(MaqCaller8192) {
            ReadBase[] rbs = buf[0 .. num_of_valid_bases];
            auto likelihood_dict = errmod.computeLikelihoods(rbs);
        } else {
            auto rbs = ReadBaseRange!(typeof(column.reads))(column.reads, minimum_base_quality);
            auto likelihood_dict = errmod.computeLikelihoods(takeExactly(rbs, num_of_valid_bases));
        }
        return GenotypeLikelihoodInfo(likelihood_dict);
    }

    /// Make call on a pileup column
    final Nullable!DiploidCall5 makeCall(C)(C column, string reference="", string sample="") {

        auto gts = genotypeLikelihoodInfo(column);

        Nullable!DiploidCall5 result;

        if (gts.count < 2) return result;

        static if (__traits(compiles, column.reference_base)) {
            auto refbase = Base5(column.reference_base);
        } else {
            auto refbase = Base5('N');
        }
        
        if (sample == "") {
            auto rg = column.reads.front["RG"];
            if (!rg.is_nothing) {
                sample = cast(string)rg;
            }
        }

        result = DiploidCall5(sample, reference, column.position,
                              refbase, gts[0].genotype,
                              gts[1].score - gts[0].score);
                
        return result;
    }

    /// main method of this class
    auto findSNPs(P)(P pileup_columns, string reference="", string sample="") {
        static assert(__traits(compiles, {pileup_columns.front.reference_base;}));

        static struct Result {
            private MaqSnpCaller _caller;
            private P _pileup;
            private DiploidCall5 _front;
            private bool _empty;
            private string _reference;
            private string _sample;

            this(MaqSnpCaller caller, P pileup, string reference, string sample) {
                _caller = caller;
                _pileup = pileup;
                _reference = reference;
                _sample = sample;
                _fetchNextSNP();
            }

            DiploidCall5 front() @property {
                return _front;
            }
           
            bool empty() @property {
                return _empty;
            }

            void popFront() {
                _pileup.popFront();
                _fetchNextSNP();
            }

            private void _fetchNextSNP() {
                while (true) {
                    if (_pileup.empty) {
                        _empty = true;
                        break;
                    }

                    auto call = _caller.makeCall(_pileup.front, _reference, _sample);
                    if (!call.isNull && call.is_variant && call.quality > _caller.minimum_call_quality) {
                        _front = call.get;
                        break;
                    } else {
                        _pileup.popFront();
                    }
                }
            }
        }

        return Result(this, pileup_columns, reference, sample);
    }
}
