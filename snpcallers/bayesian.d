module snpcallers.bayesian;

import pileuprange;

import std.algorithm;
import std.typecons;
import std.range;
import std.math;
import std.conv;

private {

    immutable float[] table;

    static this() {
        table = cast(immutable)array(take(recurrence!"a[n-1] * 0.85f"(1.0f), 512));
    }

    float f(size_t k) {
        return k >= table.length ? 0.0f : table[k];
    }

    // Phred-base quality -> log10 of probability that the base is wrong
    float normalize(int base_error) {
        if (base_error == 255) {
            return float.nan;
        }
        float minus_log10_base_error = cast(float)base_error / 10.0;
        return -minus_log10_base_error;
    }

    // Computes probability that bases corresponding to error_rates_wrong
    // are indeed wrong, and bases that correspond to error_rates_correct
    // are correct.
    //
    // error_rates_* must be normalized
    float logErrorProb(R1, R2)(R1 error_rates_wrong, R2 error_rates_correct) {

        auto log10_site_error = 0.0f;
        auto site_error_denominator = 0.0f;

        auto fix_w = array(filter!"!isnan(a)"(error_rates_wrong));
        auto fix_c = filter!"!isnan(a)"(error_rates_correct);

        fix_w.sort;

        size_t k;

        auto log_error_prob = 0.0f;
        foreach (e; fix_w) {
            log_error_prob += f(k) * e;
            site_error_denominator += f(k);
            k += 1;
        }
        
        // now k is equal to length of fix_w

        if (!fix_c.empty) {
            log10_site_error = log_error_prob / site_error_denominator;

            auto n = k + walkLength(fix_c);
            log_error_prob += coeff(n, k, log10_site_error);
        }

        return log_error_prob;
    }

    float coeff(size_t n, size_t k, float log10_base_error) {
        return 1.0f; // TODO: use precomputed table
    }
}

struct BayesianSnpInfo {
    Genotype estimated_genotype;
    float quality;
}

Nullable!BayesianSnpInfo snpInfo(C)(ref C column, ref BayesianCallerSettings settings) {

    Nullable!BayesianSnpInfo result;

    if (column.coverage < settings.minimum_coverage) {
        return result;
    }

    uint[char] base_counts;

    foreach (read; column.reads) {
        auto base = read.current_base;
        if (base != '-') {
            base_counts[base] += 1;
        }
    }

    if (base_counts.length == 0) {
        return result;
    }

    if (base_counts.length == 1) {
        // only one possible genotype

        char allele = base_counts.byKey().front;
        if (allele == column.reference_base) {
            return result;
        }

        auto error_rates = map!(r => normalize(r.current_base_quality))(column.reads);

        auto log10_site_error = logErrorProb(error_rates, cast(float[])[]);

        if (-log10_site_error < settings.minimum_quality) {
            return result;
        }

        BayesianSnpInfo info;
        info.quality = -log10_site_error;
        info.estimated_genotype = Genotype(Base(allele));

        result = info;
        return result;
    }
    /* TODO
    if (base_counts.length >= 2) {
        // select two more frequent bases
        Base[] alleles = base_counts.values;
        partialSort!((a, b) => base_counts[a] > base_counts[b])(alleles, 2);

        Genotype[3] genotypes;
        float[3] log_probs;
        genotypes[0] = Genotype(alleles[0]);
        genotypes[1] = Genotype(alleles[0], alleles[1]);
        genotypes[2] = Genotype(alleles[1]);

    }*/

    return result;
}

/// Settings for bayesian caller
struct BayesianCallerSettings {
    int minimum_coverage = 1; /// minimum number of reads at the site
    float minimum_quality = 5.0; /// minimum quality of genotype estimation
}

/// Default calling settings
BayesianCallerSettings defaultSettings;

/// Represents a base
struct Base {
    char base;
    alias base this;
}

/// Represents genotype of a diploid individual
struct Genotype {

    this(Base base) {
        bases[0] = bases[1] = base;
    }

    this(Base base1, Base base2) {
        bases[0] = base1;
        bases[1] = base2;
    }

    /// Alleles (equal to each other in case of homozygosity)
    Base[2] bases;
}

/// Column augmented with estimated genotype and its quality
struct VariantColumn(C) {
    C column;
    alias column this;

    Genotype estimated_genotype;
    float quality;
}

/// Returns range of SNPs given range of reads with 'MD' tags.
auto findSNPs(R)(R reads, ref BayesianCallerSettings settings=defaultSettings) {
    static struct Range(C) {
        this(C columns, ref BayesianCallerSettings settings) {
            _columns = columns;
            _settings = settings;
            _findNextSNP();
        }

        VariantColumn!(ElementType!C) front() @property {
            return _front;
        }
      
        bool empty() @property {
            return _empty;
        }

        void popFront() @property {
            _columns.popFront();
            _findNextSNP();
        }

        private {
            C _columns;
            BayesianCallerSettings _settings;

            bool _empty;
            VariantColumn!(ElementType!C) _front;

            void _findNextSNP() {
                while (true) {
                    if (_columns.empty) {
                        _empty = true;
                        return;
                    }

                    auto column = _columns.front;

                    auto info = snpInfo(column, _settings);

                    if (info.isNull) {
                        _columns.popFront();
                    } else {
                        _front.column = column;
                        _front.estimated_genotype = info.estimated_genotype;
                        _front.quality = info.quality;
                        return;
                    }
                }
            }
        }
    }

    static auto snpRange(C)(C columns, ref BayesianCallerSettings settings) {
        return Range!C(columns, settings);
    }

    return snpRange(pileupWithReferenceBases(reads), settings);
}
