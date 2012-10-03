module snpcallers.bayesian;

import pileuprange;
import std.algorithm;
import std.range;
import std.math;
import std.conv;

struct BayesianCallerSettings {
    int minimum_coverage = 1;
    float minimum_quality = 5.0;
}

private BayesianCallerSettings defaultSettings;

struct Base {
    char base;
    alias base this;
}

struct Genotype {

    this(Base base) {
        bases[0] = bases[1] = base;
    }

    this(Base base1, Base base2) {
        bases[0] = base1;
        bases[1] = base2;
    }

    Base[2] bases;
}

struct VariantColumn(C) {
    C column;
    alias column this;

    Genotype estimated_genotype;
    
    float quality;
}

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

            static immutable float[] table;

            static float[] generateTable() @trusted {
                float[] array;
                array ~= 1.0f;
                foreach (i; 1 .. 128) {
                    array ~= array[i - 1] * 0.85f;
                }
                return array;
            }

            static this() {
                table = cast(immutable)generateTable();
            }

            float f(uint k) {
                return k >= 128 ? 0.0f : table[k];
            }

            // Phred-base quality -> log10 of probability that the base is wrong
            float normalize(int base_error) {
                float minus_log10_base_error = cast(float)base_error / 10.0;
                return -minus_log10_base_error;
            }

            void _findNextSNP() {
                while (true) {
                    if (_columns.empty) {
                        _empty = true;
                        return;
                    }

                    auto column = _columns.front;

                    if (column.coverage < _settings.minimum_coverage) {
                        goto pop_front;
                    }

                    uint[Base] base_counts;

                    foreach (read; column.reads) {
                        auto base = Base(read.current_base);
                        if (base != '-') {
                            base_counts[base] += 1;
                        }
                    }


                    if (base_counts.length == 0) {
                        goto pop_front;
                    }

                    if (base_counts.length == 1) {
                        // only one possible genotype

                        char allele = base_counts.byKey().front;
                        if (allele == column.reference_base) {
                            goto pop_front;
                        }

                        // probability that all bases are wrong
                        float log10_site_error = 0.0f;

                        auto k = 0;
                        float denominator = 0.0f;
                        foreach (read; column.reads) {
                            if (read.current_base_quality == 255) 
                                continue;
                            auto log10_base_error = normalize(read.current_base_quality);
                            log10_site_error += f(k) * log10_base_error;
                            denominator += f(k);
                            k += 1;
                        }

                        log10_site_error /= denominator;
                        log10_site_error *= base_counts[Base(allele)];

                        _front.quality = -log10_site_error;

                        if (_front.quality < _settings.minimum_quality) {
                            goto pop_front;
                        }

                        _front.column = column;
                        _front.estimated_genotype = Genotype(Base(allele));
                        goto stay_at_current_column;
                    }

pop_front:
                    _columns.popFront();
                    continue;
stay_at_current_column:
                    break;
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
                }
            }
        }
    }

    static auto snpRange(C)(C columns, ref BayesianCallerSettings settings) {
        return Range!C(columns, settings);
    }

    return snpRange(pileupWithReferenceBases(reads), settings);
}
