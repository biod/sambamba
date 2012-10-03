module snpcallers.simple;

import pileuprange;

import std.algorithm;
import utils.algo;

struct SimpleCallerSettings {
    int minimum_coverage = 5;
    int minimum_witnesses = 2;
}

SimpleCallerSettings defaultSettings;

bool isSNP(C)(C column, ref SimpleCallerSettings settings)
{
    if (column.coverage < settings.minimum_coverage) {
        return false;
    }

    int[char] bases_count;

    foreach (read; column.reads) {
        if (read.current_base != '-') {
            bases_count[read.current_base] += 1;
        }
    }

    if (bases_count.length == 0) {
        // e.g. all overlapping reads have deletions at this location
        return false;
    }

    auto consensus = argmax!(base => bases_count[base])(bases_count.byKey());

    if (bases_count[consensus] < settings.minimum_witnesses) {
        return false;
    }

    return consensus != column.reference_base;
}

auto findSNPs(R)(R reads, ref SimpleCallerSettings settings=defaultSettings) {
    auto columns = pileupWithReferenceBases(reads);

    return filter!(column => isSNP(column, settings))(filter!"a.coverage > 0"(columns));
}
