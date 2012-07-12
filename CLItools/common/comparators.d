/* Contains comparison functions for alignment records */
module common.comparators;

import alignment;

/// Comparison function for alignments
bool compareAlignmentCoordinates(Alignment a1, Alignment a2) {
    if (a1.ref_id < a2.ref_id) return true;
    if (a1.ref_id > a2.ref_id) return false;
    if (a1.position < a2.position) return true;
    return false;
}
