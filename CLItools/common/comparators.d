/* Contains comparison functions for alignment records */
module common.comparators;

import alignment;

bool compareReadNames(Alignment a1, Alignment a2) {
    return a1.read_name < a2.read_name;
}

/// Comparison function for alignments
bool compareAlignmentCoordinates(Alignment a1, Alignment a2) {
    if (a1.ref_id == -1) return false; // unmapped reads should be last
    if (a2.ref_id == -1) return true;
    if (a1.ref_id < a2.ref_id) return true;
    if (a1.ref_id > a2.ref_id) return false;
    if (a1.position < a2.position) return true;
    return false;
}
