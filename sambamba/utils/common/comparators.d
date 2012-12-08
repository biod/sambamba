/*
    This file is part of Sambamba.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

*/
/* Contains comparison functions for alignment records */
module sambamba.utils.common.comparators;

import bio.bam.read;

bool compareReadNames(BamRead a1, BamRead a2) {
    return a1.name < a2.name;
}

/// Comparison function for alignments
bool compareAlignmentCoordinates(BamRead a1, BamRead a2) {
    if (a1.ref_id == -1) return false; // unmapped reads should be last
    if (a2.ref_id == -1) return true;
    if (a1.ref_id < a2.ref_id) return true;
    if (a1.ref_id > a2.ref_id) return false;
    if (a1.position < a2.position) return true;
    return false;
}
