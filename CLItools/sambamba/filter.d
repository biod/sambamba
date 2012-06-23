/**
  Set of filters for alignments. 
  All share a common interface and can be easily combined.
*/
module filter;

import alignment;
import tagvalue;
import validation.alignment;

/// Common interface for all filters
interface Filter {
    bool accepts(ref Alignment a) const;
}

/// Filter which accepts all alignments
final class NullFilter : Filter {
    bool accepts(ref Alignment a) const {
        return true;
    }
}

/// Filter which accepts only alignments with big enough mapping quality
final class MappingQualityFilter : Filter {

    private ubyte _min_q;

    this(ubyte minimal_quality) {
        _min_q = minimal_quality;
    }

    bool accepts(ref Alignment a) const {
        return a.mapping_quality >= _min_q;
    }
}

/// Filter which accepts only alignments with a given read group
final class ReadGroupFilter : Filter {

    private string _rg;

    this(string read_group) {
        _rg = read_group;
    }

    bool accepts(ref Alignment a) const {
        auto rg = a["RG"];
        if (!rg.is_string || (to!string(rg) != _rg)) {
            return false;
        }
        return true;
    }
}

final class ValidAlignmentFilter : Filter {
    
    bool accepts(ref Alignment a) const {
        return isValid(a);
    }
}

/// Filter which accepts only alignments accepted by both filters
final class AndFilter : Filter {

    private Filter _a, _b;

    this(Filter a, Filter b) {
        _a = a; _b = b;
    }

    bool accepts(ref Alignment a) const {
        return _a.accepts(a) && _b.accepts(a);
    }
}
