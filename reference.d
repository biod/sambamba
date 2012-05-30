module reference;

import randomaccessmanager;

import std.stream;
import std.exception;

/**
  Stores reference sequence name and length
 */
struct ReferenceSequenceInfo {
    string name;
    int length;

    /**
      Constructs the structure from input stream
     */
    this(ref Stream stream) {
        int l_name; // length of the reference name plus one
        stream.read(l_name);
        name = stream.readString(l_name)[0..$-1].idup; // strip '\0' at the end
        stream.read(length);
    }
}

/**
  Represents reference sequence.
 */
struct ReferenceSequence {
   
    /// Name of reference sequence as in BAM file
    string name() @property {
        return _info.name;
    }

    /// Length in base pairs
    int length() @property {
        return _info.length;
    }

    /// Get alignments overlapping [start, end)
    auto opSlice(int start, int end) {
        enforce(start < end, "start must be less than end");
        return _manager.getAlignments(_ref_id, start, end);
    }

    this(RandomAccessManager manager, int ref_id, ReferenceSequenceInfo info) {
        _manager = manager;
        _ref_id = ref_id;
        _info = info;
    }

private:
    RandomAccessManager _manager;
    int _ref_id;
    ReferenceSequenceInfo _info;
}
