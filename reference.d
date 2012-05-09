module reference;

import std.stream;

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
