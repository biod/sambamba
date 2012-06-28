module utils.tagstoragebuilder;

import tagvalue;
import utils.value;

import std.array;

/// Struct for building tag storage, effectively managing memory.
struct TagStorageBuilder {
    private Appender!(ubyte[]) _storage;

    /// Return tag data (little-endian, in BAM format)
    ubyte[] data() @property {
        return _storage.data();
    }

    static TagStorageBuilder create() {
        TagStorageBuilder builder;
        builder._storage = appender!(ubyte[])();
        return builder;
    }

    /// Clear underlying storage
    void clear() {
        _storage.clear();
    }

    /// Append tag value to the storage
    void put(string tag, ref Value value) {
        _storage.put(cast(ubyte[])tag);
        emplaceValue(_storage, value);
    }
}
