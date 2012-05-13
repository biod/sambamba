/**
   Contains classes implementing TagStorage interface for accessing tags. 
   The classes encapsulate tag parsing logic as well.

   The approach to parsing used by original samtools is as follows:

   1) don't parse tags eagerly
   2) store ubyte[] array in the struct
   3) change endianness in-place, if needed
        - thus we can slice strings and arrays -- no extra allocations!
   4) to retrieve value, iterate over all tags
        - micro-optimization: compare ushorts instead of 2-char strings
*/
module tagstorage;

import tagvalue;

import std.stream;
import std.conv;
import std.system;

/**
  Abstract tag storage. 

  Provides hash-like access (currently read-only) and opportunity to iterate
  storage like an associative array.
*/
abstract class TagStorage {
    static TagStorage createStorage(ubyte[] chunk); /// initializes tag storage from chunk
    abstract Value opIndex(string s); /// provides hash-like access

    abstract int opApply(int delegate(string, Value value) dg); /// iteration
}

/**
  Thrown in case of unrecognized tag type
 */
class UnknownTagTypeException : Exception {
    this(string msg) { super(msg); }
}

/**
  Eagerly parses all tags
*/
class EagerTagStorage : TagStorage {
    /**
      Takes chunk of memory which contains alignment auxiliary data.
     */
    this(ubyte[] chunk) {

        scope Stream memory_stream = new MemoryStream(chunk);
        scope Stream stream = new EndianStream(memory_stream,
                                  Endian.littleEndian);

        while (!stream.eof()) {
            char[] tag = stream.readString(2);
            char type = void;
            stream.read(type);

            if (type == 'Z' || type == 'H') {
                _tags[tag.idup] = readString(stream, type);
            } else if (type != 'B') {
                _tags[tag.idup] = readPrimitive(stream, type);
            } else {
                char elem_type = void;
                stream.read(elem_type);
                _tags[tag.idup] = readArray(stream, elem_type);
            }
        } 
        
        _tags.rehash(); // makes lookups faster
    }

    static TagStorage createStorage(ubyte[] chunk) {
        return new EagerTagStorage(chunk);
    }

    int opApply(int delegate(string, Value value) dg) {
        foreach (s, v; _tags) {
            auto result = dg(s, v);
            if (result != 0) {
                return result;
            }
        }
        return 0;
    }

    Value opIndex(string s) {
        return _tags[s];
    }

private:

    Value[string] _tags;

    Value readString(ref Stream stream, char type) {
        char[] s;
        char c;
        do {
            stream.read(c);
            s ~= c;
        } while (c != 0);
        return Value(cast(string)s[0 .. $-1]);
    }

    Value readPrimitive(ref Stream stream, char type) {
        string readPrimitiveHelper() {
            char[] cases;
            foreach (c2t; PrimitiveTagValueTypes) {
                cases ~= "case '".dup~c2t.ch~"':"~
                         "    "~c2t.ValueType.stringof~" val;"~
                         "    stream.read(val);"~
                         "    return Value(val);".dup;
            }
            return to!string("switch (type) {"~ 
                        cases~
                   "    default: throw new UnknownTagTypeException(to!string(type));"~
                   "}");
        }
        mixin(readPrimitiveHelper());
    }

    Value readArray(ref Stream stream, char elem_type) {
        string readArrayHelper() {
            char[] cases;
            foreach (c2t; ArrayElementTagValueTypes) {
                cases ~= "case '"~c2t.ch~"':"~
                         "    auto val = new "~c2t.ValueType.stringof~"[length];"~
                         "    foreach (i; 0 .. length) {"~
                         "        stream.read(val[i]);"~
                         "    }"~
                         "    return Value(val);".dup;
            }
            return to!string("switch (elem_type) {"~
                        cases~
                   "    default: throw new UnknownTagTypeException(elem_type ~ \"[]\");"~
                   "}");
        }
        uint length;
        stream.read(length);
        mixin(readArrayHelper());
    }

}
