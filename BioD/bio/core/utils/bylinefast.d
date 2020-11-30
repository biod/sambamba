/// Kudos to Juan Manuel Cabo
/// http://forum.dlang.org/post/cddkatcqmdtibcmfljff@forum.dlang.org
///
/// This piece of code is in public domain.
module bio.core.utils.bylinefast;

import std.stdio;
import std.string: indexOf;
import core.stdc.string: memmove;

/**
   Reads by line in an efficient way (10 times faster than File.byLine from std.stdio).
   This is accomplished by reading entire buffers (fgetc() is not used),
   and allocating as little as possible.

   The char \n is considered as separator, removing the previous \r if it exists.

   The \n is never returned. The \r is not returned if it was
   part of a \r\n (but it is returned if it was by itself).

   The returned string is always a substring of a temporary
   buffer, that must not be stored. If necessary, you must
   use str[] or .dup or .idup to copy to another string.

   Example:

   File f = File("file.txt");
   foreach (string line; ByLineFast(f)) {
       ...process line...
       //Make a copy:
       string copy = line[];
   }

   The file isn't closed when done iterating, unless it was the only reference to
   the file (same as std.stdio.byLine). (example: ByLineFast(File("file.txt"))).
*/
struct ByLineFast {
    File file;
    char[] line;
    bool first_call = true;
    char[] buffer;
    char[] strBuffer;

    this(File f, int bufferSize=4096) {
        assert(bufferSize > 0);
        file = f;
        buffer.length = bufferSize;
    }

    @property bool empty() const {
        //Its important to check "line !is null" instead of
        //"line.length != 0", otherwise, no empty lines can
        //be returned, the iteration would be closed.
        if (line !is null) {
            return false;
        }
        if (!file.isOpen) {
            //Clean the buffer to avoid pointer false positives:
            (cast(char[])buffer)[] = 0;
            return true;
        }

        //First read. Determine if it's empty and put the char back.
            auto mutableFP = (cast(File*) &file).getFP();
        auto c = fgetc(mutableFP);
        if (c == -1) {
            //Clean the buffer to avoid pointer false positives:
            (cast(char[])buffer)[] = 0;
            return true;
        }
        if (ungetc(c, mutableFP) != c) {
            assert(false, "Bug in cstdlib implementation");
        }
        return false;
    }

    @property char[] front() {
        if (first_call) {
            popFront();
            first_call = false;
        }
        return line;
    }

    void popFront() {
        if (strBuffer.length == 0) {
            strBuffer = file.rawRead(buffer);
            if (strBuffer.length == 0) {
                file.detach();
                line = null;
                return;
            }
        }

        long pos = strBuffer.indexOf('\n');
        if (pos != -1) {
            if (pos != 0 && strBuffer[cast(size_t)pos-1] == '\r') {
                line = strBuffer[0 .. cast(size_t)(pos-1)];
            } else {
                line = strBuffer[0 .. cast(size_t)pos];
            }
            //Pop the line, skipping the terminator:
            strBuffer = strBuffer[cast(size_t)(pos+1) .. $];
        } else {
            //More needs to be read here. Copy the tail of the buffer
            //to the beginning, and try to read with the empty part of
            //the buffer.
            //If no buffer was left, extend the size of the buffer before
            //reading. If the file has ended, then the line is the entire
            //buffer.

            if (strBuffer.ptr != buffer.ptr) {
                //Must use memmove because there might be overlap
                memmove(buffer.ptr, strBuffer.ptr,
                        strBuffer.length * char.sizeof);
            }
            auto spaceBegin = strBuffer.length;
            if (strBuffer.length == buffer.length) {
                //Must extend the buffer to keep reading.
                assumeSafeAppend(buffer);
                buffer.length = buffer.length * 2;
            }
            char[] readPart = file.rawRead(buffer[spaceBegin .. $]);
            if (readPart.length == 0) {
                //End of the file. Return whats in the buffer.
                //The next popFront() will try to read again, and then
                //mark empty condition.
                if (spaceBegin != 0 && buffer[spaceBegin-1] == '\r') {
                    line = buffer[0 .. spaceBegin-1];
                } else {
                    line = buffer[0 .. spaceBegin];
                }
                strBuffer = null;
                return;
            }
            strBuffer = buffer[0 .. spaceBegin + readPart.length];
            //Now that we have new data in strBuffer, we can go on.
            //If a line isn't found, the buffer will be extended again to read more.
            popFront();
        }
    }
}
