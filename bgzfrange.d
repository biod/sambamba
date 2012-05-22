module bgzfrange;

import std.stream;
import std.array : uninitializedArray;
import std.conv;

/**
  Structure representing BGZF block.
 */
struct BgzfBlock {
    // field types are as in the SAM/BAM specification
    // ushort ~ uint16_t, char ~ uint8_t, uint ~ uint32_t

    public ushort bsize; /// total Block SIZE minus one
    public ubyte[] compressed_data = void;
    public uint crc32;
    public uint input_size; /// size of uncompressed data
}

/// Exception type, thrown in case of encountering corrupt BGZF blocks
class BgzfException : Exception {
    this(string msg) { super(msg); }
}

/**
  Class for iterating over BGZF blocks coming from any Stream
 */
class BgzfRange {

    /**
      Constructs range from stream
     */
    this(Stream stream) {
        _stream = stream;
        loadNextBlock();
    }

    /**
        Returns: offset of the start of the current BGZF block
                 in underlying stream
     */
    @property ulong start_offset() { return _start_offset; }

    bool empty() {
        return _empty;
    }

    void popFront() {
        loadNextBlock();
    }

    BgzfBlock front() {
        return _current_block;
    }

private:
    Stream _stream;
    ulong _start_offset;

    bool _empty = false;

    BgzfBlock _current_block;

    void throwBgzfException(string msg) {
        throw new BgzfException("Error reading BGZF block starting from offset " ~
                                to!string(_start_offset) ~ ": " ~ msg);
    }

    void loadNextBlock() {
        _start_offset = _stream.position;

        if (_stream.eof()) {
            _empty = true; // indicate that range is now empty
            return;
        }

        try {
            uint bgzf_magic = void;
            _stream.read(bgzf_magic);
            if (bgzf_magic != 0x04_08_8B_1F) { // little endian
                throwBgzfException("wrong BGZF magic");
            }
            
            uint gzip_mod_time = void;
            ubyte gzip_extra_flags = void;
            ubyte gzip_os = void;
            ushort gzip_extra_length = void;

            _stream.read(gzip_mod_time);
            _stream.read(gzip_extra_flags);
            _stream.read(gzip_os);
            _stream.read(gzip_extra_length);
          
            ushort bsize = void; // total Block SIZE minus 1
            bool found_block_size = false;

            // read extra subfields
            size_t len = 0;
            while (len < gzip_extra_length) {
                ubyte si1 = void;    // Subfield Identifier1
                ubyte si2 = void;    // Subfield Identifier2
                ushort slen = void;  // Subfield LENgth
                
                _stream.read(si1);    
                _stream.read(si2);    
                _stream.read(slen);   

                if (si1 == 66 && si2 == 67) { 
                    // found 'BC' as subfield identifier
                    if (slen != 2) {
                        throwBgzfException("wrong BC subfield length: " ~ 
                                           to!string(slen) ~ "; expected 2");
                    }

                    if (found_block_size) {
                        throwBgzfException("duplicate field with block size");
                    }

                    // read block size
                    _stream.read(bsize); 
                    found_block_size = true;

                    // skip the rest
                    _stream.seekCur(slen - bsize.sizeof);
                } else {
                    // this subfield has nothing to do with block size, 
                    // just skip
                    _stream.seekCur(slen);
                }

                len += si1.sizeof + si2.sizeof + slen.sizeof + slen;
            } 

            if (len != gzip_extra_length) {
                throwBgzfException("total length of subfields in bytes (" ~ 
                                   to!string(len) ~ 
                                   ") is not equal to gzip_extra_length (" ~
                                   to!string(gzip_extra_length) ~ ")");
            }

            if (!found_block_size) {
                throwBgzfException("block size was not found in any subfield");
            }
           
            // read compressed data
            auto cdata_size = bsize - gzip_extra_length - 19;

            _current_block.bsize = bsize;
            _current_block.compressed_data = uninitializedArray!(ubyte[])(cdata_size);
            _stream.readExact(_current_block.compressed_data.ptr, cdata_size);
            
            _stream.read(_current_block.crc32);
            _stream.read(_current_block.input_size);
            
            return;

        } catch (ReadException e) {
            throwBgzfException("stream error: " ~ e.msg);
        }

        assert(0);
    }
}
