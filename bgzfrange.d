module bgzfrange;

import std.stream;

/**
  Structure representing BGZF block.
 */
struct BgzfBlock {
    // field types are as in the SAM/BAM specification
    // ushort ~ uint16_t, char ~ uint8_t, uint ~ uint32_t

    public ushort bsize; /// total Block SIZE minus one
    public char[] compressed_data = void;
    public uint crc32;
    public uint input_size;
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

    /*
        TODO: throw various exceptions instead of returning false
     */
    bool loadNextBlock() {
        _start_offset = _stream.position;

        if (_stream.eof()) {
            _empty = true; // indicate that range is now empty
            return true;
        }

        try {
            auto bgzf_magic = _stream.readString(4);
            if (bgzf_magic != x"1F 8B 08 04") {
                return false;
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
                        return false; // wrong subfield length
                    }

                    if (found_block_size) {
                        return false; // duplicate field
                    }

                    // read block size
                    _stream.read(bsize); 
                    found_block_size = true;

                    // read the rest
                    _stream.readString(slen - bsize.sizeof); 
                } else {
                    // this subfield has nothing to do with block size, 
                    // just skip
                    _stream.readString(slen);
                }

                len += si1.sizeof + si2.sizeof + slen.sizeof + slen;
            } 

            if (len != gzip_extra_length) {
                return false;
            }

            if (!found_block_size) {
                return false;
            }
           
            // read compressed data
            auto cdata_size = bsize - gzip_extra_length - 19;

            _current_block.bsize = bsize;
            _current_block.compressed_data = _stream.readString(cdata_size);
            
            _stream.read(_current_block.crc32);
            _stream.read(_current_block.input_size);
            
            return true;

        } catch (ReadException e) {
            return false; // TODO: better error handling
        }

        return false;
    }
}
