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
module bgzfrange;

public import bgzfblock;
import constants;

import std.stream;
import std.array : uninitializedArray;
import std.conv;

/// Exception type, thrown in case of encountering corrupt BGZF blocks
class BgzfException : Exception {
    this(string msg) { super(msg); }
}

/**
  Range for iterating over BGZF blocks coming from any Stream
 */
struct BgzfRange {

    // /////////////////////////////////////////////////////////////////////////
    //
    // | Here is the general picture of what happens.
    // |
    // | First of all, BgzfRange reads bytes from the stream and determines
    // | boundaries of BGZF blocks. Elements of the range are blocks of 
    // | compressed data, together with their start offsets in the file.
    // | 
    // | After that, blocks are decompressed, and another range comes into play.
    // | Start offsets are still provided together with blocks, because they are 
    // | needed for random access, to calculate virtual offsets of alignments.
    //
    // - Virtual offset is a pair of numbers which uniquely identifies 
    // - location of an individual alignment record in the file. The first is 
    // - the start offset of the bgzf block in which the alignment record begins,
    // - and the second is offset in decompressed data of that block.
    // - Blocks are required to contain no more than 65536 bytes of uncompressed
    // - data, and virtual offsets are stored as uint64_t numbers as follows:
    // - [ {block offset (48 bits)} {offset in decompressed data (16 bits)} ]
    // - 
    // - Relatively small size of BGZF blocks makes for fast random access, 
    // - still allowing good compression (about 3x).
    //
    // | Now that we have range of decompressed blocks, those blocks have to be
    // | made into a proper input stream of bytes. BGZF specification does not
    // | deal with alignment records, it deals with blocks of arbitrary data.
    // | Therefore it's possible that some alignments will get splitted, though
    // | nowadays most software which produces BAM files avoids that.
    // | 
    // | ChunkInputStream joins decompressed blocks into a stream, providing
    // | virtualTell() method, which returns current virtual offset.
    // |
    // | Then, in case of BAM file, SAM header is read first, then comes some
    // | basic information about reference sequences (names and lengths), 
    // | just so as not to duplicate it in alignment records and use integer
    // | indices instead.
    // | 
    // | After that, alignment records come. They are typically quite small,
    // | about 100-1000 bytes depending on sequence length and amount of tags.
    // |
    // | In order to avoid copying memory, and, more importantly, allocating it
    // | so frequently (GC is not well suited for 10000s of allocations/sec), 
    // | the input stream provides yet another method, readSlice, which tries to
    // | return a slice of underlying decompressed block, and only in case it's
    // | impossible it allocates memory which is very rare.
    //  
    // - There're also two possible policies for iterating alignments, 
    // - either packed together with their virtual offsets, or without them.
    // - The first one is used for random access, the second one - for serial.
    //                                                                      
    // -----------------------------------------------------------------------
    //            Picture summarizing the above description:                  
    //                                                                        
    //     Stream      BgzfRange          range of        input     range of  
    //                                  decompressed     stream    alignment  
    //                 each block       BGZF blocks,                records   
    //                  is ~20kB         each ~65kB                           
    //                                                                        
    //     -------      ---------.        ---------     ---------             
    //     |  r  |      |   1st | \       | f|d   |     |  SAM  |             
    //     |  a  |  ->  |  bgzf |\ -----> | i|e   |     | header|             
    //     |  w  |      | block | \       | r|c   |     |-------|             
    //     |     |      |_______|\ \      | s|o   |  -> | r|s   |             
    //     |  b  |      |   2nd | \ \     | t|m   |     | e|e   |             
    //     |  y  |  ->  |  bgzf |  \ ---> |   p   |     | f|q|i |             
    //     |  t  |      | block |   \     |   r   |     | e|u|n |             
    //     |  e  |      |_______|    \    |   e|b |     | r|e|f |             
    //     |  s  |      |   3rd |     \   |   s|l |  -> | e|n|o |             
    //     |  .  |  ->  |  bgzf |      \  |   s|o |     | n|c   |             
    //     |  .  |      | block |       ->|   e|c |     | c|e   |             
    //     |  .  |      |-------|         |   d|k |     | e|s   |             
    //     |  .  |  ->  |  ...  |         |       |     |-------|     --------
    //     |  .  |      |       |         |-------|  -> |records|  -> |------|
    //     |  .  |      |  ...  |         |  ...  |     |  ...  |  -> |------|
    //
    // /////////////////////////////////////////////////////////////////////////

    /**
      Constructs range from stream
     */
    this(Stream stream) {
        _stream = stream;
        _seekable = stream.seekable;
        loadNextBlock();
    }

    /**
        Returns: offset of the start of the current BGZF block
                 in underlying stream. If the stream is non-seekable, 
                 the result is always 0.
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
    bool _seekable = false;

    BgzfBlock _current_block;

    void throwBgzfException(string msg) {
        throw new BgzfException("Error reading BGZF block starting from offset " ~
                                to!string(_start_offset) ~ ": " ~ msg);
    }

    void loadNextBlock() {
        if (_seekable) {
            _start_offset = _stream.position;
        }

        if (_stream.eof()) {
            _empty = true; // indicate that range is now empty
            version(development) {
                import std.stdio;
                stderr.writeln("[info][BGZF range] EOF, current offset is ", _stream.position);
            }
            return;
        }

        try {
            uint bgzf_magic = void;
           
            // TODO: fix byte order if needed
            auto bytes_read = _stream.read((cast(ubyte*)&bgzf_magic)[0 .. 4]);

            if (bytes_read == 0) {
                _empty = true;
                version(development) {
                    import std.stdio;
                    stderr.writeln("[info][BGZF range] end of stream, current offset is ", _stream.position);
                }
                return;
                // TODO: check if last BGZF block was empty, and if not throw a warning
            }

            if (bgzf_magic != BGZF_MAGIC) { 
                throwBgzfException("wrong BGZF magic");
            }
        
            ushort gzip_extra_length = void;

            if (_seekable) {
                _stream.seekCur(uint.sizeof + 2 * ubyte.sizeof);
            } else {
                uint gzip_mod_time = void;
                ubyte gzip_extra_flags = void;
                ubyte gzip_os = void;
                _stream.read(gzip_mod_time);
                _stream.read(gzip_extra_flags);
                _stream.read(gzip_os);
            }

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

                if (si1 == BAM_SI1 && si2 == BAM_SI2) { 
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
                    if (_seekable) {
                        _stream.seekCur(slen - bsize.sizeof);
                    } else {
                        _stream.readString(slen - bsize.sizeof);
                    }
                } else {
                    // this subfield has nothing to do with block size, 
                    // just skip
                    if (_seekable) {
                        _stream.seekCur(slen);
                    } else {
                        _stream.readString(slen);
                    }
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
           
            _current_block.start_offset = start_offset;

            return;

        } catch (ReadException e) {
            throwBgzfException("stream error: " ~ e.msg);
        }

        assert(0);
    }
}
