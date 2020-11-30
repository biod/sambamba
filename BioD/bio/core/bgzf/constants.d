/*
    This file is part of BioD.
    Copyright (C) 2012-2016    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.core.bgzf.constants;

immutable ubyte[4] BGZF_MAGIC = [0x1F, 0x8B, 0x8, 0x4];

immutable ubyte[16] BLOCK_HEADER_START =
    [ 31, 139,   8,   4,  // BGZF magic
       0,   0,   0,   0,  // GZIP modification time
       0,                 // GZIP extra flags
     255,                 // GZIP OS identifier
       6,   0,            // GZIP extra length == 6 (LE)
      66,  67,            // Subfield 'BC'
       2,   0];           // Subfield length (holds 1 ushort)

// empty block
immutable ubyte[28] BGZF_EOF =
    [31, 139, 8, 4,
        0, 0, 0, 0,
                 0,
               255,
              6, 0,
            66, 67,
              2, 0,
             27, 0,
              3, 0,
        0, 0, 0, 0,
        0, 0, 0, 0];



// BGZF block header length in bytes.
// Block header holds BLOCK_HEADER_START + block size (ushort)
immutable BLOCK_HEADER_LENGTH = BLOCK_HEADER_START.length + ushort.sizeof;

// BGZF footer holds CRC32 and size of decompressed block.
immutable BLOCK_FOOTER_LENGTH = uint.sizeof + uint.sizeof;

immutable BGZF_MAX_BLOCK_SIZE = 65536;
immutable BGZF_BLOCK_SIZE = 0xFF00;
