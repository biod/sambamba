/*
    This file is part of Sambamba.
    Copyright (C) 2012-2013    Artem Tarasov <lomereiter@gmail.com>

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
module sambamba.slice;

import bio.bam.reader;
import bio.bam.writer;
import bio.bam.constants;
import bio.core.bgzf.block;
import bio.core.bgzf.compress;
import bio.core.utils.stream;
import bio.core.region;

import std.array;
import std.stream;
import std.getopt;
import std.parallelism;
import std.conv;
import std.stdio;
import std.exception;

/**
    Consider interval [beg .. end) on the reference, where 0 <= beg < end.

    Let R1 = { read from BAM file | read overlaps [beg .. end), read.position < beg },
        R2 = { read from BAM file | beg <= read.position < end },

        s1_start_offset = virtual offset of the first read that overlaps [beg .. end)
        s2_start_offset = virtual offset of the first read which position >= beg
        s2_end_offset = virtual offset of the first read which position >= end

                /\/\/\/\/\/\/\/\/\/\/\/\/#########==============  ...  =============#####      
BGZF blocks   ..........)[...........)[..........)[..........)[.  ...  ....)[.........)[.......
                .                        .                                               .     
                s1_start_offset          s2_start_offset                          s2_end_offset
                                                                                               
                /\/\/\/\/\/\/\/\/\/\/\/\/######                   ...                          
         ...)[..........)[...........)[..........)[..........)[.  ...  ....)[.........)[.......
                .                        .     .                                               
                s1_start_offset   s2_start_offset, s2_end_offset                               

    These numbers are not correctly defined in some cases, so let's extend their
    definitions.

    1) There are no reads that overlap [beg .. end) (R1 and R2 are empty)
        Set all three numbers to the virtual offset of the end of file.
    2) There are no reads that have position >= beg (R2 is empty)
        Set both s2_start_offset and s2_end_offset to the end virtual offset of the last read in R1.
    3) There are reads that have position >= beg, but no one with position >= end
        Set s2_end_offset to the end virtual offset of the last read from current reference.

    Define sets
    S1 = { read | s1_start_offset <= read.start_virtual_offset < s2_start_offset },
    S2 = { read | s2_start_offset <= read.start_virtual_offset < s2_end_offset }

    Notice that R1 is a subset of S1, and R2 = S2.

    Now we divide the algorithm into subcases.

    1) Both R1 and R2 are empty.
        
        Output BAM file with no reads.

    2) R1 is not empty, R2 is empty.

        Output BAM file with reads from R1, creating new BGZF blocks for all of them.

    3) R1 is empty, R2 is not empty.

        First of all, output BAM header and reference sequences information.

        Take first read from R2. Adjust its first BGZF block by chomping s2_start_offset.uoffset
        from the left, and output it. Set start_offset to the start file offset of the next BGZF block.
        Take last read from R2. Adjust its last BGZF block by chomping everything after the end of
        the alignment record. Set end_offset to the start file offset of this BGZF block.
        
        Output first adjusted block, then copy of file since start_offset till end_offset, then
        second adjusted block.

        (It may turn out that after chomping some of blocks are empty. In this case we skip them.)

    4) R1 and R2 are not empty.

        Since R1 and R2 are disjoint, combine approaches 2) and 3)
*/


void fetchRegion(BamReader bam, Region region, ref Stream stream)
{
    auto chr = region.reference;
    auto beg = region.beg;
    auto end = region.end;

    auto filename = bam.filename;

    auto reads1 = bam[chr][beg .. end];

    auto eof_offset = bam.eofVirtualOffset();

    VirtualOffset s1_start_offset;
    VirtualOffset s2_start_offset;
    VirtualOffset s2_end_offset = eof_offset;

    if (reads1.empty) {
        s1_start_offset = s2_start_offset = eof_offset;
    } else {
        s1_start_offset = s2_start_offset = reads1.front.start_virtual_offset;
    }

    // Set R1
    auto r1 = appender!(typeof(reads1.front)[])();

    while (!reads1.empty) {
        auto front = reads1.front;
        if (front.position < beg) {
            r1.put(front);
            s2_start_offset = front.end_virtual_offset;
            reads1.popFront();
        } else {
            break;
        }
    }

    auto reads2 = bam[chr][end .. uint.max];

    if (reads1.empty && reads2.empty) {
        // are there any reads with position >= beg?
        s2_end_offset = s2_start_offset;
    } else if (!reads1.empty && reads2.empty) { 
        // are there any reads with position >= end?
        s2_end_offset = bam[chr].endVirtualOffset();
    } else {
        foreach (read; reads2) {
            s2_end_offset = read.start_virtual_offset;
            if (read.position >= end) {
                break;
            }
        }
    }

    // write header and reference sequence information
    auto writer = new BamWriter(stream);
    writer.writeSamHeader(bam.header);
    writer.writeReferenceSequenceInfo(bam.reference_sequences);
    // write R1
    foreach (read; r1.data)
        writer.writeRecord(read);
    writer.flush();

    // R2 is non-empty
    if (s2_start_offset < s2_end_offset) {
       
        // Either R2 is fully contained in one BGZF block...
        if (s2_start_offset.coffset == s2_end_offset.coffset) {
            // write chomped block
            auto block = bam.getBgzfBlockAt(s2_start_offset.coffset);
            auto data = decompressBgzfBlock(block).decompressed_data;
            data = data[s2_start_offset.uoffset .. s2_end_offset.uoffset];
            stream.write(bgzfCompress(data, -1));
        } else { // ...or it spans several of them.

            // write left chomped block
            auto block1 = bam.getBgzfBlockAt(s2_start_offset.coffset);
            auto copy_start_offset = block1.end_offset;

            auto data1 = decompressBgzfBlock(block1).decompressed_data;
            data1 = data1[s2_start_offset.uoffset .. $];
            if (data1.length > 0) {
                stream.write(bgzfCompress(data1, -1));
            }

            auto copy_end_offset = s2_end_offset.coffset;

            // copy blocks in between

            ubyte[8192] copy_buffer;
            auto file_stream = new bio.core.utils.stream.File(filename);
            file_stream.seekSet(copy_start_offset);

            size_t curpos = cast(size_t)copy_start_offset;
            while (copy_end_offset - curpos > 8192) {
                file_stream.readExact(copy_buffer.ptr, 8192);
                stream.writeExact(copy_buffer.ptr, 8192);
                curpos += 8192;
            }

            file_stream.readExact(copy_buffer.ptr, cast(size_t)copy_end_offset - curpos);
            stream.writeExact(copy_buffer.ptr, cast(size_t)copy_end_offset - curpos);

            // write right chomped block if it's non-empty
            auto block2 = bam.getBgzfBlockAt(s2_end_offset.coffset);

            auto data2 = decompressBgzfBlock(block2).decompressed_data;
            data2 = data2[0 .. s2_end_offset.uoffset];

            if (data2.length > 0) {
                stream.write(bgzfCompress(data2, -1));
            }
        }
    }

    // write EOF
    stream.writeExact(BAM_EOF.ptr, BAM_EOF.length);
}


void printUsage()
{
    stderr.writeln("Usage: sambamba-slice [options] <input.bam> <region>");
    stderr.writeln();
    stderr.writeln("       Region is given in standard form ref:beg-end.");
    stderr.writeln("       Output is to STDOUT unless output filename is specified.");
    stderr.writeln();
    stderr.writeln("OPTIONS: -o, --output-filename=OUTPUT_FILENAME");
    stderr.writeln("            output BAM filename");
}

version(standalone) {
    int main(string[] args) {
        return slice_main(args);
    }
}

int slice_main(string[] args) {

    // at least two arguments must be presented
    string output_filename = null;

    try {
        getopt(args,
               std.getopt.config.caseSensitive,
               "output-filename|o", &output_filename);

        if (args.length < 3) {
            printUsage();
            return 0;
        }

        auto bam = new BamReader(args[1]);
        auto region = parseRegion(args[2]);

        Stream stream;
        scope(exit) stream.close();

        if (output_filename != null) {
            stream = new std.stream.BufferedFile(output_filename, FileMode.OutNew);
        } else {
            immutable BUFSIZE = 1_048_576;
            stream = new BufferedFile(stdout.fileno, FileMode.Out, BUFSIZE);
        }

        fetchRegion(bam, region, stream);

    } catch (Exception e) {
        stderr.writeln("sambamba-slice: ", e);
        version (development) {
            throw e;
        }

        return 1;
    }
    return 0;
}
