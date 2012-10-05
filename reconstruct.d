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
module reconstruct;

import alignment;
import md.core;

import std.conv;
import std.range;
import std.traits;
import std.algorithm;
import std.range;

/// Reconstruct read DNA.
/// Returns lazy sequence.
auto dna(T)(T read) 
    if (is(Unqual!T == Alignment))
{

    debug {
        /*
        import std.stdio;
        stderr.writeln("[dna] processing read ", read.read_name);
        stderr.flush();
        */
    }

    static struct QueryChunk(S) {
        S sequence;
        CigarOperation operation;
    }

    // Get read sequence chunks corresponding to query-consuming operations in read.sequence
    static auto queryChunks(ref T read) {
        static struct Result(R, S) {
            this(R ops, S seq) {
                _seq = seq;
                _ops = ops;
            }
        
            auto front() @property {
                auto op = _ops.front;
                return QueryChunk!S(_seq[0 .. op.length], op);
            }

            bool empty() @property {
                return _ops.empty;    
            }

            void popFront() {
                _seq = _seq[_ops.front.length .. _seq.length];
                _ops.popFront();
            }

            private R _ops;
            private S _seq;
        }

        static auto getResult(R, S)(S sequence, R cigar) {
            return Result!(R, S)(cigar, sequence);
        }

        return getResult(read.sequence, filter!"a.is_query_consuming"(read.cigar));
    }

    auto _read = read;

    auto query_chunks = queryChunks(_read);

    static struct Result(R, M) {
        this(ref T read, R query_sequence, M md_operations) {
            debug {
                _initial_qseq = to!string(query_sequence);
            }
            _qseq = query_sequence; 
            _md = md_operations;
            _fetchNextMdOperation();
        }

        bool empty() @property {
            return _empty;            
        }

        /*
        MD operations -> match(N)    ? consume N characters from query
                         mismatch(C) ? consume a character from query and replace it with C
                         deletion(S) ? consume S from MD
        */

        char front() @property {
            final switch (_cur_md_op.type) {
                case MdOperationType.Match:
                    return cast(char)_qseq.front;
                case MdOperationType.Mismatch:
                    return _cur_md_op.mismatch;
                case MdOperationType.Deletion:
                    return cast(char)_cur_md_op.deletion.front;
            }
        }

        private void _fetchNextMdOperation() {
            if (_md.empty) {
                _empty = true;
                return;
            }
            _cur_md_op = _md.front;
            _md.popFront();
        }

        private bool _qseqIsSuddenlyEmpty() {
            if (!_qseq.empty) {
                return false;
            }

            /* MD and CIGAR don't match */
            debug {
                import std.stdio;
                stderr.writeln("Current MD operation: ", _cur_md_op);
                stderr.writeln("Query sequence: ", _initial_qseq);
            }

            return true;
        }

        void popFront() {
            final switch (_cur_md_op.type) {
                case MdOperationType.Mismatch:
                    if (_qseqIsSuddenlyEmpty())
                        break;
                    _qseq.popFront();
                    _fetchNextMdOperation();
                    break;
                case MdOperationType.Match:
                    if (_qseqIsSuddenlyEmpty())
                        break;
                    --_cur_md_op.match;
                    _qseq.popFront();
                    if (_cur_md_op.match == 0) {
                        _fetchNextMdOperation();
                    }
                    break;
                case MdOperationType.Deletion:
                    _cur_md_op.deletion = _cur_md_op.deletion[1..$];
                    if (_cur_md_op.deletion.empty) {
                        _fetchNextMdOperation();
                    }
                    break;
            }
        }

        private {
            debug {
                string _initial_qseq;
            }
            R _qseq;
            M _md;
      
            bool _empty;
            MdOperation _cur_md_op;
        }
    }
  
    auto md = _read["MD"];
    string md_str;
    if (!md.is_nothing) {
        md_str = cast(string)_read["MD"];
    }
    
    static auto getResult(R, M)(ref T read, R query, M md_ops) {
        return Result!(R, M)(read, query, md_ops);
    }

    auto result =  getResult(_read, 
                             joiner(map!"a.sequence"(filter!"a.operation.is_reference_consuming"(query_chunks))),
                             mdOperations(md_str));

    debug {
        import std.stdio;
        if (result.empty) {
            stderr.writeln("[dna] empty DNA!");
            stderr.writeln("      read name: ", read.read_name);
            stderr.writeln("      read sequence: ", read.sequence);
            stderr.writeln("      read CIGAR: ", read.cigarString());
            stderr.writeln("      read MD tag: ", read["MD"]);
            stderr.flush();
        }
    }

    return result;
}

unittest {

    import std.stdio;
    writeln("Testing reconstruction of reference from MD tags and CIGAR");

    // Test reference reconstruction from MD and CIGAR.
    // (Tests are taken from http://davetang.org/muse/2011/01/28/perl-and-sam/)

    Alignment read;

    read = Alignment("r1",
                     "CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG",
                     [CigarOperation(36, 'M')]);
    read["MD"] = "1A0C0C0C1T0C0T27";

    assert(equal(dna(read), "CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG"));

    read = Alignment("r2",
                     "GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT",
                     [CigarOperation(6, 'M'),
                      CigarOperation(1, 'I'),
                      CigarOperation(29, 'M')]);
    read["MD"] = "0C1C0C1C0T0C27";

    assert(equal(dna(read), "CACCCCTCTGACATCCGGCCTGCTCCTTCTCACAT"));

    read = Alignment("r3",
                     "AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC",
                     [CigarOperation(9, 'M'),
                      CigarOperation(9, 'D'),
                      CigarOperation(27, 'M')]);
    read["MD"] = "2G0A5^ATGATGTCA27";
    assert(equal(dna(read), "AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC"));

    read = Alignment("r4",
                     "AGTGATGGGAGGATGTCTCGTCTGTGAGTTACAGCA",
                     [CigarOperation(2, 'M'),
                      CigarOperation(1, 'I'),
                      CigarOperation(7, 'M'),
                      CigarOperation(6, 'D'),
                      CigarOperation(26, 'M')]);
    read["MD"] = "3C3T1^GCTCAG26";
    assert(equal(dna(read), "AGGCTGGTAGCTCAGGGATGTCTCGTCTGTGAGTTACAGCA"));

}

/**
 * Returns lazy sequence of reference bases. If some bases can't be determined from reads,
 * they are replaced with 'N'.
 *
 * Reads must be a range of reads aligned to the same reference sequence, sorted by leftmost
 * coordinate.
 * Returned reference bases start from the leftmost position of the first read,
 * and end at the rightmost position of all the reads.
 */
auto dna(R)(R reads)
    if (isInputRange!R && is(Unqual!(ElementType!R) == Alignment))
{
    static struct Result(F) {
        this(F reads) {
            _reads = reads;
            if (_reads.empty) {
                _empty = true;
                return;
            }
            auto read = _reads.front;
            _chunk = dna(read);
            _reference_pos = read.position;
            _reads.popFront();
        }

        @property bool empty() {
            return _empty;
        }

        @property char front() {
            if (_bases_to_skip > 0) {
                return 'N';
            }
            return _chunk.front;
        }

        private void setSkipMode(ref Alignment read) {
            _reads.popFront();
            _chunk = dna(read);
            _bases_to_skip = read.position - _reference_pos;
        }

        void popFront() {
            _reference_pos += 1;

            if (_bases_to_skip > 0) {
                --_bases_to_skip;
                return;
            }

            _chunk.popFront();

            /*
             * If current chunk is empty, get the next one.
             *                                                                  
             * Here's the reference:                                            
             * .........................*.......................................
             *                          _reference_pos (we are here)            
             * Last chunk ended just now:                                       
             *              [..........]                                        
             * Go through subsequent reads while their leftmost position is     
             * less or equal to _reference_pos, select the one which covers     
             * more bases to the right of _reference_pos.                       
             *               [...............]                                  
             *                [....]                                            
             *                  [..........]                                    
             *                        [.........]  <- this one is the best      
             */
            if (_chunk.empty) {
                if (_reads.empty) {
                    _empty = true;
                    return;
                }
                auto next_read = _reads.front;
                if (next_read.position > _reference_pos) {
                    setSkipMode(next_read);
                    return;
                }
                auto best_read = next_read;
                // read covers half-open [position .. position + basesCovered) interval
                auto best_end_pos = best_read.basesCovered() + best_read.position;
                bool found_good = best_end_pos > _reference_pos;
                while (true) {
                    if (_reads.empty) {
                        if (!found_good) {
                            _empty = true;
                            return;
                        }
                        break;
                    }

                    auto read = _reads.front;

                    if (read.position > _reference_pos) {
                        if (!found_good) {
                            setSkipMode(read);
                            return;
                        }
                        break;
                    }

                    auto end_pos = read.basesCovered() + read.position;
                    if (end_pos > _reference_pos) {
                        found_good = true;
                        if (end_pos > best_end_pos) {
                            best_end_pos = end_pos;
                            best_read = read;
                        }
                    }
                    _reads.popFront();
                }

                // If we're here, we've found a good read.
                _chunk = dna(best_read);
                debug {
                    /*
                    import std.stdio;
                    writeln("_reference_pos = ", _reference_pos, 
                            "; best_read.position = ", best_read.position,
                            "; _chunk length = ", best_read.basesCovered());
                            */
                }
                // However, we need to strip some bases from the left.
                popFrontN(_chunk, _reference_pos - best_read.position);
            }
        }

        private size_t _bases_to_skip;
        private size_t _reference_pos;
        private typeof(dna(_reads.front)) _chunk;
        private bool _empty = false;
        private F _reads;
    }

    auto nonempty = filter!"a.basesCovered() > 0"(reads);
    return Result!(typeof(nonempty))(nonempty);
}

unittest {

    // reads are taken from HG00110.chrom20.ILLUMINA.bwa.GBR.exome.20111114.bam 

    auto r1 = Alignment("r1",
                        "AGGTTTTGTGAGTGGGACAGTTGCAGCAAAACACAACCATAGGTGCCCATCCACCAAGGCAGGCTCTCCATCTTGCTCAGAGTGGCTCTA",
                        [CigarOperation(89, 'M'),
                         CigarOperation(1, 'S')]);
    r1.position = 60246;
    r1["MD"] = "89";

    auto r2 = Alignment("r2",
                        "TGTGAGTGGGACAGTTGCAGCAAAACACAACCATAGGTGCCCATCCACCAAGGCAGGCTCTCCATCTTGCTCAGAGTGGCTCCAGCCCTT",
                        [CigarOperation(83, 'M'),
                         CigarOperation(7, 'S')]);
    r2.position = 60252;
    r2["MD"] = "82T0";
    
    auto r3 = Alignment("r3",
                        "CATAGGTGCCCATCCACCAAGGCAGGCTCTCCATCTTGCTCAGAGTGGCTCTAGCCCTTGCTGACTGCTGGGCAGGGAGAGAGCAGAGCT",
                        [CigarOperation(90, 'M')]);
    r3.position = 60283;
    r3["MD"] = "90";

    auto r4 = Alignment("r4",
                        "CCCTTGCTGACTGCTGGGCAGGGAGAGAGCAGAGCTAACTTCCTCATGGGACCTGGGTGTGTCTGATCTGTGCACACCACTATCCAACCG",
                        [CigarOperation(90, 'M')]);
    r4.position = 60337;
    r4["MD"] = "90";

    auto r5 = Alignment("r5",
                        "GAGGCTCCACCCTGGCCACTCTTGTGTGCACACAGCACAGCCTCTACTGCTACACCTGAGTACTTTGCCAGTGGCCTGGAAGCACTTTGT",
                        [CigarOperation(90, 'M')]);
    r5.position = 60432;
    r5["MD"] = "90";

    auto reads = [r1, r2, r3, r4, r5];
    assert(equal(dna(reads), "AGGTTTTGTGAGTGGGACAGTTGCAGCAAAACACAACCATAGGTGCCCATCCACCAAGGCAGGCTCTCCATCTTGCTCAGAGTGGCTCTAGCCCTTGCTGACTGCTGGGCAGGGAGAGAGCAGAGCTAACTTCCTCATGGGACCTGGGTGTGTCTGATCTGTGCACACCACTATCCAACCGNNNNNGAGGCTCCACCCTGGCCACTCTTGTGTGCACACAGCACAGCCTCTACTGCTACACCTGAGTACTTTGCCAGTGGCCTGGAAGCACTTTGT"));
}
