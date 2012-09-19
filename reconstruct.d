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
import std.conv;
import std.range;
import std.traits;

/**
 *  Tries to reconstruct reference sequence from MD tag and CIGAR.
 *  If the read has no MD tag, returns null, otherwise returns reference bases.
 */
string dna(T)(ref T read)
    if (is(Unqual!T == Alignment))
{
    // based on code from Bio::DB::Sam
    
    auto md_tag = read["MD"];
    if (md_tag.is_nothing) {
        return null;
    }

    auto md = cast(string)md_tag;

    auto cigar = read.cigar;

    auto seq = "";
    auto qseq = read.sequence;

    foreach(e; cigar)
    {
        if (e.operation == 'M') {
            seq ~= to!string(qseq[0 .. e.length]);
            qseq = qseq[e.length .. qseq.length];
        } else if (e.operation == 'S' || e.operation == 'I') {
            qseq = qseq[e.length .. qseq.length]; // WTF? $ doesn't work for user-defined types!
        }
    }

    size_t start;
    enum State {
        NOT_STARTED,
        NUMBER,
        START_DELETION,
        DELETION,
        INSERTION
    }

    State current_state = State.NOT_STARTED;
    string result;
    ulong num;
    char base;
    foreach (c; md)
    {
        switch (c) {
            case '0', '1', '2', '3', '4', '5', '6', '7', '8', '9':
                final switch(current_state) {
                    case State.NOT_STARTED:
                    case State.NUMBER:
                        num = num * 10 + cast(int)(c - '0');
                        break;
                    case State.DELETION:
                        result ~= base;
                        num = cast(int)(c - '0');
                        break;
                    case State.INSERTION:
                        result ~= base;
                        start += 1;
                        num = cast(int)(c - '0');
                        break;
                    case State.START_DELETION:
                        break; // FIXME: \^\d is an error
                }
                current_state = State.NUMBER;
                break;
            case '^':
                final switch(current_state) {
                    case State.NOT_STARTED:
                        break;
                    case State.NUMBER:
                        result ~= seq[start .. start + num];
                        start += num;
                        break;
                    case State.DELETION:
                        result ~= base;
                        break;
                    case State.INSERTION:
                        result ~= base;
                        start += 1;
                        break;
                    case State.START_DELETION:
                        break; // FIXME: \^\^ is also an error
                }
                current_state = State.START_DELETION;
                break;
            default: // assume it's a base
                final switch(current_state) {
                    case State.NOT_STARTED:
                        current_state = State.INSERTION;
                        break;
                    case State.NUMBER:
                        result ~= seq[start .. start + num];
                        start += num;
                        current_state = State.INSERTION;
                        break;
                    case State.INSERTION:
                        result ~= base;
                        start += 1;
                        break;
                    case State.START_DELETION:
                        current_state = State.DELETION;
                        break;
                    case State.DELETION:
                        result ~= base;
                        break;
                }
                base = c;
                break;
        }
    }

    final switch (current_state)
    {
        case State.NOT_STARTED:
        case State.START_DELETION:
            break;
        case State.NUMBER:
            result ~= seq[start .. start + num];
            break; 
        case State.INSERTION:
        case State.DELETION:
            result ~= base;
            break;
    }

    return result;
}

unittest {
    // Test reference reconstruction from MD and CIGAR.
    // (Tests are taken from http://davetang.org/muse/2011/01/28/perl-and-sam/)

    Alignment read;

    read = Alignment("r1",
                     "CGATACGGGGACATCCGGCCTGCTCCTTCTCACATG",
                     [CigarOperation(36, 'M')]);
    read["MD"] = "1A0C0C0C1T0C0T27";
    assert(dna(read) == "CACCCCTCTGACATCCGGCCTGCTCCTTCTCACATG");

    read = Alignment("r2",
                     "GAGACGGGGTGACATCCGGCCTGCTCCTTCTCACAT",
                     [CigarOperation(6, 'M'),
                      CigarOperation(1, 'I'),
                      CigarOperation(29, 'M')]);
    read["MD"] = "0C1C0C1C0T0C27";
    assert(dna(read) == "CACCCCTCTGACATCCGGCCTGCTCCTTCTCACAT");

    read = Alignment("r3",
                     "AGTGATGGGGGGGTTCCAGGTGGAGACGAGGACTCC",
                     [CigarOperation(9, 'M'),
                      CigarOperation(9, 'D'),
                      CigarOperation(27, 'M')]);
    read["MD"] = "2G0A5^ATGATGTCA27";
    assert(dna(read) == "AGGAATGGGATGATGTCAGGGGTTCCAGGTGGAGACGAGGACTCC");

    read = Alignment("r4",
                     "AGTGATGGGAGGATGTCTCGTCTGTGAGTTACAGCA",
                     [CigarOperation(2, 'M'),
                      CigarOperation(1, 'I'),
                      CigarOperation(7, 'M'),
                      CigarOperation(6, 'D'),
                      CigarOperation(26, 'M')]);
    read["MD"] = "3C3T1^GCTCAG26";
    assert(dna(read) == "AGGCTGGTAGCTCAGGGATGTCTCGTCTGTGAGTTACAGCA");
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
    static struct Result {
        this(R reads) {
            _reads = reads;
            while (_chunk.length == 0) {
                if (_reads.empty) {
                    _empty = true;
                    return;
                }
                auto read = _reads.front;
                _chunk = dna(read);
                _reference_pos = read.position;
                _reads.popFront();
            }
        }

        @property bool empty() const {
            return _empty;
        }

        @property char front() const {
            if (_bases_to_skip > 0) {
                return 'N';
            }
            return _chunk[0];
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

            _chunk = _chunk[1 .. $];

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
            if (_chunk.length == 0) {
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
                    import std.stdio;
                    writeln("_reference_pos = ", _reference_pos, 
                            "; best_read.position = ", best_read.position,
                            "; _chunk.length = ", _chunk.length);
                }
                // However, we need to strip some bases from the left.
                _chunk = _chunk[_reference_pos - best_read.position .. $];
            }
        }

        private size_t _bases_to_skip;
        private size_t _reference_pos;
        private string _chunk;
        private bool _empty = false;
        private R _reads;
    }

    return Result(reads);
}
