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

/**
 *  Tries to reconstruct reference sequence from MD tag and CIGAR.
 *  If the read has no MD tag, returns null, otherwise returns reference bases.
 */
string dna(ref Alignment read)
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
