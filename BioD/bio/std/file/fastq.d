/*
   This file is part of BioD.
   Copyright (C) 2016   George Githinji <biorelated@gmail.com>
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

/*
   The bio.core.fastq module is based from a question that I posted on forums.dlang.org.
   I have bundled the answer as a module to support parsing fastq files with D.
   Credit should go to Rikki Cattermole.
 */

module bio.std.file.fastq;

struct FastqRecord {
    const(char)[] id;
    const(char)[] seq;
    const(char)[] qual;

    static auto read(const(char)[] from) {
        struct Result {
            private {
                const(char)[] source;
                FastqRecord value;
                bool isEmpty;
            }

            this(const(char)[] source) {
                this.source = source;
                popFront;
            }

            @property {
                FastqRecord front() {
                    return value;
                }

                bool empty() {
                    return isEmpty;
                }
            }

            void popFront() {
                import std.string : indexOf;

                if (source is null) {
                    isEmpty = true;
                    return;
                }

                void tidyInput() {
                    foreach(i, c; source) {
                        switch(c) {
                            case 0: .. case ' ':
                                    break;
                            default:
                                    source = source[i .. $];
                                    return;
                        }
                    }
                    source = null;
                }

                tidyInput();

                if (source is null)
                    return;

                // id
                assert(source[0] == '@');

                ptrdiff_t len = source.indexOf("\n");
                assert(len > 0);

                value.id = source[1 .. len];
                if (value.id[$-1] == "\r"[0])
                    value.id = value.id[0 .. $-1];

                source = source[len + 1 .. $];

                // seq
                len = source.indexOf("\n");
                assert(len > 0);

                value.seq = source[0 .. len];
                if (value.seq[$-1] == "\r"[0])
                    value.seq = value.seq[0 .. $-1];

                source = source[len + 1 .. $];

                // +id
                len = source.indexOf("\n");
                assert(len > 0);
                source = source[len + 1 .. $];

                // qual
                len = source.indexOf("\n");
                assert(len > 0);

                value.qual = source[0 .. len];
                if (value.qual[$-1] == "\r"[0])
                    value.qual = value.qual[0 .. $-1];

                if (source.length > len + 1) {
                    source = source[len + 1 .. $];
                    tidyInput();
                } else
                    source = null;
            }
        }
        return Result(from);
    }
}

/* fails with ldc >1.11
unittest {
    string input = """
        @seq1
        TTATTTTAAT
        +
        ?+BBB/DHH@
        @seq2
        GACCCTTTGCA
        +
        ?+BHB/DIH@
        @SEQ_ID
        GATTTGGGGTTCAAAGCAGTATCGATCAAATAGTAAATCCATTTGTTCAACTCACAGTTT
        +
        !''*((((***+))%%%++)(%%%%).1***-+*''))**55CCF>>>>>>CCCCCCC65
        """[1 .. $];

    foreach(record; FastqRecord.read(input)) {
        import std.stdio;
        // stderr.writeln(record); -> should be an assert statement
    }
}

*/
