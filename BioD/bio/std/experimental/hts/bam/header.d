/*
    New style BAM reader. This file is part of Sambamba.
    Copyright (C) 2017,2018 Pjotr Prins <pjotr.prins@thebird.nl>

    Sambamba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published
    by the Free Software Foundation; either version 2 of the License,
    or (at your option) any later version.

    Sambamba is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
    02111-1307 USA

*/

// This is a complete rewrite of Artem Tarasov's original reader.

module bio.std.experimental.hts.bam.header;

/*
import std.conv;
import core.stdc.stdio: fopen, fread, fclose;
import std.typecons;
import std.bitmanip;

import bio.bam.cigar;
*/

import std.exception;
import std.file;
import std.stdio;
import std.string;

// why import this from old bio.bam
// TODO check it depends on undead. 
import bio.std.hts.bam.constants;

import bio.std.experimental.hts.bgzf;
import bio.std.experimental.hts.bgzf_writer;

// what is the difference btw these constants and the ones from bio.std.hts.bam.constants
import bio.std.experimental.hts.constants; 

struct RefSequence {
  size_d length;
  string name;
}

struct BamHeader {
  string id;
  string text;
  RefSequence[] refs;

  @disable this(this); // disable copy semantics;
}

void fetch_bam_header(ref BamHeader header, ref BgzfStream stream) {
  // stderr.writeln("Fetching BAM header");
  ubyte[4] ubyte4;
  stream.read(ubyte4);
  enforce(ubyte4 == BAM_MAGIC,"Invalid file format: expected BAM magic number");
  immutable text_size = stream.read!int();
  // stderr.writeln("Text size ",text_size.sizeof," ",text_size);
  immutable text = stream.read!string(text_size);
  header = BamHeader(BAM_MAGIC,text);
  immutable n_refs = stream.read!int();
  // stderr.writeln("Fetching ",n_refs," references");
  foreach(int n_ref; 0..n_refs) {
    immutable l_name = stream.read!int();
    // stderr.writeln("!!",l_name);
    auto ref_name = stream.read!string(l_name);
    immutable l_ref = stream.read!int(); // length of reference sequence (bps)
    // stderr.writeln(l_name," ",ref_name," ",l_ref);
    header.refs ~= RefSequence(l_ref,ref_name[0..l_name-1]); // drop zero terminator
  }
}

void write_bam_header(ref BgzfWriter bw, ref BamHeader header) {
  // stderr.writeln("Writing BAM header");
  ubyte[4] magic = cast(ubyte[])BAM_MAGIC;
  bw.write(magic);
  // stderr.writeln("Text size ",int.sizeof," ",header.text.length);
  bw.write!int(cast(int)header.text.length);
  bw.write(header.text);
  auto n_refs = cast(int)header.refs.length;
  bw.write!int(cast(int)header.refs.length);
  // stderr.writeln("Writing ",n_refs," references");
  foreach(int n_ref; 0..n_refs) {
    immutable refseq = header.refs[n_ref];
    bw.write!int(cast(int)(refseq.name.length+1));  // incl. zero terminator
    // stderr.writeln("!!",refseq.name.length+1);
    bw.write(refseq.name);
    bw.write!ubyte(cast(ubyte)'\0');
    bw.write!int(cast(int)refseq.length);
    // stderr.writeln(refseq.name.length+1," ",refseq.name," ",refseq.length);
  }
  // stderr.writeln("!!");
  bw.flush_block();
}
