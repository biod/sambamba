/*
    New style BAM writer. This file is part of Sambamba.
    Copyright (C) 2017 Pjotr Prins <pjotr.prins@thebird.nl>

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

module sambamba.bio2.bam.writer;

import std.conv;
import core.stdc.stdio: fopen, fread, fclose;
import std.exception;
import std.file;
import std.stdio;
import std.string;
import std.typecons;
import std.bitmanip;

import bio.bam.cigar;
import bio.bam.constants;

import sambamba.bio2.bgzf;
import sambamba.bio2.bgzf_writer;
import sambamba.bio2.constants;

import sambamba.bio2.bam.header;

struct BamWriter {
  BgzfWriter bgzf_writer;

  this(string fn, ref BamHeader header) {
    bgzf_writer = BgzfWriter(fn,9);
    write_bam_header(bgzf_writer,header);
  }

}
