/*
    This file is part of Sambamba.
    Copyright (C) 2012-2015    Artem Tarasov <lomereiter@gmail.com>

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
module utils.strip_bcf_header;
import undead.cstream : CFile, FileMode;
import std.exception;
import std.algorithm : min;

// used in mpileup tool to strip the common header from partial files
void stripBcfHeader(File input_file, File output_file) {
  import bio.core.bgzf.inputstream;
  import bio.core.bgzf.constants;
  import bio.core.bgzf.block;
  import bio.core.bgzf.compress;

  auto stream = new CFile(input_file.getFP(), FileMode.In);
  auto supplier = new StreamSupplier(stream);
  BgzfBlock block;
  ubyte[] tmp = new ubyte[BGZF_MAX_BLOCK_SIZE];
  ushort skip_start, skip_end;

  supplier.getNextBgzfBlock(&block, tmp.ptr, &skip_start, &skip_end);
  auto decompressed = decompressBgzfBlock(block);
  auto data = decompressed.decompressed_data;
  enforce(cast(char[])data[0 .. 3] == "BCF");

  auto header_len = *(cast(int*)(data.ptr + 5));
  size_t offset = 9 + header_len;
  header_len -= min(data.length - 9, header_len);
  while (header_len > 0) {
    supplier.getNextBgzfBlock(&block, tmp.ptr, &skip_start, &skip_end);
    decompressed = decompressBgzfBlock(block);
    data = decompressed.decompressed_data;
    offset = header_len;
    header_len -= min(data.length, header_len);
  }

  auto first_block = bgzfCompress(data[offset .. $], -1, tmp);
  output_file.rawWrite(first_block);

  while (true) {
    auto input = input_file.rawRead(tmp);
    if (input.length == 0) break;
    output_file.rawWrite(input);
  }
}

void stripUncompressedBcfHeader(File input_file, File output_file) {
  import std.exception, std.algorithm;
  ubyte[4096] buffer = void;
  auto start = input_file.rawRead(buffer[0 .. 5]);
  enforce(cast(char[])start[0 .. 3] == "BCF");

  auto header_len_ = input_file.rawRead(buffer[0 .. 4]);

  auto header_len = *(cast(int*)(buffer.ptr));
  while (header_len > 0) {
    auto bytes_to_read = min(header_len, buffer.length);
    auto chunk = input_file.rawRead(buffer[0 .. bytes_to_read]);
    header_len -= chunk.length;
  }

  while (true) {
    auto chunk = input_file.rawRead(buffer[]);
    if (chunk.length == 0)
      break;
    output_file.rawWrite(chunk);
  }
}

void stripVcfHeader(File input_file, File output_file) {
  import bio.core.utils.bylinefast;
  bool header = true;
  auto w = output_file.lockingTextWriter;
  foreach (line; ByLineFast(input_file)) {
    if (header && line.length > 0 && line[0] == '#') // header
      continue;
    header = false;
    w.put(line);
    w.put('\n');
  }
}

import std.stdio, std.conv;
int strip_bcf_header_main(string[] args) {
  if (args[1] == "--bcf")
    stripBcfHeader(stdin, stdout);
  else if (args[1] == "--ubcf")
    stripUncompressedBcfHeader(stdin, stdout);
  else if (args[1] == "--vcf")
    stripVcfHeader(stdin, stdout);
  return 0;
}
