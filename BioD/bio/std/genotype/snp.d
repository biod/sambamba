/*
    This file is part of BioD.

    Copyright (C) 2018 Pjotr Prins <pjotr.prins@thebird.nl>
*/

module bio.std.genotype.snp;

import bio.std.range.splitter;

import std.algorithm;
import std.array;
import std.conv;
import std.exception;
import std.experimental.logger;
import std.stdio;

/*

  Functions around SNPs

*/

struct SNP {
  string name;
  string chr;
  long pos;   // -1 if NA
  double cm;  // centimorgan, -1 if NA
}

immutable NullSNP = SNP("unknown","NA",-1,-1.0);

/*
  Fetch SNP annotations from tab delimited file that looks like name,
  pos, chr, cM

  rs3668922       111771071       13      65.0648
  rs13480515      17261714        10      4.72355
  rs13483034      53249416        17      30.175

  NA values can be used and the last cM column is optional.

  Note: this function should be generalized to return a named tuple, based on a
  named list
*/

SNP[] fetch_snp_annotations(string filen) {
  SNP[] list;

  bool[string] names;

  info("Parsing ",filen);
  File f = File(filen);

  try {
    foreach(line; f.byLine) {
      auto fields = array(SimpleSplitConv!(char[])(line));
      auto name = to!string(fields[0]);
      if (name in names)
        throw new Exception("SNP name "~name~" appears multiple times in "~filen);
      auto pos = (fields[1] == "NA" ? -1 : to!ulong(fields[1]));
      auto chr = to!string(fields[2]);
      double cm = -1.0;
      if (fields.length > 3)
        cm = (fields[3] == "NA" ? -1.0 : to!double(fields[3]));
      auto snp = SNP(name,chr,pos,cm);
      list ~= snp;
    }
  }
  catch(Exception e) {
    writeln(e.msg); // need to test and give proper output FIXME
    throw e;
  }

  return list;
}
