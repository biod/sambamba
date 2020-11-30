/*
    This file is part of BioD.

    Copyright (C) 2018 Pjotr Prins <pjotr.prins@thebird.nl>
*/

module bio.std.genotype.maf;

import std.algorithm;
import std.array;
import std.conv;
import std.stdio;

/*

   Functions around multi-allelic frequencies (MAF). Allele frequencies are usually
   listed as a range of values between 0.0-1.0 and 0.0-2.0.

*/

/*
   Return (multi-allelele) frequencies of values in gs as an associative array. E.g.

      double[] g2 = [1.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0];

   returns

      double[double] r = [ 0.0: 0.1, 1.0: 0.7, 2.0: 0.2];

   Note you can use any type, so this will work

      assert(maf(["AA", "AB", "BB", "AA", "BB" ]) == [ "AA": 0.4, "BB": 0.4, "AB": 0.2]);
*/

double[T] maf(T)(T[] gs) {
  uint[T] list;
  foreach (g ; gs) {
    list[g] += 1;
  }
  double[T] freq;
  foreach (k ; list.keys)
    freq[k] = 1.0 * list[k] / gs.length;
  return freq;
}

unittest {
  // List of values between 0.0 and 2.0
  double[] g2 = [1.0, 0.0, 1.0, 2.0, 2.0, 1.0, 1.0, 1.0, 1.0, 1.0];
  double[double] r = [ 0.0: 0.1, 1.0: 0.7, 2.0: 0.2];
  assert(maf(g2) == [ 0.0: 0.1, 1.0: 0.7, 2.0: 0.2]);
  // List of values between 0.0 and 1.0
  double[] g1 = array(g2.map!(a => a/2.0));
  assert(maf(g1) == [ 0.0: 0.1, 0.5: 0.7, 1.0: 0.2]);
  assert(maf(["AA", "AB", "BB", "AA", "BB" ]) == [ "AA": 0.4, "BB": 0.4, "AB": 0.2]);
}
