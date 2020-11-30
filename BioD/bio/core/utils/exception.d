/*
    Asserte functions. This file is part of BioD.
    Copyright (C) 2018 Pjotr Prins <pjotr.prins@thebird.nl>
*/

// This code is based on 'exception.d' in D's Phobos

module bio.core.utils.exception;

import std.exception;
import std.traits;

/++
    Asserts that the given value is true, but unlike standard assert
    throws an exception on error.

    Params:
        value = The value to test.
        ex = The exception to throw if the value evaluates to false.

    Returns: $(D value), if `cast(bool)value` is true. Otherwise, $(D ex) is
    thrown.

    Example:
    --------------------
    auto f = asserte(fopen("data.txt"));
    auto line = readln(f);
    asserte(line.length, new IOException); // expect a non-empty line
    --------------------
 +/
T asserte(T)(T value, lazy Throwable ex)
{
  version(assert) {
    if (!value) throw ex();
  }
  return value;
}

T asserte(T)(T value, lazy string msg = "asserte failed")
{
  version(assert) {
    if (!value) throw new Exception(msg);
  }
  return value;
}

/++
    Asserts that the given value is true, but unlike standard assert
    throws an exception on error.

    Params:
        value = The value to test.
        dg = The delegate to be called if the value evaluates to false.
        file = The source file of the caller.
        line = The line number of the caller.

    Returns: $(D value), if `cast(bool)value` is true. Otherwise, the given
    delegate is called.

    The safety and purity of this function are inferred from $(D Dg)'s safety
    and purity.
 +/
T asserte(T, Dg, string file = __FILE__, size_t line = __LINE__)
    (T value, scope Dg dg)
    if (isSomeFunction!Dg && is(typeof( dg() )) &&
        is(typeof({ if (!value) {} })))
{
  version(assert) {
    if (!value) dg();
  }
  return value;
}
