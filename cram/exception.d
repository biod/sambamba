module cram.exception;

public import std.exception;

class CramException : Exception {
    this(string msg) { super(msg); }
}
