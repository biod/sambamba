module bio.std.hts.bam.md.parse;

import bio.std.hts.bam.md.operation;
import std.ascii;
import std.array;
import std.algorithm;
import std.functional;
import std.range;
import std.conv;
import std.traits;

/// Returns bidirectional range of MD operations. Zero matches are skipped.
auto mdOperations(string md) {

    static struct Result {
        private {
            string _md = void;
            MdOperation _cached_front = void;
            MdOperation _cached_back = void;
            ubyte _rem = 255;
        }

        this(string md) {
            _md = md;
            if (!cacheFront()) {
                _rem = 0;
            } else {
                if (!cacheBack()) {
                    _cached_back = _cached_front;
                    _rem = 1;
                }
            }
        }
      
        bool empty() @property {
            return _rem == 0;
        }

        Result save() @property {
            Result res = void;
            res._md = _md;
            res._cached_front = _cached_front;
            res._cached_back = _cached_back;
            res._rem = _rem;
            return res;
        }

        ref MdOperation front() @property {
            return _cached_front;
        }

        ref MdOperation back() @property {
            return _cached_back;
        }

        void popFront() {
            if (_md.empty) {
                if (_rem == 255) {
                    _cached_front = _cached_back;
                    _rem = 1;
                } else {
                    _rem = 0;
                }
            } else {
                if (!cacheFront())
                    _rem = 0;
            }
        }

        void popBack() {
            if (_md.empty) {
                if (_rem == 255) {
                    _cached_back = _cached_front;
                    _rem = 1;
                } else {
                    _rem = 0;
                }
            } else {
                if (!cacheBack())
                    _rem = 0;
            }
        }

        private bool cacheFront() {
            if (_md.empty)
                return false;

            if (_md[0] == '^') {          // deletion, get bases
                _md = _md[1 .. $];
                auto len = countUntil!(not!isUpper)(_md);
                if (len == -1) {
                    len = _md.length;
                }
                _cached_front = MdOperation.createDeletion(_md[0 .. len]);
                _md = _md[len .. $];
            } else if (isDigit(_md[0])) { // match, get number
                auto len = countUntil!(not!isDigit)(_md);
                if (len == -1) {
                    len = _md.length;
                }
                _cached_front = MdOperation.createMatch(to!uint(_md[0 .. len]));
                _md = _md[len .. $];
            } else {                     // mismatch
                _cached_front = MdOperation.createMismatch(_md[0]);
                _md = _md[1 .. $];
            }

            return true;
        }

        private bool cacheBack() {
            if (_md.empty)
                return false;

            if (isDigit(_md[$ - 1])) { // match, get number
                auto len = countUntil!(not!isDigit)(retro(_md));
                if (len == -1) {
                    len = _md.length;
                }
                _cached_back = MdOperation.createMatch(to!uint(_md[$ - len .. $]));
                _md = _md[0 .. $ - len];
            } else {
                if (_md.length == 1 || isDigit(_md[$ - 2])) { // mismatch
                    _cached_back = MdOperation.createMismatch(_md[$ - 1]);
                    _md = _md[0 .. $ - 1];
                } else { // deletion
                    auto len = countUntil!"a == '^'"(retro(_md));
                    _cached_back = MdOperation.createDeletion(_md[$ - len .. $]);
                    _md = _md[0 .. $ - len - 1];
                }
            }

            return true;
        }
    }

    static bool isZeroMatch(MdOperation op) {
        return op.type == MdOperationType.Match &&
               op.match == 0;
    }

    return filterBidirectional!(not!isZeroMatch)(Result(md));
}

/// Alias for return type of mdOperations
alias ReturnType!mdOperations MdOperationRange;

unittest {

    import std.algorithm;

    import std.stdio;
    
    assert(equal(mdOperations("86"), 
                [MdOperation.createMatch(86)]));

    assert(equal(mdOperations("0G81"), 
                [MdOperation.createMismatch('G'), 
                 MdOperation.createMatch(81)]));

    assert(equal(mdOperations("62^T28"), 
                [MdOperation.createMatch(62), 
                 MdOperation.createDeletion("T"), 
                 MdOperation.createMatch(28)]));

    assert(equal(retro(mdOperations("3C6C0A13^A4C2")),
                 retro([MdOperation.createMatch(3),   
                        MdOperation.createMismatch('C'), 
                        MdOperation.createMatch(6),
                        MdOperation.createMismatch('C'), 
                        MdOperation.createMismatch('A'),
                        MdOperation.createMatch(13),  
                        MdOperation.createDeletion("A"), 
                        MdOperation.createMatch(4),
                        MdOperation.createMismatch('C'), 
                        MdOperation.createMatch(2)])));

    assert(equal(mdOperations("27^TTT63"),
                [MdOperation.createMatch(27), 
                 MdOperation.createDeletion("TTT"), 
                 MdOperation.createMatch(63)]));
}
