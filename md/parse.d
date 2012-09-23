module md.parse;

import md.operation;
import std.ascii;
import std.array;
import std.algorithm;
import std.functional;
import std.conv;

/// Returns range of MD operations. Zero matches are skipped.
auto mdOperations(string md) {

    struct Result {
        private {
            MdOperation _cur_md_op; // current MD operation
            bool _empty = false;
            string _md;
        }

        this(string md) {
            _md = md;
            popFront();
        }
      
        bool empty() @property {
            return _empty;
        }

        MdOperation front() @property {
            return _cur_md_op;
        }
        
        void popFront() {
            if (_md.empty) {
                _empty = true;
                return;
            }

            if (_md[0] == '^') {
                _md = _md[1 .. $];
                auto len = countUntil!(not!isUpper)(_md);
                if (len == -1) {
                    len = _md.length;
                }
                _cur_md_op = MdOperation.createDeletion(_md[0 .. len]);
                _md = _md[len .. $];
            } else if (isDigit(_md[0])) {
                auto len = countUntil!(not!isDigit)(_md);
                if (len == -1) {
                    len = _md.length;
                }
                _cur_md_op = MdOperation.createMatch(to!uint(_md[0 .. len]));
                _md = _md[len .. $];
            } else {
                _cur_md_op = MdOperation.createMismatch(_md[0]);
                _md = _md[1 .. $];
            }
        }
    }

    static bool isZeroMatch(MdOperation op) {
        return op.type == MdOperationType.Match &&
               op.match == 0;
    }

    return filter!(not!isZeroMatch)(Result(md));
}

unittest {

    import std.algorithm;

    assert(equal(mdOperations("86"), 
                [MdOperation.createMatch(86)]));

    assert(equal(mdOperations("0G81"), 
                [MdOperation.createMismatch('G'), 
                 MdOperation.createMatch(81)]));

    assert(equal(mdOperations("62^T28"), 
                [MdOperation.createMatch(62), 
                 MdOperation.createDeletion("T"), 
                 MdOperation.createMatch(28)]));

    assert(equal(mdOperations("3C6C0A13^A4C2"),
                [MdOperation.createMatch(3),   
                 MdOperation.createMismatch('C'), 
                 MdOperation.createMatch(6),
                 MdOperation.createMismatch('C'), 
                 MdOperation.createMismatch('A'),
                 MdOperation.createMatch(13),  
                 MdOperation.createDeletion("A"), 
                 MdOperation.createMatch(4),
                 MdOperation.createMismatch('C'), 
                 MdOperation.createMatch(2)]));

    assert(equal(mdOperations("27^TTT63"),
                [MdOperation.createMatch(27), 
                 MdOperation.createDeletion("TTT"), 
                 MdOperation.createMatch(63)]));
}
