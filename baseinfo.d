module baseinfo;

import BioD.Base;

import alignment;
import fz.flowcall;

import std.range;
import std.conv;
import std.traits;
import std.typetuple;

///
struct MixinArg(T, string Tag) {
    T value;
    alias value this;
    alias Tag TagName;
}

/// Wrapper for arguments to $(D basesWith) function (see below).
/// Required to distinguish to which tag each parameter refers.
MixinArg!(T, Tag) arg(string Tag, T)(T value) {
    return MixinArg!(T, Tag)(value);
}

struct PerBaseInfo(R, Tags...) {

    private alias TypeTuple!("CIGAR", Tags) Extensions;

    // /////////////////////////////////////////////////////////////////////////
    //
    // Each 'extension' is a template with name TAGbaseInfo, containing 
    // a couple of mixin templates:
    // 
    // * resultProperties
    //      These are additional properties provided by the template
    //
    // * rangeMethods
    //      These describe how to proceed to the next base.
    //      The following methods must be implemented:
    //      
    //      - void setup(Args...)(const ref R read, Args args);
    //          Gets called during range construction. All constructor
    //          arguments are forwarded, and it's this function which
    //          is responsible for getting required parameters for this
    //          particular template.
    //
    //      - void populate(Result)(ref Result r);
    //          Populates fields of the result declared in resultProperties.
    //          Should run in O(1), just copying a few variables.
    //
    //      - void update(const ref R read);
    //          Encapsulates logic of moving to the next base and updating
    //          mixin variables correspondingly.
    //
    //      - void copy(Range)(const ref Range source, ref Range target);
    //          Gets called during $(D source.save). Therefore, any ranges
    //          used in mixin templates must be saved as well at that time.
    //
    // /////////////////////////////////////////////////////////////////////////

    private static string getResultProperties(Exts...)() {
        char[] result;
        foreach (ext; Exts) 
            result ~= "mixin " ~ ext ~ "baseInfo!R.resultProperties;".dup;
        return cast(string)result;
    }

    static struct Result {
        Base base;
        alias base this;

        string opCast(T)() if (is(T == string))
        {
            return to!string(base);
        }

        mixin(getResultProperties!Extensions());
    }

    private static string getRangeMethods(Exts...)() {
        char[] result;
        foreach (ext; Exts)
            result ~= "mixin " ~ ext ~ "baseInfo!R.rangeMethods " ~ ext ~ ";".dup;
        return cast(string)result;
    }

    mixin(getRangeMethods!Extensions());

    private void setup(string tag, Args...)(R read, Args args) {
        mixin(tag ~ ".setup(read, args);");
    }

    private void populate(string tag)(ref Result r) {
        mixin(tag ~ ".populate(r);");
    }

    private void update(string tag)() {
        mixin(tag ~ ".update(_read);");
    }

    private void copy(string tag)(ref typeof(this) other) {
        mixin(tag ~ ".copy(this, other);");
    }

    this(Args...)(R read, Args args) {
        _read = read;
        _seq = read.sequence;
        _rev = read.is_reverse_strand;

        foreach (t; Extensions) {
            setup!t(read, args);
        }
    }

    bool empty() @property const {
        return _seq.empty;
    }

    Result front() @property {
        Result r = void;
        r.base = _rev ? _seq.back : _seq.front;
        foreach (t; Extensions)
            populate!t(r);
        return r;
    }

    void popFront() {
        moveToNextBase();
    }

    PerBaseInfo!(R, Tags) save() @property {
        PerBaseInfo!(R, Tags) r = void;
        r._read = _read.dup;
        r._seq = r._read.sequence;
        r._rev = _rev;
        foreach (t; Extensions)
            copy!t(r);
        return r;
    }

    private void moveToNextBase() {

        foreach (t; Extensions) {
            update!t();
        }

        if (_rev) 
            _seq.popBack();
        else 
            _seq.popFront();
    }

    private {
        R _read;
        typeof(_read.sequence) _seq;
        bool _rev;
    }
}

///
///  Collect per-base information from available tags. 
///  Use $(D arg!TagName) to pass a parameter related to a particular tag.
///
///  Example:
///
/// basesWith!"FZ"(arg!"FZ"(flow_order));
///
template basesWith(Tags...) {
    auto basesWith(R, Args...)(R read, Args args) {
        return PerBaseInfo!(R, Tags)(read, args);
    }
}

/// Provides additional property $(D flow_call).
template FZbaseInfo(R) {

    mixin template resultProperties() {
        /// Current flow call
        ReadFlowCall flow_call() @property {
            return _flow_call;
        }

        private {
            ReadFlowCall _flow_call;
        }
    }

    mixin template rangeMethods() {

        private {
            ForwardRange!ReadFlowCall _flow_calls;
            ushort _at;

            debug {
                string _read_name;
            }
        }

        void setup(Args...)(const ref R read, Args args) 
        {
            string flow_order;

            debug {
                _read_name = read.read_name.idup;
            }

            enum argExists = staticIndexOf!(MixinArg!(string, "FZ"), Args);
            static assert(argExists != -1, `Flow order must be provided via arg!"FZ"`);

            foreach (arg; args) {
                static if(is(typeof(arg) == MixinArg!(string, "FZ")))
                    flow_order = arg;
            }

            _flow_calls = readFlowCalls(read, flow_order);
            _at = 0;
        }

        void populate(Result)(ref Result result) {
            result._flow_call = _flow_calls.front;

            debug {
                if ((_rev && result.base != result._flow_call.base.complement)
                    || (!_rev && result.base != result._flow_call.base)) {
                    import std.stdio;
                    stderr.writeln("invalid flow call at ", _read_name, ": ", result.position);
                }
            }
        }

        void update(const ref R read) 
        {
            ++_at;
            if (_at == _flow_calls.front.length) {
                _flow_calls.popFront();
                _at = 0;
            }
        }

        void copy(Range)(ref Range source, ref Range target) {
            target.FZ._flow_calls = source._flow_calls.save();
            target.FZ._at = source.FZ._at;

            debug {
                target._read_name = _read_name;
            }
        }
    }
}

/// Provides additional properties
///     * position
///     * cigar_operation
template CIGARbaseInfo(R) {

    mixin template resultProperties() {
        /// Current CIGAR operation
        CigarOperation cigar_operation() @property {
            return _cigar_operation;
        }

        /// Position of the corresponding base on the reference.
        /// If current CIGAR operation is not one of 'M', '=', 'X',
        /// returns the position of the previous valid base.
        ulong position() @property {
            return _reference_position;
        }

        private {
            CigarOperation _cigar_operation;
            ulong _reference_position;
        }
    }

    mixin template rangeMethods() {

        private {
            const(CigarOperation)[] _cigar;
            long _index;
            ulong _at;
            ulong _ref_pos;
        }

        void setup(Args...)(const ref R read, Args) 
        {
            _cigar = read.cigar;

            _index = read.is_reverse_strand ? _cigar.length : -1;
            _ref_pos = read.is_reverse_strand ? (read.position + read.basesCovered() - 1)
                                              : read.position;

            _moveToNextCigarOperator(read.is_reverse_strand);
        }

        void populate(Result)(ref Result result) {
            result._cigar_operation = _cigar[_index];
            result._reference_position = _ref_pos;
        }

        void update(const ref R read) 
        {
           ++_at;
           if (_cigar[_index].is_reference_consuming) {
               _ref_pos += read.is_reverse_strand ? -1 : 1;
           }

           if (_at == _cigar[_index].length) {
               _moveToNextCigarOperator(read.is_reverse_strand);
           }
        }

        void copy(Range)(const ref Range source, ref Range target) {
            target.CIGAR._cigar = source.CIGAR._cigar;
            target.CIGAR._index = source.CIGAR._index;
            target.CIGAR._at = source.CIGAR._at;
            target.CIGAR._ref_pos = source.CIGAR._ref_pos;
        }

        private void _moveToNextCigarOperator(bool reverse) {
            _at = 0;
            if (!reverse) {
                for (++_index; _index < _cigar.length; ++_index)
                {
                    if (_cigar[_index].is_query_consuming)
                        break;
                    if (_cigar[_index].is_reference_consuming)
                        _ref_pos += _cigar[_index].length;
                }
            } else {
                for (--_index; _index >= 0; --_index) {
                    if (_cigar[_index].is_query_consuming)
                        break;
                    if (_cigar[_index].is_reference_consuming)
                        _ref_pos -= _cigar[_index].length;
                }
            }
        }
    }
}
