/*
    This file is part of BioD.
    Copyright (C) 2012    Artem Tarasov <lomereiter@gmail.com>

    Permission is hereby granted, free of charge, to any person obtaining a
    copy of this software and associated documentation files (the "Software"),
    to deal in the Software without restriction, including without limitation
    the rights to use, copy, modify, merge, publish, distribute, sublicense,
    and/or sell copies of the Software, and to permit persons to whom the
    Software is furnished to do so, subject to the following conditions:
    
    The above copyright notice and this permission notice shall be included in
    all copies or substantial portions of the Software.
    
    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
    FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
    DEALINGS IN THE SOFTWARE.

*/
module bio.std.hts.bam.baseinfo;

import bio.core.base;
import bio.core.sequence;

import bio.std.hts.bam.read;
import bio.std.hts.bam.tagvalue;
import bio.std.hts.iontorrent.flowcall;
import bio.std.hts.bam.md.core;
import bio.std.hts.bam.cigar;

import std.range;
import std.conv;
import std.traits;
import std.typecons;
import std.typetuple;

/// 
enum Option
{
    /// adds 'cigar_before' and 'cigar_after' properties
    cigarExtra, 

    /// adds 'md_operation', 'md_operation_offset' properties
    mdCurrentOp,

    /// adds 'previous_md_operation' property
    mdPreviousOp,

    /// adds 'next_md_operation' property
    mdNextOp
}

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

template staticFilter(alias P, T...)
{
    static if (T.length == 0)
        alias TypeTuple!() staticFilter;
    else static if (P!(T[0]))
        alias TypeTuple!(T[0], staticFilter!(P, T[1..$])) staticFilter;
    else
        alias staticFilter!(P, T[1..$]) staticFilter;
}

template isTag(alias argument)
{
    enum isTag = is(typeof(argument) == string);
}

template isOption(alias argument)
{
    enum isOption = is(typeof(argument) == Option);
}

struct PerBaseInfo(R, TagsAndOptions...) {

    alias staticFilter!(isTag, TagsAndOptions) Tags;
    alias staticFilter!(isOption, TagsAndOptions) Options;

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
    //          Current base of the result is updated before the call.
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
            result ~= "mixin " ~ ext ~ "baseInfo!(R, Options).resultProperties;".dup;
        return cast(string)result;
    }

    static struct Result {
        /// Actual read base, with strand taken into account.
        Base base;
        alias base this;

        string opCast(T)() if (is(T == string))
        {
            return to!string(base);
        }

        bool opEquals(T)(T base) const
            if (is(Unqual!T == Base)) 
        {
            return this.base == base;
        }

        bool opEquals(T)(T result) const
            if (is(Unqual!T == Result))
        {
            return this == result;
        }

        bool opEquals(T)(T base) const
            if (is(Unqual!T == char) || is(Unqual!T == dchar))
        {
            return this.base == base;
        }

        mixin(getResultProperties!Extensions());
    }

    private static string getRangeMethods(Exts...)() {
        char[] result;
        foreach (ext; Exts)
            result ~= "mixin " ~ ext ~ "baseInfo!(R, Options).rangeMethods " ~ ext ~ ";".dup;
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
        _rev = read.is_reverse_strand;
        _seq = reversableRange!complementBase(read.sequence, _rev);

        foreach (t; Extensions) {
            setup!t(read, args);
        }
    }

    bool empty() @property {
        return _seq.empty;
    }

    /// Allows to construct front element in-place, avoiding a copy.
    void constructFront(Result* addr)
    {
        addr.base = _seq.front;
        foreach (t; Extensions)
            populate!t(*addr);
    }

    Result front() @property {
        Result r = void;
        r.base = _seq.front;
        foreach (t; Extensions)
            populate!t(r);
        return r;
    }

    void popFront() {
        moveToNextBase();
    }

    PerBaseInfo save() @property {
        PerBaseInfo r = void;
        r._read = _read.dup;
        r._seq = _seq.save;
        r._rev = _rev;
        foreach (t; Extensions)
            copy!t(r);
        return r;
    }

    ref PerBaseInfo opAssign(PerBaseInfo other) {
        _read = other._read;
        _seq = other._seq.save;
        _rev = other._rev;
        foreach (t; Extensions)
            other.copy!t(this);
        return this;
    }

    private void moveToNextBase() {

        foreach (t; Extensions) {
            update!t();
        }

        _seq.popFront();
    }

    /// Returns true if the read is reverse strand,
    /// and false otherwise.
    bool reverse_strand() @property const {
        return _rev;
    }

    private {
        bool _rev = void;
        R _read = void;
        ReversableRange!(complementBase, typeof(_read.sequence)) _seq = void;
    }
}

///
///  Collect per-base information from available tags. 
///  Use $(D arg!TagName) to pass a parameter related to a particular tag.
///
///  Example:
///
/// basesWith!"FZ"(arg!"flowOrder"(flow_order), arg!"keySequence"(key_sequence));
///
template basesWith(TagsAndOptions...) {
    auto basesWith(R, Args...)(R read, Args args) {
        return PerBaseInfo!(R, TagsAndOptions)(read, args);
    }
}

/// Provides additional property $(D reference_base)
template MDbaseInfo(R, Options...) {

    mixin template resultProperties() {
    
        enum MdCurrentOp = staticIndexOf!(Option.mdCurrentOp, Options) != -1;
        enum MdPreviousOp = staticIndexOf!(Option.mdPreviousOp, Options) != -1;
        enum MdNextOp = staticIndexOf!(Option.mdNextOp, Options) != -1;

        /// If current CIGAR operation is reference consuming,
        /// returns reference base at this position, otherwise
        /// returns '-'.
        ///
        /// If read is on '-' strand, the result will be
        /// complementary base.
        char reference_base() @property const {
            return _ref_base;
        }

        private char _ref_base = void;

        static if (MdPreviousOp)
        {
            private Nullable!MdOperation _previous_md_operation = void;

            /// Previous MD operation
            Nullable!MdOperation previous_md_operation() @property {
                return _previous_md_operation;
            }
        }

        static if (MdCurrentOp)
        {

            private MdOperation _current_md_operation = void;
            private uint _current_md_operation_offset = void;

            /// Current MD operation
            MdOperation md_operation() @property {
                return _current_md_operation;
            }

            /// If current MD operation is match, returns how many bases
            /// have matched before the current base. Otherwise returns 0.
            uint md_operation_offset() @property const {
                return _current_md_operation_offset;
            }
        }

        static if (MdNextOp)
        {
            private Nullable!MdOperation _next_md_operation = void;
            /// Next MD operation
            Nullable!MdOperation next_md_operation() @property {
                return _next_md_operation;
            }
        }
    }

    mixin template rangeMethods() {

        enum MdCurrentOp = staticIndexOf!(Option.mdCurrentOp, Options) != -1;
        enum MdPreviousOp = staticIndexOf!(Option.mdPreviousOp, Options) != -1;
        enum MdNextOp = staticIndexOf!(Option.mdNextOp, Options) != -1;

        private {
            ReversableRange!(reverseMdOp, MdOperationRange) _md_ops = void;
            uint _match; // remaining length of current match operation
            MdOperation _md_front = void;

            static if (MdPreviousOp)
            {
                Nullable!MdOperation _previous_md_op;
                bool _md_front_is_initialized;
            }
        }

        private void updateMdFrontVariable()
        {
            static if (MdPreviousOp)
            {
                if (_md_front_is_initialized)
                    _previous_md_op = _md_front;

                _md_front_is_initialized = true;
            }

            _md_front = _md_ops.front;
            _md_ops.popFront();
        }

        void setup(Args...)(const ref R read, Args args)
        {
            auto md = read["MD"];
            auto md_str = *(cast(string*)&md);
            _md_ops = reversableRange!reverseMdOp(mdOperations(md_str),
                                                  read.is_reverse_strand);
         
            while (!_md_ops.empty)
            {
                updateMdFrontVariable();
                if (!_md_front.is_deletion) {
                    if (_md_front.is_match) {
                        _match = _md_front.match;
                    }
                    break;
                }
            }
        }

        void populate(Result)(ref Result result)
        {
            if (!current_cigar_operation.is_reference_consuming)
            {
                result._ref_base = '-';
                return;
            }

            MdOperation op = _md_front;
            if (op.is_mismatch)
                result._ref_base = op.mismatch.asCharacter;
            else if (op.is_match) {
                result._ref_base = result.base.asCharacter;
            }
            else assert(0);

            static if (MdPreviousOp)
            {
                if (_previous_md_op.isNull)
                    result._previous_md_operation.nullify();
                else
                    result._previous_md_operation = _previous_md_op.get;
            }

            static if (MdCurrentOp)
            {

                result._current_md_operation = op;
                result._current_md_operation_offset = _md_front.match - _match;
            }

            static if (MdNextOp)
            {
                if (_md_ops.empty)
                    result._next_md_operation.nullify();
                else
                    result._next_md_operation = _md_ops.front;
            }
        }

        void update(const ref R read)
        {
            if (!current_cigar_operation.is_reference_consuming)
                return;

            if (_md_front.is_mismatch)
            {
                if (_md_ops.empty)
                    return;

                updateMdFrontVariable();
            }
            else if (_md_front.is_match)
            {
                --_match;
                if (_match == 0 && !_md_ops.empty) {
                    updateMdFrontVariable();
                }
            }
            else assert(0);

            while (_md_front.is_deletion) {
                if (_md_ops.empty)
                    return;

                updateMdFrontVariable();
            }

            if (_match == 0 && _md_front.is_match)
                _match = _md_front.match;
        }

        void copy(Range)(ref Range source, ref Range target)
        {
            target.MD._md_ops = source.MD._md_ops.save;
            target.MD._md_front = source.MD._md_front;

            static if (MdPreviousOp)
            {
                if (source.MD._previous_md_op.isNull)
                    target.MD._previous_md_op.nullify();
                else
                    target.MD._previous_md_op = source.MD._previous_md_op.get;
                target.MD._md_front_is_initialized = source.MD._md_front_is_initialized;
            }
        }
    }
}

/// Provides additional property $(D flow_call).
template FZbaseInfo(R, Options...) {

    mixin template resultProperties() {
        /// Current flow call
        ReadFlowCall flow_call() @property const {
            return _flow_call;
        }

        private {
            ReadFlowCall _flow_call;
        }
    }

    mixin template rangeMethods() {

        private {
            ReadFlowCallRange!(BamRead.SequenceResult) _flow_calls = void;
            ReadFlowCall _current_flow_call = void;
            ushort _at = void;

            debug {
                string _read_name;
            }
        }

        void setup(Args...)(const ref R read, Args args) 
        {
            string flow_order = void;
            string key_sequence = void;

            debug {
                _read_name = read.name.idup;
            }

            enum flowOrderExists = staticIndexOf!(MixinArg!(string, "flowOrder"), Args);
            enum keySequenceExists = staticIndexOf!(MixinArg!(string, "keySequence"), Args);
            static assert(flowOrderExists != -1, `Flow order must be provided via arg!"flowOrder"`);
            static assert(keySequenceExists != -1, `Flow order must be provided via arg!"keySequence"`);

            foreach (arg; args) {
                static if(is(typeof(arg) == MixinArg!(string, "flowOrder")))
                    flow_order = arg;

                static if(is(typeof(arg) == MixinArg!(string, "keySequence")))
                    key_sequence = arg;
            }

            _at = 0;

            _flow_calls = readFlowCalls(read, flow_order, key_sequence);
            if (!_flow_calls.empty) {
                _current_flow_call = _flow_calls.front;
            }
        }

        void populate(Result)(ref Result result) {
            result._flow_call = _current_flow_call;

            debug {
                if (result.base != result._flow_call.base) {
                    import std.stdio;
                    stderr.writeln("invalid flow call at ", _read_name, ": ", result.position);
                }
            }
        }

        void update(const ref R read) 
        {
            ++_at;
            if (_at == _current_flow_call.length) {
                _flow_calls.popFront();
                if (!_flow_calls.empty) {
                    _current_flow_call = _flow_calls.front;
                    _at = 0;
                }
            }
        }

        void copy(Range)(ref Range source, ref Range target) {
            target.FZ._flow_calls = source._flow_calls.save();
            target.FZ._at = source.FZ._at;
            target.FZ._current_flow_call = source._current_flow_call;

            debug {
                target._read_name = _read_name;
            }
        }
    }
}

/// Retrieving flow signal intensities from ZM tags is also available.
alias FZbaseInfo ZMbaseInfo;

/// Provides additional properties
///     * position
///     * cigar_operation
///     * cigar_operation_offset
template CIGARbaseInfo(R, Options...) {

    mixin template resultProperties() {

        enum CigarExtraProperties = staticIndexOf!(Option.cigarExtra, Options) != -1;

        static if (CigarExtraProperties)
        {
            /// Current CIGAR operation
            CigarOperation cigar_operation() @property {
                return _cigar[_operation_index];
            }

            /// CIGAR operations before current one
            auto cigar_before() @property {
                return _cigar[0 .. _operation_index];
            }

            /// CIGAR operations after current one
            auto cigar_after() @property {
                return _cigar[_operation_index + 1 .. _cigar.length];
            }
        }
        else
        {
            /// Current CIGAR operation
            CigarOperation cigar_operation() @property const {
                return _current_cigar_op;
            }
        }

        /// Position of the corresponding base on the reference.
        /// If current CIGAR operation is not one of 'M', '=', 'X',
        /// returns the position of the previous mapped base.
        uint position() @property const {
            return _reference_position;
        }

        /// Offset in current CIGAR operation, starting from 0.
        uint cigar_operation_offset() @property const {
            return _cigar_operation_offset;
        }

        private {
            int _operation_index = void;
            uint _reference_position = void;
            uint _cigar_operation_offset = void;
            static if (CigarExtraProperties)
            {
                ReversableRange!(identity, const(CigarOperation)[]) _cigar = void;
            }
            else
            {
                CigarOperation _current_cigar_op;
            }
        }
    }

    mixin template rangeMethods() {

        enum CigarExtraProperties = staticIndexOf!(Option.cigarExtra, Options) != -1;

        private {
            CigarOperation _current_cigar_op = void;

            ulong _cur_cig_op_len = void;
            bool _cur_cig_op_is_ref_cons = void;

            int _index = void;
            uint _at = void;
            uint _ref_pos = void;
            ReversableRange!(identity, const(CigarOperation)[]) _cigar = void;
        }

        /// Current CIGAR operation, available to all extensions
        const(CigarOperation) current_cigar_operation() @property const {
            return _current_cigar_op;
        }

        void setup(Args...)(const ref R read, Args) 
        {
            _cigar = reversableRange(read.cigar, read.is_reverse_strand);

            _index = -1;
            _ref_pos = reverse_strand ? (read.position + read.basesCovered() - 1)
                                      : read.position;

            _moveToNextCigarOperator();
            assert(_index >= 0);
        }

        void populate(Result)(ref Result result) {
            result._reference_position = _ref_pos;
            result._cigar_operation_offset = _at;
            static if (CigarExtraProperties)
            {
                result._cigar = _cigar;
                result._operation_index = _index;
            }
            else
            {
                result._current_cigar_op = _current_cigar_op;
            }
        }

        void update(const ref R read) 
        {
           ++_at;

           if  (_cur_cig_op_is_ref_cons) {
               _ref_pos += reverse_strand ? -1 : 1;
           }

           if (_at == _cur_cig_op_len) {
               _moveToNextCigarOperator();
           }
        }

        void copy(Range)(const ref Range source, ref Range target) {
            target.CIGAR._cigar = source.CIGAR._cigar;
            target.CIGAR._index = source.CIGAR._index;
            target.CIGAR._current_cigar_op = source.CIGAR._current_cigar_op;
            target.CIGAR._cur_cig_op_len = source.CIGAR._cur_cig_op_len;
            target.CIGAR._cur_cig_op_is_ref_cons = source.CIGAR._cur_cig_op_is_ref_cons;
            target.CIGAR._at = source.CIGAR._at;
            target.CIGAR._ref_pos = source.CIGAR._ref_pos;
        }

        private void _moveToNextCigarOperator() {
            _at = 0;
            for (++_index; _index < _cigar.length; ++_index)
            {
                _current_cigar_op = _cigar[_index];
                _cur_cig_op_is_ref_cons = _current_cigar_op.is_reference_consuming;
                _cur_cig_op_len = _current_cigar_op.length;

                if (_current_cigar_op.is_query_consuming)
                    break;

                if (_cur_cig_op_is_ref_cons)
                {
                    if (reverse_strand)
                        _ref_pos -= _cur_cig_op_len;
                    else
                        _ref_pos += _cur_cig_op_len;
                }
            }
        }
    }
}
