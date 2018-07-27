module cram.wrappers;

import cram.htslib, cram.exception;

import std.typecons, std.parallelism, std.range;
debug import std.stdio;

auto zeroChecked(alias func, T...)(string err_msg, auto ref T params) {
    int ret = func(params);
    if (!ret)
        throw new CramException(err_msg);
}

struct RcPtr(T, alias Free) {
    static struct Payload {
        debug {
            static size_t payload_counter;
            size_t payload_id;
        }

        T* ptr;
        this(T* ptr) {
            this.ptr = ptr;
            debug {
                payload_id = ++payload_counter;
                // stderr.writeln("Init ", T.stringof, "* #", payload_id);
            }
        }

        ~this() {
            debug {
                // stderr.writeln("Free ", T.stringof, "* #", payload_id);
            }
            Free(ptr);
        }

        alias ptr this;

        this(this) { assert(false); }
        void opAssign(Payload rhs) { assert(false); }
    }

    alias Data = RefCounted!(Payload, RefCountedAutoInitialize.no);
    Data data;
    alias data this;

    bool isNull() @property const {
        return !data.refCountedStore.isInitialized || data.ptr is null;
    }

    this(this)
    {
        static if (is(T == cram_slice)) {
            // debug stderr.writeln("COPIED #", data.payload_id + 1);
        }
    }

    this(T* ptr) { data = Data(ptr); }
}

auto nullChecked(alias func, T...)(string err_msg, cram_fd* fd,
                                   auto ref T other_params)
{
    auto ptr = func(fd, other_params);
    if (ptr !is null)
        return ptr;
    if (fd.err == 0)
        return null;
    throw new CramException(err_msg);
}

// instances of these types must be used within the same thread!
void my_cram_close(cram_fd* fd) {
    import core.sys.posix.pthread;
    pthread_mutex_destroy(&fd.metrics_lock);
    pthread_mutex_destroy(&fd.ref_lock);
    pthread_mutex_destroy(&fd.bam_list_lock);
    cram_close(fd);
}

alias RcCramFd = RcPtr!(cram_fd, my_cram_close);
alias CramFd = RcCramFd;

CramFd openCram(string filename) {
    import std.string : toStringz;
    auto fd = cram_open(toStringz(filename), "rb");

    if (fd == null)
        throw new CramException("Can't open file " ~ filename);

    cram_set_option(fd, cram_option.CRAM_OPT_DECODE_MD);

    // initialize locks, but we will use the pool from D standard library
    // instead of the htslib implementation
    import core.sys.posix.pthread;
    pthread_mutex_init(&fd.metrics_lock, null);
    pthread_mutex_init(&fd.ref_lock, null);
    pthread_mutex_init(&fd.bam_list_lock, null);
    fd.shared_ref = 1;
    fd.own_pool = 0;
    fd.pool = null;

    return CramFd(fd);
}

enum CramFilterResult {
    skip, // skip over the object
    pass, // the filter doesn't reject the object
    stop  // skip, and further iteration won't give any results
}

alias RcCramContainer = RcPtr!(cram_container, cram_free_container);
alias CramContainer = RcCramContainer;
alias CramContainerFilter = CramFilterResult delegate(cram_container*);

alias RcCramSlice = RcPtr!(cram_slice, cram_free_slice);
alias UndecodedSliceFilter = CramFilterResult delegate(cram_slice*);

struct CramSlice {
    CramFd fd;
    RcCramContainer container;
    RcCramSlice slice;

    alias slice this;

    bool is_decoded() @property const {
        return slice.crecs !is null;
    }
}

struct CramContainerRange {
    CramContainer front;
    bool empty;

    private CramFd _fd;
    private CramContainerFilter _filter;

    // the containers that doesn't pass the filter
    // will be simply skipped without any decoding
    this(CramFd fd, CramContainerFilter f)
    {
        _fd = fd;
        _filter = f;
        popFront();
    }

    void popFront() {
        if (!front.isNull) {
            // skip remaining blocks of previous container
            while (front.curr_slice != front.max_slice) {
                ++front.curr_slice;
                auto next_slice = cram_read_slice(_fd);
                if (!next_slice)
                    throw new CramException("Failure in cram_read_slice");
                cram_free_slice(next_slice);
            }
        }

        auto err_msg = "Failed to read container header";
        while (true) {
            // read container header
            // debug stderr.writeln("cram_read_container");
            auto ptr = nullChecked!cram_read_container(err_msg, _fd);
            if (ptr is null) {
                empty = true;
                break;
            }
            front = CramContainer(ptr);
            // apply the filter
            if (_filter is null) break;
            auto res = _filter(front);

            if (res == CramFilterResult.stop) {
                empty = true;
                break;
            }

            if (res == CramFilterResult.pass) {
                break;
            }

            if (-1 == cram_seek(_fd, front.length, 1))
                throw new CramException("Failed to seek in CRAM file");
        }

        if (empty) return;

        // read compression header block
        err_msg = "Failed to read compression header block";
        front.comp_hdr_block = nullChecked!cram_read_block(err_msg, _fd);
        auto content_type = front.comp_hdr_block.content_type;
        if (content_type != cram_content_type.COMPRESSION_HEADER)
            throw new CramException(err_msg);

        err_msg = "Failed to decode compression header";
        front.comp_hdr = cram_decode_compression_header(_fd,
                                                        front.comp_hdr_block);
        if (front.comp_hdr is null)
            throw new CramException(err_msg);

        if (!(front.comp_hdr.AP_delta)) {
            // FIXME: in cram_decode.c, there's a mutex locked around this line,
            //        but it looks unnecessary
            _fd.unsorted = 1;
        }
    }
}

auto containers(CramFd fd, CramContainerFilter f) {
    return CramContainerRange(fd, f);
}

class UndecodedSliceRange {
    private {
        CramFd _fd;
        ref CramContainer _container() @property {
            assert(!_containers.empty);
            return _containers.front;
        }
        CramContainerRange _containers;
        UndecodedSliceFilter _sf;
    }

    this(CramFd fd, CramContainerFilter cf, UndecodedSliceFilter sf) {
        _fd = fd;
        _sf = sf;
        _containers = containers(_fd, cf);

        if (_containers.empty)
            empty = true;
        popFront();
    }

    CramSlice front;
    bool empty;

    // true == either front or empty is set
    private bool readNextSliceFromCurrentContainer() {
        assert(_container.curr_slice < _container.max_slice);
        _container.curr_slice++;

        auto err_msg = "Failure in cram_read_slice";
        // debug stderr.writeln("cram_read_slice (", _container.curr_slice,
        //                      "/", _container.max_slice, ")");
        auto ptr = cram_read_slice(_fd);
        if (ptr is null) {
            throw new CramException(err_msg);
        }
        assert(ptr.hdr !is null);
        ptr.last_apos = ptr.hdr.ref_seq_start;

        if (_sf is null) { setupFront(ptr); return true; }
        auto ret = _sf(ptr);
        if (ret == CramFilterResult.pass) { setupFront(ptr); return true; }
        if (ret == CramFilterResult.stop) { empty = true; return true; }
        return false;
    }

    private void setupFront(cram_slice* ptr) {
        front = CramSlice(_fd, _container, RcCramSlice(ptr));
    }

    void popFront() {
        while (!empty) {
            while (_container.curr_slice < _container.max_slice)
                if (readNextSliceFromCurrentContainer())
                    return;

            // no slice (satisfying the filter) in the current container
            _containers.popFront();
            if (_containers.empty) { empty = true; return; }
        }
    }
}

auto undecodedSlices(CramFd fd, CramContainerFilter cf,
                     UndecodedSliceFilter sf=null)
{
    return new UndecodedSliceRange(fd, cf, sf);
}

void decodeSlice(cram_fd* fd, cram_container* c, cram_slice* s) {
    auto err_msg = "Failure in cram_decode_slice";
    // debug writeln("DECODING slice #", s.id + 1);
    int ret = cram_decode_slice(fd, c, s, fd.header);
    if (ret != 0)
        throw new CramException(err_msg);
}

void decodeSlice(CramSlice slice) {
    decodeSlice(slice.fd, slice.container, slice.data.ptr);
}

import bio.core.utils.roundbuf;

struct CramSliceDecoder(R)
    if (isInputRange!R && is(ElementType!R == CramSlice))
{
    private {
        R _slices;
        TaskPool _pool;

        // FIXME: D arrays don't call element destructors when GC-d :(
        RoundBuf!CramSlice _input_queue;
        alias DecodeTask = Task!(decodeSlice,
                                 cram_fd*, cram_container*, cram_slice*)*;
        RoundBuf!DecodeTask _output_queue;

        void putNextSliceIntoQueue() {
            auto slice = _slices.front;
            _input_queue.put(slice);
            _slices.popFront();

            // don't copy it between threads (ref. counting is not atomic)
            cram_fd* fd = slice.fd;
            cram_container* c = slice.container;
            cram_slice* s = slice;
            // debug writeln("PUT slice #", s.id + 1, " into queue");
            version (serial) {
                decodeSlice(fd, c, s);
            } else {
                auto t = task!decodeSlice(fd, c, s);
                _pool.put(t);
                _output_queue.put(t);
            }
        }
    }

    this(R slices, TaskPool pool) {
        _slices = slices;
        _pool = pool;
        import std.algorithm : max;
        auto size = max(1, pool.size) * 2;
        _input_queue = RoundBuf!CramSlice(size);
        _output_queue = RoundBuf!DecodeTask(size);
        while (!_slices.empty && !_input_queue.full)
            putNextSliceIntoQueue();

        popFront();
    }

    CramSlice front;
    bool empty;

    void popFront() {
        if (_input_queue.empty) {
            empty = true;
            return;
        }

        version (serial) {} else {
            auto t = _output_queue.front;
            _output_queue.popFront();
            t.yieldForce(); // now _input_queue.front is decoded
        }

        front = _input_queue.front;
        _input_queue.popFront();
        // debug writeln("GET slice #", front.id + 1, " from queue");
        if (!_slices.empty)
            putNextSliceIntoQueue();
    }
}

auto decode(R)(R slices, std.parallelism.TaskPool pool)
    if(isInputRange!R && is(ElementType!R == CramSlice))
{
    return CramSliceDecoder!R(slices, pool);
}

auto slices(CramFd fd, CramContainerFilter cf, UndecodedSliceFilter sf,
            TaskPool pool)
{
    auto undecoded_slices = new UndecodedSliceRange(fd, cf, sf);
    return undecoded_slices.decode(pool);
}

// workaround for LDC #795 - doesn't use short-circuit evaluation
// TODO remove when 0.15.0 release comes out
auto joiner2(RoR)(RoR r)
if (isInputRange!RoR && isInputRange!(ElementType!RoR)) {
    static struct Result {
    private:
        RoR _items;
        ElementType!RoR _current;
    public:
        this(RoR r) {
            _items = r;

            _current = ElementType!RoR.init;

            while (!_items.empty) {
                _current = _items.front;
                if (!_current.empty)
                    break;
                _items.popFront();
            }
        }

        @property auto empty() { return _items.empty; }
        @property auto ref front() { assert(!empty); return _current.front; }
        void popFront()
        {
            assert(!_current.empty);
            _current.popFront();
            if (_current.empty)
            {
                assert(!_items.empty);
                _items.popFront();
                while (!_items.empty) {
                    _current = _items.front;
                    if (!_current.empty)
                        break;
                    _items.popFront();
                }
            }
        }
    }
    return Result(r);
}
