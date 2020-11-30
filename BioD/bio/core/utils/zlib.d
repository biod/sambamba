// minimized version of etc.c.zlib
module bio.core.utils.zlib;

import core.stdc.config;

extern (C) {

const char[] ZLIB_VERSION = "1.2.3";
const ZLIB_VERNUM = 0x1230;

alias void* function (void* opaque, uint items, uint size) alloc_func;
alias void  function (void* opaque, void* address) free_func;

struct z_stream
{
    ubyte    *next_in;  /* next input byte */
    uint     avail_in;  /* number of bytes available at next_in */
    c_ulong  total_in;  /* total nb of input bytes read so far */

    ubyte    *next_out; /* next output byte should be put there */
    uint     avail_out; /* remaining free space at next_out */
    c_ulong  total_out; /* total nb of bytes output so far */

    char     *msg;      /* last error message, NULL if no error */
    void*    state;     /* not visible by applications */

    alloc_func zalloc;  /* used to allocate the internal state */
    free_func  zfree;   /* used to free the internal state */
    void*      opaque;  /* private data object passed to zalloc and zfree */

    int     data_type;  /* best guess about the data type: binary or text */
    c_ulong adler;      /* adler32 value of the uncompressed data */
    c_ulong reserved;   /* reserved for future use */
}

alias z_stream* z_streamp;

/*
     gzip header information passed to and from zlib routines.  See RFC 1952
  for more details on the meanings of these fields.
*/
struct gz_header {
    int     text;       /* true if compressed data believed to be text */
    c_ulong time;       /* modification time */
    int     xflags;     /* extra flags (not used when writing a gzip file) */
    int     os;         /* operating system */
    byte    *extra;     /* pointer to extra field or Z_NULL if none */
    uint    extra_len;  /* extra field length (valid if extra != Z_NULL) */
    uint    extra_max;  /* space at extra (only when reading header) */
    byte    *name;      /* pointer to zero-terminated file name or Z_NULL */
    uint    name_max;   /* space at name (only when reading header) */
    byte    *comment;   /* pointer to zero-terminated comment or Z_NULL */
    uint    comm_max;   /* space at comment (only when reading header) */
    int     hcrc;       /* true if there was or will be a header crc */
    int     done;       /* true when done reading gzip header (not used
                           when writing a gzip file) */
}

alias gz_header* gz_headerp;

                        /* constants */

enum
{
        Z_NO_FLUSH      = 0,
        Z_PARTIAL_FLUSH = 1, /* will be removed, use Z_SYNC_FLUSH instead */
        Z_SYNC_FLUSH    = 2,
        Z_FULL_FLUSH    = 3,
        Z_FINISH        = 4,
        Z_BLOCK         = 5,
        Z_TREES         = 6,
}
/* Allowed flush values; see deflate() and inflate() below for details */

enum
{
        Z_OK            = 0,
        Z_STREAM_END    = 1,
        Z_NEED_DICT     = 2,
        Z_ERRNO         = -1,
        Z_STREAM_ERROR  = -2,
        Z_DATA_ERROR    = -3,
        Z_MEM_ERROR     = -4,
        Z_BUF_ERROR     = -5,
        Z_VERSION_ERROR = -6,
}
/* Return codes for the compression/decompression functions. Negative
 * values are errors, positive values are used for special but normal events.
 */

enum
{
        Z_NO_COMPRESSION         = 0,
        Z_BEST_SPEED             = 1,
        Z_BEST_COMPRESSION       = 9,
        Z_DEFAULT_COMPRESSION    = -1,
}
/* compression levels */

enum
{
        Z_FILTERED            = 1,
        Z_HUFFMAN_ONLY        = 2,
        Z_RLE                 = 3,
        Z_FIXED               = 4,
        Z_DEFAULT_STRATEGY    = 0,
}
/* compression strategy; see deflateInit2() below for details */

enum
{
        Z_BINARY   = 0,
        Z_TEXT     = 1,
        Z_UNKNOWN  = 2,

        Z_ASCII    = Z_TEXT
}
/* Possible values of the data_type field (though see inflate()) */

enum
{
        Z_DEFLATED   = 8,
}
/* The deflate compression method (the only one supported in this version) */

const int Z_NULL = 0;  /* for initializing zalloc, zfree, opaque */

                        /* basic functions */

char* zlibVersion();
int deflateInit(z_streamp strm, int level)
{
    return deflateInit_(strm, level, ZLIB_VERSION.ptr, z_stream.sizeof);
}

int deflate(z_streamp strm, int flush);
int deflateEnd(z_streamp strm);

int inflateInit(z_streamp strm)
{
    return inflateInit_(strm, ZLIB_VERSION.ptr, z_stream.sizeof);
}
int inflate(z_streamp strm, int flush);
int inflateEnd(z_streamp strm);

int deflateInit2(z_streamp strm,
                 int  level,
                 int  method,
                 int  windowBits,
                 int  memLevel,
                 int  strategy)
{
    return deflateInit2_(strm, level, method, windowBits, memLevel,
                         strategy, ZLIB_VERSION.ptr, z_stream.sizeof);
}


int deflateBound(z_streamp strm, size_t sourceLen);

int inflateInit2(z_streamp strm, int windowBits)
{
    return inflateInit2_(strm, windowBits, ZLIB_VERSION.ptr, z_stream.sizeof);
}

int compress(ubyte* dest,
             size_t* destLen,
             ubyte* source,
             size_t sourceLen);

int compress2(ubyte* dest,
              size_t* destLen,
              ubyte* source,
              size_t sourceLen,
              int level);

size_t compressBound(size_t sourceLen);

int uncompress(ubyte* dest,
               size_t* destLen,
               ubyte* source,
               size_t sourceLen);

alias void* gzFile;
alias int z_off_t;              // file offset

gzFile gzopen(char* path, char* mode);
gzFile gzdopen(int fd, char* mode);

int gzsetparams(gzFile file, int level, int strategy);
int gzread(gzFile file, void* buf, uint len);
int gzwrite(gzFile file, void* buf, uint len);
int gzprintf(gzFile file, char* format, ...);
int gzputs(gzFile file, char* s);
char* gzgets(gzFile file, char* buf, int len);
int gzputc(gzFile file, int c);
int    gzgetc(gzFile file);
int gzungetc(int c, gzFile file);
int gzflush(gzFile file, int flush);
z_off_t gzseek(gzFile file, z_off_t offset, int whence);
int gzrewind(gzFile file);
z_off_t  gztell(gzFile file);
int gzeof(gzFile file);
int gzdirect(gzFile file);
int gzclose(gzFile file);
char* gzerror(gzFile file, int *errnum);
void gzclearerr (gzFile file);
 uint adler32  (uint adler, ubyte *buf, uint len);
uint adler32_combine(uint adler1, uint adler2, z_off_t len2);
uint crc32(uint crc, ubyte *buf, uint len);
uint crc32_combine (uint crc1, uint crc2, z_off_t len2);

int deflateInit_(z_streamp strm,
                 int level,
                 const char* versionx,
                 int stream_size);

int inflateInit_(z_streamp strm,
                 const char* versionx,
                 int stream_size);

int deflateInit2_(z_streamp strm,
                  int level,
                  int method,
                  int windowBits,
                  int memLevel,
                  int strategy,
                  const char* versionx,
                  int stream_size);

int inflateBackInit_(z_stream* strm,
                     int windowBits,
                     ubyte* window,
                     const char* z_version,
                     int stream_size);

int inflateInit2_(z_streamp strm,
                  int windowBits,
                  const char* versionx,
                  int stream_size);

char* zError(int err);
int inflateSyncPoint(z_streamp z);
uint* get_crc_table();

}

class ZlibException : Exception
{
    this(int errnum) {
        auto msg = "[zlib] " ~ messageFromErrnum(errnum);
        super(msg);
    }

    this(string func, int errnum) {
        auto msg = "[zlib/" ~ func ~ "] " ~ messageFromErrnum(errnum);
        super(msg);
    }

    private string messageFromErrnum(int errnum) {
        switch (errnum)
        {
            case Z_STREAM_END:      msg = "stream end"; break;
            case Z_NEED_DICT:       msg = "need dict"; break;
            case Z_ERRNO:           msg = "errno"; break;
            case Z_STREAM_ERROR:    msg = "stream error"; break;
            case Z_DATA_ERROR:      msg = "data error"; break;
            case Z_MEM_ERROR:       msg = "mem error"; break;
            case Z_BUF_ERROR:       msg = "buf error"; break;
            case Z_VERSION_ERROR:   msg = "version error"; break;
            default:                msg = "unknown error";  break;
        }
        return msg;
    }
}

uint crc32(uint crc, const(void)[] buf)
{
    return crc32(crc, cast(ubyte *)buf.ptr, cast(uint)(buf.length));
}

