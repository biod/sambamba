module utils.stream;

public import std.stream;
private import core.stdc.stdio;

version(Posix){
    private import core.sys.posix.unistd;
}

final class File: std.stream.File {
    this(string filename) {
        // Issue 8528 workaround
        auto file = fopen(toStringz(filename), "rb");
        super(core.stdc.stdio.fileno(file), FileMode.In);
    }

    override ulong seek(long offset, SeekPos rel) {
        assertSeekable();
        auto hFile = handle;
        version (Windows) {
          int hi = cast(int)(offset>>32);
          uint low = SetFilePointer(hFile, cast(int)offset, &hi, rel);
          if ((low == INVALID_SET_FILE_POINTER) && (GetLastError() != 0))
            throw new SeekException("unable to move file pointer");
          ulong result = (cast(ulong)hi << 32) + low; 
        } else version (Posix) {
            // Phobos casts offset to int, leading to throwing an exception
            // on large files
            auto result = lseek(hFile, cast(off_t)offset, rel);
        }    
        if (result == cast(typeof(result))-1)
          throw new SeekException("unable to move file pointer");
        readEOF = false;
        return cast(ulong)result;
      }
}
