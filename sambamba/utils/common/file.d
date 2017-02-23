module sambamba.utils.common.file;
import undead.stream;
import std.stdio;

BufferedFile bufferedFile(string fn, undead.stream.FileMode mode, size_t buffer_size=8192) {
    if (fn == "-")
        return new BufferedFile(std.stdio.stdout.fileno, mode, buffer_size);
    else
        return new BufferedFile(fn, mode, buffer_size);
}
