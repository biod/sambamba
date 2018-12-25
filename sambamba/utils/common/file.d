module sambamba.utils.common.file;
import contrib.undead.stream;
import std.stdio;

BufferedFile bufferedFile(string fn, contrib.undead.stream.FileMode mode, size_t buffer_size=8192) {
    if (fn == "-")
        return new BufferedFile(std.stdio.stdout.fileno, mode, buffer_size);
    else
        return new BufferedFile(fn, mode, buffer_size);
}
