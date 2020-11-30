module bio.core.utils.tmpfile;

import std.path;
import std.process;

/// Absolute path of temporary file.
///
/// Params:
///     filename - base name
///     tmpdir   - temporary directory
///                                   
string tmpFile(string filename, string tmpdir=null) {
    if (tmpdir != null) {
        return buildPath(tmpdir, filename);
    }
version(Windows)
    return buildPath(std.process.getenv("TEMP"), filename);
else version(Posix)
    return "/tmp/" ~ filename;
}
