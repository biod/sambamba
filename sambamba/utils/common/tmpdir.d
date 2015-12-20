module sambamba.utils.common.tmpdir;
import std.random;
import std.path;
import std.file;
import std.process;
import std.conv;

string defaultTmpDir() {
    return environment.get("TMPDIR", "/tmp");
}

string randomSubdir(string tmpdir, string prefix="") {
    auto gen = Random(unpredictableSeed);
    char[4] subdirname;
    string result;
    do { 
      foreach (ref c; subdirname[])
          c = uniform!"[]"('a', 'z', gen);
      result = buildPath(tmpdir, "sambamba-pid" ~ thisProcessID.to!string ~ 
                                 "-" ~ prefix ~ subdirname[]);
    } while (std.file.exists(result));
    std.file.mkdirRecurse(result);
    return result;
}
