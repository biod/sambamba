module sambamba.utils.common.tmpdir;
import std.random;
import std.path;
import std.file;
import std.process;
import std.conv;

string randomSubdir(string tmpdir) {
    auto gen = Random(unpredictableSeed);
    char[4] subdirname;
    string result;
    do { 
      foreach (ref c; subdirname[])
          c = uniform!"[]"('a', 'z', gen);
      result = buildPath(tmpdir, "sambamba-pid" ~ thisProcessID.to!string ~ 
                                 "-" ~ subdirname[]);
    } while (std.file.exists(result));
    std.file.mkdirRecurse(result);
    return result;
}
