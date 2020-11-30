// Written in the D programming language

module contrib.undead.internal.file;

// Copied from std.file. undead doesn't have access to it, but some modules
// in undead used std.file.deleteme when they were in Phobos, so this gives
// them access to a version of it.
public @property string deleteme() @safe
{
    import std.conv : to;
    import std.file : tempDir;
    import std.path : buildPath;
    import std.process : thisProcessID;
    static _deleteme = "deleteme.dmd.unittest.pid";
    static _first = true;

    if(_first)
    {
        _deleteme = buildPath(tempDir(), _deleteme) ~ to!string(thisProcessID);
        _first = false;
    }

    return _deleteme;
}

