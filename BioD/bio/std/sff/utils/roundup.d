module bio.std.sff.utils.roundup;

import std.traits;
import std.conv;

/// First number
T roundup(T)(T number) 
    if (isUnsigned!T)
{
    if (number % 8 == 0)
        return number;
    return to!T(number + (8 - number % 8));
}
