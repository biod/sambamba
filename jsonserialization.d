/**
  Alignment serialization to JSON
*/
module jsonserialization;

import bamfile;
import tagvalue;
import utils.format;

import std.conv;
import std.algorithm;
import std.typecons;
import std.stdio;
import std.traits;
import std.c.stdlib;

import std.array;

/** Representation of tag value in SAM format
 
    Example:
    ----------
    Value v = 2.7;
    assert(toJson(v) == "2.7");

    v = [1, 2, 3];
    assert(toJson(v) == "[1,2,3]");
    ----------
*/
string toJson(Value v) {
    char[] buf;
    buf.reserve(16);
    jsonSerialize(v, buf);
    return cast(string)buf;
}

/// Print SAM representation to FILE* or append it to char[]/char* 
/// (in char* case it's your responsibility to allocate enough memory)
void jsonSerialize(S)(Value v, ref S stream) {

    if (v.is_numeric_array) {
        string toSamNumericArrayHelper() {
            char[] cases;
            foreach (t; ArrayElementTagValueTypes) {
                char[] printexpr = "putinteger(stream, elem);".dup;
                if (t.ch == 'f') {
                    printexpr = "append(stream, \"%g\", elem);".dup;
                }
                cases ~= `case '`~t.ch~`':` ~
                         `  putcharacter(stream, '[');`~
                         t.ValueType.stringof~`[] arr = cast(`~t.ValueType.stringof~`[])v;`~
                         `  if (arr.length != 0) { { auto elem = arr[0];`~printexpr~`}`~
                         `  foreach (elem; arr[1..$]) { putcharacter(stream, ',');`~printexpr~`}`~
                         `  }putcharacter(stream, ']');`~
                         `  return;`.dup;
            }
            return "switch (v.bam_typeid) { " ~ cases.idup ~ "default: assert(0); }";
        }
        mixin(toSamNumericArrayHelper());
    }
    if (v.is_integer) {
        switch (v.bam_typeid) {
            case 'c': putinteger(stream, to!byte(v));   return;
            case 'C': putinteger(stream, to!ubyte(v));  return; 
            case 's': putinteger(stream, to!short(v));  return; 
            case 'S': putinteger(stream, to!ushort(v)); return; 
            case 'i': putinteger(stream, to!int(v));    return; 
            case 'I': putinteger(stream, to!uint(v));   return; 
            default: assert(0);
        }
    }
    if (v.is_float) {
        append(stream, "%g", to!float(v));
        return;
    }
    switch (v.bam_typeid) {
        case 'Z', 'H':
            putcharacter(stream, '"');
            foreach (char c; cast(string)v) {
                if (c == '"' || c == '\\')
                    putcharacter(stream, '\\');
                putcharacter(stream, c);
            }
            putcharacter(stream, '"');
            return;
        case 'A': 
            putcharacter(stream, '"');
            auto c = to!char(v);
            if (c == '"' || c == '\\')
                putcharacter(stream, '\\');
            putcharacter(stream, c);
            putcharacter(stream, '"');
            return;
        default: assert(0);
    }
}

/// Get JSON representation of an alignment.
///
/// Requires providing information about reference sequences,
/// since alignment struct itself doesn't hold their names, only integer ids.
/// 
/// Example:
/// -------------
/// toJson(alignment, bam.reference_sequences);
/// -------------
string toJson(Alignment alignment, ReferenceSequenceInfo[] info) {
    char[] buf;
    buf.reserve(512);
    jsonSerialize(alignment, info, buf);
    return cast(string)buf;
}

/// Serialize $(D alignment) to FILE* or append it to char[]/char* 
/// (in char* case it's your responsibility to allocate enough memory)
void jsonSerialize(S)(Alignment alignment, ReferenceSequenceInfo[] info, ref S stream) 
    if (is(Unqual!S == FILE*) || is(Unqual!S == char*) || is(Unqual!S == char[]))
{

    // Notice: it is extremely important to exclude pointers,
    // otherwise you'll get recursion and stack overflow.
    static if (__traits(compiles, alloca(0)) && !is(Unqual!S == char*)) {

        immutable ALLOCA_THRESHOLD = 5000;

        if (alignment.size_in_bytes < ALLOCA_THRESHOLD) {

            char* buffer = cast(char*)alloca(alignment.size_in_bytes * 10 + 1000);

            if (buffer != null) {
                char* p = buffer; // this pointer will be modified
                jsonSerialize(alignment, info, p);
                putstring(stream, buffer[0 .. p - buffer]);
                return;
            } else {
                debug {
                    import std.stdio;
                    writeln("WARNING: pointer allocated with alloca was null");
                }
            }
        }
    }
   
    putstring(stream, `{"qname":"`);
    putstring(stream, alignment.read_name);
    putstring(stream, `","flag":`);
    putinteger(stream, alignment.flag);
    putstring(stream, `,"rname":"`);

    if (alignment.ref_id == -1) {
        putstring(stream, `*","pos":`);
    } else {
        putstring(stream, info[alignment.ref_id].name);
        putstring(stream, `","pos":`);
    }

    putinteger(stream, alignment.position + 1);

    putstring(stream, `,"mapq":`);
    putinteger(stream, alignment.mapping_quality);

    putstring(stream, `,"cigar":"`);

    if (alignment.cigar.length == 0) {
        putstring(stream, `*","rnext":`);
    } else {
        foreach (cigar_op; alignment.cigar) {
            putinteger(stream, cigar_op.length);
            putcharacter(stream, cigar_op.operation);
        }
        putstring(stream, `","rnext":`);
    }
    if (alignment.next_ref_id == alignment.ref_id) {
        if (alignment.next_ref_id == -1) {
            putstring(stream, `"*","pnext":`);
        } else {
            putstring(stream, `"=","pnext":`);
        }
    } else {
        if (alignment.next_ref_id == -1 ||
            info[alignment.next_ref_id].name.length == 0)
        {
            putstring(stream, `"*","pnext":`);
        } else {
            putcharacter(stream, '"');
            putstring(stream, info[alignment.next_ref_id].name);
            putstring(stream, `","pnext":`);
        }
    }

    putinteger(stream, alignment.next_pos + 1);

    putstring(stream, `,"tlen":`);
    putinteger(stream, alignment.template_length);

    putstring(stream, `,"seq":"`);
    if (alignment.raw_sequence_data.length == 0) {
        putstring(stream, `*","qual":`);
    } else {
        foreach(char c; alignment.sequence()) {
            putcharacter(stream, c);
        }
        putstring(stream, `","qual":`);
    }
    if (alignment.phred_base_quality.length == 0 || 
        alignment.phred_base_quality[0] == '\xFF')
    {
        putstring(stream, `"*","tags":{`);
    } else {
        putcharacter(stream, '"');
        foreach(char c; alignment.phred_base_quality) {
            char to_write = cast(char)(c + 33);
            if (to_write == '"')
                putstring(stream, `\"`);
            else if (to_write == '\\')
                putstring(stream, "\\\\");
            else
                putcharacter(stream, to_write);
        }
        putstring(stream, `","tags":{`);
    }
   
    bool not_first = false;
    foreach (k, v; alignment.tags) {
        assert(k.length == 2);

        if (not_first) putstring(stream, `,"`);
        else putcharacter(stream, '"');
        putstring(stream, k);
        putstring(stream, `":`);

        not_first = true;

        jsonSerialize(v, stream);
    }

    putstring(stream, `}}`);
    return;
}
