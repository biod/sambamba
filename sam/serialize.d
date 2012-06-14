/**
  Alignment and tag value serialization to SAM format
*/
module sam.serialize;

import bamfile;
import tagvalue;
import utils.format;

import std.conv;
import std.algorithm;
import std.typecons;
import std.stdio;
import std.c.stdlib;

import std.array;

/** Representation of tag value in SAM format
 
    Example:
    ----------
    Value v = 2.7;
    assert(to_sam(v) == "f:2.7");

    v = [1, 2, 3];
    assert(to_sam(v) == "B:i,1,2,3");
*/
string to_sam(Value v) {
    auto buf = appender!(ubyte[])();
    buf.reserve(16);
    serialize(v, buf);
    return cast(string)buf.data;
}

/// Prints SAM representation into FILE* or Appender!(ubyte[])
void serialize(S)(Value v, ref S stream) {
    if (v.is_numeric_array) {
        string toSamNumericArrayHelper() {
            char[] cases;
            foreach (t; ArrayElementTagValueTypes) {
                char[] loopbody = "putcharacter(stream, ',');putinteger(stream, elem);".dup;
                if (t.ch == 'f') {
                    loopbody = "append(stream, \",%g\", elem);".dup;
                }
                cases ~= `case '`~t.ch~`':` ~
                         `  putstring(stream, "B:`~t.ch~`");`~
                         t.ValueType.stringof~`[] arr = cast(`~t.ValueType.stringof~`[])v;`~
                         `  foreach (elem; arr) {`~loopbody~`}`~
                         `  return;`.dup;
            }
            return "switch (v.bam_typeid) { " ~ cases.idup ~ "default: assert(0); }";
        }
        mixin(toSamNumericArrayHelper());
    }
    if (v.is_integer) {
        putstring(stream, "i:");
        switch (v.bam_typeid) {
            case 'c': putinteger(stream, to!byte(v)); return;
            case 'C': putinteger(stream, to!ubyte(v)); return;
            case 's': putinteger(stream, to!short(v)); return;
            case 'S': putinteger(stream, to!ushort(v)); return;
            case 'i': putinteger(stream, to!int(v)); return;
            case 'I': putinteger(stream, to!uint(v)); return;
            default: assert(0);
        }
    }
    if (v.is_float) {
        append(stream, "f:%g", to!float(v));
        return;
    }
    switch (v.bam_typeid) {
        case 'Z', 'H':
            putcharacter(stream, v.bam_typeid);
            putcharacter(stream, ':');
            putstring(stream, cast(string)v);
            return;
        case 'A': 
            putstring(stream, "A:");
            putcharacter(stream, to!char(v));
            return;
        default: assert(0);
    }
}

/// Get SAM representation of an alignment.
///
/// Requires providing information about reference sequences,
/// since alignment struct itself doesn't hold their names, only integer ids.
/// 
/// Example:
/// -------------
/// to_sam(alignment, bam.reference_sequences);
///
string to_sam(Alignment alignment, ReferenceSequenceInfo[] info) {
    auto buf = appender!(ubyte[])();
    buf.reserve(512);
    serialize(alignment, info, buf);
    return cast(string)buf.data;
}

/// Serialize an alignment into FILE* or Appender!(ubyte[])
void serialize(S)(Alignment alignment, ReferenceSequenceInfo[] info, ref S stream) {

    putstring(stream, alignment.read_name);
    putcharacter(stream, '\t');

    putinteger(stream, alignment.flag);
    putcharacter(stream, '\t');

    if (alignment.ref_id == -1) {
        putstring(stream, "*\t");
    } else {
        putstring(stream, info[alignment.ref_id].name);
        putcharacter(stream, '\t');
    }

    putinteger(stream, alignment.position + 1);
    putcharacter(stream, '\t');

    putinteger(stream, alignment.mapping_quality);
    putcharacter(stream, '\t');

    if (alignment.cigar.length == 0) {
        putstring(stream, "*\t");
    } else {
        // avoid memory allocation and NOT use cigar_string()
        foreach (cigar_op; alignment.cigar) {
            putinteger(stream, cigar_op.length);
            putcharacter(stream, cigar_op.operation);
        }
        putcharacter(stream, '\t');
    }
    if (alignment.next_ref_id == alignment.ref_id) {
        if (alignment.next_ref_id == -1) {
            putstring(stream, "*\t");
        } else {
            putstring(stream, "=\t");
        }
    } else {
        if (alignment.next_ref_id == -1 ||
            info[alignment.next_ref_id].name.length == 0)
        {
            putstring(stream, "*\t");
        } else {
            putstring(stream, info[alignment.next_ref_id].name);
            putcharacter(stream, '\t');
        }
    }

    putinteger(stream, alignment.next_pos + 1);
    putcharacter(stream, '\t');

    putinteger(stream, alignment.template_length);
    putcharacter(stream, '\t');

    // use alloca to avoid overhead of appending single chars
    auto len = alignment.sequence_length;
    char* str = null;
    if (len <= 1024) {
        str = cast(char*)alloca(alignment.sequence_length);
    }

    if (alignment.raw_sequence_data.length == 0) {
        putstring(stream, "*\t");
    } else {
        if (str is null) {
            foreach(char c; alignment.sequence())
                putcharacter(stream, c);
        } else {
            size_t i = 0;
            foreach(char c; alignment.sequence()) 
                *(str + i++) = c;
            putstring(stream, str[0 .. i]);
        }
        putcharacter(stream, '\t');
    }
    if (alignment.phred_base_quality.length == 0 || 
        alignment.phred_base_quality[0] == '\xFF')
    {
        putcharacter(stream, '*');
    } else {
        if (str is null) {
            foreach(char c; alignment.phred_base_quality) {
                putcharacter(stream, cast(char)(c + 33));
            }
        } else {
            foreach(size_t i, char c; alignment.phred_base_quality)
                *(str + i) = cast(char)(c + 33);
            putstring(stream, str[0 .. alignment.sequence_length]);
        }
    }
    
    foreach (k, v; alignment.tags) {
        assert(k.length == 2);
        putcharacter(stream, '\t');
        putstring(stream, k);
        putcharacter(stream, ':');
        serialize(v, stream);
    }
}
