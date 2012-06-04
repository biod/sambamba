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

/** Representation of tag value in SAM format
 
    Example:
    ----------
    Value v = 2.7;
    assert(to_sam(v) == "f:2.7");

    v = [1, 2, 3];
    assert(to_sam(v) == "B:i,1,2,3");
*/
string to_sam(Value v) {
    char[] buf;
    serialize(v, buf);
    return cast(string)buf;
}

/// Prints SAM representation into a file or a buffer
void serialize(S)(Value v, ref S stream) {
    if (v.is_numeric_array) {
        string toSamNumericArrayHelper() {
            char[] cases;
            foreach (t; ArrayElementTagValueTypes) {
                char[] format = "%d".dup;
                if (t.ch == 'f') {
                    format = "%g".dup;
                }
                cases ~= `case '`~t.ch~`':` ~
                         `  putstring(stream, "B:`~t.ch~`");`~
                         t.ValueType.stringof~`[] arr = to!(`~t.ValueType.stringof~`[])(v);`~
                         `  foreach (elem; arr) {`~
                         `      append(stream, ",`~format~`", elem);`~
                         `  }`~
                         `  return;`.dup;
            }
            return "switch (v.bam_typeid) { " ~ cases.idup ~ "default: assert(0); }";
        }
        mixin(toSamNumericArrayHelper());
    }
    if (v.is_integer) {
        switch (v.bam_typeid) {
            case 'c': append(stream, "i:%d", to!byte(v)); return;
            case 'C': append(stream, "i:%d", to!ubyte(v)); return;
            case 's': append(stream, "i:%d", to!short(v)); return;
            case 'S': append(stream, "i:%d", to!ushort(v)); return;
            case 'i': append(stream, "i:%d", to!int(v)); return;
            case 'I': append(stream, "i:%d", to!uint(v)); return;
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
            putstring(stream, to!string(v));
            return;
        case 'A': 
            append(stream, "A:%c", to!char(v));
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
    char[] buf;
    serialize(alignment, info, buf);
    return cast(string)buf;
}

/// Serialize an alignment into a stream or file
void serialize(S)(Alignment alignment, ReferenceSequenceInfo[] info, ref S stream) {
    putstring(stream, alignment.read_name);
    append(stream, "\t%d\t", alignment.flag);
    if (alignment.ref_id == -1) {
        putcharacter(stream, '*');
    } else {
        putstring(stream, info[alignment.ref_id].name);
    }
    append(stream, "\t%d\t%d\t", alignment.position + 1, alignment.mapping_quality);
    if (alignment.cigar.length == 0) {
        putstring(stream, "*\t");
    } else {
        // avoid memory allocation and NOT use cigar_string()
        foreach (cigar_op; alignment.cigar) {
            append(stream, "%d%c", cigar_op.length, cigar_op.operation);
        }
        putcharacter(stream, '\t');
    }
    if (alignment.next_ref_id == alignment.ref_id) {
        if (alignment.next_ref_id == -1) {
            putcharacter(stream, '*');
        } else {
            putcharacter(stream, '=');
        }
    } else {
        if (alignment.next_ref_id == -1 ||
            info[alignment.next_ref_id].name.length == 0)
        {
            putcharacter(stream, '*');
        } else {
            putstring(stream, info[alignment.next_ref_id].name);
        }
    }
    append(stream, "\t%d\t%d\t", alignment.next_pos + 1, alignment.template_length);
    if (alignment.raw_sequence_data.length == 0) {
        putstring(stream, "*\t");
    } else {
        foreach(char c; alignment.sequence())
            putcharacter(stream, c);
        putcharacter(stream, '\t');
    }
    if (alignment.phred_base_quality.length == 0 || 
        alignment.phred_base_quality[0] == '\xFF')
    {
        putcharacter(stream, '*');
    } else {
        foreach (char c; alignment.phred_base_quality) {
            putcharacter(stream, cast(char)(c + 33));
        }
    }
    
    foreach (k, v; alignment.tags) {
        assert(k.length == 2);
        append(stream, "\t%c%c:", k[0], k[1]);
        serialize(v, stream);
    }
}
