import bamfile;
import tagvalue;

import std.stdio;
import std.conv;
import std.algorithm : map;

string to_sam(Value v) {
    if (v.type == typeid(string)) {
        return "Z:" ~ to!string(v);
    } else if (v.type == typeid(float)) {
        return "f:" ~ to!string(v);
    } else if (v.type == typeid(char)) {
        return "A:" ~ to!string(v);
    } else if (v.convertsTo!int()) {
        return "i:" ~ to!string(v);
    } else { // array
        return "Z:NOT_YET_IMPLEMENTED";
    }
}

void main(string[] args) {

    /* This is something like 'samtools view'.
       Unfortunately, quite slow, and the bottleneck
       seems to be usage of VariantN. 
    */

	auto bam = BamFile(args[1]);
	foreach (alignment; bam.alignments) {
        writef("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s", 
                 alignment.read_name, 
                 alignment.flag,
                 alignment.ref_id == -1 ?
                     "*" :
                     bam.reference_sequences[alignment.ref_id].name,
                 alignment.position + 1,
                 alignment.mapping_quality,
                 alignment.cigar_string.length == 0 ?
                     "*" :
                     alignment.cigar_string,
                 alignment.next_ref_id == alignment.ref_id ?
                     alignment.next_ref_id == -1 ?
                        "*" :
                        "=" :
                     alignment.next_ref_id == -1 ||
                     bam.reference_sequences[alignment.next_ref_id].name.length == 0 ?
                        "*" :
                        bam.reference_sequences[alignment.next_ref_id].name,
                 alignment.next_pos + 1,
                 alignment.template_length,
                 alignment.sequence.length == 0 ?
                     "*" :
                     alignment.sequence,
                 alignment.qual.length == 0 || alignment.qual[0] == '\xFF' ? 
                     "*" : 
                     to!string(map!"cast(char)(a+33)"(alignment.qual)));

		foreach (k, v; alignment.tags) {
            writef("\t%s:%s", k, to_sam(v));
		}

        write("\n");
	}
}
