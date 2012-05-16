import bamfile;
import tagvalue;

import std.stdio;
import std.conv;
import std.algorithm : map;

void main(string[] args) {

    /* This is something like 'samtools view'  */

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
                 alignment.phred_base_quality.length == 0 || 
                     alignment.phred_base_quality[0] == '\xFF' ? 
                     "*" : 
                     to!string(map!"cast(char)(a+33)"(alignment.phred_base_quality)));

        foreach (k, v; alignment.tags) {
            writef("\t%s:%s", k, v.to_sam);
        }

        write("\n");
    }
}
