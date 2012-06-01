import bamfile;
import tagvalue;

import std.c.stdio;
import std.conv;
import std.algorithm : map;

void main(string[] args) {

    /* This is something like 'samtools view'  */

    auto bam = BamFile(args[1]);
    foreach (alignment; bam.alignments) {
        printf("%.*s\t%u\t",
                 alignment.read_name.length, alignment.read_name.ptr,
                 alignment.flag);
        if (alignment.ref_id == -1) {
            printf("*");
        } else {
            auto refseqname = bam.reference_sequences[alignment.ref_id].name;
            printf("%.*s", refseqname.length, refseqname.ptr);
        }
        printf("\t%d\t%u\t", alignment.position + 1, alignment.mapping_quality);
        if (alignment.cigar.length == 0) {
            printf("*\t");
        } else {
            // avoid memory allocation and NOT use cigar_string()
            foreach (cigar_op; alignment.cigar) {
                printf("%i%c", cigar_op.length, cigar_op.operation);
            }
            printf("\t");
        }
        if (alignment.next_ref_id == alignment.ref_id) {
            if (alignment.next_ref_id == -1) {
                printf("*\t");
            } else {
                printf("=\t");
            }
        } else {
            if (alignment.next_ref_id == -1 ||
                bam.reference_sequences[alignment.next_ref_id].name.length == 0)
            {
                printf("*\t");
            } else {
                auto refseqname = bam.reference_sequences[alignment.next_ref_id].name;
                printf("%.*s\t", refseqname.length, refseqname.ptr);
            }
        }
        printf("%i\t%i\t", alignment.next_pos + 1, alignment.template_length);
        if (alignment.raw_sequence_data.length == 0) {
            printf("*\t");
        } else {
            foreach(c; alignment.sequence)
                printf("%c", c);
            printf("\t");
        }
        if (alignment.phred_base_quality.length == 0 || 
            alignment.phred_base_quality[0] == '\xFF')
        {
            printf("*");
        } else {
            foreach (c; alignment.phred_base_quality) {
                printf("%c", c + 33);
            }
        }

        foreach (k, v; alignment.tags) {
            auto s = v.to_sam();
            printf("\t%.*s:%.*s", k.length, k.ptr, s.length, s.ptr);
        }

        printf("\n");
    }
}
