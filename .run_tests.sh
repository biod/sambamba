./build/sambamba sort -t2 -n BioD/test/data/ex1_header.bam -o ex1_header.nsorted.bam
./build/sambamba view -t2 ex1_header.nsorted.bam | cut -f1 | LC_COLLATE=C sort -c
./build/sambamba sort -t2 ex1_header.nsorted.bam -o ex1_header.sorted.bam
./build/sambamba index -t2 ex1_header.sorted.bam
./build/sambamba slice ex1_header.sorted.bam chr1 -o /dev/null
./build/sambamba view -c ex1_header.sorted.bam chr2
./build/sambamba view -c ex1_header.sorted.bam chr1
./build/sambamba view -c ex1_header.sorted.bam '*' 
