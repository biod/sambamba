#!/bin/bash

FILENAME=benchmark.dat

echo -e "compiler\tn\tcompressed\tdecompressed\tunparsed\tparsed" > $FILENAME

echo -ne "dmd\t" >> $FILENAME
./run_benchmarks >> $FILENAME
#echo -ne "gdc\t" >> $FILENAME
#./run_benchmarks_gdc >> $FILENAME

for i in {2..8}
do
    echo -ne "dmd\t" >> $FILENAME
    ./run_benchmarks_mt $i >> $FILENAME
#    echo -ne "gdc\t" >> $FILENAME
#    ./run_benchmarks_mt_gdc $i >> $FILENAME
done
