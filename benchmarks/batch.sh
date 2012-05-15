#!/bin/bash

./run_benchmarks > serial.results

for n in {1..8} 
do
    ./run_benchmarks_mt $n > $n.results
done
