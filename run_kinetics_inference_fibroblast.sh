#!/bin/bash

mkdir -p kinetics/output_fibroblast/split

## tpm
parallel -j 2 ./txburstML.py ::: kinetics/input_fibroblast/*.csv
mv *.pkl kinetics/output_fibroblast/
parallel -j 2 ./pkl_to_tsv.py ::: kinetics/output_fibroblast/*.pkl

./split_kinetics.R kinetics/output_fibroblast/ kinetics/output_fibroblast/split/fibroblast 1,2,3
