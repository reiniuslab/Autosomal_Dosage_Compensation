#!/bin/bash

mkdir -p kinetics/output_xpress/tpm/split kinetics/output_xpress/normumi/split

## tpm
parallel -j 2 ./txburstML.py ::: kinetics/input_xpress/tpm/*.csv
mv *.pkl kinetics/output_xpress/tpm/
parallel -j 2 ./pkl_to_tsv.py ::: kinetics/output_xpress/tpm/*.pkl

./split_kinetics.R kinetics/output_xpress/tpm/ kinetics/output_xpress/tpm/split/clonal_xpress_tpm 2,3,4

## normumi
parallel -j 2 ./txburstML.py ::: kinetics/input_xpress/normumi/*.csv
mv *.pkl kinetics/output_xpress/normumi/
parallel -j 2 ./pkl_to_tsv.py ::: kinetics/output_xpress/normumi/*.pkl

./split_kinetics.R kinetics/output_xpress/normumi/ kinetics/output_xpress/normumi/split/clonal_xpress_normumi 2,3,4
