#!/bin/bash
# Copy the R reference outputs for the MiSeq SOP mock sample (forward reads)
# into test/fixtures, for the R-free regression tests in test/test_r_fixtures.jl.
#   bash make_fixtures.sh <r_out_dir> <fixtures_dir>
set -e
R=$1; F=$2; S=Mock_S280_L001
mkdir -p "$F"
cp "$R/filtered/${S}_F_filt.fastq.gz" "$F/mock_F.fastq.gz"
cp "$R/err_F.tsv"                     "$F/mock_F_err.tsv"
cp "$R/derep_F_${S}.tsv"              "$F/mock_F_derep.tsv"
cp "$R/dada_F_${S}.tsv"               "$F/mock_F_dada.tsv"
cp "$R/dadamap_F_${S}.txt"            "$F/mock_F_map.txt"
cp "$R/dadapval_F_${S}.txt"           "$F/mock_F_pval.txt"
ls -la "$F"
