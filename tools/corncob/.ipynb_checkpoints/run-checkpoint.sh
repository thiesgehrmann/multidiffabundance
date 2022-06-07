#!/usr/bin/env bash

abundance="$1"
metadata="$2"
formula="$3"
outfile="$4"

Rscript run_deseq2.R "$abundance" "$metadata" "$formula" "$outfile"