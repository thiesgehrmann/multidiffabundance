#!/usr/bin/env bash

abundance="$1"
metadata="$2"
formula="$3"
outfile="$4"

Rscript run_corncob.R "$abundance" "$metadata" "$formula" "$outfile"