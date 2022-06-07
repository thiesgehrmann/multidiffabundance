#!/usr/bin/env bash
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

abundance="$1"
metadata="$2"
formula="$3"
outfile="$4"

Rscript "$SCRIPT_DIR/run_aldex.R" "$abundance" "$metadata" "$formula" "$outfile"