#!/usr/bin/env bash
SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]:-$0}"; )" &> /dev/null && pwd 2> /dev/null; )";

abundance="$1"
metadata="$2"
formula="$3"
outprefix="$4"

Rscript "$SCRIPT_DIR/run_maaslin2.R" "$abundance" "$metadata" "$formula" "$outprefix"