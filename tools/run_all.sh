#!/usr/bin/env bash

count_data="~/repos//UA_isala/flow1_redo_revision/data/genus_level_counts.tsv"
meta_data="~/repos//UA_isala/flow1_redo_revision/data/metadata.tsv"
formula_file="~/repos//multidiffabundance/tools/list_of_functions.txt"

source "`dirname $CONDA_EXE`/../etc/profile.d/conda.sh"

conda activate limma && aldex2/run.sh "$count_data" "$meta_data" "$formula_file" "output.aldex2.tsv"
#conda activate ancom2 && ancombc/run.sh "$count_data" "$meta_data" "$formula_file" "output.ancombc.tsv"
#conda activate deseq2 && deseq2/run.sh "$count_data" "$meta_data" "$formula_file" "output.deseq2.tsv"
#conda activate limma && limma/run.sh "$count_data" "$meta_data" "$formula_file" "output.limma.tsv"
#corncob/run.sh "$count_data" "$meta_data" "$formula_file" "output.corncob.tsv"
#lmclr/run.sh "$count_data" "$meta_data" "$formula_file" "output.lmclr.tsv"
