#!/usr/bin/env bash

count_data="~/repos//UA_isala/flow1_redo_revision/data/genus_level_counts.tsv"
meta_data="~/repos//UA_isala/flow1_redo_revision/data/metadata.tsv"
formula_file="~/repos//multidiffabundance/tools/list_of_functions.txt"

source "`dirname $CONDA_EXE`/../etc/profile.d/conda.sh"

#conda activate aldex2 && aldex2/run.sh "$count_data" "$meta_data" "$formula_file" "output.aldex2"
#conda activate ancom2 && ancombc/run.sh "$count_data" "$meta_data" "$formula_file" "output.ancombc"
#conda activate corncob && corncob/run.sh "$count_data" "$meta_data" "$formula_file" "output.corncob"

#conda activate deseq2 && deseq2/run.sh "$count_data" "$meta_data" "$formula_file" "output.deseq2"
#conda activate limma && limma/run.sh "$count_data" "$meta_data" "$formula_file" "output.limma"
conda activate limma && lmclr/run.sh "$count_data" "$meta_data" "$formula_file" "output.lmclr"
conda activate maaslin2 && maaslin2/run.sh "$count_data" "$meta_data" "$formula_file" "output.maaslin"
