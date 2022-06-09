#!/usr/bin/env bash

#count_data="~/repos//UA_isala/flow1_redo_revision/data/genus_level_counts.tsv"
#meta_data="~/repos//UA_isala/flow1_redo_revision/data/metadata.tsv"
#formula="~/repos//multidiffabundance/tools/list_of_functions.txt"
#formula="~ Reproductive.haschildren + Technical.dna_conc + General.Age"
#./mda.sh "$count_data" "$meta_data" "$formula" "test_output"

source "`dirname $CONDA_EXE`/../etc/profile.d/conda.sh"

###############################################################################

function usage(){
    echo "MDA: Multi Differential Abundance tool"
    echo "Usage: $0 [options] <count_data> <meta_data> <function> <outputprefix>"
    echo "  Options:"
    echo "    -a|--aldex2   : Run ALDex2"
    echo "    -A|--ancombc  : Run ANCOMBC"
    echo "    -c|--corncob  : Run corncob"
    echo "    -d|--deseq2   : Run DESeq2"
    echo "    -l|--limma    : Run limma"
    echo "    -L|--lmclr    : Run lm CLR"
    echo "    -m|--maaslin2 : Run Maaslin2"
    echo "    --complement  : run the complement of the current selection
                              (i.e. -c --complement would NOT run corncob)"
    echo "  Inputs:"
    echo "    count_data : File path to a sample*taxa (rows*columns) tab separated file.
                           First column is sample IDs"
    echo "    meta_data  : File path to a sample*determinant (rows*columns) tab separated file.
                           First column is sample IDs"
    echo "    function   : Either a string in R formula format (e.g. ~ a + b + c), or
                           a file path to a line-separated list of these formulas"
    echo "    outputprefix : A path to a location of the outputs"

}

###############################################################################
# PROCESS INPUTS

# We use "$@" instead of $* to preserve argument-boundary information
ARGS=$(getopt -o 'aAcdlLm' --long 'aldex2,ancombc,corncob,deseq2,limma,lmclr,maaslin2,complement' -- "$@") || exit
eval "set -- $ARGS"

aldex2=0
ancombc=0
corncob=0
deseq2=0
limma=0
lmclr=0
maaslin2=0
all=1
complement=0


while true; do
    case $1 in
      (-a|--aldex2)
            all=0; aldex2=1; shift;;
      (-A|--ancombc)
            all=0; ancombc=1; shift;;
      (-c|--corncob)
            all=0; corncob=1; shift;;
      (-d|--deseq2)
            all=0; deseq2=1; shift;;
      (-l|--limma)
            all=0; limma=1; shift;;
      (-L|--lmclr)
            all=0; lmclr=1; shift;;
      (-m|--maaslin2)
            all=0; maaslin2=1; shift;;
      (--complement)
            complement=1; shift;;
      (--)  shift; break;;
      (*)   usage $0; exit 1;;           # error
    esac
done

remaining=("$@")

if [ ! "${#remaining[@]}" -eq 4 ]; then
    usage $0
    exit 1
fi

function check_if_is_formula(){
    echo "as.formula('$1')" | R --no-save &> /dev/null
    ret=$?
    echo $ret
}

count_data=`realpath ${remaining[0]/#\~/$HOME}`
meta_data=`realpath ${remaining[1]/#\~/$HOME}`
formula=${remaining[2]}
outprefix=`realpath ${remaining[3]/#\~/$HOME}`

if [ `check_if_is_formula "$formula"` -eq 1 ]; then
    formula=`realpath ${formula/#\~/$HOME}`
fi

echo "$count_data"
echo "$meta_data"
echo "$formula"
echo "$outprefix"

function check_if_exists_else_fail(){
    if [ ! -f "$1" ]; then
        echo "Error: File '$1' does not exist."
        exit $2;
    fi
}

check_if_exists_else_fail "$count_data" 2
check_if_exists_else_fail "$meta_data" 2


mkdir -p "$outprefix"


if [ $complement -eq 1 ]; then
    aldex2=$((1-aldex2))
    ancombc=$((1-ancombc))
    corncob=$((1-corncob))
    deseq2=$((1-deseq2))
    limma=$((1-limma))
    lmclr=$((1-lmclr))
    maaslin2=$((1-maaslin2))
fi

if [ $all -eq 1 ] ; then
    aldex2=1
    ancombc=1
    corncob=1
    deseq2=1
    limma=1
    lmclr=1
    maaslin2=1
fi

###############################################################################


if [ $aldex2 -eq 1 ] ; then
    echo "RUN ALDEX2"
    conda activate aldex2 && aldex2/run.sh "$count_data" "$meta_data" "$formula" "$outprefix/aldex2" > "$outprefix/aldex2.stdout" 2> "$outprefix/aldex2.stderr"
fi

if [ $ancombc -eq 1 ] ; then
    echo "RUN ANCOMBC"
    #conda activate ancom2 && ancombc/run.sh "$count_data" "$meta_data" "$formula" "$outprefix/ancombc" > "$outprefix/ancombc.stdout" 2> "$outprefix/ancombc.stderr"
fi

if [ $corncob -eq 1 ] ; then
    echo "RUN CORNCOB"
    #conda activate corncob && corncob/run.sh "$count_data" "$meta_data" "$formula" "$outprefix/corncob" > "$outprefix/corncob.stdout" 2> "$outprefix/corncob.stderr"
fi
    
if [ $deseq2 -eq 1 ] ; then
    echo "RUN DESEQ2"
    #conda activate deseq2 && deseq2/run.sh "$count_data" "$meta_data" "$formula" "$outprefix/deseq2" > "$outprefix/deseq2.stdout" 2> "$outprefix/deseq2.stderr"
fi

if [ $limma -eq 1 ] ; then
    echo "RUN LIMMA"
    #conda activate limma && limma/run.sh "$count_data" "$meta_data" "$formula" "$outprefix/limma" > "$outprefix/limma.stdout" 2> "$outprefix/limma.stderr"
fi

if [ $lmclr -eq 1 ] ; then
    echo "RUN LMCLR"
    #conda activate limma && lmclr/run.sh "$count_data" "$meta_data" "$formula" "$outprefix/lmclr" > "$outprefix/lmclr.stdout" 2> "$outprefix/lmclr.stderr"
fi

if [ $maaslin2 -eq 1 ] ; then
    echo "RUN MAASLIN2"
    #conda activate maaslin2 && maaslin2/run.sh "$count_data" "$meta_data" "$formula" "$outprefix/maaslin2" > "$outprefix/maaslin2.stdout" 2> "$outprefix/maaslin2.stderr"
fi
