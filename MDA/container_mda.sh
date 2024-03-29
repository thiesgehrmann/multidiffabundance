#!/usr/bin/env bash

###############################################################################

function usage(){
    echo "MDA: Multi Differential Abundance tool (Docker/singularity wrapper)"
    echo "Usage: $0 [options] <count_data> <meta_data> <function> <outputprefix>"
    echo "  Options:"
    echo "    -a|--aldex2   : Run ALDex2"
    echo "    -A|--ancombc  : Run ANCOMBC"
    echo "    -c|--corncob  : Run corncob"
    echo "    -d|--deseq2   : Run DESeq2"
    echo "    -l|--limma    : Run limma"
    echo "    -L|--lmclr    : Run lm CLR"
    echo "    -m|--maaslin2 : Run Maaslin2"
    echo "    --alpha       : Run alpha diversity test"
    echo "    --beta        : Run beta diversity test (permanova/adonis2)"
    echo "    --complement  : run the complement of the current selection"
    echo "                    (i.e. -c --complement would NOT run corncob)"
    echo "    -o|--options  : Arguments to the mda functions. Comma seperated."
    echo "                    (e.g. -o 'beta.permutations=10000')"
    echo "    --nocache     : Do not use a cache"
    echo "    --recache     : Re-do cached analyses"
    echo "    --help        : Display this help message"
    echo "    --docker      : Use docker instead of singularity"
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


ARGS=$(getopt -o 'aAcdlLmo:' --long 'alpha,beta,aldex2,ancombc,corncob,deseq2,limma,lmclr,maaslin2,complement,options:,docker,nocache' -- "$@") || exit
eval "set -- $ARGS"

singularity=1

args=""

while true; do
    case $1 in
      (-a|--aldex2)
            args="$args -a"; shift;;
      (-A|--ancombc)
            args="$args -A"; shift;;
      (-c|--corncob)
            args="$args -c"; shift;;
      (-d|--deseq2)
            args="$args -d"; shift;;
      (-l|--limma)
            args="$args -l"; shift;;
      (-L|--lmclr)
            args="$args -L"; shift;;
      (-m|--maaslin2)
            args="$args -m"; shift;;
      (--alpha)
            args="$args --alpha"; shift;;
      (--beta)
            args="$args --beta"; shift;;
      (--complement)
            args="$args --complement"; shift;;
      (-o|--options)
            args="$args -o '$2'"; shift 2;;
      (--nocache)
            args="$args --nocache"; shift;;
      (--recache)
            args="$args --recache"; shift;;
      (--help)
            usage $0; exit 1;;
      (--docker)
            singularity=0; shift;;
      (--)  shift; break;;
      (*)   usage $0; exit 1;;           # error
    esac
done

#############################
# Deal with file inputs

remaining=("$@")

if [ ! "${#remaining[@]}" -eq 4 ]; then
    usage $0
    exit 1
fi

count_data=`realpath ${remaining[0]/#\~/$HOME}`
meta_data=`realpath ${remaining[1]/#\~/$HOME}`
formula=${remaining[2]}
outprefix=`realpath ${remaining[3]/#\~/$HOME}`

#############################
# Prepare docker bind mount & soft link files to there
bindm_source="$outprefix"
bindm_target="/mda_inputs/"
mkdir -p "$bindm_source"

function check_if_exists_else_fail(){
    if [ ! -f "$1" ]; then
        echo "Error: File '$1' does not exist."
        exit $2;
    fi
}

check_if_exists_else_fail "$count_data" 2
rm -f "$bindm_source/"`basename "$count_data"`
ln "$count_data" "$bindm_source"
bindm_count_data="$bindm_target"`basename "$count_data"`

check_if_exists_else_fail "$meta_data" 2
rm -f "$bindm_source/"`basename "$meta_data"`
ln "$meta_data" "$bindm_source"
bindm_meta_data="$bindm_target"`basename "$meta_data"`


formula_check=`realpath "${formula/#\~/$HOME}"`
if [ -f "$formula_check" ] ; then
    rm -f "$bindm_source/"`basename "$formula_check"`
    ln "$formula_check" "$bindm_source"
    bindm_formula="$bindm_target"`basename "$formula_check"`
else
    bindm_formula="$formula"
fi


#############################
# Run things

docker_tag="thiesgehrmann/multidiffabundance:1"

if [ "$singularity" -eq 1 ]; then
    singularity exec \
        --bind "$bindm_source/":"$bindm_target" \
        --containall \
        docker://"${docker_tag}" \
        /usr/bin/mda \
        $args \
        "$bindm_count_data" \
        "$bindm_meta_data" \
        "$bindm_formula" \
        "$bindm_target/" | sed -e "s%$bindm_target%$bindm_source/%g"

else 
    docker run \
        --mount type=bind,source="$bindm_source/",target="$bindm_target" \
        "$docker_tag" \
        "/usr/bin/mda" \
        $args \
        "$bindm_count_data" \
        "$bindm_meta_data" \
        "$bindm_formula" \
        "$bindm_target/" | sed -e "s%$bindm_target%$bindm_source/%g"
fi
