#!/usr/bin/env bash

for env in maaslin3; do # aldex2 aldex3 ancombc2 corncob deseq2 limma maaslin2 maaslin3 zicoseq mda; do
    echo ""
    echo "#######################################################"
    echo "#  BUILDING $env"
    echo "#######################################################"
    docker build --no-cache --progress=plain -f multidiffabundance_${env}.Dockerfile -t thiesgehrmann/multidiffabundance_${env}:1 /home/thies/repos/devel/multidiffabundance
    docker push thiesgehrmann/multidiffabundance_${env}:1
done
