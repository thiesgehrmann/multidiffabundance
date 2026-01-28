#!/usr/bin/env bash

for env in ancombc2; do # aldex2 aldex3 ancombc2 deseq2 limma maaslin2 maaslin3 mda; do
    echo ""
    echo "#######################################################"
    echo "#  BUILDING $env"
    echo "#######################################################"
    docker build --no-cache --progress=plain -f multidiffabundance_${env}.Dockerfile -t thiesgehrmann/multidiffabundance_${env}:1 ../
    docker push thiesgehrmann/multidiffabundance_${env}:1
done
