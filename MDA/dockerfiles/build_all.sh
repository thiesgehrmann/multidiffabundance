#!/usr/bin/env bash

for env in mda; do # aldex2 ancombc2 deseq2 limma maaslin2; do
    echo ""
    echo "#######################################################"
    echo "#  BUILDING $env"
    echo "#######################################################"
    docker build --no-cache -f multidiffabundance_${env}.Dockerfile -t thiesgehrmann/multidiffabundance_${env}:1 ../
    docker push thiesgehrmann/multidiffabundance_${env}:1
done