#!/bin/bash
# Evaluate N(gen events)
if [[ "$1" == "" ]]; then
    exit
fi

FILE=`samweb list-files "dh.dataset=${1} and availability:anylocation" | head -n 1`
if [[ "FILE" == "" ]]; then
    echo "File for dataset ${1} not found"
    exit
fi
if [[ "$2" != "" ]]; then
    samweb get-metadata ${FILE}
fi
NFILES=`samweb list-files "dh.dataset=${1} and availability:anylocation" | wc |  awk '{print $1}'`
NGEN=`samweb get-metadata ${FILE} | awk '{if($1 == "dh.gencount:") print $2}'`
if [[ "${NGEN}" == "" ]]; then
    echo "No gen count field found for file ${FILE}"
    exit
fi
echo "N(files)=${NFILES}, N(gen-per-file)=${NGEN}, N(gen)=$((NFILES * NGEN))"
