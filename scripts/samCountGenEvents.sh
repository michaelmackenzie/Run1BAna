#!/bin/bash
# Evaluate N(gen events)
if [[ "$1" == "" ]]; then
    exit
fi

NTESTS=$3
if [[ "${NTESTS}" == "" ]]; then
    NTESTS=10
fi

setup dhtools
NFILES=`samweb list-files "dh.dataset=${1} and availability:anylocation" | wc -l`
if [[ "${NFILES}" == "" ]] || [[ ${NFILES} -eq 0 ]]; then
    echo "No dataset files found!"
    exit
fi

# Ensure there are enough files for the tests requested
if [ ${NTESTS} -gt ${NFILES} ]; then
    echo "Only ${NFILES} in the dataset --> Reducing to this many files to check"
    NTESTS=${NFILES}
fi

NGEN=0
for ((i = 1; i <= NTESTS; i++))
do
    FILE=`samweb list-files "dh.dataset=${1} and availability:anylocation" | head -n $i | tail -n 1`
    if [[ "FILE" == "" ]]; then
        echo "File for dataset ${1} not found"
        exit
    fi
    if [[ "$2" != "" ]]; then
        samweb get-metadata ${FILE}
    fi
    NFILES=`samweb list-files "dh.dataset=${1} and availability:anylocation" | wc |  awk '{print $1}'`
    NGEN_I=`samweb get-metadata ${FILE} | awk '{if($1 == "dh.gencount:") print $2}'`
    if [[ "${NGEN_I}" == "" ]]; then
        echo "No gen count field found for file ${FILE}"
        exit
    fi
    NGEN=$((NGEN + NGEN_I))
done

echo "N(files)=${NFILES}, N(gen per $NTESTS file(s))=${NGEN}, N(gen)=$((NFILES * NGEN / NTESTS))"
