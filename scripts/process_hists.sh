#! /bin/bash

DATASETS="$1"
RECOVERSION=$2
TAG=$3
NMAX=$4
DRYRUN=$5

if [[ "${DATASETS}" == "" ]]; then
    DATASETS="RMC RPC DIO COSMIC CE PILEUP NEUTRON PROTON POLY"
fi
if [[ "${RECOVERSION}" == "" ]]; then
    RECOVERSION="v40"
fi
HISTTAG=""
if [[ "${TAG}" != "" ]]; then
    HISTTAG="-${TAG}"
fi
if [[ "${NMAX}" == "" ]]; then
    NMAX="2e6"
fi

echo "DATASET=${DATASETS}"
echo "RECOVERSION=${RECOVERSION}"
echo "TAG=${TAG}"

# Script + dataset inputs
SCRIPT="Run1BAna/scripts/hist_run1bana_tree_v2.C"

# Version with 1.75 cm degrader target
if [[ "${RECOVERSION}" == "v40" ]]; then
    RMC="fgam0b1s51r0004"
    # PILEUP="mnbs0b1s51r0004" # Unfiltered
    PILEUP="mnbs1b1s51r0004" # Filtered
    COSMIC="csms0b1s51r0004"
    RPC="rpce0b1s51r0004"
    CE="cele0b1s51r0004"
    DIO="fele0b1s51r0004"
    NEUTRON="neut0b1s51r0004"
    PROTON="prot0b1s51r0004"
    POLY="pgamcb1s51r0004"
else
    echo "Unknown reco version ${RECOVERSION}"
    exit 1
fi

HEAD=""
if [[ "${DRYRUN}" != "" ]]; then
    HEAD="echo"
    echo "Performing a dry run!"
fi

# Pileup histogram
if [[ "${DATASETS}" == *"PILEUP"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PILEUP}/nts.mmackenz.${PILEUP}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${PILEUP}${HISTTAG}.hist"
    # Always process all of the pileup due to low acceptance
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", -1)"
fi

# DIO histogram
if [[ "${DATASETS}" == *"DIO"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${DIO}/nts.mmackenz.${DIO}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${DIO}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi

# CE histogram
if [[ "${DATASETS}" == *"CE"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${CE}/nts.mmackenz.${CE}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${CE}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi

# Cosmic histogram
if [[ "${DATASETS}" == *"COSMIC"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${COSMIC}/nts.mmackenz.${COSMIC}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${COSMIC}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi

# RMC histogram
if [[ "${DATASETS}" == *"RMC"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RMC}/nts.mmackenz.${RMC}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${RMC}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi

# Neutron histogram
if [[ "${DATASETS}" == *"NEUTRON"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${NEUTRON}/nts.mmackenz.${NEUTRON}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${NEUTRON}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi

# Proton histogram
if [[ "${DATASETS}" == *"PROTON"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PROTON}/nts.mmackenz.${PROTON}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${PROTON}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi

# Poly RMC histogram
if [[ "${DATASETS}" == *"POLY"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${POLY}/nts.mmackenz.${POLY}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${POLY}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi

# RPC histogram
if [[ "${DATASETS}" == *"RPC"* ]]; then
    INDATA="nts.owner.${RPC}.Run1BAna.sequencer.root"
    OUTDATA="Run1BAna.${RPC}${HISTTAG}.hist"
    ${HEAD} root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", ${NMAX})"
fi
