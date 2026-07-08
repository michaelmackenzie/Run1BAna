#! /bin/bash

DATASETS="$1"
DORPC=$2
NOPLOT=$3

if [[ "${DATASETS}" == "" ]]; then
    DATASETS="CE PILEUP"
fi

# Script + dataset inputs
SCRIPT="Run1BAna/scripts/hist_run1bana_tree.C"
PILEUP="mnbs0b1s51r0003"
DIO="diob0b1s51r0004"
COSMIC="csms0b1s51r0003"
CE="cele0b1s51r0003"
TAG="v40"


# Pileup histogram
if [[ "${DATASETS}" == *"PILEUP"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PILEUP}/nts.mmackenz.${PILEUP}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${PILEUP}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", -1)"
fi

# DIO histogram
if [[ "${DATASETS}" == *"DIO"* ]] && [[ "${DORPC}" == "" ]] && [[ "${TAG}" != "v40" ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${DIO}/nts.mmackenz.diobb1s51r0000.Run1BAna.*.root"
    OUTDATA="Run1BAna.${DIO}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# CE histogram
if [[ "${DATASETS}" == *"CE"* ]] && [[ "${DORPC}" == "" ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${CE}/nts.mmackenz.${CE}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${CE}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Cosmic histogram
if [[ "${DATASETS}" == *"COSMIC"* ]] && [[ "${DORPC}" == "" ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${COSMIC}/nts.mmackenz.${COSMIC}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${COSMIC}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Make plots
if [[ "${NOPLOT}" == "" ]]; then
    SCRIPT="Run1BAna/scripts/plotCEvsBkgFromNtuple.C(\"${TAG}\")"
    root -l -q -b "${SCRIPT}"
fi
