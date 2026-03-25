#! /bin/bash

DATASETS="$1"
TAG=$2
NOPLOT=$3

if [[ "${DATASETS}" == "" ]]; then
    DATASETS="RMC DIO RPC COSMIC CE PILEUP"
fi
if [[ "${TAG}" == "" ]]; then
    TAG="v04"
fi

# Script + dataset inputs
SCRIPT="Run1BAna/scripts/hist_run1bana_tree.C"
PILEUP="mnbs4b1s51r0002"
RMC="fgam4b1s51r0002"
RPC="rpce4b0s51r0002"
DIO="diob4b1s51r0002"
COSMIC="csms4b0s51r0002"
CE="cele4b1s51r0001"

# Version with 2 cm target + 10 cm poly
if [[ "${TAG}" == "v06" ]]; then
    RMC="fgam6b0s51r0002"
    PILEUP="mnbs6b1s51r0002"
    COSMIC="csms6b0s51r0002"
fi

# Pileup histogram
if [[ "${DATASETS}" == *"PILEUP"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PILEUP}/nts.mmackenz.${PILEUP}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${PILEUP}.hist"
    [ -f ${PILEUP}.log ] && rm ${PILEUP}.log
    # root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")" | tee ${PILEUP}.log
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", -1)"
fi

# DIO histogram
if [[ "${DATASETS}" == *"DIO"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${DIO}/nts.mmackenz.${DIO}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${DIO}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# CE histogram
if [[ "${DATASETS}" == *"CE"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${CE}/nts.mmackenz.${CE}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${CE}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Cosmic histogram
if [[ "${DATASETS}" == *"COSMIC"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${COSMIC}/nts.mmackenz.${COSMIC}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${COSMIC}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# RMC histogram
if [[ "${DATASETS}" == *"RMC"* ]]; then
    if [[ "${TAG}" == "v06" ]]; then
        INDATA="nts.mmackenz.fgam6b0s51r0002.Run1BAna.root"
    else
        INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RMC}/nts.mmackenz.${RMC}.Run1BAna.*.root"
    fi
    OUTDATA="Run1BAna.${RMC}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# RPC histogram
if [[ "${DATASETS}" == *"RPC"* ]]; then
    # INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RPC}/nts.mmackenz.${RPC}.Run1BAna.*.root"
    INDATA="nts.owner.rpce4b0s51r0002.Run1BAna.sequencer.root"
    OUTDATA="Run1BAna.${RPC}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Make plots
if [[ "${NOPLOT}" == "" ]]; then
    SCRIPT="Run1BAna/scripts/plotRMCvsBkgFromNtuple.C(\"${TAG}\")"
    root -l -q -b "${SCRIPT}"
fi
