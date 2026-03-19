#! /bin/bash

DATASETS="$1"
DORPC=$2
NOPLOT=$3

if [[ "${DATASETS}" == "" ]]; then
    DATASETS="RMC DIO RPC PILEUP"
fi

# Script + dataset inputs
SCRIPT="Run1BAna/scripts/hist_run1bana_tree.C"
PILEUP="mnbs4b1s51r0000"
RMC="fgam4b1s51r0000"
RPC="rpce4b0s51r0001"
DIO="diob4b1s51r0000"
TAG="v04"

if [[ "${DORPC}" != "" ]]; then
    echo "Performing RPC processing"
    PILEUP="mnbs4b1s51r0001"
    TAG="v05"
fi

# Pileup histogram
if [[ "${DATASETS}" == *"PILEUP"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PILEUP}/nts.mmackenz.${PILEUP}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${PILEUP}.hist"
    [ -f ${PILEUP}.log ] && rm ${PILEUP}.log
    # root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")" | tee ${PILEUP}.log
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# DIO histogram
if [[ "${DATASETS}" == *"DIO"* ]] && [[ "${DORPC}" == "" ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${DIO}/nts.mmackenz.diobb1s51r0000.Run1BAna.*.root"
    OUTDATA="Run1BAna.${DIO}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# RMC histogram
if [[ "${DATASETS}" == *"RMC"* ]] && [[ "${DORPC}" == "" ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RMC}/nts.mmackenz.${RMC}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${RMC}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# RPC histogram
if [[ "${DATASETS}" == *"RPC"* ]] && [[ "${DORPC}" != "" ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RPC}/nts.mmackenz.${RPC}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${RPC}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Make plots
if [[ "${NOPLOT}" == "" ]]; then
    if  [[ "${DORPC}" == "" ]]; then
        SCRIPT="Run1BAna/scripts/plotRMCvsBkgFromNtuple.C(\"Run1BAna.${RMC}.hist\", \"Run1BAna.${PILEUP}.hist\", \"${TAG}\")"
        root -l -q -b "${SCRIPT}"
    else
        SCRIPT="Run1BAna/scripts/plotRPCvsBkgFromNtuple.C(\"${TAG}\")"
        root -l -q -b "${SCRIPT}"
    fi
fi
