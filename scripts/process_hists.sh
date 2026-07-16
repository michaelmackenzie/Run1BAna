#! /bin/bash

DATASETS="$1"
RECOVERSION=$2
TAG=$3
NOPLOT=$4

if [[ "${DATASETS}" == "" ]]; then
    DATASETS="RMC DIO COSMIC CE PILEUP NEUTRON PROTON POLY"
fi
if [[ "${RECOVERSION}" == "" ]]; then
    RECOVERSION="v40"
fi
HISTTAG=""
if [[ "${TAG}" != "" ]]; then
    HISTTAG="-${TAG}"
fi

echo "DATASET=${DATASETS}"
echo "RECOVERSION=${RECOVERSION}"
echo "TAG=${TAG}"

# Script + dataset inputs
SCRIPT="Run1BAna/scripts/hist_run1bana_tree.C"
PILEUP="mnbs4b1s51r0002"
RMC="fgam4b1s51r0002"
RPC="rpce4b0s51r0002"
DIO="diob4b1s51r0002"
COSMIC="csms4b0s51r0002"
CE="cele4b1s51r0001"
NEUTRON="neut0b1s51r0003"
PROTON="prot0b1s51r0003"
POLY="pgamcb1s51r0003"

# Version with 2 cm target + 10 cm poly
if [[ "${RECOVERSION}" == "v06" ]]; then
    RMC="fgam6b0s51r0002"
    PILEUP="mnbs6b1s51r0002"
    COSMIC="csms6b0s51r0002"
    RPC="rpce6b0s51r0002"
fi

# Version with 2 cm target + 10 cm poly, early time gate
if [[ "${RECOVERSION}" == "v07" ]]; then
    PILEUP="mnbs7b1s51r0002"
    COSMIC="csms7b0s51r0002"
    RPC="rpce7b1s51r0002"
fi


# Version with 1.75 cm degrader target
if [[ "${RECOVERSION}" == "v40" ]]; then
    RMC="fgam0b1s51r0003"
    PILEUP="mnbs0b1s51r0003"
    COSMIC="csms0b1s51r0003"
    RPC="rpce0b1s51r0003"
    CE="cele0b1s51r0003"
    DIO="fele0b1s51r0003"
fi

# Pileup histogram
if [[ "${DATASETS}" == *"PILEUP"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PILEUP}/nts.mmackenz.${PILEUP}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${PILEUP}${HISTTAG}.hist"
    [ -f ${PILEUP}.log ] && rm ${PILEUP}.log
    # root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")" | tee ${PILEUP}.log
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\", -1)"
fi

# DIO histogram
if [[ "${DATASETS}" == *"DIO"* ]]; then
    # INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${DIO}/nts.mmackenz.${DIO}.Run1BAna.*.root"
    INDATA="nts.owner.${DIO}.Run1BAna.sequencer.root"
    OUTDATA="Run1BAna.${DIO}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# CE histogram
if [[ "${DATASETS}" == *"CE"* ]]; then
    INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${CE}/nts.mmackenz.${CE}.Run1BAna.*.root"
    OUTDATA="Run1BAna.${CE}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Cosmic histogram
if [[ "${DATASETS}" == *"COSMIC"* ]]; then
    if [[ "${RECOVERSION}" == "v07" ]]; then
        INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/csms6b0s51r0002/nts.mmackenz.csms6b0s51r0002.Run1BAna.*.root"
    else
        INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${COSMIC}/nts.mmackenz.${COSMIC}.Run1BAna.*.root"
    fi
    OUTDATA="Run1BAna.${COSMIC}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# RMC histogram
if [[ "${DATASETS}" == *"RMC"* ]]; then
    if [[ "${RECOVERSION}" == "v06" ]]; then
        INDATA="nts.mmackenz.fgam6b0s51r0002.Run1BAna.root"
    else
        INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RMC}/nts.mmackenz.${RMC}.Run1BAna.*.root"
    fi
    OUTDATA="Run1BAna.${RMC}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Neutron histogram
if [[ "${DATASETS}" == *"NEUTRON"* ]] && [[ "${DORPC}" == "" ]]; then
    # INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${NEUTRON}/nts.mmackenz.${NEUTRON}.Run1BAna.*.root"
    INDATA="nts.owner.${NEUTRON}.Run1BAna.sequencer.root"
    OUTDATA="Run1BAna.${NEUTRON}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Proton histogram
if [[ "${DATASETS}" == *"PROTON"* ]] && [[ "${DORPC}" == "" ]]; then
    # INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PROTON}/nts.mmackenz.${PROTON}.Run1BAna.*.root"
    INDATA="nts.owner.${PROTON}.Run1BAna.sequencer.root"
    OUTDATA="Run1BAna.${PROTON}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Poly RMC histogram
if [[ "${DATASETS}" == *"POLY"* ]]; then
    # INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RMC}/nts.mmackenz.${POLY}.Run1BAna.*.root"
    INDATA="nts.owner.${POLY}.Run1BAna.sequencer.root"
    OUTDATA="Run1BAna.${POLY}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# RPC histogram
if [[ "${DATASETS}" == *"RPC"* ]] && [[ "${RECOVERSION}" != "v40" ]]; then
    # INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RPC}/nts.mmackenz.${RPC}.Run1BAna.*.root"
    INDATA="nts.owner.rpce4b0s51r0002.Run1BAna.sequencer.root"
    OUTDATA="Run1BAna.${RPC}${HISTTAG}.hist"
    root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"
fi

# Make plots
if [[ "${NOPLOT}" == "" ]]; then
    if [[ "${DATASETS}" != "CE" ]]; then
        SCRIPT="Run1BAna/scripts/plotRMCvsBkgFromNtuple.C(\"${RECOVERSION}\", \"${TAG}\")"
        root -l -q -b "${SCRIPT}"
    fi
    if [[ "${DATASETS}" != "RMC" ]]; then
        SCRIPT="Run1BAna/scripts/plotCEvsBkgFromNtuple.C(\"${RECOVERSION}\", \"${TAG}\")"
        root -l -q -b "${SCRIPT}"
    fi
fi
