#! /bin/bash

# Script + dataset inputs
SCRIPT="Run1BAna/scripts/hist_run1bana_tree.C"
PILEUP="mnbs4b1s51r0000"
RMC="fgam4b1s51r0000"
DIO="diob4b1s51r0000"
TAG="v04"

# Pileup histogram
INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${PILEUP}/nts.mmackenz.${PILEUP}.Run1BAna.*.root"
OUTDATA="Run1BAna.${PILEUP}.hist"
[ -f ${PILEUP}.log ] && rm ${PILEUP}.log
root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")" | tee ${PILEUP}.log

# DIO histogram
INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${DIO}/nts.mmackenz.diobb1s51r0000.Run1BAna.*.root"
OUTDATA="Run1BAna.${DIO}.hist"
root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"

# Signal histogram
INDATA="/exp/mu2e/data/users/mmackenz/run1b/data/${RMC}/nts.mmackenz.${RMC}.Run1BAna.*.root"
OUTDATA="Run1BAna.${RMC}.hist"
root -l -q -b "${SCRIPT}(\"${INDATA}\", \"${OUTDATA}\")"

# Make plots
SCRIPT="Run1BAna/scripts/plotRMCvsBkgFromNtuple.C(\"Run1BAna.${RMC}.hist\", \"Run1BAna.${PILEUP}.hist\", \"${TAG}\")"
root -l -q -b "${SCRIPT}"
