#! /bin/bash

DATASETS="mcs.mu2e.CeEndpointMixLowTriggerable-KL.Run1Bab3_best_v1_2.art    \
          mcs.mu2e.DIOtail90_infMixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
          mcs.mu2e.DIOtail80_90MixLowTriggerable-KL.Run1Bab3_best_v1_2.art  \
          mcs.mu2e.DIOtail60_80MixLowTriggerable-KL.Run1Bab3_best_v1_2.art  \
          mcs.mu2e.DIOtail0_60MixLowTriggerable-KL.Run1Bab3_best_v1_2.art   \
          mcs.mu2e.FlatGammaMixLowTriggerable-KL.Run1Bab3_best_v1_2.art     \
          mcs.mu2e.FlateMinusMixLowTriggerable-KL.Run1Bab3_best_v1_2.art    \
          mcs.mu2e.NoPrimaryMixLowTriggerable-KL.Run1Bab3_best_v1_2.art"

setup dhtools; setup mu2efiletools
for DATASET in ${DATASETS}; do
    echo "${DATASET}:"
    ./Run1BAna/scripts/samCountGenEvents.sh ${DATASET}
    NEVENTS=`./Run1BAna/scripts/samCountEvents.sh ${DATASET}`
    echo "N(events in dataset)=${NEVENTS}"
done
