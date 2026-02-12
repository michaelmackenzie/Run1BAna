#! /bin/bash
# Process mcs -> histogram files


DATASETS=$1
MAXEVENTS=$2
DRYRUN=$3

if [[ "${DATASETS}" == "" ]]; then
    #v01: Required a line fit on output, pileup generation was wrong
    # DATASETS="mcs.mu2e.CeEndpointMixLowTriggerable-KL.Run1Bab2_best_v1_2.art \
    #         mcs.mu2e.DIOtail90MixLowTriggerable-KL.Run1Bab2_best_v1_2.art \
    #         mcs.mu2e.DIOtail80MixLowTriggerable-KL.Run1Bab2_best_v1_2.art \
    #         mcs.mu2e.DIOtail50MixLowTriggerable-KL.Run1Bab2_best_v1_2.art \
    #         mcs.mu2e.NoPrimaryMixLowTriggerable-KL.Run1Bab2_best_v1_2.art"

    #v02: No line fit requirement, 1% 1BB
    DATASETS="mcs.mu2e.CeEndpointMixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
              mcs.mu2e.DIOtail90_infMixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
              mcs.mu2e.DIOtail80_90MixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
              mcs.mu2e.DIOtail60_80MixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
              mcs.mu2e.DIOtail0_60MixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
              mcs.mu2e.FlatGammaMixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
              mcs.mu2e.FlateMinusMixLowTriggerable-KL.Run1Bab3_best_v1_2.art \
              mcs.mu2e.NoPrimaryMixLowTriggerable-KL.Run1Bab3_best_v1_2.art"
fi

if [[ "${MAXEVENTS}" == "" ]]; then
    MAXEVENTS="100000"
fi


for DATASET in ${DATASETS}; do
    OUTNAME=`echo ${DATASET} | sed 's/mcs./nts./g' | sed 's/.art/.root/g'`
    echo "${DATASET} --> ${OUTNAME}"
    FILELIST="Run1BAna/file_lists/${DATASET}.files"
    if [ ! -f ${FILELIST} ]; then
        echo "Can't find file list ${FILELIST}!"
        continue
    fi
    FCL="ana.fcl"
    if [[ "${DATASET}" == *"DIOtail50"* ]]; then
        FCL="ana_dio_50.fcl"
    elif [[ "${DATASET}" == *"DIOtail80"* ]]; then
        FCL="ana_dio_80.fcl"
    fi
    echo mu2e -c "Run1BAna/fcl/${FCL}" -n ${MAXEVENTS} -S ${FILELIST} -T ${OUTNAME}
    if [[ "${DRYRUN}" == "" ]]; then
        mu2e -c "Run1BAna/fcl/${FCL}" -n ${MAXEVENTS} -S ${FILELIST} -T ${OUTNAME}
        RESULT=$?
        if [ ${RESULT} -ne 0 ]; then
            echo "mu2e job exited with code ${RESULT}!"
            exit 1
        fi
        if [ ! -f ${OUTNAME} ]; then
            echo "Failed to create ${OUTNAME}!"
        fi
    fi
done

echo "Finished processing"
