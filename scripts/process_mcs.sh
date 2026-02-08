#! /bin/bash
# Process mcs -> histogram files


DATASETS=$1
MAXEVENTS=$2

if [[ "${DATASETS}" == "" ]]; then
    DATASETS="mcs.mu2e.CeEndpointMixLowTriggerable-KL.Run1Bab2_best_v1_2.art mcs.mu2e.DIOtail90MixLowTriggerable-KL.Run1Bab2_best_v1_2.art mcs.mu2e.DIOtail80MixLowTriggerable-KL.Run1Bab2_best_v1_2.art mcs.mu2e.DIOtail50MixLowTriggerable-KL.Run1Bab2_best_v1_2.art mcs.mu2e.NoPrimaryMixLowTriggerable-KL.Run1Bab2_best_v1_2.art"
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
         mu2e -c "Run1BAna/fcl/${FCL}" -n ${MAXEVENTS} -S ${FILELIST} -T ${OUTNAME}
done

echo "Finished processing"
