# Create a new dataset directory

Help() {
    echo "Create a new dataset directory"
    echo "Usage: ./init_dataset.sh <name> <reference> [optional \"d\" for dryrun]"
    echo "Example: ./init_dataset cele0b2 datasets/cele0b1"
}

DATASET=$1
REFERENCE=$2
DRYRUN=$3

if [[ "${DATASET}" == "" ]] || [[ "${REFERENCE}" == "" ]]; then
    echo "Missing inputs!"
    Help
    exit
fi

if [ ! -d ${REFERENCE} ]; then
    echo "Reference ${REFERENCE} not found"
    Help
    exit
fi

DIR=`dirname ${REFERENCE}`
if [ -d ${DIR}/${DATASET} ]; then
    echo "Dataset ${DATASET} already exists, remove before recreating"
    Help
    exit
fi

HEAD=""
if [[ "${DRYRUN}" != "" ]]; then
    HEAD="echo"
fi

BASE=`basename ${REFERENCE}`
${HEAD} cp -r ${REFERENCE} ${DIR}/${DATASET}

${HEAD} sed -i "s/${BASE}/${DATASET}/g" ${DIR}/${DATASET}/*.*
${HEAD} rename ${BASE} ${DATASET} ${DIR}/${DATASET}/*.*
${HEAD} rename ${BASE} ${DATASET} ${DIR}/${DATASET}/catalog/*.files
