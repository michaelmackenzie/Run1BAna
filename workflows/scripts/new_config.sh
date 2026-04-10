# Create a new configuration, starting from config_v03 as a reference
#! /bin/bash

#----------------------------------------------------------------------------------------
Help() {
    echo "Create a new configuration area"
    echo "new_config.sh <config name, e.g. config_v10> <run number, e.g. 1480> <reference, default = config_v03> <reference run, default = 1450>"
}


#----------------------------------------------------------------------------------------
# Get and validate the inputs
#----------------------------------------------------------------------------------------
CONFIG=$1
RUN=$2
REFERENCE=$3
REFRUN=$4


if [[ "$1" == "" ]] || [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    Help
    exit
fi

if [[ "${RUN}" == "" ]]; then
    echo "No run number given!"
    Help
    exit 1
fi

if ! [[ "${RUN}" =~ ^[0-9]+$ ]]; then
    echo "Run number ${RUN} is not a number!"
    Help
    exit 1
fi

if [ -d ${CONFIG} ]; then
    echo "Config directory ${CONFIG} already exists!"
    Help
    exit 1
fi

if [[ "${REFERENCE}" == "" ]]; then
    REFERENCE="config_v03"
    REFRUN=1450
elif [[ "${REFRUN}" == "" ]]; then
    echo "No run number given for reference ${REFERENCE}"
    Help
    exit 1
elif ! [[ "${REFRUN}" =~ ^[0-9]+$ ]]; then
    echo "Reference run number ${REFRUN} is not a number!"
    Help
    exit 1
fi

#----------------------------------------------------------------------------------------
# Build the new area
#----------------------------------------------------------------------------------------

# Offset by 100 for Run 1A
RUNA=$((RUN + 100))
REFRUNA=$((REFRUN + 100))

mkdir ${CONFIG}
cp -r ${REFERENCE}/* ${CONFIG}/
sed -i "s/${REFRUN}/${RUN}/g" ${CONFIG}/run1b_beam/epilog_run.fcl
sed -i "s/${REFRUNA}/${RUNA}/g" ${CONFIG}/run1a_beam/epilog_run.fcl
sed -i "s/${REFERENCE}/${CONFIG}/g" ${CONFIG}/*.fcl
sed -i "s/${REFERENCE}/${CONFIG}/g" ${CONFIG}/*/*.fcl
sed -i "s/${REFERENCE}/${CONFIG}/g" ${CONFIG}/*/*.txt

echo "Built config area ${CONFIG}, please update ${CONFIG}/run1b_beam/geom.txt as needed"
