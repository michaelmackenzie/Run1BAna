# Copy over ntuple job results
#! /bin/bash

Help() {
    echo "Copy over ntuple job outputs"
    echo "Usage: copy_ntuples.sh <INDATA> <OUTDATA> <JOBID> <JOB)> <STAGE> <PROJECT>"
}

#---------------------------------------------------------------------------
# Get the inputs
#---------------------------------------------------------------------------

INDATA=$1
OUTDATA=$2
JOBID=$3
JOB=$4
STAGE=$5
PROJECT=$6

if [[ "${INDATA}" == "" ]] || [[ "${INDATA}" == "-h" ]] || [[ "${INDATA}" == "--help" ]]; then
    Help
    exit 1
fi

if [[ "${OUTDATA}" == "" ]] || [[ "${JOBID}" == "" ]] || [[ "${JOB}" == "" ]] || [[ "${STAGE}" == "" ]] || [[ "${PROJECT}" == "" ]]; then
    echo "Missing inputs!"
    Help
    exit 1
fi


#---------------------------------------------------------------------------
# Ensure permissions/access are set up
#---------------------------------------------------------------------------

if ! klist -s; then
    echo "Kerberos ticket is not valid!"
    exit 2
fi
kinit -R

getToken
if [ $? -ne 0 ]; then
    echo "Error getting token!"
    exit 2
fi

#---------------------------------------------------------------------------
# Copy over the ntuples/data
#---------------------------------------------------------------------------

OUTDIR="/exp/mu2e/data/users/${USER}/run1b/data/${OUTDATA}/"
[ ! -d ${OUTDIR} ] && mkdir -p ${OUTDIR}

INDIR=`ls -d /pnfs/mu2e/scratch/users/${USER}/workflow/${PROJECT}.${INDATA}.${STAGE}_${JOB}/outstage/${JOBID}/ 2>/dev/null | head -n 1`
if [[ "${INDIR}" == "" ]] || [ ! -d ${INDIR} ]; then
    echo "No input directory found!"
    exit 3
fi

# Copy over the data
setup ifdhc
for dir in `ls -d ${INDIR}*/*/`; do
    LOGFILE=`ls -d ${dir}/*.log 2>/dev/null | head -n 1`
    if [[ "${LOGFILE}" == "" ]] || [ ! -f $LOGFILE ]; then
        echo ">>> Directory ${dir} has no log file!"
        continue
    fi
    # Check for errors
    RETURNLINE=`grep "Art has completed" ${LOGFILE} | head -n 1`
    if [[ "${RETURNLINE}" == "" ]]; then
        echo ">>> No art return code: ${dir}"
        continue;
    fi
    if [[ "${RETURNLINE}" != *"status 0"* ]]; then
        echo ">>> Error: ${RETURNLINE}: ${dir}"
        FILEOPENERROR=`grep "FileOpenError" ${LOGFILE} | head -n 1`
        if [[ "${FILEOPENERROR}" != "" ]]; then
            echo "--> Continuing processing given a file open error"
        else
            continue;
        fi
    fi
    # Copy over data
    for f in `ls -d ${dir}*.root* ${dir}jsn*.json 2>/dev/null`; do
        bn=`basename $f`
        if [ -f ${OUTDIR}${bn} ]; then
            continue
        fi
        echo "ifdh cp ${f} ${OUTDIR}"
        ifdh cp ${f} ${OUTDIR}
    done
done

echo "Finished copying over files"
