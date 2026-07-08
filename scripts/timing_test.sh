#! /bin/bash

Help() {
    echo "Run a trigger timing test"
    echo "  Arguments:"
    echo "    1: Paths, default = calo_cluster_50,apr_TC_calo"
    echo "    2: Data tag, default = Run1Ban_best_v1_4-000"
    echo "    3: Output tag, default = Run1Ban"
    echo "    4: N(events) to test, default = 20000"
    echo "    5: Optional skip processing step"
}

if [[ "$1" == "-h" ]] || [[ "$1" == "--help" ]]; then
    Help
    exit
fi

PATHS="$1"
DATATAG="$2"
NAME="$3"
NEVENTS="$4"
SKIPPROC="$5"

if [[ "${PATHS}" == "" ]]; then
    PATHS="calo_cluster_50,apr_TC_calo"
fi
if [[ "${DATATAG}" == "" ]]; then
    DATATAG="Run1Ban_best_v1_4-000"
fi
if [[ "${NAME}" == "" ]]; then
    NAME="Run1Ban"
fi
if [[ "${NEVENTS}" == "" ]]; then
    NEVENTS="20000"
fi

TAG="1BB"
DIR="timing/${NAME}"
DATAFILE="Run1BAna/file_lists/dig.mu2e.NoPrimaryMix1BB.${DATATAG}.art.files"

if [ ! -f ${DATAFILE} ]; then
    echo "Data file list ${DATAFILE} not found!"
    Help
    exit 1
fi

[ ! -d ${DIR} ] && mkdir -p ${DIR}

if [[ "${SKIPPROC}" == "" ]]; then
    python mu2e-trig-config/python/genTimingFcl.py -n ${NAME} -p ${PATHS}
    if [ ! -f ${NAME}.fcl ]; then
        echo "Failed to generate timing fcl!"
        exit 1
    fi

    # Set the services for Run 1B
    echo 'services.DbService.purpose: "Sim_best"' >> ${NAME}.fcl
    echo 'services.DbService.version: "v1_5"' >> ${NAME}.fcl
    echo 'services.DbService.nearestMatch: "true"' >> ${NAME}.fcl

    mu2e -c ${NAME}.fcl --debug-config ${NAME}_debug.fcl
    if [ ! -f ${NAME}_debug.fcl ]; then
        echo "Failed to generate debug fcl!"
        exit 1
    fi
    cp ${NAME}.fcl ${DIR}/
    cp ${NAME}_debug.fcl ${DIR}/
    /usr/bin/time -v mu2e -c ${NAME}.fcl -S ${DATAFILE} -n ${NEVENTS}
    RC=$?
    if [ $RC -ne 0 ]; then
        echo "art exited with return code ${RC} --> exiting!"
        exit 1
    fi
    if [ ! -f triggerTiming.db ]; then
        echo "Failed to generate timing database!"
        exit 1
    fi
    cp triggerTiming.db ${DIR}/triggerTiming.db
fi
source mu2eTimingPlotsMaker/bash/ProcessSQL.sh ${DIR} ${TAG}
python mu2eTimingPlotsMaker/python/merge_timing_files.py -i ${DIR}/csv_${TAG}/
root.exe -q -b "mu2eTimingPlotsMaker/scripts/runAllTiming.C(\"${DIR}/csv_${TAG}\")"
python mu2eTimingPlotsMaker/python/plot_filters_time_root.py -s ${DIR}/csv_${TAG}/timing_plots.root -o ${DIR}/csv_${TAG}/plots/
