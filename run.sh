#!/bin/sh
module purge
module load gcc/8.3.1 hpcx/2.13.0-mt

# IPM
USE_IPM=0
if [[ $USE_IPM -eq 1 ]]; then
    export IPM_KEYFILE=$HPCX_IPM_DIR/etc/ipm_key_mpi                             
    export IPM_REPORT=full                                           
    export IPM_LOG=full                                              
    export IPM_STATS=all                                             
    export LD_PRELOAD=$HPCX_IPM_LIB
fi

# Score-P
USE_SCOREP=0
if [[ $USE_SCOREP -eq 1 ]]; then
    module load ectrans/1.2.0-scorep
    export SCOREP_TIMER=clock_gettime
    export SCOREP_ENABLE_PROFILING=true
    export SCOREP_ENABLE_TRACING=true
    export SCOREP_EXPERIMENT_DIRECTORY=scorep
    export SCOREP_PROFILING_MAX_CALLPATH_DEPTH=2048
    #export SCOREP_FILTERING_FILE=$SLURM_SUBMIT_DIR/filter
    export SCOREP_TOTAL_MEMORY=64MB
    #export SCOREP_METRIC_PAPI=PAPI_L1_TCA,PAPI_L1_TCM,PAPI_L2_TCA,PAPI_L2_TCM,PAPI_L3_TCA,PAPI_L3_TCM,PAPI_MEM_RCY,PAPI_MEM_SCY,PAPI_MEM_WCY,PAPI_PRF_DM,PAPI_RES_STL
else
    module load ectrans/1.2.0
fi

# ECTRANS
SMALL=1279
MEDIUM=2559
LARGE=3839
CHECK=1
ITER=1
LEVEL=137
NPROMA=32
PRECISION="sp"
#MEMINFO="--meminfo"

SLURM_JOB_CPUS_PER_NODE=32
if [[ -n SLURM_JOB_NUM_NODES ]]; then
    N=$SLURM_JOB_NUM_NODES
else
    N=1
fi
if [[ -n $SLURM_NTASKS_PER_NODE ]]; then
    PPN=$SLURM_NTASKS_PER_NODE
elif [[ -n SLURM_JOB_CPUS_PER_NODE ]]; then
    PPN=$SLURM_JOB_CPUS_PER_NODE
else
    PPN=2
fi
NP=$((N*PPN))
MPIRUN_CMD="\
    mpirun \
        -np $NP \
        -map-by ppr:$((PPN/2)):socket \
        -mca pml ucx \
            -x UCX_NET_DEVICES=mlx5_0:1 \
            -x UCX_LOG_LEVEL=error \
        ectrans-benchmark-$PRECISION -t $SMALL -c 100 -n $ITER -l $LEVEL --vordiv --scders --uvders --nproma $NPROMA --norms $MEMINFO
"

# Report
LOGFILE="run.log"
echo "Benchmark: ECTRANS" | tee $LOGFILE
echo "Partition: $SLURM_JOB_PARTITION" | tee -a $LOGFILE
echo "Nodelist : $SLURM_JOB_NODELIST" | tee -a $LOGFILE
echo "Node #   : $SLURM_JOB_NUM_NODES" | tee -a $LOGFILE
echo "PPN      : $SLURM_JOB_CPUS_PER_NODE" | tee -a $LOGFILE
echo "NP       : $((SLURM_JOB_NUM_NODES * SLURM_JOB_CPUS_PER_NODE))" | tee -a $LOGFILE
echo "MPIRUN   : "$(echo ${MPIRUN_CMD} | tr -s ' ') | tee -a $LOGFILE

# RUN
$MPIRUN_CMD 2>&1 | tee -a $LOGFILE
