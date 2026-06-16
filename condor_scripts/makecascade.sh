#!/bin/bash
ProcId=$2
NEVENTS=1000000
SEED=$((ProcId+5000000))

echo $SEED

sleep $ProcId

SNDLHC_mymaster=/afs/cern.ch/work/a/aiulian/public/SNDLHCBuild
export ALIBUILD_WORK_DIR=$SNDLHC_mymaster/sw #for alienv

echo "SETUP"
source /cvmfs/sndlhc.cern.ch/SNDLHC-2025/Jan30/setUp.sh
eval `alienv load sndsw/latest`

OUTPUTPATH=/afs/cern.ch/work/a/aiulian/public/sim_snd/cascadeproduction_Wtarget_2026
python /afs/cern.ch/work/a/aiulian/public/SHiP/FairShip/macro/makeCascade.py -n $NEVENTS -s $SEED -t $OUTPUTPATH/Run_${ProcId}_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root
