#!/bin/bash
ProcId=$2
#normalization to 5e+13 with all 1000 runs
POTNORM=5e+10

SNDLHC_mymaster=/afs/cern.ch/work/a/aiulian/public/SNDLHCBuild
export ALIBUILD_WORK_DIR=$SNDLHC_mymaster/sw #for alienv

echo "SETUP"
source /cvmfs/sndlhc.cern.ch/SNDLHC-2025/Jan30/setUp.sh
eval `alienv load sndsw/latest`

INPUTPATH=/afs/cern.ch/work/a/aiulian/public/sim_snd/cascadeproduction_Wtarget_2026
OUTPUTPATH=/afs/cern.ch/work/a/aiulian/public/sim_snd/decayproduction_Wtarget_2026_Ecut

cd $OUTPUTPATH
python /afs/cern.ch/work/a/aiulian/public/SHiP/FairShip/macro/makeDecay.py -f $INPUTPATH/Run_${ProcId}_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root -p $POTNORM
