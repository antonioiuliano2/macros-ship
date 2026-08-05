#!/bin/bash
ProcId=$2
NEVENTS=1000000
SEED=$((ProcId+5000000))

echo $SEED

sleep $ProcId

SHIP_mymaster=/afs/cern.ch/work/a/aiuliano/public/ECN3_SHIPBuild
export ALIBUILD_WORK_DIR=$SHIP_mymaster/sw #for alienv

echo "SETUP"
source /cvmfs/ship.cern.ch/26.06/setUp.sh
source /afs/cern.ch/work/a/aiuliano/public/sim_ecn3ship/condor_sims/2026_08_05_Cascade_EduardSpreadUpdate/updatecascadeEduard_26_06.env
set -o nounset

OUTPUTPATH=/afs/cern.ch/work/a/aiuliano/public/sim_ecn3ship/condor_sims/2026_08_05_Cascade_EduardSpreadUpdate
python /afs/cern.ch/work/a/aiuliano/public/ECN3_SHIPBuild/FairShip/macro/makeCascade.py --cascade-lambda 12.0 -n $NEVENTS -s $SEED -t $OUTPUTPATH/Run_${ProcId}_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root
