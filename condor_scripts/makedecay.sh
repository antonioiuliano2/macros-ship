#!/bin/bash
ProcId=$2
#normalization to 5e+13 with 1000 runs
POTNORM=5e+10

SHIP_mymaster=/afs/cern.ch/work/a/aiuliano/public/ECN3_SHIPBuild
export ALIBUILD_WORK_DIR=$SHIP_mymaster/sw #for alienv

echo "SETUP"
source /cvmfs/ship.cern.ch/26.06/setUp.sh
source /afs/cern.ch/work/a/aiuliano/public/sim_ecn3ship/condor_sims/2026_08_05_Decay_EduardSpreadUpdate/updatecascadeEduard_26_06.env
set -o nounset

INPUTPATH=/eos/experiment/ship/user/aiuliano/2026_08_05_Cascade_EduardSpreadUpdate
OUTPUTPATH=/afs/cern.ch/work/a/aiuliano/public/sim_ecn3ship/condor_sims/2026_08_05_Decay_EduardSpreadUpdate

cd $OUTPUTPATH
python /afs/cern.ch/work/a/aiuliano/public/ECN3_SHIPBuild/FairShip/macro/makeDecay.py -f $INPUTPATH/Run_${ProcId}_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root -p $POTNORM
