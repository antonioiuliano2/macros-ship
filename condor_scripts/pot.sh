#!/bin/bash
#BSUB -q 1nd
ProcId=$2
NEVENTS=13500 #approximatively npot per spill
echo "SETUP"

SHIPBUILD_mymaster=/afs/cern.ch/work/a/aiuliano/public/SHIPBuild
export ALIBUILD_WORK_DIR=$SHIPBUILD_mymaster/sw #for alienv
source /cvmfs/ship.cern.ch/SHiP-2020/latest/setUp.sh
eval `alienv load FairShip/latest`

echo "Starting script from ", $STARTEVENT

OUTPUTFOLDER=/afs/cern.ch/work/a/aiuliano/public/sim_charm/condor_sims/CH2_pot_10_02_20/simulation

python $FAIRSHIP/muonShieldOptimization/run_MufluxfixedTarget.py --CharmdetSetup 1 --CharmTarget 2 -G -e 0.1 -n $NEVENTS -o $OUTPUTFOLDER/$ProcId
