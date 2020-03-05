#!/bin/bash
#BSUB -q 1nd
ProcId=$2
NEVENTS=2000 #approximatively npot per spill
STARTEVENT=$((ProcId * $NEVENTS))
echo "SETUP"

SHIPBUILD_mymaster=/afs/cern.ch/work/a/aiuliano/public/SHIPBuild
export ALIBUILD_WORK_DIR=$SHIPBUILD_mymaster/sw #for alienv
source /cvmfs/ship.cern.ch/SHiP-2020/latest/setUp.sh
eval `alienv load FairShip/latest`

INPUTFILE=/afs/cern.ch/work/a/aiuliano/public/Cascadefiles/Cascade_500k_Lead_withprimaries.root

OUTPUTDIR=/afs/cern.ch/work/a/aiuliano/public/sim_charm/condor_sims/CH2_charmcascade_14_02_20

echo "Starting script from ", $STARTEVENT

python /afs/cern.ch/work/a/aiuliano/public/SHIPBuild/FairShip/macro/run_simScript.py --charm 1 -A charmonly --CharmdetSetup 1 -f $INPUTFILE --CharmTarget 2 --TrackingCharm -n $NEVENTS -i $STARTEVENT -o  $OUTPUTDIR/simulation/$ProcId




