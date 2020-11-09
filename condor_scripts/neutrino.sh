#!/bin/bash
#launched with deepcopy option in run_simScript.py
ProcId=$2
LSB_JOBINDEX=$((ProcId+1))
echo $LSB_JOBINDEX

SHIPBUILD_mymaster=/afs/cern.ch/work/a/aiuliano/public/SHIPBuild
export ALIBUILD_WORK_DIR=$SHIPBUILD_mymaster/sw #for alienv

echo "SETUP"
source /cvmfs/ship.cern.ch/SHiP-2020/latest/setUp.sh
eval `alienv load FairShip/latest`

echo "Starting script."
INPUTFILES=/afs/cern.ch/work/a/aiuliano/public/GenieEvents/condor_production_ccdis/seed_$LSB_JOBINDEX/
echo $INPUTFILES

OUTPUTDIR=/eos/user/a/aiuliano/public/sims_FairShip/sim_nutau/neutrinos2019_v2/
python $SHIPBUILD_mymaster/FairShip/macro/run_simScript.py --Genie -f $INPUTFILES/nu_e/genie-nu_e.root -n 10000 -o $OUTPUTDIR/$LSB_JOBINDEX/nu_e/

python $SHIPBUILD_mymaster/FairShip/macro/run_simScript.py --Genie -f $INPUTFILES/nu_mu/genie-nu_mu.root -n 10000 -o $OUTPUTDIR/$LSB_JOBINDEX/nu_mu/

python $SHIPBUILD_mymaster/FairShip/macro/run_simScript.py --Genie -f $INPUTFILES/nu_tau/genie-nu_tau.root -n 10000 -o $OUTPUTDIR/$LSB_JOBINDEX/nu_tau/
