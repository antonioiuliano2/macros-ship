#!/bin/bash
#launched with deepcopy option in run_simScript.py
ProcId=$2
Run=$3
NEVENTS=10000
STARTEVENT=$((Run*NEVENTS))
LSB_JOBINDEX=$((Run+1))
echo $LSB_JOBINDEX

sleep $ProcId

SHIPBUILD_mymaster=/afs/cern.ch/work/a/aiuliano/public/SHIPBuild
export ALIBUILD_WORK_DIR=$SHIPBUILD_mymaster/sw #for alienv

echo "SETUP"
source /cvmfs/ship.cern.ch/SHiP-2021/latest/setUp.sh
eval `alienv load FairShip/latest`

echo "Starting script, from neutrino number "
echo $STARTEVENT 
INPUTFILES=/afs/cern.ch/work/a/aiuliano/public/GenieEvents_SHIP/GenieEventsTungsten_12_2020
echo "From file "
echo $INPUTFILES

OUTPUTDIR=/afs/cern.ch/work/a/aiuliano/public/sim_nutau/nutau_CCDIS_tungsten_20_January_2021

/cvmfs/ship.cern.ch/SHiP-2021/2020/November/16/sw/slc7_x86-64/Python/v3.6.8-1/bin/python3 $SHIPBUILD_mymaster/FairShip/macro/run_simScript.py --Genie -f $INPUTFILES/CCDIS/nu_tau/genie-nu_tau.root -i $STARTEVENT -n $NEVENTS -o $OUTPUTDIR/$LSB_JOBINDEX
