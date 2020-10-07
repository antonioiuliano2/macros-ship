#!/bin/bash

echo "SETUP"

quarter=$3

SHIPBUILD_mymaster=/afs/cern.ch/work/a/aiuliano/public/SHIPBuild
source /cvmfs/ship.cern.ch/SHiP-2020/latest/setUp.sh
eval `alienv load FairShip/latest`

source /afs/cern.ch/work/a/aiuliano/public/fedra/setup_new.sh

echo "Starting script for quarter ", $quarter

cd /afs/cern.ch/work/a/aiuliano/public/Charmdata/condor_reconstruction_emulsion/CH2R5/$quarter

root -l /afs/cern.ch/work/a/aiuliano/public/Charmdata/condor_reconstruction_emulsion/CH2R5/$quarter/charm_vertexing.C


