#!/bin/bash

for imuon in {0..3398}
 do
  echo "process file $imuon"
  root -l -q /afs/cern.ch/work/a/aiuliano/public/macros-ship/Muons_SHiP/Muons_macros/Extract_EduardScoringPoints2024_Muons.C\($imuon\)
 done
 
