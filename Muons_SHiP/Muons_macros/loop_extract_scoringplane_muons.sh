#!/bin/bash

for imuon in {0..100}
 do
  echo "process file $imuon"
  root -l -q /home/utente/Lavoro/macros-ship/Muons_SHiP/Muons_macros/Extract_EduardScoringPoints2024_Muons.C\($imuon\)
 done
 