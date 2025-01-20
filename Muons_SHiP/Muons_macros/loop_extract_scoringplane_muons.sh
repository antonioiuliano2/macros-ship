#!/bin/bash

for imuon in {0..66}
 do
  echo "process file $imuon"
  root -l -q /home/utente/Lavoro/macros-ship/Muons_SHiP/Muons_macros/ExtractVetoPointMuons.C\($imuon\)
 done
 