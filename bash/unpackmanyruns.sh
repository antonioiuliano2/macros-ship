#!/bin/bash 
#folder indicated by Maria Elena in date 26 June
for i in $(cat /afs/cern.ch/work/a/aiuliano/public/Charmdata/DATA_9000_test/runlist.txt);

 do
  echo "now processing run ", $i
  mkdir /eos/experiment/ship/data/charmxsec/DATA_0900/rootdata/RUN_0900_$i
  source unpackscifi.sh $i
 done
