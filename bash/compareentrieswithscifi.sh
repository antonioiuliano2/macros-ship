#!/bin/bash
NRUN=$1

SPILLCODELIST=/afs/cern.ch/work/a/aiuliano/public/Charmdata/spillcodes/spillcodes_run$NRUN.csv
DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/rootdata/RUN_8000_$NRUN
SCIFIDATADIR=/eos/experiment/ship/data/charmxsec/DATA_0900/rootdata/RUN_0900_$NRUN/

echo "Entries Global vs SciFi"
for SPILLCODE in $(cat $SPILLCODELIST);

do
 echo "Now at spillcode "$SPILLCODE
 FILE1=$(ls $DATADIR/*$SPILLCODE*.root)
 FILE2=$(ls $SCIFIDATADIR/*$SPILLCODE*.root)
 #python compareentries.py $FILE1 $FILE2
 python /afs/cern.ch/work/a/aiuliano/public/macros-ship/bash/comparetimestamps.py $FILE1 $FILE2 $SPILLCODE
done

#unique name for file from each run
mv entrycomparison.root entrycomparison$NRUN.root
