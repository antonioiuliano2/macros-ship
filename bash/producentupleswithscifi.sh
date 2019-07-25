#!/bin/bash
NRUN=$1

SPILLCODELIST=/afs/cern.ch/work/a/aiuliano/public/Charmdata/spillcodes/spillcodes_run$NRUN.csv
DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/rootdata/RUN_8000_$NRUN
SCIFIDATADIR=/eos/experiment/ship/data/charmxsec/DATA_0900/rootdata/RUN_0900_$NRUN/
OUTPUTDIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/ntuples/test_250719/

for SPILLCODE in $(cat $SPILLCODELIST);

do
 echo "Now at spillcode "$SPILLCODE
 echo "NOW PROCESSING "$(ls $DATADIR/*$SPILLCODE*.root)" AND "$(ls $SCIFIDATADIR/*$SPILLCODE*.root)
 python $FAIRSHIP/charmdet/CharmdetHitPositions.py -f $DATADIR/*$SPILLCODE*.root -s $SCIFIDATADIR/*$SPILLCODE*.root -w -o $OUTPUTDIR/ntuples_$SPILLCODE.root
done
