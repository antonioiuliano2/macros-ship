#!/bin/bash
DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/
#starting loop on folders
for RUNFOLDER in $(ls $DATADIR/rawdata);
do
 RUNNUMBER=${RUNFOLDER:(-4)} #string:-4 retrieves the last four characters of a string
 echo "processing run " $RUNNUMBER

 for FILENAME in $(ls $DATADIR/rootdata/$RUNFOLDER/);

 do
  echo "NOW PROCESSING "$FILENAME
  python $FAIRSHIP/charmdet/CharmdetHitPositions.py -f $DATADIR/rootdata/$RUNFOLDER/$FILENAME -w -o $DATADIR/ntuples/$RUNFOLDER/ntuples_$FILENAME 
 done 
done
