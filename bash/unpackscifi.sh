#!/bin/bash
#DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/
DATADIR=/eos/experiment/ship/data/charmxsec
#need to define ONLINEFOLDER, for example 

ONLINEFOLDER=/afs/cern.ch/work/a/aiuliano/public/oliver_online/
#run unpacker on all raw data files (after Event Builder), source unpacktest.sh 2781
for FILENAME in $(ls $DATADIR/sciFiData/190625_updates/RUN_0900_$1/);

do
 echo "unpacking file "$FILENAME 
 $ONLINEFOLDER/unpack.py -f $DATADIR/sciFiData/190625_updates/RUN_0900_$1/$FILENAME -n $1 --charm -o $DATADIR/DATA_0900/rootdata/RUN_0900_$1/$FILENAME.root
#$DATADIR/rootdata/RUN_8000_$1/$NOEXTENSIONFILENAME.root
done
