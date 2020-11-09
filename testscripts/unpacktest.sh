#!/bin/bash
DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/
#need to define ONLINEFOLDER, for example 

ONLINEFOLDER=/afs/cern.ch/work/a/aiuliano/public/oliver_online/
#run unpacker on all raw data files (after Event Builder), source unpacktest.sh 2781
for FILENAME in $(ls /eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/rawdata/RUN_8000_$1/);

do
 NOEXTENSIONFILENAME=${FILENAME/.raw/}
 echo $NOEXTENSIONFILENAME.root
 echo "unpacking file "$FILENAME 
 $ONLINEFOLDER/unpack.py -f $DATADIR/rawdata/RUN_8000_$1/$FILENAME -n $1 --charm -o $DATADIR/rootdata/RUN_8000_$1/$NOEXTENSIONFILENAME.root
done
