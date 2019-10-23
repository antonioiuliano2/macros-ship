#!/bin/bash
#environment definition
DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest
ONLINEFOLDER=/afs/cern.ch/work/a/aiuliano/public/oliver_online

#loop on folders

for RUNFOLDER in $(ls $DATADIR/rawdata);
 do
 RUNNUMBER=${RUNFOLDER:(-4)} #string:-4 retrieves the last four characters of a string
 echo "processing run " $RUNNUMBER
 #loop on files for each run
 for FILENAME in $(ls $DATADIR/rawdata/$RUNFOLDER);
  do
   NOEXTENSIONFILENAME=${FILENAME/.raw/}

   INPUTFILE=$DATADIR/rawdata/RUN_8000_$RUNNUMBER/$FILENAME
   OUTPUTFILE=$DATADIR/rootdata/RUN_8000_$RUNNUMBER/$NOEXTENSIONFILENAME.root 

   echo "unpacking file "$INPUTFILE
   echo "into "  $OUTPUTFILE
   #note, if outputfile already exists, will be recreated by FairRoot->be careful to check if the paths are set correctly before launching, expecially for outputfiles!
   $ONLINEFOLDER/unpack.py -f $INPUTFILE -n $RUNNUMBER --charm -o $OUTPUTFILE
  done

 done
