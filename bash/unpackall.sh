#!/bin/bash
#environment definition
DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/
ONLINEFOLDER=/afs/cern.ch/work/a/aiuliano/public/SHIPBuild/oliver_online/

#loop on folders

for RUNFOLDER in $(ls $DATADIR/rawdata);
 do
 RUNNUMBER=${RUNFOLDER:(-4)} #string:-4 retrieves the last four characters of a string
 echo "processing run " $RUNNUMBER
 #loop on files for each run
 for FILENAME in $(ls $DATADIR/rawdata/$RUNFOLDER);
  do
   NOEXTENSIONFILENAME=${FILENAME/.raw/}
   echo "unpacking file "$FILENAME "into "  $NOEXTENSIONFILENAME.root
   #note, if outputfile already exists, will be recreated by FairRoot->be careful to check if the paths are set correctly, expecially for outputfiles!
   #$ONLINEFOLDER/unpack.py -f $DATADIR/rawdata/RUN_8000_$1/$FILENAME -n $1 --   charm -o $DATADIR/rootdata/RUN_8000_$1/$NOEXTENSIONFILENAME.root 
   fi
  done

 done
