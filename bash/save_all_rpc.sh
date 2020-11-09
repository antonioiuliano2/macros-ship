#!/bin/bash
#input CHARM0, CHARM1, CHARM2.....
DATAPATH=/eos/experiment/ship/data/rpc_charm/with_hits

cat $DATAPATH/$1/*/file_pairs.txt | while read word myword

do 
 RPCTRACKSPATH=$(ls $DATAPATH/$1/*/$word)
 OTHERHITSPATH=$(ls $DATAPATH/$1/*/$myword)
 echo "NOW PROCESSING ", $RPCTRACKSPATH, $OTHERHITSPATH
 root -l -q /afs/cern.ch/work/a/aiuliano/public/macros-ship/analisi_charmdata/rpc_charm/Save_Tracks_Hits.C\(\"$RPCTRACKSPATH\",\"$OTHERHITSPATH\"\)
done
