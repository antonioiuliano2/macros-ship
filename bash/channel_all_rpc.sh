#input: CHARM0,... CHARM1, CHARM2....CHARM6
DATADIR=/eos/experiment/ship/data/rpc_charm/

for FILENAME in $(ls $DATADIR$1/*/*.root); #when using *root, FILENAME acquires the full path
do
 FILENAME=${FILENAME//$DATADIR/}
 #FILENAME=${FILENAME//$1"/"/}
 echo "NOW PROCESSING ", $FILENAME
 root -l -q /afs/cern.ch/work/a/aiuliano/public/macros-ship/analisi_charmdata/rpc_charm/channel_rpc_clusters.C\(\"$FILENAME\"\)
done

