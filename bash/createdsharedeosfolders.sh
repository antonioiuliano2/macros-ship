#create a set of folders with the same names as input main folder (I prefer to do it in a separate script for safety)

INPUTPATH=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/rootdata
DESTPATH=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/ntuples

for RUNFOLDER in $(ls $INPUTPATH);

do
 echo "creating folder " $DESTPATH/$RUNFOLDER
 mkdir $DESTPATH/$RUNFOLDER
done
