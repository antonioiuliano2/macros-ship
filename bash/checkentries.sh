#!/bin/bash
DATADIR=/eos/experiment/ship/data/charmxsec/

for FILENAME in $(ls $DATADIR/DATA_8000_charmTest/rootdata/RUN_8000_$1/);

do
 root -l -q checkentries.C\(\"$DATADIR/DATA_8000_charmTest/rootdata/RUN_8000_$1/$FILENAME\"\)
 echo " for file ",$FILENAME
 #python $FAIRSHIP/charmdet/CharmdetHitPositions.py -f $DATADIR/rootdata/RUN_8000_$1/$FILENAME -w -o $DATADIR/ntuples/RUN_8000_$1/ntuples_$FILENAME
done

echo "start SciFi files "

for FILENAME in $(ls $DATADIR/DATA_0900/rootdata/RUN_0900_0$1/);

do
 root -l -q checkentries.C\(\"$DATADIR/DATA_0900/rootdata/RUN_0900_0$1/$FILENAME\"\)
 echo " for file ",$FILENAME
 #python $FAIRSHIP/charmdet/CharmdetHitPositions.py -f $DATADIR/rootdata/RUN_8000_$1/$FILENAME -w -o $DATADIR/ntuples/RUN_8000_$1/ntuples_$FILENAME
done
