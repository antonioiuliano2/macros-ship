#!/bin/bash
DATADIR=/eos/experiment/ship/data/charmxsec/DATA_8000_charmTest/

for FILENAME in $(ls $DATADIR/rootdata/RUN_8000_$1/);

do
 echo "NOW PROCESSING "$FILENAME
 python $FAIRSHIP/charmdet/CharmdetHitPositions.py -f $DATADIR/rootdata/RUN_8000_$1/$FILENAME -w -o $DATADIR/ntuples/RUN_8000_$1/ntuples_$FILENAME
done
