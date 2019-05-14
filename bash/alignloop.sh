#!/bin/bash

cp firstalign.rootrc align.rootrc
makescanset -set=1.0.0.0 -dzbase=195 -dz=-1300 -from_plate=$1 -to_plate=$2 -v=2

echo "Starting pre-align on a small area"

emalign -set=1.0.0.0 -new -v=2

cp b000001.0.0.0.align.pdf b000001.0.0.0.prealign.pdf

cp secondalign.rootrc align.rootrc
makescanset -set=1.0.0.0 -dzbase=195 -dz=-1300 -from_plate=$1 -to_plate=$2 -v=2

echo "Starting true align on the full area"

emalign -set=1.0.0.0 -new -v=2
