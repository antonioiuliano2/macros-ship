#!/bin/bash
makescanset -set=1.0.0.0 -dzbase=195 -dz=-1300 -fromplate=$1 -to_plate=$2

emlink -set=1.0.0.0 -new -v=2
