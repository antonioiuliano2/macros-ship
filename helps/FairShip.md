#General guide to FairShip and my simulations  (last change on 7 December 2018)
## Local installation

Follow Thomas guide to install (12th Collaboration Meeting).

After installation, to redo build (local modification, global updates, ...):

`> cd SHIPBuild`
`> FairShip/aliBuild.sh`

Set environment variables with the command

`>alibuild/alienv enter --shellrc FairShip/latest`

The computing group name on lxplus is ship-cg. It can be set from the CERN resource page.

On Lxplus launch
`'cd SHIPBuild/FairShip'`
export SHIPBUILD=/cvmfs/ship.cern.ch/SHIPBuild
`./localBuild.sh`
source ../FairShipRun/config.sh

## Run FairShip simulation

Standard run_simScript.py simulation, with --Genie option and no options about the detector. Usage:

python $FairShip/macro/run_simScript.py --Genie -f genie_nu.root -n nevents -o outputfile.

Important: we have to provide Genie with the positions where we want to generate these neutrinos. Look for 

Geniegen.SetPositions()

usually we want to set it from the start to the end of the neutrino target.


## Using c++ macros
(Updated: This can be avoided by using TTreeReader and removing the include lines, thanks to cling (clang?) compiling and loading the shared libraries itself)

root -l  
`>>#include "FairMCPoint.h"`  
`>>#include "TMCProcess.h`  
`>>.x mymacro()`   


## Checking geometry output  Event display can be launched in the following way:
`python -i $FAIRSHIP/macro/eventDisplay.py -f simulationfile.root -g
geofile.root`

(actual names of `simulationfile.root` and `geofile.root` depend on the
launched simulation)

Positions and dimensions of volumes can be checked in the following way:
`python $FAIRSHIP/macro/getGeoInformation.py -g geofile.root`
 
Useful options:

* `-v`: name of the volume to expand (see list of volume daughters)
* `-l`: 'depth' level of the subnode expansion (how many daughters are showed)  
