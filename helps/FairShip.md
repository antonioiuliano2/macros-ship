#General guide to FairShip and my simulations  (last change on 7 December 2018)
## Local installation

Follow Thomas guide to install (12th Collaboration Meeting).

After installation, to redo build (local modification, global updates, ...):

`> cd SHIPBuild`
`> FairShip/aliBuild.sh`

Set environment variables with the command:

`>alibuild/alienv enter --shellrc FairShip/latest`

# Simulation for Charm cross section and Muon flux measurements

Charmdet folder contains the geometry classes for the detectors used in charm cross section and Muon flux measurements.

## Simulations for charm cross section measurement

In general:

Option `--CharmdetSetup 1` activates charm cross section geometry, while
`--CharmdetSetup 0` activates muon flux geometry.

Only the most useful options have been explained here. For the complete list of
available options please refer to the related scripts.

Charm production simulations are done from `macro/run_simScript.py`. Example
syntax:    
`python $FAIRSHIP/macro/run_simScript.py --charm 1 -A charmonly
--CharmdetSetup 1 -f Cascadefile -n 1000 -o outputfolder`

Useful options:
* `--charm 1`: activates `charmdet` configuration instead of SHiP standard (both
  for charm cross section and muon flux measurements)
* `-A charmonly`: activates
  charm production simulation
* `-f`: input file with charm vertices (if you have Kerberos configured (e.g.
  by default on `lxplus`), this will be taken directly
  from `/eos/ship/data/Charm/Cascade-parp16-MSTP82-1-MSEL4-978Bpot.root` by
  default)
* `-n`: number of events
* `-o`: output of folder where geometry and output of simulation will be saved  

General POT simulations are done from
`muonShieldOptimization/run_MufluxfixedTarget.py`. Example of syntax:
 
 `python $FAIRSHIP/muonShieldOptimization/run_MufluxfixedTarget.py
 --CharmdetSetup 1 -G -e 0.001 -n 1000 -o outputfolder`
 
 It is a derivation of the fixed target simulation used in SHiP, applied to
 `charmdet` geometry.
 
 Useful options:
 
* `-e`: Energy cut for adding tracks to `Geant4` propagation (choosing a high cut
  allows to save memory for larger simulations)
* `-n`: Number of events
* `-o`: output of folder where geometry and output of simulation will be saved
* `-r`: number of run (can be used as folder naming if `-o` option is not used)
* `-f` : force overwriting of directory (DANGEROUS: if used in a wrong
    directory, it will delete it. DO NOT USE IT together with `-o` option)  


Different options for proton generation:

* `-V`: default one, proton interactions generated with `Pythia` and `EvtGen` is
  used for decays
* `-P`: both proton interactions and decays handled with `Pythia`
* `-G`: most basic simulation: one 400 GeV proton directly sent to `Geant4`

Details can be found here: <https://cds.cern.ch/record/2280572>.

All simulations use `Geant4` for propagation. IMPORTANT: Both `-V` and `-P`
generate all interactions in target. Due to small dimensions of target used in
`charmdet` measurement, many protons pass through without interacting. To
correctly simulate surviving protons and their tracks in detectors, use `-G`
option.

For any question or doubt about these simulations, contact Thomas Ruf
(<mailto:thomas.ruf@cern.ch>) or Antonio Iuliano
(<mailto:antonio.iuliano@cern.ch>)  

## Simulations for muon flux measurements

* For full simulation, proton (400 GeV) on SHiP muflux target, plus detector setup:
 `python  run_MufluxfixedTarget.py -e ecut -P --CharmdetSetup 0`
  more options available (boosting di-muon BR, di-muon cross sections,
  charmonium, charm, ...), see `--help`
* For fast simulation, with muons from external file:
`python  run_simScript.py --MuonBack -n 1000 --charm=1 --CharmdetSetup=0 -f
inputFile`, where `inputFile`:
`/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_10.0_withCharmandBeauty0_mu.root`
(+ 66 more files)
* For digitization:
`python runMufluxDigi.py -f ship.conical.FixedTarget-TGeant4.root -g
geofile_full.conical.FixedTarget-TGeant4.root`→
`ship.conical.FixedTarget-TGeant4_dig.root` Track reconstruction & analysis use
`drifttubeMonitoring.py`, will automatically detect if it is MC data.

# Simulation for neutrino interactions

Neutrino interactions are simulated with GENIE MonteCarlo software, then the produced particles are sent to Geant4 for propagation.  
Starting from the 2018 spectra produced by Thomas, they can be found in:  

* /eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_charm_nu_1.0.root  

* /eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_charm_nu_10.0.root  

Thomas suggested to use the 1 GeV file untile 10 GeV, then the 10 GeV files. So I created a neutino_merged.root file with these features.  

In this file you can find 1D histogram of neutrino momentum, along with 2D histograms of pt vs p. The 1D histograms are the input for GENIE, while the 2D ones will be used by FairShip generator to generate in the transverse plane.  

## Genie simulation

Most of this material here is credit to Annarita (<3)

First (and important!), we need to set Genie User Decay Settings, to tell Genie not to make tau leptons and charmed hadrons decay, because we want to let Geant4 decay them.

This is set in UserPhysicsOption.xml, setting the option DecayParticlewithCode = pdgcode, as false, with pdg codes for charmed hadrons:  

* ± 421    
* ± 411  
* ± 431
* 4122  
* 4112  
* 4212  
* 4222  

and for tau leptons:

* ± 16  

Export GXMLPATH variable with the path of the modified UserPhysicsOptions.xml will tell Genie your request.

If you do not have splines, you have to create them by launching 

`gmkspl -p #nu -t 1000822040[0.014], 1000822060[0.241], 1000822070[0.221],
1000822080[0.524] -n 500 -e 400 -o out_file_name.xml`  

I am currently using Annarita splines. Actually Genie manual does not recommend that, but to use official splines if no particular physical requests are asked. I did not see any difference last time I checked with the official splines. I will probably use official splines for future generations (not current ones)  

Launch Genie simulations with:

`gevgen -n Nevents -p neutrinocode -t targetcode -e  0.5,350  --run nrun -f neutrinofile,histoname
    --cross-sections splines --message-thresholds $GENIE/config/Messenger_laconic.xml --seed nseed`

Convert the output in gst format:

`gntpc -i gntp.0.ghep.root -f gst --message-thresholds $GENIE/config/Messenger_laconic.xml"`

All these operations are written in macro/makeGenieEvents.py 
Usage: 
\> python -i
\> makeEvent(100)
\> makeNtuples()

## FairShip simulation

Standard run_simScript.py simulation, with --Genie option and no options about the detector. Usage:

python $FairShip/macro/run_simScript.py --Genie -f genie_nu.root -n nevents -o outputfile.

Important: we have to provide Genie with the positions where we want to generate these neutrinos. Look for 

Geniegen.SetPositions()

usually we want to set it from the start to the end of the neutrino target.


## Using c++ macros

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
