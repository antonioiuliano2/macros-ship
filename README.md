# Personal repo of SHiP related macros

Reorganized to have clear (well, a bit more clear) folders for each project:

* shipcharm: not working on it anymore, kept mostly for archive;

* Muons_SHiP: input muons, especially to evaluate noise to emulsions and other detectors:

* muon_background: background evaluation from muon interactions. For now, muonDIS

* neutrinos_old, neutrions_2018, neutrinos_2025: input neutrino fluxes, and related codes. To be merged in a single folder

* condor_scripts: input for HTCondor submissions

* bash: bash related code, may be also inserted in the related project(s)

* nu_ship_simulations: all code to analyze neutrino event simulations after the FairShip/Geant4 simulation. Special case, efficiency scripts

* FEDRA: all FEDRA related stuff, may overlap a lot with macros-snd;

* general_utils: any small general script which can work on any project. I.e. drawrootvolume.

* testscripts: small templates used to test functions before adding them to a project;

* tutorials: scripts designed especially for some tutorials sessions;

* helps: some markdown notes, written for documentation before my gitbook repo;

## To do

Arrange more clear structure within the project subfolders.

I would avoid deleting old code, but I should avoid also being overwhelmed by it and unable to find useful code anymore.