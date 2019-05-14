# How I am doing FairShip2FEDRA conversion

In order to reproduce our reconstruction procedure for emulsion data, I am following these steps:
* Launch a FairShip simulation (pot or charm, is the same format);
* Pass the 'BoxPoints' into 'EdbSegP couples' (x,y,tx,ty), with a dummy weight
* Add MC true information (MCEvt,MCTrackID)
* Launch FEDRA reconstruction (use -suff=cp.root in makescanset, since raw.root files are not present)