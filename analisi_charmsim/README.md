# FORMAT of FairShip simulation for pot and charm

## Charm format

Each event represents the decay of two charmed hadrons. The event contains the decay of the two hadrons, along with all the particles produced in the primary proton interactions:

Description of TrackID and MotherId (TrackID is not present as an attribute of MCTrack branch, but works as default index):

* The first track is always the mother proton (TrackID 0);

* The two charmed hadrons have then always MotherId equal to 0. Also, one charmed hadron has TrackID 1. This means that we always have charm daughters with MotherId == 1 (but they are not the only ones in the event);

* The other involved particles have motherID -1;

Event information is contained in the `MCEventHeader.fEventId`: warning, this starts from 1, whereas indeces `cbmsim->GetEntry(eventid)` start from 0.
