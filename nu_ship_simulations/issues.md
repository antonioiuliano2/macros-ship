# List of issues encountered during simulations

* Production Z not in target:
  Solution: remember to check GenieGen.SetPositions() before launching a new simulation;

* Overlaps:
  Always check for overlaps between volumes, since they usually screw up hit registration;

* DetectorID in SHiPRpcPoint:
  Ecountered issue due to TGeoAssemblyVolumes, most of IDs were 1. Solved with pull request;
