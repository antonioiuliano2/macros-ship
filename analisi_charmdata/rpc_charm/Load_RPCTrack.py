#reads locally reconstructed tracks by Bari group and builds RPCTrack objects to be stored in FairShip

import ROOT  as r

localRPCfile = r.TFile.Open("root:://eospublic.cern.ch//eos/experiment/ship/data/rpc_charm/RPC_RecoTracks_run2793_s1f22ade3.root")
sTree = localRPCfile.Get("RPC_RecoTracks")

nevents = sTree.GetEntries()

for i in range(nevents):
 print (i)
 sTree.GetEntry(i)
 nclusters = sTree.nclusters

 localRPCTrack = r.RPCTrack(sTree.trk_teta, sTree.trk_phi, sTree.trk_slopexz, sTree.trk_slopeyz)
 #adding elements to the container
 localRPCTrack.SetTrackID(sTree.id_track)
 localRPCTrack.SetRunNumber(sTree.id_run)

 spillname = r.TString(sTree.id_spill) #converting the char* in a TString
 localRPCTrack.SetSpillName(spillname) #currently not working due to char* to string conversion
 localRPCTrack.SetTrigger(sTree.trigger)

 print " Track ", localRPCTrack.GetTrackID(), " for run ", localRPCTrack.GetRunNumber(), ", spill ", localRPCTrack.GetSpillName().Data(), "and trigger ", localRPCTrack.GetTrigger() #<<" and trigger ",trigger<<endl;
 print " Teta and phi angles: ",localRPCTrack.GetTheta()," , ", localRPCTrack.GetPhi()
 print " Slopes (xz, yz): ( ",localRPCTrack.GetSlopeXZ(),", ", localRPCTrack.GetSlopeYZ()
 print " # clusters = ",nclusters;

 cl_x_list = sTree.cl_x
 cl_y_list = sTree.cl_y
 cl_z_list = sTree.cl_z
 cl_dir_list = sTree.cl_dir
 cl_rpc_list = sTree.cl_rpc

 for icluster in range (nclusters):
  #adding the cluster
  localRPCTrack.AddCluster(cl_x_list[icluster], cl_y_list[icluster], cl_z_list[icluster], cl_dir_list[icluster], cl_rpc_list[icluster])
  #getting information for print
  vector = localRPCTrack.GetClusterPos(icluster)
  print "Cluster in rpc ", localRPCTrack.GetClusterStation(icluster), "in plane ", localRPCTrack.GetClusterDir(icluster), "with coordinates (x,y,z) = ( ", vector[0], ", ", vector[1], ", ",vector[2]," )"
  #cout<<"Cluster in rpc "<<cl_rpc->at(j)<<" in plane "<<cl_dir->at(j)<<" with coordinates (x,y,z) = ( "<< cl_x->at(j)<<", "<<cl_y->at(j)<<", "<<cl_z->at(j)<<" )"<<endl;
