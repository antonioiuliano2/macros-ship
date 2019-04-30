// *** Read ROOT Tree containing tracks reconstructed in the muon tagger (SHiP-charm TB) ***
//
// .L Read_RPC_Tracks.C
// Read_Tracks("./RPC_RecoTracks_run2793_s1f22ae29.root");//CHARM1-RUN6, spill1f22ae29, no ECC on the beamline
//
//

#include <TTree.h>
//#/home/antonio/Lavoro/Analisi/macros-ship/analisi_charmdata/rpc_charm
//-------------------------------------------------------------
void Read_Tracks(const char *filein){ //last update 14/11/18, A.P.

  if(!filein){
    cout<<"Input file not exists! Exiting Read_RPCTracks .."<<endl; return; }
  else{
    cout<<"*************************************************************"<<endl;
    cout<<"*                READ RPC Reconstructed Tracks              *"<<endl;
    cout<<"*************************************************************"<<endl;
    cout<<endl<<" RPC plane 0/1 = Horizontal/Vertical direction (for strips)"<<endl<<endl;
    cout<<"_____________________________________________________________"<<endl<<endl;
  }
  
  TFile *fin = new TFile(filein);
  TTree *RPC_Trks = (TTree*)fin->Get("RPC_RecoTracks");

  int id_run = -1, trigger = -1, nclusters = -1, id_track = -1;
  float trk_teta = -99., trk_phi = -99., trk_slopexz = -99., trk_slopeyz = -99.;
  char id_spill[10];
  std::vector<int> *cl_rpc = NULL;
  std::vector<float> *cl_x = NULL;
  std::vector<float> *cl_y = NULL;
  std::vector<float> *cl_z = NULL;
  std::vector<int> *cl_dir = NULL;//0 for H, 1 for V
  
  RPC_Trks->SetBranchAddress("id_run",&id_run);
  RPC_Trks->SetBranchAddress("trigger",&trigger);
  RPC_Trks->SetBranchAddress("nclusters",&nclusters);
  RPC_Trks->SetBranchAddress("id_track",&id_track);
  RPC_Trks->SetBranchAddress("id_spill",&id_spill);
  RPC_Trks->SetBranchAddress("trk_teta",&trk_teta);
  RPC_Trks->SetBranchAddress("trk_phi",&trk_phi);
  RPC_Trks->SetBranchAddress("trk_slopexz",&trk_slopexz);
  RPC_Trks->SetBranchAddress("trk_slopeyz",&trk_slopeyz);

  RPC_Trks->SetBranchAddress("cl_x",&cl_x);
  RPC_Trks->SetBranchAddress("cl_y",&cl_y);
  RPC_Trks->SetBranchAddress("cl_z",&cl_z);
  RPC_Trks->SetBranchAddress("cl_dir",&cl_dir);
  RPC_Trks->SetBranchAddress("cl_rpc",&cl_rpc);

  int ntracks = RPC_Trks->GetEntries();
  cout<<"RPC Tracks Tree with "<<ntracks<<" entries "<<endl;

  for(int i = 0; i<ntracks; i++){

    RPC_Trks->GetEntry(i);

    cout<<" Track "<< id_track <<" for run "<< id_run << ", spill "<< id_spill <<" and trigger "<<trigger<<endl;
    cout<<" Teta and phi angles: "<<trk_teta<<" , "<<trk_phi<<endl;
    cout<<" Slopes (xz, yz): ( "<<trk_slopexz<<", "<<trk_slopeyz<<" )"<<endl;
    cout<<" # clusters = "<<nclusters<<endl;

    for(int j = 0; j<nclusters; j++){

      cout<<"Cluster in rpc "<<cl_rpc->at(j)<<" in plane "<<cl_dir->at(j)<<" with coordinates (x,y,z) = ( "<< cl_x->at(j)<<", "<<cl_y->at(j)<<", "<<cl_z->at(j)<<" )"<<endl;

    }//for j

  }  //for i
  return;
}

void Load_Tracks_Hits(const char *filein, const char *fileout){
    if(!filein){
    cout<<"Input file not exists! Exiting Read_RPCTracks .."<<endl; return; }

  TFile *fin = TFile::Open(filein);
  TTree *RPC_Trks = (TTree*)fin->Get("RPC_RecoTracks");

  int id_run = -1, trigger = -1, nclusters = -1, id_track = -1;
  float trk_teta = -99., trk_phi = -99., trk_slopexz = -99., trk_slopeyz = -99.;
  char id_spill[10];
  std::vector<int> *cl_rpc = NULL;
  std::vector<float> *cl_x = NULL;
  std::vector<float> *cl_y = NULL;
  std::vector<float> *cl_z = NULL;
  std::vector<int> *cl_channel = NULL;//0 for H, 1 for V
  std::vector<int> *cl_dir = NULL;//0 for H, 1 for V
  
  //branches to read the tree
  RPC_Trks->SetBranchAddress("id_run",&id_run);
  RPC_Trks->SetBranchAddress("trigger",&trigger);
  RPC_Trks->SetBranchAddress("nclusters",&nclusters);
  RPC_Trks->SetBranchAddress("id_track",&id_track);
  RPC_Trks->SetBranchAddress("id_spill",&id_spill);
  RPC_Trks->SetBranchAddress("trk_teta",&trk_teta);
  RPC_Trks->SetBranchAddress("trk_phi",&trk_phi);
  RPC_Trks->SetBranchAddress("trk_slopexz",&trk_slopexz);
  RPC_Trks->SetBranchAddress("trk_slopeyz",&trk_slopeyz);

  RPC_Trks->SetBranchAddress("cl_channel",&cl_channel);
  RPC_Trks->SetBranchAddress("cl_x",&cl_x);
  RPC_Trks->SetBranchAddress("cl_y",&cl_y);
  RPC_Trks->SetBranchAddress("cl_z",&cl_z);
  RPC_Trks->SetBranchAddress("cl_dir",&cl_dir);
  RPC_Trks->SetBranchAddress("cl_rpc",&cl_rpc);

  int ntracks = RPC_Trks->GetEntries();
  cout<<"RPC Tracks Tree with "<<ntracks<<" entries "<<endl;

  //output file where data will be converted
  TFile *fout = new TFile(fileout, "RECREATE");
  TTree *FairShip_RPC_Trks = new TTree("cbmsim","Hits associated to reconstructed tracks");
  TClonesArray rpcarray("MuonTaggerHit",100);
  
  FairShip_RPC_Trks->Branch("MuonTaggerHit", &rpcarray);
  //FairShip_RPC_Trks->Branch("runID",&id_run);
  //FairShip_RPC_Trks->Branch("spillID",&id_spill);

  MuonTaggerHit *recocluster;
  for(int itrk = 0; itrk<ntracks; itrk++){ //main loop on tracks
    int iarray = 0;
    rpcarray.Clear();
    RPC_Trks->GetEntry(itrk); //getting the entry
    for (int icluster = 0; icluster < nclusters; icluster++){
      int detid = cl_dir->at(icluster)*1000 + cl_rpc->at(icluster) * 10000 + cl_channel->at(icluster); //coding used by FairShip to identify the RPC channel
      new (rpcarray[iarray]) MuonTaggerHit(detid,0.);
      iarray++;
      //rpcarray->Add(recocluster);      
     }
     FairShip_RPC_Trks->Fill();
    }
    FairShip_RPC_Trks->Write();
    fout->Close();
}
/*
void Load_Tracks(const char *filein, const char *fileout){ //last update 06/12/18, A.I.  
  if(!filein){
    cout<<"Input file not exists! Exiting Read_RPCTracks .."<<endl; return; }

  TFile *fin = TFile::Open(filein);
  TTree *RPC_Trks = (TTree*)fin->Get("RPC_RecoTracks");

  int id_run = -1, trigger = -1, nclusters = -1, id_track = -1;
  float trk_teta = -99., trk_phi = -99., trk_slopexz = -99., trk_slopeyz = -99.;
  char id_spill[10];
  std::vector<int> *cl_rpc = NULL;
  std::vector<float> *cl_x = NULL;
  std::vector<float> *cl_y = NULL;
  std::vector<float> *cl_z = NULL;
  std::vector<int> *cl_dir = NULL;//0 for H, 1 for V
  
  //branches to read the tree
  RPC_Trks->SetBranchAddress("id_run",&id_run);
  RPC_Trks->SetBranchAddress("trigger",&trigger);
  RPC_Trks->SetBranchAddress("nclusters",&nclusters);
  RPC_Trks->SetBranchAddress("id_track",&id_track);
  RPC_Trks->SetBranchAddress("id_spill",&id_spill);
  RPC_Trks->SetBranchAddress("trk_teta",&trk_teta);
  RPC_Trks->SetBranchAddress("trk_phi",&trk_phi);
  RPC_Trks->SetBranchAddress("trk_slopexz",&trk_slopexz);
  RPC_Trks->SetBranchAddress("trk_slopeyz",&trk_slopeyz);

  RPC_Trks->SetBranchAddress("cl_x",&cl_x);
  RPC_Trks->SetBranchAddress("cl_y",&cl_y);
  RPC_Trks->SetBranchAddress("cl_z",&cl_z);
  RPC_Trks->SetBranchAddress("cl_dir",&cl_dir);
  RPC_Trks->SetBranchAddress("cl_rpc",&cl_rpc);

  int ntracks = RPC_Trks->GetEntries();
  cout<<"RPC Tracks Tree with "<<ntracks<<" entries "<<endl;

  //output file where data will be converted
  TFile *fout = new TFile(fileout, "RECREATE");
  TTree *FairShip_RPC_Trks = new TTree("SHIPcharm_recodata","Reconstructed tracks in SHiP-charm detectors");
  TClonesArray *rpcarray = new TClonesArray("RPCTrack");

  FairShip_RPC_Trks->Branch("RpcTrack", &rpcarray);
  FairShip_RPC_Trks->Branch("runID",&id_run);
  //FairShip_RPC_Trks->Branch("spillID",&id_spill);

  RPCTrack *recotrack;

  int ntracks_perrun = 0;
  for(int itrk = 0; itrk<ntracks; itrk++){ //main loop on tracks

    RPC_Trks->GetEntry(itrk); //getting the entry
    if ((id_track == 1) && (itrk > 0)){ //if trackID is back to 1, get Event and reset all counters
      FairShip_RPC_Trks->Fill();
      rpcarray->Clear(); 
      ntracks_perrun = 0;      
    }

    new((*rpcarray)[ntracks_perrun]) RPCTrack(trk_teta, trk_phi);    

    for(int j = 0; j<nclusters; j++){ //loop on clusters associated to a track
          //rpcarray[ntracks_perrun].AddCluster(cl_x, cl_y, cl_z, cl_dir, cl_rpc);
    }

    ntracks_perrun++;
    delete recotrack;
  }
  return;
}*/
