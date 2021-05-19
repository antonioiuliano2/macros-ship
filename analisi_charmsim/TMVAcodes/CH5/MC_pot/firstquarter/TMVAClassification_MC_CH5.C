// V. GENTILE  LAST UPDATE ON (2020/02/11)

#include "TClonesArray.h"
#include <TChain.h>
#include <TFile.h>
#include <math.h>

using namespace std;


// Impact parameter of a track with respect to a vertex
float IPtoVertex(TVector3 vertexpos, TVector3 trackstartpos, float tracktx, float trackty){
 
 float dz = -vertexpos(2) + trackstartpos(2);
 float x0 = -tracktx * dz + trackstartpos(0) - vertexpos(0);
 float y0 = -trackty * dz + trackstartpos(1) - vertexpos(1);
 float m = (trackty/tracktx);
 float ip = TMath::Abs(-(y0-m*x0))/TMath::Sqrt(1+m*m);
 //cout << x0 << " " << y0 << " " << m << endl; 
 //float ip = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));

 return ip;
}

// IPmax over IPrms calculated on segments of a track with respet to a vertex
float SegIPtoVertex2(TVector3 vertexpos, int nseg, float* seg_x, float* seg_tx, float* seg_y, float* seg_ty, float* seg_z){
 
  float ipseg[nseg];
  float delta_ipseg[nseg-1];
  float ipseg_nomax[nseg-2];  
  for (int i = 0; i < nseg; i++){   
    float dz = -vertexpos(2) + seg_z[i];
    float ipx = -seg_tx[i] * dz + seg_x[i] - vertexpos(0);
    float ipy = -seg_ty[i] * dz + seg_y[i] - vertexpos(1);
    ipseg[i] = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));
  }

  for (int i = 0; i < nseg-1; i++){
    delta_ipseg[i] = TMath::Abs(ipseg[i+1] - ipseg[i]);
    //cout << i << " " << delta_ipseg[i] << endl;
  }

  float ipsegrms = TMath::RMS(nseg-1,delta_ipseg);
  return ipsegrms;
}

// Impact parameter of segments with respect to a vertex
void SegIPtoVertex3(TVector3 vertexpos, int nseg, float* seg_x, float* seg_tx, float* seg_y, float* seg_ty, float* seg_z, float trk_ip[3]){
 
  float ipseg[nseg];
  float delta_ipseg[nseg-1];
  float ipseg_nomax[nseg-2];
  
  for (int i = 0; i < nseg; i++){   
    float dz = -vertexpos(2) + seg_z[i];
    float ipx = -seg_tx[i] * dz + seg_x[i] - vertexpos(0);
    float ipy = -seg_ty[i] * dz + seg_y[i] - vertexpos(1);
    ipseg[i] = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));
  }

  for (int i = 0; i < nseg-1; i++){
    delta_ipseg[i] = TMath::Abs(ipseg[i+1] - ipseg[i]);
    //cout << i << " " << delta_ipseg[i] << endl;
  }

  //float trk_ip[3]={};
  trk_ip[0] = TMath::RMS(nseg-1,delta_ipseg);
  trk_ip[1] = TMath::Mean(nseg-1,delta_ipseg);
  trk_ip[2] = TMath::MaxElement(nseg-1,delta_ipseg);
  //return trk_ip;
}

//// Kinkmax over Kinkrms calculated on segments of a track with respet to a vertex
float FedraTrackKink2(int nseg, float* seg_tx, float* seg_ty){
  
  float kinkangles[nseg-1];
  float kinkangles_nomax[nseg-2];
  //loop on subsequent segment pairs of the tracks
  for (int i = 0; i < nseg-1; i++){
    kinkangles[i]=TMath::Sqrt(pow(seg_tx[i+1]-seg_tx[i],2)+pow(seg_ty[i+1]-seg_ty[i],2));
  }
  float deltathetarms = TMath::RMS(nseg-1, kinkangles);
  return deltathetarms;
}

// Kink angles between consecutive track segments
void FedraTrackKink3(int nseg, float* seg_tx, float* seg_ty, float trk_kink[3]){

  int dim = (nseg*(nseg-1))/2.;
  int idelta=0;
  float kinkangles[nseg-1];
  //loop on random segment pairs of the tracks
  for (int i = 0; i < nseg-1; i++){
    for (int j = (i+1); j < nseg-1; j++){
      kinkangles[idelta]=TMath::Sqrt(pow(seg_tx[j]-seg_tx[i],2)+pow(seg_ty[j]-seg_ty[i],2));
      idelta++;
    }
  }
  //float trk_kink[3]={};
  trk_kink[0] = TMath::RMS(nseg-1, kinkangles);
  trk_kink[1] = TMath::Mean(nseg-1, kinkangles);
  trk_kink[2] = TMath::MaxElement(nseg-1, kinkangles);
  //return trk_kink;
  }


// Mutual kink angles among track semgents
void FedraTrackKinkRandom(int nseg, float* seg_tx, float* seg_ty, float trk_kink[3]){

  int dim = (nseg*(nseg-1))/2.;
  int idelta=0;
  float kinkangles[dim];
  //loop on random segment pairs of the tracks
  for (int i = 0; i < nseg; i++){
    for (int j =(i+1); j < nseg; j++){
      kinkangles[idelta]=TMath::Sqrt(pow(seg_tx[j]-seg_tx[i],2)+pow(seg_ty[j]-seg_ty[i],2));
      //cout <<"kink 3 "<< dim << " " << nseg << " " << idelta << " " << kinkangles[idelta] << endl;
      idelta++;
    }
  }
  //float trk_kink[3]={};
  trk_kink[0] = TMath::RMS(dim, kinkangles);
  trk_kink[1] = TMath::Mean(dim, kinkangles);
  trk_kink[2] = TMath::MaxElement(dim, kinkangles);
  //return trk_kink;
}

// BDT INPUT FOR INTERACTION VERTICES
void prepareTMVAtree(){

   int choice=0;
  
   cout << "Scegli il dataset (1: CHARM simulation; 2: BACKGROUND simulation\t";
   cin >> choice;
   cout << endl;

   TFile *f2(0);
   TFile *f(0);
   TString f_name="";
   TString f2_name="";
   switch(choice){
   case 1:
     f_name = "../vtx_MC_analysis.root";
     f = TFile::Open( f_name );
     if(!f){
       std::cout << "ERROR: could not open data file" << std::endl;
       exit(1);
     }
     break;
   case 2:
     f_name = "vtx_MC_analysis_POT.root";
     f = TFile::Open( f_name );
     if(!f) {
       std::cout << "ERROR: could not open data file" << std::endl;
       exit(1);
     }
     break;
   }
   
   f->cd();   
   TTreeReader vtxreader("vtx",f);

   TTreeReaderArray<EdbSegP> segments(vtxreader,"s");
   TTreeReaderArray<EdbSegP> tracks(vtxreader,"t.");
   TTreeReaderArray<float> IP(vtxreader,"impactparameter");
   TTreeReaderArray<int> nseg(vtxreader,"nseg");
   TTreeReaderValue<float> maxap(vtxreader,"maxaperture");
   TTreeReaderValue<float> prob(vtxreader,"probability");
   TTreeReaderValue<int> n(vtxreader,"n");
   TTreeReaderValue<int> Flag(vtxreader,"flag");
   TTreeReaderValue<int> vid(vtxreader,"vID");
   TTreeReaderArray<int> Incoming(vtxreader,"incoming");
   TTreeReaderValue<float> v_x(vtxreader,"vx");
   TTreeReaderValue<float> v_y(vtxreader,"vy");
   TTreeReaderValue<float> v_z(vtxreader,"vz");
   TTreeReaderValue<char> vgood(vtxreader,"vtx_good");
   TTreeReaderValue<int> mppdg(vtxreader,"mp_pdgID");
   TTreeReaderValue<int> mpmother(vtxreader,"mp_motherID");
   TTreeReaderValue<int> mpevent(vtxreader,"mp_eventID");
   TTreeReaderValue<int> n_mothers(vtxreader,"n_mothers");
   TTreeReaderValue<int> n_events(vtxreader,"n_events");
   TTreeReaderValue<int> n_electrons(vtxreader,"n_electrons");
   TTreeReaderValue<int> max_freq(vtxreader,"max_freq");
   TTreeReaderValue<float> mean_freq(vtxreader,"mean_freq");
   TTreeReaderValue<bool> charmdaug(vtxreader,"charm_daug");    //no POT
  
  
  
  
  Char_t vtx_good;
  Int_t vID, mp_pdgID,ntracks, mp_motherID,flag, mp_eventID, nelectrons, maxfreq, nmothers, nevents;
  Float_t meanTX, meanTY,meannseg, maxaperture, probability,vx,vy,vz,meanfreq;
  Float_t maxIP, meanIP, vtx_node, vtx_holes, vtx_gap, vtx_fill;
  Bool_t charm_daug;
  
  const Int_t maxtracks = 1000; //I do not expect a larger number of tracks than this
  const Int_t maxseg=30;
  Int_t incoming[maxtracks], trk_num_holes[maxtracks], trk_max_gap[maxtracks];
  Float_t t_eTX[maxtracks], t_eTY[maxtracks], t_eX[maxtracks], t_eY[maxtracks], t_eZ[maxtracks], trackfill[maxtracks],impactparameter[maxtracks];
  Int_t s_ePID[maxseg], t__, t_eNseg[maxtracks];
  vector<float> s_eX, s_eY, s_eZ;
  vector<int> s_eTID;


  TString outfileName("");
  switch(choice){
  case 1:
    outfileName = "tmva_input_vertices.root";
    break;
  case 2:
    outfileName = "tmva_input_vertices_POT.root";
    break;
  }
  
  TFile * outputfile = new TFile(outfileName,"RECREATE");
  
  TTree *outputvtx = new TTree("vertices","Input file for TMVA analysis");
  
  outputvtx->Branch("vtx_good",&vtx_good,"vtx_good/B");
  outputvtx->Branch("mp_eventID",&mp_eventID,"mp_eventID/I");
  outputvtx->Branch("mp_pdgID",&mp_pdgID,"mp_pdgID/I");
  outputvtx->Branch("mp_motherID",&mp_motherID,"mp_motherID/I");
  outputvtx->Branch("charm_daug",&charm_daug,"charm_daug/O");
  outputvtx->Branch("nevents",&nevents,"nevents/I");
  outputvtx->Branch("nmothers",&nmothers,"nmothers/I");
  outputvtx->Branch("nelectrons",&nelectrons,"nelectrons/I");
  outputvtx->Branch("maxfreq",&maxfreq,"maxfreq/I");
  outputvtx->Branch("meanfreq",&meanfreq,"meanfreq/F");
   outputvtx->Branch("vID",&vID,"vID/I");
   outputvtx->Branch("vx",&vx,"vx/F");
   outputvtx->Branch("vy",&vy,"vy/F");
   outputvtx->Branch("vz",&vz,"vz/F");
   outputvtx->Branch("ntracks",&ntracks,"ntracks/I");
   outputvtx->Branch("probability",&probability,"probability/F");
   outputvtx->Branch("maxaperture",&maxaperture,"maxaperture/F");
   outputvtx->Branch("meannseg",&meannseg,"meannseg/F");
   outputvtx->Branch("meanTX",&meanTX,"meanTX/F");
   outputvtx->Branch("meanTY",&meanTY,"meanTY/F");
   outputvtx->Branch("maxIP",&maxIP,"maxIP/F");
   outputvtx->Branch("meanIP",&meanIP,"meanIP/F");
   outputvtx->Branch("vtx_fill",&vtx_fill,"vtx_fill/F");
   outputvtx->Branch("t_eNseg",t_eNseg,"t_eNseg[ntracks]/I");
   outputvtx->Branch("t_eX",t_eX,"t_eX[ntracks]/F");
   outputvtx->Branch("t_eY",t_eY,"t_eY[ntracks]/F");
   outputvtx->Branch("t_eZ",t_eZ,"t_eZ[ntracks]/F");
   outputvtx->Branch("t_eTX",t_eTX,"t_eTX[ntracks]/F");
   outputvtx->Branch("t_eTY",t_eTY,"t_eTY[ntracks]/F");
   outputvtx->Branch("s_eTID",&s_eTID);
   outputvtx->Branch("s_eX",&s_eX);
   outputvtx->Branch("s_eY",&s_eY);
   outputvtx->Branch("s_eZ",&s_eZ);
  const Int_t nvertices = vtxreader.GetEntries();
  //cout << nvertices << endl;
  
  for (int ivtx=0;ivtx<nvertices;ivtx++){
    
    vtxreader.Next();
    
    meannseg = 0.;
    meanTX = 0.;
    meanTY = 0.;
    maxIP = 0.;
    meanIP = 0.;
    
    vtx_node=0.;
    vtx_holes=0.;
    vtx_gap=0.;
    vtx_fill=0.;
    
    vID=*vid;
    vx=*v_x;
    vy=*v_y;
    vz=*v_z;
    ntracks = *n;
    flag =*Flag;
    vtx_good = *vgood;
    mp_eventID = *mpevent;
    mp_pdgID = *mppdg;
    mp_motherID = *mpmother;
    charm_daug = *charmdaug;   //no POT
    nevents = *n_events;
    nmothers = *n_mothers;
    nelectrons = *n_electrons;
    maxfreq = *max_freq;
    meanfreq = *mean_freq;
    maxaperture = *maxap;
    probability = *prob;

    
    if(flag!=2 && flag!=5 && vz<-3000 && ntracks>=4){ // TOPOLOGIC CUT
      
      int index_tracks=0;
      int sum_seg=0;
      
      for (int itrk = 0; itrk < ntracks; itrk++){

	t_eTX[itrk] = tracks[itrk].TX();
	t_eTY[itrk] = tracks[itrk].TY();
	t_eX[itrk] = tracks[itrk].X();
	t_eY[itrk] = tracks[itrk].Y();
	t_eZ[itrk] = tracks[itrk].Z();
	t_eNseg[itrk] = nseg[itrk];

      
	for(int u=sum_seg;u<(sum_seg+nseg[itrk]);u++){
	  s_eTID.push_back(itrk);
	  s_eX.push_back(segments[u].X());
	  s_eY.push_back(segments[u].Y());
	  s_eZ.push_back(segments[u].Z());
	}

	float trk_fill=0;
	meannseg = meannseg + nseg[itrk];
	meanTX = meanTX + tracks[itrk].TX();
	meanTY = meanTY + tracks[itrk].TY();
	meanIP = meanIP + IP[itrk];
	if (maxIP < IP[itrk]) maxIP = IP[itrk];
	
	// FILL FACTOR
	float plate_step=1300;
	int remaining_plates = floor(TMath::Abs(vz)/plate_step) + 1;
	if(remaining_plates!=0 && Incoming[itrk]==1){
	  trk_fill = (float) (nseg[itrk]) /(float) (remaining_plates);
	  index_tracks++;
	}
	else trk_fill=0;
	
	vtx_fill += trk_fill; 
	sum_seg += nseg[itrk];
      }
      
      meannseg = (Float_t)(meannseg/(ntracks));
      meanTX = (Float_t)(meanTX/(ntracks));
      meanTY = (Float_t)(meanTY/(ntracks));
      meanIP = (Float_t)(meanIP/(ntracks));      
      if(index_tracks!=0)vtx_fill = (Float_t)(vtx_fill/(index_tracks));
      else vtx_fill=0;
      
      outputvtx->Fill();
    } // flag statement

    s_eTID.clear();
    s_eX.clear();
    s_eY.clear();
    s_eZ.clear();
  }
  //f->Close();
  outputfile->Write();
  outputfile->Close();
}



// BDT INPUT FOR DECAY VERTICES 
void prepareTMVAtree2nd(){

   int choice=0;
  
   cout << "Scegli il dataset (1: CHARM simulation; 2: BACKGROUND simulation)\t";
   cin >> choice;
   cout << endl;

   TFile *f2(0);
   TFile *f(0);
   TString f_name="";
   TString f2_name="";
   switch(choice){
   case 1:
     f2_name = "vtx_BDT_data_evaluated.root";     // tutto
     f_name = "vtx_MC_analysis.root";
     f = TFile::Open( f_name );
     f2 = TFile::Open( f2_name );
     if(!f || !f2) {
       std::cout << "ERROR: could not open data file" << std::endl;
       exit(1);
     }
     break;
   case 2:
     f2_name = "vtx_BDT_data_evaluated_POT.root";     // tutto
     f_name = "vtx_MC_analysis_POT.root";
     if (!gSystem->AccessPathName(f_name )) {
       f = TFile::Open( f_name );
       f2 = TFile::Open( f2_name );
     }
     if(!f || !f2) {
       std::cout << "ERROR: could not open data file" << std::endl;
       exit(1);
     }
     break;
   }
     
   f->cd();   
   TTreeReader vtxreader("vtx",f);
   f2->cd();   
   TTreeReader mckreader("bdt",f2);
   
   TTreeReaderArray<EdbSegP> segments(vtxreader,"s");
   TTreeReaderArray<EdbSegP> tracks(vtxreader,"t.");
   TTreeReaderArray<float> IP(vtxreader,"impactparameter");
   TTreeReaderArray<int> nseg(vtxreader,"nseg");
   TTreeReaderValue<float> maxap(vtxreader,"maxaperture");
   TTreeReaderValue<float> prob(vtxreader,"probability");
   TTreeReaderValue<int> n(vtxreader,"n");
   TTreeReaderValue<int> vid(vtxreader,"vID");
   TTreeReaderArray<int> Incoming(vtxreader,"incoming");
   TTreeReaderValue<int> Flag(vtxreader,"flag");
   TTreeReaderValue<float> v_x(vtxreader,"vx");
   TTreeReaderValue<float> v_y(vtxreader,"vy");
   TTreeReaderValue<float> v_z(vtxreader,"vz");
   TTreeReaderValue<char> vgood(vtxreader,"vtx_good");
   TTreeReaderValue<int> mpevent(vtxreader,"mp_eventID");
   TTreeReaderValue<int> mppdg(vtxreader,"mp_pdgID");
   TTreeReaderValue<int> mpmother(vtxreader,"mp_motherID");
   TTreeReaderValue<bool> charmdaug(vtxreader,"charm_daug");    //no POT
   TTreeReaderValue<float> bdt_value(mckreader,"bdt_value");
   TTreeReaderValue<int> bdt_vID(mckreader,"bdt_vID");
  
   Char_t vtx_good;
   Int_t vID, mp_pdgID,ntracks,mp_motherID,flag, mp_eventID, bdt_vtxID;
   Float_t bdt_val;
   Float_t meanTX, meanTY,meannseg, maxaperture, probability,vx,vy,vz;
   Float_t maxIP, meanIP, vtx_node, vtx_holes, vtx_gap, vtx_fill;
   float mean_rms_ip=0;
   float mean_rms_kink=0;
   
   const Int_t maxtracks = 1000; //I do not expect a larger number of tracks than this
   const Int_t maxseg=30;
   Int_t incoming[maxtracks], trk_num_holes[maxtracks], trk_max_gap[maxtracks];
   Float_t t_eTX[maxtracks], t_eTY[maxtracks],  trackfill[maxtracks],impactparameter[maxtracks];
   Int_t s_ePID[maxseg], t__;
   Bool_t charm_daug;
   
   float cut_bdt_val;
   cout << "Inserire cut BDT interaction ";
   cin >> cut_bdt_val;
   cout << endl;
   
   TString outfileName("");
   switch(choice){
   case 1:
     outfileName = "tmva_input_vertices2nd.root";
     break;
   case 2:
     outfileName = "tmva_input_vertices2nd_POT.root";
     break;
   }
   
   TFile * outputfile = new TFile(outfileName,"RECREATE");
  TTree *outputvtx = new TTree("vertices","Input file for TMVA analysis");
  
  outputvtx->Branch("vtx_good",&vtx_good,"vtx_good/B");
  outputvtx->Branch("mp_eventID",&mp_eventID,"mp_eventID/I");
  outputvtx->Branch("mp_pdgID",&mp_pdgID,"mp_pdgID/I");
  outputvtx->Branch("mp_motherID",&mp_motherID,"mp_motherID/I");
  outputvtx->Branch("charm_daug",&charm_daug,"charm_daug/O");
  outputvtx->Branch("vID",&vID,"vID/I");
  outputvtx->Branch("vx",&vx,"vx/F");
  outputvtx->Branch("vy",&vy,"vy/F");
  outputvtx->Branch("vz",&vz,"vz/F");
  outputvtx->Branch("ntracks",&ntracks,"ntracks/I");
  outputvtx->Branch("probability",&probability,"probability/F");
  outputvtx->Branch("maxaperture",&maxaperture,"maxaperture/F");
  outputvtx->Branch("meannseg",&meannseg,"meannseg/F");
  outputvtx->Branch("meanTX",&meanTX,"meanTX/F");
  outputvtx->Branch("meanTY",&meanTY,"meanTY/F");
  outputvtx->Branch("maxIP",&maxIP,"maxIP/F");
  outputvtx->Branch("meanIP",&meanIP,"meanIP/F");
  outputvtx->Branch("mean_rms_kink",&mean_rms_kink,"mean_rms_kink/F");
  outputvtx->Branch("mean_rms_ip",&mean_rms_ip,"mean_rms_ip/F");
  outputvtx->Branch("vtx_fill",&vtx_fill,"vtx_fill/F");
  outputvtx->Branch("bdt_value",&bdt_val,"bdt_value/F");
  
  const Int_t nvertices = vtxreader.GetEntries();
  const Int_t nvtx = vtxreader.GetTree()->GetMaximum("vID");

  float BDT1_val[nvtx];
  for (int ivtx=0;ivtx<nvtx;ivtx++){
    BDT1_val[ivtx]=-1;   
  }

  const Int_t nbdt1 = mckreader.GetEntries();
  
  for (int ivtx=0;ivtx<nbdt1;ivtx++){
    mckreader.Next();
    bdt_val = *bdt_value;
    bdt_vtxID = *bdt_vID;
    //cout << bdt_vtxID << " " << bdt_val << endl;
    BDT1_val[bdt_vtxID]=bdt_val;
  }
  
  for (int ivtx=0;ivtx<nvertices;ivtx++){

    vtxreader.Next();
    //mckreader.Next();
    
    meanTX = 0.;
    meanTY = 0.;
    maxIP = 0.;
    meanIP = 0.;
    
    vtx_node=0.;
    vtx_holes=0.;
    vtx_gap=0.;
    vtx_fill=0.;

    vID=*vid;
    vx=*v_x;
    vy=*v_y;
    vz=*v_z;
    ntracks = *n;
    flag =*Flag;
    
    vtx_good = *vgood;
    mp_eventID = *mpevent;
    mp_pdgID = *mppdg;
    mp_motherID = *mpmother;
    charm_daug = *charmdaug;  // no POT
    
    maxaperture = *maxap;
    probability = *prob;
    //bdt_val = *bdt_value;
    bdt_val = BDT1_val[vID];
    
    int sum_seg=0;
    int index_tracks=0;
    int index_trk=0;
    int index_seg=0;

    mean_rms_ip=0;
    mean_rms_kink=0;

    TVector3 vtx_pos;
    vtx_pos.SetXYZ(vx,vy,vz);
    if((bdt_val<cut_bdt_val || (bdt_val>cut_bdt_val && ntracks<=8)) && flag!=2 && flag!=5 && vz<-3000){  // decay vertices dataset

      for (int itrk = 0; itrk < ntracks; itrk++){
	
	float trk_kink_random[3]={};
	float trk_fill=0;
	meannseg = meannseg + nseg[itrk];
	//cout << meannseg << " " << nseg[itrk] << endl;
	meanTX = meanTX + tracks[itrk].TX();
	meanTY = meanTY + tracks[itrk].TY();
	meanIP = meanIP + IP[itrk];
	if (maxIP < IP[itrk]) maxIP = IP[itrk];
	
	
	// segments
	if(nseg[itrk]>=3){
	  
	  float seg_x[nseg[itrk]];
	  float seg_y[nseg[itrk]];
	  float seg_z[nseg[itrk]];
	  float seg_tx[nseg[itrk]];
	  float seg_ty[nseg[itrk]];
	  float seg_mc_trk[nseg[itrk]];
	  float seg_mc_evt[nseg[itrk]];
	  
	  for(int u=sum_seg;u<(sum_seg+nseg[itrk]);u++){
	    seg_x[u-sum_seg]=segments[u].X();
	    seg_y[u-sum_seg]=segments[u].Y();
	    seg_z[u-sum_seg]=segments[u].Z();
	    seg_tx[u-sum_seg]=segments[u].TX();
	    seg_ty[u-sum_seg]=segments[u].TY();
	  }
	  
	  FedraTrackKinkRandom(nseg[itrk],seg_tx,seg_ty,trk_kink_random);
	  mean_rms_kink += trk_kink_random[0]*nseg[itrk]; 	 
	  mean_rms_ip += SegIPtoVertex2(vtx_pos, nseg[itrk], seg_x, seg_tx, seg_y, seg_ty, seg_z)*nseg[itrk];
	  index_trk++;
	  index_seg += nseg[itrk];
	}
	else{
	  trk_kink_random[0]=0;
	  trk_kink_random[1]=0;
	  trk_kink_random[2]=0;
	
	}

	// FILL FACTOR
	float plate_step=1300;
	int remaining_plates = floor(TMath::Abs(vz)/plate_step) + 1;
	if(remaining_plates!=0 && Incoming[itrk]==1){
	  trk_fill = (float) (nseg[itrk]) /(float) (remaining_plates);
	  index_tracks++;
	}
	else trk_fill=0;
	
	vtx_fill += trk_fill;  
	sum_seg += nseg[itrk];
	
      }
      
      meannseg = (Float_t)(meannseg/(ntracks));
      meanTX = (Float_t)(meanTX/(ntracks));
      meanTY = (Float_t)(meanTY/(ntracks));
      meanIP = (Float_t)(meanIP/(ntracks));
      vtx_fill = (Float_t)(vtx_fill/(index_tracks));

      // MEAN RMS IP AND KINK WEIGHTED FOR THE NUMBER OF SEGMENTS
      if(index_trk!=0){
      mean_rms_kink = (Float_t)(mean_rms_kink/(index_seg));
      mean_rms_ip = (Float_t)(mean_rms_ip/(index_seg));
      }
      else {
	mean_rms_kink = -1;
	mean_rms_ip = -1;
      }
      
      outputvtx->Fill();
    }
  }
  
  outputfile->Write();
  outputfile->Close();
  
}



void prepareTMVAtreeDS(){

  //using T = int;
  //using V = std::vector<T>;
  //using VC = std::vector<T,ROOT::Detail::VecOps::RAdoptAllocator<T>>;

 
  TFile *f = TFile::Open("annotated_ds_data_result.root");
  TTreeReader dsreader("ds",f);
  //TTreeReader dsshower("dsshower",f);
 
  TTreeReaderValue<int> Vtx_id1(dsreader,"vtx_fe_id");
  TTreeReaderValue<int> Vtx_mc_ev(dsreader,"vtx_mc_ev");
  TTreeReaderArray<int> Vtx2_mc_ev(dsreader,"dsvtx_vtx2_mc_ev");
  TTreeReaderArray<int> Vtx_id2(dsreader,"dsvtx_vtx2_vid");
  TTreeReaderArray<int> Trk_mc_id(dsreader,"dsvtx_vtx2_mc_tid");
  TTreeReaderArray<int> Trk_id(dsreader,"dsvtx_vtx2_tid");
  TTreeReaderArray<int> Trk_pid(dsreader,"dsvtx_vtx2_mc_pid");
  TTreeReaderValue<int> Event(dsreader,"vtx_topology");
  TTreeReaderArray<float> ka(dsreader,"dsvtx_vtx2_vka");
  TTreeReaderArray<float> ip(dsreader,"dsvtx_vtx2_vip");
  TTreeReaderArray<float> Trk_pms(dsreader,"dsvtx_vtx2_trk_pms");
  TTreeReaderArray<float> Tka_rms(dsreader,"dsvtx_vtx2_tka_rms");
  TTreeReaderArray<float> Tka_mean(dsreader,"dsvtx_vtx2_tka_mean");
  TTreeReaderArray<float> Tka_max(dsreader,"dsvtx_vtx2_tka_max");
  TTreeReaderArray<float> TX(dsreader,"dsvtx_vtx2_tx");
  TTreeReaderArray<float> TY(dsreader,"dsvtx_vtx2_ty");
  TTreeReaderArray<float> Trk_x(dsreader,"dsvtx_vtx2_xt");
  TTreeReaderArray<float> Trk_y(dsreader,"dsvtx_vtx2_yt");
  TTreeReaderArray<float> Trk_z(dsreader,"dsvtx_vtx2_zt");
  TTreeReaderArray<int> Trk_incoming(dsreader,"dsvtx_vtx2_incoming");
  TTreeReaderArray<int> Nseg(dsreader,"dsvtx_vtx2_tnseg");  
  TTreeReaderArray<int> Zpositive(dsreader,"dsvtx_vtx2_positivedz"); 
  TTreeReaderArray<int> Samevent(dsreader,"dsvtx_vtx2_samevent");
  TTreeReaderArray<int> Charmdaug(dsreader,"dsvtx_vtx2_charmdaughter");
  TTreeReaderArray<float> Decaylength(dsreader,"dsvtx_vtx2_dl");
  TTreeReaderArray<float> Dphi(dsreader,"dsvtx_vtx2_phidifference");
  TTreeReaderArray<float> Pointing(dsreader,"dsvtx_vtx2_endingdeltaphi");
  TTreeReaderArray<float> Meanphi(dsreader,"dsvtx_vtx2_meanphi");
  TTreeReaderArray<float> Meanlife(dsreader,"dsvtx_vtx2_tau");
  TTreeReaderValue<vector<int>> Vtx2_ntrk(dsreader,"dsvtx_vtx2_ntrk");
  TTreeReaderValue<float> Vtx_x(dsreader,"vtx_x");
  TTreeReaderValue<float> Vtx_y(dsreader,"vtx_y");
  TTreeReaderValue<float> Vtx_z(dsreader,"vtx_z");
  TTreeReaderArray<float> Vtx2_x(dsreader,"dsvtx_vtx2_vx");
  TTreeReaderArray<float> Vtx2_y(dsreader,"dsvtx_vtx2_vy");
  TTreeReaderArray<float> Vtx2_z(dsreader,"dsvtx_vtx2_vz");
  TTreeReaderArray<float> Vtx2_prob(dsreader,"dsvtx_vtx2_prob");
  TTreeReaderArray<float> Vtx2_bdt2(dsreader,"dsvtx_vtx2_bdt2");
  TTreeReaderValue<float> Vtx_bdt1(dsreader,"vtx_bdt1");
  TTreeReaderValue<int> Vtx_ntrk(dsreader,"vtx_ntrk");
  TTreeReaderValue<float> Vtx_maxaperture(dsreader,"vtx_maxap");
  //TTreeReaderArray<float> Out15(dsshower,"dsvtx_vtx2_trk_output15");
  //TTreeReaderArray<float> Out30(dsshower,"dsvtx_vtx2_trk_output30");
  //TTreeReaderArray<int> Nshower(dsshower,"dsvtx_vtx2_trk_nshower");
  //TTreeReaderArray<int> Sizeb(dsshower,"dsvtx_vtx2_trk_sizeb");
  
  int nseg, ntracks, npl, nholes, maxgap, event, charm1, charm2, vtx_ntrk, starting2, incoming;
  int zpositive, samevent,charmdaughter, vtx_id, vtx2_id, trk_id, vtx_mc_ev, vtx2_mc_ev, trk_pid, trk_mc_id, goodvtx, vtx2_ntrk, trk_nshower, trk_sizeb, trk_shower;
  float kink, impact, tka_rms, tka_mean, tka_max, tip_rms, tip_mean, vtx_maxap, trk_x, trk_y, trk_z, trk_tx, trk_ty, tip_vtx2, vtx_bdt1, trk_pms, trk_fill, vtx2_maxap;
  float tip_max, decaylength, dist_tr, vtx_x, vtx_y, vtx_z, vtx2_x, vtx2_y, vtx2_z, vtx2_prob, meanlife, vtx2_bdt2, dphi, dphi2, meannseg, meanfill, pointing, trk_shower_val15, trk_shower_val30, trk_shower_mean, trk_shower_rms, trk_shower_min, trk_shower_max;
  
  TFile * outputfile = new TFile("tmva_input_verticesDS.root","RECREATE");
  TTree *outputvtx = new TTree("vertices","Input file for TMVA analysis");

  outputvtx->Branch("vtx_mc_ev",&vtx_mc_ev,"vtx_mc_ev/I");
  outputvtx->Branch("vtx2_mc_ev",&vtx2_mc_ev,"vtx2_mc_ev/I");
  outputvtx->Branch("vtx_id",&vtx_id,"vtx_id/I");
  outputvtx->Branch("vtx2_id",&vtx2_id,"vtx2_id/I");
  outputvtx->Branch("vtx_ntrk",&vtx_ntrk,"vtx_ntrk/I");
  outputvtx->Branch("vtx2_ntrk",&vtx2_ntrk,"vtx2_ntrk/I");
  outputvtx->Branch("vtx_maxap",&vtx_maxap,"vtx_maxap/F");
  outputvtx->Branch("vtx_bdt1",&vtx_bdt1,"vtx_bdt1/F");
  outputvtx->Branch("vtx2_bdt2",&vtx2_bdt2,"vtx2_bdt2/F");
  outputvtx->Branch("charm1",&charm1,"charm1/I");
  outputvtx->Branch("charm2",&charm2,"charm2/I");
  outputvtx->Branch("event",&event,"event/I");
  outputvtx->Branch("meannseg",&meannseg,"meannseg/F");
  outputvtx->Branch("meanfill",&meanfill,"meanfill/F");
  outputvtx->Branch("zpositive",&zpositive,"zpositive/I");
  outputvtx->Branch("samevent",&samevent,"samevent/I");
  outputvtx->Branch("charmdaughter",&charmdaughter,"charmdaughter/I");
  outputvtx->Branch("goodvtx",&goodvtx,"goodvtx/I");
  outputvtx->Branch("starting2",&starting2,"starting2/I");
  outputvtx->Branch("pointing",&pointing,"pointing/F");
  outputvtx->Branch("kink",&kink,"kink/F");
  outputvtx->Branch("impact",&impact,"impact/F");
  outputvtx->Branch("decaylength",&decaylength,"decaylength/F");
  outputvtx->Branch("dist_tr",&dist_tr,"dist_tr/F");
  outputvtx->Branch("vtx_x",&vtx_x,"vtx_x/F");
  outputvtx->Branch("vtx_y",&vtx_y,"vtx_y/F");
  outputvtx->Branch("vtx_z",&vtx_z,"vtx_z/F");
  outputvtx->Branch("vtx2_x",&vtx2_x,"vtx2_x/F");
  outputvtx->Branch("vtx2_y",&vtx2_y,"vtx2_y/F");
  outputvtx->Branch("vtx2_z",&vtx2_z,"vtx2_z/F");
  outputvtx->Branch("vtx2_prob",&vtx2_prob,"vtx2_prob/F");
  outputvtx->Branch("vtx2_maxap",&vtx2_maxap,"vtx2_maxap/F");
  outputvtx->Branch("tka_rms",&tka_rms,"tka_rms/F");
  outputvtx->Branch("tka_mean",&tka_mean,"tka_mean/F");
  outputvtx->Branch("tka_max",&tka_max,"tka_max/F");
  outputvtx->Branch("trk_pms",&trk_pms,"trk_pms/F");
  outputvtx->Branch("tkr_sizeb",&trk_sizeb,"trk_sizeb/I");
  /*outputvtx->Branch("tkr_nshower",&trk_nshower,"trk_nshower/I");
  outputvtx->Branch("tkr_shower_val15",&trk_shower_val15,"trk_shower_val15/F");
  outputvtx->Branch("tkr_shower_val30",&trk_shower_val30,"trk_shower_val30/F");
  outputvtx->Branch("tkr_shower_mean",&trk_shower_mean,"trk_shower_mean/F");
  outputvtx->Branch("tkr_shower_rms",&trk_shower_rms,"trk_shower_rms/F");
  outputvtx->Branch("tkr_shower_min",&trk_shower_min,"trk_shower_min/F");
  outputvtx->Branch("tkr_shower_max",&trk_shower_max,"trk_shower_max/F");
  outputvtx->Branch("tkr_shower",&trk_shower,"trk_shower/I");*/
  

  const Int_t nvtx = dsreader.GetEntries();
  //cout << nvertices << endl;
  for (int ivtx=0;ivtx<nvtx;ivtx++){
    dsreader.Next();
    //dsshower.Next();
    event = *Event;
    vtx_mc_ev = *Vtx_mc_ev;
    vtx_ntrk = *Vtx_ntrk;
    vtx_maxap = *Vtx_maxaperture;
    vtx_id = *Vtx_id1;
    vtx_x = *Vtx_x;
    vtx_y = *Vtx_y;
    vtx_z = *Vtx_z;
    vtx_bdt1 = *Vtx_bdt1;
    
    vector<int> ntrk = *Vtx2_ntrk;
    int sum_trk=0;
    
    
    
    for(int ivtx2=0;ivtx2<ntrk.size();ivtx2++){

      charm1=0;
      charm2=0;
      goodvtx=0;
      kink=0;
      impact=0;
      tka_rms=0;
      trk_pms=0;
      tka_mean=0;
      tka_max=0;
      meannseg=0;
      meanfill=0;
      vtx2_maxap=0;
      /*
      trk_shower_val15=0;
      trk_shower_val30=0;
      trk_nshower=0;
      trk_sizeb=0;

      trk_shower_mean=0;
      trk_shower_rms=0;
      trk_shower_min=0;
      trk_shower_max=0;
      trk_shower=0
      
      */

      float shower_info[3]={};
      float shower_val15[ntrk[ivtx2]];
      
      //float meanfillmin=1;
      //float meanfillmax=0;
      int trk_before=0;
      int ntracks=0;
      trk_shower=0;
      int index_shower15=0;
      int index_shower30=0;
      
      dist_tr = TMath::Sqrt(TMath::Power(Vtx2_x[ivtx2]-vtx_x,2)+TMath::Power(Vtx2_y[ivtx2]-vtx_y,2));

      decaylength = Decaylength[ivtx2];
      vtx2_id = Vtx_id2[ivtx2];
      vtx2_bdt2 = Vtx2_bdt2[ivtx2];
      vtx2_mc_ev = Vtx2_mc_ev[ivtx2];
      vtx2_ntrk = ntrk[ivtx2];
      vtx2_x = Vtx2_x[ivtx2];
      vtx2_y = Vtx2_y[ivtx2];
      vtx2_z = Vtx2_z[ivtx2];
      vtx2_prob = Vtx2_prob[ivtx2];
      pointing = Pointing[ivtx2];
      
      zpositive = Zpositive[ivtx2];
      samevent = Samevent[ivtx2];
      charmdaughter = Charmdaug[ivtx2];
      //goodvtx = Goodvtx[ivtx2];
      
      
      for(int k=0;k<ntrk[ivtx2];k++){
	
	nseg = Nseg[sum_trk+k];
	meannseg +=nseg;
	
	kink += ka[sum_trk+k];
	impact += ip[sum_trk+k];
	
	//cout << kink << " " << impact <<  endl;
	
	tka_rms += Tka_rms[sum_trk+k]*nseg;
	tka_mean += Tka_mean[sum_trk+k]*nseg;
	tka_max += Tka_max[sum_trk+k]*nseg;

	/*
	// SHOWER
	if(Out15[sum_trk+k]!=-10){
	  trk_shower_val15 += Out15[sum_trk+k];
	  trk_nshower += Nshower[sum_trk+k];
	  trk_sizeb += Sizeb[sum_trk+k];
	  shower_val15[index_shower15] = Out15[sum_trk+k];
	  index_shower15++;
	}
	
	if(Out30[sum_trk+k]!=-10){
	  trk_shower_val30 += Out30[sum_trk+k];
	  index_shower30++;
	}
	
	if(trk_shower_val15>0.5)trk_shower++;
	*/
	
	trk_z = Trk_z[sum_trk+k];
	if((trk_z-vtx_z)>0)goodvtx=1;
	float vtrk_pms = Trk_pms[sum_trk+k];
	if(vtrk_pms==-99)vtrk_pms=10;
	if(vtrk_pms==-1)vtrk_pms=0;
	trk_pms += vtrk_pms;
	if(Trk_pid[sum_trk+k]==1)charm1=1;
	if(Trk_pid[sum_trk+k]==2)charm2=1;
	
	float plate_step=1300;
	int remaining_plates = floor(TMath::Abs(vtx2_z)/plate_step) + 1;
	if(remaining_plates!=0){
	  trk_fill = (float) (nseg) /(float) (remaining_plates);
	}
	else trk_fill=0;
	
	meanfill += trk_fill;
	
	// MAX APERTURE
	for(int j=(k+1);j<ntrk[ivtx2];j++){
	  float dtx= TX[sum_trk+j] - TX[sum_trk+k];
	  float dty= TY[sum_trk+j] - TY[sum_trk+k];
	  
	  float tmp_ap = TMath::Sqrt(dtx*dtx+dty*dty);
	  if(tmp_ap>vtx2_maxap)vtx2_maxap=tmp_ap;
	}
	
	
	ntracks++;
      }
      
      kink /= ntracks;
      impact /= ntracks;
      meanfill /= ntracks;

            
      tka_rms /= meannseg;
      tka_mean /= meannseg;
      tka_max /= meannseg;

      meannseg /=vtx2_ntrk;

      /*
      // SHOWER
      if(index_shower15!=0){
	trk_shower_val15 /= index_shower15;
	trk_nshower /= index_shower15;
	trk_sizeb /= index_shower15;
      }
      else {
	trk_shower_val15=0;
	trk_nshower=0;
	trk_sizeb=0;
      }
      
      if(index_shower30!=0)trk_shower_val30 /= index_shower30;
      else trk_shower_val30=0;

      
      trk_shower_mean = TMath::Mean(index_shower15,shower_val15);
      trk_shower_rms = TMath::RMS(index_shower15,shower_val15);
      trk_shower_min = TMath::MinElement(index_shower15,shower_val15);
      trk_shower_max = TMath::MaxElement(index_shower15,shower_val15);
      */
      outputvtx->Fill();

      sum_trk += vtx2_ntrk;
      
    }
    
  }
  
  outputfile->Write();
  outputfile->Close();
  
}



/// \file
/// \ingroup tutorial_tmva
/// \notebook -nodraw
/// This macro provides examples for the training and testing of the
/// TMVA classifiers.
///
/// As input data is used a toy-MC sample consisting of four Gaussian-distributed
/// and linearly correlated input variables.
/// The methods to be used can be switched on and off by means of booleans, or
/// via the prompt command, for example:
///
///     root -l ./TMVAClassification.C\(\"BDT,Likelihood\"\)
///
/// (note that the backslashes are mandatory)
/// If no method given, a default set of classifiers is used.
/// The output file "TMVA.root" can be analysed with the use of dedicated
/// macros (simply say: root -l <macro.C>), which can be conveniently
/// invoked through a GUI that will appear at the end of the run of this macro.
/// Launch the GUI via the command:
///
///     root -l ./TMVAGui.C
///
/// You can also compile and run the example with the following commands
///
///     make
///     ./TMVAClassification <Methods>
///
/// where: `<Methods> = "method1 method2"` are the TMVA classifier names
/// example:
///
///     ./TMVAClassification Fisher LikelihoodPCA BDT
///
/// If no method given, a default set is of classifiers is used
///
/// - Project   : TMVA - a ROOT-integrated toolkit for multivariate data analysis
/// - Package   : TMVA
/// - Root Macro: TMVAClassification
///
/// \macro_output
/// \macro_code
/// \author Andreas Hoecker


#include <cstdlib>
#include <iostream>
#include <map>
#include <string>

#include "TChain.h"
#include "TFile.h"
#include "TTree.h"
#include "TString.h"
#include "TObjString.h"
#include "TSystem.h"
#include "TROOT.h"

#include "TMVA/Factory.h"
#include "TMVA/DataLoader.h"
#include "TMVA/Tools.h"
#include "TMVA/TMVAGui.h"

int TMVAClassification_MC_CH5( TString myMethodList = "" )
{
  // The explicit loading of the shared libTMVA is done in TMVAlogon.C, defined in .rootrc
  // if you use your private .rootrc, or run from a different directory, please copy the
  // corresponding lines from .rootrc

  // Methods to be processed can be given as an argument; use format:
  //
  //     mylinux~> root -l TMVAClassification.C\(\"myMethod1,myMethod2,myMethod3\"\)

  //---------------------------------------------------------------
  // This loads the library
  TMVA::Tools::Instance();

  // Default MVA methods to be trained + tested
  std::map<std::string,int> Use;

  // Cut optimisation
  Use["Cuts"]            = 1;
  Use["CutsD"]           = 1;
  Use["CutsPCA"]         = 0;
  Use["CutsGA"]          = 0;
  Use["CutsSA"]          = 0;
  //
  // 1-dimensional likelihood ("naive Bayes estimator")
  Use["Likelihood"]      = 1;
  Use["LikelihoodD"]     = 0; // the "D" extension indicates decorrelated input variables (see option strings)
  Use["LikelihoodPCA"]   = 1; // the "PCA" extension indicates PCA-transformed input variables (see option strings)
  Use["LikelihoodKDE"]   = 0;
  Use["LikelihoodMIX"]   = 0;
  //
  // Mutidimensional likelihood and Nearest-Neighbour methods
  Use["PDERS"]           = 1;
  Use["PDERSD"]          = 0;
  Use["PDERSPCA"]        = 0;
  Use["PDEFoam"]         = 1;
  Use["PDEFoamBoost"]    = 0; // uses generalised MVA method boosting
  Use["KNN"]             = 1; // k-nearest neighbour method
  //
  // Linear Discriminant Analysis
  Use["LD"]              = 1; // Linear Discriminant identical to Fisher
  Use["Fisher"]          = 0;
  Use["FisherG"]         = 0;
  Use["BoostedFisher"]   = 0; // uses generalised MVA method boosting
  Use["HMatrix"]         = 0;
  //
  // Function Discriminant analysis
  Use["FDA_GA"]          = 1; // minimisation of user-defined function using Genetics Algorithm
  Use["FDA_SA"]          = 0;
  Use["FDA_MC"]          = 0;
  Use["FDA_MT"]          = 0;
  Use["FDA_GAMT"]        = 0;
  Use["FDA_MCMT"]        = 0;
  //
  // Neural Networks (all are feed-forward Multilayer Perceptrons)
  Use["MLP"]             = 0; // Recommended ANN
  Use["MLPBFGS"]         = 0; // Recommended ANN with optional training method
  Use["MLPBNN"]          = 1; // Recommended ANN with BFGS training method and bayesian regulator
  Use["CFMlpANN"]        = 0; // Depreciated ANN from ALEPH
  Use["TMlpANN"]         = 0; // ROOT's own ANN
  Use["DNN_GPU"]         = 0; // CUDA-accelerated DNN training.
  Use["DNN_CPU"]         = 0; // Multi-core accelerated DNN.
  //
  // Support Vector Machine
  Use["SVM"]             = 1;
  //
  // Boosted Decision Trees
  Use["BDT"]             = 1; // uses Adaptive Boost
  Use["BDTG"]            = 0; // uses Gradient Boost
  Use["BDTB"]            = 0; // uses Bagging
  Use["BDTD"]            = 0; // decorrelation + Adaptive Boost
  Use["BDTF"]            = 0; // allow usage of fisher discriminant for node splitting
  //
  // Friedman's RuleFit method, ie, an optimised series of cuts ("rules")
  Use["RuleFit"]         = 1;
  // ---------------------------------------------------------------

  std::cout << std::endl;
  std::cout << "==> Start TMVAClassification" << std::endl;

  // Select methods (don't look at this code - not of interest)
  if (myMethodList != "") {
    for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) it->second = 0;

    std::vector<TString> mlist = TMVA::gTools().SplitString( myMethodList, ',' );
    for (UInt_t i=0; i<mlist.size(); i++) {
      std::string regMethod(mlist[i]);

      if (Use.find(regMethod) == Use.end()) {
	std::cout << "Method \"" << regMethod << "\" not known in TMVA under this name. Choose among the following:" << std::endl;
	for (std::map<std::string,int>::iterator it = Use.begin(); it != Use.end(); it++) std::cout << it->first << " ";
	std::cout << std::endl;
	return 1;
      }
      Use[regMethod] = 1;
    }
  }

  // --------------------------------------------------------------------------------------------------


   // BDT CHOICE

  int choice=0;
  
   cout << "Scegli il tipo di BDT (1: interaction vtx; 2: decay vtx; 3: DS vtx) \t";
   cin >> choice;
   cout << endl;
  
  // Here the preparation phase begins

  // Read training and test data
  // (it is also possible to use ASCII format as input -> see TMVA Users Guide)
   TFile *input(0);
   TFile *input2(0);
   TString fname="";
   TString fname2="";
   TString fname_dataset="";
   switch(choice){
   case 1:
     fname = "tmva_input_vertices.root";     
     // fname = "tmva_input_vertices_POT.root"; 
     fname_dataset = "dataset_1st";     
     break;
   case 2:
     fname = "tmva_input_vertices2nd.root";    // esclusi i primari
     fname2 = "tmva_input_vertices2nd_POT.root";    // POT
     fname_dataset = "dataset_2nd";     
     break;
   case 3:
     fname = "tmva_input_verticesDS.root";     // Decay search DS
     fname_dataset = "dataset_DS";     
     break;
   }

   if (!gSystem->AccessPathName( fname )) {
     input = TFile::Open( fname ); // check if file in local directory exists
     if(choice==2)input2 = TFile::Open( fname2 ); // check if file in local directory exists
   }
   else {
     TFile::SetCacheFileDir(".");
     input = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD");
     if(choice==2)input2 = TFile::Open("http://root.cern.ch/files/tmva_class_example.root", "CACHEREAD");
   }
   if (!input) {
     std::cout << "ERROR: could not open data file" << std::endl;
     exit(1);
   }
   std::cout << "--- TMVAClassification       : Using input file: " << input->GetName() << std::endl;
   if(choice==2){
     if (!input2) {
       std::cout << "ERROR: could not open data file" << std::endl;
       exit(1);
     }
     std::cout << "--- TMVAClassification       : Using input file: " << input2->GetName() << std::endl;

   }
   
  
  // Register the training and test trees

  TTree * sampleTree = dynamic_cast<TTree*>(input->Get("vertices"));
  TTree * sampleTree_POT =0;
  if(choice==2)sampleTree_POT = dynamic_cast<TTree*>(input2->Get("vertices"));   // DS
  TTree * signalTree =0;
  TTree * backgroundTree=0;
  TString outfileName("");

  switch(choice){
  case 1:
    outfileName = "TMVA_1st.root";
    break;
  case 2:
    outfileName = "TMVA_2nd.root";
    break;
  case 3:
    outfileName = "TMVA_DS.root";
  }
  
  // Create a ROOT output file where TMVA will store ntuples, histograms, etc.
  //TString outfileName( "TMVA.root" );
  TFile* outputFile = TFile::Open(outfileName, "RECREATE" );
 
 
  // solo primari tutto
  switch(choice){
  case 1:
    //signalTree = sampleTree->CopyTree("mp_motherID==-1");
    //backgroundTree = sampleTree->CopyTree("!(mp_motherID==-1)");
    int cut_ntrk;
    float cut_vz;
    cout << "Inserisci cut numero minimo di tracce ";
    cin >> cut_ntrk;
    cout << "Inserisci cut superiore in Vz ";
    cin >> cut_vz;
    cout << endl;
    outputFile->SetName(outfileName);    
    signalTree = sampleTree->CopyTree(Form("(((ntracks-nelectrons)>=8 || (mp_motherID==0 && ntracks>=8)) && maxfreq!=1 && (ntracks-nevents-1)>=8 && ((ntracks-maxfreq)/maxfreq<0.75 || ntracks>10))"));
    backgroundTree = sampleTree->CopyTree(Form("!(((ntracks-nelectrons)>=8 || (mp_motherID==0 && ntracks>=8)) && maxfreq!=1 && (ntracks-nevents-1)>=8 && ((ntracks-maxfreq)/maxfreq<0.75 || ntracks>10))"));
    break;
  case 2:
    // escludo primari
    float cut_bdt;
    cout << "Inserisci cut bdt di interazione ";
    cin >> cut_bdt;
    cout << endl;
    signalTree = sampleTree->CopyTree(Form("(((mp_motherID==1 || mp_motherID==2) && (bdt_value>%.2f && ntracks<=8)) || (bdt_value<%.2f && charm_daug)) && mean_rms_kink!=-1",cut_bdt,cut_bdt));
    //backgroundTree = sampleTree->CopyTree("(!(mp_motherID==1 || mp_motherID==2) && (bdt_value>0.04 && ntracks<=8) || (bdt_value<0.04 && !charm_daug)) && mean_rms_kink!=-1");
    backgroundTree = sampleTree_POT->CopyTree(Form("((bdt_value>%.2f && ntracks<=8) || (bdt_value<%.2f)) && mean_rms_kink!=-1",cut_bdt,cut_bdt));
    break;
  case 3:
    //DS
    signalTree = sampleTree->CopyTree("(zpositive && samevent && charmdaughter) && (vtx2_z-vtx_z)>500");
    backgroundTree = sampleTree->CopyTree("!(zpositive && samevent && charmdaughter) && (vtx2_z-vtx_z)>500");
    break;
  }
  
  

     

  // Create the factory object. Later you can choose the methods
  // whose performance you'd like to investigate. The factory is
  // the only TMVA object you have to interact with
  //
  // The first argument is the base of the name of all the
  // weightfiles in the directory weight/
  //
  // The second argument is the output file for the training results
  // All TMVA output can be suppressed by removing the "!" (not) in
  // front of the "Silent" argument in the option string
  TMVA::Factory *factory = new TMVA::Factory( "TMVAClassification", outputFile,
					      "!V:!Silent:Color:DrawProgressBar:Transformations=I;D;P;G,D:AnalysisType=Classification" );
  //"!V:!Silent:Color:DrawProgressBar:Transformations=I;P;G,D:AnalysisType=Classification" );

  TMVA::DataLoader *dataloader=new TMVA::DataLoader(fname_dataset);
  // If you wish to modify default settings
  // (please check "src/Config.h" to see all available global options)
  //
  //    (TMVA::gConfig().GetVariablePlotting()).fTimesRMS = 8.0;
  //    (TMVA::gConfig().GetIONames()).fWeightFileDir = "myWeightDirectory";

  // Define the input variables that shall be used for the MVA training
  // note that you may also use variable expressions, such as: "3*var1/var2*abs(var3)"
  // [all types of expressions that can also be parsed by TTree::Draw( "expression" )]
  //dataloader->AddVariable( "myvar1 := var1+var2", 'F' );
  //dataloader->AddVariable( "myvar2 := var1-var2", "Expression 2", "", 'F' );
  //dataloader->AddVariable( "var3",                "Variable 3", "units", 'F' );
  //dataloader->AddVariable( "var4",                "Variable 4", "units", 'F' );

  switch(choice){
  case 1:
    // tutto
    dataloader->AddVariable("probability", 'F');
    dataloader->AddVariable("maxaperture", 'F');   // tutto
    //dataloader->AddVariable("IP_var := meanIP/maxIP",'F');
    dataloader->AddVariable("Fill_var := vtx_fill/ntracks",'F');
    //dataloader->AddVariable("meannseg",'F');
    //dataloader->AddVariable("meanTX",'F');
    //dataloader->AddVariable("ntracks",'I');
    //dataloader->AddVariable("vtx_fill",'F');
    //dataloader->AddVariable("meanTY",'F');
    dataloader->AddVariable("maxIP",'F');
    dataloader->AddVariable("meanIP",'F');
    break;
    //
  case 2:
     // escluso primari
    dataloader->AddVariable("meannseg", 'F');    // escluso primari
    dataloader->AddVariable("bdt_value", 'F');    // escluso primari
    dataloader->AddVariable("Rms_kink :=mean_rms_kink",'F');
    dataloader->AddVariable("probability",'F');
    dataloader->AddVariable("fill:=vtx_fill*ntracks",'F');
    dataloader->AddVariable("maxaperture", 'F');   // tutto
    break;
    //
  case 3:
    //Decay search DS
    //
    dataloader->AddVariable("tka_rms", 'F');
    dataloader->AddVariable("kink",'F');
    dataloader->AddVariable("impact",'F');
    dataloader->AddVariable("decaylength",'F');
    //dataloader->AddVariable("(decaylength-dist_tr)/decaylength",'F');
    dataloader->AddVariable("(dist_tr)",'F');
    dataloader->AddVariable("vtx2_bdt2",'F');
    dataloader->AddVariable("trk_pms", 'F');
    //dataloader->AddVariable("decaylength/trk_pms", 'F');
    dataloader->AddVariable("vtx2_maxap", 'F');
    break;
    //
  }
    
    // You can add so-called "Spectator variables", which are not used in the MVA training,
    // but will appear in the final "TestTree" produced by TMVA. This TestTree will contain the
    // input variables, the response values of all trained MVAs, and the spectator variables


  if(choice!=3){
    dataloader->AddSpectator( "vID",  "Index of the vertex", "units", 'I' );
    dataloader->AddSpectator( "meanTX",  "meanTX", "units", 'F' );
    dataloader->AddSpectator( "meanTY", "meanTY","units",'F');
  }

  // global event weights per tree (see below for setting event-wise weights)
  Double_t signalWeight     = 1.0;
  Double_t backgroundWeight = 1.0;

  // You can add an arbitrary number of signal or background trees
  dataloader->AddSignalTree    ( signalTree,     signalWeight );
  dataloader->AddBackgroundTree( backgroundTree, backgroundWeight );

  // To give different trees for training and testing, do as follows:
  //
  //     dataloader->AddSignalTree( signalTrainingTree, signalTrainWeight, "Training" );
  //     dataloader->AddSignalTree( signalTestTree,     signalTestWeight,  "Test" );

  // Use the following code instead of the above two or four lines to add signal and background
  // training and test events "by hand"
  // NOTE that in this case one should not give expressions (such as "var1+var2") in the input
  //      variable definition, but simply compute the expression before adding the event
  // ```cpp
  // // --- begin ----------------------------------------------------------
  // std::vector<Double_t> vars( 4 ); // vector has size of number of input variables
  // Float_t  treevars[4], weight;
  //
  // // Signal
  // for (UInt_t ivar=0; ivar<4; ivar++) signalTree->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
  // for (UInt_t i=0; i<signalTree->GetEntries(); i++) {
  //    signalTree->GetEntry(i);
  //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
  //    // add training and test events; here: first half is training, second is testing
  //    // note that the weight can also be event-wise
  //    if (i < signalTree->GetEntries()/2.0) dataloader->AddSignalTrainingEvent( vars, signalWeight );
  //    else                              dataloader->AddSignalTestEvent    ( vars, signalWeight );
  // }
  //
  // // Background (has event weights)
  // background->SetBranchAddress( "weight", &weight );
  // for (UInt_t ivar=0; ivar<4; ivar++) background->SetBranchAddress( Form( "var%i", ivar+1 ), &(treevars[ivar]) );
  // for (UInt_t i=0; i<background->GetEntries(); i++) {
  //    background->GetEntry(i);
  //    for (UInt_t ivar=0; ivar<4; ivar++) vars[ivar] = treevars[ivar];
  //    // add training and test events; here: first half is training, second is testing
  //    // note that the weight can also be event-wise
  //    if (i < background->GetEntries()/2) dataloader->AddBackgroundTrainingEvent( vars, backgroundWeight*weight );
  //    else                                dataloader->AddBackgroundTestEvent    ( vars, backgroundWeight*weight );
  // }
  // // --- end ------------------------------------------------------------
  // ```
  // End of tree registration

  // Set individual event weights (the variables must exist in the original TTree)
  // -  for signal    : `dataloader->SetSignalWeightExpression    ("weight1*weight2");`
  // -  for background: `dataloader->SetBackgroundWeightExpression("weight1*weight2");`
  //dataloader->SetBackgroundWeightExpression( "weight" );

  // Apply additional cuts on the signal and background samples (can be different)
  TCut mycuts = ""; // for example: TCut mycuts = "abs(var1)<0.5 && abs(var2-0.5)<1";
  TCut mycutb = ""; // for example: TCut mycutb = "abs(var1)<0.5";

  // Tell the dataloader how to use the training and testing events
  //
  // If no numbers of events are given, half of the events in the tree are used
  // for training, and the other half for testing:
  //
  //    dataloader->PrepareTrainingAndTestTree( mycut, "SplitMode=random:!V" );
  //
  // To also specify the number of testing events, use:
  //
  //    dataloader->PrepareTrainingAndTestTree( mycut,
  //         "NSigTrain=3000:NBkgTrain=3000:NSigTest=3000:NBkgTest=3000:SplitMode=Random:!V" );

  switch(choice){
  case 1:
    // CUT TUTTO
    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
					    //"nTrain_Signal=1088:nTrain_Background=14426:SplitMode=Random:NormMode=NumEvents:!V" );
					    "nTrain_Signal=1154:nTrain_Background=16601:SplitMode=Random:NormMode=NumEvents:!V" );
    break;//(1042,12156)
  case 2:
  // CUT ESCLUSI PRIMARI
    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
					    "nTrain_Signal=1103:nTrain_Background=12100:SplitMode=Random:NormMode=NumEvents:!V" );
    break;
  case 3:
    
    // CUT DS
    dataloader->PrepareTrainingAndTestTree( mycuts, mycutb,
					   "nTrain_Signal=476:nTrain_Background=5038:SplitMode=Random:NormMode=NumEvents:!V" );
    break;
  }
  
  //
  // ### Book MVA methods
  //
  // Please lookup the various method configuration options in the corresponding cxx files, eg:
  // src/MethoCuts.cxx, etc, or here: http://tmva.sourceforge.net/optionRef.html
  // it is possible to preset ranges in the option string in which the cut optimisation should be done:
  // "...:CutRangeMin[2]=-1:CutRangeMax[2]=1"...", where [2] is the third input variable

  // Cut optimisation
  if (Use["Cuts"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "Cuts",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart" );

  if (Use["CutsD"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsD",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=Decorrelate" );

  if (Use["CutsPCA"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsPCA",
			 "!H:!V:FitMethod=MC:EffSel:SampleSize=200000:VarProp=FSmart:VarTransform=PCA" );

  if (Use["CutsGA"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsGA",
			 "H:!V:FitMethod=GA:CutRangeMin[0]=-10:CutRangeMax[0]=10:VarProp[1]=FMax:EffSel:Steps=30:Cycles=3:PopSize=400:SC_steps=10:SC_rate=5:SC_factor=0.95" );

  if (Use["CutsSA"])
    factory->BookMethod( dataloader, TMVA::Types::kCuts, "CutsSA",
			 "!H:!V:FitMethod=SA:EffSel:MaxCalls=150000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  // Likelihood ("naive Bayes estimator")
  if (Use["Likelihood"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "Likelihood",
			 "H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmoothBkg[1]=10:NSmooth=1:NAvEvtPerBin=50" );

  // Decorrelated likelihood
  if (Use["LikelihoodD"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodD",
			 "!H:!V:TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=Decorrelate" );

  // PCA-transformed likelihood
  if (Use["LikelihoodPCA"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodPCA",
			 "!H:!V:!TransformOutput:PDFInterpol=Spline2:NSmoothSig[0]=20:NSmoothBkg[0]=20:NSmooth=5:NAvEvtPerBin=50:VarTransform=PCA" );

  // Use a kernel density estimator to approximate the PDFs
  if (Use["LikelihoodKDE"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodKDE",
			 "!H:!V:!TransformOutput:PDFInterpol=KDE:KDEtype=Gauss:KDEiter=Adaptive:KDEFineFactor=0.3:KDEborder=None:NAvEvtPerBin=50" );

  // Use a variable-dependent mix of splines and kernel density estimator
  if (Use["LikelihoodMIX"])
    factory->BookMethod( dataloader, TMVA::Types::kLikelihood, "LikelihoodMIX",
			 "!H:!V:!TransformOutput:PDFInterpolSig[0]=KDE:PDFInterpolBkg[0]=KDE:PDFInterpolSig[1]=KDE:PDFInterpolBkg[1]=KDE:PDFInterpolSig[2]=Spline2:PDFInterpolBkg[2]=Spline2:PDFInterpolSig[3]=Spline2:PDFInterpolBkg[3]=Spline2:KDEtype=Gauss:KDEiter=Nonadaptive:KDEborder=None:NAvEvtPerBin=50" );

  // Test the multi-dimensional probability density estimator
  // here are the options strings for the MinMax and RMS methods, respectively:
  //
  //      "!H:!V:VolumeRangeMode=MinMax:DeltaFrac=0.2:KernelEstimator=Gauss:GaussSigma=0.3" );
  //      "!H:!V:VolumeRangeMode=RMS:DeltaFrac=3:KernelEstimator=Gauss:GaussSigma=0.3" );
  if (Use["PDERS"])
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERS",
			 "!H:!V:NormTree=T:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600" );

  if (Use["PDERSD"])
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSD",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=Decorrelate" );

  if (Use["PDERSPCA"])
    factory->BookMethod( dataloader, TMVA::Types::kPDERS, "PDERSPCA",
			 "!H:!V:VolumeRangeMode=Adaptive:KernelEstimator=Gauss:GaussSigma=0.3:NEventsMin=400:NEventsMax=600:VarTransform=PCA" );

  // Multi-dimensional likelihood estimator using self-adapting phase-space binning
  if (Use["PDEFoam"])
    factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoam",
			 "!H:!V:SigBgSeparate=F:TailCut=0.001:VolFrac=0.0666:nActiveCells=500:nSampl=2000:nBin=5:Nmin=100:Kernel=None:Compress=T" );

  if (Use["PDEFoamBoost"])
    factory->BookMethod( dataloader, TMVA::Types::kPDEFoam, "PDEFoamBoost",
			 "!H:!V:Boost_Num=30:Boost_Transform=linear:SigBgSeparate=F:MaxDepth=4:UseYesNoCell=T:DTLogic=MisClassificationError:FillFoamWithOrigWeights=F:TailCut=0:nActiveCells=500:nBin=20:Nmin=400:Kernel=None:Compress=T" );

  // K-Nearest Neighbour classifier (KNN)
  if (Use["KNN"])
    factory->BookMethod( dataloader, TMVA::Types::kKNN, "KNN",
			 "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

  // H-Matrix (chi2-squared) method
  if (Use["HMatrix"])
    factory->BookMethod( dataloader, TMVA::Types::kHMatrix, "HMatrix", "!H:!V:VarTransform=None" );

  // Linear discriminant (same as Fisher discriminant)
  if (Use["LD"])
    factory->BookMethod( dataloader, TMVA::Types::kLD, "LD", "H:!V:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher discriminant (same as LD)
  if (Use["Fisher"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "Fisher", "H:!V:Fisher:VarTransform=None:CreateMVAPdfs:PDFInterpolMVAPdf=Spline2:NbinsMVAPdf=50:NsmoothMVAPdf=10" );

  // Fisher with Gauss-transformed input variables
  if (Use["FisherG"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "FisherG", "H:!V:VarTransform=Gauss" );

  // Composite classifier: ensemble (tree) of boosted Fisher classifiers
  if (Use["BoostedFisher"])
    factory->BookMethod( dataloader, TMVA::Types::kFisher, "BoostedFisher",
			 "H:!V:Boost_Num=20:Boost_Transform=log:Boost_Type=AdaBoost:Boost_AdaBoostBeta=0.2:!Boost_DetailedMonitoring" );

  // Function discrimination analysis (FDA) -- test of various fitters - the recommended one is Minuit (or GA or SA)
  if (Use["FDA_MC"])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MC",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:SampleSize=100000:Sigma=0.1" );

  if (Use["FDA_GA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:PopSize=100:Cycles=2:Steps=5:Trim=True:SaveBestGen=1" );

  if (Use["FDA_SA"]) // can also use Simulated Annealing (SA) algorithm (see Cuts_SA options])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_SA",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=SA:MaxCalls=15000:KernelTemp=IncAdaptive:InitialTemp=1e+6:MinTemp=1e-6:Eps=1e-10:UseDefaultScale" );

  if (Use["FDA_MT"])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=2:UseImprove:UseMinos:SetBatch" );

  if (Use["FDA_GAMT"])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_GAMT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=GA:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:Cycles=1:PopSize=5:Steps=5:Trim" );

  if (Use["FDA_MCMT"])
    factory->BookMethod( dataloader, TMVA::Types::kFDA, "FDA_MCMT",
			 "H:!V:Formula=(0)+(1)*x0+(2)*x1+(3)*x2+(4)*x3:ParRanges=(-1,1);(-10,10);(-10,10);(-10,10);(-10,10):FitMethod=MC:Converger=MINUIT:ErrorLevel=1:PrintLevel=-1:FitStrategy=0:!UseImprove:!UseMinos:SetBatch:SampleSize=20" );

  // TMVA ANN: MLP (recommended ANN) -- all ANNs in TMVA are Multilayer Perceptrons
  if (Use["MLP"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

  if (Use["MLPBFGS"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

  if (Use["MLPBNN"])
    factory->BookMethod( dataloader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators


  // Multi-architecture DNN implementation.
  if (Use["DNN_CPU"] or Use["DNN_GPU"]) {
    // General layout.
    TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

    // Training strategies.
    TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
		      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
		      "WeightDecay=1e-4,Regularization=L2,"
		      "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
    TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
		      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
		      "WeightDecay=1e-4,Regularization=L2,"
		      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
		      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
		      "WeightDecay=1e-4,Regularization=L2,"
		      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString trainingStrategyString ("TrainingStrategy=");
    trainingStrategyString += training0 + "|" + training1 + "|" + training2;

    // General Options.
    TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
			"WeightInitialization=XAVIERUNIFORM");
    dnnOptions.Append (":"); dnnOptions.Append (layoutString);
    dnnOptions.Append (":"); dnnOptions.Append (trainingStrategyString);

    // Cuda implementation.
    if (Use["DNN_GPU"]) {
      TString gpuOptions = dnnOptions + ":Architecture=GPU";
      factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_GPU", gpuOptions);
    }
    // Multi-core CPU implementation.
    if (Use["DNN_CPU"]) {
      TString cpuOptions = dnnOptions + ":Architecture=CPU";
      factory->BookMethod(dataloader, TMVA::Types::kDNN, "DNN_CPU", cpuOptions);
    }
  }

  // CF(Clermont-Ferrand)ANN
  if (Use["CFMlpANN"])
    factory->BookMethod( dataloader, TMVA::Types::kCFMlpANN, "CFMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N"  ); // n_cycles:#nodes:#nodes:...

  // Tmlp(Root)ANN
  if (Use["TMlpANN"])
    factory->BookMethod( dataloader, TMVA::Types::kTMlpANN, "TMlpANN", "!H:!V:NCycles=200:HiddenLayers=N+1,N:LearningMethod=BFGS:ValidationFraction=0.3"  ); // n_cycles:#nodes:#nodes:...

  // Support Vector Machine
  if (Use["SVM"])
    factory->BookMethod( dataloader, TMVA::Types::kSVM, "SVM", "Gamma=0.25:Tol=0.001:VarTransform=Norm" );

  // Boosted Decision Trees
  if (Use["BDTG"]) // Gradient Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTG",
			 "!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

  if (Use["BDT"])  // Adaptive Boost
      factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
  		 "!H:!V:NTrees=100:MinNodeSize=2.5%:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=40" );

  // Boosted Decision Trees with adaptive boosting
  //if (Use["BDT"])
  //factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDT",
  //			 "!H:!V:NTrees=400:nEventsMin=0.05:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:PruneMethod=NoPruning" );
  
  if (Use["BDTB"]) // Bagging
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTB",
			 "!H:!V:NTrees=400:BoostType=Bagging:SeparationType=GiniIndex:nCuts=20" );

  if (Use["BDTD"]) // Decorrelation + Adaptive Boost
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTD",
			 "!H:!V:NTrees=400:MinNodeSize=2%:MaxDepth=3:BoostType=AdaBoost:SeparationType=GiniIndex:nCuts=20:VarTransform=Decorrelate" );

  if (Use["BDTF"])  // Allow Using Fisher discriminant in node splitting for (strong) linearly correlated variables
    factory->BookMethod( dataloader, TMVA::Types::kBDT, "BDTF",
			 "!H:!V:NTrees=400:MinNodeSize=2.5%:UseFisherCuts:MaxDepth=3:BoostType=AdaBoost:AdaBoostBeta=0.5:SeparationType=GiniIndex:nCuts=20" );

  // RuleFit -- TMVA implementation of Friedman's method
  if (Use["RuleFit"])
    factory->BookMethod( dataloader, TMVA::Types::kRuleFit, "RuleFit",
			 "H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

  // For an example of the category classifier usage, see: TMVAClassificationCategory
  //
  // --------------------------------------------------------------------------------------------------
  //  Now you can optimize the setting (configuration) of the MVAs using the set of training events
  // STILL EXPERIMENTAL and only implemented for BDT's !
  //
  //     factory->OptimizeAllMethods("SigEffAt001","Scan");
  //     factory->OptimizeAllMethods("ROCIntegral","FitGA");
  //
  // --------------------------------------------------------------------------------------------------

  // Now you can tell the factory to train, test, and evaluate the MVAs
  //
  // Train MVAs using the set of training events
  factory->TrainAllMethods();

  // Evaluate all MVAs using the set of test events
  factory->TestAllMethods();

  // Evaluate and compare performance of all configured MVAs
  factory->EvaluateAllMethods();

  // --------------------------------------------------------------

  // Save the output
  outputFile->Close();

  std::cout << "==> Wrote root file: " << outputFile->GetName() << std::endl;
  std::cout << "==> TMVAClassification is done!" << std::endl;

  delete factory;
  delete dataloader;
  // Launch the GUI for the root macros
  if (!gROOT->IsBatch()) TMVA::TMVAGui( outfileName );

  return 0;
}

int main( int argc, char** argv )
{
  // Select methods (don't look at this code - not of interest)
  TString methodList;
  for (int i=1; i<argc; i++) {
    TString regMethod(argv[i]);
    if(regMethod=="-b" || regMethod=="--batch") continue;
    if (!methodList.IsNull()) methodList += TString(",");
    methodList += regMethod;
  }
  return TMVAClassification_MC_CH5(methodList);
}
