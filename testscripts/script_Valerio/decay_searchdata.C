// Classificazione dei vertici simulati ricostruiti in FEDRA
// Training per analisi multivariata
// V. Gentile 2019

//#include "GiuliDecaySearch.h"
#include "Definitions.h"
#include <TH2.h>
#include <TH3F.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>
#include <iomanip>
#include <locale>
#include <sstream>
#include <string>
#include <fstream>
#include <algorithm>
#include <numeric>
#include <stdio.h>
#include <TF1.h>
#include <TStyle.h>

#include <TGraph.h>
#include <TGraph2D.h>
#include <TMultiGraph.h>
#include <TMath.h>
#include <TRandom3.h>
#include <TSystem.h>
#include <TLine.h>
#include <TGraph.h>
#include <TF1.h>
#include <TTree.h>
#include <TFile.h>
#include <TTreeReader.h>
#include <TTreeReaderArray.h>
#include <TObjArray.h>

#include <TParticlePDG.h>
#include <TPDGCode.h>
#include <TDatabasePDG.h>

#include <EdbVertex.h>

using namespace std;


struct Vertex
{
  Int_t event;   // index of the tree
  Int_t flag;    // vertex type: -1 fake; 1 primary; 0 charm; 2 secondary;
  Int_t fe_id;  // fedra_id
  Int_t mc_id;  // monte carlo id
  Int_t mc_pdg;  // monte carlo pdg
  Int_t fe_pproton_id;  // fedra primary proton vertex id linked to the mc id 
  Int_t ntrk;   // number of tracks attached to the vertex
  Int_t mc_ch_id1; // monte carlo charm id 1
  Int_t mc_ch_id2; // monte carlo charm id 2
  Int_t mc_nd[2] ; // number montecarlo charged daugthers 
  Double_t mc_dz[2];   // delta z between charm and primary vertex
  Double_t mc_dl[2];   // decay length between charm and primary vertex
};
  
struct Track
{
  vector<int> icharm;
  vector<int> fe_itrk;              // number fedra charged daugthers 
  vector<int>  fe_id;
  vector<int>  mc_id;
  vector<int>  pdg;               // track pdg
  vector<double>  fe_dz;          // starting coordinates of the track
  vector<double> fe_dl;           // 2d angle projections in xz and yz planes
  vector<double> ip;              // impact parameter
  vector<double> ka;              // kink angle
  vector<int> iseg_rmax;          // i segment of max kink over rms kink of a track
  vector<float> rmax;             // max kink over rms kink of a track
  vector<int> iseg_ipmax;         // i segment of max ip over rms ip of a track
  vector<float> ipmax;            // max ip over rms ip of a track
  vector<int> nseg;               // number of segments
  vector<int> type;               // 0 - track alone, 1 - in primary vtx, 2 - in other vtx same event, 3 - in other vtx, 4 - fake vertex
  vector<bool> split;             // true if there are two tracks with the same mc_trk_id and same mc_event_id
  vector<bool> pardaug;           // true if a track has segments with different mc_trk_id
};

Vertex vtx;
Track trk;

void FillVertex(int ievent, int iflag, int ife_id, int imc_id, int imc_pdg, int ipproton, int intrk, int icharm_id1, int icharm_id2, int imc_nd1, int imc_nd2, double imc_dz1, double imc_dz2, double imc_dl1, double imc_dl2){
  vtx.event = ievent;
  vtx.flag = iflag;
  vtx.fe_id = ife_id;
  vtx.mc_id = imc_id;
  vtx.mc_pdg = imc_pdg;
  vtx.fe_pproton_id= ipproton;
  vtx.ntrk = intrk;
  vtx.mc_ch_id1 = icharm_id1;
  vtx.mc_ch_id2 = icharm_id2;
  vtx.mc_nd[0] = imc_nd1;
  vtx.mc_nd[1] = imc_nd2;
  vtx.mc_dz[0] = imc_dz1;
  vtx.mc_dz[1] = imc_dz2;
  vtx.mc_dl[0] = imc_dl1;
  vtx.mc_dl[1] = imc_dl2;
}

void FillTrack(int ife_itrk, int iicharm, int ife_id, int imc_id, int ipdg, int inseg,  int itype, double ife_dz, double ife_dl, double iip, double ika, float isrmax, float irmax, float isipmax, float iipmax, bool isplit, bool ipardaug){
  trk.fe_itrk.push_back(ife_itrk);
  trk.icharm.push_back(iicharm);
  trk.fe_id.push_back(ife_id);
  trk.mc_id.push_back(imc_id);
  trk.pdg.push_back(ipdg);
  trk.fe_dz.push_back(ife_dz);
  trk.fe_dl.push_back(ife_dl);
  trk.ip.push_back(iip);
  trk.ka.push_back(ika);
  trk.iseg_rmax.push_back(isrmax);
  trk.rmax.push_back(irmax);
  trk.iseg_ipmax.push_back(isipmax);
  trk.ipmax.push_back(iipmax);
  trk.pardaug.push_back(ipardaug);
  trk.nseg.push_back(inseg);
  trk.type.push_back(itype);
  trk.split.push_back(isplit);
}

void Set0() {
  vtx.event=0, vtx.flag=0, vtx.fe_id=0, vtx.mc_id=0, vtx.mc_pdg=0; vtx.ntrk=0, vtx.mc_ch_id1=0, vtx.mc_ch_id2=0, vtx.mc_nd[0]=0, vtx.mc_nd[1]=0, vtx.mc_dz[0]=0, vtx.mc_dz[1]=0, vtx.mc_dl[0]=0, vtx.mc_dl[1]=0, vtx.fe_pproton_id=0;
  //vtx.mc_nd=0, vtx.mc_dz1=0, vtx.mc_dz2=0, vtx.mc_dl1=0, vtx.mc_dl2=0;
  trk.fe_itrk.clear(), trk.icharm.clear(), trk.fe_id.clear(), trk.mc_id.clear(), trk.pdg.clear(), trk.fe_dz.clear(), trk.fe_dl.clear(), trk.ip.clear(), trk.ka.clear(), trk.nseg.clear(), trk.type.clear(), trk.iseg_rmax.clear(), trk.rmax.clear(), trk.iseg_ipmax.clear(), trk.ipmax.clear(), trk.split.clear(), trk.pardaug.clear();
}

float IPtoVertex(TVector3 vertexpos, TVector3 trackstartpos, float tracktx, float trackty){
  //'''Impact parameter of track with respect to primary vertex'''
 
  float dz = -vertexpos(2) + trackstartpos(2);
  float ipx = -tracktx * dz + trackstartpos(0) - vertexpos(0);
  float ipy = -trackty * dz + trackstartpos(1) - vertexpos(1);

  float ip = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));

  //cout <<"pos "<< vertexpos(2) << " " << trackstartpos(2) << " " << tracktx << " " << trackty << endl;
  //cout << trackstartpos(0) << " " << trackstartpos(1) << endl;
  //cout <<"ip "<< dz << " " << ipx << " " << ipy << " " << ip << endl;

  return ip;
}


void SegIPtoVertex(TVector3 vertexpos, EdbTrackP *mytrack, float IpSeg[2]){
  //'''Impact parameter of track with respect to primary vertex'''

  int nseg = mytrack->N();
  float ipseg[nseg];
  float delta_ipseg[nseg-1];
  float ipseg_nomax[nseg-2];
  for (int i = 0; i < nseg; i++){
    EdbSegP *aseg = (EdbSegP*) mytrack->GetSegment(i);
    
    float dz = -vertexpos(2) + aseg->Z();
    float ipx = -aseg->TX() * dz + aseg->X() - vertexpos(0);
    float ipy = -aseg->TY() * dz + aseg->Y() - vertexpos(1);
    ipseg[i] = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));
  }

  for (int i = 0; i < nseg-1; i++){
    delta_ipseg[i] = ipseg[i+1] - ipseg[i];
  }
  
  float ipsegmax = TMath::MaxElement(nseg-1, delta_ipseg);
  int index=0;
  for (int i = 0; i < nseg-2; i++){
    if(delta_ipseg[i]!=ipsegmax){
      ipseg_nomax[index]=delta_ipseg[i];
      index++;
    }
    else IpSeg[0]=i;
  }
  IpSeg[0]++;
  float ipsegrms = TMath::RMS(index, ipseg_nomax);
  IpSeg[1] = ipsegmax/ipsegrms;
}

void SegIPtoVertex2(TVector3 vertexpos, int nseg, double* seg_x, double* seg_tx, double* seg_y, double* seg_ty, double* seg_z, float IpSeg[2], int trk_id){
  //'''Impact parameter of track with respect to primary vertex'''

  float ipseg[nseg];
  float delta_ipseg[nseg-1];
  float ipseg_nomax[nseg-2];
  for (int i = 0; i < nseg; i++){   
    float dz = -vertexpos(2) + seg_z[i];
    float ipx = -seg_tx[i] * dz + seg_x[i] - vertexpos(0);
    float ipy = -seg_ty[i] * dz + seg_y[i] - vertexpos(1);
    ipseg[i] = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));
    if(trk_id==34432)cout <<"ip "<< nseg << " " <<  i <<" " <<  ipseg[i] << endl;
  }

  for (int i = 0; i < nseg-1; i++){
    delta_ipseg[i] = ipseg[i+1] - ipseg[i];
    if(trk_id==34432)cout <<"delta ip "<< nseg << " " <<  i <<" " <<  delta_ipseg[i] << endl;
  }

  float ipsegmax = TMath::MaxElement(nseg-1, delta_ipseg);
  int index=0;
  IpSeg[0]=0;
  for (int i = 0; i < nseg-2; i++){
    if(delta_ipseg[i]!=ipsegmax){
      ipseg_nomax[index]=delta_ipseg[i];
      index++;
    }
    else IpSeg[0]=i;
    if(trk_id==34432)cout <<"nomax ip "<< nseg << " " <<  i <<" " <<  delta_ipseg[i] << " " << ipsegmax << " " <<index << " " <<  ipseg_nomax[index-1] << " " << IpSeg[0] << endl;
  }
  IpSeg[0]++;
  float ipsegrms = TMath::RMS(index, ipseg_nomax);
  IpSeg[1] = ipsegmax/ipsegrms;
  if(trk_id==34432)cout << "finale " << IpSeg[0] << " " << IpSeg[1] << endl;
}

float AverageKinkAngle(float parenttx, float parentty,const float* daughterstx, const float* daughtersty, int ndaughters){
  float kink=0.;

  //loop on daughters
  for (int i = 0; i < ndaughters; i++){
 
    float daughtertx = daughterstx[i];
    float daughterty = daughtersty[i];

    kink += TMath::Sqrt(pow(parenttx - daughtertx,2) + pow(parentty - daughterty,2));
  }
 
  return kink/ndaughters;
}


float KinkAngle(float parenttx, float parentty, float daughtertx, float daughterty){
  //only one daughter, average kinkangle equal to kink angle
  return AverageKinkAngle(parenttx,parentty,&daughtertx,&daughterty,1);

}

void FedraTrackKink(EdbTrackP* mytrack, float dtheta[2]){

  int nseg = mytrack->N();
  float kinkangles[nseg-1];
  float kinkangles_nomax[nseg-2];
  //loop on subsequent segment pairs of the tracks
  //cout << "num seg " << nseg << endl;
  for (int i = 0; i < nseg-1; i++){
    EdbSegP *firstseg = (EdbSegP*) mytrack->GetSegment(i);
    EdbSegP *secondseg = (EdbSegP*) mytrack->GetSegment(i+1);
    kinkangles[i]=TMath::Sqrt(pow(secondseg->TX()-firstseg->TX(),2)+pow(secondseg->TY()-firstseg->TY(),2));
    // cout <<"kink "<< nseg << " " <<  i <<" " <<  kinkangles[i] << endl;
  }
  //getting maximum and rms
  float deltathetamax = TMath::MaxElement(nseg-1, kinkangles);
  int index=0;
  for (int i = 0; i < nseg-2; i++){
    if(kinkangles[i]!=deltathetamax){
      kinkangles_nomax[index]=kinkangles[i];
      index++;
    }
    else dtheta[0]=i;
  }
  dtheta[0]++;
  float deltathetarms = TMath::RMS(index, kinkangles_nomax);
  dtheta[1] = deltathetamax/deltathetarms;
  //return rmax;

}

void FedraTrackKink2(int nseg, double* seg_tx, double* seg_ty, float dtheta[2],int trk_id){

  float kinkangles[nseg-1];
  float kinkangles_nomax[nseg-2];
  //loop on subsequent segment pairs of the tracks
  //cout << "num seg " << nseg << endl;
  for (int i = 0; i < nseg-1; i++){
    kinkangles[i]=TMath::Sqrt(pow(seg_tx[i+1]-seg_tx[i],2)+pow(seg_ty[i+1]-seg_ty[i],2));
    if(trk_id==34433)cout <<"kink "<< nseg << " " << i << " " << kinkangles[i] << endl;
  }
  //getting maximum and rms
  float deltathetamax = TMath::MaxElement(nseg-1, kinkangles);
  int index=0;
  for (int i = 0; i < nseg-2; i++){
    if(kinkangles[i]!=deltathetamax){
      kinkangles_nomax[index]=kinkangles[i];
      index++;
    }
    else dtheta[0]=i;
  }
  dtheta[0]++;
  float deltathetarms = TMath::RMS(index, kinkangles_nomax);
  dtheta[1] = deltathetamax/deltathetarms;
  if(trk_id==34433) cout << "final "<<deltathetamax << " " << deltathetarms << " " << dtheta[1] << endl;

}


bool ParentDaugther(EdbTrackP* mytrack){

  int nseg = mytrack->N();
  //float kinkangles[nseg-1];
  //loop on subsequent segment pairs of the tracks
  //cout << "num seg " << nseg << endl;

  bool pardaug=false;
 
  for (int i = 0; i < nseg-1; i++){
    EdbSegP *firstseg = (EdbSegP*) mytrack->GetSegment(i);
    EdbSegP *secondseg = (EdbSegP*) mytrack->GetSegment(i+1);

    int pre_trkid = firstseg->MCTrack();
    int post_trkid = secondseg->MCTrack();
    int pre_evid = firstseg->MCEvt();
    int post_evid = secondseg->MCEvt();

    if(pre_evid==post_evid && pre_trkid!=post_trkid)pardaug=true;
    //cout << pre_trkid << " " << post_trkid << " " << pre_evid << " " << post_evid << " " << pardaug << endl;
  }

  return pardaug;

}

bool ParentDaugther2(int nseg, double* seg_mc_trk, double* seg_mc_evt){

  bool pardaugh=false;
 
  for (int i = 0; i < nseg-1; i++){
   
    int pre_trkid = seg_mc_trk[i];
    int post_trkid = seg_mc_trk[i+1];
    int pre_evid = seg_mc_evt[i];
    int post_evid = seg_mc_evt[i+1];

    if(pre_evid==post_evid && pre_trkid!=post_trkid)pardaugh=true;
  
  }

  return pardaugh;

}



void Loop()
{
  
  ofstream log_vtx;
  log_vtx.open("log_vtx_search.txt",ios::out);
  log_vtx << "//----V. GENTILE 10/03/2019 ------// \n VERTEX SEARCH ON MC SIMULATION" << endl;

  ofstream log_decay;
  log_decay.open("log_decay_search.txt",ios::out);
  log_decay << "//----V. GENTILE 01/06/2019 ------// \n DECAY SEARCH ON MC SIMULATION" << endl;
  
  
  TFile * simulationfile = TFile::Open("../ship.conical.Pythia8CharmOnly-TGeant4.root");  
  TFile * fedrafile = TFile::Open("vertices_MC_modified.root");  // full
  TFile * trackfile = TFile::Open("b000001.0.0.0.trk.root");  // full

  EdbVertexRec *vertexrec = (EdbVertexRec*)fedrafile->Get("EdbVertexRec");
  EdbVertex *vertexobject = 0;
  EdbVertex *vertexobject_tmp = 0;
  EdbTrackP *ftrack = 0;
  
  //int nv = vertexrec->Nvtx();
  //cout << nv << endl;
  
  //TTreeReader tracks_ttreader("tracks",fedrafile);
  //TTreeReaderValue<float> trackX(tracks_ttreader, "t.eX");
  //TTreeReaderArray<int> trk_mID(tracks_ttreader, "t.eAid[2]");
  //tracks_ttreader.Next();
  //cout<<"PROVA READER: "<<*trackX<<endl;

 
  
  TTree *simtree = (TTree*) simulationfile->Get("cbmsim");
  TTree *vtx_fedratree = (TTree*) fedrafile->Get("vtx");
  TTree *trk_fedratree = (TTree*) trackfile->Get("tracks");
  //TTree *trk_fedratree = (TTree*) fedrafile->Get("tracks");

  if(vtx_fedratree == 0) return;
  cout << "HELLO WORLD" << endl;
  if(trk_fedratree == 0) return;
  if(simtree == 0) return;

  cout << "HELLO WORLD" << endl;
    
  vtx_reader_Fedra(vtx_fedratree);
  tracks_reader_Fedra(trk_fedratree);
  cbmsim_reader_MC(simtree);

  Long64_t nentries_vtx_fedra = vtx_fedratree->GetEntriesFast();
  log_decay<<"tot entries vtx fedra "<< nentries_vtx_fedra << endl;
  log_vtx<<"tot entries vtx fedra "<< nentries_vtx_fedra << endl;
  
  Long64_t nentries_trk_fedra = trk_fedratree->GetEntriesFast();
  log_decay<<"tot entries trk fedra "<< nentries_trk_fedra << endl;
  
  Long64_t nentries_MC = simtree->GetEntriesFast();
  log_decay<<"tot entries MC "<< nentries_MC << endl;
  log_vtx<<"tot entries MC "<< nentries_MC << endl;

  vtx_fedratree->BuildIndex("vID");
  //simtree->BuildIndex("MCTrack.fStartZ");

  const int MCevents=2500;
  
  int fa,fb,fc,fe,ff,fg=0;
  double fd,fh,fi,fl;
  
  std::vector<std::vector<int> > fEventId;
  std::vector<std::vector<int> > fPdg;
  std::vector<std::vector<int> > fTrkId;
  std::vector<std::vector<int> > fMPdg;
  std::vector<std::vector<int> > fMTrkId;
  std::vector<std::vector<int> > fMTrkId_true;
  std::vector<std::vector<double> > fMom;
  std::vector<std::vector<double> > fSX;
  std::vector<std::vector<double> > fSY;
  std::vector<std::vector<double> > fSZ;
  std::vector<std::vector<double> > c1momMC;
  std::vector<std::vector<double> > c2momMC;

  std::vector<double> pvxMC;
  std::vector<double> pvyMC;
  std::vector<double> pvzMC;
  std::vector<double> c1vxMC;
  std::vector<double> c1vyMC;
  std::vector<double> c1vzMC;
  std::vector<double> c2vxMC;
  std::vector<double> c2vyMC;
  std::vector<double> c2vzMC;
  //std::vector<double> c1momMC;
  //std::vector<double> c2momMC;
  std::vector<double> charm_1_z;
  std::vector<double> charm_2_z;
  std::vector<double> charm_1_3d;
  std::vector<double> charm_2_3d;

  std::vector<std::vector<int> > charm_1_id;
  std::vector<std::vector<int> > charm_2_id;
  std::vector<std::vector<int> > daug_1_id;
  std::vector<std::vector<int> > daug_2_id;
 
  pvxMC.resize(MCevents);
  pvyMC.resize(MCevents);
  pvzMC.resize(MCevents);
  c1vxMC.resize(MCevents);
  c1vyMC.resize(MCevents);
  c1vzMC.resize(MCevents);
  c2vxMC.resize(MCevents);
  c2vyMC.resize(MCevents);
  c2vzMC.resize(MCevents);
  c1momMC.resize(MCevents);
  c2momMC.resize(MCevents);
  charm_1_z.resize(MCevents);
  charm_2_z.resize(MCevents);
  charm_1_3d.resize(MCevents);
  charm_2_3d.resize(MCevents);
  charm_1_id.resize(MCevents);
  charm_2_id.resize(MCevents);
  daug_1_id.resize(MCevents);
  daug_2_id.resize(MCevents);
  
  fEventId.resize(MCevents);
  fPdg.resize(MCevents);
  fTrkId.resize(MCevents);
  fMPdg.resize(MCevents);
  fMTrkId.resize(MCevents);
  fMTrkId_true.resize(MCevents);
  fMom.resize(MCevents);
  fSX.resize(MCevents);
  fSY.resize(MCevents);
  fSZ.resize(MCevents);
			    
  
  ifstream input;
  int iline=0;
  int tmp_fa=0;
  bool found_primary_vtx=false;
  bool found_charm_vtx1=false;
  bool found_charm_vtx2=false;

  TH1D * hcharm1_z = new TH1D("hcharm1_z","",100,0,20000);
  TH1D * hcharm2_z = new TH1D("hcharm2_z","",100,0,20000);
  TH1D * hcharm_z = new TH1D("hcharm_z","",100,0,20000);

  TH1D * hcharm1_3d = new TH1D("hcharm1_3d","",100,0,20000);
  TH1D * hcharm2_3d = new TH1D("hcharm2_3d","",100,0,20000);
  TH1D * hcharm_3d = new TH1D("hcharm_3d","",100,0,20000);

  TH1D * hcharm1_mom = new TH1D("hcharm1_mom","",500,0,500);
  TH1D * hcharm2_mom = new TH1D("hcharm2_mom","",500,0,500);


  // MC FILES
  
  input.open("dumpinfo.txt",ios::in);
  while(1){
    //input >> fa >> fb >> fc >> fd >> fe >> ff >> fg >> fh >> fi >> fl;
    input >> fa >> fb >> fc >> fd >> fe >> ff  >> fg >> fh >> fi >> fl;
    if(input.eof())break;
    if(!input.good())break;
    if(tmp_fa!=fa){
      iline=0;
      found_primary_vtx=false;
      found_charm_vtx1=false;
      found_charm_vtx2=false;
    }
    
    fEventId[fa].push_back(fa);
    fTrkId[fa].push_back(fb);
    fPdg[fa].push_back(fc);    
    fMom[fa].push_back(fd);
    fMPdg[fa].push_back(fe);
    fMTrkId[fa].push_back(ff);
    fMTrkId_true[fa].push_back(fg);
    fSX[fa].push_back(fh);
    fSY[fa].push_back(fi);
    fSZ[fa].push_back(fl);
    tmp_fa=fa;
    
    if(!found_primary_vtx && fPdg[fa][iline]==2212 && fTrkId[fa][iline]==0){
      found_primary_vtx=true;
      pvxMC[fa]=fSX[fa][iline]; //um
      pvyMC[fa]=fSY[fa][iline]; //um
      pvzMC[fa]=fSZ[fa][iline]; //um
      
    }


    if(fMTrkId_true[fa][iline]==1){
      
      if(!found_charm_vtx1){
	found_charm_vtx1=true;
	c1vxMC[fa]=fSX[fa][iline]; //um
	c1vyMC[fa]=fSY[fa][iline]; //um
	c1vzMC[fa]=fSZ[fa][iline]; //um
	// cout <<"charm1 " << fa << " " << c1momMC[fa] << " " << fTrkId[fa][iline] << endl;  
      }
      
      simtree->GetEntry(fa);
      int daug_index=0;
      for (Long64_t ientry=0; ientry<MCTrack_;ientry++) {
	TParticlePDG *pdgPart = (TDatabasePDG::Instance())->GetParticle(fc);
	Double_t charge  = pdgPart->Charge();
	int mycharge = abs(int(charge))/3.;
	if(ientry==fb && abs(fc)<10000 && charge!=0){
	  if(fa==521)cout << fa << " " << fb << " " << fc << " " << charge << endl;
	  charm_1_id[fa].push_back(MCTrack_fMotherId[ientry]);
	  daug_1_id[fa].push_back(fb);
	  if(fa==521)cout << "charm1 moth " << " " << fb << " " << charm_1_id[fa][daug_index] <<  endl;
	  c1momMC[fa].push_back(fMom[fa][iline]);
	  daug_index++;
	}
      }
    }
    
    if(fMPdg[fa][iline]!=1 && fMTrkId[fa][iline]>1 &&  (abs(fMPdg[fa][iline])==411 || abs(fMPdg[fa][iline])==421 || abs(fMPdg[fa][iline])==431 || abs(fMPdg[fa][iline]==4122))){
      if(!found_charm_vtx2){
	found_charm_vtx2=true;
	c2vxMC[fa]=fSX[fa][iline]; //um
	c2vyMC[fa]=fSY[fa][iline]; //um
	c2vzMC[fa]=fSZ[fa][iline]; //um
	
	
	charm_1_z[fa] = c1vzMC[fa]-pvzMC[fa];
	charm_2_z[fa] = c2vzMC[fa]-pvzMC[fa];
	charm_1_3d[fa] = sqrt(pow(c1vxMC[fa]-pvxMC[fa],2)+pow(c1vyMC[fa]-pvyMC[fa],2)+pow(c1vzMC[fa]-pvzMC[fa],2));
	charm_2_3d[fa] = sqrt(pow(c2vxMC[fa]-pvxMC[fa],2)+pow(c2vyMC[fa]-pvyMC[fa],2)+pow(c2vzMC[fa]-pvzMC[fa],2));
	
	hcharm1_z->Fill(charm_1_z[fa]);
	hcharm2_z->Fill(charm_2_z[fa]);
	hcharm1_3d->Fill(charm_1_3d[fa]);
	hcharm2_3d->Fill(charm_2_3d[fa]);
	hcharm_z->Fill(charm_1_z[fa]);
	hcharm_z->Fill(charm_2_z[fa]);
	hcharm_3d->Fill(charm_1_3d[fa]);
	hcharm_3d->Fill(charm_2_3d[fa]);
	
	//cout <<"charm2 " << fa << " " << c2momMC[fa] << " " << fTrkId[fa][iline] << endl;  	
      }
      simtree->GetEntry(fa);
      int daug_index=0;
      for (Long64_t ientry=0; ientry<MCTrack_;ientry++) {

	TParticlePDG *pdgPart = (TDatabasePDG::Instance())->GetParticle(fc);
	Double_t charge  = pdgPart->Charge();
	int mycharge = abs(int(charge))/3.;
	
	if(ientry==fb && abs(fc)<10000 && charge!=0){
	  charm_2_id[fa].push_back(MCTrack_fMotherId[ientry]);
	  daug_2_id[fa].push_back(fb);
	  if(fa==521)cout << "charm2 moth " << " " << fb << " " << charm_2_id[fa][daug_index] << endl;
	  c2momMC[fa].push_back(fMom[fa][iline]);
	  daug_index++;
	}
      }      
    }
    //if(fa==521)cout << fa << " " <<  fb << " " <<  fc << " " <<  fd << " " <<  fe << " " <<  ff << " " <<  fg << " " <<  fh << " " <<  fi << endl;
    iline++;
  }


  
  std::vector<bool> event_primary_proton;
  std::vector<int> vtx_trackmother;
  std::vector<int> vtx_trackprocess;
  std::vector<int> vtx_trackpdgcode;
  std::vector<int> vtx_motherpdgcode;
  std::vector<int> vtx_trackid;
  std::vector<double> vtx_trackstartX;
  std::vector<double> vtx_trackstartY;
  std::vector<double> vtx_trackstartZ;
  std::vector<double> vtx_trackeventId;
  std::vector<double> vtx_trackmom;
  
  float mean_freq=0;
  int max_freq=0;
  float mean_ev_freq=0;
  int max_ev_freq=0;
  int mp_motherID=0;
  int mp_eventID=0;
  int mp_procID=0;
  int mp_pdgID=0;
  int mp_ev_pdgID=0;
  int n_electrons=0;

  bool flag=0;
  bool flag_ev=0;
  //bool evt_good=false;
  //bool daughter_proton=false;
  bool vtx_good=false;
  

  float mp_vx=0;
  float mp_vy=0;
  float mp_vz=0;
  float mp_vdx=-1;
  float mp_vdy=-1;
  float mp_vdz=-1;



  TH1F *hoccurrence = new TH1F("occurrence / ntracks","",50,0,1);

  int vbinx=13;
  int vbiny=6;
  int vbinz=4;
  
  TH3D * hvertex = new TH3D("hv","hv",vbinx,0,130000,vbiny,20000,80000,vbinz,-37000,3000);

  int width=hvertex->GetNbinsX();
  int height=hvertex->GetNbinsY();
  int depth=hvertex->GetNbinsZ();

  int xmin = hvertex->GetXaxis()->GetXmin();
  int xmax = hvertex->GetXaxis()->GetXmax();
  int ymin = hvertex->GetYaxis()->GetXmin();
  int ymax = hvertex->GetYaxis()->GetXmax();
  int zmin = hvertex->GetZaxis()->GetXmin();
  int zmax = hvertex->GetZaxis()->GetXmax();
  float xrange = (xmax-xmin)/width;
  float yrange = (ymax-ymin)/height;
  float zrange = (zmax-zmin)/depth; 
  //int counts=3;
  
  cout << width << " " << height << " " << depth << endl;
  cout << xmin << " " << xmax << " " << ymin << " " << ymax << " " << zmin << " " << zmax << endl;
  cout << xrange << " " << yrange << " " << zrange << endl;
  
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > vtx_list;
  //std::vector<std::vector<std::vector<std::vector<unsigned int> > > > vtx_ntrk;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > trk_list;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > trk_entry;   
  vtx_list.resize(width);
  //vtx_ntrk.resize(width);
  trk_list.resize(width);
  trk_entry.resize(width);
  for(int i=0;i<width;i++){
    //y axis size
    vtx_list[i].resize(height);
    //vtx_ntrk[i].resize(height);
    trk_list[i].resize(height);
    trk_entry[i].resize(height);
    for(int j=0;j<height;j++){
      //z axis size
      vtx_list[i][j].resize(depth);
      //vtx_ntrk[i][j].resize(depth);
      trk_list[i][j].resize(depth);
      trk_entry[i][j].resize(depth);
    }
  }
  

  
  //  std::vector <std::vector<double> > cl_el;
  //  ev_x_pos2.clear();
  TFile * f_out = new TFile("vtx_MC_analysis_27_05_19.root","RECREATE");
  //fedratree->CloneTree()->Write();
  TTree * Tree_out = new TTree();
  //Tree_out =(TTree*)fedratree->CopyTree("");
  vtx_fedratree->CloneTree()->Write();
  Tree_out->SetName("ana");
  Tree_out->Branch("vtx_good",&vtx_good,"vtx_good/B");
  //Tree_out->Branch("evt_",&primary_proton,"primary_proton/B");
  //Tree_out->Branch("daughter_proton",&daughter_proton,"daughter_proton/B");
  Tree_out->Branch("mp_eventID",&mp_eventID,"mp_eventID/I");
  //Tree_out->Branch("mp_ev_pdgID",&mp_ev_pdgID,"mp_ev_pdgID/I");
  Tree_out->Branch("mp_motherID",&mp_motherID,"mp_motherID/I");
  Tree_out->Branch("mp_procID",&mp_procID,"mp_procID/I");
  Tree_out->Branch("mp_pdgID",&mp_pdgID,"mp_pdgID/I");
  Tree_out->Branch("mp_vx",&mp_vx,"mp_vx/F");
  Tree_out->Branch("mp_vy",&mp_vy,"mp_vy/F");
  Tree_out->Branch("mp_vz",&mp_vz,"mp_vz/F");
  Tree_out->Branch("mp_vdx",&mp_vdx,"mp_vdx/F");
  Tree_out->Branch("mp_vdy",&mp_vdy,"mp_vdy/F");
  Tree_out->Branch("mp_vdz",&mp_vdz,"mp_vdz/F");
  Tree_out->Branch("max_freq",&max_freq,"max_freq/I");
  Tree_out->Branch("mean_freq",&mean_freq,"mean_freq/F");
  Tree_out->Branch("max_ev_freq",&max_ev_freq,"max_ev_freq/I");
  Tree_out->Branch("mean_ev_freq",&mean_ev_freq,"mean_ev_freq/F");
  Tree_out->Branch("n_electrons",&n_electrons,"n_electrons/I");
  Tree_out->AddFriend(vtx_fedratree);
     
  
  std::vector<int> vtx_good_vec(nentries_vtx_fedra);
  std::vector<int> vtx_mp_eventID(nentries_vtx_fedra);
  std::vector<int> vtx_mp_motherID(nentries_vtx_fedra);
  std::vector<int> vtx_mp_pdgID(nentries_vtx_fedra);
  std::vector<int> vtx_mp_flag(nentries_vtx_fedra);
  std::vector<int> vtx_ch_id1(nentries_vtx_fedra);
  std::vector<int> vtx_ch_id2(nentries_vtx_fedra);
  std::vector<int> vtx_mc_id_primary_proton(MCevents);
  std::fill(vtx_mc_id_primary_proton.begin(), vtx_mc_id_primary_proton.end(), -1);
  


  // VERTEX STUDY

  //int first_entry=0;
  //int last_entry=nentries_vtx_fedra;
   
  for (Long64_t ientry=0;ientry<nentries_vtx_fedra;ientry++) { // Vertex index

    if(ientry%100==0)cout << "VID " << ientry << endl;
    log_vtx << "\nFedra Vertex number "<< ientry << endl;
    log_vtx << "MC_EventID / MCtrkID / MCMotherID / MCMotherpdgCode / MCpdgCode / MCProcessId" << endl;
    vtx_fedratree->GetEntry(ientry);
    vtx_good_vec[ientry]=-1;
    vtx_mp_eventID[ientry]=-1;
    vtx_mp_motherID[ientry]=-1;
    vtx_mp_pdgID[ientry]=-1;
    vtx_mp_flag[ientry]=-1;
    vtx_ch_id1[ientry]=-1;
    vtx_ch_id2[ientry]=-1;
		
    if(vx>xmin && vx<xmax && vy>ymin && vy<ymax && vz>zmin && vz<zmax ){

      vtx_trackmother.clear();
      vtx_trackprocess.clear();
      vtx_trackpdgcode.clear();
      vtx_motherpdgcode.clear();
      vtx_trackid.clear();
      vtx_trackstartX.clear();
      vtx_trackstartY.clear();
      vtx_trackstartZ.clear();
      vtx_trackeventId.clear();
    


      int ix = floor((vx-xmin)/xrange);
      int iy = floor((vy-ymin)/yrange);
      int iz = floor((vz-zmin)/zrange);
      //cout << vx << " " << vy <<" " << vz<< endl;
      //cout << ix << " " << iy << " " << iz << endl;
     
      vtx_list.at(ix).at(iy).at(iz).push_back(vID);
      //vtx_ntrk.at(ix).at(iy).at(iz).push_back(n);

      //cout <<  "size " << vtx_list.at(ix).at(iy).at(iz).size() << endl;
     
      for (int itrk=0; itrk<n;itrk++) {  // Track index
       
	//vtx_trackeventId.push_back(MCEventID[itrk]);
	simtree->GetEntry(MCEventID[itrk]);
	//cout << "tracks " << ientry << " " << n << " " << MCTrackID[itrk] <<  endl;
       
       
	//std::vector<int> vtx_motherid_event(MCTrack_);
	//cout << MCTrack_ << endl;
       
	for (Long64_t jn=0; jn<MCTrack_;jn++) {
	 
	  if(jn==MCTrackID[itrk]){
	    vtx_trackeventId.push_back(MCEventID[itrk]);
	    vtx_trackpdgcode.push_back(MCTrack_fPdgCode[jn]);
	    vtx_trackmother.push_back(MCTrack_fMotherId[jn]);
	    vtx_trackprocess.push_back(MCTrack_fProcID[jn]);
	    vtx_trackid.push_back(jn);
	    //cout << ientry << " " <<  jn << " " << MCTrack_fPdgCode[jn] << endl;
	    vtx_trackstartX.push_back(MCTrack_fStartX[jn]);
	    vtx_trackstartY.push_back(MCTrack_fStartY[jn]);
	    vtx_trackstartZ.push_back(MCTrack_fStartZ[jn]);
	  }
	}
      }
     
     
      for(int h=0; h<vtx_trackmother.size(); h++){
	int tmp_jn=-1;
	if(vtx_trackmother.at(h) == -1) vtx_motherpdgcode.push_back(2212);
	else{
	  bool found_mother=false;
	  simtree->GetEntry(vtx_trackeventId.at(h));
	  for (Long64_t jn=0; jn<MCTrack_;jn++) {
	    if(vtx_trackmother.at(h) == jn && !found_mother){
	      found_mother=true;
	      tmp_jn = jn;
	      vtx_motherpdgcode.push_back(MCTrack_fPdgCode[jn]);
	    }
	  }
	}
	//cout <<"LOG "<< ientry << " " << vtx_trackeventId.at(h)  << " " << h << " " << vtx_trackid.at(h) << " " << vtx_trackmother.at(h) << " " << vtx_trackpdgcode.at(h) << " " << vtx_motherpdgcode.at(h) << " " << tmp_jn << endl;
	if(abs(vtx_motherpdgcode.at(h))==411)log_vtx << "Figlio del D+ " << endl;
	if(abs(vtx_motherpdgcode.at(h))==421)log_vtx << "Figlio del D0 " << endl;
	if(abs(vtx_motherpdgcode.at(h))==431)log_vtx << "Figlio del Ds " << endl;
	if(abs(vtx_motherpdgcode.at(h))==441)log_vtx << "Figlio dell'Eta_c " << endl;
	if(abs(vtx_motherpdgcode.at(h))==4122)log_vtx << "Figlio del Lambda_c+ " << endl;
	if(abs(vtx_motherpdgcode.at(h))==4132)log_vtx << "Figlio della Xsi_c0 " << endl;
	if(abs(vtx_motherpdgcode.at(h))==4232)log_vtx << "Figlio della Xsi_c+ " << endl;
	if(abs(vtx_motherpdgcode.at(h))==4332)log_vtx << "Figlio dell'Omega_c0 " << endl;
	if(abs(vtx_trackpdgcode.at(h))==411)log_vtx << "Charm D+ " << endl;
	if(abs(vtx_trackpdgcode.at(h))==421)log_vtx << "Charm D0 " << endl;
	if(abs(vtx_trackpdgcode.at(h))==431)log_vtx << "Charm Ds " << endl;
	if(abs(vtx_trackpdgcode.at(h))==441)log_vtx << "Charm Eta_c " << endl;
	if(abs(vtx_trackpdgcode.at(h))==4122)log_vtx << "Charm Lambda_c+ " << endl;
	if(abs(vtx_trackpdgcode.at(h))==4132)log_vtx << "Charm Xsi_c0 " << endl;
	if(abs(vtx_trackpdgcode.at(h))==4232)log_vtx << "Charm Xsi_c+ " << endl;
	if(abs(vtx_trackpdgcode.at(h))==4332)log_vtx << "Charm Omega_c0 " << endl;
	if(vtx_trackmother.at(h)==1)vtx_ch_id1[ientry]=1;
	if(vtx_trackmother.at(h)!=1 && (abs(vtx_motherpdgcode.at(h))==411 || abs(vtx_motherpdgcode.at(h))==421 || abs(vtx_motherpdgcode.at(h))==431 || abs(vtx_motherpdgcode.at(h))==4122))vtx_ch_id2[ientry] = vtx_trackmother.at(h);   
	log_vtx << vtx_trackeventId.at(h) <<  " "  << vtx_trackid.at(h) << " " << vtx_trackmother.at(h) << " " << vtx_trackpdgcode.at(h) << " " << vtx_motherpdgcode.at(h) << endl;
       
      }
     
     
      std::vector<int> vtx_ntracks_same_mother(n);
      //vtx_ntracks_same_mother.clear();
      std::vector<int> vtx_ntracks_same_event(n);
      //vtx_ntracks_same_event.clear();
     
      /// CONTA LE OCCORRENZE PER MOTHER
      for(int h=0; h<vtx_trackmother.size(); h++){
	flag=0;
	for (int i=0;i<h;i++){
	  if (vtx_trackmother[h]==vtx_trackmother[i]){// controlla che quel numero non sia già stato considerato in precedenza
	    flag=1;
	    break;
	  }
	  //cout << "before "<< ientry << " " << h << " " << i << " " << flag << " " << vtx_trackmother[h] << " " << vtx_ntracks_same_mother[h] << endl;
	}
	if (flag==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	  vtx_ntracks_same_mother[h]=1; // assegna 1 alla mother della prima traccia che incontra
	  for (int i=h+1;i<vtx_trackmother.size();i++){		 
	    if(vtx_trackmother[h] == vtx_trackmother[i] && vtx_trackid[h]!=vtx_trackid[i] && vtx_trackeventId[h]==vtx_trackeventId[i]) vtx_ntracks_same_mother[h]++;
	    //  cout << "after "<< ientry << " " <<  h << " " << i << " " << flag << " " <<" " << vtx_trackmother[i] << " " <<  vtx_trackmother[h] << " " << vtx_ntracks_same_mother[h] << endl;
	  }
	  //cout << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackmother[h] << "trk Frequency " << vtx_ntracks_same_mother[h] << " times" << endl;
	  //log << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackmother[h] << "trk Frequency " << vtx_ntracks_same_mother[h] << " times" << endl;
	}
	//cout << "ntracks "<< ientry << " " <<  h << " " <<  vtx_trackmother[h] << endl;
      }

     
      /// CONTA LE OCCORRENZE PER EVENTO
      for(int h=0; h<vtx_trackmother.size(); h++){
	flag_ev=0;
	for (int i=0;i<h;i++){
	  if (vtx_trackeventId[h]==vtx_trackeventId[i]){// controlla che quel numero non sia già stato considerato in precedenza
	    flag_ev=1;
	    break;
	  }
	}
	if (flag_ev==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	  vtx_ntracks_same_event[h]=1;
	  for (int i=h+1;i<vtx_trackmother.size();i++){		 
	    if(vtx_trackeventId[h] == vtx_trackeventId[i] && vtx_trackid[h]!=vtx_trackid[i]) vtx_ntracks_same_event[h]++;
	  }
	  //cout << ientry << " " << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackid[h] << "trk Frequency " << vtx_ntracks_same_event[h] << " times" << endl;
	}
      }
     
     
      vector<int>::iterator max_mID = max_element(vtx_ntracks_same_mother.begin(), vtx_ntracks_same_mother.end()); // c++11
      auto max_occurrence = *max_element(vtx_ntracks_same_mother.begin(), vtx_ntracks_same_mother.end()); // c++11    
      int iel = distance(vtx_ntracks_same_mother.begin(), max_mID);
     
      vector<int>::iterator max_evID = max_element(vtx_ntracks_same_event.begin(), vtx_ntracks_same_event.end()); // c++11
      auto max_ev_occurrence = *max_element(vtx_ntracks_same_event.begin(), vtx_ntracks_same_event.end()); // c++11
      int iel_ev = distance(vtx_ntracks_same_event.begin(), max_evID);
     
      vector <float> vtx_x;
      vector <float> vtx_y;
      vector <float> vtx_z;
      vtx_x.clear();
      vtx_y.clear();
      vtx_z.clear();

      int max_motherId = *max_mID;
      //auto max_eventId = *max_evID;

      //cout << "max occ "<<ientry << " " << iel << " " << max_occurrence << endl;
      //cout << vtx_trackeventId.at(iel_ev) << endl;
      n_electrons=0;
      //daughter_proton=false;
     
      if(max_ev_occurrence==1){
	mp_eventID=-1;
	//primary_proton=false;
      }
      if(max_ev_occurrence>1) {
	//cout << "hello " << max_ev_occurrence <<" " << max_occurrence  << endl;
	mp_eventID = vtx_trackeventId.at(iel_ev);
	//cout << ientry << " " << mp_eventID << endl;
	//primary_proton=true;
	for(int i=0; i<vtx_trackpdgcode.size(); i++){
	  //if(abs(vtx_trackpdgcode.at(i))==2212 && vtx_motherpdgcode.at(i)==-1)daughter_proton=true;
	  if(abs(vtx_trackpdgcode.at(i))==11 && vtx_trackeventId.at(i)==mp_eventID)n_electrons++;
	}
      }
     
      if(max_occurrence==1){
	mp_motherID=-2;
	mp_procID=-2;
	mp_pdgID=-2;
	mp_vx=0;
	mp_vy=0;
	mp_vz=0;
	mp_vdx=-1;
	mp_vdy=-1;
	mp_vdz=-1;
      }
      if(max_occurrence>1) {
              
	int trk_index=0;

       
	for(unsigned int itrk=0;itrk<vtx_trackmother.size();itrk++){
	  //cout <<"crash "<< ientry << " " << vtx_trackmother.size() << " " <<  max_occurrence << " " << itrk << " " << iel  << endl;
	 
	  if(vtx_trackmother.at(itrk)==vtx_trackmother.at(iel)){
	    vtx_x.push_back(vtx_trackstartX.at(itrk));
	    vtx_y.push_back(vtx_trackstartY.at(itrk));
	    vtx_z.push_back(vtx_trackstartZ.at(itrk));	   
	    //cout << ientry << " " << trk_index << " " << vtx_x.at(trk_index) << " " << vtx_y.at(trk_index) << " " << vtx_z.at(trk_index) << endl;
	    trk_index++;
	  }
	 
	}
       
	mp_vx = accumulate(vtx_x.begin(), vtx_x.end(), 0.0)/vtx_x.size();
	mp_vy = accumulate(vtx_y.begin(), vtx_y.end(), 0.0)/vtx_y.size();
	mp_vz = accumulate(vtx_z.begin(), vtx_z.end(), 0.0)/vtx_z.size();
      
	mp_vx = mp_vx*pow(10,4) + 62500;
	mp_vy = mp_vy*pow(10,4) + 49500;
	mp_vz = (mp_vz-125.6)*pow(10,4);       

	auto max_vtx_x = *max_element(vtx_x.begin(), vtx_x.end()); 
	auto max_vtx_y = *max_element(vtx_y.begin(), vtx_y.end()); 
	auto max_vtx_z = *max_element(vtx_z.begin(), vtx_z.end());
	auto min_vtx_x = *min_element(vtx_x.begin(), vtx_x.end()); 
	auto min_vtx_y = *min_element(vtx_y.begin(), vtx_y.end()); 
	auto min_vtx_z = *min_element(vtx_z.begin(), vtx_z.end());        

	mp_vdx = max_vtx_x - min_vtx_x;
	mp_vdy = max_vtx_y - min_vtx_y;
	mp_vdz = max_vtx_z - min_vtx_z;
                    
	mp_motherID = vtx_trackmother.at(iel);
	mp_procID = vtx_trackprocess.at(iel);
	mp_pdgID = vtx_motherpdgcode.at(iel);

	//cout << ientry << " " << mp_vx << " " << vx << " " << mp_vy << " " << vy << " " << mp_vz << " " << vz << endl;
	//cout << ientry <<" " <<  mp_vdx << " " << mp_vdy << " " << mp_vdz << endl;
        
	//cout << ientry << " " << vtx_trackmother.at(iel) << " " << endl;
       
      }
      mean_freq=(1.0*max_occurrence)/n;
      max_freq=max_occurrence;

      mean_ev_freq=(1.0*max_ev_occurrence)/n;
      max_ev_freq=max_ev_occurrence;

      log_vtx << "The max number of tracks with same mother ID "<< mp_motherID << " in the Vertex " << ientry << " is " << max_occurrence << endl;
      log_vtx << "The max number of tracks with same event ID "<< mp_eventID << " in the Vertex " << ientry << " is " << max_ev_occurrence << endl;


      //------ DESCRIPTION  -----------//
     
      //if(primary_proton){
      if(mp_pdgID!=-2){
	if(mp_motherID==-1){
	  log_vtx << "It's a primary proton vertex with pdg code " << mp_pdgID << endl;
	  vtx_good=true;
	  vtx_good_vec[ientry]=0;
	  vtx_mp_eventID[ientry]=mp_eventID;
	  vtx_mp_motherID[ientry]=mp_motherID;
	  vtx_mp_pdgID[ientry]=mp_pdgID;
	  vtx_mp_flag[ientry]=1;
	  vtx_mc_id_primary_proton[mp_eventID]=ientry;
	}
	if(mp_motherID==0 || (abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122)){
	  log_vtx << "It's a charm vertex with pdg code " << mp_pdgID << endl;
	  vtx_good=true;
	  vtx_good_vec[ientry]=1;
	  vtx_mp_eventID[ientry]=mp_eventID;
	  vtx_mp_motherID[ientry]=mp_motherID;
	  vtx_mp_pdgID[ientry]=mp_pdgID;
	  vtx_mp_flag[ientry]=0;
	}
	if( mp_motherID>0 && !(abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122)){
	  log_vtx << "It's a  secondary vertex with pdg code " << mp_pdgID << endl;
	  vtx_good=true;
	  vtx_good_vec[ientry]=2;
	  vtx_mp_eventID[ientry]=mp_eventID;
	  vtx_mp_motherID[ientry]=mp_motherID;
	  vtx_mp_pdgID[ientry]=mp_pdgID;
	  vtx_mp_flag[ientry]=2;
	}	 
      }
      else {
	log_vtx << "It's a fake\n";
	vtx_mp_flag[ientry]=-1;
	vtx_good=false;
      }

      //------------------------- //
     
      //if(daughter_proton)log << "Warning! A primary proton is a daughter!\n";
      //if(mp_motherID==0 && mp_pdgID==2212)log << "It's a primary proton\n";
      //if(mp_motherID>0 && mp_pdgID==2212)log << "It's a secondary proton\n";
      //if(mp_motherID==-1 && mp_pdgID==2212)log << "Warning! A primary proton is a daughter!\n";

      log_vtx << "VID / MostProbable Vertex X (Y,Z) from MC / Vertex X (Y,Z) from Fedra\n";
      log_vtx << ientry << " " << mp_vx << " " << vx << " " << mp_vy << " " << vy << " " << mp_vz << " " << vz << endl;
      log_vtx << "Gap between MC and Fedra Vertex position (X,Y,Z)\n";
      log_vtx << ientry <<" " <<  mp_vdx << " " << mp_vdy << " " << mp_vdz << endl;
     
      hoccurrence->Fill(mean_freq);
      Tree_out->Fill();
       
      //cout << vx << " " << vy << " " << vz << endl;
      hvertex->Fill(vx,vy,vz);
     
    }
    //if(ientry==99)ientry=18800;
  }
   
  log_vtx << "VERTEXING LOG TERMINED" << endl;
  TCanvas *c1 = new TCanvas();
  hoccurrence->Draw();
  TCanvas *c2 = new TCanvas();
  hvertex->Draw("");
  log_vtx.close();
  Tree_out->Write();
  f_out->Close();
  // TRACKS


  // TRACKS NOT IN A VERTEX
   
  for (Long64_t ientry=0; ientry<nentries_trk_fedra;ientry++) { // Tracks index

    if(ientry%10000==0)cout << "TID " << ientry << endl;
    //log << "\nFedra Track number "<< ientry << endl;
    //log << "MC_EventID / MCtrkID / MCMotherID / MCMotherpdgCode / MCpdgCode / MCProcessId" << endl;
    trk_fedratree->GetEntry(ientry);

    if(t_->X()>xmin && t_->X()<xmax && t_->Y()>ymin && t_->Y()<ymax && t_->Z()>zmin && t_->Z()<zmax ){
     
      int ix = floor((t_->X()-xmin)/xrange);
      int iy = floor((t_->Y()-ymin)/yrange);
      int iz = floor((t_->Z()-zmin)/zrange);
      //cout << ientry << " " << t_->X() << " " << t_->Y() << " " << t_->Z() <<  endl;
      //cout << ix << " " << iy << " " << iz << endl;
     
      trk_list.at(ix).at(iy).at(iz).push_back(trid);
      trk_entry.at(ix).at(iy).at(iz).push_back(ientry);

      //cout << ientry << " " << t_->MCEvt() << " " << t_->Aid(0) << endl;
      //cout << ientry << " " << t_->X() << " " << t_->Y() << " " << t_->Z() <<  endl;
    }
  }

   
  TFile * f_ds = new TFile("decay_search_result.root","RECREATE");
  TTree *TreeDS = new TTree("ds","decay search");
  TreeDS->Branch("vtx",&vtx,"event/I:flag/I:fe_id/I:mc_id/I:mc_pdg/I:fe_pproton_id/I:ntrk/I:mc_ch_id1/I:mc_ch_id2/I:mc_nd[2]/I:mc_dz[2]/D:mc_dl[2]/D");
  //TreeDS->Branch("trk",&trk,"trk.fe_itrk:trk.icharm");

  TreeDS->Branch("trk.fe_itrk",&trk.fe_itrk);
  TreeDS->Branch("trk.icharm",&trk.icharm);
  TreeDS->Branch("trk.fe_id",&trk.fe_id);
  TreeDS->Branch("trk.mc_id",&trk.mc_id);
  TreeDS->Branch("trk.pdg",&trk.pdg);
  TreeDS->Branch("trk.nseg",&trk.nseg);
  TreeDS->Branch("trk.type",&trk.type);
  TreeDS->Branch("trk.split",&trk.split);
  TreeDS->Branch("trk.fe_dz",&trk.fe_dz);
  TreeDS->Branch("trk.fe_dl",&trk.fe_dl);
  TreeDS->Branch("trk.ip",&trk.ip);
  TreeDS->Branch("trk.ka",&trk.ka);
  TreeDS->Branch("trk.iseg_rmax",&trk.iseg_rmax);
  TreeDS->Branch("trk.rmax",&trk.rmax);
  TreeDS->Branch("trk.iseg_ipmax",&trk.iseg_ipmax);
  TreeDS->Branch("trk.ipmax",&trk.ipmax);
  TreeDS->Branch("trk.pardaug",&trk.pardaug);
  
   

  // DECAY SEARCH

  int first_entry=0;
  int last_entry=nentries_vtx_fedra;
  int event = 0;

   
  //EdbSegP *seg=0;
   
  for (Long64_t ientry=first_entry; ientry<last_entry;ientry++) { // Vertex index
     
    if(ientry%1==0)cout << "VID DS " << ientry << endl;

     
    int mp_event =  vtx_mp_eventID[ientry];
     
    //cout <<"primary_proton " << ientry << " " <<  mp_event << " " << primary_proton_vtx <<endl;
    int mp_mother =  vtx_mp_motherID[ientry];
    int mp_pdgID = vtx_mp_pdgID[ientry];
    int mp_flag = vtx_mp_flag[ientry];
    int ch_id1 = vtx_ch_id1[ientry];
    int ch_id2 = vtx_ch_id2[ientry];
    //cout << ch_id1 << " " << ch_id2 << endl;
    bool is_a_charm = false;
     
    if(vtx_good_vec.at(ientry)!=-1 && (abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122))is_a_charm=true;

    if(vtx_good_vec.at(ientry)!=-1 && (mp_mother==-1 || (abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122))){

      Set0();
      int nfound=0;
       
      int tmp_i=-1;
      int tmp_j=-1;
      int tmp_k=-1;
      int tmp_l=-1;
      int tmp_i_p1=-1;
      int tmp_j_p1=-1;
      int tmp_k_p1=-1;
      int tmp_i_m1=-1;
      int tmp_j_m1=-1;
      int tmp_k_m1=-1;
      int tmp_vID=-1;
      double tmp_vx=-1;
      double tmp_vy=-1;
      double tmp_vz=-1;

      double mc_tx=0;
      double mc_ty=0;

      TVector3 vtx_pos;
      TVector3 vtx_primary_proton_pos;
      int primary_proton_vtx = vtx_mc_id_primary_proton[mp_event];
      if(vtx_mc_id_primary_proton[mp_event]!=-1){
	vtx_fedratree->GetEntry(vtx_mc_id_primary_proton[mp_event]);
	vtx_primary_proton_pos.SetXYZ(vx,vy,vz);
      }
      else vtx_primary_proton_pos.SetXYZ(-1,-1,-1);
      TVector3 trk_pos;

      vector<int> trk_mc_id;
      vector<int> trk_z;
      trk_mc_id.clear();
      trk_z.clear();

      bool splitted_track=false;
      bool pardaug=false;
       
      double fdz = 0;
      double fdl = 0;
      double ip = 0;
      double ka = 0;
      float rmax=0;
      int ty = 0;
      int mc_nd1=0;
      int mc_nd2=0;
      int mc_pdg=0;
      int fnd=0;
      int icharm=-1;

      float dtheta[2]={-1,-1};
      float ipseg[2]={-1,-1};
       
      vtx_fedratree->GetEntry(ientry);

      // RICERCA DEL CUBO IESIMO DEL VERTICE IESIMO
       
      for(int i=0;i<width;i++){
	for(int j=0;j<height;j++){
	  for(int k=0;k<depth;k++){
	    //cout << "size " << vtx_list.at(i).at(j).at(k).size() << endl;
	    for(int l=0;l<vtx_list.at(i).at(j).at(k).size();l++){
	      if(vtx_list.at(i).at(j).at(k).size()>0){
		//log2 << i << " " << j << " " << k << " " << l <<  " " << vtx_list.at(i).at(j).at(k).at(l) << " " << "0" <<  endl;
		if(vtx_list.at(i).at(j).at(k).at(l)==vID){
		  log_decay << "\nVtx " << vtx_list.at(i).at(j).at(k).at(l) << " found in coord " << i << " " << j << " " << k << " " << l << " of MCEvent " << mp_event << " with Fedra xyz " << vx << " " << vy << " " << vz <<  endl;
		  log_decay << "MC DZ CHARM1 FROM PRIMARY VTX  " << charm_1_z[mp_event] << " - MC DZ CHARM2 FROM PRIMARY VTX " << charm_2_z[mp_event] << endl;
		  if(is_a_charm==true)log_decay << "CHARM VERTEX" <<  endl;

		  int ndaug1 = charm_1_id.at(mp_event).size();
		  int ndaug2 = charm_2_id.at(mp_event).size();
		  //mnd = ndaug1 + ndaug2;
		  
		  FillVertex(event,mp_flag,vID,mp_event,mp_pdgID,primary_proton_vtx,n,ch_id1,ch_id2,ndaug1,ndaug2,charm_1_z[mp_event],charm_2_z[mp_event],charm_1_3d[mp_event],charm_2_3d[mp_event]);

		  // inizializzo alcune variabili
		  vtx_pos.SetXYZ(vx,vy,vz);
		  tmp_vID=vID;
		  tmp_vx=vx;
		  tmp_vy=vy;
		  tmp_vz=vz;
		  tmp_i=i;
		  tmp_j=j;
		  tmp_k=k;
		  tmp_l=l;

		  //cout << vtx_pos(0) << " " << vtx_primary_proton_pos(0) << endl;

		  vertexobject = (EdbVertex *)(vertexrec->eVTX->At(vID));

		  // loop sulle tracce del vertice
		  for(int itrk=0;itrk<n;itrk++){
		    icharm=-1;
		    rmax=0;
		    pardaug=false;
		     
		    simtree->GetEntry(MCEventID[itrk]);

		    // vedo se la traccia è un charm
		    for (Long64_t jn=0; jn<MCTrack_;jn++) {			
		      if(jn==MCTrackID[itrk]){
			mc_pdg=MCTrack_fPdgCode[jn];
			if(jn==1)icharm=1;
			if(jn!=1 && (abs(mc_pdg)==411 || abs(mc_pdg)==421 || abs(mc_pdg)==431 || abs(mc_pdg)==4122))icharm=2;
		      }
		    } // ---------//
		     
		    ftrack = vertexobject->GetTrack(itrk);
		    //cout << ientry << " " << itrk << " " << ftrack->N() << endl;

		    // SOLO CHARM
		    if(icharm==1 || icharm==2){
		       
		      trk_mc_id.push_back(MCTrackID[itrk]);
		      trk_z.push_back(ftrack->Zstart());
		      //cout <<"size "<< (trk_mc_id.size()-1) << endl;
		      for(int s=0;s<(trk_mc_id.size()-1);s++){
			//cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
			if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){ // 2 tracce in Fedra hanno lo stesso ID monte carlo
			  //cout << "hello"<<endl;
			  splitted_track=true;
			   
			  if(ftrack->Zstart()<trk_z.at(s))break;
			  else {
			    trk.fe_dz[s]=ftrack->Zstart()-tmp_vz;
			    //cout << "eccomi3"<< endl;
			    break;
			  }
			}
		      }
		       
		      //if(!splitted_track){
		      if(icharm==1)log_decay << "Linked to the vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Charm1 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		      if(icharm==2)log_decay << "Linked to the vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Charm2 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		      log_decay << "MCEvt " <<  MCEventID[itrk] << " - MCTrkId " << MCTrackID[itrk] << " MotherID " << MCMotherID[itrk] << endl;
		      if(vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==411 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==421 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==431 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==4122)log_decay << "Charm vertex with pdg " << vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]  << endl;;
		      //(abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122)
		      //trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
		      trk_pos.SetXYZ(ftrack->GetSegmentFirst()->X(),ftrack->GetSegmentFirst()->Y(),ftrack->Zstart());
		       
		       
		      // tipo di traccia
		      if(vID==tmp_vID)ty=1;
		      else{
			if(vtx_good_vec.at(ientry)==-1)ty=4;
			else{
			  if(MCEventID[itrk]==mp_event)ty=2;
			  if(MCEventID[itrk]!=mp_event)ty=3;
			}
		      } //---------//
		       
		       
		      fdz = trk_pos(2) - tmp_vz;
		      fdl = sqrt(pow(trk_pos(0) - tmp_vx,2)+pow(trk_pos(1) - tmp_vy,2)+pow(trk_pos(2) - tmp_vz,2));
		      ip = IPtoVertex(vtx_pos, trk_pos, TX[itrk], TY[itrk]);
		      // cout << "ip_vtx " << ip << endl;
		      //cout << "ip_vtx " << ip << " " << vtx_pos(0) <<" " << vtx_pos(1) << vtx_pos(2) << " " << trk_pos(0) << " " << trk_pos(1)<< " " << trk_pos(2) << " " << TX[itrk] << " " << TY[itrk] <<   endl;
		      ka = sqrt(pow(mc_tx - TX[itrk],2) + pow(mc_ty - TY[itrk],2));
		      FedraTrackKink(ftrack,dtheta);
		      SegIPtoVertex(vtx_primary_proton_pos, ftrack, ipseg);
		      pardaug = ParentDaugther(ftrack);
		      FillTrack(nfound,icharm,0,MCTrackID[itrk],mc_pdg,nseg[itrk],ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
		      nfound++; 
		    } // fine solo charm

		    else{
		      int ndaug1 = charm_1_id.at(mp_event).size();
		      mc_tx = (c1vxMC[mp_event]-pvxMC[mp_event])/charm_1_z[mp_event];
		      mc_ty = (c1vyMC[mp_event]-pvyMC[mp_event])/charm_1_z[mp_event];
		     

		      splitted_track=false;
		      for (int idaug = 0; idaug < ndaug1; idaug++) {
		     
			if(MCTrackID[itrk]==daug_1_id[mp_event][idaug] && MCEventID[itrk]==fEventId[mp_event][idaug] && MCMotherID[itrk]==charm_1_id[mp_event][idaug]){
			  icharm=11;
			  trk_mc_id.push_back(MCTrackID[itrk]);
			  trk_z.push_back(ftrack->Zstart());
			  //cout <<"size "<< (trk_mc_id.size()-1) << endl;
			  for(int s=0;s<(trk_mc_id.size()-1);s++){
			    // cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
			    if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){
			      log_decay << "Splitted track"<<endl;
			      splitted_track=true;
			      if(ftrack->Zstart()<trk_z.at(s))break;
			      else {
				trk.fe_dz[s]=ftrack->Zstart()-tmp_vz;
				//cout << "eccomi1" << endl;
				break;
			      }
			    }
			  }
			  //cout <<"2 "<< ftrack->GetSegmentFirst()->Y() << endl;
			  //if(!splitted_track){
			  log_decay << "Linked to the vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Daugther charm1 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
			  log_decay << "MCEvt " <<  MCEventID[itrk] << " - MCTrkId " << MCTrackID[itrk] << " MotherID " << MCMotherID[itrk] << endl;
			  if(vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==411 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==421 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==431 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==4122)log_decay << "Charm vertex with pdg " << vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]  << endl;;
			  //(abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122)
			  //trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
			   
			  trk_pos.SetXYZ(ftrack->GetSegmentFirst()->X(),ftrack->GetSegmentFirst()->Y(),ftrack->Zstart());
			   
			   
			  if(vID==tmp_vID)ty=1;
			  else{
			    if(vtx_good_vec.at(ientry)==-1)ty=4;
			    else{
			      if(MCEventID[itrk]==mp_event)ty=2;
			      if(MCEventID[itrk]!=mp_event)ty=3;
			    }
			  }
			   
			   
			  fdz = trk_pos(2) - tmp_vz;
			  fdl = sqrt(pow(trk_pos(0) - tmp_vx,2)+pow(trk_pos(1) - tmp_vy,2)+pow(trk_pos(2) - tmp_vz,2));
			  ip = IPtoVertex(vtx_pos, trk_pos, TX[itrk], TY[itrk]);
			  // cout << "ip_vtx " << ip << endl;
			  //cout << "ip_vtx " << ip << " " << vtx_pos(0) <<" " << vtx_pos(1) << vtx_pos(2) << " " << trk_pos(0) << " " << trk_pos(1)<< " " << trk_pos(2) << " " << TX[itrk] << " " << TY[itrk] <<   endl;
			  ka = sqrt(pow(mc_tx - TX[itrk],2) + pow(mc_ty - TY[itrk],2));
			  //rmax = FedraTrackKink(ftrack);
			  FedraTrackKink(ftrack,dtheta);
			  SegIPtoVertex(vtx_primary_proton_pos, ftrack, ipseg);
			  pardaug = ParentDaugther(ftrack);
			  FillTrack(nfound,icharm,0,MCTrackID[itrk],mc_pdg,nseg[itrk],ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
			  nfound++;
			  //cout << "SOLO CHARM" <<endl;
			}
			//cout << "DAUG1" <<endl;
		      }
		   
		      int ndaug2 = charm_2_id.at(mp_event).size();
		      mc_tx = (c2vxMC[mp_event]-pvxMC[mp_event])/charm_2_z[mp_event];
		      mc_ty = (c2vyMC[mp_event]-pvyMC[mp_event])/charm_2_z[mp_event];
		 
		   
		      splitted_track=false;
		      for (int idaug = 0; idaug < ndaug2; idaug++) {   
			//if(MCEventID[itrk]==765)cout << "linked "<<idaug << " " << MCTrackID[itrk] << " " << daug_1_id[mp_event][idaug] << endl;
			if(MCTrackID[itrk]==daug_2_id[mp_event][idaug] && MCEventID[itrk]==fEventId[mp_event][idaug] && MCMotherID[itrk]==charm_2_id[mp_event][idaug]){

			  icharm=22;
			  trk_mc_id.push_back(MCTrackID[itrk]);
			  trk_z.push_back(ftrack->Zstart());
			  // cout <<"size "<< (trk_mc_id.size()-1) << endl;
			  for(int s=0;s<(trk_mc_id.size()-1);s++){
			    //cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
			    if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){
			      //cout << "hello"<<endl;
			      log_decay << "Splitted track"<<endl;
			      splitted_track=true;
			      if(ftrack->Zstart()<trk_z.at(s))break;
			      else {
				trk.fe_dz[s] = ftrack->Zstart()-tmp_vz;
				//cout << "eccomi2 " << trk.fe_dz[s] << endl;
				break;
			      }
			    }
			  }

			  //if(!splitted_track){ 
			  log_decay << "Linked to the vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Daugther charm2 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
			  log_decay << "MCEvt " << MCEventID[itrk] << " - MCTrkId " << MCTrackID[itrk] << " MotherID " << MCMotherID[itrk] << endl;

			  if(vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==411 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==421 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==431 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==4122)log_decay << "Nearby charm vertex with pdg " << vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]  << endl;;
		       
			  //trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
			  trk_pos.SetXYZ(ftrack->GetSegmentFirst()->X(),ftrack->GetSegmentFirst()->Y(),ftrack->Zstart());
		       
			  if(vID==tmp_vID)ty=1;
			  else{
			    if(vtx_good_vec.at(ientry)==-1)ty=4;
			    else{
			      if(MCEventID[itrk]==mp_event)ty=2;
			      if(MCEventID[itrk]!=mp_event)ty=3;
			    }
			  }
		       
		       
			  fdz = trk_pos(2) - tmp_vz;
			  fdl = sqrt(pow(trk_pos(0) - tmp_vx,2)+pow(trk_pos(1) - tmp_vy,2)+pow(trk_pos(2) - tmp_vz,2));
			  ip = IPtoVertex(vtx_pos, trk_pos, TX[itrk], TY[itrk]);
			  ka = sqrt(pow(mc_tx - TX[itrk],2) + pow(mc_ty - TY[itrk],2));
			  //rmax = FedraTrackKink(ftrack);
			  FedraTrackKink(ftrack,dtheta);
			  SegIPtoVertex(vtx_primary_proton_pos, ftrack, ipseg);
			  pardaug = ParentDaugther(ftrack);
			  FillTrack(nfound,icharm,0,MCTrackID[itrk],0,nseg[itrk],ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
		      
			  nfound++;
			  // }
			}
			//cout << "DAUG2" <<endl;
		      }
		    }
		  }

		   
		   
		  event++;
		  break;
		}
	      }
	    }	     
	  }
	}
      }

       
       
      // RICERCA SUI PRIMI VICINI
       

      if(tmp_i==0){
	tmp_i_m1=0;
	tmp_i_p1=tmp_i+1;
      }
      if(tmp_j==0){
	tmp_j_m1=0;
	tmp_j_p1=tmp_j+1;
      }
      if(tmp_k==0){
	tmp_k_m1=0;
	tmp_k_p1=tmp_k+1;
      }
      if(tmp_i==(vbinx-1)){
	tmp_i_m1=tmp_i-1;
	tmp_i_p1=vbinx;
      }
      if(tmp_j==(vbiny-1)){
	tmp_j_m1=tmp_j-1;
	tmp_j_p1=vbiny;
      }
      if(tmp_k==(vbinz-1)){
	tmp_k_m1=tmp_k-1;
	tmp_k_p1=vbinz;
      }

      if(tmp_i!=0 && tmp_i!=(vbinx-1)){
	tmp_i_m1=tmp_i-1;
	tmp_i_p1=tmp_i+2;
      }

      if(tmp_j!=0 && tmp_j!=(vbiny-1)){
	tmp_j_m1=tmp_j-1;
	tmp_j_p1=tmp_j+2;
      }

      if(tmp_k!=0 && tmp_k!=(vbinz-1)){
	tmp_k_m1=tmp_k-1;
	tmp_k_p1=tmp_k+2;
      }

      //cout << tmp_i << " " << tmp_i_m1 << " " << tmp_i_p1 << endl;
      //cout << tmp_j << " " << tmp_j_m1 << " " << tmp_j_p1 << endl;
      //cout << tmp_k << " " << tmp_k_m1 << " " << tmp_k_p1 << endl;


      for(int i=tmp_i_m1;i<tmp_i_p1;i++){
	for(int j=tmp_j_m1;j<tmp_j_p1;j++){
	  for(int k=tmp_k_m1;k<tmp_k_p1;k++){

	    // VERTICI VICINI
	    for(int l=0;l< vtx_list.at(i).at(j).at(k).size();l++){
	      
	      int tmp_index = vtx_fedratree->GetEntryNumberWithIndex(vtx_list.at(i).at(j).at(k).at(l));		 
		 
	      vtx_fedratree->GetEntry(tmp_index);
	      
	      //cout << "view ID "<< vID << " " << tmp_vID << endl;
	      if(vID!=tmp_vID && mp_event==vtx_mp_motherID[tmp_index]){ // solo se l'altro vertice appartiene allo stesso evento montecarlo
	
		//cout << vID << " " << tmp_vID << endl;
		vertexobject_tmp = (EdbVertex *)(vertexrec->eVTX->At(vID));
	
		for(int itrk=0;itrk<n;itrk++){

		  ftrack = vertexobject_tmp->GetTrack(itrk);
		  
		  icharm=-1;
		  rmax=0;
		  pardaug=false;

		  simtree->GetEntry(MCEventID[itrk]);
		 
		  for (Long64_t jn=0; jn<MCTrack_;jn++) {			
		    if(jn==MCTrackID[itrk]){
		      mc_pdg=MCTrack_fPdgCode[jn];
		      if(jn==1)icharm=1;
		      if(jn!=1 && (abs(mc_pdg)==411 || abs(mc_pdg)==421 || abs(mc_pdg)==431 || abs(mc_pdg)==4122))icharm=2;
		    }
		  }

		  //cout <<"1 "<< itrk << " " << MCTrackID[itrk] << " " << icharm << endl; // icharm =-1
		  //cout <<"2 "<< ftrack->GetSegmentFirst()->Y() << endl;
		  //cout <<"2 "<< ftrack->X() << endl;
		  // SOLO CHARM
		  
		  if(icharm==1 || icharm==2){
		    //cout <<"1 "<< itrk << " " << MCTrackID[itrk] <<endl;
		    //cout <<"2 "<< ftrack->X() << endl;
		    trk_mc_id.push_back(MCTrackID[itrk]);
		    trk_z.push_back(ftrack->Zstart());
		    
		    //cout <<"size "<< (trk_mc_id.size()-1) << endl;
		    for(int s=0;s<(trk_mc_id.size()-1);s++){
		      //cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
		      if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){ // 2 tracce in Fedra hanno lo stesso ID monte carlo
			//cout << "hello"<<endl;
			splitted_track=true;
			if(ftrack->Zstart()<trk_z.at(s))break;
			else {
			  trk.fe_dz[s]=ftrack->Zstart()-tmp_vz;
			  //cout << "eccomi3"<< endl;
			  break;
			}
		      }
		    }
		    
		    if(icharm==1)log_decay << "Linked to another vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Charm1 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		    if(icharm==2)log_decay << "Linked to another vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Charm2 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		    log_decay << "MCEvt " <<  MCEventID[itrk] << " - MCTrkId " << MCTrackID[itrk] << " MotherID " << MCMotherID[itrk] << endl;
		    if(vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==411 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==421 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==431 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==4122)log_decay << "Nearby charm vertex with pdg " << vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]  << endl;;
		    //(abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122)
		    //trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
		    trk_pos.SetXYZ(ftrack->GetSegmentFirst()->X(),ftrack->GetSegmentFirst()->Y(),ftrack->Zstart());
		     
		       
		    if(vID==tmp_vID)ty=1;
		    else{
		      if(vtx_good_vec.at(tmp_index)==-1)ty=4;
		      else{
			if(MCEventID[itrk]==mp_event)ty=2;
			if(MCEventID[itrk]!=mp_event)ty=3;
		      }
		    }
		       
		    
		    fdz = trk_pos(2) - tmp_vz;
		    fdl = sqrt(pow(trk_pos(0) - tmp_vx,2)+pow(trk_pos(1) - tmp_vy,2)+pow(trk_pos(2) - tmp_vz,2));
		    ip = IPtoVertex(vtx_pos, trk_pos, TX[itrk], TY[itrk]);
		    // cout << "ip_vtx " << ip << endl;
		    //cout << "ip_vtx " << ip << " " << vtx_pos(0) <<" " << vtx_pos(1) << vtx_pos(2) << " " << trk_pos(0) << " " << trk_pos(1)<< " " << trk_pos(2) << " " << TX[itrk] << " " << TY[itrk] <<   endl;
		    ka = sqrt(pow(mc_tx - TX[itrk],2) + pow(mc_ty - TY[itrk],2));
		    //rmax = FedraTrackKink(ftrack);
		    FedraTrackKink(ftrack,dtheta);
		    SegIPtoVertex(vtx_primary_proton_pos, ftrack, ipseg);
		    pardaug = ParentDaugther(ftrack);
		    FillTrack(nfound,icharm,0,MCTrackID[itrk],mc_pdg,nseg[itrk],ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
		    nfound++;
		   
		    }
		  // fine solo charm
		  else{
		    int ndaug1 = charm_1_id.at(mp_event).size();
		    mc_tx = (c1vxMC[mp_event]-pvxMC[mp_event])/charm_1_z[mp_event];
		    mc_ty = (c1vyMC[mp_event]-pvyMC[mp_event])/charm_1_z[mp_event];
		    
		    
		    splitted_track=false;
		    for (int idaug = 0; idaug < ndaug1; idaug++) {
		     
		      if(MCTrackID[itrk]==daug_1_id[mp_event][idaug] && MCEventID[itrk]==fEventId[mp_event][idaug] && MCMotherID[itrk]==charm_1_id[mp_event][idaug]){
			
			icharm=11;
			trk_mc_id.push_back(MCTrackID[itrk]);
			trk_z.push_back(ftrack->Zstart());
			//cout <<"size "<< (trk_mc_id.size()-1) << endl;
			for(int s=0;s<(trk_mc_id.size()-1);s++){
			  // cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
			  if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){
			    //cout << "hello"<<endl;
			    log_decay << "Splitted track"<<endl;
			    splitted_track=true;
			    if(ftrack->Zstart()<trk_z.at(s))break;
			    else {
			      trk.fe_dz[s]=ftrack->Zstart()-tmp_vz;
			      //cout << "eccomi1" << endl;
			      break;
			    }
			  }
			}
		
			//if(!splitted_track){
			log_decay << "Linked to another vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Daugther charm1 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
			log_decay << "MCEvt " <<  MCEventID[itrk] << " - MCTrkId " << MCTrackID[itrk] << " MotherID " << MCMotherID[itrk] << endl;
			if(vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==411 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==421 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==431 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==4122)log_decay << "Nearby charm vertex with pdg " << vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]  << endl;;
			//(abs(mp_pdgID)==411 || abs(mp_pdgID)==421 || abs(mp_pdgID)==431 || abs(mp_pdgID)==4122)
			//trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
		        trk_pos.SetXYZ(ftrack->GetSegmentFirst()->X(),ftrack->GetSegmentFirst()->Y(),ftrack->Zstart());

		      
		       
			if(vID==tmp_vID)ty=1;
			else{
			  if(vtx_good_vec.at(tmp_index)==-1)ty=4;
			  else{
			    if(MCEventID[itrk]==mp_event)ty=2;
			    if(MCEventID[itrk]!=mp_event)ty=3;
			  }
			}
		       
			
			fdz = trk_pos(2) - tmp_vz;
			fdl = sqrt(pow(trk_pos(0) - tmp_vx,2)+pow(trk_pos(1) - tmp_vy,2)+pow(trk_pos(2) - tmp_vz,2));
			ip = IPtoVertex(vtx_pos, trk_pos, TX[itrk], TY[itrk]);
			// cout << "ip_vtx " << ip << endl;
			//cout << "ip_vtx " << ip << " " << vtx_pos(0) <<" " << vtx_pos(1) << vtx_pos(2) << " " << trk_pos(0) << " " << trk_pos(1)<< " " << trk_pos(2) << " " << TX[itrk] << " " << TY[itrk] <<   endl;
			ka = sqrt(pow(mc_tx - TX[itrk],2) + pow(mc_ty - TY[itrk],2));
			//rmax = FedraTrackKink(ftrack);
			FedraTrackKink(ftrack,dtheta);
			SegIPtoVertex(vtx_primary_proton_pos, ftrack, ipseg);
			pardaug = ParentDaugther(ftrack);
			FillTrack(nfound,icharm,0,MCTrackID[itrk],mc_pdg,nseg[itrk],ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
			nfound++;
			// }
		      }
		      
		    }
		    
		    int ndaug2 = charm_2_id.at(mp_event).size();
		    mc_tx = (c2vxMC[mp_event]-pvxMC[mp_event])/charm_2_z[mp_event];
		    mc_ty = (c2vyMC[mp_event]-pvyMC[mp_event])/charm_2_z[mp_event];
		  
		    
		    splitted_track=false;
		    for (int idaug = 0; idaug < ndaug2; idaug++) {   
		      //if(MCEventID[itrk]==765)cout << "linked "<<idaug << " " << MCTrackID[itrk] << " " << daug_1_id[mp_event][idaug] << endl;
		      if(MCTrackID[itrk]==daug_2_id[mp_event][idaug] && MCEventID[itrk]==fEventId[mp_event][idaug] && MCMotherID[itrk]==charm_2_id[mp_event][idaug]){

			icharm=22;
			trk_mc_id.push_back(MCTrackID[itrk]);
			trk_z.push_back(ftrack->Zstart());
			// cout <<"size "<< (trk_mc_id.size()-1) << endl;
			for(int s=0;s<(trk_mc_id.size()-1);s++){
			  //cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
			  if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){
			    //cout << "hello"<<endl;
			    log_decay << "Splitted track"<<endl;
			    splitted_track=true;
			    if(ftrack->Zstart()<trk_z.at(s))break;
			    else {
			      trk.fe_dz[s] = ftrack->Zstart()-tmp_vz;
			      //cout << "eccomi2 " << trk.fe_dz[s] << endl;
			      break;
			    }
			  }
			}
		
			//if(!splitted_track){ 
			log_decay << "Linked to another vertex " << vtx_list.at(i).at(j).at(k).at(l)  << " Daugther charm2 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
			log_decay << "MCEvt " << MCEventID[itrk] << " - MCTrkId " << MCTrackID[itrk] << " MotherID " << MCMotherID[itrk] << endl;

			if(vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==411 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==421 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==431 || vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]==4122)log_decay << "Nearby charm vertex with pdg " << vtx_mp_pdgID[vtx_list.at(i).at(j).at(k).at(l)]  << endl;;
		       
			//trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
		        trk_pos.SetXYZ(ftrack->GetSegmentFirst()->X(),ftrack->GetSegmentFirst()->Y(),ftrack->Zstart());
			if(vID==tmp_vID)ty=1;
			else{
			  if(vtx_good_vec.at(tmp_index)==-1)ty=4;
			  else{
			    if(MCEventID[itrk]==mp_event)ty=2;
			    if(MCEventID[itrk]!=mp_event)ty=3;
			  }
			}
			
			
			fdz = trk_pos(2) - tmp_vz;
			fdl = sqrt(pow(trk_pos(0) - tmp_vx,2)+pow(trk_pos(1) - tmp_vy,2)+pow(trk_pos(2) - tmp_vz,2));
			ip = IPtoVertex(vtx_pos, trk_pos, TX[itrk], TY[itrk]);
			ka = sqrt(pow(mc_tx - TX[itrk],2) + pow(mc_ty - TY[itrk],2));
			//rmax = FedraTrackKink(ftrack);
			FedraTrackKink(ftrack,dtheta);
			SegIPtoVertex(vtx_primary_proton_pos, ftrack, ipseg);
			pardaug = ParentDaugther(ftrack);
			FillTrack(nfound,icharm,0,MCTrackID[itrk],mc_pdg,nseg[itrk],ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
		
			nfound++;
		      }
		    }
		  }
		}
	      }
	    }
	     
	     
	    // TRACCE ISOLATE
	    
	    vtx_fedratree->GetEntry(ientry);
	     
	    for(int l=0;l< trk_list.at(i).at(j).at(k).size();l++){
	      trk_fedratree->GetEntry(trk_entry.at(i).at(j).at(k).at(l));

	      icharm=-1;
	      rmax=0;
	      pardaug=false;
	      // info dalla simulazione
	      if(t_->MCEvt()==mp_event){
	       
		simtree->GetEntry(t_->MCEvt());
	       
		for (Long64_t jn=0; jn<MCTrack_;jn++) {			
		  if(jn==t_->MCTrack()){
		    mc_pdg=MCTrack_fPdgCode[jn];
		    if(jn==1)icharm=1;
		    if(jn!=1 && (abs(mc_pdg)==411 || abs(mc_pdg)==421 || abs(mc_pdg)==431 || abs(mc_pdg)==4122))icharm=2;
		  }
		}

		///

		// SOLO CHARM
		if(icharm==1 || icharm==2){
		 
		  trk_mc_id.push_back(t_->MCTrack());
		  trk_z.push_back(t_->Z());
		  //cout <<"size "<< (trk_mc_id.size()-1) << endl;
		  for(int s=0;s<(trk_mc_id.size()-1);s++){
		    //cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
		    if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){ // 2 tracce in Fedra hanno lo stesso ID monte carlo
		      //cout << "hello"<<endl;
		      splitted_track=true;
		      if(t_->Z()<trk_z.at(s))break;
		      else {
			trk.fe_dz[s]=t_->Z()-tmp_vz;
			//cout << "eccomi3"<< endl;
			break;
		      }
		    }
		  }


		  if(icharm==1)log_decay << "Charm1 Track  in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		  if(icharm==2)log_decay << "Charm2 Track  in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		  log_decay << "MCEvt " <<  t_->MCEvt() << " - MCTrkId " << t_->MCTrack() << " MotherID " << t_->Aid(0) << endl;
		 
		  //trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
		  trk_pos.SetXYZ(t_->X(),t_->Y(),t_->Z());
		  ty=0;		    
		  fdz = t_->Z() - tmp_vz;
		  fdl = sqrt(pow(t_->X() - tmp_vx,2)+pow(t_->Y() - tmp_vy,2)+pow(t_->Z() - tmp_vz,2));
		  ip = IPtoVertex(vtx_pos, trk_pos, t_->TX(), t_->TY());
		  //cout << "ip_trk " << ip << " " << vtx_pos(0) <<" " << vtx_pos(1) << vtx_pos(2) << " " << trk_pos(0) << " " << trk_pos(1)<< " " << trk_pos(2) << " " << t_->TX() << " " << t_->TY() <<  endl;
		  ka = sqrt(pow(mc_tx - t_->TX(),2) + pow(mc_ty - t_->TY(),2));
		  //cout << ka << " " << pow(mc_tx - t_->TX(),2) << " " << pow(mc_ty - t_->TY(),2) << endl;
		  // cout << "nseg_fuori "<< tnseg << endl;
		  double seg_x[tnseg];
		  double seg_y[tnseg];
		  double seg_z[tnseg];
		  double seg_tx[tnseg];
		  double seg_ty[tnseg];
		  double seg_mc_trk[tnseg];
		  double seg_mc_evt[tnseg];
		  int trk_id = t_->ID();
		 
		  for(int k=0;k<s_;k++){
		    seg_x[k]=s_eX[k];
		    seg_y[k]=s_eY[k];
		    seg_z[k]=s_eZ[k];
		    seg_tx[k]=s_eTX[k];
		    seg_ty[k]=s_eTY[k];
		    seg_mc_trk[k]=s_eMCTrack[k];
		    seg_mc_evt[k]=s_eMCEvt[k];
		    //cout << k << " " << s_eTrack[k] << " " << tnseg << " " << s_eTX[k] << " " <<  s_eTY[k] << endl;
		  }
		  //rmax = FedraTrackKink2(tnseg,seg_tx,seg_ty);
		  FedraTrackKink2(tnseg,seg_tx,seg_ty,dtheta,trk_id);
		  SegIPtoVertex2(vtx_primary_proton_pos, tnseg, seg_x, seg_tx, seg_y, seg_ty, seg_z, ipseg, trk_id);
		  pardaug = ParentDaugther2(tnseg,seg_mc_trk, seg_mc_evt);
		  FillTrack(nfound,icharm,trid,t_->MCTrack(),mc_pdg,tnseg,ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
		  nfound++;		 
		}
		// fine solo charm
	       
		else{
		  int ndaug1 = charm_1_id.at(mp_event).size();
		  mc_tx = (c1vxMC[mp_event]-pvxMC[mp_event])/charm_1_z[mp_event];
		  mc_ty = (c1vyMC[mp_event]-pvyMC[mp_event])/charm_1_z[mp_event];
		
	       
		  splitted_track=false;
		  for (int idaug = 0; idaug < ndaug1; idaug++) {   
		    if(t_->MCTrack()==daug_1_id[mp_event][idaug] && t_->MCEvt()==fEventId[mp_event][idaug] && t_->Aid(0)==charm_1_id[mp_event][idaug]){
		   
		      icharm=11;
		      trk_mc_id.push_back(t_->MCTrack());
		      trk_z.push_back(t_->Z());
		      //cout <<"size "<< (trk_mc_id.size()-1) << endl;
		      for(int s=0;s<(trk_mc_id.size()-1);s++){
			//cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
			if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){ // 2 tracce in Fedra hanno lo stesso ID monte carlo
			  //cout << "hello"<<endl;
			  splitted_track=true;
			  if(t_->Z()<trk_z.at(s))break;
			  else {
			    trk.fe_dz[s]=t_->Z()-tmp_vz;
			    //cout << "eccomi3"<< endl;
			    break;
			  }
			}
		      }
		    
		    
		      //if(!splitted_track){
		      log_decay << "Daugther charm1 Track  in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		      log_decay << "MCEvt " <<  t_->MCEvt() << " - MCTrkId " << t_->MCTrack() << " MotherID " << t_->Aid(0) << endl;
		    
		      //trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
		      trk_pos.SetXYZ(t_->X(),t_->Y(),t_->Z());
		      ty=0;
		    
		      //cout << mc_tx << " " << mc_ty << " " << t_->TX() << " " << t_->TY() << endl;
		      //cout << "vtx " << vtx_pos(0) << " " << vtx_pos(1) << " " << vtx_pos(2);
		      //cout << "trk " << t_->X() << " " << t_->Y() << " " << t_->Z() << " " <<  t_->SX() << " " << t_->SY() << " " << t_->SZ() << endl;
		  
		    
		      fdz = t_->Z() - tmp_vz;
		      fdl = sqrt(pow(t_->X() - tmp_vx,2)+pow(t_->Y() - tmp_vy,2)+pow(t_->Z() - tmp_vz,2));
		      ip = IPtoVertex(vtx_pos, trk_pos, t_->TX(), t_->TY());
		      //cout << "ip_trk " << ip << " " << vtx_pos(0) <<" " << vtx_pos(1) << vtx_pos(2) << " " << trk_pos(0) << " " << trk_pos(1)<< " " << trk_pos(2) << " " << t_->TX() << " " << t_->TY() <<  endl;
		      ka = sqrt(pow(mc_tx - t_->TX(),2) + pow(mc_ty - t_->TY(),2));
		      //cout << ka << " " << pow(mc_tx - t_->TX(),2) << " " << pow(mc_ty - t_->TY(),2) << endl;
		      // cout << "nseg_fuori "<< tnseg << endl;
		      double seg_x[tnseg];
		      double seg_y[tnseg];
		      double seg_z[tnseg];
		      double seg_tx[tnseg];
		      double seg_ty[tnseg];
		      double seg_mc_trk[tnseg];
		      double seg_mc_evt[tnseg];
		      int trk_id = t_->ID();

		      //cout << "nseg " << tnseg << " " << s_ << endl;
		    
		      for(int k=0;k<s_;k++){
			seg_x[k]=s_eX[k];
			seg_y[k]=s_eY[k];
			seg_z[k]=s_eZ[k];
			seg_tx[k]=s_eTX[k];
			seg_ty[k]=s_eTY[k];
			seg_mc_trk[k]=s_eMCTrack[k];
			seg_mc_evt[k]=s_eMCEvt[k];
			//cout << k << " " << s_eTrack[k] << " " << tnseg << " " << s_eTX[k] << " " <<  s_eTY[k] << endl;
		      }
		      //rmax = FedraTrackKink2(tnseg,seg_tx,seg_ty);
		      FedraTrackKink2(tnseg,seg_tx,seg_ty,dtheta, trk_id);
		      SegIPtoVertex2(vtx_primary_proton_pos, tnseg, seg_x, seg_tx, seg_y, seg_ty, seg_z, ipseg, trk_id);
		      pardaug = ParentDaugther2(tnseg,seg_mc_trk, seg_mc_evt);
		      FillTrack(nfound,icharm,trid,t_->MCTrack(),mc_pdg,tnseg,ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
		      nfound++;
		      // }
		    }
		  }
	       
		
		  int ndaug2 = charm_2_id.at(mp_event).size();
		  mc_tx = (c2vxMC[mp_event]-pvxMC[mp_event])/charm_2_z[mp_event];
		  mc_ty = (c2vyMC[mp_event]-pvyMC[mp_event])/charm_2_z[mp_event];
	       

	       
		  splitted_track=false;
		  for (int idaug = 0; idaug < ndaug2; idaug++) {   
		    if(t_->MCTrack()==daug_2_id[mp_event][idaug] && t_->MCEvt()==fEventId[mp_event][idaug] && t_->Aid(0)==charm_2_id[mp_event][idaug]){

		      icharm=22;
		      trk_mc_id.push_back(t_->MCTrack());
		      trk_z.push_back(t_->Z());
		      //cout <<"size "<< (trk_mc_id.size()-1) << endl;
		      for(int s=0;s<(trk_mc_id.size()-1);s++){
			//cout << s << " " << trk_mc_id.at(s) << " " << trk_mc_id.at(trk_mc_id.size()-1) << " " << splitted_track << endl;
			if(trk_mc_id.at(s)==trk_mc_id.at(trk_mc_id.size()-1)){
			  //cout << "hello"<<endl;
			  splitted_track=true;
			  if(t_->Z()<trk_z.at(s))break;
			  else {
			    trk.fe_dz[s]=t_->Z()-tmp_vz;
			    //cout << "eccomi4"<<endl;
			    break;
			  }
			}
		      }
		   
		      //if(!splitted_track){
		      log_decay << "Daugther charm2 Track in "<< i << " " << j << " " << k << " " <<   l <<  " with TrkID " << trk_list.at(i).at(j).at(k).at(l) << endl;
		      log_decay << "MCEvt " <<  t_->MCEvt() << " - MCTrkId " << t_->MCTrack() << " MotherID " << t_->Aid(0) << endl;
		     
		      //trk_pos.SetXYZ(ftrack->X(),ftrack->Y(),ftrack->Z());
		      trk_pos.SetXYZ(t_->X(),t_->Y(),t_->Z());
		      ty=0;
		     
		      fdz = t_->Z() - tmp_vz;
		      fdl = sqrt(pow(t_->X() - tmp_vx,2)+pow(t_->Y() - tmp_vy,2)+pow(t_->Z() - tmp_vz,2));
		      ip = IPtoVertex(vtx_pos, trk_pos, t_->TX(), t_->TY());
		      ka = sqrt(pow(mc_tx - t_->TX(),2) + pow(mc_ty - t_->TY(),2));

		      double seg_x[tnseg];
		      double seg_y[tnseg];
		      double seg_z[tnseg];
		      double seg_tx[tnseg];
		      double seg_ty[tnseg];
		      double seg_mc_trk[tnseg];
		      double seg_mc_evt[tnseg];
		      int trk_id = t_->ID();
		     
		      for(int k=0;k<s_;k++){
			seg_x[k]=s_eX[k];
			seg_y[k]=s_eY[k];
			seg_z[k]=s_eZ[k];
			seg_tx[k]=s_eTX[k];
			seg_ty[k]=s_eTY[k];
			seg_mc_trk[k]=s_eMCTrack[k];
			seg_mc_evt[k]=s_eMCEvt[k];
			//cout << k << " " << s_eTrack[k] << " " << tnseg << " " << s_eTX[k] << " " <<  s_eTY[k] << endl;
		      }
		      //rmax = FedraTrackKink2(tnseg,seg_tx,seg_ty);
		      FedraTrackKink2(tnseg,seg_tx,seg_ty,dtheta,trk_id);
		      SegIPtoVertex2(vtx_primary_proton_pos, tnseg, seg_x, seg_tx, seg_y, seg_ty, seg_z, ipseg, trk_id);
		      pardaug = ParentDaugther2(tnseg,seg_mc_trk, seg_mc_evt);
		      FillTrack(nfound,icharm,trid,t_->MCTrack(),mc_pdg,tnseg,ty,fdz,fdl,ip,ka,dtheta[0],dtheta[1],ipseg[0],ipseg[1],splitted_track,pardaug);
		      nfound++;
		      // }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
      TreeDS->Fill();
    }
  }
   
  TreeDS->Write();
  f_ds->Close();
  log_decay << "DECAY SEARCH LOG TERMINED" << endl;
  log_decay.close();
   
  ofstream log2;
  log2.open("list_vertex_and_tracks.txt",ios::out);
  log2 << "//----V. GENTILE ------// \n LIST OF VERTEX (6th 0) AND TRACKS (6th 1)" << endl;
   
  for(int i=0;i<width;i++){
    for(int j=0;j<height;j++){
      for(int k=0;k<depth;k++){
	//cout << "size " << vtx_list.at(i).at(l).at(k).size() << endl;
	for(int l=0;l<vtx_list.at(i).at(j).at(k).size();l++){
	  if(vtx_list.at(i).at(j).at(k).size()>0){
	    //cout << i << " " << j << " " << k << " " << l <<  endl;
	    log2 << i << " " << j << " " << k << " " << l <<  " " << vtx_list.at(i).at(j).at(k).at(l) << " " << "0" <<  endl;
	  }
	}
	for(int l=0;l< trk_list.at(i).at(j).at(k).size();l++){
	  log2 << i << " " << j << " " << k << " " << l <<  " " << trk_list.at(i).at(j).at(k).at(l) << " " << "1" <<  endl;
	}
      }
    }
  }

  log2.close();
  
   
}

int myrun(){
  Loop();
  return 0;
}
