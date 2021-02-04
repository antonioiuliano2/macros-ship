// Classificazione dei vertici simulati ricostruiti in FEDRA
// Training per analisi multivariata
// V. Gentile 2019

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

#include <algorithm>    // std::equal
#include <vector>       // std::vector

using namespace std;


void Loop()
{

  int cut_ntrk;
  cout << "Insert cut on ntracks ";
  cin >> cut_ntrk;
  cout << endl;

  int cut_vz_min;
  cout << "Insert cut on vz_min ";
  cin >> cut_vz_min;
  cout << endl;

  ofstream log_vtx;
  log_vtx.open("log_vtx_search.txt",ios::out);
  log_vtx << "//----V. GENTILE 10/03/2019 ------// \n VERTEX SEARCH ON MC SIMULATION" << endl;
  
  ofstream log_decay;
  log_decay.open("log_data_search.txt",ios::out);
  log_decay << "//----V. GENTILE ------// \n VERTEX SEARCH ON DATA" << endl;

  TFile * fedrafile = TFile::Open("vertextree.root");  // full
  //TFile * simulationfile = TFile::Open("../ship.conical.Pythia8CharmOnly-TGeant4_dig.root");
  TFile * simulationfile = TFile::Open("../../pythia8_Geant4_1000_0.1.root");

  //TFile * bdtfile = TFile::Open("vtx_BDT_data_evaluated.root");  // full
  //TTree * bdttree = (TTree*) bdtfile->Get("bdt");
  //TTree * vtxbdttree = (TTree*) bdtfile->Get("vtx");

  //bdt_reader(bdttree);
  //vtxbdt_reader(vtxbdttree);
  
  TTree *simtree = (TTree*) simulationfile->Get("cbmsim");     
  TTree *vtx_fedratree = (TTree*) fedrafile->Get("vtx");

  if(vtx_fedratree == 0) return;
      
  vtx_reader_Fedra(vtx_fedratree);
  cbmsim_reader_MC(simtree);

  
  Long64_t nentries_vtx_fedra = vtx_fedratree->GetEntriesFast();
  log_decay<<"tot entries vtx fedra "<< nentries_vtx_fedra << endl;
   
  vtx_fedratree->BuildIndex("vID");


  int vbinx=13;//52;
  int vbiny=10;//24;
  int vbinz=8;//16;

  double epsilon=1; //micron
  double epsilon_phi=0.1; //rad
  int min_points=1;
  
  TH3D * hvertex = new TH3D("hv","hv",vbinx,0,130000,vbiny,0,100000,vbinz,-77000,3000); //WARNING

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

  log_decay << "\nxrange [" << xmin << " , " << xmax << "] ; yrange [" << ymin << " , " << ymax << "] ; zrange [" << zmin << " , " << zmax << "]" << endl;
  log_decay << "Cube dimensions XYZ (" << xrange << " , " << yrange << " , " << zrange << ")" << endl;
  log_decay << "N cubes in XYZ "<< width << " " << height << " " << depth << endl;
  

  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > vtx_list;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > vtx_list_nearby;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > trk_list;
  std::vector<std::vector<std::vector<std::vector<unsigned int> > > > trk_entry;   
  vtx_list.resize(width);
  vtx_list_nearby.resize(width);
  trk_list.resize(width);
  trk_entry.resize(width);
  for(int i=0;i<width;i++){
    //y axis size
    vtx_list[i].resize(height);
    vtx_list_nearby[i].resize(height);
    trk_list[i].resize(height);
    trk_entry[i].resize(height);
    for(int j=0;j<height;j++){
      //z axis size
      vtx_list[i][j].resize(depth);
      vtx_list_nearby[i][j].resize(depth);
      trk_list[i][j].resize(depth);
      trk_entry[i][j].resize(depth);
    }
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
  int charm_daug=0;
  int mp_ev_pdgID=0;
  int n_electrons=0;
  int n_primaries=0;
  int n_charmdaug=0;
  int n_events=0;
  int n_mothers=0;

  bool flag_trk=0;
  bool flag_ev=0;
  bool vtx_good=false;
  

  float mp_vx=0;
  float mp_vy=0;
  float mp_vz=0;
  float mp_vdx=-1;
  float mp_vdy=-1;
  float mp_vdz=-1;
  
  TFile * f_out = new TFile("vtx_MC_analysis.root","RECREATE");
  TTree * Tree_out = new TTree();
  Tree_out = vtx_fedratree->CloneTree(0);
  Tree_out->Branch("vID",&vID,"vID/I");
  Tree_out->Branch("vtx_good",&vtx_good,"vtx_good/B");
  Tree_out->Branch("mp_eventID",&mp_eventID,"mp_eventID/I");
  Tree_out->Branch("mp_motherID",&mp_motherID,"mp_motherID/I");
  Tree_out->Branch("mp_procID",&mp_procID,"mp_procID/I");
  Tree_out->Branch("mp_pdgID",&mp_pdgID,"mp_pdgID/I");
  Tree_out->Branch("charm_daug",&charm_daug,"charm_daug/O");
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
  Tree_out->Branch("n_primaries",&n_primaries,"n_primaries/I");
  Tree_out->Branch("n_charmdaug",&n_charmdaug,"n_charmgaug/I");
  Tree_out->Branch("n_events",&n_events,"n_events/I");
  Tree_out->Branch("n_mothers",&n_mothers,"n_mothers/I");
  Tree_out->AddFriend(vtx_fedratree);
  
  
   std::vector<int> vtx_good_vec(nentries_vtx_fedra);
   std::vector<int> vtx_mp_eventID(nentries_vtx_fedra);
   std::vector<int> vtx_mp_motherID(nentries_vtx_fedra);


   //bdttree->BuildIndex("bdt_vID");

   // VERTEX STUDY
   
   for (Long64_t ientry=0; ientry<nentries_vtx_fedra;ientry++) { // Vertex index
     //cout << ientry << " " << flag << endl;
     if(ientry>=0){ // for debugging
       if(ientry%100==0)cout << "VID " << ientry << endl;
       log_vtx << "\nFedra Vertex number "<< ientry << endl;
       log_vtx << "MC_EventID / MCtrkID / MCMotherID / MCpdgCode / MCMotherpdgCode" << endl;
       vtx_fedratree->GetEntry(ientry);
       if(flag!=2 && flag!=5){
	 vtx_good_vec[ientry]=-1;
	 vtx_mp_eventID[ientry]=-1;
	 vtx_mp_motherID[ientry]=-1;
	 charm_daug=false;
	 	 
	 // ONLY FOR INTERACTION VTX IDENTIFICATION
	 /*int tmp_index = -1;
	 float tmp_bdt = -1;
	 tmp_index = bdttree->GetEntryNumberWithIndex(vID); // indice sul vertice secondario
	 if(tmp_index!=-1){
	   bdttree->GetEntry(tmp_index);
	   tmp_bdt = bdt_value;
	   cout << "vID " << vID << " " << ientry << " " << tmp_index << " " << tmp_bdt << endl;
	 }*/
	 
	 vtx_fedratree->GetEntry(ientry);	 
	 /////
	 
	 if(vx>xmin && vx<xmax && vy>ymin && vy<ymax && vz>zmin && vz<zmax && n>=cut_ntrk && vz<-4000 && vz>cut_vz_min){ // LAST 4 CUTS ONLY FOR INTERACTION VTX IDENTIFICATION
	   
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
	   
	   vtx_list.at(ix).at(iy).at(iz).push_back(vID);
	   
	   for (int itrk=0; itrk<n;itrk++) {  // Track index
	  
	     if(MCMotherID[itrk]==1 || MCMotherID[itrk]==2)charm_daug=true;
	  
	     simtree->GetEntry(MCEventID[itrk]);
	     
	     for (Long64_t jn=0; jn<MCTrack_;jn++) {
	       
	       if(jn==MCTrackID[itrk]){
		 vtx_trackeventId.push_back(MCEventID[itrk]);
		 vtx_trackpdgcode.push_back(MCTrack_fPdgCode[jn]);
		 vtx_trackmother.push_back(MCTrack_fMotherId[jn]);
		 vtx_trackprocess.push_back(MCTrack_fProcID[jn]);
		 vtx_trackid.push_back(jn);
		 vtx_trackstartX.push_back(MCTrack_fStartX[jn]);
		 vtx_trackstartY.push_back(MCTrack_fStartY[jn]);
		 vtx_trackstartZ.push_back(MCTrack_fStartZ[jn]);
	       }
	     }
	     
	     for (Long64_t jn=0; jn<EmuBaseTrks_;jn++) {	 
	       if(jn==MCTrackID[itrk]){
		 vtx_trackstartX.push_back(EmuBaseTrks_fX[jn]);
		 vtx_trackstartY.push_back(EmuBaseTrks_fY[jn]);
	       }
	     }       
	   }
	   
	   
	   for(uint h=0; h<vtx_trackmother.size(); h++){
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
	     
	     if(abs(vtx_motherpdgcode.at(h))==411)log_vtx << "Figlio del D+ " << endl;
	     if(abs(vtx_motherpdgcode.at(h))==421)log_vtx << "Figlio del D0 " << endl;
	     if(abs(vtx_motherpdgcode.at(h))==431)log_vtx << "Figlio del Ds " << endl;
	     if(abs(vtx_motherpdgcode.at(h))==441)log_vtx << "Figlio dell'Eta_c " << endl;
	     if(abs(vtx_motherpdgcode.at(h))==4122)log_vtx << "Figlio del Lambda_c+ " << endl;
	     if(abs(vtx_motherpdgcode.at(h))==4132)log_vtx << "Figlio della Xsi_c0 " << endl;
	     if(abs(vtx_motherpdgcode.at(h))==4232)log_vtx << "Figlio della Xsi_c+ " << endl;
	     if(abs(vtx_motherpdgcode.at(h))==4332)log_vtx << "Figlio dell'Omega_c0 " << endl;
	     log_vtx << vtx_trackeventId.at(h) <<  " "  << vtx_trackid.at(h) << " " << vtx_trackmother.at(h) << " " << vtx_trackpdgcode.at(h) << " " << vtx_motherpdgcode.at(h) << endl;
	     
	   }
     
     
	   std::vector<int> vtx_ntracks_same_mother(n);
	   //vtx_ntracks_same_mother.clear();
	   std::vector<int> vtx_ntracks_same_event(n);
	   //vtx_ntracks_same_event.clear();
	   
	   /// CONTA LE OCCORRENZE PER MOTHER
	   for(uint h=0; h<vtx_trackmother.size(); h++){
	     flag_trk=0;
	     for (uint i=0;i<h;i++){
	       if (vtx_trackmother[h]==vtx_trackmother[i]){// controlla che quel numero non sia già stato considerato in precedenza
		 flag_trk=1;
		 break;
	       }
	       //cout << "before "<< ientry << " " << h << " " << i << " " << flag << " " << vtx_trackmother[h] << " " << vtx_ntracks_same_mother[h] << endl;
	     }
	     if (flag_trk==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	       vtx_ntracks_same_mother[h]=1; // assegna 1 alla mother della prima traccia che incontra
	       for (uint i=h+1;i<vtx_trackmother.size();i++){		 
		 if(vtx_trackmother[h] == vtx_trackmother[i] && vtx_trackid[h]!=vtx_trackid[i] && vtx_trackeventId[h]==vtx_trackeventId[i]) vtx_ntracks_same_mother[h]++;
	       }
	       
	     }
	   }
	   
     
	   /// CONTA LE OCCORRENZE PER EVENTO
	   for(uint h=0; h<vtx_trackmother.size(); h++){
	     flag_ev=0;
	     for (uint i=0;i<h;i++){
	       if (vtx_trackeventId[h]==vtx_trackeventId[i]){// controlla che quel numero non sia già stato considerato in precedenza
		 flag_ev=1;
		 break;
	       }
	     }
	     if (flag_ev==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	       vtx_ntracks_same_event[h]=1;
	       for (uint i=h+1;i<vtx_trackmother.size();i++){		 
		 if(vtx_trackeventId[h] == vtx_trackeventId[i] && vtx_trackid[h]!=vtx_trackid[i]) vtx_ntracks_same_event[h]++;
	       }
	     }
	   }
	   
	   
	   
	   n_electrons=0;
	   n_primaries=0;
	   n_charmdaug=0;
	   n_events=0;
	   n_mothers=0;

	   
	   // counts the number of events linked to the vertex
	   for(uint i=0;i<vtx_ntracks_same_event.size();i++){
	     if(vtx_ntracks_same_event[i]>0)n_events++;
	   }
	   
	   // counts the number of mothers linked to the vertex
	   for(uint i=0;i<vtx_ntracks_same_mother.size();i++){
	     if(vtx_ntracks_same_mother[i]>0)n_mothers++;
	   }
	   
	   for(uint i=0; i<vtx_trackpdgcode.size(); i++){
	     if(abs(vtx_trackpdgcode.at(i))==11)n_electrons++;
	     if(vtx_trackmother.at(i)==-1)n_primaries++;
	     if((vtx_trackmother.at(i)==1 || vtx_trackmother.at(i)==2))n_charmdaug++;
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
	   
	   
	   if(max_ev_occurrence<=1){
	     mp_eventID=-1;
	   }
	   else {
	     //cout << "iel_ev " << max_occurrence << " " << iel_ev << endl;
	     mp_eventID = vtx_trackeventId.at(iel_ev);
	     //for(uint i=0; i<vtx_trackpdgcode.size(); i++){
	     //if(abs(vtx_trackpdgcode.at(i))==11 /*&& vtx_trackeventId.at(i)==mp_eventID*/)n_electrons++;
	     //if(vtx_trackmother.at(i)==-1 && vtx_trackeventId.at(i)==mp_eventID)n_primaries++;
	     //if((vtx_trackmother.at(i)==1 || vtx_trackmother.at(i)==2) && vtx_trackeventId.at(i)==mp_eventID)n_charmdaug++;
	     //}
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
	   if(max_occurrence>1){
	     
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
	     
	   }
	   mean_freq=(1.0*max_occurrence)/n;
	   max_freq=max_occurrence;
	   
	   mean_ev_freq=(1.0*max_ev_occurrence)/n;
	   max_ev_freq=max_ev_occurrence;
	   
	   log_vtx << "The max number of tracks with same mother ID "<< mp_motherID << " in the Vertex " << ientry << " is " << max_occurrence << endl;
	   log_vtx << "The max number of tracks with same event ID "<< mp_eventID << " in the Vertex " << ientry << " is " << max_ev_occurrence << endl;
	   
	   
	   //------ DESCRIPTION  -----------//
	   
  
	   if(mp_pdgID!=-2){
	     if(mp_motherID==-1){
	       log_vtx << "It's a primary proton vertex with pdg code " << mp_pdgID << endl;
	       vtx_good=true;
	       vtx_good_vec[ientry]=0;
	       vtx_mp_eventID[ientry]=mp_eventID;
	       vtx_mp_motherID[ientry]=mp_motherID;
	     }
	     if(mp_motherID==0){
	       log_vtx << "It's a charm vertex with pdg code " << mp_pdgID << endl;
	       vtx_good=true;
	       vtx_good_vec[ientry]=1;
	       vtx_mp_eventID[ientry]=mp_eventID;
	       vtx_mp_motherID[ientry]=mp_motherID;
	     }
	     if( mp_motherID>0){
	       log_vtx << "It's a  secondary vertex with pdg code " << mp_pdgID << endl;
	       vtx_good=true;
	       vtx_good_vec[ientry]=2;
	       vtx_mp_eventID[ientry]=mp_eventID;
	       vtx_mp_motherID[ientry]=mp_motherID;
	     }
	     
	   }
	   else {
	     log_vtx << "It's a fake\n";
	     vtx_good=false;
	   }
	   
	   //------------------------- //
	   
	   log_vtx << "VID / MostProbable Vertex X (Y,Z) from MC / Vertex X (Y,Z) from Fedra\n";
	   log_vtx << ientry << " " << mp_vx << " " << vx << " " << mp_vy << " " << vy << " " << mp_vz << " " << vz << endl;
	   log_vtx << "Gap between MC and Fedra Vertex position (X,Y,Z)\n";
	   log_vtx << ientry <<" " <<  mp_vdx << " " << mp_vdy << " " << mp_vdz << endl;
	   if(charm_daug==true)log_vtx << "A charm daughter is present" << endl;
	   if(n_electrons>0)log_vtx << "The number of electron tracks is " <<n_electrons <<  endl;
	   if(n_primaries>0)log_vtx << "The number of primary tracks is " <<n_primaries <<  endl;
	   if(n_charmdaug>0)log_vtx << "The number of daughter charm tracks is " <<n_charmdaug <<  endl;
	   if(n_events>1)log_vtx << "Warning: the number of linked montecarlo events are " <<n_events <<  endl;
	   if(n_mothers>1)log_vtx << "Warning: the number of linked montecarlo mothers are " <<n_mothers <<  endl;
	   
	   Tree_out->Fill();
	   if(ientry%1000==0)Tree_out->AutoSave("SaveSelf");
	   
	 }
       } // flag
     } // only for special vertex
   }
   log_vtx << "VERTEXING LOG TERMINED" << endl;
   log_vtx.close();
   Tree_out->Write();
   f_out->Close();
   
   
}

int myrun(){
  Loop();
  return 0;
}
