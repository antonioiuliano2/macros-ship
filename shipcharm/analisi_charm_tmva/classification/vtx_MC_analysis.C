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

using namespace std;

void Loop()
{

  int fa,fb,fc;
  vector <int> fMolt;
  vector <int> fTrkId;
  vector <int> fEventId;

  vector <int> fsMolt;
  vector <int> fsTrkId;
  vector <int> fsEventId;
  
  ifstream input;
  input.open("interacting.txt",ios::in);
  while(1){
    input >> fa >> fb >> fc;
    if(input.eof())break;
    if(!input.good())break;
    fMolt.push_back(fa);
    fTrkId.push_back(fb);
    fEventId.push_back(fc);
  }

  TH1I *h = new TH1I("","",200,0,100);
  ifstream input2;
  input2.open("secondaries_halfspill.txt",ios::in);
  while(1){
    input2 >> fa >> fb >> fc;
    if(input2.eof())break;
    if(!input2.good())break;
    h->Fill(fa);
    fsMolt.push_back(fa);
    fsTrkId.push_back(fb);
    fsEventId.push_back(fc);
  }
  //h->Draw("");
  
  ofstream log;
  log.open("log_vtx.txt",ios::out);
  log << "//----V. GENTILE 10/03/2019 ------// \n FEDRA VERTEXING ON MC SIMULATION" << endl;
  
  TFile * simulationfile = TFile::Open("pythia8_Geant4_1000_0.5.root");
  //TFile * fedrafile = TFile::Open("vtx_MC_root6_small_1190ev.root"); // small file
  TFile * fedrafile = TFile::Open("vtx_MC_16_04_19.root");  // full

  TTree *simtree = (TTree*) simulationfile->Get("cbmsim");
  TTree *fedratree = (TTree*) fedrafile->Get("vtx");

  
  std::vector<bool> event_primary_proton;
  std::vector<int> vtx_trackmother;
  std::vector<int> vtx_trackprocess;
  std::vector<int> vtx_trackpdgcode;
  std::vector<int> vtx_motherpdgcode;
  std::vector<int> vtx_trackid;
  std::vector<float> vtx_trackstartX;
  std::vector<float> vtx_trackstartY;
  std::vector<float> vtx_trackstartZ;
  std::vector<float> vtx_trackeventId;

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
  bool primary_proton=false;
  bool daughter_proton=false;
  bool vtx_good=false;
  

  float mp_vx=0;
  float mp_vy=0;
  float mp_vz=0;
  float mp_vdx=-1;
  float mp_vdy=-1;
  float mp_vdz=-1;

  TH1F *hoccurrence = new TH1F("occurrence / ntracks","",50,0,1);
  
  //  std::vector <std::vector<double> > cl_el;
  //  ev_x_pos2.clear();
  TFile * f_out = new TFile("vtx_MC_analysis_15_04_19.root","RECREATE");
  //fedratree->CloneTree()->Write();
  TTree * Tree_out = new TTree();
  //Tree_out =(TTree*)fedratree->CopyTree("");
  fedratree->CloneTree()->Write();
  Tree_out->SetName("ana");
  Tree_out->Branch("vtx_good",&vtx_good,"vtx_good/B");
  Tree_out->Branch("primary_proton",&primary_proton,"primary_proton/B");
  Tree_out->Branch("daughter_proton",&daughter_proton,"daughter_proton/B");
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
  Tree_out->AddFriend(fedratree);

  vtx_reader_Fedra(fedratree);
  vtx_reader_MC(simtree);
  
  if(fedratree == 0) return;
  if(simtree == 0) return;

   Long64_t nentries_fedra = fedratree->GetEntriesFast();
   log<<"tot entries fedra "<< nentries_fedra << endl;

   Long64_t nentries_MC = simtree->GetEntriesFast();
   log<<"tot entries MC "<< nentries_MC << endl;

   /// Inizializzo a false tutti gli eventi primary proton
   for (Long64_t ientry=0; ientry<nentries_MC;ientry++) { // SIM index
     event_primary_proton.push_back(false);
   }

   // Sono true gli eventi nella lista interacting.txt
   for (Long64_t ientry=0; ientry<fEventId.size();ientry++) { // file index
     event_primary_proton[fEventId.at(ientry)]=true; 
   }

   /*
   // Sono true gli eventi nella lista secondaries.txt
   for (Long64_t ientry=0; ientry<fsEventId.size();ientry++) { // file index
     simtree->GetEntry(fsEventId.at(ientry));
     if(MCTrack_fPdgCode[fsTrkId.at(ientry)]!=2112)cout << MCTrack_fPdgCode[fsTrkId.at(ientry)] << endl;
     }*/

   
   for (Long64_t ientry=0; ientry<nentries_fedra;ientry++) { // Vertex index

     if(ientry%100==0)cout << "VID " << ientry << endl;
     log << "\nFedra Vertex number "<< ientry << endl;
     log << "MC_EventID / MCtrkID / MCMotherID / MCpdgCode / MCProcessId" << endl;
     fedratree->GetEntry(ientry);

     vtx_trackmother.clear();
     vtx_trackprocess.clear();
     vtx_trackpdgcode.clear();
     vtx_motherpdgcode.clear();
     vtx_trackid.clear();
     vtx_trackstartX.clear();
     vtx_trackstartY.clear();
     vtx_trackstartZ.clear();
     vtx_trackeventId.clear();

     for (int itrk=0; itrk<n;itrk++) {  // Track index

       vtx_trackeventId.push_back(MCEventID[itrk]);
       simtree->GetEntry(MCEventID[itrk]);
       
       for (Long64_t jn=0; jn<MCTrack_;jn++) {
	 
	 if(jn==MCTrackID[itrk]){
	   //cout << MCEventID[itrk]  << " "  << jn << " " << MCTrack_fMotherId[jn] << " " << MCTrack_fPdgCode[jn] << " " << MCTrack_fProcID[jn] << endl;
	   log << MCEventID[itrk]  << " "  << jn << " " << MCTrack_fMotherId[jn] << " " << MCTrack_fPdgCode[jn] << " " << MCTrack_fProcID[jn] << endl;
	   vtx_trackpdgcode.push_back(MCTrack_fPdgCode[jn]);
	   vtx_motherpdgcode.push_back(MCTrack_fPdgCode[MCTrack_fMotherId[jn]]);
	   vtx_trackmother.push_back(MCTrack_fMotherId[jn]);
	   vtx_trackprocess.push_back(MCTrack_fProcID[jn]);
	   vtx_trackid.push_back(jn);
	   vtx_trackstartX.push_back(MCTrack_fStartX[jn]);
	   vtx_trackstartY.push_back(MCTrack_fStartY[jn]);
	   vtx_trackstartZ.push_back(MCTrack_fStartZ[jn]);	   	   
	 }
       }
     }

     std::vector<int> vtx_ntracks_same_mother(n);
     //vtx_ntracks_same_mother.clear();
     std::vector<int> vtx_ntracks_same_event(n);
     //vtx_ntracks_same_event.clear();
     
     /// CONTA LE OCCORRENZE PER MOTHER
     for(int h=0; h<n; h++){
       flag=0;
       for (int i=0;i<h;i++){
	 if (vtx_trackmother[h]==vtx_trackmother[i]){// controlla che quel numero non sia già stato considerato in precedenza
	   flag=1;
	   break;
	 }
       }
	 if (flag==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	   vtx_ntracks_same_mother[h]=1;
	   for (int i=h+1;i<n;i++){		 
	     if(vtx_trackmother[h] == vtx_trackmother[i] && vtx_trackid[h]!=vtx_trackid[i] && vtx_trackeventId[h]==vtx_trackeventId[i]) vtx_ntracks_same_mother[h]++;
	   }
	   //cout << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackmother[h] << "trk Frequency " << vtx_ntracks_same_mother[h] << " times" << endl;
	   log << "MC_EventID " << vtx_trackeventId[h] << "ev MCmotherID " << vtx_trackmother[h] << "trk Frequency " << vtx_ntracks_same_mother[h] << " times" << endl;
	   }
	 //cout << h << " " << vtx_ntracks_same_mother[h] << endl;
     }


     /// CONTA LE OCCORRENZE PER EVENTO
     for(int h=0; h<n; h++){
       flag_ev=0;
       for (int i=0;i<h;i++){
	 if (vtx_trackeventId[h]==vtx_trackeventId[i]){// controlla che quel numero non sia già stato considerato in precedenza
	   flag_ev=1;
	   break;
	 }
       }
       if (flag_ev==false){ // se il numero non è stato considerato in precedenza, allora conta quante volte compare nel vettore e visualizza il numero di comparse
	   vtx_ntracks_same_event[h]=1;
	   for (int i=h+1;i<n;i++){		 
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


     n_electrons=0;
     daughter_proton=false;
    
     if(max_ev_occurrence==1) {
       mp_eventID=-1;
       primary_proton=false;
     }
     else {
       mp_eventID = vtx_trackeventId.at(iel_ev);
       //cout << ientry << " " << mp_eventID << endl;
       if(event_primary_proton.at(mp_eventID)==true){
	 primary_proton=true;
	 for(int i=0; i<vtx_trackpdgcode.size(); i++){
	   if(abs(vtx_trackpdgcode.at(i))==2212 && vtx_motherpdgcode.at(i)==-1)daughter_proton=true;
	   if(abs(vtx_trackpdgcode.at(i))==11 && vtx_trackeventId.at(i)==mp_eventID)n_electrons++;
	 }
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
     else {
              
       int trk_index=0;
             
       for(unsigned int itrk=0;itrk<vtx_trackmother.size();itrk++){
	 
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

     log << "The max number of tracks with same mother ID "<< mp_motherID << " in the Vertex " << ientry << " is " << max_occurrence << endl;
     log << "The max number of tracks with same event ID "<< mp_eventID << " in the Vertex " << ientry << " is " << max_ev_occurrence << endl;


     //------ DESCRIPTION  -----------//
     
     if(primary_proton){
       if(mp_pdgID!=-2){
	 if(abs(mp_pdgID)==2212 && mp_motherID==0){
	   log << "It's a primary proton\n";
	   vtx_good=true;
	 }
	 else{
	   if(!daughter_proton){
	     log << "Is it a secondary with pdg code " << mp_pdgID << "?\n";
	     vtx_good=true;
	   }
	   else {
	     log << "Warning! A primary proton is a daughter!\n";
	     vtx_good=false;
	   }
	 }
       }
       else {	 	 
	 if(!daughter_proton) log << "Is it a fake or a shower with " << n_electrons << " electrons?\n";
	 else log << "Warning! A primary proton is a daughter!\n";
	 vtx_good=false;
       }
     }
     else {
       log << "It's a fake\n";
       vtx_good=false;
     }

     //------------------------- //
     
       //if(daughter_proton)log << "Warning! A primary proton is a daughter!\n";
     //if(mp_motherID==0 && mp_pdgID==2212)log << "It's a primary proton\n";
     //if(mp_motherID>0 && mp_pdgID==2212)log << "It's a secondary proton\n";
     //if(mp_motherID==-1 && mp_pdgID==2212)log << "Warning! A primary proton is a daughter!\n";

     log << "VID / MostProbable Vertex X (Y,Z) from MC / Vertex X (Y,Z) from Fedra\n";
     log << ientry << " " << mp_vx << " " << vx << " " << mp_vy << " " << vy << " " << mp_vz << " " << vz << endl;
     log << "Gap between MC and Fedra Vertex position (X,Y,Z)\n";
     log << ientry <<" " <<  mp_vdx << " " << mp_vdy << " " << mp_vdz << endl;
     
     hoccurrence->Fill(mean_freq);
     Tree_out->Fill();
   }


   log << "VERTEXING LOG TERMINED" << endl;
   hoccurrence->Draw();
   log.close();
   Tree_out->Write();
   f_out->Close();
   
}

int myrun(){
  Loop();
  return 0;
}
