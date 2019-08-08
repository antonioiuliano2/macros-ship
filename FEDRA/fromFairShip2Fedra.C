//tool to load hit from emulsion to FEDRA (13 aprile 2018) //updated to work with the new TTreeReader structure (6 March)
//to use it, go in a directory and create the folders b000001/p001 to b000001/p029
//then launch it from the directory mother of b000001

#include <stdio.h>
#include <TROOT.h>
#include "TRandom.h"

void fromFairShip2Fedra(TString filename);
void smearing (Float_t &TX, Float_t &TY, const float angres);
bool efficiency(const float emuefficiency);
bool efficiency(const float tantheta, TH1D * emuefficiency);

//start script
void fromFairShip2Fedra(){
 
 fromFairShip2Fedra("ship.conical.Pythia8CharmOnly-TGeant4.root");
}

TF1 angularresolution(){
  //estimate angular resolution from OPERA data (Nuclear Instruments and Methods in Physics Research A 554 (2005) 247–254, M.De Serio et al.)
 const int npoints = 3;
 float angles[npoints] = {0.05, 0.2, 0.3}; //angles of tracks from which measurements have been performed
 float accuracy[npoints] = {0.0004, 0.00064,0.00094}; //angular accuracy
  
 TF1 fres = TF1 ("fres","pol1",0,0.3); //resolution function to interpolate
 
 TGraph *resgraph = new TGraph(npoints, angles, accuracy);
  
 TCanvas *rescanvas = new TCanvas();
 resgraph->Draw("AP*");
 resgraph->Fit(&fres);
 rescanvas->Print("resfunction.root");
 return fres;
}

//#include "/home/utente/fedra/include/EdbCouplesTree.h"
using namespace TMath;
TRandom *grandom = new TRandom3(); //creating every time a TRandom3 is a bad idea
TFile *file = NULL;
TH1D *heff = NULL ; //efficiency at different angles
/*
void EmuTracksfromFairShip2Fedra(TString filename){ //directly read emutracks, no corrections here because they are done directly at digitization level
 const int nplates = 29;
 int nbrick = 1;

 //**********************OPENING INPUT FILE***************************
 TFile * inputfile = TFile::Open(filename.Data());
 //TFile * inputfile = TFile::Open("/eos/user/a/aiuliano/sims_FairShip/sim_charm/pot/uniformonespill_onelayer_ch1_03_03_19/pythia8_Geant4_1000_0.5.root");
 if (!inputfile) return;

 //getting tree and arrays
 TTreeReader reader("cbmsim",inputfile);
 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<EmuBaseTrk> emulsionhits(reader,"EmuBaseTrks");
 
 //TTree* cbmsim = (TTree*)inputfile->Get("cbmsim");
 Float_t tx = 0, ty=0, xem= 0, yem = 0;
 const Int_t nevents = reader.GetTree()->GetEntries();
 int ihit = 0, ievent = 0;
 int nfilmhit = 0;
 float tantheta;
 int trackID = 0, motherID = 0, pdgcode = 0;
 //***********************CREATING FEDRA TREES**************************
 gInterpreter->AddIncludePath("/afs/cern.ch/work/a/aiuliano/public/fedra/include");
 EdbCouplesTree *ect[nplates];
 for (int i = 1; i <= nplates; i++){
  ect[i-1] = new EdbCouplesTree();
  if (i <10) ect[i-1]->InitCouplesTree("couples",Form("b00000%i/p00%i/1.%i.0.0.cp.root",nbrick,i,i),"RECREATE");
  else ect[i-1]->InitCouplesTree("couples",Form("b00000%i/p0%i/1.%i.0.0.cp.root",nbrick,i,i),"RECREATE");
 }
 Int_t Flag = 1;
 cout<<"Start processing nevents: "<<nevents<<endl;  
 //************************STARTING LOOP ON SIMULATION******************  
// while (reader.Next()){
 for (int i = 0; i < nevents; i++){
  if (ievent%1000==0) cout<<"processing event "<<ievent<<" out of "<<nevents<<endl;
  reader.Next();
   for (const EmuBaseTrk& emupoint:emulsionhits){   
     bool savehit = true; //by default I save all hits
//no you don't want to do this//     if (j % 2 == 0) continue;
     trackID = emupoint.GetMCTrackID();

     if (trackID >= 0) motherID = tracks[trackID].GetMotherId();
     xem = emupoint.GetX()* 1E+4 + 62500;
     yem = emupoint.GetY()* 1E+4 + 49500;
     tx = emupoint.GetTX();
     ty = emupoint.GetTY();

     //**************SAVING HIT IN FEDRA BASE-TRACKS****************
             
     ect[nfilmhit-1]->eS->Set(ihit,xem,yem,tx,ty,1,Flag);
     ect[nfilmhit-1]->eS->SetMC(ievent, trackID); //objects used to store MC true information
     ect[nfilmhit-1]->eS->SetAid(motherID, 0); //forcing areaID member to store mother MC track information
     ect[nfilmhit-1]->eS->SetW(ngrains); //need a high weight to do tracking
     ect[nfilmhit-1]->eS->SetVid(pdgcode,0); //forcing viewID[0] member to store pdgcode information
     ect[nfilmhit-1]->Fill();
     ihit++; //hit entry, increasing as the tree is filled        
     
     }//end of loop on emulsion points
    ievent++;
   } //end of loop on tree
  for (int iplate = 0; iplate < nplates; iplate++){
   ect[iplate]->Write();  
   ect[iplate]->Close();  
 }
}
*/
void fromFairShip2Fedra(TString filename){
 const int nplates = 29;
 int nbrick = 1; // to set b00000%i number

 const bool useefficiencymap = false; //use the map instead of the constant value down
 const bool dosmearing = true; //gaussian smearing or not

 if (useefficiencymap){ 
  file = TFile::Open("efficiency_alltracks.root");
  heff = (TH1D*) file->Get("heff");
 }
 const float emuefficiency = 1.0; // i am efficient to everything
 const float angres = 0.003; // 3 milliradians
 const float ngrains = 70; //the same number for all the couples, so they have the same weigth.
 //**********************OPENING INPUT FILE***************************
 TFile * inputfile = TFile::Open(filename.Data());
 //TFile * inputfile = TFile::Open("/eos/user/a/aiuliano/sims_FairShip/sim_charm/pot/uniformonespill_onelayer_ch1_03_03_19/pythia8_Geant4_1000_0.5.root");
 if (!inputfile) return;

 //getting tree and arrays
 TTreeReader reader("cbmsim",inputfile);
 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<BoxPoint> emulsionhits(reader,"BoxPoint");
 
 //TTree* cbmsim = (TTree*)inputfile->Get("cbmsim");
 Float_t tx = 0, ty=0, xem= 0, yem = 0;
 const Int_t nevents = reader.GetTree()->GetEntries();
 int ihit = 0, ievent = 0;
 int nfilmhit = 0;
 float tantheta;
 int trackID = 0, motherID = 0, pdgcode = 0;
 //***********************CREATING FEDRA TREES**************************
 gInterpreter->AddIncludePath("/afs/cern.ch/work/a/aiuliano/public/fedra/include");
 EdbCouplesTree *ect[nplates];
 for (int i = 1; i <= nplates; i++){
  ect[i-1] = new EdbCouplesTree();
  if (i <10) ect[i-1]->InitCouplesTree("couples",Form("b00000%i/p00%i/2.%i.0.0.cp.root",nbrick,i,i),"RECREATE");
  else ect[i-1]->InitCouplesTree("couples",Form("b00000%i/p0%i/2.%i.0.0.cp.root",nbrick,i,i),"RECREATE");
 }
 Int_t Flag = 1;
 cout<<"Start processing nevents: "<<nevents<<endl;  
 //************************STARTING LOOP ON SIMULATION******************  
// while (reader.Next()){
 for (int i = 0; i < nevents; i++){
  if (ievent%1000==0) cout<<"processing event "<<ievent<<" out of "<<nevents<<endl;
  reader.Next();
   for (const BoxPoint& emupoint:emulsionhits){   
     bool savehit = true; //by default I save all hits
//no you don't want to do this//     if (j % 2 == 0) continue;
     pdgcode = emupoint.PdgCode();
     trackID = emupoint.GetTrackID();

     if (trackID >= 0) motherID = tracks[trackID].GetMotherId();
     xem = emupoint.GetX()* 1E+4 + 62500;
     yem = emupoint.GetY()* 1E+4 + 49500;
     tx = emupoint.GetPx()/emupoint.GetPz();
     ty = emupoint.GetPy()/emupoint.GetPz();

     tantheta = pow(pow(tx,2) + pow(ty,2),0.5);
     double charge;        

     if ((TDatabasePDG::Instance()->GetParticle(pdgcode))!=NULL) charge = TDatabasePDG::Instance()->GetParticle(pdgcode)->Charge();
     else charge = 0.;
     nfilmhit = emupoint.GetDetectorID(); //getting number of the film
     //*************EXCLUDE HITS FROM BEING SAVED*******************
     if (nfilmhit > 1000) savehit = false;
     if (tantheta > 1) savehit = false; //we scan from theta 0 to a maximum of 1 rad
     if(charge == 0.) savehit = false; //we do not track neutral particles
     //saving the hits for a plate in the corresponding couples (only one layer saved, the other has ID + 10000)             
 
     if (useefficiencymap){ //efficiency map with angle
     if(!efficiency(tantheta, heff)) savehit = false; //inserting some holes due to emulsion inefficiency           
     }
     //constant value of efficiency
     else if(!efficiency(emuefficiency)) savehit = false;
     if (dosmearing) smearing(tx,ty,angres);
     //**************SAVING HIT IN FEDRA BASE-TRACKS****************
     if (savehit){        
      ect[nfilmhit-1]->eS->Set(ihit,xem,yem,tx,ty,1,Flag);
      ect[nfilmhit-1]->eS->SetMC(ievent, trackID); //objects used to store MC true information
      ect[nfilmhit-1]->eS->SetAid(motherID, 0); //forcing areaID member to store mother MC track information
      ect[nfilmhit-1]->eS->SetVid(pdgcode,0); //forcing viewID[0] member to store pdgcode information
      ect[nfilmhit-1]->eS->SetW(ngrains); //need a high weight to do tracking
      ect[nfilmhit-1]->Fill();
      ihit++; //hit entry, increasing as the tree is filled        
      }
     }//end of loop on emulsion points
    ievent++;
   } //end of loop on tree
  for (int iplate = 0; iplate < nplates; iplate++){
   ect[iplate]->Write();  
   ect[iplate]->Close();  
 }
}

void smearing (Float_t &TX, Float_t &TY, const float angres){
 float deltaTX = grandom->Gaus(0,angres); //angular resolution, adding a gaussian offset to TX and TY
 float deltaTY = grandom->Gaus(0,angres);
 //cout<<TX<<endl;
 TX = TX + deltaTX;
 TY = TY + deltaTY;
}

bool efficiency(const float emuefficiency){ //for now, just a constant, to be replaced with an efficiency map (probably with the angle)

 float prob = grandom->Uniform(0,1);
 if (prob < emuefficiency) return true; //efficiency larger than probability, we take the event
 else return false;
}

bool efficiency(const float tantheta, TH1D * emuefficiency){ //for now, just a constant, to be replaced with an efficiency map (probably with the angle)

 float prob = grandom->Uniform(0,1);
 int ibin = emuefficiency->FindBin(tantheta); 
 const float efficiency = emuefficiency->GetBinContent(ibin);
 if (prob < efficiency) return true; //efficiency larger than probability, we take the event
 else return false;
}
/*
script ispirato a void ReadGEMdata da parte di Annarita
*/
