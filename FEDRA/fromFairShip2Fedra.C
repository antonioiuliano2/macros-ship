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
 
 fromFairShip2Fedra("ship.conical.Pythia8CharmOnly-TGeant4_dig.root");
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

void set_default(TEnv &cenv){ //setting default parameters, if not presents from file
 cenv.SetValue("FairShip2Fedra.nbrick",1);//to set b00000%i number
 cenv.SetValue("FairShip2Fedra.nplates",29);
 cenv.SetValue("FairShip2Fedra.nevents",10000); // number of events to be passed to FEDRA
 cenv.SetValue("FairShip2Fedra.neventsxspill",1000); // number of events to be passed to FEDRA
 cenv.SetValue("FairShip2Fedra.useefficiencymap",0);
 cenv.SetValue("FairShip2Fedra.emuefficiency",0.85); //only if useefficiency map is set to false
 cenv.SetValue("FairShip2Fedra.dosmearing",1);
 cenv.SetValue("FairShip2Fedra.maxtheta",1); //angular max of scanning
 cenv.SetValue("FairShip2Fedra.minkinenergy",0.1); //do not pass particles beyond this value, track ID would be -2
 cenv.SetValue("FairShip2Fedra.ngrains",70); // to set weight
 cenv.SetValue("FairShip2Fedra.angres",0.003);//used for smearing, if dosmearing = true

}

//#include "/home/utente/fedra/include/EdbCouplesTree.h"
using namespace TMath;
TRandom *grandom = new TRandom3(); //creating every time a TRandom3 is a bad idea
TFile *file = NULL;
TH1D *heff = NULL ; //efficiency at different angles
void fromFairShip2Fedra(TString filename){

 TEnv cenv("FairShip2Fedra");
 set_default(cenv);
 cenv.ReadFile("FairShip2Fedra.rootrc" ,kEnvLocal);
 //getting options from file
 const Int_t nevents = cenv.GetValue("FairShip2Fedra.nevents",10000);
 const int nplates = cenv.GetValue("FairShip2Fedra.nplates",29);
 int nbrick = cenv.GetValue("FairShip2Fedra.nbrick",1); // to set b00000%i number

 float angres = cenv.GetValue("FairShip2Fedra.angres",0.003); //Used cases: 3, 5milliradians. Constant value overwritten if useresfunction=true
 float minkinE = cenv.GetValue("FairShip2Fedra.minkinenergy",0.1);
 float maxtheta = cenv.GetValue("FairShip2Fedra.maxtheta",1);

 const float ngrains = cenv.GetValue("FairShip2Fedra.ngrains",70) ; //the same number for all the couples, so they have the same weigth.
 const float emuefficiency = cenv.GetValue("FairShip2Fedra.emuefficiency",0.85); // flat value
 
 const bool useefficiencymap = cenv.GetValue("FairShip2Fedra.useefficiencymap",0); //use the map instead of the constant value down
 const bool dosmearing = cenv.GetValue("FairShip2Fedra.dosmearing",1); //gaussian smearing or not
 const bool useresfunction = false; //use resfunction from operadata instead of constant value
 
 cout<<"Starting conversion with efficiency "<<emuefficiency<<" maxtheta "<<maxtheta<<" and min kin E "<<minkinE<<endl;
 //if not performed digitization
 const bool donedigi = false;
 int neventsxspill = cenv.GetValue("FairShip2Fedra.neventsxspill",1000);
 int ntotspills = nevents/neventsxspill;
 float spilldy = 10./ntotspills;
 cout<<"Generating "<<ntotspills<<" spills with dy "<<spilldy<<endl;
 int nspill = 0;
 float pottime = 0.;
 float targetmoverspeed = 2.6; //speed of Target Mover

 if (useefficiencymap){ 
  file = TFile::Open("efficiency_alltracks.root");
  heff = (TH1D*) file->Get("heff");
 }
 TF1 resfunction;
 if (useresfunction) resfunction = angularresolution();
 // **********************OPENING INPUT FILE***************************
 TFile * inputfile = TFile::Open(filename.Data());

 if (!inputfile) return;

 //getting tree and arrays
 TTreeReader reader("cbmsim",inputfile);
 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
// TTreeReaderArray<BoxPoint> emulsionhits(reader,"EmuBaseTrks");
 TTreeReaderArray<BoxPoint> emulsionhits(reader,"BoxPoint");
 
 Float_t tx = 0, ty=0, xem= 0, yem = 0;
 //const Int_t nevents = reader.GetTree()->GetEntries();
 int ihit = 0, ievent = 0;
 int nfilmhit = 0;
 float tantheta, momentum;
 int trackID = 0, motherID = 0, pdgcode = 0;
 // ***********************CREATING FEDRA TREES**************************
 gInterpreter->AddIncludePath("/afs/cern.ch/work/a/aiuliano/public/fedra/include");
 EdbCouplesTree *ect[nplates];
 for (int i = 1; i <= nplates; i++){
  ect[i-1] = new EdbCouplesTree();
  if (i <10) ect[i-1]->InitCouplesTree("couples",Form("b00000%i/p00%i/%i.%i.0.0.cp.root",nbrick,i,nbrick,i),"RECREATE");
  else ect[i-1]->InitCouplesTree("couples",Form("b00000%i/p0%i/%i.%i.0.0.cp.root",nbrick,i,nbrick,i),"RECREATE");
 }
 Int_t Flag = 1;
 cout<<"Start processing nevents: "<<nevents<<endl;  
 // ************************STARTING LOOP ON SIMULATION******************  
// while (reader.Next()){
 for (int i = 0; i < nevents; i++){
  if (i%1000==0) cout<<"processing event "<<i<<" out of "<<nevents<<endl;
  reader.Next();
  pottime = gRandom->Uniform()*4.8;
  nspill = i/neventsxspill;
  for (const BoxPoint& emupoint:emulsionhits){   
     bool savehit = true; //by default I save all hits
//no you don't want to do this//     if (j % 2 == 0) continue;
     momentum = TMath::Sqrt(pow(emupoint.GetPx(),2) + pow(emupoint.GetPy(),2) + pow(emupoint.GetPz(),2));
     pdgcode = emupoint.PdgCode();
     trackID = emupoint.GetTrackID();
     bool emubasetrackformat = true;
     
     if (trackID >= 0) motherID = tracks[trackID].GetMotherId();
     else motherID = -2; //hope I do not see them

     if (!donedigi){    

      xem = emupoint.GetX() -12.5/2. + pottime * targetmoverspeed;
      yem = emupoint.GetY() - 9.9/2. + nspill * spilldy + 0.5;
      } 
    
      else{

       xem = emupoint.GetX();
       yem = emupoint.GetY();

      }

      xem = xem* 1E+4 + 62500;
      yem = yem* 1E+4 + 49500;         
     
      tx = emupoint.GetPx()/emupoint.GetPz();
      ty = emupoint.GetPy()/emupoint.GetPz();  
      tantheta = pow(pow(tx,2) + pow(ty,2),0.5);

     double charge,mass;        

     if ((TDatabasePDG::Instance()->GetParticle(pdgcode))!=NULL){ 
      charge = TDatabasePDG::Instance()->GetParticle(pdgcode)->Charge();
      mass = TDatabasePDG::Instance()->GetParticle(pdgcode)->Mass();
      }
     else{ 
      charge = 0.;
      mass = 0.;
      }
     nfilmhit = emupoint.GetDetectorID(); //getting number of the film
     double kinenergy = TMath::Sqrt(pow(mass,2)+pow(momentum,2)) - mass;
     // *************EXCLUDE HITS FROM BEING SAVED*******************
     if (nfilmhit > 1000) savehit = false;
     if (tantheta > TMath::Tan(maxtheta)) savehit = false; //we scan from theta 0 to a maximum of 1 rad
     if(charge == 0.) savehit = false; //we do not track neutral particles
     if(kinenergy < minkinE) savehit = false; //particles with too low kin energy 
     //saving the hits for a plate in the corresponding couples (only one layer saved, the other has ID + 10000)             
 
     if (useefficiencymap){ //efficiency map with angle
     if(!efficiency(tantheta, heff)) savehit = false; //inserting some holes due to emulsion inefficiency           
     }
     //constant value of efficiency
     else if(!efficiency(emuefficiency)) savehit = false;
     if (dosmearing){ 
       if (useresfunction) angres = resfunction.Eval(TMath::ATan(tantheta));
      smearing(tx,ty,angres);
     }
     // **************SAVING HIT IN FEDRA BASE-TRACKS****************
     if (savehit){        
      ect[nfilmhit-1]->eS->Set(ihit,xem,yem,tx,ty,1,Flag);
      ect[nfilmhit-1]->eS->SetMC(ievent, trackID); //objects used to store MC true information
      ect[nfilmhit-1]->eS->SetP(momentum); //storing true momentum of the hit
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
 cout<<"end of script, saving rootrc wih used parameters"<<endl;
 cenv.WriteFile("FairShip2Fedra.save.rootrc");
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
