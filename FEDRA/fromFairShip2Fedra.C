//tool to load hit from emulsion to FEDRA (13 aprile 2018) //updated to work with the new TTreeReader structure (6 March)

#include <stdio.h>
#include <TROOT.h>
#include "TRandom.h"

//#include "/home/utente/fedra/include/EdbCouplesTree.h"
using namespace TMath;

void smearing (Float_t &TX, Float_t &TY, const float angres);
bool efficiency(const float emuefficiency);
void fromFairShip2Fedra(int nplate);

void fromFairShip2Fedra(){
	const int nplates = 29;
	for (int iplate = 1; iplate <= 29; iplate++){
		fromFairShip2Fedra(iplate);
	}
}

void fromFairShip2Fedra(int nplate){
 const float emuefficiency = 0.9;
 const float angres = 0.003; // 3 milliradians
 //TFile * inputfile = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/sim_charm/pot/Charm1_uniform_19_11_2018/pythia8_evtgen_Geant4_1000_0.001.root");
 TFile * inputfile = TFile::Open("/eos/user/a/aiuliano/sims_FairShip/sim_charm/pot/uniformonespill_onelayer_ch1_03_03_19/pythia8_Geant4_1000_0.5.root");
 if (!inputfile) return;

 //getting tree and arrays
 TTreeReader reader("cbmsim",inputfile);
 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<BoxPoint> emulsionhits(reader,"BoxPoint");
 
 //TTree* cbmsim = (TTree*)inputfile->Get("cbmsim");
 Float_t tx = 0, ty=0, xem= 0, yem = 0;
 const Int_t nevents = reader.GetTree()->GetEntries();
 int ihit = 0, ievent = 0;
 int trackID = 0, motherID = 0;
 gInterpreter->AddIncludePath("/afs/cern.ch/work/a/aiuliano/public/fedra/include");
 EdbCouplesTree ect;
 if (nplate <10) ect.InitCouplesTree("couples",Form("b000001/p00%i/1.%i.0.0.cp.root",nplate,nplate),"RECREATE");
 else ect.InitCouplesTree("couples",Form("b000001/p0%i/1.%i.0.0.cp.root",nplate,nplate),"RECREATE");
 Int_t Flag = 1;
   
   
 while (reader.Next()){
   for (const BoxPoint& emupoint:emulsionhits){
     bool savehit = true; //by default I save all hits
//no you don't want to do this//     if (j % 2 == 0) continue;
     trackID = emupoint.GetTrackID();
     motherID = emupoint.GetTrackID();
     xem = emupoint.GetX()* 1E+4 + 62500;
     yem = emupoint.GetY()* 1E+4 + 49500;
     tx = emupoint.GetPx()/emupoint.GetPz();
     ty = emupoint.GetPy()/emupoint.GetPz();
     //saving the hits for a plate in the corresponding couples (only one layer saved, the other has ID + 10000)    
     if (emupoint.GetDetectorID() == nplate){
       //Inserting real data effects in the simulation
       if(!efficiency(emuefficiency)) savehit = false; //inserting some holes due to emulsion inefficiency
       smearing(tx,ty,angres);
       //Saving the hit in FEDRA
       if (savehit){        
        ect.eS->Set(ihit,xem,yem,tx,ty,1,Flag);
        ect.eS->SetMC(ievent, trackID) //objects used to store MC true information
        ect.eS->SetAid(motherID, 0); //forcing areaID member to store mother MC track information
        ect.eS->SetW(70.); //need a high weight to do tracking
     //if(do_invert) ect.eS->Transform(&aff_invert);
        ect.Fill();
        ihit++; //hit entry, increasing as the tree is filled
        }
       }
     }//end of loop on emulsion points
    ievent++;
   } //end of loop on tree

  ect.Close();  
}

void smearing (Float_t &TX, Float_t &TY, const float angres){
 TRandom *grandom = new TRandom3();
 float deltaTX = grandom->Gaus(0,angres); //angular resolution, adding a gaussian offset to TX and TY
 float deltaTY = grandom->Gaus(0,angres);
 TX = TX + deltaTX;
 TY = TY + deltaTY;
 delete grandom;
}

bool efficiency(const float emuefficiency){ //for now, just a constant, to be replaced with an efficiency map (probably with the angle)

 TRandom3 *grandom = new TRandom3();
 float prob = grandom->Uniform(0,1);
 delete grandom;
 if (prob > emuefficiency) return true;
 else return false;
}
/*
script ispirato a void ReadGEMdata da parte di Annarita
*/
