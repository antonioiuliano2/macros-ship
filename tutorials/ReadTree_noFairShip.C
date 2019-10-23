#define ReadTree_noFairShip_cxx
#include "ReadTree_noFairShip.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void ReadTree_noFairShip::Loop()
{
//   In a ROOT session, you can do:
//      root> .L ReadTree_noFairShip.C
//      root> ReadTree_noFairShip t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   //*************HISTOGRAMS AND VARIABLES************
   TH2D *hbeamxy = new TH2D("hbeamxy", "Simulated profile of the neutrino beam", 80, -40 ,40, 80, -40, 40);
   TH2D *hxy = new TH2D("hxy", "2D distribution for first Rpc plane",200,-100,100,150,-50,100);

   //*************START LOOP*****************
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      //*******loop on tracks****
      Double32_t neutrinostartx = MCTrack_fStartX[0]; //(itrk 0 is the first neutrino)
      Double32_t neutrinostarty = MCTrack_fStartY[0];
        
      hbeamxy->Fill(neutrinostartx, neutrinostarty);
      for(int itrk = 0; itrk < MCTrack_;itrk++){
	//do whatever you need to do with the tracks
      }
      //******loop on RPC hits*****
      for(int ihit = 0; ihit < ShipRpcPoint_; ihit++){
        Double32_t rpcx = ShipRpcPoint_fX[ihit];
        Double32_t rpcy = ShipRpcPoint_fY[ihit]; 

        if (ShipRpcPoint_fDetectorID[ihit] == 10000) hxy->Fill(rpcx,rpcy); //detectorIDs for ShipRpcs: 10000, 10001... 
      }
      // if (Cut(ientry) < 0) continue;
   }
   //**************END LOOP******************
   //drawing histograms
   TCanvas * c1 = new TCanvas();
   hbeamxy->Draw("COLZ");
   hbeamxy->GetXaxis()->SetName("x[cm]");
   hbeamxy->GetYaxis()->SetName("y[cm]");
   TCanvas * c2 = new TCanvas();
   hxy->Draw("COLZ");
   hxy->GetXaxis()->SetName("x[cm]");
   hxy->GetYaxis()->SetName("y[cm]");
}
