//where high number of electrons are generated?
#include <assert.h>     /* assert */
void muonrpc_loop(){
 TFile *inputfile = TFile::Open("ship.conical.MuonBack-TGeant4.root");
 TTreeReader reader("cbmsim",inputfile); //reading file loaded before executing the script

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");
 
 //histograms to be instantiated
 TFile *histofile = new TFile("origin_hitsrpc.root","RECREATE");
 TH1I *hprocessID = new TH1I("hprocessID","IDs of the prcoess generating the electrons/positrons registered in RPC station 0",50,0,50);
 TH1I *hmotherpdg = new TH1I("hmotherpdg","PdgCode of the mother of the electrons/positrons registered in RPC station 0",50,-25,25);

 TH1F *hstartz = new TH1F("hstartz", "production position of the electron/positrons registered in RPC station 0",4000,-6000,-2000);
 TH2F *hstartxy = new TH2F("hstartxy", "production position of the electron/positrons registered in RPC station 0", 100,-500,500,100, -500, 500);
 //loop on the events
 const int nevents = (reader.GetTree())->GetEntries();
 cout<<"number of entries: "<<nevents<<endl;
 int ievent = 0;
 while (reader.Next()){
  if (ievent % 100000 == 0) cout<<"now processing event "<<ievent<<" over "<<nevents<<endl;
  // if (ievent >= 100000) break;
  int irpc = 0; //to keep track of rpcindex
  //loop on rpc points
  for (const ShipRpcPoint &rpcpoint: rpcpoints){ 
   int trackID = rpcpoint.GetTrackID();
   int hitpdg = rpcpoint.PdgCode();
   int detID = rpcpoint.GetDetectorID();
   int nplane = detID - 10000;

   assert(trackID>0); //we can access track information for this hit
   int procID = tracks[trackID].GetProcID();
   //xyz start position of the track
   Float_t startx = tracks[trackID].GetStartX();
   Float_t starty = tracks[trackID].GetStartY();
   Float_t startz = tracks[trackID].GetStartZ();
   Float_t momentum = tracks[trackID].GetP();
   Float_t weight = tracks[trackID].GetWeight();   
   //we can also access information about the mother
   int motherID = tracks[trackID].GetMotherId();
   Float_t mothermomentum = tracks[motherID].GetP();
   int motherpdgcode = tracks[motherID].GetPdgCode();
   if ((nplane == 0) && (TMath::Abs(hitpdg) == 11)){ //CUT for histogram filling (electrons in RPC number 0
    cout<<"rpcpos: "<<rpcpoint.GetX()<<" "<<rpcpoint.GetY()<<" "<<rpcpoint.GetZ()<<endl;
    hprocessID->Fill(procID,weight);
    hstartxy->Fill(startx, starty,weight);
    hstartz->Fill(startz,weight);
    hmotherpdg->Fill(motherpdgcode,weight);

   } //end fill 
   
   irpc++;
  } //end loop on hits
 ievent++;
 } //end of the loop over the vents
 inputfile->Close();
 //**********************drawing of the histograms
 TCanvas *cstartpos = new TCanvas();
 cstartpos->Divide(1,2);
 cstartpos->cd(1);
 hstartz->GetXaxis()->SetTitle("z[cm]");
 hstartz->Draw("histo");
 hstartz->Write();
 cstartpos->cd(2);
 hstartxy->GetXaxis()->SetTitle("x[cm]");
 hstartxy->GetYaxis()->SetTitle("y[cm]");
 hstartxy->Draw("COLZ");
 hstartxy->Write();

 TCanvas *cmotherpdg = new TCanvas();
 hmotherpdg->GetXaxis()->SetTitle("mumPDGcode");
 hmotherpdg->Draw("histo");
 hmotherpdg->Write();
 TCanvas *cprocID = new TCanvas();
 hprocessID->GetXaxis()->SetTitle("procID");
 hprocessID->Draw("histo");
 hprocessID->Write();

 histofile->Close();
}
