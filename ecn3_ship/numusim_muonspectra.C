//reading numu simulation, AFTER they have cut within target acceptance
void numusim_muonspectra(){
 
 TString prefix("root:://eosuser.cern.ch/");
 TString simpath_mu_bar("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_04_14_numu_bar_CHARMCCDIS_spectro_1_2T/");
 TFile *inputfile = TFile::Open((prefix+simpath_mu_bar+TString("inECC_ship.conical.Genie-TGeant4.root")).Data());
 TTree *simtree = (TTree*) inputfile->Get("cbmsim");
 
 TTreeReader reader(simtree);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");

 TH1D * hmuE = new TH1D("hmuE","Muon energy;E[GeV]",400,0,400);
 

 TFile *outputfile = new TFile("muoninfo.root","RECREATE");
 TNtuple * muoninfo = new TNtuple("muoninfo","Information of primary muons from numu CCDIS","ievent:px:py:pz:x:y:z:weight");
 

 const int nentries = reader.GetEntries(true);
 cout<<"Number of events"<<nentries<<endl;

 for(int ientry = 0;ientry<nentries;ientry++){   
  reader.SetEntry(ientry);

  if (TMath::Abs(tracks[1].GetPdgCode())==13){ //saving initial muon information
    Double_t weight = tracks[1].GetWeight();

    hmuE->Fill(tracks[1].GetEnergy(),weight);
    //storing muon information
    muoninfo->Fill(ientry, tracks[1].GetPx(), tracks[1].GetPy(), tracks[1].GetPz(), tracks[1].GetStartX(), tracks[1].GetStartY(), tracks[1].GetStartZ(), weight);
   } 
 } 
 TCanvas *c = new TCanvas();
 hmuE->Scale(1./hmuE->Integral());
 hmuE->Draw("histo");

 //storing ntuple and closing outputfile
 outputfile->cd();
 muoninfo->Write();
 outputfile->Close();

 


}
