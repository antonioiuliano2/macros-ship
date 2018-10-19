void writeweights(){
TFile *file = TFile::Open("/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_10.0_withCharmandBeauty50000_mu.root");
TTree *tree = (TTree*) file->Get("cbmsim");

tree->Scan("MCTrack.fW","MCTrack.fPdgCode == 13 || MCTrack.fPdgCode == -13");

}
