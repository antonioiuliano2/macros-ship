//simply create spill infotree from csv file
void spillinfotree(){
  TFile *outputfile = new TFile("charm_spills.root","recreate");
  TTree *spilltree = new TTree("spilltree","Scaler information from SHiP-Charm July 2018 testbeam");

  TString columns = TString("runcode/I:name/I:pot[10]/I:spilltime[10]/I");
  spilltree->ReadFile("/eos/experiment/ship/data/charmxsec/bookkeeping/spillinfo.csv",columns.Data());

  outputfile->Write();
  outputfile->Close();
}
