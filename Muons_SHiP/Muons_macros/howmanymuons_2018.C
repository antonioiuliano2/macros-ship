void howmanymuons_2018(){
Int_t nmuons = 0;
for (int i = 0; i < 67; i++){
  if (i == 53) continue; //broken file
  TFile *file = TFile::Open(Form("/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_10.0_withCharmandBeauty%i_mu.root",i*1000));
  TTree *tree = (TTree*) file->Get("cbmsim");
  cout<<"In file: "<<i*1000<<" there are: "<<tree->GetEntries()<<" entries "<<endl;
  nmuons += tree->GetEntries();
 }
cout<<"Number of total entries: "<<nmuons<<endl;
}
