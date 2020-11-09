void checkentries(TString filename){
 gErrorIgnoreLevel = kWarning; //avoid printing ROOT standard messages
 TFile *myfile = TFile::Open(filename.Data());
 TTree *tree = (TTree*) myfile->Get("cbmsim");
 
 cout<<tree->GetEntries();
}
