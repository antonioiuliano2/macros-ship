//copy just some entries in a new tree

void copyentry_cbmsim(){
 TFile * origfile = TFile::Open("pythia8_Geant4_1000_0.1_dig.root","READ");
 TTree *origtree = (TTree*) origfile->Get("cbmsim"); 
 TList * origlist = (TList*) origfile->Get("BranchList"); //we need this list for the display too
 TList * origtimelist = (TList*) origfile->Get("TimeBasedBranchList"); //we need this list for the display too

 TFolder *cbmfolder = (TFolder*) origfile->Get("cbmroot");
 FairFileHeader *fileheader = (FairFileHeader*) origfile->Get("FileHeader");

 TFile * newfile = new TFile("selectedtree.root","RECREATE"); 
 TTree * newtree = origtree->CloneTree(0);


 origtree->GetEntry(29584);
 newtree->Fill();

 //saving output
 origlist->Write("BranchList",1);
 origtimelist->Write("TimeBasedBranchList",1);

 cbmfolder->Write("cbmroot",1);
 fileheader->Write("FileHeader",1);

 newtree->Write();

 origfile->Close();
 newfile->Close();
}
