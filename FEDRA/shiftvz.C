void shiftvz(){
 TFile *vertexfile = TFile::Open("vertextree.root");

 TTree *vtx = (TTree*) vertexfile->Get("vtx");

 float vz, newvz;
 const float vzshift = 36820;
 vtx->SetBranchAddress("vz",&vz);

 TFile *newvertexfile = new TFile("vertextree_correctvz.root","RECREATE");
 vtx->SetBranchStatus("vz",0); //do not copy it
 TTree *newvtx = vtx->CloneTree(0);
 newvtx->Branch("vz",&newvz,"vz/F");
 const int nvertices = vtx->GetEntries();

 vtx->SetBranchStatus("vz",1);//now I need to read it
 for (int ivtx = 0; ivtx < nvertices; ivtx++){
  vtx->GetEntry(ivtx);
  newvz = vz + vzshift;
  newvtx->Fill();
 }
 newvertexfile->cd();
 newvtx->Write();
 newvertexfile->Close();
}
