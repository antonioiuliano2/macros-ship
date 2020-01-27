using namespace ROOT;
void nsegments_firstplate(){

  //getting tree and setting branches
  TFile *showerfile = TFile::Open("/home/antonio/cernbox/Synched/Charmsim_Showreco/Showreco_ds_23_01/Shower_60000_0.root");
  TTree *showertree = (TTree*) showerfile->Get("treebranch");

  int sizeb;
  const int maxsize = 1000;
  float zb[maxsize];

  for (int i = 0; i< maxsize;i++) zb[i] = 100000.;

  showertree->SetBranchAddress("sizeb",&sizeb);
  showertree->SetBranchAddress("zb",&zb);

  //histogram to check shower molteplicity at first plate
  TH1I *hnseg_min = new TH1I("hnseg_min","Shower molteplicity at first plate;nsegments",100,0,100);
  //loop within shower
 
  const int nshowers = showertree->GetEntries();
   cout<<"Looping over "<<nshowers<<" showers "<<endl;
  for (int ishower = 0; ishower<nshowers; ishower++){
     showertree->GetEntry(ishower);
     RVec<float> zsegments(zb, maxsize); //converting to RVec for easier operations

     //find minimum;
     float zmin = Min(zsegments);
     //find how many segments are at minimum
     hnseg_min->Fill(zsegments[zsegments==zmin].size());

  }
 hnseg_min->Draw();
}