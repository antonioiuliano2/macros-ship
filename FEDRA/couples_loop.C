//basic script to show how to loop in couples (written in 23 May 2019)
void couples_loop(){

 TH2D *htxty = new TH2D("htxty","Histogram for test of reading the tree",300,-1.5,1.5,300,-1.5,1.5);

 TFile *couplesfile = TFile::Open("/ship/CHARM2018/CH3-R3/b000001/p001/1.1.0.0.cp.root");
 TTree *couples = (TTree*) couplesfile->Get("couples"); //just for knowing the entries
 const int nentries = couples->GetEntries();

 EdbCouplesTree *mytree = new EdbCouplesTree();
 mytree->InitCouplesTree("couples","/ship/CHARM2018/CH3-R3/b000001/p001/1.1.0.0.cp.root","READ");

 for (int i = 0; i < nentries; i++){
  mytree->GetEntry(i);

  EdbSegP *myseg = mytree->eS;
  htxty->Fill(myseg->TX(),myseg->TY());
 }
 htxty->Draw("COLZ");
}
