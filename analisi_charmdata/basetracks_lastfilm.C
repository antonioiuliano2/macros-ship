//a small script to remember which script I used for the beam plot in my PhD thesis (6 February 2021), CH1-R2 and CH1-R6
void basetracks_lastfilm(){
 TFile *cpfile = TFile::Open("1.29.0.0.cp.root");
 TTree *couples = (TTree*) cpfile->Get("couples");

 couples->Draw("s.eY*1e-3:s.eX*1e-3>>hxy(125,0,125,100,0,100)","s.Theta()<0.1","COLZ");
 TH2D *hxy = (TH2D*) gDirectory->FindObject("hxy");

 hxy->Draw("COLZ");
 hxy->SetTitle(";x[mm];y[mm]");
 gStyle->SetOptStat("n");
}