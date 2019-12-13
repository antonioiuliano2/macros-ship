void plot_thomasds(){
 
 TFile *thomasfile = TFile::Open("/eos/experiment/ship/data/Mbias/background-prod-2018/Cascade-run0-19-parp16-MSTP82-1-MSEL4-10Bpot.root");
 TTree *cascadetree = (TTree*) thomasfile->Get("pythia6");

 TH1D *hE_dsminus_prim = new TH1D("hE_dsminus_prim","Histogram of energy dsminus;E[GeV]",400,0,400);
 TH1D *hE_dsplus_prim = new TH1D("hE_dsplus_prim","Histogram of energy dsplus;E[GeV]",400,0,400);

 TH1D *hE_dsminus_other = new TH1D("hE_dsminus_other","Histogram of energy dsminus;E[GeV]",400,0,400);
 TH1D *hE_dsplus_other = new TH1D("hE_dsplus_other","Histogram of energy dsplus;E[GeV]",400,0,400);

 TCanvas *cds_comparison_primary = new TCanvas();
 cascadetree->Draw("E>>hE_dsplus_prim","id==431&&mE>398");
 cascadetree->Draw("E>>hE_dsminus_prim","id==-431&&mE>398");

 hE_dsplus_prim->Draw();

 hE_dsplus_prim->SetLineColor(kRed);
 hE_dsminus_prim->Draw("SAMES");

 cds_comparison_primary->BuildLegend();
 hE_dsplus_prim->SetTitle("Mother energy larger than 398 GeV");
 cds_comparison_primary->Draw();

 TCanvas *cds_comparison_other = new TCanvas();
 cascadetree->Draw("E>>hE_dsplus_other","id==431&&mE<398");
 cascadetree->Draw("E>>hE_dsminus_other","id==-431&&mE<398");

 hE_dsplus_other->Draw();

 hE_dsplus_other->SetLineColor(kRed);
 hE_dsminus_other->Draw("SAMES");

 hE_dsplus_other->SetTitle("Mother energy less than 398 GeV");
 cds_comparison_other->BuildLegend();
 cds_comparison_other->Draw();

}
