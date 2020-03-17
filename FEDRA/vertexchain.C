void vertexchain(){
 TString prepath("/eos/user/a/aiuliano/public/sims_fedra/CH3_pot_05_03_2020/b000001/reconstruction_output/");
 TChain vtxchain("vtx");

 //adding files for various quarters
 vtxchain.Add((prepath+TString("firstquarter/vertextree_firstquarter.root")).Data());
 vtxchain.Add((prepath+TString("secondquarter/vertextree_secondquarter.root")).Data());
 vtxchain.Add((prepath+TString("thirdquarter/vertextree_thirdquarter.root")).Data());
 vtxchain.Add((prepath+TString("fourthquarter/vertextree_fourthquarter.root")).Data());
 
 TCut selection("flag!=2 && flag!=5");
 
 //drawing vertex distributions
 vtxchain.Draw("vz>>hz",selection);
 vtxchain.Draw("vy:vx>>hxy",selection);
 vtxchain.Draw("n>>hn",selection);
 //getting histograms from memory
 TH1D *hz = (TH1D*) gDirectory->FindObject("hz");
 TH1I *hn = (TH1I*) gDirectory->FindObject("hn");
 TH2D *hxy = (TH2D*) gDirectory->FindObject("hxy");
 //setting title and drawing them
 gStyle->SetStatX(0.5);
 gStyle->SetStatY(0.9);  

 TCanvas *cz = new TCanvas();
 hz->SetTitle("vz distribution;z[#mum]");
 hz->Draw();

 TCanvas *cxy = new TCanvas();
 hxy->SetTitle("vxy distribution;x[#mum];y[#mum]");
 hxy->Draw("COLZ"); 

 TCanvas *cn = new TCanvas();
 hn->SetTitle("ntracks;Ntrk");
 hn->Draw();
}
