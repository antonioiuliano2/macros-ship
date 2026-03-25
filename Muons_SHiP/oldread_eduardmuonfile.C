//OLD VERSION; WITH TTREE DRAW! WEIGHTS ARE NOT CORRECT; DO NOT USE!
void read_eduardmuonfile(){
    //first, accessing the file from EOS
    TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/snd_ship_ecn3/muons_3_planes.root");
    TTree *cbmsim = (TTree*) inputfile->Get("cbmsim");
    cout<<"Reading a tree with "<<cbmsim->GetEntries()<<endl;

    TFile *outputfile = new TFile("muon_snd_ship_ecn3_plots.root","RECREATE");

    TCanvas *cxy = new TCanvas("cxy","xy distribution",800,800);
    cbmsim->Draw("sco_1Point.fY:sco_1Point.fX>>hxy(140,-600,800,140,-600,800)","MCTrack[0].fW","COLZ");
    TH2D *hxy = (TH2D*) gDirectory->Get("hxy");
    hxy->SetTitle("xy in scoring plane 1;x[cm];y[cm]");
    cxy->Write();
    
    TCanvas *cxy_zoomed = new TCanvas("cxy_zoomed","xy distribution zoomed",800,800);
    cbmsim->Draw("sco_1Point.fY:sco_1Point.fX>>hxy_zoomed(400,-200,200,400,-200,200)","MCTrack[0].fW","COLZ");
    TH2D *hxy_zoomed = (TH2D*) gDirectory->Get("hxy_zoomed");
    hxy_zoomed->SetTitle("xy zoomed in scoring plane 1;x[cm];y[cm]");
    cxy_zoomed->Write();

    TCanvas *cxy_zoomed_text = new TCanvas("cxy_zoomed_text","xy distribution zoomed_withtext",800,800);
    cbmsim->Draw("sco_1Point.fY:sco_1Point.fX>>hxy_zoomed_text(40,-200,200,40,-200,200)","MCTrack[0].fW","COLZ && TEXT");
    TH2D *hxy_zoomed_text = (TH2D*) gDirectory->Get("hxy_zoomed_text");
    hxy_zoomed_text->SetTitle("xy zoomed in scoring plane 1;x[cm];y[cm]");   
    cxy_zoomed_text->Write(); 

    TCanvas *cP = new TCanvas("cP","P distribution");
    cbmsim->Draw("TMath::Sqrt(sco_1Point.fPx*sco_1Point.fPx+sco_1Point.fPy*sco_1Point.fPy+sco_1Point.fPz*sco_1Point.fPz)>>hP","MCTrack[0].fW");
    TH1D *hP = (TH1D*) gDirectory->Get("hP");
    hP->SetTitle("momentum in scoring plane 1;P[GeV/c]");
    cP->Write();

    TCanvas *cpdgcode = new TCanvas("cpdgcode","PdgCode");
    cbmsim->Draw("sco_1Point.fPdgCode>>hPdgCode","MCTrack[0].fW");
    TH1D *hPdgCode = (TH1D*) gDirectory->Get("hPdgCode");
    hPdgCode->SetTitle("PdgCode in scoring plane 1;PdgCode");
    cpdgcode->Write();
    outputfile->Close();
    
}