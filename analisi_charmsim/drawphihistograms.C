void drawphihistograms(){
    //opening files
    const double phimin = 2.2; //minimum phi to accept selection
    TFile *signalfile = TFile::Open("$HOME/cernbox/Synched/Archivio_cronologico/Marzo_2021/testdeltaphi_sig.root");
    TCanvas *csignal = (TCanvas*) signalfile->Get("c1_n12");
    TFile *bkgfile = TFile::Open("$HOME/cernbox/Synched/Archivio_cronologico/Marzo_2021/testdeltaphi_bkg.root");
    TCanvas *cbkg = (TCanvas*) bkgfile->Get("c1_n12");
    //retrieving histograms
    csignal->Draw();
    TH1D* hphi = (TH1D*) csignal->GetPrimitive("hphi");
    TH1D* hphisignal = new TH1D(*hphi); //I need to copy, otherwise I will lose itwhen I get the other histo
    hphisignal->SetName("hphisignal");
    hphisignal->SetTitle("Phi distribution for CCDIS tau neutrino signal");
    cbkg->Draw();
    hphi = (TH1D*) cbkg->GetPrimitive("hphi");
    TH1D* hphibkg = new TH1D(*hphi); //I need to copy, otherwise I will lose itwhen I get the other histo
    hphibkg->SetName("hphibackground");
    hphibkg->SetTitle("Phi distribution for CharmCCDIS muon neutrino background");
    
    //setting range properly, scaling histograms

    hphibkg->GetXaxis()->SetRangeUser(0,3.2);
    hphibkg->Scale(1./hphibkg->Integral());
    hphisignal->GetXaxis()->SetRangeUser(0,3.2);
    hphisignal->Scale(1./hphisignal->Integral());
    TCanvas *c1 = new TCanvas();

    hphisignal->Draw("histo");
    hphibkg->SetLineColor(kRed);
    hphibkg->Draw("SAMES&&histo");
    c1->BuildLegend();
    TLine *l = new TLine(phimin, 0, phimin, hphisignal->GetMaximum());
    l->SetLineColor(kRed);
    l->Draw("SAME");

}