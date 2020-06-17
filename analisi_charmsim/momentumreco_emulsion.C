//study momentum reconstruction in emulsion from Valerio's file with kinematic variables

void momentumreco_emulsion(){
    TFile *inputfile = TFile::Open("annotated_ds_data_result_wihkinematics.root");
    TTreeReader bdtreader("ds",inputfile);

    //defining variables to study
    TTreeReaderArray<Int_t> eventID(bdtreader,"dsvtx_vtx2_mc_ev");
    TTreeReaderArray<Int_t> trackID(bdtreader,"dsvtx_vtx2_mc_tid");

    TTreeReaderArray<Int_t> charmdaughter(bdtreader,"dsvtx_vtx2_trk_charmdaughter");
    TTreeReaderArray<Int_t> samevent(bdtreader,"dsvtx_vtx2_trk_samevent");
    TTreeReaderArray<Int_t> positivedz(bdtreader,"dsvtx_vtx2_trk_positivedz");

    TTreeReaderArray<Float_t> pmsang(bdtreader,"dsvtx_vtx2_trk_pms_ang");
    TTreeReaderArray<Float_t> momtrue(bdtreader,"dsvtx_vtx2_trk_mc_mom");

    TH2F *hpeff = new TH2F("hpeff","Relative momentum resolution;P_{true};(P_{reco}-P_{true})/P_{true}",100,0,100,55,-1,10);
    TH2F *hpreco_ptrue = new TH2F("hpreco_ptrue","P_reco vs P_true;P_true[GeV];P_rec[GeV]",1000,0,100,200,0,20);
    TH2F *hpwasfound = new TH2F("hpwasfound","Was it found downsream?",100,0,100,2,0,2);
    
    //file with daughter distribution information
    TFile *distfile = TFile::Open("../distributions_mctrue_withdaughters.root");
    TTreeReader distreader("charmdecays",distfile);

    TTreeReaderArray<Int_t> daugh_trackID(distreader,"daugh_trackid");
    TTreeReaderArray<Int_t> daugh_founddownstream(distreader,"daugh_founddownstream");
    
    const int nentries = bdtreader.GetEntries();
    int ntrks;
    int wasfound;

    bool signaltrack;

    for (int ientry = 0; ientry < nentries; ientry++){
        bdtreader.SetEntry(ientry);
        ntrks = pmsang.GetSize();
        for (int itrk = 0; itrk < ntrks; itrk++){
            //looking at distributions for that eventID
            distreader.SetEntry(eventID[itrk]);
            //looping over all daughters to find the one with same ID -> was it found in SciFi?
            signaltrack = false;
            if (charmdaughter[itrk] && samevent[itrk] && positivedz[itrk]) signaltrack = true;
            wasfound = 0;
            for (int idaughter = 0; idaughter < daugh_trackID.GetSize(); idaughter++){
               if(daugh_trackID[idaughter] == trackID[itrk]) wasfound = daugh_founddownstream[idaughter];
            }

            if (pmsang[itrk] > 0.&&signaltrack){ 
                hpeff->Fill(momtrue[itrk],(pmsang[itrk] - momtrue[itrk])/momtrue[itrk]);
                hpwasfound->Fill(momtrue[itrk],wasfound);
                hpreco_ptrue->Fill(momtrue[itrk],pmsang[itrk]);
            }
        }
    }
    
    hpeff->Draw("COLZ");
    TH1D *hpres = hpeff->ProjectionY("hpres",1,5);
    TCanvas *cres = new TCanvas();
    hpres->Draw();
    TCanvas *ccomp = new TCanvas();
    hpreco_ptrue->Draw("COLZ");
    TCanvas *cwasfound = new TCanvas();
    hpwasfound->Draw("COLZ");

}

void read_dataframe(){
    TFile *inputfile = TFile::Open("annotated_ds_data_result_2.root");
    TTree *inputtree = (TTree*) inputfile->Get("ds");
    
    string lowmomentumbackground = "dsvtx_vtx2_alllowmomentum && !dsvtx_vtx2_charmdaughtersamevent && dsvtx_vtx2_prec  > 0. && dsvtx_vtx2_ntrk > 1";
    string lowmomentumsignal = "dsvtx_vtx2_charmdaughtersamevent && dsvtx_vtx2_alllowmomentum && dsvtx_vtx2_prec > 0. && dsvtx_vtx2_ntrk > 1";
    
    TCanvas *cprec = new TCanvas();
    inputtree->Draw("dsvtx_vtx2_prec>>hprecbkg(25)",lowmomentumbackground.data());
    inputtree->Draw("dsvtx_vtx2_prec>>hprecsig(25)",lowmomentumsignal.data());

    TH1D * hprecbkg = (TH1D*) gDirectory->FindObject("hprecbkg");
    TH1D * hprecsig = (TH1D*) gDirectory->FindObject("hprecsig");

    hprecbkg->Scale(1./hprecbkg->Integral());
    hprecsig->Scale(1./hprecsig->Integral());
    
    hprecsig->SetLineColor(kRed);
    hprecbkg->SetTitle("Background vertices");
    hprecsig->SetTitle("Signal vertices");

    hprecbkg->Draw("hist");
    hprecsig->Draw("hist&&SAMES");
    cprec->BuildLegend();

    hprecbkg->SetTitle("Reconstructed momentum for low momentum vertices;Prec[GeV/c]");

    TCanvas *cmassinv = new TCanvas();
    inputtree->Draw("dsvtx_vtx2_mass_inv>>hmassinvbkg(50)",lowmomentumbackground.data());
    inputtree->Draw("dsvtx_vtx2_mass_inv>>hmassinvsig(50)",lowmomentumsignal.data());
    
    TH1D * hmassinvbkg = (TH1D*) gDirectory->FindObject("hmassinvbkg");
    TH1D * hmassinvsig = (TH1D*) gDirectory->FindObject("hmassinvsig");
    
    hmassinvbkg->Scale(1./hmassinvbkg->Integral());
    hmassinvsig->Scale(1./hmassinvsig->Integral());
    
    hmassinvsig->SetLineColor(kRed);
    hmassinvbkg->SetTitle("Background vertices");
    hmassinvsig->SetTitle("Signal vertices");

    hmassinvbkg->Draw("hist");
    hmassinvsig->Draw("hist&&SAMES");
    cmassinv->BuildLegend();

    hmassinvbkg->SetTitle("Reconstructed invariant mass for low momentum vertices;Mass Inv[GeV]");

    
    TCanvas *cminmassinv = new TCanvas();
    inputtree->Draw("dsvtx_vtx2_min_mass_inv>>hminmassinvbkg(50)",lowmomentumbackground.data());
    inputtree->Draw("dsvtx_vtx2_min_mass_inv>>hminmassinvsig(50)",lowmomentumsignal.data());

    TH1D * hminmassinvbkg = (TH1D*) gDirectory->FindObject("hminmassinvbkg");
    TH1D * hminmassinvsig = (TH1D*) gDirectory->FindObject("hminmassinvsig");

    hminmassinvbkg->Scale(1./hminmassinvbkg->Integral());
    hminmassinvsig->Scale(1./hminmassinvsig->Integral());
    
    hminmassinvsig->SetLineColor(kRed);
    hminmassinvbkg->SetTitle("Background vertices");
    hminmassinvsig->SetTitle("Signal vertices");

    hminmassinvbkg->Draw("hist");
    hminmassinvsig->Draw("hist&&SAMES");
    cminmassinv->BuildLegend();

    hminmassinvbkg->SetTitle("Reconstructed minimum invariant mass for low momentum vertices;Min Mass Inv[Gev]");
}
