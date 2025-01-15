double GetParticleCharge (int pdgcode, TDatabasePDG *pdg){
  //from PDG, get charge
  double charge = 0.;
  if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
  else if (pdgcode > 1e+8) charge = 1.; //test storing heavy nuclei
  return charge;
}

void mudis_normalization(){
    TFile *ouputputfile = new TFile("SND_MuonShield_Scoringplanes_Histograms_1000files_weightedscale.root","RECREATE");

    TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles

    TString prefix("root:://eosuser.cern.ch/");
    TString simpath_mu("/eos/user/a/aiuliano/public/sims_FairShip/sim_muon_background/sim_MuDIS/2024_12_09_file2000_1000jobs/");
    TString sigmapath("/home/utente/Simulations/muon_background_production/file_2000_sigmadata/");

    TH1D *hw1 = new TH1D("hw1","Weight 1",100,0,800);
    TH1D *hw2 = new TH1D("hw2","Weight 2",100,0,4000);
    TH1D *hsigma = new TH1D("hsigma","DIS Sigma",1000,0,0.001);
    TH1D *hnorm = new TH1D("hnorm","Normalization_Rate",100,0,1000);

    //muon starting positions
    TH2D *hxy = new TH2D("hxy","xy of muon starting points;x[cm];y[cm]",120,-600.,600.,100,-500.,500.);

    //histograms for all 10 scoring planes, from 0 to 9    
    const int nscoringplanes = 10;

    double scoringplane_surface = 20. * 40.; //only positive x are stored in the scoring planes
    
    TH2D *hxy_sc[nscoringplanes];
    for (int i = 0; i < nscoringplanes; i++){
        hxy_sc[i] = new TH2D(Form("hxy_sc%i",i),Form("xy of hits in scoring plane #%i;x[cm];y[cm]",i),40,-20.,20.,40,-20.,20.);
    }

    TH2D *hxy_sc0_xpos = new TH2D("xy_sc0_xpos","xy of hits in scoring plane 0 positive x;x[cm];y[cm]",20,0.,20.,40,-20.,20.);

    //binning parameters must be the same for all momentum histograms, since I want to plot them together
    const int nbins_p = 100;
    const double p_min = 0.;
    const double p_max = 10.;

    TH1D *hp_sc0 = new TH1D("hp_sc0","Momentum in last scoring plane;p[GeV/c]",nbins_p ,p_min,p_max);
    TH1D *hp_sc0_charged = new TH1D("hp_sc0_charged","Charged particles",nbins_p ,p_min,p_max);

    TH1D *hp_sc0_neutrons = new TH1D("hp_sc0_neutrons","Neutrons",nbins_p,p_min,p_max);
    TH1D *hp_sc0_gamma = new TH1D("hp_sc0_gamma","Gamma",nbins_p ,p_min,p_max);
    TH1D *hp_sc0_pions_n = new TH1D("hp_sc0_pions_n","Pi0",nbins_p,p_min,p_max);
    
    TH1D *hp_sc0_pions_ch = new TH1D("hp_sc0_pions_ch","Charged Pions",nbins_p ,p_min,p_max);
    TH1D *hp_sc0_protons = new TH1D("hp_sc0_protons","Charged Protons",nbins_p ,p_min,p_max);

    TH1D *hp_sc0_electrons = new TH1D("hp_sc0_electrons","Electrons",nbins_p ,p_min,p_max);
    TH1D *hp_sc0_muons = new TH1D("hp_sc0_muons","Muons",nbins_p ,p_min,p_max);

    const int startfile = 0;
    const int endfile = 1000;
    //preparing TChain
    TChain *simchain = new TChain("cbmsim");
    for (int ifile = startfile; ifile <= endfile; ifile++){
        simchain->Add((prefix+simpath_mu+TString(Form("%i/ship.conical.muonDIS-TGeant4.root",ifile))).Data()); //missing files are considered automatically in tree number, and skipped
    }

    const double spill_weight = 15694520448.0; //total weight from the pkl middle file
    const int nmuons_spill = 463588111; // about 400M muons

    TTreeReader reader(simchain);
    TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
    TTreeReaderArray<UpstreamTaggerPoint> scoringpoints(reader,"UpstreamTaggerPoint");
    const int nentries = reader.GetEntries(true);
    
    cout<<"Start loop over number of entries "<<nentries<<endl;
    
    fstream sigmafile;

    const double N_avogadro = 6.022e+23;

    int treenumber = -1;
    for(int ientry = 0;ientry<nentries;ientry++){
        reader.SetEntry(ientry);        
        if (simchain->GetTreeNumber()>treenumber){
            treenumber = simchain->GetTreeNumber();
            cout<<"Started new tree at entry: "<<ientry<<" "<<treenumber<<endl;            

            if (sigmafile.is_open()) sigmafile.close();
            sigmafile.open((sigmapath+TString(Form("sigmadata_%i.txt",treenumber))).Data(),ios::in);
        }
        double w1 = tracks[0].GetWeight();
        double w2 = tracks[2].GetWeight();

        double sigma;
        sigmafile >> sigma; //taking cross section from input file
        sigma = sigma *1e-27; //from mb to cm2
        
        hw1->Fill(w1);
        hw2->Fill(w2);
        hsigma->Fill(sigma);
        double normalization_rate = sigma * w1 * w2 * N_avogadro;
        hnorm->Fill(normalization_rate);

        double xmuon = tracks[0].GetStartX();
        double ymuon = tracks[0].GetStartY();
        //xy distribution of input muons
        hxy->Fill(xmuon,ymuon,normalization_rate);

        //loop in scoringplanes
         for (const UpstreamTaggerPoint& scoringpoint:scoringpoints){

            int pdgcode = scoringpoint.PdgCode();
            double charge = GetParticleCharge(pdgcode, pdg);

            double px = scoringpoint.GetPx();
            double py = scoringpoint.GetPy();
            double pz = scoringpoint.GetPz();

            double momentum = TMath::Sqrt(px * px + py * py + pz * pz);
            int detID = scoringpoint.GetDetectorID();
            //fill corresponding scoring plane
            if (detID>=0 && detID < 10) hxy_sc[detID]->Fill(scoringpoint.GetX(),scoringpoint.GetY(),normalization_rate);
            //studying last scoring plane in more detail
            if (detID==0.){
                if (scoringpoint.GetX() >= 0.) hxy_sc0_xpos->Fill(scoringpoint.GetX(),scoringpoint.GetY(),normalization_rate);
                //momentum histogram
                hp_sc0->Fill(momentum,normalization_rate);
                if (TMath::Abs(charge)>0.) hp_sc0_charged->Fill(momentum,normalization_rate);
                //according to pdgcode, fill different histograms
                switch (TMath::Abs(pdgcode))
                {
                 case 22: hp_sc0_gamma->Fill(momentum,normalization_rate);
                 break;

                 case 11: hp_sc0_electrons->Fill(momentum,normalization_rate);
                 break;

                 case 13: hp_sc0_muons->Fill(momentum,normalization_rate);
                 break;

                 case 2212: hp_sc0_protons->Fill(momentum,normalization_rate);
                 break;

                 case 2112: hp_sc0_neutrons->Fill(momentum,normalization_rate);
                 break;

                 case 211: hp_sc0_pions_ch->Fill(momentum,normalization_rate);
                 break;
                 
                 case 111: hp_sc0_pions_n->Fill(momentum,normalization_rate);
                 break;
                }//end switch pdgcode
                

            } //end if last scoring plane
         }//Upstream Tagger Point positions loop
    } //event loop

    //computing normalization factor
    double norm_factor = spill_weight/hw1->Integral();
    //double norm_factor = nmuons_spill/nentries;
    cout<<"Normalization factor for one spill "<< norm_factor<<endl;
    
    //plotting results, saving histograms
    ouputputfile->cd();
    TCanvas *c1 = new TCanvas("c1","Weight 1");
    hw1->Draw();
    TCanvas *c2 = new TCanvas("c2","Weight 2");
    hw2->Draw();
    TCanvas *c3 = new TCanvas("c3","Cross section from Pythia6");
    hsigma->Draw();
    TCanvas *c4 = new TCanvas("c4","Normalization Rate");
    hnorm->Draw();
    TCanvas *c5 = new TCanvas("c5","Muons XY Starting Positions",800,800);
    hxy->Scale(norm_factor);
    hxy->Draw("COLZ");
    hxy->Write();
    TCanvas *c6 = new TCanvas("c6","MC XY Scoring Plane 0",800,800);
    hxy_sc[0]->Scale(norm_factor);
    hxy_sc[0]->Draw("COLZ");
    hxy_sc[0]->Write();
    cout<<"Hits in scoring plane "<<0<<" are "<<hxy_sc[0]->Integral()<<" density "<<hxy_sc[0]->Integral()/scoringplane_surface<<endl;

    TCanvas *c6_pos = new TCanvas("c6_pos","MC XY Scoring Plane 0 X>=0",800,800);
    hxy_sc0_xpos->Scale(norm_factor);
    hxy_sc0_xpos->Draw("COLZ");
    hxy_sc0_xpos->Write();
    cout<<"Hits in scoring plane "<<0<<" with positive x are "<<hxy_sc0_xpos->Integral()<<" density "<<hxy_sc0_xpos->Integral()/scoringplane_surface<<endl;

    TCanvas *c_sc = new TCanvas("c_sc","Scoring Planes distributions XY",800,800);
    c_sc->Divide(3,3);
    for (int i=1; i < nscoringplanes; i++){
        c_sc->cd(i);
        hxy_sc[i]->Scale(norm_factor);
        hxy_sc[i]->Draw("COLZ");
        hxy_sc[i]->Write();

        cout<<"Hits in scoring plane "<<i<<" are "<<hxy_sc[i]->Integral()<<" density "<<hxy_sc[i]->Integral()/scoringplane_surface<<endl;
    }

    //momentum plot
    TCanvas *c9 = new TCanvas("c9","Momentum particles at scoring plane 0");

    hp_sc0->SetLineColor(kBlue);
    hp_sc0->Scale(norm_factor);
    hp_sc0->Draw("hist");
    hp_sc0->Write();

    hp_sc0_charged->SetLineColor(kRed);
    hp_sc0_charged->Scale(norm_factor);
    hp_sc0_charged->Draw("hist&&SAMES");
    hp_sc0_charged->Write();

    hp_sc0_electrons->SetLineColor(kCyan);
    hp_sc0_electrons->Scale(norm_factor);
    hp_sc0_electrons->Draw("hist&&SAMES");
    hp_sc0_electrons->Write();

    hp_sc0_muons->SetLineColor(kTeal);
    hp_sc0_muons->Scale(norm_factor);
    hp_sc0_muons->Draw("hist&&SAMES");
    hp_sc0_muons->Write();
    
    hp_sc0_pions_ch->SetLineColor(kMagenta);
    hp_sc0_pions_ch->Scale(norm_factor);
    hp_sc0_pions_ch->Draw("hist&&SAMES");
    hp_sc0_pions_ch->Write();

    hp_sc0_gamma->SetLineColor(kYellow);
    hp_sc0_gamma->Scale(norm_factor);
    hp_sc0_gamma->Draw("hist&&SAMES");
    hp_sc0_gamma->Write();

    hp_sc0_pions_n->SetLineColor(kOrange);
    hp_sc0_pions_n->Scale(norm_factor);
    hp_sc0_pions_n->Draw("hist&&SAMES");    
    hp_sc0_pions_n->Write();
    
    hp_sc0_protons->SetLineColor(kBlack);
    hp_sc0_protons->Scale(norm_factor);
    hp_sc0_protons->Draw("hist&&SAMES");
    hp_sc0_protons->Write();

    hp_sc0_neutrons->SetLineColor(kGreen);
    hp_sc0_neutrons->Scale(norm_factor);
    hp_sc0_neutrons->Draw("hist&&SAMES");
    hp_sc0_neutrons->Write();

    c9->BuildLegend();

}