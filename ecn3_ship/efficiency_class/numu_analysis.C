//estimate efficiency for numu CCDIS
void numu_analysis(){
    //efficiency estimation parameters
    const double offsetxy = 0.1; const int Nminplates = 4; //parameters for geometrical acceptance
    const double maxtheta = 1.; const double minmomentum = 1.; // parameters for primary vertex visibility
    const double posres = 100.*1e-4; const double sagittares = 0.02122; //parameters for magnetic spectrometer

    //histograms
    TH1D *hnuP = new TH1D("hnuP","Neutrino momentum;P[GeV/C]",400,0,400);
    TH2D *hq2_x = new TH2D("hq2_x",";log10(x);log10(Q2)",20,-4.,1.,10,-1.5,3.5);

    //histograms_passed
    TH1D *hnuP_passed = new TH1D("hnuP_passed","Neutrino momentum;P[GeV/C]",400,0,400);
    TH2D *hq2_x_passed = new TH2D("hq2_x_passed",";log10(x);log10(Q2)",20,-4.,1.,10,-1.5,3.5);

    TString prefix("root:://eosuser.cern.ch/");
    TString simfile("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_06_17_numu_CCDIS_ECN3geom/inECC_ship.conical.Genie-TGeant4.root");
    TString geofile("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_06_17_numu_CCDIS_ECN3geom/1/geofile_full.conical.Genie-TGeant4.root");
    EfficiencyCut effcut = EfficiencyCut((prefix+simfile).Data(),(prefix+geofile).Data());
    const int nentries = effcut.GetEntries();
    //const int nentries = 1000;
    int ninbrick = 0;

    double ngeomok = 0, nvisible = 0, nspectro = 0;

    double totalweight = 0;

    cout<<"looping over "<<nentries<<" entries "<<endl;
    for (int ientry = 0; ientry < nentries; ientry++){//loop over entries
     effcut.GetEntry(ientry);
     //cout<<ientry<<endl;
     //apply geometrical efficiency
     int geomok = effcut.GeometricalEfficiency(offsetxy, Nminplates);
     if (geomok < 0){ 
        //cout<<"No brick found!, please check event "<<ientry<<endl;
        continue; //it is not in the brick
     }
     totalweight += effcut.GetEventWeight();
     ninbrick++;
     //apply visiblevertex location
     bool isvisible = effcut.VisibleVertexLocation(maxtheta, minmomentum);
     //apply momentum measurement
     int checkmcs = effcut.MCSmeasurement(); //-1 if failed, 0 if exited normally
     bool checkspectrometer = effcut.SpectrometerAcceptance(posres, sagittares);
     //if(checkmcs==0) cout<<effcut.measured_momentum_nudaughterID[1]<<endl;
     //cout<<endl;
     //increasing counters:
     if (geomok == 1){
        ngeomok += effcut.GetEventWeight();
        if (isvisible){
            nvisible += effcut.GetEventWeight();
            if (checkspectrometer){
                nspectro += effcut.GetEventWeight();
                effcut.FillHistograms(hnuP_passed, hq2_x_passed);
            }
        }
     }

     effcut.FillHistograms(hnuP, hq2_x);
    }
    cout<<"fraction geom ok "<<ngeomok/totalweight<<" "<<ngeomok<<endl;
    cout<<"fraction visible ok "<<nvisible/totalweight<<" "<<nvisible<<endl;
    cout<<"fraction spectro ok "<<nspectro/totalweight<<" "<<nspectro<<endl;
    cout<<"total weight"<<totalweight<<" total events considered "<<ninbrick<<endl;
    //drawing histograms, store the passed q2_x into file
    const double normalization = 5.33e+05; //numuCCDIS per year
    const int nyears = 15; //number of years to process
    
    TCanvas *cq2x = new TCanvas();
    hq2_x->Scale(1./totalweight * normalization * nyears);
    hq2_x->Draw("COLZ");
    cq2x->SetLogz();

    TCanvas *cq2x_passed = new TCanvas();
    hq2_x_passed->Scale(1./totalweight * normalization * nyears);
    hq2_x_passed->Draw("COLZ");
    cq2x_passed->SetLogz();
    //looping over bins
    int nbinsx = hq2_x_passed->GetNbinsX();
    int nbinsy = hq2_x_passed->GetNbinsY();
    
    std::ofstream ofs;
    ofs.open ("numu_q2x_passed_histogram.txt", std::ofstream::out);
    ofs<<"x q2 #nu"<<endl;
    
    for (int ibinx = 1; ibinx<= nbinsx; ibinx++)    {
     for (int ibiny = 1; ibiny<= nbinsy; ibiny++){
        double xcenter = pow(10,hq2_x_passed->GetXaxis()->GetBinCenter(ibinx));
        double ycenter = pow(10,hq2_x_passed->GetYaxis()->GetBinCenter(ibiny));
        double nreco = hq2_x->GetBinContent(ibinx,ibiny);
        ofs<<xcenter<<" "<<ycenter<<" "<<nreco<<endl;
     }
    }
    ofs.close();

    TCanvas *cnuP = new TCanvas();
    hnuP->Scale(1./totalweight);
    hnuP->Draw();
    hnuP_passed->Scale(1./totalweight);
    hnuP_passed->Draw("SAMES");

}