//estimate efficiency for numu CCDIS
void numu_CHARMCCDIS_analysis(){

    cout<<"Here spectrometer acceptance cuts DO NOT REQUIRE OTHER CUTS"<<endl;
    cout<<"MuPassed requires also not DS passed"<<endl;
    //efficiency estimation parameters
    const double offsetxy = 0.1; const int Nminplates = 4; //parameters for geometrical acceptance
    const double maxtheta = 1.; const double minmomentum = 1.; // parameters for primary vertex visibility
    const double posres = 100.*1e-4; const double sagittares = 0.02122; //parameters for magnetic spectrometer
    const double dsposres = 100.*1e-4; const double dssagittares = 0.02122; //parameters for magnetic spectrometer
    //decay seach parameters
    const double maxdl = 0.4; //max decay length
    const double minkinkangle = 0.02; //minimum angle between decaying particle and daughter
    const double minip = 10e-4; //minimum impact parameter of daughter with respect to primary vertex
    const double mindaumomentum = 0.1; //min momentum of decay daughter
    const double maxdautantheta = 1.; //max tan theta decay daughter

    //histograms
    TH1D *hnuP = new TH1D("hnuP","Neutrino momentum;P[GeV/C]",400,0,400);
    TH2D *hq2_x = new TH2D("hq2_x",";log10(x);log10(Q2)",20,-4.,1.,10,-1.5,3.5);
    TH1D *hmuP = new TH1D("hmuP","All muons;P[GeV/C]",400,0,400);

    //histograms_passed
    TH1D *hnuP_passed = new TH1D("hnuP_passed","Neutrino momentum;P[GeV/C]",400,0,400);
    TH2D *hq2_x_passed = new TH2D("hq2_x_passed",";log10(x);log10(Q2)",20,-4.,1.,10,-1.5,3.5);
    TH1D *hmuP_passed = new TH1D("hmuP_passed","Muons measured in SND Spectrometer but NOT entering in HS spectrometer;P[GeV/C]",400,0,400);

   //histograms_passed (in SHiP Decay Spectrometer instead of SND Muon Spectrometer)
    TH1D *hnuP_dspassed = new TH1D("hnuP_dspassed","Accepted events in DS;P[GeV/C]",400,0,400);
    TH2D *hq2_x_dspassed = new TH2D("hq2_x_dspassed",";log10(x);log10(Q2)",20,-4.,1.,10,-1.5,3.5);
    TH1D *hmuP_dspassed = new TH1D("hmuP_dspassed","Muons entering in HS Spectrometer;P[GeV/C]",400,0,400);


    TString prefix("root:://eosuser.cern.ch/");
    TString simfile("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_07_09_numu_CHARMCCDIS_ECN3geom/inECC_ship.conical.Genie-TGeant4.root");
    TString geofile("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_07_09_numu_CHARMCCDIS_ECN3geom/geofile_full.conical.Genie-TGeant4.root");
    EfficiencyCut effcut = EfficiencyCut((prefix+simfile).Data(),(prefix+geofile).Data());
    const int nentries = effcut.GetEntries();
    //const int nentries = 1000;
    int ninbrick = 0;

    double ngeomok = 0., nvisible = 0., nspectro = 0., ndecayspectro = 0, ndecaysearch = 0.;

    double totalweight = 0.;

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
     bool checkdecaysearch = effcut.DecaySearch(false,maxdl,minkinkangle,minip,mindaumomentum,maxdautantheta);
     bool checkspectrometer = effcut.SpectrometerAcceptance(posres, sagittares);
     bool checkdecayspectrometer = effcut.DecaySpectrometerAcceptance(dsposres, dssagittares);
     //if(checkmcs==0) cout<<effcut.measured_momentum_nudaughterID[1]<<endl;
     //cout<<endl;
     //increasing counters:
     if (geomok == 1){
        ngeomok += effcut.GetEventWeight();
        if (isvisible){
            nvisible += effcut.GetEventWeight();
            if (checkdecaysearch){
             ndecaysearch += effcut.GetEventWeight();
            }            
        }
     }
    
     if ((checkspectrometer)&&(!checkdecayspectrometer)){ //the spectro measurement requires all the ones before
                nspectro += effcut.GetEventWeight();            
                effcut.FillHistograms(hnuP_passed, hq2_x_passed, hmuP_passed);
             }

     if (checkdecayspectrometer){ //entering HS alone!
                ndecayspectro += effcut.GetEventWeight();
                effcut.FillHistograms(hnuP_dspassed, hq2_x_dspassed, hmuP_dspassed);
     }

     effcut.FillHistograms(hnuP, hq2_x, hmuP);
    }
    cout<<"fraction geom ok "<<ngeomok/totalweight<<" "<<ngeomok<<endl;
    cout<<"fraction visible ok "<<nvisible/totalweight<<" "<<nvisible<<endl;
    cout<<"fraction decay search ok "<<ndecaysearch/totalweight<<" "<<ndecaysearch<<endl;
    cout<<"fraction spectro ok "<<nspectro/totalweight<<" "<<nspectro<<endl;
    cout<<"fraction spectro in ds ok "<<ndecayspectro/totalweight<<" "<<ndecayspectro<<endl;    
    cout<<"total weight"<<totalweight<<" total events considered "<<ninbrick<<endl;
    //drawing histograms, store the passed q2_x into file
    const double normalization = 2.34e+04; //numuCHARMCCDIS per year
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
    ofs.open ("numucharmccdis_q2x_passed_histogram.txt", std::ofstream::out);
    ofs<<"x q2 Enu #nu"<<endl;
    
    for (int ibinx = 1; ibinx<= nbinsx; ibinx++)    {
     for (int ibiny = 1; ibiny<= nbinsy; ibiny++){
        double xcenter = pow(10,hq2_x_passed->GetXaxis()->GetBinCenter(ibinx));
        double ycenter = pow(10,hq2_x_passed->GetYaxis()->GetBinCenter(ibiny));
        double nreco = hq2_x->GetBinContent(ibinx,ibiny);
        //global bin number, used for map
        int whichbin = hq2_x_passed->GetBin(ibinx,ibiny);
        double meanenergy = ROOT::VecOps::Mean(effcut.energies_q2x[whichbin]);
        ofs<<xcenter<<" "<<ycenter<<" "<<meanenergy<<" "<<nreco<<endl;   
     }
    }
    ofs.close();

    TCanvas *cnuP = new TCanvas();
    hnuP->Scale(1./totalweight);
    hnuP->Draw("histo");
    hnuP_passed->SetFillColor(kBlue);
    hnuP_passed->SetFillStyle(3305);
    hnuP_passed->Scale(1./totalweight);
    hnuP_passed->Draw("SAMES&&histo");
    hnuP_dspassed->SetLineColor(kRed);
    hnuP_dspassed->SetFillColor(kRed);
    hnuP_dspassed->SetFillStyle(3315);
    hnuP_dspassed->Scale(1./totalweight);
    hnuP_dspassed->Draw("SAMES&&histo");
    cnuP->BuildLegend();

    TCanvas *cmuP = new TCanvas();
    auto hmuStack = new THStack("hmuStack","Muons entering DS vs not entering but measured in SND");
    hmuP->Scale(1./totalweight);
    hmuP->Draw("histo");
    hmuP_passed->SetFillColor(kBlue);
    hmuP_passed->SetFillStyle(3305);
    hmuP_passed->Scale(1./totalweight);
    //hmuP_passed->Draw("SAMES&&histo");
    hmuP_dspassed->SetLineColor(kRed);
    hmuP_dspassed->SetFillColor(kRed);
    hmuP_dspassed->SetFillStyle(3315);
    hmuP_dspassed->Scale(1./totalweight);
    //hmuP_dspassed->Draw("SAMES&&histo");
    hmuStack->Add(hmuP_dspassed);
    hmuStack->Add(hmuP_passed);
    hmuStack->Draw("SAMES&&histo");
    cmuP->BuildLegend();

}