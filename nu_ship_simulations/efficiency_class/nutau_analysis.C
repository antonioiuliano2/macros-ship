//estimate efficiency for numu CCDIS
void nutau_analysis(){
    //efficiency estimation parameters
    const double offsetxy = 0.1; const int Nminplates = 4; //parameters for geometrical acceptance
    const double maxtheta = 1.; const double minmomentum = 1.; // parameters for primary vertex visibility
    const double posres = 100.*1e-4; const double sagittares = 0.02122; //parameters for magnetic spectrometer
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

    TString prefix("root:://eosuser.cern.ch/");
    TString simfile("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_07_04_nutau_bar_CCDIS_ECN3geom/inECC_ship.conical.Genie-TGeant4.root");
    TString geofile("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_07_04_nutau_bar_CCDIS_ECN3geom/1/geofile_full.conical.Genie-TGeant4.root");
    EfficiencyCut effcut = EfficiencyCut((prefix+simfile).Data(),(prefix+geofile).Data());
    const int nentries = effcut.GetEntries();
    //const int nentries = 1000;
 
    const int ndecaychannels = 4;
 
    int ninbrick[ndecaychannels] = {0,0,0,0};
    double ngeomok[ndecaychannels] = {0.,0.,0.,0.};
    double nvisible[ndecaychannels] = {0.,0.,0.,0.};
    double nspectro[ndecaychannels] = {0.,0.,0.,0.};
    double ndecaysearch[ndecaychannels] = {0.,0.,0.,0.};

    double totalweight[ndecaychannels] = {0.,0.,0.,0.};

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
     //apply visiblevertex location
     bool isvisible = effcut.VisibleVertexLocation(maxtheta, minmomentum);
     //apply momentum measurement
     int checkmcs = effcut.MCSmeasurement(); //-1 if failed, 0 if exited normally
     bool checkdecaysearch = effcut.DecaySearch(true,maxdl,minkinkangle,minip,mindaumomentum,maxdautantheta);
     int decaychannel = effcut.DecayChannel();
     bool checkspectrometer = effcut.SpectrometerAcceptance(posres, sagittares);
     if (decaychannel==0){
        cout<<"ERROR: Unidentified decay channel in event "<<ientry<<endl;;
        continue;
     }
     ninbrick[decaychannel-1]++;
     totalweight[decaychannel-1] += effcut.GetEventWeight();
     //if(checkmcs==0) cout<<effcut.measured_momentum_nudaughterID[1]<<endl;
     //cout<<endl;
     //increasing counters:
     if (geomok == 1){
        ngeomok[decaychannel-1] += effcut.GetEventWeight();
        if (isvisible){
            nvisible[decaychannel-1] += effcut.GetEventWeight();
            if (checkdecaysearch){
             ndecaysearch[decaychannel-1] += effcut.GetEventWeight();
             if (decaychannel==1){ //only for muon channel
              if (checkspectrometer){
                nspectro[decaychannel-1] += effcut.GetEventWeight();                          
                effcut.FillHistograms(hnuP_passed, hq2_x_passed, hmuP_passed);
              }
             } //end of decay channel check
             else effcut.FillHistograms(hnuP_passed, hq2_x_passed, hmuP_passed); //for other decay channels the previous steps are enough
            } //end of decay search check
        } // end of visibility check
     } //end of geometrical acceptance

     effcut.FillHistograms(hnuP, hq2_x, hmuP);
    }
    double alltotalweight = 0.;
    for (int idecaychannel = 0; idecaychannel < ndecaychannels; idecaychannel++){
     cout<<"fraction geom ok "<<ngeomok[idecaychannel]/totalweight[idecaychannel]<<" "<<ngeomok[idecaychannel]<<endl;
     cout<<"fraction visible ok "<<nvisible[idecaychannel]/totalweight[idecaychannel]<<" "<<nvisible[idecaychannel]<<endl;
     cout<<"fraction decay search ok "<<ndecaysearch[idecaychannel]/totalweight[idecaychannel]<<" "<<ndecaysearch[idecaychannel]<<endl;
     cout<<"fraction spectro ok "<<nspectro[idecaychannel]/totalweight[idecaychannel]<<" "<<nspectro[idecaychannel]<<endl;
     cout<<"total weight"<<totalweight[idecaychannel]<<" total events considered "<<ninbrick[idecaychannel]<<endl;

     alltotalweight += totalweight[idecaychannel];
    }
    //drawing histograms, store the passed q2_x into file
    const double normalization = 2.34e+04; //numuCHARMCCDIS per year
    const int nyears = 15; //number of years to process
    
    TCanvas *cq2x = new TCanvas();
    hq2_x->Scale(1./alltotalweight * normalization * nyears);
    hq2_x->Draw("COLZ");
    cq2x->SetLogz();

    TCanvas *cq2x_passed = new TCanvas();
    hq2_x_passed->Scale(1./alltotalweight * normalization * nyears);
    hq2_x_passed->Draw("COLZ");
    cq2x_passed->SetLogz();
    //looping over bins
    int nbinsx = hq2_x_passed->GetNbinsX();
    int nbinsy = hq2_x_passed->GetNbinsY();
    
    std::ofstream ofs;
    ofs.open ("nutauccdis_q2x_passed_histogram.txt", std::ofstream::out);
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
    hnuP->Scale(1./alltotalweight);
    hnuP->Draw();
    hnuP_passed->Scale(1./alltotalweight);
    hnuP_passed->Draw("SAMES");

}