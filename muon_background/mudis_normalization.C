void mudis_normalization(){

    TString prefix("root:://eosuser.cern.ch/");
    TString simpath_mu("/eos/user/a/aiuliano/public/sims_FairShip/sim_muon_background/sim_MuDIS/2024_12_09_file2000_1000jobs/");
    TString sigmapath("/home/utente/Simulations/muon_background_production/file_2000_sigmadata/");

    TH1D *hw1 = new TH1D("hw1","Weight 1",100,0,800);
    TH1D *hw2 = new TH1D("hw2","Weight 2",100,0,4000);
    TH1D *hsigma = new TH1D("hsigma","DIS Sigma",1000,0,0.001);

    const int startfile = 0;
    const int endfile = 1;
    //preparing TChain
    TChain *simchain = new TChain("cbmsim");
    for (int ifile = startfile; ifile <= endfile; ifile++){
        simchain->Add((prefix+simpath_mu+TString(Form("%i/ship.conical.muonDIS-TGeant4.root",ifile))).Data());
    }

    TTreeReader reader(simchain);
    TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
    const int nentries = reader.GetEntries(true);
    
    cout<<"Start loop over number of entries "<<nentries<<endl;
    
    fstream sigmafile;

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
        hw1->Fill(w1);
        hw2->Fill(w2);
        hsigma->Fill(sigma);
        double normalization_rate = sigma * w1 * w2;
    }
    //plotting results
    TCanvas *c1 = new TCanvas();
    hw1->Draw();
    TCanvas *c2 = new TCanvas();
    hw2->Draw();
    TCanvas *c3 = new TCanvas();
    hsigma->Draw();
}