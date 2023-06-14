//estimate efficiency for numu CCDIS
void numu_analysis(){
    EfficiencyCut effcut = EfficiencyCut(
        "root:://eosuser.cern.ch//eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_05_26_numu_CCDIS_spectrosagitta/inECC_ship.conical.Genie-TGeant4.root","root:://eosuser.cern.ch//eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_05_26_numu_CCDIS_spectrosagitta/1/geofile_full.conical.Genie-TGeant4.root");
    int nentries = effcut.GetEntries();
    cout<<"looping over "<<nentries<<" entries "<<endl;
    for (int ientry = 0; ientry < nentries; ientry++){//loop over entries
     effcut.GetEntry(ientry);
     cout<<ientry<<endl;
     cout<<effcut.GeometricalEfficiency()<<endl;
     cout<<effcut.VisibleVertexLocation()<<endl;
     int checkmcs = effcut.MCSmeasurement(); //-1 if failed, 0 if exited normally
     if(checkmcs==0) cout<<effcut.measured_momentum_nudaughterID[1]<<endl;
     cout<<endl;
    }



}