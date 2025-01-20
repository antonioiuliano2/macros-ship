//Extract points from muons in scoring plane after hadron stopper (vetopoint) (Created by A. Iuliano 20 January 2025)

void ExtractVetoPointMuons(int ifile){
    TString prefix("root:://eospublic.cern.ch/");//for ROOTXD

    //outputfile
    TFile *outputfile = new TFile(Form("/home/utente/Simulations/background-prod-2018/vetopoints_muons/vetopoints_scoringplane_muons_withCharmandBeauty%i.root",ifile*1000),"RECREATE");
    TNtuple * muon_mcpoints = new TNtuple("muon_mcpoints","MC Point in scoring plane after hadron stopper","trackID:pdgcode:x:y:z:px:py:pz:weight");

    //open input files for reading as a TChain
    TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
    TChain *simchain = new TChain("cbmsim");

    //const int nfiles = 67; //number of input files
    simchain->Add((prefix+TString(Form("/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_10.0_withCharmandBeauty%i_mu.root",ifile*1000))).Data());
    //const int nfiles = 1;
    //for (int i = 0; i < nfiles; i++){
    //    simchain->Add((prefix+TString(Form("/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_10.0_withCharmandBeauty%i_mu.root",i*1000))).Data());
 //}

    TTreeReader reader(simchain);
    TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
    TTreeReaderArray<vetoPoint> vetopoints(reader,"vetoPoint");
    
    //start loop
    const int nentries = reader.GetEntries(true);
    cout<<"Number of events"<<nentries<<endl;
    int pdgcode, trackID;
    double weight;
    double momentum;
    for (int ientry = 0; ientry < nentries; ientry++){
        if (ientry % (int)1e+6 == 0) cout<<"Arrived at entry "<<ientry<<endl;
        reader.SetEntry(ientry); //reading the entry
        for (const vetoPoint& vetopoint: vetopoints){
            pdgcode = vetopoint.PdgCode();
            if (TMath::Abs(pdgcode)==13){                
                trackID = vetopoint.GetTrackID();
                weight = tracks[trackID].GetWeight();
                muon_mcpoints->Fill(trackID,pdgcode,vetopoint.GetX(), vetopoint.GetY(), vetopoint.GetZ(), vetopoint.GetPx(), vetopoint.GetPy(), vetopoint.GetPz(), weight);    
            }
        } //end vetoPoint loop
    } //end event loop

    //store output file
    outputfile->cd();
    muon_mcpoints->Write();
    outputfile->Close();

}