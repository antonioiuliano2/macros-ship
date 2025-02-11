//Extract points from muons in scoring planes in MuonShield (sco_NPoint) (Created by A. Iuliano 11 February 2025)
void Extract_EduardScoringPoints2024_Muons(int ifile){
    TString prefix("root:://eospublic.cern.ch/");//for ROOTXD

    //outputfile
    TFile *outputfile = new TFile(Form("/home/utente/Simulations/scoringpoints_muons/scoringplane_muons_%i.root",ifile),"RECREATE");
    TNtuple * muon_mcpoints = new TNtuple("muon_mcpoints","MC Point in scoring plane in Muon Shield","trackID:pdgcode:x:y:z:px:py:pz:weight");

    //open input files for reading as a TChain
    TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
    TChain *simchain = new TChain("cbmsim");

    //const int nfiles = 67; //number of input files
    simchain->Add((prefix+TString(Form("/eos/experiment/ship/user/edursov/pycondor_out/muons_in_ms_sc_unrolled_spill_full_info_for_warm_part/11/%i/ship.conical.MuonBack-TGeant4.root",ifile))).Data());

    TTreeReader reader(simchain);
    TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
    TTreeReaderArray<vetoPoint> scopoints(reader,"sco_66Point");
    
    //start loop
    const int nentries = reader.GetEntries(true);
    cout<<"Number of events"<<nentries<<endl;
    int pdgcode, trackID;
    double weight;
    double momentum;
    for (int ientry = 0; ientry < nentries; ientry++){
        if (ientry % (int)1e+6 == 0) cout<<"Arrived at entry "<<ientry<<endl;
        reader.SetEntry(ientry); //reading the entry
        for (const vetoPoint& scopoint: scopoints){
            pdgcode = scopoint.PdgCode();          
            trackID = scopoint.GetTrackID();
            weight = tracks[trackID].GetWeight();
            muon_mcpoints->Fill(trackID,pdgcode,scopoint.GetX(), scopoint.GetY(), scopoint.GetZ(), scopoint.GetPx(), scopoint.GetPy(), scopoint.GetPz(), weight);                
        } //end scopoint loop
    } //end event loop

    //store output file
    outputfile->cd();
    muon_mcpoints->Write();
    outputfile->Close();

}