//Compute xy spatial distributions of muons after hadron stopper (Created by A. Iuliano 20 January 2025)
double GetCharge(float PdgCode){
 TDatabasePDG *pdgdata = TDatabasePDG::Instance();
 double charge = (pdgdata->GetParticle(PdgCode)->Charge())/3.; //pdgdatabase returns charge in quark units (i.e. e- is -3.)
 return charge;
}

void SpatialDistributions(){
    TString prefix("root:://eospublic.cern.ch/");//for ROOTXD

    //outputfile
    TFile *outputfile = new TFile("vetopoints_hadronstopper_muons.root","RECREATE"); //actually only muons! Not worked

    //open input files for reading as a TChain
    TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
    TChain *simchain = new TChain("cbmsim");

    //const int nfiles = 67; //number of input files
    
    const int nfiles = 1;
    for (int i = 0; i < nfiles; i++){
        simchain->Add((prefix+TString(Form("/eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_10.0_withCharmandBeauty%i_mu.root",i*1000))).Data());
 }

    TTreeReader reader(simchain);
    TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
    TTreeReaderArray<vetoPoint> vetopoints(reader,"vetoPoint");
    
    //start loop
    const int nentries = reader.GetEntries(true);
    cout<<"Number of events"<<nentries<<endl;
    int pdgcode, trackID;
    double weight;
    double momentum;

    //initialize histograms
    TH2D * hxy = new TH2D("hxy", "Charged xy distribution;x[cm];y[cm]",200,-100,100,200,-100,100);
    TH1D * hz = new TH1D("hz", "Charged z distribution;z[cm]",1000,-7100,-6100);
    TH1D * hp = new TH1D("hp", "Charged momentum;p[GeV/c]",400,0,400);

    for (int ientry = 0; ientry < nentries; ientry++){
        if (ientry % (int)1e+6 == 0) cout<<"Arrived at entry "<<ientry<<endl;
        reader.SetEntry(ientry); //reading the entry
        for (const vetoPoint& vetopoint: vetopoints){
            pdgcode = vetopoint.PdgCode();
            double charge = GetCharge(pdgcode);
            if (TMath::Abs(pdgcode)==13){                 //we only have muons, actually
                trackID = vetopoint.GetTrackID();
                weight = tracks[trackID].GetWeight();

                hxy->Fill(vetopoint.GetX(),vetopoint.GetY(), weight);
                momentum = TMath::Sqrt(pow(vetopoint.GetPx(),2) + pow(vetopoint.GetPy(),2) + pow(vetopoint.GetPz(),2));  
                hp->Fill(momentum, weight);    
            }
        } //end vetoPoint loop
    } //end event loop

    //draw histograms
    TCanvas *cxy = new TCanvas("cxy","XY Distribution",800,800);
    hxy->Draw("COLZ");

    TCanvas *cp = new TCanvas("cp", "Momentum distribution");
    hp->Draw();

    outputfile->cd();
    hxy->Write();
    hp->Write();

}