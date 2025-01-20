//Compute xy spatial distributions of muons after hadron stopper (Created by A. Iuliano 20 January 2025)

void SpatialDistributions(){
    //open input files for reading as a TChain
    TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
    TChain *simchain = new TChain("muon_mcpoints");

    //const int nfiles = 67; //number of input files
    const int nfiles = 1;
    cout<<"reading file "<<endl;
    for (int i = 0; i < nfiles; i++){
        simchain->Add(Form("/home/utente/Simulations/background-prod-2018/vetopoints_muons/vetopoints_scoringplane_muons_withCharmandBeauty%i.root",i*1000));
 }
    const int nentries = simchain->GetEntries();
    cout<<"Number of events"<<nentries<<endl;

    float px, py, pz;
    float x, y, z;
    float trackID, pdgcode;
    float weight;

    //setting branch addresses
    simchain->SetBranchAddress("trackID",&trackID);
    simchain->SetBranchAddress("pdgcode",&pdgcode);

    simchain->SetBranchAddress("x",&x);
    simchain->SetBranchAddress("y",&y);
    simchain->SetBranchAddress("z",&z);
    
    simchain->SetBranchAddress("px",&px);
    simchain->SetBranchAddress("py",&py);
    simchain->SetBranchAddress("pz",&pz);

    simchain->SetBranchAddress("weight",&weight);
    //initialize histograms
    TH2D * hxy = new TH2D("hxy", "Muons xy distribution;x[cm];y[cm]",200,-100,100,200,-100,100);
    TH1D * hz = new TH1D("hz", "Muons z distribution;z[cm]",1000,-7100,-6100);
    TH1D * hp = new TH1D("hp", "Muons momentum;p[GeV/c]",400,0,400);
    //start loop
    double momentum;
    for (int ientry = 0; ientry < nentries; ientry++){
        simchain->GetEntry(ientry); //reading the entry
            
        hxy->Fill(x,y, weight);
                
        momentum = TMath::Sqrt(pow(px,2) + pow(py,2) + pow(pz,2));
        hp->Fill(momentum, weight);    
    } //end event loop

    //draw histograms
    TCanvas *cxy = new TCanvas("cxy","XY Distribution",800,800);
    hxy->Draw("COLZ");

    TCanvas *cp = new TCanvas("cp", "Momentum distribution");
    hp->Draw();

}