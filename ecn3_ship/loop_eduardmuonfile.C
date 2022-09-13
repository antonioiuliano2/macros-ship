//we actually need to loop to correctly assign weights
void loop_eduardmuonfile(){
 //first, opening file, setting ttreereader
 TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/snd_ship_ecn3/muons_3_planes.root");
 TTree *simtree = (TTree*)inputfile->Get("cbmsim");
 TTreeReader reader(simtree);

 //setting tracks and hits branches;
 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<vetoPoint> sco1points(reader,"sco_1Point");

 const int nentries = reader.GetEntries();
 //const int nentries = 100000;

 cout<<"Starting looping over entries "<<nentries<<endl;

 TH2D *hxy = new TH2D("hxy","xy distribution of mcpoints in scoring plane 1;x[cm];y[cm]",40,-200,200,40,-200,200);

 for(int ientry = 0;ientry<nentries;ientry++){    
    reader.SetEntry(ientry);
    //looping over hits
    for (const vetoPoint& sco1point: sco1points){
        //accessing track to get weight
        int trackID = sco1point.GetTrackID();
        double weight = tracks[trackID].GetWeight();

        hxy->Fill(sco1point.GetX(),sco1point.GetY(),weight);

    }
 }
    TCanvas *cxy = new TCanvas("cxy","xy zoomed over region");
    hxy->Draw("COLZ && TEXT");

}