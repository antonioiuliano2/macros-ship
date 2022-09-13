using namespace ROOT;

RVec<double> GetX(TClonesArray &pointarray){
    RVec<double> myx;
    const int npoints = pointarray.GetEntries();
    for (int i=0; i < npoints;i++){
        FairMCPoint *point = (FairMCPoint*) pointarray.At(i);
        myx.push_back(point->GetX());
    }
    return myx;
}

RVec<double> GetY(TClonesArray &pointarray){
    RVec<double> myy;
    const int npoints = pointarray.GetEntries();
    for (int i=0; i < npoints;i++){
        FairMCPoint *point = (FairMCPoint*) pointarray.At(i);
        myy.push_back(point->GetY());
    }
    return myy;
}

RVec<double> GetWeight(TClonesArray &trackarray,TClonesArray &pointarray){

    RVec<double> sco1_weights;
    const int npoints = pointarray.GetEntries();
    for (int i=0; i < npoints;i++){
        FairMCPoint *point = (FairMCPoint*) pointarray.At(i);
        //get corresponding track
        ShipMCTrack* mctrack = (ShipMCTrack*) trackarray.At(point->GetTrackID());
        sco1_weights.push_back(mctrack->GetWeight());
    }

    
    return sco1_weights;
}

void read_eduardmuonfile(){

    EnableImplicitMT(); //we need to go parallel
    //first, accessing the file from EOS
    //TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/snd_ship_ecn3/muons_3_planes.root");
    TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/ntuples_for_dis/combi_ecn3/merged_combi_ecn3_spill_11.root");
    //build rdataframe
    RDataFrame df("cbmsim",inputfile);

    //const char* detname = "sco_1Point";
    const char* detname = "UpstreamTaggerPoint";
    auto df1 = df.Define("sco1X",GetX,{detname})
                .Define("sco1Y",GetY,{detname})
                .Define("sco1Weight",GetWeight,{"MCTrack",detname});

    cout<<"Starting rdataframe analysis"<<endl;
    auto hxy = df1.Histo2D({"hxy","xy distribution of muons;x[cm];y[cm]",140,-600,800,140,-600,800},"sco1X","sco1Y","sco1Weight");
    auto hxy_zoomed = df1.Histo2D({"hxy_zoomed","xy distribution of muons;x[cm];y[cm]",400,-200,200,400,-200,200},"sco1X","sco1Y","sco1Weight");
    auto hxy_zoomed_text = df1.Histo2D({"hxy_zoomed_text","xy distribution of muons;x[cm];y[cm]",40,-200,200,40,-200,200},"sco1X","sco1Y","sco1Weight");

    TFile *outputfile = new TFile("plots_muonshitseduard_upstreamBIG.root","RECREATE");
    TCanvas *cxy = new TCanvas("cxy","xy distribution",800,800);
    hxy->DrawClone("COLZ");
    cxy->Write();
    TCanvas *cxy_zoomed = new TCanvas("cxy_zoomed","xy distribution",800,800);
    hxy_zoomed->DrawClone("COLZ");
    cxy_zoomed->Write();
    TCanvas *cxy_zoomed_text = new TCanvas("cxy_zoomed_text","xy distribution with text in bins",800,800);
    hxy_zoomed_text->DrawClone("COLZ&&TEXT");    
    cxy_zoomed_text->Write();


}