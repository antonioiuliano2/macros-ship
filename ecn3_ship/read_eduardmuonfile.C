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

RVec<double> GetP(TClonesArray &pointarray){
    RVec<double> myP;
    const int npoints = pointarray.GetEntries();
    for (int i=0; i < npoints; i++){
        FairMCPoint *point = (FairMCPoint*) pointarray.At(i);
        myP.push_back(TMath::Sqrt(point->GetPx()* point->GetPx() + point->GetPy() * point->GetPy() + point->GetPz() * point->GetPz()));
    }

    return myP;
}

void read_eduardmuonfile(){

    EnableImplicitMT(); //we need to go parallel
    //first, accessing the file from EOS
    //TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/snd_ship_ecn3/muons_3_planes.root");
    //TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/ntuples_for_dis/combi_ecn3/merged_combi_ecn3_spill_11.root");
    TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/ntuples_for_dis/combi_ecn3_fixed_overlaps_new_field/merged_combi_ecn3_spill_2.root");
    //build rdataframe
    RDataFrame df("cbmsim",inputfile);

    //const char* detname = "sco_1Point";
    const char* detname = "UpstreamTaggerPoint";
    auto df1 = df.Define("sco1X",GetX,{detname})
                .Define("sco1Y",GetY,{detname})
                .Define("sco1P",GetP,{detname})
                .Define("sco1Weight",GetWeight,{"MCTrack",detname});

    cout<<"Starting rdataframe analysis"<<endl;
    auto hxy = df1.Histo2D({"hxy","xy distribution of muons;x[cm];y[cm]",140,-600,800,140,-600,800},"sco1X","sco1Y","sco1Weight");
    auto hxy_zoomed = df1.Histo2D({"hxy_zoomed","xy distribution of muons;x[cm];y[cm]",400,-200,200,400,-200,200},"sco1X","sco1Y","sco1Weight");
    auto hxy_zoomed_text = df1.Histo2D({"hxy_zoomed_text","xy distribution of muons;x[cm];y[cm]",40,-200,200,40,-200,200},"sco1X","sco1Y","sco1Weight");

    auto hPmap = df1.Profile2D({"hPmap","Momentum map of muons;x[cm];y[cm];P[GeV/c]",20,-100,100,20,-100,100,0,14000},"sco1X","sco1Y","sco1P");

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

    TCanvas *cP_map = new TCanvas("cP_map","Momentum map of muons",800,800);
    hPmap->DrawClone("COLZ");
    cP_map->Write();


}

void plothistograms(){
    TFile *histfile = TFile::Open("plots_muonshitseduard_upstreamBIG.root");
    TCanvas *cxy = (TCanvas*) histfile->Get("cxy_zoomed_text");

    TH2D *hxy_zoomed_text = (TH2D*) cxy->GetPrimitive("hxy_zoomed_text");

    TCanvas *cnew = new TCanvas("cnew","Muon occupancy in one spill",800,800);
    hxy_zoomed_text->GetXaxis()->SetRangeUser(-50,50);
    hxy_zoomed_text->GetYaxis()->SetRangeUser(-50,50);
    hxy_zoomed_text->Draw("COLZ&&TEXT");

    TLine * lxmin = new TLine(-20,-20,-20,+20);
    TLine * lxmax = new TLine(20,-20,20,+20);
    TLine * lymin = new TLine(-20,-20,20,-20);
    TLine * lymax = new TLine(-20,20,20,20);

    lxmin->SetLineColor(kRed);
    lxmax->SetLineColor(kRed);
    lymin->SetLineColor(kRed);
    lymax->SetLineColor(kRed);

    lxmin->Draw("SAME");
    lxmax->Draw("SAME");
    lymin->Draw("SAME");
    lymax->Draw("SAME");
    //computing integral
    TAxis *xaxis = hxy_zoomed_text->GetXaxis();
    TAxis *yaxis = hxy_zoomed_text->GetYaxis();

    double muonsinarea = hxy_zoomed_text->Integral(
            xaxis->FindBin(-19.9),xaxis->FindBin(19.9),yaxis->FindBin(-19.9),yaxis->FindBin(19.9));
    
    cout<<"Integral in area: "<<muonsinarea<<endl;


}