//still cannot use rrdataframe in cbmsim to apply selections easily, better go back to TTreeReader loops (28 September 2022)
void create_tree(){
 //first, opening file, setting ttreereader
 //TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/snd_ship_ecn3/muons_3_planes.root");
 //TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/ntuples_for_dis/combi_ecn3_fixed_overlaps_new_field/merged_combi_ecn3_spill_2.root");
 //TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/ntuples_for_dis/combi_ecn3_fixed_snd_field_off/merge_combi_ecn3_spill_snd_scor_planes_fixed_fields_11102022.root");
 TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/ntuples_for_dis/optimized_18102022/merged_optimized_18102022.root");
 TTree *simtree = (TTree*)inputfile->Get("cbmsim");
 TTreeReader reader(simtree);

 //setting tracks and hits branches;
 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 //TTreeReaderArray<vetoPoint> sco1points(reader,"sco_1Point");
 TTreeReaderArray<UpstreamTaggerPoint> points(reader,"UpstreamTaggerPoint");

 const int nentries = reader.GetEntries();
 
 TFile *outputfile = new TFile("mcpoints_upstreamtagger_optimized_18102022.root","RECREATE");
 TNtuple *simplepointstree = new TNtuple("mcpoints","MCPoints from UpstreamTagger","MCEventID:MCTrackID:MCMotherId:PdgCode:Px:Py:Pz:X:Y:Z:Weight");
 //const int nentries = 100000;

 cout<<"Starting looping over entries "<<nentries<<endl;

 for(int ientry = 0;ientry<nentries;ientry++){    
    reader.SetEntry(ientry);
    //looping over hits
    for (const UpstreamTaggerPoint& point: points){
        //accessing track to get weight
        int trackID = point.GetTrackID();
        int motherId = tracks[trackID].GetMotherId();

        //separating positive and negative muons
        int pdgcode = point.PdgCode();
        
        double px = point.GetPx();
        double py = point.GetPy();
        double pz = point.GetPz();
        
        double x = point.GetX();
        double y = point.GetY();
        double z = point.GetZ();
        
        double weight = tracks[trackID].GetWeight();
        
        simplepointstree->Fill(ientry,trackID,motherId,pdgcode,px,py,pz,x,y,z,weight);

    }
 }
 outputfile->cd();
 simplepointstree->Write();
 outputfile->Close(); 
}

double GetCharge(float PdgCode){
 TDatabasePDG *pdgdata = TDatabasePDG::Instance();
 double charge = (pdgdata->GetParticle(PdgCode)->Charge())/3.; //pdgdatabase returns charge in quark units (i.e. e- is -3.)
 return charge;
}

//reading produced ntuple with rdataframe
void read_tree(){
 //getting file and dataframe
 TFile *inputfile = TFile::Open("mcpoints_upstreamtagger_optimized_18102022.root");
 ROOT::RDataFrame df("mcpoints",inputfile);
 
 //derivative variables, (tri-momentum, energy, mass, charge...);

 auto df0 = df.Define("P","TMath::Sqrt(Px*Px+Py*Py+Pz*Pz)").Define("Theta","TMath::ACos(Pz/P)").Define("Charge",GetCharge,{"PdgCode"}); 
 //auto df1 = df0.Filter("Theta < 0.3");
 double xmin = -20.;
 double xmax = 20.;
 double ymin = -20.;
 double ymax = 20.;
 //auto df1 = df0.Filter("1");
 auto df1 = df0.Filter(Form("X>=%f && X<=%f && Y>=%f && Y<=%f",xmin,xmax,ymin,ymax));
 //filling histograms
 auto hxy = df1.Histo2D({"hxy","xy distribution of muons;x[cm];y[cm]",140,-600,800,140,-600,800},"X","Y","Weight");
 auto hxy_zoomed = df1.Histo2D({"hxy_zoomed","xy distribution of muons;x[cm];y[cm]",400,-200,200,400,-200,200},"X","Y","Weight");
 auto hxy_zoomed_text = df1.Histo2D({"hxy_zoomed_text","xy distribution of muons;x[cm];y[cm]",40,-200,200,40,-200,200},"X","Y","Weight");

 auto hP = df1.Histo1D({"hp",";P[GeV/c]",60,0,300},"P","Weight");
 auto hPmap = df1.Profile2D({"hPmap","Momentum map of muons;x[cm];y[cm];P[GeV/c]",20,-100,100,20,-100,100,0,14000},"X","Y","P","Weight");
 
 auto hChargemap = df1.Profile2D({"hChargemap","Charge map of muons;x[cm];y[cm];Charge",20,-100,100,20,-100,100,-1,1},"X","Y","Charge","Weight");
 
 auto hTheta = df1.Histo1D({"hTheta",";#theta[rad]",60,0,3.},"Theta","Weight");
 auto hThetaStera = df1.Histo1D({"hThetaStera",";#theta[sr];tracks/cm^{2}",40,0,1.},"Theta","Weight");
 auto hThetamap = df1.Profile2D({"hThetamap","Theta map of muons;x[cm];y[cm];#theta[rad]",20,-100,100,20,-100,100,0,3.},"X","Y","Theta","Weight");
 
 //drawing histograms and saving them to a ROOT plot file
 TFile *outputfile = new TFile("plots_muonshitseduard_upstreamBIG_optimized_onlyemuarea.root","RECREATE");
 TCanvas *cxy = new TCanvas("cxy","xy distribution",800,800);
 hxy->DrawClone("COLZ");
 cxy->Write();
 
 TCanvas *cxy_zoomed = new TCanvas("cxy_zoomed","xy distribution",800,800);
 hxy_zoomed->DrawClone("COLZ");
 cxy_zoomed->Write();
 TCanvas *cxy_zoomed_text = new TCanvas("cxy_zoomed_text","xy distribution with text in bins",800,800);
 hxy_zoomed_text->DrawClone("COLZ&&TEXT");    
 cxy_zoomed_text->Write();
 
 TCanvas *cTheta = new TCanvas("cTheta","Theta angle");
 hTheta->DrawClone();
 cTheta->Write();

 TCanvas *cThetaStera = new TCanvas("cThetaStera","Theta angle in steradians");
 //theta in steradians
 const int nbins = hThetaStera->GetNbinsX();
 double surface = (xmax - xmin) * (ymax - ymin);
 cout<<"Evaluate tracks/cm2 in area of "<<surface<<endl;
 for (int ibin = 1; ibin<=nbins;ibin++){
    double binvalue = hThetaStera->GetBinContent(ibin);
    double thetavalue = hThetaStera->GetBinCenter(ibin);
    double updatedbinvalue = binvalue/(1600*thetavalue*2*TMath::Pi()); //according to Komatsu Suggestion;
    hThetaStera->SetBinContent(ibin,updatedbinvalue);
 }
 hThetaStera->DrawClone("hist");
 hThetaStera->Write();

 TCanvas *cTheta_map = new TCanvas("cTheta_map","Angular map of muons",800,800);
 hThetamap->DrawClone("COLZ");
 cTheta_map->Write();
 
 
 TCanvas *cP = new TCanvas("cP","Momentum");
 hP->DrawClone();
 cP->SetLogy();
 cP->Write();

 TCanvas *cP_map = new TCanvas("cP_map","Momentum map of muons",800,800);
 hPmap->DrawClone("COLZ");
 cP_map->Write();

 TCanvas *cCharge_map = new TCanvas("cCharge_map","Charge map of muons",800,800);
 hChargemap->DrawClone("COLZ");
 cCharge_map->Write();
}

void plotmuondensityhistogram(){
    //TFile *histfile = TFile::Open("plots_muonshitseduard_upstreamBIG.root");
    TFile *histfile = TFile::Open("plots_muonshitseduard_upstreamBIG_lowangle.root");
    TCanvas *cxy = (TCanvas*) histfile->Get("cxy_zoomed_text");

    TH2D *hxy_zoomed_text = (TH2D*) cxy->GetPrimitive("hxy_zoomed_text");

    hxy_zoomed_text = (TH2D*) hxy_zoomed_text->RebinX(2);
    hxy_zoomed_text = (TH2D*) hxy_zoomed_text->RebinY(2);

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