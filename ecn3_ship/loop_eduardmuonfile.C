//still cannot use rrdataframe in cbmsim to apply selections easily, better go back to TTreeReader loops (28 September 2022)
void create_tree(){
 //first, opening file, setting ttreereader
 //TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/snd_ship_ecn3/muons_3_planes.root");
 TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/ntuples_for_dis/combi_ecn3_fixed_overlaps_new_field/merged_combi_ecn3_spill_2.root");
 TTree *simtree = (TTree*)inputfile->Get("cbmsim");
 TTreeReader reader(simtree);

 //setting tracks and hits branches;
 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 //TTreeReaderArray<vetoPoint> sco1points(reader,"sco_1Point");
 TTreeReaderArray<UpstreamTaggerPoint> points(reader,"UpstreamTaggerPoint");

 const int nentries = reader.GetEntries();
 
 TFile *outputfile = new TFile("mcpoints_upstreamtagger_merged_combi_ecn3_spill_2.root","RECREATE");
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
 TFile *inputfile = TFile::Open("mcpoints_upstreamtagger_merged_combi_ecn3_spill_2.root");
 ROOT::RDataFrame df("mcpoints",inputfile);
 
 //derivative variables, (tri-momentum, energy, mass, charge...);
 auto df1 = df.Define("P","TMath::Sqrt(Px*Px+Py*Py+Pz*Pz)").Define("Charge",GetCharge,{"PdgCode"}); 
 
 //filling histograms
 auto hxy = df1.Histo2D({"hxy","xy distribution of muons;x[cm];y[cm]",140,-600,800,140,-600,800},"X","Y","Weight");
 auto hxy_zoomed = df1.Histo2D({"hxy_zoomed","xy distribution of muons;x[cm];y[cm]",400,-200,200,400,-200,200},"X","Y","Weight");
 auto hxy_zoomed_text = df1.Histo2D({"hxy_zoomed_text","xy distribution of muons;x[cm];y[cm]",40,-200,200,40,-200,200},"X","Y","Weight");

 auto hP = df1.Histo1D({"hp",";P[GeV/c]",60,0,300},"P","Weight");
 auto hPmap = df1.Profile2D({"hPmap","Momentum map of muons;x[cm];y[cm];P[GeV/c]",20,-100,100,20,-100,100,0,14000},"X","Y","P","Weight");
 
 auto hChargemap = df1.Profile2D({"hChargemap","Charge map of muons;x[cm];y[cm];Charge",20,-100,100,20,-100,100,-1,1},"X","Y","Charge","Weight");
 
 //drawing histograms and saving them to a ROOT plot file
 TFile *outputfile = new TFile("plots_negativemuons_upstreamBIG.root","RECREATE");
 TCanvas *cxy = new TCanvas("cxy","xy distribution",800,800);
 hxy->DrawClone("COLZ");
 cxy->Write();
 
 TCanvas *cxy_zoomed = new TCanvas("cxy_zoomed","xy distribution",800,800);
 hxy_zoomed->DrawClone("COLZ");
 cxy_zoomed->Write();
 TCanvas *cxy_zoomed_text = new TCanvas("cxy_zoomed_text","xy distribution with text in bins",800,800);
 hxy_zoomed_text->DrawClone("COLZ&&TEXT");    
 cxy_zoomed_text->Write();
 
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
