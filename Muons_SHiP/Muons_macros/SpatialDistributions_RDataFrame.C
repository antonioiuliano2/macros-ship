//once again, with speed
using namespace ROOT;
void SpatialDistributions_RDataFrame(){

 //chain preparation
 TString prefix("root:://eosuser.cern.ch/");//for ROOTXD
 
 TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
 TChain *simchain = new TChain("muon_mcpoints");

 const int nfiles = 67; //number of input files
 cout<<"reading file "<<endl;
 for (int i = 0; i < nfiles; i++){
    simchain->Add((prefix+TString(Form("/eos/user/a/aiuliano/public/sims_FairShip/sim_muon_background/vetopoints_muons2018/vetopoints_scoringplane_muons_withCharmandBeauty%i.root",i*1000))).Data());
 }
 const int nentries = simchain->GetEntries();
 cout<<"Number of events"<<nentries<<endl;

 //build dataframe
 auto df = RDataFrame(*simchain);
 
 auto df1 = df.Define("p","TMath::Sqrt(px * px + py * py + pz * pz)");
 //make plots
 auto hxy = df1.Histo2D({"hxy","XY of muons at end of hadron stopper;x[cm];y[cm]",1000,-500,500,1000,-500,500},"x","y","weight");
 auto hp = df1.Histo1D({"hp","Momentum of muons at end of hadron stopper;p[GeV/c]",400,0,400},"p","weight");
 
 //draw histograms
 TCanvas *cxy = new TCanvas("cxy","XY Distribution",800,800);
 hxy->DrawClone("COLZ");

 TCanvas *cp = new TCanvas("cp", "Momentum distribution");
 hp->DrawClone();
 
}
