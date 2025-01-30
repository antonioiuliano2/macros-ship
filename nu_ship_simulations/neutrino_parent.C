//reading old neutrino production (2017) with information about neutrino parent (AI: 30 January 2025)
bool ischarm(Float_t PdgCode){
  ROOT::RVec<int> charmpdglist = {421,411,431,4122,4232,4132,4332};
      return ROOT::VecOps::Any(charmpdglist==TMath::Abs(PdgCode));
      }

void neutrino_parent(){
 ROOT::RDataFrame df("pythia8-Geant4","pythia8_Geant4-withCharm_onlyNeutrinos.root"); // /home/utente/Simulations/pythia8_Geant4-withCharm_onlyNeutrinos.root
 //define momentum from px,py,pz
 auto df_p = df.Define("p","TMath::Sqrt(px*px+py*py+pz*pz)");
 //filtering to neutrino flavours
 auto df_nue = df_p.Filter("TMath::Abs(id)==12");
 auto df_numu = df_p.Filter("TMath::Abs(id)==14");
 auto df_nutau = df_p.Filter("TMath::Abs(id)==16");

 auto df_numu_charmparent = df_numu.Filter(ischarm,{"parentid"});
 auto hp_charmnumu = df_numu_charmparent.Histo1D({"hp","Momentum of muons at end of hadron stopper;p[GeV/c]",400,0,400},"p","w");

 TCanvas *cp = new TCanvas("cp","Momentum");
 hp_charmnumu->DrawClone();
}