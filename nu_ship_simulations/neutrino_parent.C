//reading old neutrino production (2017) with information about neutrino parent (AI: 30 January 2025)
bool ischarm(Float_t PdgCode){
  ROOT::RVec<int> charmpdglist = {421,411,431,4122,4232,4132,4332};
  return ROOT::VecOps::Any(charmpdglist==TMath::Abs(PdgCode));
      }

bool isnotcharm(Float_t PdgCode){
  ROOT::RVec<int> charmpdglist = {421,411,431,4122,4232,4132,4332};
  return !ROOT::VecOps::Any(charmpdglist==TMath::Abs(PdgCode));
      }

bool iskaon(Float_t PdgCode){
  ROOT::RVec<int> kaonpdglist = {321,130,310};
  return ROOT::VecOps::Any(kaonpdglist==TMath::Abs(PdgCode));
}

void neutrino_parent(){
 ROOT::RDataFrame df("pythia8-Geant4","pythia8_Geant4-withCharm_onlyNeutrinos.root"); // /home/utente/Simulations/pythia8_Geant4-withCharm_onlyNeutrinos.root
 //define momentum from px,py,pz
 auto df_p = df.Define("p","TMath::Sqrt(px*px+py*py+pz*pz)");
 //filtering to neutrino flavours

 auto df_nue = df_p.Filter("TMath::Abs(id)==12");
 auto df_numu = df_p.Filter("TMath::Abs(id)==14");
 auto df_nutau = df_p.Filter("TMath::Abs(id)==16");

 //parent selection
 auto df_numu_charmparent = df_numu.Filter(ischarm,{"parentid"}); //for charm production, parent is stored in parentid
 auto df_numu_notcharmparent = df_numu.Filter(isnotcharm,{"parentid"}); //for charm production, parent is stored in parentid

 auto df_nue_charmparent = df_nue.Filter(ischarm,{"parentid"}); //for charm production, parent is stored in parentid
 auto df_nue_notcharmparent = df_nue.Filter(isnotcharm,{"parentid"}); //for charm production, parent is stored in parentid

 auto df_nutau_charmparent = df_nutau.Filter(ischarm,{"parentid"}); //for charm production, parent is stored in parentid
 auto df_nutau_notcharmparent = df_nutau.Filter(isnotcharm,{"parentid"}); //for charm production, parent is stored in parentid

 //auto df_numu_pionparent = df_numu.Filter("TMath::Abs(pythiaid)==211"); //for not charm, parent is stored in pythiaid
 //auto df_numu_kaonparent = df_numu.Filter(iskaon,{"pythiaid"}); //for not charm, parent is stored in pythiaid

 auto hp_nue = df_nue.Histo1D(
  {"hp_nue","All;p[GeV/c]",400,0,400},"p","w");
 auto hp_charmnue = df_nue_charmparent.Histo1D(
  {"hp_charmnue","From charm decay;p[GeV/c]",400,0,400},"p","w");
 auto hp_notcharmnue = df_nue_notcharmparent.Histo1D(
  {"hp_notcharmnue","Not from charm decay;p[GeV/c]",400,0,400},"p","w");
 //normalize to unity (old target, absolute integrals can be misleading)
 hp_charmnue->Scale(1./hp_nue->Integral());
 hp_notcharmnue->Scale(1./hp_nue->Integral());

 auto hp_numu = df_numu.Histo1D(
  {"hp_numu","All;p[GeV/c]",400,0,400},"p","w");
 auto hp_charmnumu = df_numu_charmparent.Histo1D(
  {"hp_charmnumu","From charm decay;p[GeV/c]",400,0,400},"p","w");
 auto hp_notcharmnumu = df_numu_notcharmparent.Histo1D(
  {"hp_notcharmnumu","Not from charm decay;p[GeV/c]",400,0,400},"p","w");
 //normalize to unity (old target, absolute integrals can be misleading)
 hp_charmnumu->Scale(1./hp_numu->Integral());
 hp_notcharmnumu->Scale(1./hp_numu->Integral());

 auto hp_nutau = df_nutau.Histo1D(
  {"hp_nutau","All;p[GeV/c]",400,0,400},"p","w");
  auto hp_charmnutau = df_nutau_charmparent.Histo1D(
  {"hp_charmnutau","From charm decay;p[GeV/c]",400,0,400},"p","w");
 auto hp_notcharmnutau = df_nutau_notcharmparent.Histo1D(
  {"hp_notcharmnutau","Not from charm decay;p[GeV/c]",400,0,400},"p","w");
 //normalize to unity (old target, absolute integrals can be misleading)
 hp_charmnutau->Scale(1./hp_nutau->Integral());
 hp_notcharmnutau->Scale(1./hp_nutau->Integral());
 
gStyle->SetOptStat("nei");

 auto dtest = df_nutau_notcharmparent.Display({"parentid"}, 128);

 TCanvas *cp_nue = new TCanvas("cp_nue","Momentum electron neutrinos");
 hp_notcharmnue->DrawClone();
 hp_charmnue->SetLineColor(kRed);
 hp_charmnue->DrawClone("SAMES");
 cp_nue->SetLogy();
 cp_nue->BuildLegend();


 TCanvas *cp_numu = new TCanvas("cp_numu","Momentum muon neutrinos");
 hp_notcharmnumu->DrawClone();
 hp_charmnumu->SetLineColor(kRed);
 hp_charmnumu->DrawClone("SAMES");
 cp_numu->SetLogy();
 cp_numu->BuildLegend();

 TCanvas *cp_nutau = new TCanvas("cp_nutau","Momentum tau neutrinos");
 hp_charmnutau->SetLineColor(kRed);
 hp_charmnutau->DrawClone();
 hp_notcharmnutau->DrawClone("SAMES");
 cp_nutau->SetLogy();
 cp_nutau->BuildLegend();

 cout<<"Nutau not from charm, who is the parent? "<<endl;
 dtest->Print();
}