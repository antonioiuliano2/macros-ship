//macro per l'input FLUKA, semplificare studiofondo per ridurlo a ottenere l'energia cinetica dei primari (creato il 9 Maggio 2018 :( )

#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
//#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/latest/include/FairMCEventHeader.h"
bool isunknown(Int_t PdgCode);
bool isquark(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);
void primari(){
  
  //TFile *inputfile = TFile::Open("/home/utente/Simulations/sim_charmdet/pot/Lead6ECCs_nofieldnogaps_G4only_uniform/pythia8_Geant4_1_0.001.root");
  TFile *inputfile = TFile::Open("/home/utente/Simulations/sim_charmdet/pot/charm2_pot_distributions/pythia8_evtgen_Geant4_1000_0.001.root");
  //TFile *inputfile = TFile::Open("ship.10.0.Pythia8_NoCharm-TGeant4_centralbeam_10000.root");
  //Opnening sim tree and getting desired branches
  TTree* cbmsim = (TTree*)inputfile->Get("cbmsim");
  const Int_t neventi = cbmsim->GetEntries();

  Double_t pcut = 0.1; //cut for visibility in emulsion
  Double_t Tmax = 0.6; //both for TX and TY separately
  Double_t nmin = 2; //almost 2 tracks for a vertex

  TClonesArray *arr0 = new TClonesArray("ShipMCTrack",1000);
  cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
  cbmsim->SetBranchAddress("MCTrack",&arr0);
  cbmsim->SetBranchStatus("*",1);

  TDatabasePDG *pdg = TDatabasePDG::Instance();
  pdg->AddParticle("F0P0    ", " ", 0.9960, kFALSE, 0.0, 0, "meson",  9010221);//dalla pagina di AliRoot
  const Double_t kAu2Gev=0.9314943228;
  const Double_t khSlash = 1.0545726663e-27;
  const Double_t kErg2Gev = 1/1.6021773349e-3;
  const Double_t khShGev = khSlash*kErg2Gev;
  const Double_t kYear2Sec = 3600*24*365.25;
  pdg->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,0,3,"Ion",1000010020); //dal GitHub di FairRootGroup
  pdg->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,0,6,"Ion",1000020030); 
  pdg->AddParticle("Triton","Triton",3*kAu2Gev+14.931e-3,kFALSE, khShGev/(12.33*kYear2Sec),3,"Ion",1000010030);

  //histograms setting
  TH1D *hEkin_pion = new TH1D("hEkin_pion", "Energia pioni primari",400,0,400);
  TH1D *hmomentum = new TH1D("hmomentum_proton", "Impulso figli di protone", 400, 0, 400);
  TH1I *hmolteplicity = new TH1I("hmolteplicity", "Molteplicity", 40, 0, 40);
  TH1D *htheta = new TH1D("htheta", "Angle Theta", 3500, 0, 3.5);
  TH1D *hvz = new TH1D("hvz","Z position", 400, -40000, 0);

  //ofstream outputfile("./studiofondo/histopionkinenergy_lead.txt", std::ofstream::out);
  ofstream daughterfile("temp_outputs/proton_daughters.txt",std::ofstream::out);
  for (int i  = 0; i < neventi; i++){
    arr0->Clear();
    cbmsim->GetEntry(i);
    Double_t Vz = ((ShipMCTrack*) arr0->At(0))->GetStartZ();
    Int_t molteplicity = 0;
    hvz->Fill((Vz - 124.94)*10000); //in order to compare with MC
    for (int j = 2; j < arr0->GetEntriesFast();j++){
      ShipMCTrack* track = (ShipMCTrack*) arr0->At(j);
      Int_t mumID = track->GetMotherId();
      Int_t pdgcode = track->GetPdgCode();
      Double_t momentum = track->GetP();
      Double_t mass = track->GetMass();
      Double_t TX = track->GetPx()/track->GetPz();
      Double_t TY = track->GetPy()/track->GetPz();      
      Double_t theta = TMath::ATan(TMath::Sqrt(TX * TY + TY * TY));

      //escludo le particelle non rivelabili
      if (!isquark(pdgcode) && !isintermediate(pdgcode) && !isunknown(pdgcode)){
     // cout<<pdgcode<<endl;
      //if (mumID == 0){ //figlie del protone iniziale
       if ((track->GetStartZ() == Vz)){ //exactly started at Vz
       //if ((track->GetStartZ() >= Vz) && (track->GetStartZ() <= (Vz + 0.0200))){ //check in a range
       
        if (TMath::Abs(pdg->GetParticle(pdgcode)->Charge()) > 0){
         molteplicity++;

         if ((track->GetP() > pcut) && (TX < Tmax) && (TY < Tmax)){
          hmomentum->Fill(track->GetP()); //CUTS ON MOMENTUM AND TRACK ANGLE
          htheta->Fill(theta);
          }
         daughterfile<<i<<" "<<pdg->GetParticle(pdgcode)->GetName()<<" "<<track->GetPdgCode()<<" "<<track->GetP()<<" "<<track->GetStartZ()<<endl;
        }        
	
        Double_t kinen = pow(pow(momentum,2) + pow(mass,2),0.5) - mass ;
	if (abs(track->GetPdgCode()) == 211) hEkin_pion->Fill(kinen);
      } 
     }      
    }    
  daughterfile<<endl;
  if (molteplicity >= nmin) hmolteplicity->Fill(molteplicity);
  }//fine loop sugli eventi
  hEkin_pion->Draw();
  hEkin_pion->GetXaxis()->SetTitle("GeV");
  hEkin_pion->Scale(1./hEkin_pion->Integral());
  TCanvas *c = new TCanvas();
  hmolteplicity->Draw();
  TCanvas *c2 = new TCanvas();
  hmomentum->Draw();
  TCanvas *c3 = new TCanvas();
  htheta->Draw();
  TCanvas *c4 = new TCanvas();
  hvz->Draw();
  //for (int k = 1; k < (hEkin_pion->GetNbinsX() + 1); k++){
   //outputfile<<hEkin_pion->GetXaxis()->GetBinUpEdge(k)<<" "<<hEkin_pion->GetBinContent(k)<<endl;
// }
  //outputfile.close();
  daughterfile.close();
}

bool isunknown(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 9902110) || (abs(PdgCode) == 9902210) || (abs(PdgCode) == 990) || (abs(PdgCode) == 1000060120) || (abs(PdgCode) == 1000350810)) check = true;
  return check;
}

bool isquark(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 1) || (abs(PdgCode) == 2) || (abs(PdgCode) == 21) || (abs(PdgCode) == 2103) || (abs(PdgCode) == 2101) || (abs(PdgCode) == 2203)) check = true;
  return check;
}

bool isintermediate(Int_t PdgCode){ //risonanze forti a vita molto breve
  bool check = false;
  if ((abs(PdgCode) == 223) || (abs(PdgCode) == 3224) || (abs(PdgCode) == 331) || (abs(PdgCode) == 221) || (abs(PdgCode) == 20213) || (abs(PdgCode) == 3212) ||(abs(PdgCode)==213) || (abs(PdgCode) == 113) || (abs(PdgCode) == 2224) || (abs(PdgCode) == 323)) check =  true;
  return check;
  }
