//macro per l'input FLUKA, semplificare studiofondo per ridurlo a ottenere l'energia cinetica dei primari (creato il 9 Maggio 2018 :( )

#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
//#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

void primari(){
  
  //TFile *inputfile = TFile::Open("/home/utente/Simulations/sim_charmdet/pot/Lead6ECCs_nofieldnogaps_G4only_uniform/pythia8_Geant4_1_0.001.root");
  TFile *inputfile = TFile::Open("ship.10.0.Pythia8_NoCharm-TGeant4_centralbeam_10000.root");
  //Opnening sim tree and getting desired branches
  TTree* cbmsim = (TTree*)inputfile->Get("cbmsim");
  const Int_t neventi = cbmsim->GetEntries();

  TClonesArray *arr0 = new TClonesArray("ShipMCTrack",1000);
  cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
  cbmsim->SetBranchAddress("MCTrack",&arr0);


  cbmsim->SetBranchStatus("*",1);

  //histograms setting
  TH1D *hEkin_pion = new TH1D("hEkin_pion", "Energia pioni primari",400,0,400);

  ofstream outputfile("./studiofondo/histopionkinenergy_lead.txt", std::ofstream::out);
  
  for (int i  = 0; i < neventi; i++){
    arr0->Clear();
    cbmsim->GetEntry(i);
    Double_t Vz = ((ShipMCTrack*) arr0->At(0))->GetStartZ();
    for (int j = 0; j < arr0->GetEntriesFast();j++){
      ShipMCTrack* track = (ShipMCTrack*) arr0->At(j);
      Int_t mumID = track->GetMotherId();
      Double_t momentum = track->GetP();
      Double_t mass = track->GetMass();
      
      //if (mumID == 0){ //figlie del protone iniziale
      if (track->GetStartZ() == Vz){
        cout<<i<<" "<<Vz<<" "<<j<<" "<<track->GetPdgCode()<<" "<<track->GetMotherId()<<endl;
	Double_t kinen = pow(pow(momentum,2) + pow(mass,2),0.5) - mass ;
	if (abs(track->GetPdgCode()) == 211) hEkin_pion->Fill(kinen);
      }      
    }    
  cout<<endl;
  }//fine loop sugli eventi
  hEkin_pion->Draw();
  hEkin_pion->GetXaxis()->SetTitle("GeV");
  hEkin_pion->Scale(1./hEkin_pion->Integral());
  for (int k = 1; k < (hEkin_pion->GetNbinsX() + 1); k++){
   outputfile<<hEkin_pion->GetXaxis()->GetBinUpEdge(k)<<" "<<hEkin_pion->GetBinContent(k)<<endl;
 }
  outputfile.close();
}
