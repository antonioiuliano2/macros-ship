#include <set>

#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
//#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

void studiofondo();
void conteggioelastici();
bool isintermediate(Int_t PdgCode);
void stimafondo();

void studiofondo(){
  stimafondo();
  //conteggioelastici();
}

void stimafondo(){
  ofstream outputfile("./studiofondo/histopionkinenergy_prova.txt", std::ofstream::out);
  
  //studio sul numero di tracce di noise nei film di emulsione.
  //TFile *file1 = TFile::Open("ship.10.0.Pythia8_NoCharm-TGeant4_inclusive.root");
  TFile *file1 = TFile::Open("ship.10.0.Pythia8_NoCharm-TGeant4_centralbeam_10000.root");
  TTree* cbmsim = (TTree*)file1->Get("cbmsim");

  TRandom3 *randomgen = new TRandom3;
 
  TFile *probfile = TFile::Open("./studiofondo/Out_Charm_1m.root");
  TH1D * hprob = (TH1D*) probfile->Get("PART/PRIMARY/hmyE_passed");
  
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  pdg->AddParticle("F0P0    ", " ", 0.9960, kFALSE, 0.0, 0, "meson",  9010221);//dalla pagina di AliRoot
  
  FairMCEventHeader *primario;
  cbmsim->GetBranch("MCEventHeader.");
  cbmsim->SetBranchAddress("MCEventHeader.", &primario);
  
  TClonesArray *arr0 = new TClonesArray("ShipMCTrack",1000);
  cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
  cbmsim->SetBranchAddress("MCTrack",&arr0);
  
  cbmsim->SetBranchStatus("*",1);
  
  TH1I *hnprimari = new TH1I("hnprimari", "Numero particelle al primario", 10, 0, 30);
  TH1D *hl = new TH1D("hl", "Lunghezza di volo primari", 40, 0, 100);
  TH1D *hp1 = new TH1D("hp1", "Impulso primari",100,0,110);
  TH1D *hE1 = new TH1D("hE1", "Energia pioni primari",400,0,400);
  TH1D *hEkin_pion = new TH1D("hEkin_pion", "Energia pioni primari",400,0,400);
  TH1D *hp2 = new TH1D("hp2", "Impulso figlie di primari", 100, 0, 10);
  TH1D *hp_T = new TH1D("hp_T", "Impulso trasverso", 25, 0, 2.5);
  TH1D *hkink = new TH1D("hkink", "Angolo kink", 100, 0, 2);

  //istogrammi spettri differenziati per particella
  TH1D *hprotonspectrum = new TH1D("hprotonspectrum", "Proton momentum", 4000, 0, 400);
  TH1D *hneutronspectrum = new TH1D("hneutronspectrum", "Neutron momentum", 4000, 0, 400);
  TH1D *hpionspectrum = new TH1D("hpionspectrum","Pion momentum", 4000, 0, 400);
  TH1D *hkaonspectrum = new TH1D("hkaonspectrum","Kaon momentum", 4000, 0, 400);
  TH1D *hpion0spectrum = new TH1D("hpion0spectrum","Neutral Pion momentum", 4000, 0, 400);
  TH1D *hkaon0spectrum = new TH1D("hkaon0spectrum","Neutral Kaon momentum", 4000, 0, 400);
  
  Int_t mumID;
  Int_t pdg1, pdg2;
  Double_t px1, py1, pz1;
  Double_t px2, py2, pz2;
  Double_t differenzax, differenzay;
  
  Double_t length,kink, impulso, energiacin;
  Double_t Vz;
  
  bool initialscattering;
  Int_t nkink = 0; //numero figli di primari che interagiscono in 6 mm;
  Int_t nshort = 0; //numero figli di primari che interagiscono in 6 mm;
  
  Int_t nevent = 0; //numero figli di primari che interagiscono in 6 mm in un singolo evento;
  Int_t nlength = 0;
  Int_t nprimari = 0; //numero particelle al primario;
  Int_t nbackgroundtracks = 0;
  Int_t nbackgroundevents = 0;
  
  set<Int_t> primari;
  set<Int_t> primari_visti;
  const Int_t neventi = cbmsim->GetEntries();
  for (int i = 0; i < neventi; i++){
  if (i%100 == 0) cout<<i<<endl;
  cbmsim->GetEntry(i);
  Vz = primario->GetZ();

  nbackgroundtracks = 0;
  
  nevent = 0;
  
  initialscattering = false;

   
   ShipMCTrack *check = (ShipMCTrack*) arr0->At(2);
   if (check->GetPz() > 390.){ //individuo eventuali scattering elastici, nel caso individuo i figli della successiva interazione di protone come i figli dell'interazione della particella 2
  initialscattering = true;
     }
   for (int j = 0; j < arr0->GetEntries(); j++){
     ShipMCTrack *trk2 = (ShipMCTrack*) arr0->At(j); 
     mumID = trk2->GetMotherId();
     pdg2 = trk2->GetPdgCode();
     if (j < 2) continue; //le prime due tracce sono il protone primario e la traccia di sistema
     if(pdg2 > 100000) continue;
     if (pdg2 == 990) continue;    
     px2 = trk2->GetPx();
     py2 = trk2->GetPy();
     pz2 = trk2->GetPz();
     if (trk2->GetStartZ() == Vz){
     //  if ((trk2->GetStartZ() == Vz) || ((trk2->GetMotherId() == 2) && (initialscattering == true))){
     //   if (((trk2->GetStartZ() == Vz) || ((trk2->GetMotherId() == 2) && (initialscattering == true))) && (trk2->GetPz() < 390.)){ //individuo i primari
       //if((pdg2 != 9902210) && (abs(pdg2) != 3) && (abs(pdg2) != 2203) && (abs(pdg2) != 2101) && (abs(pdg2) != 2103) && (abs(pdg->GetParticle(pdg2)->Charge()) > 0) && (isintermediate(pdg2) == false )) { //bisogna escludere gli stati intermedi e 'anomali', come i diquark
         if((pdg2 != 9902210) && (abs(pdg2) != 3) && (abs(pdg2) != 2203) && (abs(pdg2) != 2101) && (abs(pdg2) != 2103) && (isintermediate(pdg2) == false )) { //inserisco anche le particelle neutre
	 nprimari++;
	 //cout<<pdg2<<endl;	 
	 primari.insert(j);
	 impulso = TMath::Sqrt(pow(trk2->GetPx(),2) + pow(trk2->GetPy(),2) + pow(trk2->GetPz(),2));

	 switch(abs(pdg2))
         {
         case (321): hkaonspectrum->Fill(impulso);
         break;
         case (211): hpionspectrum->Fill(impulso);
         break;
         case (2212): hprotonspectrum->Fill(impulso);
         break;
         case (2112): hneutronspectrum->Fill(impulso);
         break;        
         case (111): hpion0spectrum->Fill(impulso);
         break;
         case (311): hkaon0spectrum->Fill(impulso);
         break;
         //default: cout<<"Particella non identificata:"<<pdg2<<endl;
	 // break;
         }
	 energiacin = TMath::Sqrt(pow(impulso,2) + pow(pdg->GetParticle(pdg2)->Mass(),2)) - pdg->GetParticle(pdg2)->Mass();
	 if (abs(pdg2) == 211) hEkin_pion->Fill(energiacin);
	 Double_t prob = hprob->GetBinContent(hprob->FindBin(energiacin));
	 if (prob > randomgen->Uniform(0,1)) nbackgroundtracks++;
	 //cout<<pdg->GetParticle(pdg2)->GetName()<<" "<<pdg2<<endl;	  
       }
     }
     
     if (mumID > 1){
       ShipMCTrack *trk1 = (ShipMCTrack*) arr0->At(mumID); 
       pdg1 = trk1->GetPdgCode();
       px1 = trk1->GetPx();
       py1 = trk1->GetPy();
       pz1 = trk1->GetPz();
       //cout<<pdg2<<endl;
       if ((primari.find(mumID) != primari.end()) && (abs(pdg->GetParticle(pdg2)->Charge()) > 0)){ //ora pdg2 è il figlio del primario, pdg1 è il primario
       //if (trk1->GetStartZ() == Vz){

       //if (pz1 > 390.) cout<<"protone scatterato:"<<i<<" "<<mumID<<endl;
	 if (pz1 < 390.){ //evito protoni scatterati
	   length = trk2->GetStartZ() - trk1->GetStartZ();

	 
	   
	   if (length > 0.){ //evito stati intermedi
	     nlength++;
	        if (length < 0.6){
		  nshort++;
		}
	     if (primari_visti.find(mumID) == primari_visti.end()){ //la length la calcolo una sola volta
	       hl->Fill(length);
	  
	     }
 
	     differenzax = px2/pz2 - px1/pz1;
	     differenzay = py2/pz2 - py1/pz1;
	     kink = TMath::ATan(TMath::Sqrt(pow(differenzax,2) + pow(differenzay,2)));
	     
	     hkink->Fill(kink);
	     if (kink > 0.02){
	       nkink++;           
	     }

	    
	     if ((length < 0.6) && (kink > 0.02)) nevent++;

	     if (primari_visti.find(mumID) == primari_visti.end()){
	     impulso = TMath::Sqrt(pow(px1,2) + pow(py1,2) + pow(pz1,2));
	     hp1->Fill(impulso);
	     energiacin = TMath::Sqrt(pow(impulso,2) + pow(pdg->GetParticle(pdg1)->Mass(),2)) - pdg->GetParticle(pdg1)->Mass();
	     if (abs(pdg1) == 211) hE1->Fill(energiacin); //riempo con energia del primario
	       }
	     //prob = hprob->GetBinContent(hprob->FindBin(energiacin));
	     //if (prob > randomgen->Uniform(0,1)) nbackgroundtracks++;             
	     impulso = TMath::Sqrt(pow(px2,2) + pow(py2,2) + pow(pz2,2));
	     hp2->Fill(impulso);
	     hp_T->Fill(impulso * TMath::Sin(kink));
	      primari_visti.insert(mumID);
	     
	      }
	 }
	
       }//chiude l'if per l'associazione al primario
       
     }//chiude l'if sul mumID
     // if (i == 4) cout<<pdg2<<" "<<trk2->GetStartZ()<<" "<<Vz<<endl; 
   }//chiude il ciclo sulle tracce
   //cout<<primari_visti.size()<<" "<<primari.size()<<endl;
   hnprimari->Fill(primari_visti.size());
   if (nprimari == 0) cout<<"No primari: "<<i<<endl;
   //cout<<i<<" "<<nevent<<endl;

   //cout<<"Numero primari:"<<primari.size()<<" "<<i<<endl;
   primari.clear();
   primari_visti.clear();

   if (nbackgroundtracks >= 2) nbackgroundevents++;
 }//chiude il ciclo sugli eventi
 //TCanvas *c1 = new TCanvas(); 
 hnprimari->Draw();
 TCanvas *c2 = new TCanvas();
 hl->Draw();
 //cout<<nprimari<<endl;
 //cout<<"Frazione di particelle al primario che interagiscono in 6 mm: "<<(double) nshort/nlength<<endl;
 //cout<<"Frazione con angolo di kink maggiore di 0.02 rad: "<<(double) nkink/hkink->GetEntries()<<endl;

 cout<<"Numero eventi di fondo previsti: "<<nbackgroundevents<<endl;
 
 TCanvas *c3 = new TCanvas();
 hkink->Draw();
 TCanvas *c4 = new TCanvas();
 hp2->Draw();
 TCanvas *c5 = new TCanvas();
 hp1->Scale(1./hp1->Integral());
 hp1->Draw();
 hp1->SetTitle("");
 TCanvas *c6 = new TCanvas();
 hp_T->Draw();
 TCanvas *c7 = new TCanvas();
 hE1->Scale(1./hE1->Integral());
 hE1->Draw();
 TCanvas *c8 = new TCanvas();
 hEkin_pion->Scale(1./hEkin_pion->Integral());
 hEkin_pion->Draw();
 
 /*TCanvas *cspectra = new TCanvas();
 cspectra->Divide(2,2);
 cspectra->cd(1);
 hpionspectrum->Draw();
 hpionspectrum->GetXaxis()->SetTitle("GeV/c");
 cspectra->cd(2);
 hprotonspectrum->Draw();
 hprotonspectrum->GetXaxis()->SetTitle("GeV/c");
 cspectra->cd(3);
 hneutronspectrum->Draw();
 hneutronspectrum->GetXaxis()->SetTitle("GeV/c");
 cspectra->cd(4);
 hkaonspectrum->Draw();
 hkaonspectrum->GetXaxis()->SetTitle("GeV/c");*/
 
 TCanvas *cfitproton = new TCanvas();
 hprotonspectrum->Draw();
 hprotonspectrum->GetXaxis()->SetTitle("GeV/c");

 TCanvas *cfitneutron = new TCanvas();
 hneutronspectrum->Draw();
 hneutronspectrum->GetXaxis()->SetTitle("GeV/c");

 TCanvas *cfitkaon = new TCanvas();
 hkaonspectrum->Draw();
 hkaonspectrum->GetXaxis()->SetTitle("GeV/c");

 TCanvas *cfitpion = new TCanvas();
 hpionspectrum->Draw();
 hpionspectrum->GetXaxis()->SetTitle("GeV/c");
 
 TFile *histofile = new TFile("spectrum_histograms.root","RECREATE");
 hprotonspectrum->Write();
 hpionspectrum->Write();
 hkaonspectrum->Write();
 hneutronspectrum->Write();
 histofile->Close();
 

 TCanvas *cspectra2 = new TCanvas();
 cspectra2->Divide(2,1);
 cspectra2->cd(1);
 hpion0spectrum->Draw();
 hpion0spectrum->GetXaxis()->SetTitle("GeV/c");
 cspectra2->cd(2);
 hkaon0spectrum->Draw();
 hkaon0spectrum->GetXaxis()->SetTitle("GeV/c");

 for (int k = 0; k < (hE1->GetNbinsX() + 1); k++){
   outputfile<<hE1->GetXaxis()->GetBinUpEdge(k)<<" "<<hE1->GetBinContent(k)<<endl;
 }
 outputfile.close();
}


bool isintermediate(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 223) || (abs(PdgCode) == 3222) || (abs(PdgCode) == 3224) || (abs(PdgCode) == 331) || (abs(PdgCode) == 221) || (abs(PdgCode) == 20213) || (abs(PdgCode) == 3212) ||(abs(PdgCode)==213) || (abs(PdgCode) == 113) || (abs(PdgCode) == 2224)|| (abs(PdgCode) == 323) || (abs(PdgCode) == 3122)) check =  true;
  return check;
  }

void conteggioelastici(){
  //studio sulla frazione di interazioni elastiche
  TFile *file1 = TFile::Open("ship.10.0.Pythia8_NoCharm-TGeant4_TTmedium.root"); 
  TTree* cbmsim = (TTree*)file1->Get("cbmsim");

  TRandom3 *randomgen = new TRandom3;
 
  TFile *probfile = TFile::Open("./studiofondo/Out_Charm_1m.root");
  TH1D * hprob = (TH1D*) probfile->Get("PART/PRIMARY/hmyE_passed");
  
  TDatabasePDG *pdg = TDatabasePDG::Instance();
  pdg->AddParticle("F0P0    ", " ", 0.9960, kFALSE, 0.0, 0, "meson",  9010221);//dalla pagina di AliRoot
  
  FairMCEventHeader *primario;
  cbmsim->GetBranch("MCEventHeader.");
  cbmsim->SetBranchAddress("MCEventHeader.", &primario);
  
  TClonesArray *arr0 = new TClonesArray("ShipMCTrack",1000);
  cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
  cbmsim->SetBranchAddress("MCTrack",&arr0);
  
  cbmsim->SetBranchStatus("*",1);
  
  Int_t mumID;
  Int_t pdgcode;
  Double_t Vz;

  Int_t nelastici = 0;
  Int_t nelasticiinteragentifuori = 0; //interazioni elastiche, in cui i protoni interagiscono anelasticamente fuori dal bersaglio (dove sta il Muon Tagger)
  Int_t ntotali = 0;
  
  const Int_t neventi = cbmsim->GetEntries();
  cout<<"Numero eventi: "<<neventi<<endl; 
  for (int i = 0; i < neventi; i++){
    cbmsim->GetEntry(i);
    ShipMCTrack *trk = (ShipMCTrack*) arr0->At(2); //prima traccia dell'evento (0 è il sistema, 1 è il protone padre)
    //cout<<arr0->GetEntries()<<" "<<trk->GetPz()<<endl;
    Vz = primario->GetZ();
    if (Vz < 0){
      ntotali++;
      if ((trk->GetPdgCode() == 2212) && (trk->GetPz() > 399.)){
	if (arr0->GetEntries() < 4) nelastici++; //un evento con solo sistema, protone incidente e protone diffuso è sicuramente elastico
	else{
	  ShipMCTrack *trk = (ShipMCTrack*) arr0->At(3); //passo alla traccia successiva, vedo se viene prodotta dopo
	  if (trk->GetStartZ() > Vz) nelastici++;
	  if (trk->GetStartZ() > 0 && trk->GetStartZ() < 750){
	    cout<<i<<" "<<trk->GetPdgCode()<<" "<<trk->GetStartZ()<<" "<<trk->GetMotherId()<<endl;
	    nelasticiinteragentifuori++;
	  }
	}
      }
    }
  }
    cout<<"Numero interazioni elastiche generate"<<nelastici<<"su "<<ntotali<<endl;
    cout<<"Interazioni elastiche, con successiva interazione anelastica fuori"<<nelasticiinteragentifuori<<endl;
}
