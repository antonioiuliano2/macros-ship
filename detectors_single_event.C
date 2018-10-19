//studio di quanto avviene nello spettrometro magnetico e nell'RPC evento per evento (macro creata il 17 Marzo 2017)
//I display sulle emulsioni devono avere i vertici generati in xy in modo uniforme sul bersaglio (bordi a parte, vecchia sim. 8 x 8 cm quadri), ma quelli sullo spettrometro e il rivelatore per muoni no, devono averi i vertici generati in modo gaussiano con sigma 0.5 cm (bersaglio mobile, rivelatori a valle fissi!), tranne di fatto il primo spettrometro, quindi salvo solo quello
#include<vector>
#include "/home/utente/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
#include "/home/utente/SHIPBuild/FairShip/charmdet/MuonTaggerPoint.h"
#include "/home/utente/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
#include "/home/utente/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

void emulsion();
void spectrometer();
void RPC();
void RPCgraphs();
bool ischarm (Int_t PdgCode);
bool isintermediate(Int_t PdgCode);

void detectors_single_event(){
  //emulsion();
  spectrometer(); //disegna i display per i 4 piani dello spettrometro magnetico
  //RPC(); //disegna i display per gli RPC 
  //RPCgraphs(); //grafici sugli RPC (prima calcolava anche le distanze dei muoni dal centro)
}
void RPCgraphs(){
 TFile *file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_accettanzarivelatori.root");
 //TFile *file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_isolamentomuoni_10000.root");
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");

 const bool charm_event = true; //sto studiando gli eventi di charm o gli eventi di rumore ('fondo strumentale')?
 Int_t muonsfromcharm = 0;

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

 FairMCEventHeader *header;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.",&header);
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);

 TClonesArray *arr2 = new TClonesArray("SpectrometerPoint",1000);
 cbmsim->GetBranch("SpectrometerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("SpectrometerPoint", &arr2);
 
 TClonesArray *arr3 = new TClonesArray("MuonTaggerPoint",1000);
 cbmsim->GetBranch("MuonTaggerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MuonTaggerPoint", &arr3);

 cbmsim->SetBranchStatus("*",1);

 Double_t xPoint, yPoint, zPoint, impulso, angle;
 Double_t Vx, Vy, Vz; //coordinate del vertice primario
 Int_t nplane, pdgcode;
 
 const int nevents = cbmsim->GetEntries();

 TH1I *hcharmoccupancy = new TH1I("hcharmoccupancy", "Charm and other particles in same cell",5,0,5);

 TH2D *hspectropos1big = new TH2D("hspectropos1big", "Track position in first spectrometer station", 1000, -500, 500, 1000, -500, 500);
 TH2D *hspectropos2big = new TH2D("hspectropos2big", "Track position in second spectrometer station", 1000, -500, 500, 1000, -500, 500);
 TH2D *hspectropos3big = new TH2D("hspectropos3big", "Track position in third spectrometer station", 1000, -500, 500, 1000, -500, 500);
 TH2D *hspectropos4big = new TH2D("hspectropos4big", "Track position in fourth spectrometer station", 1000, -500, 500, 1000, -500, 500);
 
 TH2D *hspectropos1 = new TH2D("hspectropos1", "Track position in first spectrometer station", 100, -6.5, 6.5, 100, -5.5, 5.5);
 TH2D *hspectropos2 = new TH2D("hspectropos2", "Track position in second spectrometer station", 100, -10, 10, 100, -10, 10);
 TH2D *hspectropos3 = new TH2D("hspectropos3", "Track position in third spectrometer station", 200, -100, 100, 150, -75, 75); //avevo fatto a 130 cm per errore prima
 TH2D *hspectropos4 = new TH2D("hspectropos4", "Track position in fourth spectrometer station", 200, -100, 100, 150, -75, 75);


 TH2D * hrpcpos1 = new TH2D("hrpcpos1","Track position in first RPC", 200, -100, 100, 130, -65, 65);
 TH2D * hrpcpos2 = new TH2D("hrpcpos2","Track position in second RPC", 200, -100, 100, 130, -65, 65);
 TH2D * hrpcpos3 = new TH2D("hrpcpos3","Track position in third RPC", 200, -100, 100,130, -65, 65);
 TH2D * hrpcpos4 = new TH2D("hrpcpos4","Track position in fourth RPC", 200, -100, 100, 130, -65, 65);
 TH2D * hrpcpos5 = new TH2D("hrpcpos5","Track position in fifth RPC", 200, -100, 100, 130, -65, 65);
 TH2D * hrpcpos6 = new TH2D("hrpcpos6","Track position in sixth RPC", 200, -100, 100, 130, -65, 65);

 TH2D * hrpcpossingle[6];
 hrpcpossingle[0] = new TH2D("hrpcpos1single","Track position in first RPC", 200, -100, 100, 130, -65, 65);
 hrpcpossingle[1] = new TH2D("hrpcpos2single","Track position in second RPC", 200, -100, 100, 130, -65, 65);
 hrpcpossingle[2] = new TH2D("hrpcpos3single","Track position in third RPC", 200, -100, 100,130, -65, 65);
 hrpcpossingle[3] = new TH2D("hrpcpos4single","Track position in fourth RPC", 200, -100, 100, 130, -65, 65);
 hrpcpossingle[4] = new TH2D("hrpcpos5single","Track position in fifth RPC", 200, -100, 100, 130, -65, 65);
 hrpcpossingle[5] = new TH2D("hrpcpos6single","Track position in sixth RPC", 200, -100, 100, 130, -65, 65);
 /*
 TH2D * hrpcpos1single = new TH2D("hrpcpos1single","Track position in first RPC", 500, -100, 100, 130, -65, 65);
 TH2D * hrpcpos2single = new TH2D("hrpcpos2single","Track position in second RPC", 500, -100, 100, 130, -65, 65);
 TH2D * hrpcpos3single = new TH2D("hrpcpos3single","Track position in third RPC", 500, -100, 100,130, -65, 65);
 TH2D * hrpcpos4single = new TH2D("hrpcpos4single","Track position in fourth RPC", 500, -100, 100, 130, -65, 65);
 TH2D * hrpcpos5single = new TH2D("hrpcpos5single","Track position in fifth RPC", 500, -100, 100, 130, -65, 65);
 */
 
 TH2D * hrpcpos1big = new TH2D("hrpcpos1big","Track position in first RPC", 500, -250, 250, 500, -250, 250);
 TH2D * hrpcpos2big = new TH2D("hrpcpos2big","Track position in second RPC", 500, -250, 250, 500, -250, 250);
 TH2D * hrpcpos3big = new TH2D("hrpcpos3big","Track position in third RPC", 500, -250, 250, 500, -250, 250);
 TH2D * hrpcpos4big = new TH2D("hrpcpos4big","Track position in fourth RPC", 500, -250, 250, 500, -250, 250);
 TH2D * hrpcpos5big = new TH2D("hrpcpos5big","Track position in fifth RPC", 500, -250, 250, 500, -250, 250);
 TH2D * hrpcpos6big = new TH2D("hrpcpos6big","Track position in sixth RPC", 500, -250, 250, 500, -250, 250); 
 
 TH1D *hx1 = new TH1D("hx1", "X position first RPC", 200, -100, 100);
 TH1D *hy1 = new TH1D("hy1", "Y position first RPC", 130, -65, 65);
 TH1D *hx2 = new TH1D("hx2", "X position second RPC", 200, -100, 100);
 TH1D *hy2 = new TH1D("hy2", "Y position second RPC", 130, -65, 65);
 TH1D *hx3 = new TH1D("hx3", "X position third RPC", 200, -100, 100);
 TH1D *hy3 = new TH1D("hy3", "Y position third RPC", 130, -65, 65);
 TH1D *hx4 = new TH1D("hx4", "X position fourth RPC", 200, -100, 100);
 TH1D *hy4 = new TH1D("hy4", "Y position fourth RPC", 130, -65, 65);
 TH1D *hx5 = new TH1D("hx5", "X position fifth RPC", 200, -100, 100);
 TH1D *hy5 = new TH1D("hy5", "Y position fifth RPC", 130, -65, 65);
 TH1D *hx6 = new TH1D("hx6", "X position sixth RPC", 200, -100, 100);
 TH1D *hy6 = new TH1D("hy6", "Y position sixth RPC", 130, -65, 65);

 TH1I *hnrpcmuons = new TH1I ("hnrpcmuons", "Number of RPCs transversed by muons from charm decay", 6,0.5,6.5);
 TH1I *hnrpcisolated = new TH1I ("hnrpcisolated", "Number of RPCs where the muon is isolated both in x and y", 7,-0.5,6.5);
 
 Int_t mumID, mumpdg, intermediatepdg, trackID;

 bool muonfromcharm;
 const int nbuiltRPCs = 6; 
 Double_t xmuon[nbuiltRPCs], ymuon[nbuiltRPCs];
 Int_t nhitx[nbuiltRPCs], nhity[nbuiltRPCs]; //numero di hit adiacenti in x o y ai muoni
 
 Int_t nmuonsfromcharm = 0; //numero muoni totali figli di charm nell'evento
 Int_t nRPCs = 0; // numero di piani di RPC attraversati dalla traccia
 Int_t nRPCisolated;
 
 Double_t distanza; //distanza muoni dal centro

 Int_t ninterazioni; //numero interazioni nel target
  for (int i = 0; i < nevents; i++){
    //ripulisco gli istogrammi per il singolo evento.
    nmuonsfromcharm = 0;
    trackID = 0;
    muonfromcharm = false;
    for (int k = 0; k < nbuiltRPCs; k++){
      hrpcpossingle[k]->Reset();
      xmuon[k] = 0.;
      ymuon[k] = 0.;
    }
    
    
    if (i % 100 == 0) cout<<i<<endl;

    arr1->Clear();
    arr2->Clear();
    arr3->Clear();
    cbmsim->GetEntry(i);    

    ShipMCTrack *firsttrack = (ShipMCTrack*) arr1->At(0);
    Vx= firsttrack->GetStartX();
    Vy = firsttrack->GetStartY();
    Vz = firsttrack->GetStartZ();
    
    ninterazioni++;
    //Spectrometer Stations
    for (int j = 0; j < arr2->GetEntriesFast(); j++){
      SpectrometerPoint *spectro =(SpectrometerPoint*) arr2->At(j);
      xPoint = spectro->GetX();
      yPoint = spectro->GetY();
      zPoint = spectro->GetZ();
      pdgcode = spectro->PdgCode();
      
      bool checktrack = false; //voglio salvare la posizione di figlie di charm e tracce al vertice primario per lo studio dell'efficienza
      
      bool charge = false;
      if (pdgcode > 100000.) charge = true;
      else if (abs(pdg->GetParticle(pdgcode)->Charge()) > 0) charge = true;
      if (charge == false) continue;
      
      impulso = pow((pow(spectro->GetPx(),2) + pow(spectro->GetPy(),2) + pow(spectro->GetPz(),2)),0.5);
      angle = TMath::ACos(spectro->GetPz()/impulso);
      nplane = 0;
      if (zPoint < 0.2) nplane = 1; //quando ci sono più piani per stazione, considero il piano più a monte
      else if ((zPoint > 15) && (zPoint < 20.75)) nplane = 2;
      else if ((zPoint > 100.) &&  (zPoint < 490.)) nplane = 3; 
      else if (zPoint > 490.) nplane = 4;

      if (spectro->GetTrackID() > 0){
	ShipMCTrack * track = (ShipMCTrack*) arr1->At(spectro->GetTrackID()); //devo aver simulato con -F per poter leggere le tracce di tutti gli hit      
	Double_t startz = track->GetStartZ();
	if (startz == Vz) checktrack = true;
	mumID = track->GetMotherId();
	//controllo se è una figlia di charm
	if (mumID>0){
	  ShipMCTrack* trk = (ShipMCTrack*)arr1->At(mumID);
	  mumpdg = trk->GetPdgCode();
	  while ((isintermediate(mumpdg) == true) && (mumID>0)) //se la madre è uno stato intermedio controllo se lo stato intermedio è stato prodotto dal charm
	  {
	    intermediatepdg = mumpdg;
	    mumID = trk->GetMotherId();
	    if (mumID <=0 ) break; //evita crash se mumID <0
	    ShipMCTrack* trk = (ShipMCTrack*) arr1->At(mumID);	  
	    mumpdg = trk->GetPdgCode();
	  }     
	  if (ischarm(mumpdg) == true) checktrack = true;
	}
      }
      if (checktrack==true){
	if (nplane == 1){ //mi concentro sul primo RPC
	  if (checktrack==true) hspectropos1big->Fill(xPoint, yPoint); 
	  if (checktrack==true) hspectropos1->Fill(xPoint, yPoint); 
	}
	if (nplane == 2){
	  if (checktrack==true) hspectropos2big->Fill(xPoint, yPoint); 
	  if (checktrack==true) hspectropos2->Fill(xPoint,yPoint);
	}
	if (nplane == 3){
	  if (checktrack==true) hspectropos3big->Fill(xPoint, yPoint); 
	  if (checktrack==true) hspectropos3->Fill(xPoint,yPoint);
	}	 
	if (nplane == 4){
	  if (checktrack==true) hspectropos4big->Fill(xPoint, yPoint); 
	  if (checktrack==true) hspectropos4->Fill(xPoint,yPoint);
	}
      }
    } //fine del loop sugli hit nello spettrometro
    //Muon Tagger
    for (int j = 0; j < arr3->GetEntriesFast(); j++){
      MuonTaggerPoint *muon =(MuonTaggerPoint*) arr3->At(j);
      xPoint = muon->GetX();
      yPoint = muon->GetY();
      zPoint = muon->GetZ();
      pdgcode = muon->PdgCode();     

      if ((muon->GetTrackID() != trackID) || (j == arr3->GetEntries() - 1)){ //sfrutto il fatto che gli hit della stessa traccia sono contigui. Quando il trackID cambia, vuol dire che la traccia si è fermata. In alternativa, può essere l'ultima traccia nell'evento.
	if (muonfromcharm){
	  if (j == arr3->GetEntries() - 1) nRPCs += 1; //se è l'ultima, bisogna sommare un nRPCs, relativo all'ultimo hit attraversato.
	  cout<<"Piani attraversati dal muone figlio di charm "<<nRPCs<<endl;
	  if (nRPCs > 0) hnrpcmuons->Fill(nRPCs);
	}
	nRPCs = 0;
      }
      
      nRPCs++;      
      trackID = muon->GetTrackID();
      muonfromcharm = false;
      bool checktrack = false; //voglio salvare la posizione di figlie di charm e tracce al vertice primario per lo studio dell'efficienza      
     

      bool charge = false;
      if (pdgcode > 100000.) charge = true;
      else if (abs(pdg->GetParticle(pdgcode)->Charge()) > 0) charge = true;
      if (charge == false) continue;
      impulso = pow((pow(muon->GetPx(),2) + pow(muon->GetPy(),2) + pow(muon->GetPz(),2)),0.5);
      angle = TMath::ACos(muon->GetPz()/impulso);
      
     if (zPoint < 550.) nplane = 1;
     else if (zPoint < 650.) nplane = 2;
     else if (zPoint < 700.) nplane = 3;
     else if (zPoint < 750.) nplane = 4;
     else if (zPoint< 800.) nplane = 5;
     else nplane = 6;
      /*if (zPoint < 680.) nplane = 1;
      else if (zPoint < 800.) nplane = 2;
      else if (zPoint < 1000.) nplane = 3;
      else if (zPoint < 1050.) nplane = 4;
      else nplane = 5;*/

      if (muon->GetTrackID() > 0){
	ShipMCTrack * track = (ShipMCTrack*) arr1->At(muon->GetTrackID());
	Double_t startz = track->GetStartZ();
	if (startz == Vz) checktrack = true;
	Int_t mumID = track->GetMotherId();
	//controllo se è una figlia di charm
	if (mumID>0){
	  ShipMCTrack* trk = (ShipMCTrack*)arr1->At(mumID);
	  mumpdg = trk->GetPdgCode();
	  while ((isintermediate(mumpdg) == true) && (mumID>0)) //se la madre è uno stato intermedio controllo se lo stato intermedio è stato prodotto dal charm
	    {
	      intermediatepdg = mumpdg;
	      mumID = trk->GetMotherId();
	      if (mumID <=0 ) break; //evita crash se mumID <0
	      ShipMCTrack* trk = (ShipMCTrack*) arr1->At(mumID);	  
	      mumpdg = trk->GetPdgCode();
	    }     
	  if (ischarm(mumpdg) == true){
	    checktrack = true;	  
	    if (abs(pdgcode) == 13){
	      cout<<"Muone figlio di charm nell'evento "<<i<<" nell'RPC numero "<<nplane<<"id muone"<<muon->GetTrackID()<<endl;
	      nmuonsfromcharm++;
	      muonfromcharm = true;
	      xmuon[nplane-1] = xPoint;
	      ymuon[nplane-1] = yPoint;	    
	    }
	  }
	}
      }
      
      //  if (checktrack == true){
	if (nplane == 1){ //mi concentro sul primo RPC
	  hx1->Fill(xPoint);
	  hy1->Fill(yPoint);
	  if (checktrack == true) hrpcpos1->Fill(xPoint, yPoint);
	  //hrpcpos1single->Fill(xPoint, yPoint);
	  hrpcpossingle[0]->Fill(xPoint, yPoint);
	  if (checktrack == true) hrpcpos1big->Fill(xPoint, yPoint);
	}//fine richiesta sul primo piano
	if (nplane == 2){       
	  hx2->Fill(xPoint);
	  hy2->Fill(yPoint);
	  if (checktrack == true) hrpcpos2->Fill(xPoint,yPoint);
	  if (checktrack == true) hrpcpos2big->Fill(xPoint, yPoint);
	  hrpcpossingle[1]->Fill(xPoint, yPoint);
	}
	if (nplane == 3){
	  hx3->Fill(xPoint);
	  hy3->Fill(yPoint);
	  if (checktrack == true) hrpcpos3->Fill(xPoint,yPoint);
	  if (checktrack == true) hrpcpos3big->Fill(xPoint, yPoint);
	  hrpcpossingle[2]->Fill(xPoint, yPoint);
	}	 
	if (nplane == 4){
	  hx4->Fill(xPoint);
	  hy4->Fill(yPoint);
	  if (checktrack == true) hrpcpos4->Fill(xPoint,yPoint);
	  if (checktrack == true) hrpcpos4big->Fill(xPoint, yPoint);
	  hrpcpossingle[3]->Fill(xPoint, yPoint);
       }     
	if (nplane == 5){
	  hx5->Fill(xPoint);
	  hy5->Fill(yPoint);
	  if (checktrack == true) hrpcpos5->Fill(xPoint,yPoint);
	  if (checktrack == true) hrpcpos5big->Fill(xPoint, yPoint);
	  hrpcpossingle[4]->Fill(xPoint, yPoint);
	}
	if (nplane == 6){
	  hx6->Fill(xPoint);
	  hy6->Fill(yPoint);
	  if (checktrack == true) hrpcpos6->Fill(xPoint,yPoint);
	  if (checktrack == true) hrpcpos6big->Fill(xPoint, yPoint);
          hrpcpossingle[5]->Fill(xPoint, yPoint);
	  //hrpcpos6single->Fill(xPoint, yPoint);
	}
	// }
    }//fine ciclo sugli hit    

    //devo ora confrontare la posizione dei muoni con quella delle altre tracce.
    // cout<<"Evento "<<i<<" NUMERO HIT NEL PRIMO RPC "<<hrpcpossingle[0]->GetEntries()<<endl;
    //cout<<"Evento "<<i<<" NUMERO HIT NEL QUINTO RPC "<<hrpcpossingle[4]->GetEntries()<<endl;
    nRPCisolated = 0;
    for (int k = 0; k < nbuiltRPCs; k++){
      nhitx[k] = 0;
      nhity[k] = 0;
      if (xmuon[k] != 0.){	
	nhitx[k] = hrpcpossingle[k]->ProjectionY("prova", hrpcpossingle[k]->GetXaxis()->FindBin(xmuon[k])-1, hrpcpossingle[k]->GetXaxis()->FindBin(xmuon[k]) + 1)->GetEntries();
	nhity[k] = hrpcpossingle[k]->ProjectionX("prova", hrpcpossingle[k]->GetYaxis()->FindBin(ymuon[k])-1, hrpcpossingle[k]->GetYaxis()->FindBin(ymuon[k]) + 1)->GetEntries();

	if (nhitx[k] == 1 && nhity[k] == 1) nRPCisolated ++;
      }
    }
    if(xmuon[0] != 0.) hnrpcisolated->Fill(nRPCisolated);
  }//fine ciclo sugli eventi
  
  cout<<"NUMERO MUONI DA CHARM CHE ARRIVANO AGLI RPC:" <<nmuonsfromcharm<<endl;
  //grafici distribuzioni spaziali sulle stazioni dello spettrometro
  TCanvas *cspectropos1 = new TCanvas();
  cspectropos1->Divide(1,2);
  cspectropos1->cd(1);
  hspectropos1big->GetXaxis()->SetTitle("cm");
  hspectropos1big->GetYaxis()->SetTitle("cm");
  hspectropos1big->Draw("COLZ");
  cspectropos1->cd(2);
  hspectropos1->GetXaxis()->SetTitle("cm");
  hspectropos1->GetYaxis()->SetTitle("cm");
  hspectropos1->Draw("COLZ");
  cspectropos1->Print("studionoise/dimensioni_maggiori_rivelatori/T1_positions.png");
  cspectropos1->Print("studionoise/dimensioni_maggiori_rivelatori/T1_positions.root");  
  TCanvas *cspectropos2 = new TCanvas();
  cspectropos2->Divide(1,2);
  cspectropos2->cd(1);
  hspectropos2big->GetXaxis()->SetTitle("cm");
  hspectropos2big->GetYaxis()->SetTitle("cm");
  hspectropos2big->Draw("COLZ");
  cspectropos2->cd(2);
  hspectropos2->GetXaxis()->SetTitle("cm");
  hspectropos2->GetYaxis()->SetTitle("cm");
  hspectropos2->Draw("COLZ");
  cspectropos2->Print("studionoise/dimensioni_maggiori_rivelatori/T2_positions.png");
  cspectropos2->Print("studionoise/dimensioni_maggiori_rivelatori/T2_positions.root");
  TCanvas *cspectropos3 = new TCanvas();
  cspectropos3->Divide(1,2);
  cspectropos3->cd(1);
  hspectropos3big->GetXaxis()->SetTitle("cm");
  hspectropos3big->GetYaxis()->SetTitle("cm");
  hspectropos3big->Draw("COLZ");
  cspectropos3->cd(2);
  hspectropos3->GetXaxis()->SetTitle("cm");
  hspectropos3->GetYaxis()->SetTitle("cm");
  hspectropos3->Draw("COLZ");
  cspectropos3->Print("studionoise/dimensioni_maggiori_rivelatori/T3_positions.png");
  cspectropos3->Print("studionoise/dimensioni_maggiori_rivelatori/T3_positions.root");
  TCanvas *cspectropos4 = new TCanvas();
  cspectropos4->Divide(1,2);
  cspectropos4->cd(1);
  hspectropos4big->GetXaxis()->SetTitle("cm");
  hspectropos4big->GetYaxis()->SetTitle("cm");
  hspectropos4big->Draw("COLZ");
  hspectropos4->GetXaxis()->SetTitle("cm");
  hspectropos4->GetYaxis()->SetTitle("cm");
  cspectropos4->cd(2);
  hspectropos4->Draw("COLZ");
  cspectropos4->Print("studionoise/dimensioni_maggiori_rivelatori/T4_positions.png");
  cspectropos4->Print("studionoise/dimensioni_maggiori_rivelatori/T4_positions.root");  
  //grafici distribuzioni spaziali sugli RPC
  TCanvas *crpcpos1 = new TCanvas();
  crpcpos1->Divide(1,2);
  crpcpos1->cd(1);
  hrpcpos1big->GetXaxis()->SetTitle("cm");
  hrpcpos1big->GetYaxis()->SetTitle("cm");
  hrpcpos1big->Draw("COLZ");
  crpcpos1->cd(2);
  hrpcpos1->GetXaxis()->SetTitle("cm");
  hrpcpos1->GetYaxis()->SetTitle("cm");
  hrpcpos1->Draw("COLZ");
  crpcpos1->Print("studionoise/dimensioni_maggiori_rivelatori/RPC1_positions.png");
  crpcpos1->Print("studionoise/dimensioni_maggiori_rivelatori/RPC1_positions.root");
  TCanvas *crpcpos2 = new TCanvas();
  crpcpos2->Divide(1,2);
  crpcpos2->cd(1);
  hrpcpos2big->GetXaxis()->SetTitle("cm");
  hrpcpos2big->GetYaxis()->SetTitle("cm");
  hrpcpos2big->Draw("COLZ");
  crpcpos2->cd(2);
  hrpcpos2->GetXaxis()->SetTitle("cm");
  hrpcpos2->GetYaxis()->SetTitle("cm");
  hrpcpos2->Draw("COLZ");
  crpcpos2->Print("studionoise/dimensioni_maggiori_rivelatori/RPC2_positions.png");
  crpcpos2->Print("studionoise/dimensioni_maggiori_rivelatori/RPC2_positions.root");
  TCanvas *crpcpos3 = new TCanvas();
  crpcpos3->Divide(1,2);
  crpcpos3->cd(1);
  hrpcpos3big->GetXaxis()->SetTitle("cm");
  hrpcpos3big->GetYaxis()->SetTitle("cm");
  hrpcpos3big->Draw("COLZ");
  crpcpos3->cd(2);
  hrpcpos3->GetXaxis()->SetTitle("cm");
  hrpcpos3->GetYaxis()->SetTitle("cm");
  hrpcpos3->Draw("COLZ");
  crpcpos3->Print("studionoise/dimensioni_maggiori_rivelatori/RPC3_positions.png");
  crpcpos3->Print("studionoise/dimensioni_maggiori_rivelatori/RPC3_positions.root");
  TCanvas *crpcpos4 = new TCanvas();
  crpcpos4->Divide(1,2);
  crpcpos4->cd(1);
  hrpcpos4big->GetXaxis()->SetTitle("cm");
  hrpcpos4big->GetYaxis()->SetTitle("cm");
  hrpcpos4big->Draw("COLZ");
  crpcpos4->cd(2);
  hrpcpos4->GetXaxis()->SetTitle("cm");
  hrpcpos4->GetYaxis()->SetTitle("cm");
  hrpcpos4->Draw("COLZ");
  crpcpos4->Print("studionoise/dimensioni_maggiori_rivelatori/RPC4_positions.png");
  crpcpos4->Print("studionoise/dimensioni_maggiori_rivelatori/RPC4_positions.root");
  TCanvas *crpcpos5 = new TCanvas();
  crpcpos5->Divide(1,2);
  crpcpos5->cd(1);
  hrpcpos5big->GetXaxis()->SetTitle("cm");
  hrpcpos5big->GetYaxis()->SetTitle("cm");
  hrpcpos5big->Draw("COLZ");
  crpcpos5->cd(2);
  hrpcpos5->GetXaxis()->SetTitle("cm");
  hrpcpos5->GetYaxis()->SetTitle("cm");
  hrpcpos5->Draw("COLZ");
  crpcpos5->Print("studionoise/dimensioni_maggiori_rivelatori/RPC5_positions.png");
  crpcpos5->Print("studionoise/dimensioni_maggiori_rivelatori/RPC5_positions.root");  
  TCanvas *crpcpos6 = new TCanvas();
  crpcpos6->Divide(1,2);
  crpcpos6->cd(1);
  hrpcpos6big->GetXaxis()->SetTitle("cm");
  hrpcpos6big->GetYaxis()->SetTitle("cm");
  hrpcpos6big->Draw("COLZ");
  crpcpos6->cd(2);
  hrpcpos6->GetXaxis()->SetTitle("cm");
  hrpcpos6->GetYaxis()->SetTitle("cm");
  hrpcpos6->Draw("COLZ");
  crpcpos6->Print("studionoise/dimensioni_maggiori_rivelatori/RPC6_positions.png");
  crpcpos6->Print("studionoise/dimensioni_maggiori_rivelatori/RPC6_positions.root");  

  cout<<"Fraction of entries in range of histogram spectrometer 1: "<<(Double_t) hspectropos1->Integral() / hspectropos1->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram spectrometer 2: "<<(Double_t) hspectropos2->Integral() / hspectropos2->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram spectrometer 3: "<<(Double_t) hspectropos3->Integral() / hspectropos3->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram spectrometer 4: "<<(Double_t) hspectropos4->Integral() / hspectropos4->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram rpc 1: "<<(Double_t) hrpcpos1->Integral() / hrpcpos1->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram rpc 2: "<<(Double_t) hrpcpos2->Integral() / hrpcpos2->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram rpc 3: "<<(Double_t) hrpcpos3->Integral() / hrpcpos3->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram rpc 4: "<<(Double_t) hrpcpos4->Integral() / hrpcpos4->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram rpc 5: "<<(Double_t) hrpcpos5->Integral() / hrpcpos5->GetEntries()<<endl;
  cout<<"Fraction of entries in range of histogram rpc 6: "<<(Double_t) hrpcpos6->Integral() / hrpcpos6->GetEntries()<<endl;
  cout<<endl;
  cout<<"Fraction muons isolated in almost two RPCS"<<(Double_t) hnrpcisolated->Integral(3,6)/hnrpcisolated->GetEntries()<<endl;
  cout<<"Fraction muons isolated in almost three RPCS"<<(Double_t) hnrpcisolated->Integral(4,6)/hnrpcisolated->GetEntries()<<endl;
  
  hx1->Scale(1./ninterazioni);
  hy1->Scale(1./ninterazioni);
  hx2->Scale(1./ninterazioni);
  hy2->Scale(1./ninterazioni);
  hx3->Scale(1./ninterazioni);
  hy3->Scale(1./ninterazioni);
  hx4->Scale(1./ninterazioni);
  hy4->Scale(1./ninterazioni);
  hx5->Scale(1./ninterazioni);
  hy5->Scale(1./ninterazioni);
  
  TCanvas *cpos1 = new TCanvas();
  cpos1->Divide(1,2);
  cpos1->cd(1);
  hx1->GetXaxis()->SetTitle("cm");
  hx1->Draw();
  cpos1->cd(2);
  hy1->Draw();
  hy1->GetXaxis()->SetTitle("cm");
  TCanvas *cpos2 = new TCanvas();
  cpos2->Divide(1,2);
  cpos2->cd(1);
  hx2->Draw();
  hx2->GetXaxis()->SetTitle("cm");
  cpos2->cd(2);
  hy2->Draw();
  hy2->GetXaxis()->SetTitle("cm");
  TCanvas *cpos3 = new TCanvas();
  cpos3->Divide(1,2);
  cpos3->cd(1);
  hx3->Draw();
  hx3->GetXaxis()->SetTitle("cm");
  cpos3->cd(2);
  hy3->Draw();
  hy3->GetXaxis()->SetTitle("cm");
  TCanvas *cpos4 = new TCanvas();
  cpos4->Divide(1,2);
  cpos4->cd(1);
  hx4->Draw();
  hx4->GetXaxis()->SetTitle("cm");
  cpos4->cd(2);
  hy4->Draw();
  hy4->GetXaxis()->SetTitle("cm");
  TCanvas *cpos5 = new TCanvas();
  cpos5->Divide(1,2);
  cpos5->cd(1);
  hx5->Draw();
  hx5->GetXaxis()->SetTitle("cm");
  cpos5->cd(2);
  hy5->Draw();
  hy5->GetXaxis()->SetTitle("cm");

  TCanvas *c6 = new TCanvas();
  hnrpcmuons->Draw();

  TCanvas *c7 = new TCanvas();
  hnrpcisolated->Draw();
}

void emulsion(){
 TFile *file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4.root");
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");

 Int_t whatfilm = 45; //numero del film di cui si vuole salvare i display
 
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

 FairMCEventHeader *prova;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.", &prova);
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);
 
 TClonesArray *arr2 = new TClonesArray("BoxPoint",1000);
 cbmsim->GetBranch("BoxPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("BoxPoint", &arr2);
 
 cbmsim->SetBranchStatus("*",1);

 Double_t xPoint, yPoint, zPoint, impulso;

 //istogramma film di emulsione
 TH2D *hxy = new TH2D("hxy", "Posizione tracce nell'evento", 100, -5 , 5, 100, -5, 5); //primo piano
 
 const int nevents = cbmsim->GetEntries();

 gROOT->SetBatch(kTRUE);


 ifstream input; 
 //input.open("posizioni_lastrine_new3.txt");
 input.open("posizioni_lastrine_cut_rate.txt"); //apro il file con salvate le posizioni delle lastre
 vector<Double_t> zlastre;
 Double_t tempvalue;

 Int_t nfilm;
 
 std::vector<Double_t>::iterator low,up;
 
 while (input.good()){
   input>>tempvalue;
   zlastre.push_back(tempvalue);
 }
 sort(zlastre.begin(), zlastre.end());
 input.close();
    //preparo le canvas per i display

 Int_t ntracks = 0;
 Int_t pdgcode = 0;
 
 for (int i = 0; i < 1000; i++){

   arr1->Clear();
   arr2->Clear();
   cbmsim->GetEntry(i);

   if (i % 10 == 0) cout<<"Inizio evento: "<<i<<endl;
   //prendo la posizione del vertice.
   ShipMCTrack *first = (ShipMCTrack*) arr1->At(0);
   Double_t Vx = first->GetStartX();
   Double_t Vy = first->GetStartY();
   cout<<Vx<<" "<<Vy<<" "<<i<<endl;

   hxy->Reset();
   //inizio il ciclo sugli hit
   for (int j = 0; j < arr2->GetEntriesFast(); j++){
     BoxPoint *box =(BoxPoint*) arr2->At(j);
     xPoint = box->GetX();
     yPoint = box->GetY();
     zPoint = box->GetZ();
     impulso = pow((pow(box->GetPx(),2) + pow(box->GetPy(),2) + pow(box->GetPz(),2)),0.5);
     pdgcode = box->PdgCode();
	
     bool charge = false;
     if (pdgcode > 10000.) charge = true; //all the hits with pdg higher than 10000 are charged
     else if (abs(pdg->GetParticle(pdgcode)->Charge()) >0) charge = true;
     if (charge == false) continue;
     
     low = lower_bound(zlastre.begin(), zlastre.end(), box->GetZ());
     up = lower_bound(zlastre.begin(), zlastre.end(), box->GetZ());
     low--;

     if (TMath::Abs(*(low) - box->GetZ()) < TMath::Abs(*(up) - box->GetZ())) nfilm = low - zlastre.begin() + 1;
     else nfilm = up-zlastre.begin() + 1;
     if (up == zlastre.end()) nfilm = low - zlastre.begin(); //per risolvere il problema dell'ultimo film

     if (pdgcode < 10000){
       if (abs(pdg->GetParticle(pdgcode)->Charge()) > 0){	 	
	 if (nfilm == whatfilm){ //mi trovo nella lastra numero...
	   hxy->Fill(xPoint,yPoint);
	 }
       }
     }
   }
   
   TCanvas *c1 = new TCanvas();
   hxy->GetXaxis()->SetTitle("cm");
   hxy->GetYaxis()->SetTitle("cm");
   hxy->Draw();
   ntracks+=hxy->GetEntries();
   
   TMarker * marker = new TMarker(Vx, Vy, 29);
   marker->SetMarkerColor(kGreen);
   marker->Draw("SAME");
   c1->Print(Form("./displays/displays_mobile/evento%i/Emu_%i.png",i+1, whatfilm));

   delete c1;
     } //fine ciclo sugli eventi
 cout<<"NUMERO TOTALE DI TRACCE NEL FILM NUMERO "<<whatfilm<<" PER 1000 EVENTI: "<< ntracks<<endl; //vediamo che numerone esce :(
}
void spectrometer(){
// TFile *file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_isolamentomuoni_10000.root");
  //TFile *file1 = TFile::Open("./studionoise/studio_precampomagnetico/ship.10.0.Pythia8CharmOnly-TGeant4_8cm.root");
 TFile *file1 = TFile::Open("/home/utente/SHIPBuild/sim_charmdet/charm_events/RUN3_nofield_1cmprespectrometer_1int/ship.conical.Pythia8CharmOnly-TGeant4.root");
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");


 
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

 FairMCEventHeader *prova = NULL;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.", &prova);
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);
 
 
 TClonesArray *arr2 = new TClonesArray("SpectrometerPoint",1000);
 cbmsim->GetBranch("SpectrometerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("SpectrometerPoint", &arr2);

 TClonesArray *arr3 = new TClonesArray("MuonTaggerPoint",1000);
 cbmsim->GetBranch("MuonTaggerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MuonTaggerPoint", &arr3);

 cbmsim->SetBranchStatus("*",1);

 // Double_t Vz;
 Double_t xPoint, yPoint, zPoint, impulso;
 int nplane; //numero di piano per l'hit (0-4 per lo spettrometro, 5-10 per il filtro per muoni)

 //istogrammi spettrometro (bin da 200 micron, x e y per T1 e T2)
/* TH1D *hx1 = new TH1D("hx1", "Strip x primo piano spettrometro", 700, -7, 7); //primo piano
 TH1D *hy1 = new TH1D("hy1", "Strip y primo piano spettrometro", 600, -6, 6);
 TH1D *hx2 = new TH1D("hx2", "Strip x secondo piano spettrometro", 1000, -10, 10); //secondo piano
 TH1D *hy2 = new TH1D("hy2", "Strip y secondo piano spettrometro", 1000, -10, 10); 
 TH1D *hx3 = new TH1D("hx3", "Strip x terzo piano spettrometro", 1000, -100, 100); //terzo piano
 TH1D *hy3 = new TH1D("hy3", "Strip y terzo piano spettrometro", 800, -80, 80);
 TH1D *hx4 = new TH1D("hx4", "Strip x quarto piano spettrometro", 1000, -100, 100); //quarto piano
 TH1D *hy4 = new TH1D("hy5", "Strip y quarto piano spettrometro", 800, -80, 80);*/
 //istogrammi spettroemtero (il pitch utilizzato è 50 x 250 micron per T1 e T2, 200 x 200 micron per T3 e T4)
 /*TH2D *hxy1 = new TH2D("hxy1", "Posizione xy primo piano spettrometro", 700, -7, 7, 600, -6, 6);
 TH2D *hxy2 = new TH2D("hxy2", "Posizione xy secondo piano spettrometro", 1000, -10, 10, 1000, -10, 10);
 TH2D *hxy3 = new TH2D("hxy3", "Posizione xy terzo piano spettrometro", 10000, -100, 100, 800, -80, 80);
 TH2D *hxy4 = new TH2D("hxy4", "Posizione xy quarto piano spettrometro", 10000, -100, 100, 800, -80, 80);*/

 Double_t DeltaTheta = 0.003;
 Double_t Deltaxstart = 1.;
 Double_t Deltax = Deltaxstart * DeltaTheta;

 Int_t Nbins = (Int_t) 4./(3*Deltax);
 TH1D *hcharmoccupancy = new TH1D("hcharmoccupancy", "Charm and other particles in same cell",5,0,5);
 TH1D *hoccupancy = new TH1D("hoccupancy", "Occupancy di ogni bin dell'istogramma",10,0,10);


 TH1D *hx1 = new TH1D("hx1", "Strip x primo piano spettrometro", Nbins, -2, 2); //primo piano
 TH1D *hy1 = new TH1D("hy1", "Strip y primo piano spettrometro", Nbins, -2, 2);
 TH1D *hx2 = new TH1D("hx2", "Strip x secondo piano spettrometro", 800, -10, 10); //secondo piano
 TH1D *hy2 = new TH1D("hy2", "Strip y secondo piano spettrometro", 4000, -10, 10); 
 TH1D *hx3 = new TH1D("hx3", "Strip x terzo piano spettrometro", 8000, -100, 100); //terzo piano
 TH1D *hy3 = new TH1D("hy3", "Strip y terzo piano spettrometro", 6400, -80, 80);
 TH1D *hx4 = new TH1D("hx4", "Strip x quarto piano spettrometro", 8000, -100, 100); //quarto piano
 TH1D *hy4 = new TH1D("hy5", "Strip y quarto piano spettrometro", 6400, -80, 80);

 TH2D *hxy1 = new TH2D("hxy1", "Posizione xy primo piano spettrometro", Nbins, -2, 2, Nbins, -2, 2);
 TH2D *hxy2 = new TH2D("hxy2", "Posizione xy secondo piano spettrometro", 800, -10, 10, 4000, -10, 10);
 TH2D *hxy3 = new TH2D("hxy3", "Posizione xy terzo piano spettrometro", 8000, -100, 100, 6400, -80, 80);
 TH2D *hxy4 = new TH2D("hxy4", "Posizione xy quarto piano spettrometro", 8000, -100, 100, 6400, -80, 80);

 // numero di strip accese
 TH1D *nx1 = new TH1D("nx1","Numero strip accese", 560, 0, 560);
 TH1D *ny1 = new TH1D("ny1","Numero strip accese", 480, 0, 480);
 TH1D *nx2 = new TH1D("nx2","Numero strip accese", 800, 0, 800);
 TH1D *ny2 = new TH1D("ny2","Numero strip accese", 800, 0, 800);
 TH1D *nx3 = new TH1D("nx3","Numero strip accese", 8000, 0, 8000);
 TH1D *ny3 = new TH1D("ny3","Numero strip accese", 6400, 0, 6400);
 TH1D *nx4 = new TH1D("nx4","Numero strip accese", 8000, 0, 8000);
 TH1D *ny4 = new TH1D("ny4","Numero strip accese", 6400, 0, 6400);
 
 int nevents;
 if (cbmsim->GetEntries() < 1000) nevents = cbmsim->GetEntries();
 else nevents = 1000;

 nevents = cbmsim->GetEntries();
 //preparazione contatori per i 1000 eventi

 Int_t nbins1, nbins2, nbins3, nbins4; //numero TOTALE di bin non vuoti
 Int_t ninbins1, ninbins2, ninbins3, ninbins4; //numero TOTALE di entries nei bin non vuoti

 Int_t nbinsy1, nbinsy2, nbinsy3, nbinsy4; //numero di bin non vuoti lungo y
 Int_t nbinsx1, nbinsx2, nbinsx3, nbinsx4; //numero di bin non vuoti lungo x
 Int_t ninbinsx1, ninbinsx2, ninbinsx3, ninbinsx4; //numero di entries nei bin non vuoti lungo x
 
 Int_t nentries1, nentries2, nentries3, nentries4; //numero di entries
 Int_t mumid, mumpdg, intermediatepdg;
 Double_t rmsx1, rmsx2, rmsx3, rmsx4; //calcolo l'rms medio lungo x
 Double_t rmsy1, rmsy2, rmsy3, rmsy4; //calcolo l'rms medio lungo y

 vector<Double_t> charmdaughterx1, charmdaughterx2, charmdaughterx3, charmdaughterx4;
 vector<Double_t> charmdaughtery1, charmdaughtery2, charmdaughtery3, charmdaughtery4;

 Int_t ncharmdaughters, ncharmdaughterscontigued; //numero figlie di charm nel primo piano dello spettrometro
 Int_t ntotcharmdaughters = 0;
 Double_t fracdaughterscontigued = 0;
 
 Int_t nhits, ncontigued;
 Int_t neventscontigued = 0;
 bool almosttwocontigued, charmcontigued1, charmcontigued2, charmcontigued3, charmcontigued4; //segnala la presenza di almeno due hit contigued in an event 
 Int_t neventscharmcontigued1 = 0, neventscharmcontigued2 = 0, neventscharmcontigued3 = 0, neventscharmcontigued4 = 0; 

 rmsx1 = 0.;
 rmsy1 = 0.;
 rmsx2 = 0.;
 rmsy2 = 0.;
 rmsx3 = 0.;
 rmsy3 = 0.;
 rmsx4 = 0.;
 rmsy4 = 0.;

 nbins1 = 0; 
 ninbins1 = 0;
 nbins2 = 0;
 ninbins2 = 0;
 nbins3 = 0;
 ninbins3 = 0;
 nbins4 = 0;
 ninbins4 = 0;



 nentries1 = 0;
 nentries2 = 0;
 nentries3 = 0;
 nentries4 = 0;
 

 nhits = 0;
 ncontigued = 0;
 ncharmdaughters = 0;
 
 gROOT->SetBatch(kTRUE);
 
 
 
 cout<<"Numero totale di eventi da salvare:"<<nevents<<endl;


 TFile *outfile = new TFile("displays/spectrometer_display_plots/spectrometer_display_plots.root","RECREATE");
// INIZIO DEL CICLO SUGLI EVENTI
 
for (int i = 0; i < nevents; i++){

   arr1->Clear();
   arr2->Clear();
   arr3->Clear();
   cbmsim->GetEntry(i);
   
   //prendo la posizione del vertice.
   ShipMCTrack *first = (ShipMCTrack*) arr1->At(0);
   Double_t Vx = first->GetStartX();
   Double_t Vy = first->GetStartY();     
   //reset degli istogrammi e del numero di bin non nulli lungo x e y
   
   hx1->Reset();
   hy1->Reset();
   hxy1->Reset();
   hx2->Reset();
   hy2->Reset();
   hxy2->Reset();
   hx3->Reset();
   hy3->Reset();
   hxy3->Reset();
   hx4->Reset();
   hy4->Reset();
   hxy4->Reset();

   nbinsx1 = 0; 
   ninbinsx1 = 0;
   nbinsx2 = 0;
   ninbinsx2 = 0;
   nbinsx3 = 0;
   ninbinsx3 = 0;
   nbinsx4 = 0;
   ninbinsx4 = 0;
   
   nbinsy1 = 0;
   nbinsy2 = 0;
   nbinsy3 = 0;
   nbinsy4 = 0;

   ncharmdaughters = 0;
   almosttwocontigued = false; 
   charmcontigued1 = false;
   charmcontigued2 = false;
   charmcontigued3 = false;
   charmcontigued4 = false;
   
   cout<<i<<endl;

   //preparazione delle canvas
   
   TCanvas *c1_1 = new TCanvas(); 
   c1_1->Range(-7,-6,7,6);
   TH1F* h1_1 = c1_1->DrawFrame(-7,-6,7,6);

   TCanvas *c1_2 = new TCanvas();
   c1_2->Range(-7,-6,7,6);
   TH1F* h1_2 = c1_2->DrawFrame(-7,-6,7,6);

   TCanvas *c2_1 = new TCanvas();
   c2_1->Range(-10,-10,10,10);
   TH1F* h2_1 = c2_1->DrawFrame(-10,-10,10,10);

   TCanvas *c2_2 = new TCanvas();
   c2_2->Range(-10,-10,10,10);
   TH1F* h2_2 = c2_2->DrawFrame(-10,-10,10,10);

   TCanvas * c3_1 = new TCanvas();
   c3_1->Range(-100,-80,100, 80);
   TH1F* h3_1 = c3_1->DrawFrame(-100,-80,100,80);

   TCanvas *c3_2 = new TCanvas();
   c3_2->Range(-100,-80,100,80);
   TH1F* h3_2 = c3_2->DrawFrame(-100,-80,100,80);
   
   TCanvas *c4_1 = new TCanvas();
   c4_1->Range(-100,-80,100,80);
   TH1F* h4_1 = c4_1->DrawFrame(-100,-80,100,80);

   TCanvas *c4_2 = new TCanvas();
   c4_2->Range(-100,-80,100,80);
   TH1F* h4_2 = c4_2->DrawFrame(-100,-80,100,80);
   //disegno un punto verde per il vertice su tutti e quattro i piani
   c1_1->cd();
   TMarker* marker = new TMarker(Vx, Vy, 29); //faccio nuovi new, di norma andrebbe fatto delete, ma temo che cancelli anche il disegno effettuato, essendo la canvas interattiva. Dopo aver disegnato il punto, comunque dovrei cancellarlo.
   marker->SetMarkerColor(kGreen);
   marker->Draw();
   c2_1->cd();
   marker = new TMarker(Vx, Vy, 29);
   marker->SetMarkerColor(kGreen);
   marker->Draw();
   c3_1->cd();
   marker = new TMarker(Vx, Vy, 29);
   marker->SetMarkerColor(kGreen);
   marker->Draw();
   c4_1->cd();
   marker = new TMarker(Vx, Vy, 29);
   marker->SetMarkerColor(kGreen);
   marker->Draw();


   
   //inizio il ciclo sugli hit
   for (int j = 0; j < arr2->GetEntriesFast(); j++){
     SpectrometerPoint *spectro =(SpectrometerPoint*) arr2->At(j);
     xPoint = spectro->GetX();
     yPoint = spectro->GetY();
     zPoint = spectro->GetZ(); 

     Int_t pdgcode = spectro->PdgCode();     
     bool charge = false;    
     if (pdgcode > 10000.) charge = true; //all the hits with pdg higher than 10000 are charged
     else if (abs(pdg->GetParticle(pdgcode)->Charge()) >0) charge = true;
     if (charge == false) continue;
     impulso = pow((pow(spectro->GetPx(),2) + pow(spectro->GetPy(),2) + pow(spectro->GetPz(),2)),0.5);
   
     nplane = 0;
/*     if (zPoint < 0.2) nplane = 1; //prendere il primo o l'ultimo piano quando ho più piani per stazione?
     else if ((zPoint > 10) && (zPoint < 20.75)) nplane = 2;
     else if ((zPoint > 100) && (zPoint < 490)) nplane = 3;
     else if (zPoint > 490) nplane = 4;*/
//nuova configurazione
     if (spectro->GetDetectorID() == 101) nplane = 1;
     if (spectro->GetDetectorID() == 103) nplane = 2;
     if (spectro->GetDetectorID() == 3) nplane = 3;
     if (spectro->GetDetectorID() == 4) nplane = 4;

     //studio l'indice della traccia per capire se è figlia di charm    
     if (spectro->GetTrackID() >= 0){ //hit con meno di 100 MeV hanno TrackID -1 e non sono dunque tracciabili
	   ShipMCTrack * track = (ShipMCTrack*) arr1->At(spectro->GetTrackID());	 
	   mumid = track->GetMotherId();
	   
	   //while(mumid > 0){ //fino al charm o alle particelle prodotte con esso	  
	   if(mumid >= 0){	   
	    ShipMCTrack *mother = (ShipMCTrack*) arr1->At(mumid);
	    mumpdg = mother->GetPdgCode();

	    while ((isintermediate(mumpdg) == true) && (mumid>0)){	
	      intermediatepdg = mumpdg;
	      mumid = mother->GetMotherId();
	      if (mumid > 0) mother = (ShipMCTrack*) arr1->At(mumid);	  
	      if (mumid > 0) mumpdg = mother->GetPdgCode();
	    }	   
	   
	   if (ischarm(mumpdg) == true){ //figlia di uno dei due charm            
	     if (nplane == 1){
	       ncharmdaughters++;
	       ntotcharmdaughters++;
	       charmdaughterx1.push_back(xPoint);
	       charmdaughtery1.push_back(yPoint);
	       cout<<i<<" "<<track->GetPdgCode()<<" "<<mumpdg<<" "<<xPoint<<" "<<yPoint<<endl;
	     }
	     if (nplane == 2){
	       charmdaughterx2.push_back(xPoint);
	       charmdaughtery2.push_back(yPoint);
	     }
	     if (nplane == 3){
	       charmdaughterx3.push_back(xPoint);
	       charmdaughtery3.push_back(yPoint);
	     }
	     if (nplane == 4){
	       charmdaughterx4.push_back(xPoint);
	       charmdaughtery4.push_back(yPoint);
	     }
	     
	   }
	 }
}  
     
	
     if (nplane == 1){ //mi trovo nel primo piano
       hx1->Fill(xPoint);
       hy1->Fill(yPoint);
       hxy1->Fill(xPoint,yPoint);
       c1_1->cd();        
       
       TBox * box = new TBox(hx1->GetBinLowEdge(hx1->FindBin(xPoint)), hy1->GetBinLowEdge(hy1->FindBin(yPoint)), hx1->GetBinLowEdge(hx1->FindBin(xPoint)+1), hy1->GetBinLowEdge(hy1->FindBin(yPoint)+1));
       if (abs(spectro->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(spectro->PdgCode()) == 11) box->SetLineColor(kOrange);	 
	 /* if (spectro->GetTrackID() >= 0){
	   ShipMCTrack * track = (ShipMCTrack*) arr1->At(spectro->GetTrackID());	 
	   Double_t mumid = track->GetMotherId();
	   
	   while(mumid > 0){ //fino al charm o alle particelle prodotte con esso
	    ShipMCTrack *mother = (ShipMCTrack*) arr1->At(mumid);
	    if (mumid == mother->GetMotherId()){
	      ShipMCTrack *mother = (ShipMCTrack*) arr1->At(mumid - 1); //per il baco sulla figlia del secondo charm
	    }
	    mumid = mother->GetMotherId();
	   }
	   if (mumid == 0){ //uno dei due charm
	     box->SetLineColor(kYellow);
	   }
	   }*/	   
       box->Draw("l");
     }

     if (nplane == 2){ //mi trovo nel secondo piano
       hx2->Fill(xPoint);
       hy2->Fill(yPoint);
       hxy2->Fill(xPoint,yPoint);
       c2_1->cd();
       TBox * box = new TBox(hx2->GetBinLowEdge(hx2->FindBin(xPoint)), hy2->GetBinLowEdge(hy2->FindBin(yPoint)), hx2->GetBinLowEdge(hx2->FindBin(xPoint)+1), hy2->GetBinLowEdge(hy2->FindBin(yPoint)+1));
       if (abs(spectro->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(spectro->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }
     
     if (nplane == 3){ //mi trovo nel terzo piano
       hx3->Fill(xPoint);
       hy3->Fill(yPoint);
       hxy3->Fill(xPoint,yPoint);
       c3_1->cd();
       TBox * box = new TBox(hx3->GetBinLowEdge(hx3->FindBin(xPoint)), hy3->GetBinLowEdge(hy3->FindBin(yPoint)), hx3->GetBinLowEdge(hx3->FindBin(xPoint)+1), hy3->GetBinLowEdge(hy3->FindBin(yPoint)+1));
       if (abs(spectro->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(spectro->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }

     if (nplane == 4){ //mi trovo nel quarto piano
       hx4->Fill(xPoint);
       hy4->Fill(yPoint);
       hxy4->Fill(xPoint,yPoint);
       c4_1->cd();
       TBox * box = new TBox(hx4->GetBinLowEdge(hx4->FindBin(xPoint)), hy4->GetBinLowEdge(hy4->FindBin(yPoint)), hx4->GetBinLowEdge(hx4->FindBin(xPoint)+1), hy4->GetBinLowEdge(hy4->FindBin(yPoint)+1));
       if (abs(spectro->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(spectro->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }
     
     
     
   }
  
   if (i % 100 == 0) cout<<"Fine evento numero :"<<i<<endl;

   nentries1 += hxy1->GetEntries();
   nentries2 += hxy2->GetEntries();
   nentries3 += hxy3->GetEntries();
   nentries4 += hxy4->GetEntries();
   
   TCanvas *c1_3 = new TCanvas();
   hxy1->SetMarkerStyle(7);
   hxy1->Draw("COLZ");
   marker->Draw("SAME");
   for (int imarker = 0; imarker < ncharmdaughters; imarker++){
     Int_t charmdaughterbinx = hxy1->GetXaxis()->FindBin(charmdaughterx1[imarker]);
     Int_t charmdaughterbiny = hxy1->GetYaxis()->FindBin(charmdaughtery1[imarker]);
    // cout<<"Integral in prossimità del charm: "<<hxy1->Integral(charmdaughterbinx - 1, charmdaughterbinx + 1, charmdaughterbiny - 1, charmdaughterbiny + 1)<<endl;
     //cout<<"Bin Content figlio charm: "<<hxy1->GetBinContent(charmdaughterbinx, charmdaughterbiny)<<endl;
     if(hxy1->Integral(charmdaughterbinx, charmdaughterbinx, charmdaughterbiny, charmdaughterbiny) > 1){     
     //if(hxy1->GetBinContent(charmdaughterbinx, charmdaughterbiny) > 1){
         marker = new TMarker(Vx, Vy, 29);
         marker->SetMarkerColor(kGreen);
         marker->Draw("SAME");
	 ncharmdaughterscontigued++;
         charmcontigued1 = true;
    }
    hcharmoccupancy->Fill(hxy1->Integral(charmdaughterbinx, charmdaughterbinx, charmdaughterbiny, charmdaughterbiny));
     //cout<<"Entries nel punto del charm in T1: "<<hxy1->GetBinContent(hxy1->GetBin(charmdaughterbinx, charmdaughterbiny))<<endl;
     /*TMarker * charmmarker = new TMarker(charmdaughterx[imarker], charmdaughtery[imarker],29);
     charmmarker->SetMarkerColor(kRed);
     charmmarker->Draw("SAME");*/
   }
   for (int nbin = 1; nbin <= hxy1->GetSize();nbin++) hoccupancy->Fill(hxy1->GetBinContent(nbin));
   cout<<"Numero di hit totali figlie di charm: "<<ntotcharmdaughters<<" size vector: "<<charmdaughterx1.size()<<endl;
   charmdaughterx1.clear();
   charmdaughtery1.clear();
   TCanvas *c2_3 = new TCanvas();
   hxy2->SetMarkerStyle(7);
   hxy2->Draw("COLZ");
   marker->Draw("SAME");
/*
   for (int idaughter = 0; idaughter < charmdaughterx2.size(); idaughter++){
     Int_t charmdaughterbinx = hxy2->GetXaxis()->FindBin(charmdaughterx2[idaughter]);
     Int_t charmdaughterbiny = hxy2->GetYaxis()->FindBin(charmdaughtery2[idaughter]);
     if(hxy2->Integral(charmdaughterbinx - 1, charmdaughterbinx + 1, charmdaughterbiny - 1, charmdaughterbiny + 1) > 1){       	 
         charmcontigued2 = true;
    }
   }  
   charmdaughterx2.clear();
   charmdaughtery2.clear();
   
   TCanvas *c3_3 = new TCanvas();
   hxy3->SetMarkerStyle(7);
   hxy3->Draw("COLZ");
   marker->Draw("SAME");

     for (int idaughter = 0; idaughter < charmdaughterx3.size(); idaughter++){
     Int_t charmdaughterbinx = hxy3->GetXaxis()->FindBin(charmdaughterx3[idaughter]);
     Int_t charmdaughterbiny = hxy3->GetYaxis()->FindBin(charmdaughtery3[idaughter]);
     if(hxy3->Integral(charmdaughterbinx - 1, charmdaughterbinx + 1, charmdaughterbiny - 1, charmdaughterbiny + 1) > 1){       	 
         charmcontigued3 = true;
    }
   }  
   charmdaughterx3.clear();
   charmdaughtery3.clear();
   TCanvas *c4_3 = new TCanvas();
   hxy4->SetMarkerStyle(7);
   hxy4->Draw("COLZ");
   marker->Draw("SAME"); 
   // cout<<"PROVA "<<i<<" size "<<charmdaughterx4.size()<<endl;
     for (int idaughter = 0; idaughter < charmdaughterx4.size(); idaughter++){
     Int_t charmdaughterbinx = hxy4->GetXaxis()->FindBin(charmdaughterx2[idaughter]);
     Int_t charmdaughterbiny = hxy4->GetYaxis()->FindBin(charmdaughtery2[idaughter]);
     if(hxy4->Integral(charmdaughterbinx - 1, charmdaughterbinx + 1, charmdaughterbiny - 1, charmdaughterbiny + 1) > 1){       	 
         charmcontigued4 = true;
    }
   }  */ 
   charmdaughterx4.clear();
   charmdaughtery4.clear();
   

   c1_3->Write();
   
   delete c1_1;
   delete c1_2;
   delete c1_3;
   delete c2_1;
   delete c2_2;
   delete c2_3;
   delete c3_1;
   delete c3_2;
   //delete c3_3;
   delete c4_1;
   delete c4_2;
  // delete c4_3;

   //riempo i grafici con il numero di strip accese lungo x e lungo y
   nx1->Fill(nbinsx1);
   ny1->Fill(nbinsy1);
   nx2->Fill(nbinsx2);
   ny2->Fill(nbinsy2);
   nx3->Fill(nbinsx3);
   ny3->Fill(nbinsy3);
   nx4->Fill(nbinsx4);
   ny4->Fill(nbinsy4);

   //cout<<"Numero hit contigui "<<ncontigued<<" su un numero di hit totali: "<<nhits<<endl;
   if (charmcontigued1) neventscharmcontigued1++;
   if (charmcontigued2) neventscharmcontigued2++;
   if (charmcontigued3) neventscharmcontigued3++;
   if (charmcontigued4) neventscharmcontigued4++;
   cout<<"Almeno una figlia di charm contigua: "<<neventscharmcontigued1<<endl;
   cout<<"Numero eventi con almeno una figlia di charm contigua nella stazione 2: "<<neventscharmcontigued2<<endl;
   cout<<"Numero eventi con almeno una figlia di charm contigua nella stazione 3: "<<neventscharmcontigued3<<endl;
   cout<<"Numero eventi con almeno una figlia di charm contigua nella stazione 4: "<<neventscharmcontigued4<<endl;
     } //fine ciclo sugli eventi
   TCanvas *coccupancycharm = new TCanvas();
   hcharmoccupancy->Draw();
   hcharmoccupancy->Write();
   TCanvas *coccupancy = new TCanvas();
   hoccupancy->Draw();
   hoccupancy->Write();
   outfile->Close();
}

//--------------------------------------------------------------------------------------------------------------------------------

//parte relativa agli RPC del rivelatore per muoni (6 piani, dimensioni 200 x 150 cm, 2 cm circa di bin)
void RPC(){
  // TFile *file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_centralbeam_10000.root");
  // TFile *file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_muontaggerfe.root");
  TFile *file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_isolamentomuoni_10000.root");
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");
 
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

 FairMCEventHeader *prova;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.", &prova);
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);

 TClonesArray *arr2 = new TClonesArray("MuonTaggerPoint",1000);
 cbmsim->GetBranch("MuonTaggerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MuonTaggerPoint", &arr2);

 cbmsim->SetBranchStatus("*",1);

 // Double_t Vz;
 Double_t xPoint, yPoint, zPoint, impulso;
 int nplane; //numero di piano per l'hit (0-4 per lo spettrometro, 5-10 per il filtro per muoni)

 //istogrammi rpc (bin da 2 cm o 2mm, x e y)
 TH1D *hx1 = new TH1D("hx1", "Strip x primo RPC", 1000, -100, 100); //primo piano
 TH1D *hy1 = new TH1D("hy1", "Strip y primo RPC", 800, -80, 80);
 TH1D *hx2 = new TH1D("hx2", "Strip x secondo RPC", 1000, -100, 100); //secondo piano
 TH1D *hy2 = new TH1D("hy2", "Strip y secondo RPC", 800, -80, 80);
 TH1D *hx3 = new TH1D("hx3", "Strip x terzo RPC", 1000, -100, 100); //terzo piano
 TH1D *hy3 = new TH1D("hy3", "Strip y terzo RPC", 800, -80, 80);
 TH1D *hx4 = new TH1D("hx4", "Strip x quarto RPC", 1000, -100, 100); //quarto piano
 TH1D *hy4 = new TH1D("hy4", "Strip y quarto RPC", 800, -80, 80);
 TH1D *hx5 = new TH1D("hx5", "Strip x quinto RPC", 1000, -100, 100); //quinto piano
 TH1D *hy5 = new TH1D("hy5", "Strip y quinto RPC", 800, -80, 80);
 TH1D *hx6 = new TH1D("hx6", "Strip x sesto RPC", 1000, -100, 100); //sesto piano
 TH1D *hy6 = new TH1D("hy6", "Strip y sesto RPC", 800, -80, 80);

 //bin da 2cm o 2mm, cambiare opportunamente come da recente richiesta
 
 TH2D *hxy1 = new TH2D("hxy1", "xy first RPC", 1000, -100 , 100 , 800, -80, 80);
 TH2D *hxy2 = new TH2D("hxy2", "xy second RPC", 1000, -100 , 100 , 800, -80, 80);
 TH2D *hxy3 = new TH2D("hxy3", "xy third RPC", 1000, -100 , 100 , 800, -80, 80);
 TH2D *hxy4 = new TH2D("hxy4", "xy fourth RPC", 1000, -100 , 100 , 800, -80, 80);
 TH2D *hxy5 = new TH2D("hxy5", "xy fifth RPC", 1000, -100 , 100 , 800, -80, 80);
 TH2D *hxy6 = new TH2D("hxy6", "xy sixth RPC", 1000, -100 , 100 , 800, -80, 80);

  // numero di strip accese
 TH1D *nx1 = new TH1D("nx1","Numero strip accese", 100, 0, 100);
 TH1D *ny1 = new TH1D("ny1","Numero strip accese", 80, 0, 80);
 TH1D *nx2 = new TH1D("nx2","Numero strip accese", 100, 0, 100);
 TH1D *ny2 = new TH1D("ny2","Numero strip accese", 80, 0, 80);
 TH1D *nx3 = new TH1D("nx3","Numero strip accese", 100, 0, 100);
 TH1D *ny3 = new TH1D("ny3","Numero strip accese", 80, 0, 80);
 TH1D *nx4 = new TH1D("nx4","Numero strip accese", 100, 0, 100);
 TH1D *ny4 = new TH1D("ny4","Numero strip accese", 80, 0, 80);
 TH1D *nx5 = new TH1D("nx5","Numero strip accese", 100, 0, 100);
 TH1D *ny5 = new TH1D("ny5","Numero strip accese", 80, 0, 80);
 TH1D *nx6 = new TH1D("nx6","Numero strip accese", 100, 0, 100);
 TH1D *ny6 = new TH1D("ny6","Numero strip accese", 80, 0, 80);
 
 Int_t nevents;
 if (cbmsim->GetEntries() < 1000) nevents = cbmsim->GetEntries();
 else nevents = 1000;

 Double_t muonx[6], muony[6];
 
  //contatori per i 1000 eventi

 Int_t nbins1, nbins2, nbins3, nbins4, nbins5, nbins6; //numero di bin non vuoti
 Int_t ninbins1, ninbins2, ninbins3, ninbins4, ninbins5, ninbins6; //numero di entries nei bin non vuoti

 Int_t nbinsx1, nbinsx2, nbinsx3, nbinsx4, nbinsx5, nbinsx6;
 Int_t nbinsy1, nbinsy2, nbinsy3, nbinsy4, nbinsy5, nbinsy6;
 
 Int_t nentries1, nentries2, nentries3, nentries4, nentries5, nentries6;
 
 Double_t rmsx1, rmsx2, rmsx3, rmsx4, rmsx5, rmsx6; //calcolo l'rms medio lungo x
 Double_t rmsy1, rmsy2, rmsy3, rmsy4, rmsy5, rmsy6; //calcolo l'rms medio lungo y
 
 rmsx1 = 0.;
 rmsy1 = 0.;
 rmsx2 = 0.;
 rmsy2 = 0.;
 rmsx3 = 0.;
 rmsy3 = 0.;
 rmsx4 = 0.;
 rmsy4 = 0.;
 rmsx5 = 0.;
 rmsy5 = 0.;
 rmsx6 = 0.;
 rmsy6 = 0.;

 nbins1 = 0; 
 ninbins1 = 0;
 nbins2 = 0;
 ninbins2 = 0;
 nbins3 = 0;
 ninbins3 = 0;
 nbins4 = 0;
 ninbins4 = 0;
 nbins5 = 0;
 ninbins5 = 0;
 nbins6 = 0;
 ninbins6 = 0;

 nentries1 = 0;
 nentries2 = 0;
 nentries3 = 0;
 nentries4 = 0;
 nentries5 = 0;
 nentries6 = 0;
 
 gROOT->SetBatch(kTRUE);

 for (int i = 0; i < nevents; i++){

   arr1->Clear();
   arr2->Clear();
   cbmsim->GetEntry(i);

   //resetto gli istogrammi
   hx1->Reset();
   hy1->Reset();
   hx2->Reset();
   hy2->Reset();
   hx3->Reset();
   hy3->Reset();
   hx4->Reset();
   hy4->Reset();
   hx5->Reset();
   hy5->Reset();
   hx6->Reset();
   hy6->Reset();

   hxy1->Reset();
   hxy2->Reset();
   hxy3->Reset();
   hxy4->Reset();
   hxy5->Reset();
   hxy6->Reset();

   for (int k = 0; k < 6; k++){
     muonx[k] = 0.;
     muony[k] = 0.;
   }

   //azzero il contatore del numero di bin non vuoti lungo x e lungo y per evento
   
   nbinsx1 = 0;
   nbinsy1 = 0;
   nbinsx2 = 0;
   nbinsy2 = 0;
   nbinsx3 = 0;
   nbinsy3 = 0;
   nbinsx4 = 0;
   nbinsy4 = 0;
   nbinsx5 = 0;
   nbinsy5 = 0;
   nbinsx6 = 0;
   nbinsy6 = 0;
   
   //preparo le canvas per i display
   TCanvas * c1_1 = new TCanvas(); 
   c1_1->Range(-100,-80,100, 80);
   c1_1->DrawFrame(-100,-80,100,80);
   TLatex t(-60, 90, "Orange for electrons, red for muons");
   t.Draw();
   
   TCanvas * c1_2 = new TCanvas(); //nota: anche se per il momento non inserisco il codice per il riempimento di queste canvas, le preparo se serviranno in futuro.
   c1_2->Range(-100,-80,100, 80);
   c1_2->DrawFrame(-100,-80,100,80);

   TCanvas *c2_1 = new TCanvas();
   c2_1->Range(-100,-80,100, 80);
   c2_1->DrawFrame(-100,-80,100,80);
   t.Draw();

   TCanvas * c2_2 = new TCanvas();
   c2_2->Range(-100,-80,100, 80);
   c2_2->DrawFrame(-100,-80,100,80);

   TCanvas *c3_1 = new TCanvas();
   c3_1->Range(-100,-80,100, 80);
   c3_1->DrawFrame(-100,-80,100,80);
   t.Draw();
   
   TCanvas * c3_2 = new TCanvas();
   c3_2->Range(-100,-80,100,80);
   c3_2->DrawFrame(-100,-80,100,80);
   
   TCanvas *c4_1 = new TCanvas();
   c4_1->Range(-100,-80,100,80);
   c4_1->DrawFrame(-100,-80,100,80);
   t.Draw();
   
   TCanvas *c4_2 = new TCanvas();
   c4_2->Range(-100,-80,100,80);
   c4_2->DrawFrame(-100,-80,100,80);
   
   TCanvas * c5_1 = new TCanvas();
   c5_1->Range(-100,-80,100,80);
   c5_1->DrawFrame(-100,-80,100,80);
   t.Draw();
   
   TCanvas *c5_2 = new TCanvas();
   c5_2->Range(-100,-80,100,80);
   c5_2->DrawFrame(-100,-80,100,80);

   TCanvas * c6_1 = new TCanvas();
   c6_1->Range(-100,-80,100,80);
   c6_1->DrawFrame(-100,-80,100,80);
   t.Draw();
   
   TCanvas *c6_2 = new TCanvas();
   c6_2->Range(-100,-80,100,80);
   c6_2->DrawFrame(-100,-80,100,80);

   for (int j = 0; j < arr2->GetEntriesFast(); j++){
     MuonTaggerPoint *muon =(MuonTaggerPoint*) arr2->At(j);
     xPoint = muon->GetX();
     yPoint = muon->GetY();
     zPoint = muon->GetZ();
     Int_t pdgcode = muon->PdgCode();
     bool charge = false;
     if (pdgcode > 10000.) charge = true;
     else if (abs(pdg->GetParticle(pdgcode)->Charge()) >0) charge = true;
     if (charge == false) continue;
     
     impulso = pow((pow(muon->GetPx(),2) + pow(muon->GetPy(),2) + pow(muon->GetPz(),2)),0.5);
     
     
     
     if (zPoint < 550.) nplane = 1;
     else if (zPoint < 650.) nplane = 2;
     else if (zPoint < 700.) nplane = 3;
     else if (zPoint < 750.) nplane = 4;
     else if (zPoint< 800.) nplane = 5;
     else nplane = 6;

     /*if (zPoint < 680.) nplane = 1;
     else if (zPoint < 720.) nplane = 2;
     else if (zPoint < 760.) nplane = 3;
     else if (zPoint < 800.) nplane = 4;
     else if (zPoint< 820.) nplane = 5;
     else nplane = 6;*/

     if (abs(muon->PdgCode()) == 13){
       muonx[nplane-1] = xPoint;
       muony[nplane-1] = yPoint;
       cout<<endl;
       cout<<"MUONE NELL'EVENTO "<<i+1<<endl;
       cout<<endl;
     }
     //cout<<nplane<<" "<<muon->PdgCode()<<endl;
	
     if (nplane == 1){ //mi trovo nel primo piano
       hx1->Fill(xPoint);
       hy1->Fill(yPoint);
       hxy1->Fill(xPoint, yPoint);
       
       c1_1->cd();
       TBox * box = new TBox(hx1->GetBinLowEdge(hx1->FindBin(xPoint)), hy1->GetBinLowEdge(hy1->FindBin(yPoint)), hx1->GetBinLowEdge(hx1->FindBin(xPoint)+1), hy1->GetBinLowEdge(hy1->FindBin(yPoint)+1));
       if (abs(muon->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(muon->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }

     if (nplane == 2){ //mi trovo nel secondo piano
       hx2->Fill(xPoint);
       hy2->Fill(yPoint);
       hxy2->Fill(xPoint,yPoint);

       c2_1->cd();
       TBox * box = new TBox(hx2->GetBinLowEdge(hx2->FindBin(xPoint)), hy2->GetBinLowEdge(hy2->FindBin(yPoint)), hx2->GetBinLowEdge(hx2->FindBin(xPoint)+1), hy2->GetBinLowEdge(hy2->FindBin(yPoint)+1));
       if (abs(muon->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(muon->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }
     
     if (nplane == 3){ //mi trovo nel terzo piano
       hx3->Fill(xPoint);
       hy3->Fill(yPoint);
       hxy3->Fill(xPoint, yPoint);
	      
       c3_1->cd();
       TBox * box = new TBox(hx3->GetBinLowEdge(hx3->FindBin(xPoint)), hy3->GetBinLowEdge(hy3->FindBin(yPoint)), hx3->GetBinLowEdge(hx3->FindBin(xPoint)+1), hy3->GetBinLowEdge(hy3->FindBin(yPoint)+1));
       if (abs(muon->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(muon->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }

     if (nplane == 4){ //mi trovo nel quarto piano
       hx4->Fill(xPoint);
       hy4->Fill(yPoint);
       hxy4->Fill(xPoint, yPoint);
       
       c4_1->cd();
       TBox * box = new TBox(hx4->GetBinLowEdge(hx4->FindBin(xPoint)), hy4->GetBinLowEdge(hy4->FindBin(yPoint)), hx4->GetBinLowEdge(hx4->FindBin(xPoint)+1./2.), hy4->GetBinLowEdge(hy4->FindBin(yPoint)+1./2.));
       if (abs(muon->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(muon->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }

     if (nplane == 5){ //mi trovo nel quinto piano
       hx5->Fill(xPoint);
       hy5->Fill(yPoint);
       hxy5->Fill(xPoint, yPoint);
       
       c5_1->cd();
       TBox * box = new TBox(hx5->GetBinLowEdge(hx5->FindBin(xPoint)), hy5->GetBinLowEdge(hy5->FindBin(yPoint)), hx5->GetBinLowEdge(hx5->FindBin(xPoint)+1), hy5->GetBinLowEdge(hy5->FindBin(yPoint)+1));
       if (abs(muon->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(muon->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }

     if (nplane == 6){ //mi trovo nel sesto piano
       hx6->Fill(xPoint);
       hy6->Fill(yPoint);
       hxy6->Fill(xPoint, yPoint);
       
       c6_1->cd();
       TBox * box = new TBox(hx6->GetBinLowEdge(hx6->FindBin(xPoint)), hy4->GetBinLowEdge(hy6->FindBin(yPoint)), hx4->GetBinLowEdge(hx6->FindBin(xPoint)+1), hy6->GetBinLowEdge(hy6->FindBin(yPoint)+1));
       if (abs(muon->PdgCode()) == 13) box->SetLineColor(kRed);
       else if (abs(muon->PdgCode()) == 11) box->SetLineColor(kOrange);
       box->Draw("l");
     }

     //inizio ciclo sui bin degli istogrammi
     for (int k = 0; k < hxy1->GetNbinsX(); k++){ //faccio il ciclo sulle x, quindi faccio il ciclo sulle y solo dove il bin in x non è nullo
       if (hx1->GetBinContent(k+1) > 0){
	 nbinsx1++;
	 for (int ky = 0; ky < hxy1->GetNbinsY(); ky++){
	 if (hxy1->GetBinContent(k+1,ky+1) > 0){
	   nbins1++;
	   ninbins1 += hxy1->GetBinContent(k+1,ky+1);
	 }
       }
	 c1_2->cd();
	 TLine *l = new TLine(hx1->GetBinLowEdge(k+1), -80, hx1->GetBinLowEdge(k+1), 80); //lo so, è memory leak che è un piacere, però io non le voglio più riprendere dopo il draw, quindi anche se non posso più accedere all'indirizzo di memoria fa nulla. Qualunque programmatore serio è libero di picchiarmi.
	 l->SetLineColor(kBlue);
	 l->Draw();
       }
       if (hy1->GetBinContent(k+1) > 0){
	 c1_2->cd();
	 TLine *l = new TLine(-100, hy1->GetBinLowEdge(k+1), 100,hy1->GetBinLowEdge(k+1));
	 l->SetLineColor(kRed);
	 l->Draw();
       }
     }
     for (int k = 0; k < hxy2->GetNbinsX(); k++){ //conto i bin con almeno una entry (strip accesa)
       if (hx2->GetBinContent(k+1) > 0){
	 nbinsx2++;
	 for (int ky = 0; ky < hxy2->GetNbinsY(); ky++){
	 if (hxy2->GetBinContent(k+1,ky+1) > 0){
	   nbins2++;
	   ninbins2 += hxy2->GetBinContent(k+1,ky+1);
	 }
       }
	 c2_2->cd();
	 TLine *l = new TLine(hx2->GetBinLowEdge(k+1), -80, hx2->GetBinLowEdge(k+1), 80); 
	 l->SetLineColor(kBlue);
	 l->Draw();
       }
       if (hy2->GetBinContent(k+1) > 0){
	 c2_2->cd();
	 TLine *l = new TLine(-100, hy2->GetBinLowEdge(k+1), 100,hy2->GetBinLowEdge(k+1));
	 l->SetLineColor(kRed);
	 l->Draw();
       }       
     }
     for (int k = 0; k < hxy3->GetNbinsX(); k++){ //conto i bin con almeno una entry (strip accesa)
       if (hx3->GetBinContent(k+1) > 0){
	 nbinsx3++;
	 for (int ky = 0; ky < hxy3->GetNbinsY(); ky++){
	 if (hxy3->GetBinContent(k+1,ky+1) > 0){
	   nbins3++;
	   ninbins3 += hxy3->GetBinContent(k+1,ky+1);
	 }
       }
	 c3_2->cd();
	 TLine *l = new TLine(hx3->GetBinLowEdge(k+1), -80, hx3->GetBinLowEdge(k+1), 80);
	 l->SetLineColor(kBlue);
	 l->Draw();
       }
       if (hy3->GetBinContent(k+1) > 0){
	 c3_2->cd();
	 TLine *l = new TLine(-100, hy3->GetBinLowEdge(k+1), 100,hy3->GetBinLowEdge(k+1));
	 l->SetLineColor(kRed);
	 l->Draw();
       }
     }
     
     for (int k = 0; k < hxy4->GetNbinsX(); k++){ //conto i bin con almeno una entry (strip accesa)
       if (hx4->GetBinContent(k+1) > 0){
	 nbinsx4++;
	 for (int ky = 0; ky < hxy4->GetNbinsY(); ky++){
	 if (hxy4->GetBinContent(k+1,ky+1) > 0){
	   nbins4++;
	   ninbins4 += hxy4->GetBinContent(k+1,ky+1);
	 }
       }
	 c4_2->cd();
	 TLine *l = new TLine(hx4->GetBinLowEdge(k+1), -80, hx4->GetBinLowEdge(k+1), 80);
	 l->SetLineColor(kBlue);
	 l->Draw();
       }
       if (hy4->GetBinContent(k+1) > 0){
	 c4_2->cd();
	 TLine *l = new TLine(-100, hy4->GetBinLowEdge(k+1), 100,hy4->GetBinLowEdge(k+1));
	 l->SetLineColor(kRed);
	 l->Draw();
       }
     }
     for (int k = 0; k < hxy5->GetNbinsX(); k++){ //conto i bin con almeno una entry (strip accesa)
       if (hx5->GetBinContent(k+1) > 0){
	 nbinsx5++;
	 for (int ky = 0; ky < hxy5->GetNbinsY(); ky++){
	   if (hxy5->GetBinContent(k+1,ky+1) > 0){
	   nbins5++;
	   ninbins5 += hxy5->GetBinContent(k+1,ky+1);
	 }
       }
	 c5_2->cd();
	 TLine *l = new TLine(hx5->GetBinLowEdge(k+1), -80, hx5->GetBinLowEdge(k+1), 80);
	 l->SetLineColor(kBlue);
	 l->Draw();
       }
       if (hy5->GetBinContent(k+1) > 0){
	 c5_2->cd();
	 TLine *l = new TLine(-100, hy5->GetBinLowEdge(k+1), 100,hy5->GetBinLowEdge(k+1));
	 l->SetLineColor(kRed);
	 l->Draw();
       }
     }

     for (int k = 0; k < hxy6->GetNbinsX(); k++){ //conto i bin con almeno una entry (strip accesa)
       if (hx6->GetBinContent(k+1) > 0){
	 nbinsx6++;
	 for (int ky = 0; ky < hxy6->GetNbinsY(); ky++){
	 if (hxy6->GetBinContent(k+1,ky+1) > 0){
	   nbins6++;
	   ninbins6 += hxy6->GetBinContent(k+1,ky+1);
	 }
       }
	 c6_2->cd();
	 TLine *l = new TLine(hx6->GetBinLowEdge(k+1), -80, hx6->GetBinLowEdge(k+1), 80); //lo so, è memory leak che è un piacere, però io non le voglio più riprendere dopo il draw, quindi anche se non posso più accedere all'indirizzo di memoria fa nulla. Qualunque programmatore serio è libero di picchiarmi.
	 l->SetLineColor(kBlue);
	 l->Draw();
       //delete l;
       }
       if (hy6->GetBinContent(k+1) > 0){
	 c6_2->cd();
	 TLine *l = new TLine(-100, hy6->GetBinLowEdge(k+1), 100,hy6->GetBinLowEdge(k+1));
	 l->SetLineColor(kRed);
	 l->Draw();
       }
     }
   }//fine ciclo sugli hit
  
   cout<<"Fine evento numero "<<i+1<<endl;
   c1_1->Print(Form("./displays/displays_centrale/evento%i/RPC_1.png",i+1));
   c1_2->Print(Form("./displays/displays_centrale/evento%i/RPC_1_2.png",i+1));
   c2_1->Print(Form("./displays/displays_centrale/evento%i/RPC_2.png",i+1));
   c2_2->Print(Form("./displays/displays_centrale/evento%i/RPC_2_2.png",i+1));
   c3_1->Print(Form("./displays/displays_centrale/evento%i/RPC_3.png",i+1));
   c3_2->Print(Form("./displays/displays_centrale/evento%i/RPC_3_2.png",i+1));
   c4_1->Print(Form("./displays/displays_centrale/evento%i/RPC_4.png",i+1));
   c4_2->Print(Form("./displays/displays_centrale/evento%i/RPC_4_2.png",i+1));
   c5_1->Print(Form("./displays/displays_centrale/evento%i/RPC_5.png",i+1));
   c5_2->Print(Form("./displays/displays_centrale/evento%i/RPC_5_2.png",i+1));
   c6_1->Print(Form("./displays/displays_centrale/evento%i/RPC_6.png",i+1));
   c6_2->Print(Form("./displays/displays_centrale/evento%i/RPC_6_2.png",i+1));

   c1_1->cd();
   hxy1->SetStats(false);
   hxy1->Draw("COLZ");
   hxy1->GetXaxis()->SetTitle("x[cm]");
   hxy1->GetYaxis()->SetTitle("y[cm]");
   /*c1_1->Update();
   TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
   stats->SetName("h1stats");
   stats->SetX1NDC(.5);
   stats->SetX2NDC(.7);*/
   TMarker * muonposition = new TMarker(muonx[0], muony[0], 5);
   muonposition->SetMarkerColor(kRed);
   muonposition->SetMarkerSize(0.5);
   muonposition->Draw();
   gSystem->Unlink("./displays_centrale/evento%i/event.gif"); //delete old file, otherwise it will add frames at the end
   c1_1->Print(Form("./displays/displays_centrale/evento%i/event.gif+50",i+1));
   c1_1->Print(Form("./displays/displays_centrale/evento%i/RPC_1.colz.png",i+1));
   c1_1->Print(Form("./displays/displays_centrale/evento%i/RPC_1.colz.root",i+1));
   c2_1->cd();
   hxy2->SetStats(false);
   hxy2->Draw("COLZ");
   hxy2->GetXaxis()->SetTitle("x[cm]");
   hxy2->GetYaxis()->SetTitle("y[cm]");
   /*c2_1->Update();
   TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
   stats->SetName("h1stats");
   stats->SetX1NDC(.5);
   stats->SetX2NDC(.7);*/
   muonposition = new TMarker(muonx[1], muony[1], 5);
   muonposition->SetMarkerColor(kRed);
   muonposition->SetMarkerSize(0.5);
   muonposition->Draw();
   c2_1->Print(Form("./displays/displays_centrale/evento%i/event.gif+50",i+1));
   c2_1->Print(Form("./displays/displays_centrale/evento%i/RPC_2.colz.png",i+1));
   c2_1->Print(Form("./displays/displays_centrale/evento%i/RPC_2.colz.root",i+1));
   c3_1->cd();
   hxy3->SetStats(false);
   hxy3->Draw("COLZ");
   hxy3->GetXaxis()->SetTitle("x[cm]");
   hxy3->GetYaxis()->SetTitle("y[cm]");
   /* c3_1->Update();
   TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
   stats->SetName("h1stats");
   stats->SetX1NDC(.5);
   stats->SetX2NDC(.7);*/
   muonposition = new TMarker(muonx[2], muony[2],5);
   muonposition->SetMarkerColor(kRed);
   muonposition->SetMarkerSize(0.5);
   muonposition->Draw();
   c3_1->Print(Form("./displays/displays_centrale/evento%i/event.gif+50",i+1));
   c3_1->Print(Form("./displays/displays_centrale/evento%i/RPC_3.colz.png",i+1));
   c3_1->Print(Form("./displays/displays_centrale/evento%i/RPC_3.colz.root",i+1));
   c4_1->cd();
   hxy4->SetStats(false);
   hxy4->Draw("COLZ");
   hxy4->GetXaxis()->SetTitle("x[cm]");
   hxy4->GetYaxis()->SetTitle("y[cm]");
   /*c4_1->Update();
   TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
   stats->SetName("h1stats");
   stats->SetX1NDC(.5);
   stats->SetX2NDC(.7);*/
   muonposition = new TMarker(muonx[3], muony[3],5);
   muonposition->SetMarkerColor(kRed);
   muonposition->SetMarkerSize(0.25);
   muonposition->Draw();
   c4_1->Print(Form("./displays/displays_centrale/evento%i/event.gif+50",i+1));
   c4_1->Print(Form("./displays/displays_centrale/evento%i/RPC_4.colz.png",i+1));
   c4_1->Print(Form("./displays/displays_centrale/evento%i/RPC_4.colz.root",i+1));
   c5_1->cd();
   hxy5->SetStats(false);
   hxy5->Draw("COLZ");
   hxy5->GetXaxis()->SetTitle("x[cm]");
   hxy5->GetYaxis()->SetTitle("y[cm]");
   /*c5_1->Update();
   TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
   stats->SetName("h1stats");
   stats->SetX1NDC(.5);
   stats->SetX2NDC(.7);*/
   muonposition = new TMarker(muonx[4], muony[4],5);
   muonposition->SetMarkerColor(kRed);
   muonposition->SetMarkerSize(0.5);
   muonposition->Draw();
   c5_1->Print(Form("./displays/displays_centrale/evento%i/event.gif+50",i+1));
   c5_1->Print(Form("./displays/displays_centrale/evento%i/RPC_5.colz.png",i+1));
   c5_1->Print(Form("./displays/displays_centrale/evento%i/RPC_5.colz.root",i+1));
   c6_1->cd();
   hxy6->SetStats(false);
   hxy6->Draw("COLZ");
   hxy6->GetXaxis()->SetTitle("x[cm]");
   hxy6->GetYaxis()->SetTitle("y[cm]");
   /*c6_1->Update();
   TPaveStats *stats = (TPaveStats*)gPad->GetPrimitive("stats");
   stats->SetName("h1stats");
   stats->SetX1NDC(.5);
   stats->SetX2NDC(.7);*/
   muonposition = new TMarker(muonx[5], muony[5],5);
   muonposition->SetMarkerColor(kRed);
   muonposition->SetMarkerSize(0.5);
   muonposition->Draw();
   c6_1->Print(Form("./displays/displays_centrale/evento%i/event.gif+50",i+1));
   c6_1->Print(Form("./displays/displays_centrale/evento%i/RPC_6.colz.png",i+1));
   c6_1->Print(Form("./displays/displays_centrale/evento%i/RPC_6.colz.root",i+1));
   //ricevo gli rms degli istogrammi
   rmsx1 += hxy1->GetRMS(1);
   rmsy1 += hxy1->GetRMS(2);
   rmsx2 += hxy2->GetRMS(1);
   rmsy2 += hxy2->GetRMS(2);
   rmsx3 += hxy3->GetRMS(1);
   rmsy3 += hxy3->GetRMS(2);
   rmsx4 += hxy4->GetRMS(1);
   rmsy4 += hxy4->GetRMS(2);
   rmsx5 += hxy5->GetRMS(1);
   rmsy5 += hxy5->GetRMS(2);
   rmsx6 += hxy6->GetRMS(1);
   rmsy6 += hxy6->GetRMS(2);

   nentries1 += hxy1->GetEntries();
   nentries2 += hxy2->GetEntries();
   nentries3 += hxy3->GetEntries();
   nentries4 += hxy4->GetEntries();
   nentries5 += hxy5->GetEntries();
   nentries6 += hxy6->GetEntries();
   
   delete c1_1;
   delete c2_1;
   delete c3_1;
   delete c4_1;
   delete c5_1;
   delete c6_1;
   
 }//fine ciclo sugli eventi
 cout<<endl;
 cout<<"RMS medio nel primo RPC lungo X: "<<(double) rmsx1/nevents<<" cm e lungo Y: "<<(double) rmsy1/nevents<< " cm "<<endl;
 cout<<"Numero medio di hits nei bin non nulli nel primo RPC: "<<(double) ninbins1/nbins1<<endl;
 cout<<"Numero medio di tracce nel primo piano: "<<(double) nentries1/nevents<<endl;
 cout<<endl;
 cout<<"RMS medio nel secondo RPC lungo X: "<<(double) rmsx2/nevents<<" cm e lungo Y: "<<(double) rmsy2/nevents<< " cm "<<endl;
 cout<<"Numero medio di hits nei bin non nulli nel secondo RPC: "<<(double) ninbins2/nbins2<<endl;
 cout<<"Numero medio di tracce nel secondo piano: "<<(double) nentries2/nevents<<endl;
 cout<<endl;
 cout<<"RMS medio nel terzo RPC lungo X: "<<(double) rmsx3/nevents<<" cm e lungo Y: "<<(double) rmsy3/nevents<< " cm "<<endl;
 cout<<"Numero medio di hits nei bin non nulli nel terzo RPC: "<<(double) ninbins3/nbins3<<endl;
 cout<<"Numero medio di tracce nel terzo piano: "<<(double) nentries3/nevents<<endl;
 cout<<endl;
 cout<<"RMS medio nel quarto RPC lungo X: "<<(double) rmsx4/nevents<<" cm e lungo Y: "<<(double) rmsy4/nevents<< " cm "<<endl;
 cout<<"Numero medio di hits nei bin non nulli nel quarto RPC: "<<(double) ninbins4/nbins4<<endl;
 cout<<"Numero medio di tracce nel quarto piano: "<<(double) nentries4/nevents<<endl;
 cout<<endl;
 cout<<"RMS medio nel quinto RPC lungo X: "<<(double) rmsx4/nevents<<" cm e lungo Y: "<<(double) rmsy5/nevents<< " cm "<<endl;
 cout<<"Numero medio di hits nei bin non nulli nel quinto RPC: "<<(double) ninbins5/nbins5<<endl;
 cout<<"Numero medio di tracce nel quinto piano: "<<(double) nentries5/nevents<<endl;
 cout<<endl;
 cout<<"RMS medio nel sesto RPC lungo X: "<<(double) rmsx6/nevents<<" cm e lungo Y: "<<(double) rmsy6/nevents<< " cm "<<endl;
 cout<<"Numero medio di hits nei bin non nulli nel sesto RPC: "<<(double) ninbins6/nbins6<<endl;
 cout<<"Numero medio di tracce nel sesto piano: "<<(double) nentries6/nevents<<endl;
 cout<<endl;
 
}

 bool ischarm(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 431) || (abs(PdgCode) == 411) || (abs(PdgCode) == 4122)  || (abs(PdgCode) == 421) || (abs(PdgCode) == 4132) || (abs(PdgCode) == 4232)||(abs(PdgCode) == 4332)||(abs(PdgCode) == 442))  check = true; 
      return check;
      }
//lista mesoni charmati: D0(bar), D+(-), Ds0(bar), Ds+(-), Lambdac+(-), Csi+(-), Csi0(bar),Omega_co(bar),eta_c
 
 bool isintermediate(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 223) || (abs(PdgCode) == 3224) || (abs(PdgCode) == 331) || (abs(PdgCode) == 221) || (abs(PdgCode) == 20213) || (abs(PdgCode) == 3212) ||(abs(PdgCode)==213) || (abs(PdgCode) == 113) || (abs(PdgCode) == 2224) || (abs(PdgCode) == 323)) check =  true;
  return check;
  }
 //lista intermedi: omega, sigma*, eta', eta, a1,pho+-,pho0, K*+-, delta++, sigma0


