//studio sulle efficienze di rivelazione di charm tenendo conto delle variabili misurabili sperimentalmente (creato il 1 Giugno 2016). Aggiornato il 22 Maggio 2017, per tenere conto del nuovo formato del file di Thomas (l'event header non corrisponde più al vertice primario, ma è posto sempre a (0,0,0), ora bisogna risalire al vertice dallo (StartX, StartY, StartZ) della traccia 0. 
bool ischarm(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);
#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/BoxPoint.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
//#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-13-02-ship-1/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

int efficienze(){
//void efficienze(TString filename = "", bool single = false){
//TFile *file1 = NULL;
 /*if (filename != "") file1 = TFile::Open(filename);
 else{
 filename = "./ship.10.0.Pythia8CharmOnly-TGeant4_1A10000events.root";
 cout<<"Apro il file di nome "<<filename<<endl;
 //TFile *file1 = TFile::Open("./studionoise/studio_precampomagnetico/ship.10.0.Pythia8CharmOnly-TGeant4_12cm.root");
 //TFile *file1 = TFile::Open("./ship.10.0.Pythia8CharmOnly-TGeant4_centralbeam_10000.root");
// TFile *file1 = TFile::Open("./ship.10.0.Pythia8CharmOnly-TGeant4_1A10000events.root");
 file1 = TFile::Open(filename);
 }*/
 TFile *file1 = TFile::Open("/home/utente/Simulations/sim_charmdet/prova/ship.conical.Pythia8CharmOnly-TGeant4.root");
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");

 //obtaining positions of target directly from geometry
 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import("/home/utente/Simulations/sim_charmdet/prova/geofile_full.conical.Pythia8CharmOnly-TGeant4.root"); 

 TGeoVolume *top = gGeoManager->GetTopVolume(); 
 TGeoNode *targetnode = top->GetNode("volTarget_1");
 if (targetnode == NULL){ cout<<"Error: Target volume not found"<<endl;
   return -1;
 }
 TGeoMatrix *transmatrix = targetnode->GetMatrix();
 Double_t zpostarget = transmatrix->GetTranslation()[2];
 Double_t zdimtarget = (((TGeoBBox*)targetnode->GetVolume()->GetShape())->GetDZ())*2.;
 cout<<"z position of center of target: "<<zpostarget<<"dim z: "<<zdimtarget<<endl;
 
 TGeoNode *firstemulsionnode = targetnode->GetVolume()->GetNode("Emulsion2_1"); //start of ECC (active part of target)
 TGeoNode *lastemulsionnode = targetnode->GetVolume()->GetNode("Emulsion_20"); //end of target
 if (firstemulsionnode == NULL || lastemulsionnode == NULL){
   cout<<"Error:almost one of the emulsion volumes not found"<<endl;
   return -1;
 }

 transmatrix = firstemulsionnode->GetMatrix();
 Double_t zstartactive = transmatrix->GetTranslation()[2] + zpostarget;
 transmatrix = lastemulsionnode->GetMatrix();
 Double_t zendactive = transmatrix->GetTranslation()[2] + zpostarget;
 cout<<"start of active part of target: "<<zstartactive<<" end of target: "<<zendactive<<endl;

 
 bool single = false; //se true calcola l'efficienza di single charm, altrimenti calcola l'efficienze di double charm
 if (single) cout<<"Single charm efficiency setted"<<endl;
 else cout<<"Double charm efficiency setted"<<endl;
 const Double32_t MolTh = 0.3;
 
 TDatabasePDG *pdg = TDatabasePDG::Instance();
 pdg->AddParticle("F0P0    ", " ", 0.9960, kFALSE, 0.0, 0, "meson",  9010221);//dalla pagina di AliRoot
 const Double_t kAu2Gev=0.9314943228;
 pdg->AddParticle("Deuteron","Deuteron",2*kAu2Gev+8.071e-3,kTRUE,0,3,"Ion",1000010020); //dal GitHub di FairRootGroup
 pdg->AddParticle("HE3","HE3",3*kAu2Gev+14.931e-3,kFALSE,0,6,"Ion",1000020030);
 //pdg->AddParticle("O17","O17",17 * 0.9314943228,kFALSE,0,8,"Ion",1000080170); //ossigeno????


 // const Double32_t check = 0; //fuori dalla prima lunghezza di interazione
 // const Double32_t check1 = -1.1; //ad almeno tre blocchi di Mo

 const Double32_t check = zendactive;
 const Double32_t check1 = zendactive - 3 * 0.3; 
 
 const int nrun = 1;

 
 Double32_t checkmin;
 switch(nrun){
 case 1: checkmin = -1000;//const Double32_t checkmin = -7.6;//per la seconda metà, se voglio tutto const Double32_t checkmin = -1000; 
   break;
 case 2: checkmin = -15.8; //posizione di inizio della porzione sensibile di target (run 2);
   break;
 case 3: checkmin = -14.6; //posizione di inizio della porzione sensibile di target (run 3);
   break;
 case 4: checkmin = -18.2; //posizione di inizio della porzione sensibile di target (run 4);
   break;
 case 5: checkmin = -14.6; //posizione di inizio della porzione sensibile di target (run 5);
   break;
 case 6: checkmin = -11.6; //posizione di inizio della porzione sensibile di target (run 6);
   break;
 case 7: checkmin = -9.8; //posizione di inizio della porzione sensibile di target (run 7);
   break;
 case 8: checkmin = -9.5; //posizione di inizio della porzione sensibile di target (run 8);
   break;
 case 9: checkmin = -9.7; //posizione di inizio della porzione sensibile di target (run 9);
   break;
 case 10:checkmin = -9.7; //posizione di inizio della porzione sensibile di target (run 10);
   break;
 }

 checkmin = zstartactive;
 
 TH1D *hkinklong = new TH1D ("hkinklong", "Angolo di kink figlie charm",100,0,0.5);
 TH1D *hkinkshort = new TH1D ("hkinkshort", "Angolo di kink figlie charm",100,0,0.5);
 TH1D * isto = new TH1D("isto", "Numero figlie di charm cariche",5,0,5);
 TH1D *hipexp = new TH1D ("hipexp", "Parametro di impatto particella figlia rispetto al vertice primario", 100, 0, 1);
 TH1D *hvis  = new TH1D("hvis", "Numero tracce cariche visibili", 30, 0, 30);
 TH1D *hpen = new TH1D ("hpen", "Numero tracce cariche penetranti", 30, 0, 30);
 TH1D *hcharmvis = new TH1D("hcharmvis", "Numero figlie charm cariche visibili", 5, 0, 5);
 TH1D *hl = new TH1D("hl", "Lunghezza volo charm", 100,0,7);

 TH2D *hxy = new TH2D("istoposizioni1", "Posizioni figlie di charm sul primo brick", 200, -10, 10, 200, -10, 10);
 	    
 FairMCEventHeader *primario = NULL;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.", &primario);
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);
 
 
 TClonesArray *arr2 = new TClonesArray("BoxPoint",10000);
 cbmsim->GetBranch("BoxPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("BoxPoint", &arr2);
 cbmsim->SetBranchStatus("*",1);

 Double32_t Vx = 0; //coordinate vertice primario
 Double32_t Vy = 0;
 Double32_t Vz = 0;

 bool preselezione;
 bool selezione;
 bool selezionecharm[2];
 const Int_t nevents = cbmsim->GetEntries();
 //const Int_t nevents = 1000; //prova
 Int_t ncut = 0; //numero di eventi tagliati perchè il vertice secondario è fuori dalla zona di check
 Int_t ntotalefigliecariche = 0; //numero totale di figlie cariche (non nel singolo evento)

 Int_t startpdg = 0; //pdg della particella figlia di charm
 Double_t startx = 0; //z del vertice secondario
 Double_t starty = 0;
 Double_t startz = 0;
 Double32_t px, py, pz;
 
 
 Int_t intermediatepdg = 0;
 Int_t mumID = 0;
 Int_t mumpdg = 0;
 Int_t hitID = 0;
 Int_t pdgcharm[2];
 Double_t charmpx[2];
 Double_t charmpy[2];
 Double_t charmpz[2];
 Int_t numerofigliecariche[2];
 Int_t figliecharmvisibili[2];
 Double_t pxtotfiglia[2];
 Double_t pytotfiglia[2];
 Double_t pztotfiglia[2];
 
 Double_t impulsototfiglia[2];
 Double_t impulsotrasversofiglia[2];
 Int_t ncaricapiu = 0; //numero di figlie di charm positive in cui deltathetax supera la soglia di 3 sigma per la misura di carica
 Int_t ncaricameno = 0; //numero di figlie di charm negative come sopra

 Int_t traccevisibili = 0;
 Int_t traccepenetranti = 0;

 Int_t nonelbrick = 0;
 Int_t nolocalizzati = 0;
 Int_t nodecay = 0;
 Int_t nok = 0;
 
 bool controllo[2]; //il check long-short va fatto una sola volta per charm (poiché un charm può avere più figlie e il controllo parte dalle figlie)

 Double32_t ipexp; 
 Double32_t ipmax[2]; //un massimo parametro di impatto e angolo di kink per ogni charm
 Double32_t delta;
 Double32_t impulso;
 Double32_t lunghezzavolo[2];
 bool checklunghezzavolo[2];
 
 Double32_t tanx;
 Double32_t tany;
 Double32_t tanxcharm[2];
 Double32_t tanycharm[2];
 Double32_t impulsocharm[2];
 
 Double32_t carica;
 Double32_t angolo;
 Double32_t angoloprimario;
 Double32_t angolocharm[2];
 Double32_t kinkminimo[2];

 Double32_t kinkmax[2];
 
 bool decay_check1[2];
 bool decay_check2[2];
 Int_t flong = 0;
 Int_t fshort = 0;

 bool off[2];
  Int_t nulli = 0;
 bool checklong[2]; //check if the charm is long or short;

 Double_t signalkink, signalip, signallength, signalmomentum, totalmomentum, totalpT; //variabili da salvare nel tree di output

 TFile *outtreefile = new TFile("signaltree.root","RECREATE");
 TTree *outtree1 = new TTree("oneprong", "signal variables 1 prong");
 TTree *outtree2 = new TTree("twoprong", "signal variables 2 prong");
 TTree *outtree3 = new TTree("threeprong", "signal variables 3 prong");
 
 //Branches list
 outtree1->Branch("kink", &signalkink, "kink/D");
 outtree1->Branch("ip", &signalip, "ip/D");
 outtree1->Branch("length", &signallength, "length/D");
 outtree1->Branch("pionmomentum", &signalmomentum, "pionmomentum/D");
 outtree1->Branch("totalmomentum", &totalmomentum, "totalmomentum/D");
 outtree1->Branch("totalpT", &totalpT, "totalpT/D");

 outtree2->Branch("kink", &signalkink, "kink/D");
 outtree2->Branch("ip", &signalip, "ip/D");
 outtree2->Branch("length", &signallength, "length/D");
 outtree2->Branch("pionmomentum", &signalmomentum, "pionmomentum/D");
 outtree2->Branch("totalmomentum", &totalmomentum, "totalmomentum/D");
 outtree2->Branch("totalpT", &totalpT, "totalpT/D");

 outtree3->Branch("kink", &signalkink, "kink/D");
 outtree3->Branch("ip", &signalip, "ip/D");
 outtree3->Branch("length", &signallength, "length/D");
 outtree3->Branch("pionmmomentum", &signalmomentum, "pionmomentum/D");
 outtree3->Branch("totalmomentum", &totalmomentum, "totalmomentum/D");
 outtree3->Branch("totalpT", &totalpT, "totalpT/D");
 cout<<"Numero di eventi: "<< nevents<<endl;
 bool isprimary;
 Int_t primaryin, primaryok;
 for (int i = 0; i < nevents; i++){
   isprimary = false;

   arr1->Clear();
   arr2->Clear();
   cbmsim->GetEntry(i);

   if (i%10 == 0) cout<<i<<endl;
   selezione = true;
   ShipMCTrack *trk = (ShipMCTrack*) arr1->At(0);
   Vx = trk->GetStartX();
   Vy = trk->GetStartY();
   Vz = trk->GetStartZ();
   if ((Vz > check) || (Vz < checkmin)){
     nulli++;
     continue;
   }
   if ((trk->GetPz() > 390.) && (trk->GetPdgCode() == 2212)){
     primaryin++;
     isprimary = true;
   }

   lunghezzavolo[0] = 0.;
   lunghezzavolo[1] = 0.;
   checklunghezzavolo[0] = false; //Fare il fill una sola volta per charm
   checklunghezzavolo[1] = false;
   traccevisibili = 0;
   traccepenetranti = 0;
   figliecharmvisibili[0] = 0;
   figliecharmvisibili[1] = 0;
   numerofigliecariche[0] = 0;
   numerofigliecariche[1] = 0;
   pxtotfiglia[0] = 0.;
   pytotfiglia[0] = 0.;
   pztotfiglia[0] = 0.;
   pxtotfiglia[1] = 0.;
   pytotfiglia[1] = 0.;
   pztotfiglia[1] = 0.;
   pdgcharm[0] = 0;
   pdgcharm[1] = 0;

   ipmax[0] = 0.;
   ipmax[1] = 0.;
   kinkminimo[0] = 0.;
   kinkminimo[1] = 0.;
   kinkmax[0] = 0.;
   kinkmax[1] = 0.;

   decay_check1[0] = false;
   decay_check1[1] = false;
   decay_check2[0] = false;
   decay_check2[1] = false;

   checklong[0] = false;
   checklong[1] = false;

   off[0] = false;
   off[1] = false;
   
   for (int j = 1; j < arr1->GetEntriesFast();j++){

     ShipMCTrack *trk = (ShipMCTrack*) arr1->At(j);
     if (ischarm(trk->GetPdgCode()) && (pdgcharm[0]==0)){  //memorizzo gli indici dei due charm
       pdgcharm[0] = trk->GetPdgCode();
      // charmpx[0] = trk->GetPx();
      // charmpy[0] = trk->GetPy();
      // charmpz[0] = trk->GetPz();
       tanxcharm[0] = trk->GetPx()/trk->GetPz(); 
       tanycharm[0] = trk->GetPy()/trk->GetPz();
       angolocharm[0] = TMath::ATan(TMath::Sqrt(pow(tanxcharm[0],2) + pow(tanycharm[0],2)));
       impulsocharm[0] = pow(pow(trk->GetPx(),2) + pow(trk->GetPy(),2) + pow(trk->GetPz(),2),0.5);
       controllo[0] = false; //azzero il controllo successivo sul decadimento long-short
     }
     if (ischarm(trk->GetPdgCode()) && (trk->GetPdgCode() != pdgcharm[0]) && (pdgcharm[1]==0)){ //il secondo lo memorizza solo se diverso dal primo
       pdgcharm[1] = trk->GetPdgCode();
      // charmpx[1] = trk->GetPx();
      // charmpy[1] = trk->GetPy();
      // charmpz[1] = trk->GetPz();
       tanxcharm[1] = trk->GetPx()/trk->GetPz(); 
       tanycharm[1] = trk->GetPy()/trk->GetPz();
       angolocharm[1] = TMath::ATan(TMath::Sqrt(pow(tanxcharm[1],2) + pow(tanycharm[1],2)));
       impulsocharm[1] = pow(pow(trk->GetPx(),2) + pow(trk->GetPy(),2) + pow(trk->GetPz(),2),0.5);
       controllo[1] = false;
     }

     if (trk->GetStartZ() == Vz){ //associo le particelle al vertice primario
	   tanx = trk->GetPx()/trk->GetPz();
	   tany = trk->GetPy()/trk->GetPz();
	  
	   impulso = TMath::Sqrt(pow(trk->GetPx(),2) + pow(trk->GetPy(),2) + pow(trk->GetPz(),2));
	    
	   angolo = TMath::ATan(TMath::Sqrt(pow(tanx,2) + pow(tany,2)));

	   carica = (pdg->GetParticle(trk->GetPdgCode())->Charge())/3.; //in unità di 'e'
	 	  
	   if ((angolo < 1.) && (impulso > 0.1) && (abs(carica)>0) && (ischarm(trk->GetPdgCode()) == false)) traccevisibili++; //visibile: angolo minore di 1 rad e impulso maggiore di 100 MeV
	   if ((angolo < 1.) && (impulso > 1.) && (abs(carica)>0) && (ischarm(trk->GetPdgCode()) == false)) traccepenetranti++; //penetranti: angolo minore di 1 rad e impulso maggiore di 1 GeV.
	   }


     
     mumID = trk->GetMotherId(); //prendo l'id della madre di ogni traccia (diretta, non madre della catena)
     mumpdg = trk->GetPdgCode();
     startpdg = trk->GetPdgCode();
     startz = trk->GetStartZ();
     startx = trk->GetStartX();
     starty = trk->GetStartY();
     px = trk->GetPx();
     py = trk->GetPy();
     pz = trk->GetPz();
    if (mumID>0){
       trk = (ShipMCTrack*)arr1->At(mumID);
       mumID = trk->GetMotherId();
       mumpdg = trk->GetPdgCode();
      while ((isintermediate(mumpdg) == true) && (mumID>0)) //se la madre è uno stato intermedio controllo se lo stato intermedio è stato prodotto dal charm
	 {
	   intermediatepdg = mumpdg;
	   
	   trk = (ShipMCTrack*) arr1->At(mumID);
	   mumID = trk->GetMotherId();
	   mumpdg = trk->GetPdgCode();
	 }
      //salto il passaggio in cui lo stato intermedio è ricondotto alla madre, avendo contato i suoi prodotti
      for (int charmindex = 0; charmindex < 2; charmindex++){
	if ((mumpdg==pdgcharm[charmindex]) && (isintermediate(startpdg) == false)){
	  lunghezzavolo[charmindex] = TMath::Sqrt(pow((startz-Vz),2) + pow((startx-Vx),2) + pow((starty-Vy),2));

	  if ((lunghezzavolo[0] != 0.) &&(checklunghezzavolo[0]==false)){
	    hl->Fill(lunghezzavolo[0]);
	    checklunghezzavolo[0] = true;
	  }
	  if ((lunghezzavolo[1] != 0.) && (checklunghezzavolo[1]==false)){
	      hl->Fill(lunghezzavolo[1]);
	      checklunghezzavolo[1] = true;
	    }
	}
	if ((mumpdg==pdgcharm[charmindex]) && (isintermediate(startpdg) == false)){
	  /*if (startpdg > 100000){
	    cout<<"Frammento figlio di charm! "<<startpdg<<" nell'evento "<<i<<endl;
	    continue; //possibili particelle non riconosciute (studiare meglio, comunque eventuali frammenti non possono essere inclusi nella ricerca delle figlie di charm)
	  }*/
	  if (abs((pdg->GetParticle(startpdg))->Charge()) > 0.){ //prima chiedo se è figlia di charm, poi accedo al pdg per la carica (diminuisco la probabilità di trovare particelle 'sconosciute' con conseguente crash)
	    carica = (pdg->GetParticle(startpdg)->Charge())/3.; //in unità di 'e'
	  
	  if (off[charmindex] == true) continue;
	  if (startz > check1){ //controllo la posizione del vertice secondario, è all'interno della regione di controllo?
	 off[charmindex] = true;
	 continue;
       }
	  int ntagged = 0; //per ogni figlia di charm conto quanti piani sono stati attraversati

	  //taglio il ciclo sui box points per ottenere le altre distribuzioni nel modo più rapido possibile
	  //kinkminimo[charmindex] = 100.;
	  //  ipexp = 100.;
	  delta = 0.;
	  ipexp = 0.;
	  kinkminimo[charmindex] = 0;
          for (int l = 0; l <arr2->GetEntriesFast();l++){
                   BoxPoint *target = (BoxPoint*) arr2->At(l);
		     if (target->GetTrackID() == j){		       
		       hxy->Fill(target->GetX(), target->GetY());
		       impulso = TMath::Sqrt(pow(target->GetPx(),2) + pow(target->GetPy(),2) + pow(target->GetPz(),2)); //parte del parametro di impatto
		         delta = 0.;
		         ipexp = 0.;
		         delta += (Vx- target->GetX()) * target->GetPx()/impulso;
			 delta += (Vy - target->GetY()) * target->GetPy()/impulso;
			 delta += (Vz - target->GetZ()) * target->GetPz()/impulso;
			 ipexp += pow((Vx - target->GetX() - delta * target->GetPx()/impulso),2);
			 ipexp += pow((Vy - target->GetY() - delta * target->GetPy()/impulso),2);
			 ipexp += pow((Vz - target->GetZ() - delta * target->GetPz()/impulso),2);
			 ipexp = TMath::Sqrt(ipexp);
			 if (ipexp > ipmax[charmindex]) ipmax[charmindex] = ipexp;
			 //hipexp->Fill(ipexp*10); //passo da centimetri a millimetri

			 
			
			 
			 //	 if (controllo[charmindex] == false){ //parte del check long-short
			 if ((target->GetZ() - Vz) > MolTh){
			     if (checklong[charmindex] == false) flong++;
			     checklong[charmindex] = true;
			     //tanx = target->GetPx()/target->GetPz(); //parte dell'angolo di kink minimo (ipotesi long)
			     //tany = target->GetPy()/target->GetPz();
			     tanx = px/pz; //parte dell'angolo di kink minimo (ipotesi long)
			     tany = py/pz;
			     angolo = TMath::ATan(TMath::Sqrt(pow(tanx,2) + pow(tany,2)));

			     tanx = tanx - tanxcharm[charmindex];
			     tany = tany - tanycharm[charmindex];

			     //Double_t kink = TMath::ATan(TMath::Sqrt((tanx)**2 + (tany)**2));
			     kinkminimo[charmindex] = TMath::ATan(TMath::Sqrt(pow((tanx),2) + pow((tany),2)));
			      
			     //kinkminimo = TMath::Abs(angolo - angolocharm[charmindex]); //se il charm lascia una traccia, uso l'impulso del charm
			     if (kinkminimo[charmindex] > kinkmax[charmindex]) kinkmax[charmindex] = kinkminimo[charmindex];
			      //hkinklong->Fill(kinkminimo);
			     controllo[charmindex] = true;
			   }
			   else{
			     if (checklong[charmindex] == false) fshort++;
			     checklong[charmindex] = true;
			     tanx = px/pz; //parte dell'angolo di kink minimo (ipotesi short)
			     tany = py/pz;
			     angolo = TMath::ATan(TMath::Sqrt(pow(tanx,2) + pow(tany,2)));
			     // cout<<"Evento short, angoli effettivi: "<<tanxcharm[charmindex]<<" "<<tanycharm[charmindex]<<endl;
			     tanxcharm[charmindex] = (target->GetX() - Vx)/(target->GetZ()-Vz);
			     tanycharm[charmindex] = (target->GetY() - Vy)/(target->GetZ()-Vz);
			     //cout<<"Evento short, angoli stimati: "<<tanxcharm[charmindex]<<" "<<tanycharm[charmindex]<<endl;
			     angoloprimario = TMath::ATan(TMath::Sqrt(pow(tanxcharm[charmindex],2) + pow(tanycharm[charmindex],2)));

			     tanx = tanx - tanxcharm[charmindex];
			     tany = tany - tanycharm[charmindex];

			     // Double_t kink = TMath::ATan(TMath::Sqrt((tanx)**2 + (tany)**2));
			     
			     kinkminimo[charmindex] = TMath::ATan(TMath::Sqrt(pow((tanx),2) + pow((tany),2)));
			     
			     if (kinkminimo[charmindex] > kinkmax[charmindex]) kinkmax[charmindex] = kinkminimo[charmindex];
			     // if (kinkminimo[charmindex] > kinkmax[charmindex]) kinkmax[charmindex] = kinkminimo;
			     //kinkminimo = TMath::Abs(angolo - angoloprimario);
			     //hkinkshort->Fill(kinkminimo);
			     controllo[charmindex] = true;
			       }
			   // }
			 break;//prendo solo il primo hit
		       }
		   
		     
	       }//fine ciclo sul box point
	  if ((angolo < 1.) && (impulso > 0.1) && (TMath::Abs(carica)>0)){
		 figliecharmvisibili[charmindex]++;
		 pxtotfiglia[charmindex] += px;
		 pytotfiglia[charmindex] += py;
		 pztotfiglia[charmindex] += pz;
		 
	       }//else cout<<angolo<<" "<<impulso<<" "<<(TMath::Abs(carica))<<endl;
               if (kinkminimo[charmindex] > 0.02) decay_check1[charmindex] = true;
	       if (ipexp > 0.001) decay_check2[charmindex] = true;	      
	  }
	  
	
	}//fine ciclo sui due charm
      }
    }
    }//fine ciclo sulle tracce
   //calcolo impulso totale e impulso trasverso totale
   impulsototfiglia[0] = pow((pow(pxtotfiglia[0],2) + pow(pytotfiglia[0],2) + pow(pztotfiglia[0],2)),0.5);
   impulsototfiglia[1] = pow((pow(pxtotfiglia[1],2) + pow(pytotfiglia[1],2) + pow(pztotfiglia[1],2)),0.5);
   
   impulsotrasversofiglia[0] = impulsototfiglia[0] * TMath::Sin(TMath::ACos((charmpx[0] * pxtotfiglia[0] + charmpy[0] * pytotfiglia[0] + charmpz[0] * pztotfiglia[0])/(impulsocharm[0] * impulsototfiglia[0])));
   impulsotrasversofiglia[1] = impulsototfiglia[1] * TMath::Sin(TMath::ACos((charmpx[1] * pxtotfiglia[1] + charmpy[1] * pytotfiglia[1] + charmpz[1] * pztotfiglia[1])/(impulsocharm[1] * impulsototfiglia[1])));

   if (single){
      for (int charmindex = 0; charmindex<2; charmindex++){
     preselezione = false;
     if (off[charmindex] == false) preselezione = true;
     if (preselezione == true){
       selezione = false;
       if ((traccevisibili >=2) && (traccepenetranti >=1)) selezione = true; //localizzazione
       if (selezione == true){

	 //salvo gli eventi nel tree
	 signallength = lunghezzavolo[charmindex];
	 signalip = ipmax[charmindex];
	 signalkink = kinkmax[charmindex];
	 signalmomentum = impulsocharm[charmindex];
	 totalmomentum = impulsototfiglia[charmindex];
	 totalpT = impulsotrasversofiglia[charmindex];
	 // cout<<kinkmax[charmindex]<<endl;

	 // if (figliecharmvisibili[charmindex] > 0 && totalpT == 0) cout<<(charmpx[1] * pxtotfiglia[1] + charmpy[1] * pytotfiglia[1] + charmpz[1] * pztotfiglia[1])/(impulsocharm[1] * impulsototfiglia[1])<<endl;
	 if (figliecharmvisibili[charmindex] == 1) outtree1->Fill();
	 if (figliecharmvisibili[charmindex] == 2) outtree2->Fill();
	 if (figliecharmvisibili[charmindex] == 3) outtree3->Fill();	    	
	 if (ipmax[charmindex] != 0) hipexp->Fill(ipmax[charmindex]*10); //la condizione di 0 serve a evitare fill nel caso di nessuna figlia di charm carica.
	 if (kinkmax[charmindex] != 0) hkinklong->Fill(kinkmax[charmindex]);
	 
	 selezionecharm[charmindex] = false;
	 
	 if ((decay_check1[charmindex] == true) &&  (decay_check2[charmindex] == true) && (figliecharmvisibili[charmindex] >= 1) && (lunghezzavolo[charmindex]<0.6)) selezionecharm[charmindex] = true;
	 if (selezionecharm[charmindex] == true) nok++;
	 else nodecay++;
       }else nolocalizzati++;
     }else nonelbrick++;
   }
   if (i%100 == 0) cout <<i<<endl;
 }
 else{
   preselezione = false;
   //cout<<pdgcharm[0]<<" "<<pdgcharm[1]<<endl;
   if ((off[0] == false) && (off[1] == false)) preselezione = true; //decadimenti nel brick
   if (preselezione == true){
     selezione = false;
     if ((traccevisibili >=2) && (traccepenetranti >=1)) selezione = true; //localizzazione vertice primario
     if (selezione == true){
   for (int charmindex = 0; charmindex<2; charmindex++){

     //salvo gli eventi nel tree
     signallength = lunghezzavolo[charmindex];
     signalip = ipmax[charmindex];
     signalkink = kinkmax[charmindex];
     signalmomentum = impulsocharm[charmindex];
     totalmomentum = impulsototfiglia[charmindex];
     totalpT = impulsotrasversofiglia[charmindex];	
     if (figliecharmvisibili[charmindex] == 1) outtree1->Fill();
     if (figliecharmvisibili[charmindex] == 2) outtree2->Fill();
     if (figliecharmvisibili[charmindex] == 3) outtree3->Fill();	       
     hcharmvis->Fill(figliecharmvisibili[charmindex]);
     selezionecharm[charmindex] = false;
     if ((decay_check1[charmindex] == true) && (decay_check2[charmindex] == true) && (figliecharmvisibili[charmindex] >= 1) && (lunghezzavolo[charmindex]<0.6)) selezionecharm[charmindex] = true; //decaysearch
   }
   
   if ((selezionecharm[0] == true) && (selezionecharm[1] == true)){
     nok++;
     if (isprimary == true) primaryok++;
   }
   else nodecay++;
 }else nolocalizzati++;
   
 }else nonelbrick++;
 }
}//fine ciclo sugli eventi
 cout<<endl;
 cout<<"Eventi long: "<<flong<<endl;
 cout<<"Eventi short: "<<fshort<<endl;
 cout<<"Protoni interagiuti fuori dal target" <<nulli<<endl;
 cout<<"No nel brick: "<<nonelbrick<<endl;
 cout<<"No localizzati: "<< nolocalizzati<<endl;
 cout<<"Nodecay: "<<nodecay<<endl;
 cout<<"Passano tutti i check" << nok<<endl;
 if (single) cout<<"Efficienze di single:" <<(double)nok/(nevents*2 - nulli*2)<<endl;
 else cout<<"Efficienza di double: "<< (double)nok/(nevents - nulli)<<endl;
 cout<<"Numero eventi primari che passano la selezione: "<<primaryok<<"su "<<primaryin<<endl;
 TCanvas *c1 = new TCanvas();
 hkinklong->Draw();
 hkinklong->SetTitle("");
 hkinklong->GetXaxis()->SetTitle("rad");
 c1->Print("./temporanei/efficienze/kink_angle.root","root");
 TCanvas *cip = new TCanvas();
 hipexp->Draw();
 hipexp->SetTitle("");
 hipexp->GetXaxis()->SetTitle("mm");
 cip->Print("./temporanei/efficienze/impact_parameter.root","root");
 TCanvas *c2 = new TCanvas();
 hl->Draw();
 TCanvas *c3 = new TCanvas();
 hcharmvis->Draw();
 TCanvas *c4 = new TCanvas();
 hxy->Draw();


 outtreefile->Write();
 outtreefile->Close();

 return 0;
}

bool ischarm(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 431) || (abs(PdgCode) == 411) || (abs(PdgCode) == 4122)  || (abs(PdgCode) == 421) || (abs(PdgCode) == 4132) || (abs(PdgCode) == 4232) ||(abs(PdgCode) == 4332))  check = true; 
      return check;
      }
//lista mesoni charmati: D0(bar), D+(-), Ds0(bar), Ds+(-), Lambdac+(-), Csi+(-), Csi0(bar), Omega_co(bar)

bool isintermediate(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 223) || (abs(PdgCode) == 3224) || (abs(PdgCode) == 331) || (abs(PdgCode) == 221) || (abs(PdgCode) == 20213) || (abs(PdgCode) == 3212) ||(abs(PdgCode)==213) || (abs(PdgCode) == 113) || (abs(PdgCode) == 2224) || (abs(PdgCode) == 323)) check =  true;
  return check;
  }

//lista intermedi: omega, sigma*, eta', eta, a1,pho+-,pho0, K*+-, delta++, sigma0 (includere (abs(PdgCode) == 3222)?)
