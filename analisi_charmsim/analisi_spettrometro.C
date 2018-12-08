//studio sui prodotti di decadimento dei charm nello spettrometro (creato il 13 Maggio, nuova versione con interazioni di protoni). T1 e T2 non sono più presenti, sostituiti dai 12 piani di PIXEL. La stima dell'angolo con la differenza in posizione non ha più senso, per il momento usare gli impulsi degli spectrometerpoint.
#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/MufluxSpectrometerPoint.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-13-02-ship-1/include/TDatabasePDG.h"
//#include "/home/utente/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

#include "TGeoManager.h"
#include "TGeoVolume.h"
#include "TGeoNode.h"
#include "TGeoMatrix.h"
bool ischarm(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);

int analisi_spettrometro(){
// TFile *file1 = TFile::Open("./ship.10.0.Pythia8CharmOnly-TGeant4_siliconTT.root");
 // TFile *file1 = TFile::Open("./ship.10.0.Pythia8CharmOnly-TGeant4_checkmagneticfield.root");
  //TFile *file1 = TFile::Open("./10run/ship.10.0.Pythia8CharmOnly-TGeant4_secondo.root");
// file1 = TFile::Open("/home/utente/SHIPBuild/sim_charm/OPERAECC/muontaggeracceptance/ship.conical.Pythia8CharmOnly-TGeant4.root");
 //TFile *file1 = TFile::Open("spectrohits_EPFL.root");
 TFile *file1 = TFile::Open("/home/utente/Simulations/sim_charmdet/charm_events/lowerfieldmap/ship.conical.Pythia8CharmOnly-TGeant4.root");
 TTree* cbmsim = (TTree*)file1->Get("cbmsim"); 

 //ottenere posizione e dimensione bersaglio di neutrini tau
 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import("/home/utente/Simulations/sim_charmdet/charm_events/lowerfieldmap/geofile_full.conical.Pythia8CharmOnly-TGeant4.root"); 

 TGeoVolume *top = gGeoManager->GetTopVolume(); 
 TGeoNode *targetnode = top->GetNode("volTarget_1");
 if (targetnode == NULL){ cout<<"Error: Target volume not found"<<endl;
   return -1;
 }
 TGeoMatrix *transmatrix = targetnode->GetMatrix();
 Double_t zpostarget = transmatrix->GetTranslation()[2];
 Double_t zdimtarget = (((TGeoBBox*)targetnode->GetVolume()->GetShape())->GetDZ())*2.;
 cout<<"z position of center of target: "<<zpostarget<<"dim z: "<<zdimtarget<<endl;

 /*TGeoNode *firstemulsionnode = targetnode->GetVolume()->GetNode("Emulsion2_1"); //start of ECC (active part of target)
 TGeoNode *lastemulsionnode = targetnode->GetVolume()->GetNode("Emulsion_20"); //end of target
 if (firstemulsionnode == NULL || lastemulsionnode == NULL){
   cout<<"Error:almost one of the emulsion volumes not found"<<endl;
   return -1;
 }*/
 //transmatrix = firstemulsionnode->GetMatrix();
 //Double_t zstartactive = transmatrix->GetTranslation()[2] + zpostarget;
 //transmatrix = lastemulsionnode->GetMatrix();
 //Double_t zendactive = transmatrix->GetTranslation()[2] + zpostarget;
 Double_t zendactive = zpostarget + zdimtarget/2.;
 Double_t zstartactive = zendactive - 5.6;
 cout<<"start of active part of target: "<<zstartactive<<" end of target: "<<zendactive<<endl;

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
 
 //DEFININING HISTOGRAMS
 TH1D *isto = new TH1D("isto","Numeri di prodotti carichi del charm", 7, 0, 7);
 TH1D *hl = new TH1D("hl", "Lunghezza volo charm", 100,0,7);
 TH1D* hz = new TH1D("hz", "Posizione in z tracce figlie di charm", 100,-22,-2);
 TH1D* hz1 = new TH1D("hz1", "Posizione in z tracce figlie di charm non arrivate all'ultimo", 100,-22,-2);
 
 TH1D *hefflocal = new TH1D("hefflocal", "Accettanza Goliath in funzione dell'impulso", 20, 0, 100);
 TH1D *hp = new TH1D("hp", "Impulso muoni figli di charm che arrivano alla fine dello spettrometro", 20, 0, 100);
 TH1D *hpin = new TH1D("hpin", "Impulso muoni figli di charm entranti", 20, 0, 100);
 hp->Sumw2(); //voglio usarli per un calcolo di efficienza, è fondamentale salvarsi gli errori affinchè il calcolo sia fatto bene (Divide).
 hpin->Sumw2();

 TH1D *hangle = new TH1D("hangle", "Angolo figlie di charm che arrivano alla fine dello spettrometro", 100, 0, 1);
 TH1D *hanglewrong = new TH1D("hanglewrong", "Angolo figlie di charm che non arrivano alla fine dello spettrometro", 100, 0, 1);
 TH1D *istodeltaxpiu = new TH1D("istodeltaxpiu","Differenze fra le componenti x degli angoli", 100, -0.8, 0.8);
 TH1D *istodeltaxmeno = new TH1D("istodeltaxmeno","Differenze fra le componenti x degli angoli", 100, -0.8, 0.8);
 TH1D *istodeltax = new TH1D("istodeltax","Differenze fra le componenti x degli angoli", 100, -0.9, 0.9);
 TH1D *istodeltay = new TH1D("istodeltay","Differenze fra le componenti y degli angoli", 100, -0.03, 0.03);
 TH1D *istodp = new TH1D("istodp", "Deltap su p", 100, -1, 1);
  
 TH1I *hnpixel = new TH1I("hnpixel", "Number of transversed pixel planes", 7,0,7);
 
 TH2D *istoposizioni1 = new TH2D("istoposizioni1", "Posizioni figlie di charm sul primo piano dello spettrometro", 40, -2, 2, 40, -2, 2);
 TH2D *istoposizioni2 = new TH2D("istoposizioni2", "Posizioni figlie di charm sul secondo piano dello spettrometro", 40, -2, 2, 40, -2, 2);
 TH2D *istoposizioni3 = new TH2D("istoposizioni3", "Posizioni figlie di charm sul terzo piano dello spettrometro", 200, -100, 100, 200, -100, 100);
 TH2D *istoposizioni4 = new TH2D("istoposizioni4", "Posizioni figlie di charm sul quarto piano dello spettrometro", 200, -100, 100, 200, -100, 100);


 // TObjArray * branchTrackID = (TObjArray*) cbmsim->GetBranch("BoxPoint.fUniqueID");
 Int_t nevents = cbmsim->GetEntries(); //nota: le entries di tree sono il doppio ci quelle di cbmsim (ogni evento->due charm)
 Int_t ntracks;


// const Double_t check = 0; //fuori dal target prima lunghezza interazione
// const Double_t check1 = check; //charm che decadono dopo il target non li devo considerare

  const Double_t check = zendactive;
  const Double_t check1 = check;

  Double_t checkmin = zstartactive;  
  //const Double_t checkmin = -6.0;//ultimi 60 cm
 //const Double_t checkmin = -15.8; //posizione di inizio della porzione sensibile di target (run 2);
 //const Double_t checkmin = -19.5; //posizione di inizio della porzione sensibile di target (run 2);
 //const Double_t checkmin = -14.6; //posizione di inizio della porzione sensibile di target (run 3);
 //const Double_t checkmin = -18.2; // (run 4)
 //const Double_t checkmin = -14.6; //posizione di inizio della porzione sensibile di target (run 5);

 Int_t ncut = 0;
 
 FairMCEventHeader *prova = NULL;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.", &prova);
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);
 
 
 TClonesArray *arr2 = new TClonesArray("SpectrometerPoint",1000);
 cbmsim->GetBranch("SpectrometerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("SpectrometerPoint", &arr2);

 TClonesArray *arr3 = new TClonesArray("MufluxSpectrometerPoint",1000);
 cbmsim->GetBranch("MufluxSpectrometerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MufluxSpectrometerPoint", &arr3);
 /*
 TClonesArray *arr3 = new TClonesArray("BoxPoint",10000);
 cbmsim->GetBranch("BoxPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("BoxPoint", &arr3);*/
 
 cbmsim->SetBranchStatus("*",1);
 Double_t px = 0;
 Double_t py = 0;
 Double_t pz = 0;
 Double_t Vz = 0;
 Double_t Vy = 0;
 Double_t Vx = 0;
 Double_t startz, starty, startx; //posizioni vertici secondari

 Double_t SpectrometerPX[4];
 Double_t SpectrometerPY[4];
 Double_t SpectrometerPZ[4];
 Double_t SpectrometerX[4]; //posizioni hit nello spettrometro
 Double_t SpectrometerY[4];
 Double_t SpectrometerZ[4];
 Double_t Theta1x,Theta2x,Theta1y,Theta2y; //angoli fra i primi due piani e i secondi due piani, componente x e componente y;
 Double_t Deltathetax, Deltathetay;
 
 bool off[2]; //il charm è decaduto fuori dalla regione di check?
 bool selezione;
 bool lunghezzavolocharm[2]; //la lunghezza di volo del charm va calcolata una volta per evento

 Int_t tutteoutspettrometrosingle = 0;
 Int_t tutteoutspettrometrodouble = 0;
 Int_t tutteinspettrometrosingle = 0;
 Int_t tutteinspettrometrodouble = 0;

 Int_t ntotalefigliecariche = 0; //numero totale di figlie cariche (non nel singolo evento)
 Int_t rivelate = 0;
 Int_t rivelate2 = 0;
 Int_t rivelate3 = 0;
 Int_t uscite = 0;
 Int_t startpdg = 0;
 Int_t intermediatepdg = 0;
 Int_t mumID = 0;
 Int_t mumpdg = 0;
 Int_t hitID = 0;
 Int_t numerofigliecariche[2];
 Int_t nfiglieinspettrometro[2];
 Int_t nfiglieoutspettrometro[2];
 Int_t pdgcharm[2];
 Int_t ncaricapiu = 0; //numero di figlie di charm positive in cui deltathetax supera la soglia di 3 sigma per la misura di carica
 Int_t ncaricameno = 0; //numero di figlie di charm negative come sopra
 Int_t nmisid = 0;
 Double_t momentum;

 Int_t scatteringID;  //variabili per il controllo dell'inefficienza
 Int_t scatteringmumID;
 Int_t suspectID;
 Int_t nemulsions; //numero di emulsioni attraversate dalle figlie di charm
 Double_t angolox;
 Double_t angoloy;
 Double_t proiezionex;
 Double_t proiezioney;
 Double_t zmax; //massimo z rivelato nel target
 Int_t nulli = 0;

 TRandom3 *g  = new TRandom3(); //generatore casuale per la risoluzione.
 const Double_t resolution = 0.01;
 TF1 * f1 = new TF1("f1", "gaus", -0.003, 0.003); //gaussiana per il fit

 const Double_t sigmatheta = 0.0012; //risoluzione ottenuta dal grafico sulle y


 TGraph *cali = new TGraph(); //curva di calibrazione per l'impulso

 TF1 * f2 = new TF1("f2", "[0] + [1] * x", 0,60);

 //f2->SetParameters(0.502,0.511); //vecchia configurazione
 f2->SetParameters(0.762,0.496);
 f2->SetParameters(0.15,0.53);
 TF1 * f2mono = new TF1("f2mono", "[0] * x",0,60);
 f2mono->SetParameter(0,0.536);
 TF1 * f3 = new TF1("f3", "gaus", -0.2, 0.4);

 bool isprimary;
 Int_t primaryok = 0;

 Double_t impulso1;
 Double_t impulso2;
 Double_t impulso3;
 Double_t impulso4;

 Int_t npixel; //piani di pixel attraversati
 bool transversedstation[6];
 for (int i = 0; i < nevents; i++){

   arr1->Clear();
   arr2->Clear();
   cbmsim->GetEntry(i);

   if (i% 100 == 0) cout<<i<<endl;
   selezione = true;
   isprimary = false;
   ShipMCTrack *track = (ShipMCTrack*) arr1->At(0);
   Vx = track->GetStartX();
   Vy = track->GetStartY();
   Vz = track->GetStartZ();
  // cout<<i<<" "<<Vx<<" "<<Vy<<" "<<Vz<<endl;
   /*Vz = prova->GetZ();
   Vx = prova->GetX();
   Vy = prova->GetY();*/
   if ((Vz > check) || (Vz < checkmin)){
     nulli++;
     continue;
     }
   if ((track->GetPz() > 390.) && (track->GetPdgCode() == 2212)) isprimary = true;
   numerofigliecariche[0] = 0;
   numerofigliecariche[1] = 0;
   nfiglieinspettrometro[0] = 0;
   nfiglieinspettrometro[1] = 0;
   nfiglieoutspettrometro[0] = 0;
   nfiglieoutspettrometro[1] = 0;
   pdgcharm[0] = 0;
   pdgcharm[1] = 0;
   off[0] = false;
   off[1] = false;
   lunghezzavolocharm[0] = false;
   lunghezzavolocharm[1] = false;      
   
   for (int j = 1; j < arr1->GetEntriesFast();j++){
     ShipMCTrack *trk = (ShipMCTrack*) arr1->At(j);
     //if (trk->GetPdgCode() == 9010221) continue; //questa particella non è contenuta nel TDatabasePDG di root
     momentum = TMath::Sqrt(pow(trk->GetPx(),2)+ pow(trk->GetPy(),2)+pow(trk->GetPz(),2));

     if (ischarm(trk->GetPdgCode()) && (pdgcharm[0]==0)){  //memorizzo gli indici dei due charm
       pdgcharm[0] = trk->GetPdgCode();
  
     }
     if (ischarm(trk->GetPdgCode()) && (trk->GetPdgCode() != pdgcharm[0]) && (pdgcharm[1]==0)){ //il secondo lo memorizza solo se diverso dal primo
       pdgcharm[1] = trk->GetPdgCode();
  
     }
     mumID = trk->GetMotherId(); //prendo l'id della madre di ogni traccia (diretta, non madre della catena)
     mumpdg = trk->GetPdgCode();
     startpdg = trk->GetPdgCode();
     startz = trk->GetStartZ();
     startx = trk->GetStartX();
     starty = trk->GetStartY();
     pz = trk->GetPz();
     px = trk->GetPx();
     py = trk->GetPy();
     
    
    if (mumID>0){
       trk = (ShipMCTrack*)arr1->At(mumID);
       mumpdg = trk->GetPdgCode();
       
       //if (i == 26) cout <<"ID traccia: "<<j<<" di pdg "<<startpdg<<" figlia di "<<mumpdg<<" id madre "<<mumID<<" con impulso "<<momentum<<endl;
       //if ((i == 26) && (mumID == 10)) cout<<"Bingo: "<<startz<<endl;
       
       while ((isintermediate(mumpdg) == true) && (mumID>0)) //se la madre è uno stato intermedio controllo se lo stato intermedio è stato prodotto dal charm
	 {
	   intermediatepdg = mumpdg;
	   mumID = trk->GetMotherId();
	   if (mumID <=0 ) break; //evita crash se mumID <0
	   trk = (ShipMCTrack*) arr1->At(mumID);
	  
	   mumpdg = trk->GetPdgCode();
	   }
      //salto il passaggio in cui lo stato intermedio è ricondotto alla madre, avendo contato i suoi prodotti
       
        for (int charmindex = 0; charmindex < 2; charmindex++){
	  //figlie di charm
        //if (true){
	if ((mumpdg==pdgcharm[charmindex]) && (isintermediate(startpdg) == false) && (abs((pdg->GetParticle(startpdg))->Charge()) > 0.) ){	
	   if (startz > check1){ //controllo la posizione del vertice secondario, è all'interno della regione di controllo? Se sì, aumenta il contatore ncut
	     if (off[charmindex] == false) ncut++; //la somma la faccio una volta per charm
	     off[charmindex] = true;
	  } //questa parte va adattata alla posizione del bersaglio finale (nota del 7 ottobre)

	   numerofigliecariche[charmindex]++; //questo lo uso per il numero di figlie cariche di ogni charm
	   ntotalefigliecariche++; //conteggio particelle cariche, va fatto solo per eventi 'realistici', in cui ha interagito nel target.
	   
	   for (int p = 0; p < 4; p++){ //reset delle posizioni
	     SpectrometerX[p] = -1000.; //valore di reset mai presente nelle distribuzioni
	     SpectrometerY[p] = -1000.;
	     SpectrometerZ[p] = -1000.; 
	  }
	   ///////////////////////////////////////////////PIXEL ANALYSIS PART/////////////////////////////7
	   //reset pixel counter
	   npixel = 0;
	     for (int istation = 0; istation < 6; istation++){
	     transversedstation[istation] = false;
	   }
	   //start loop
	   for (int l = 0; l < arr2->GetEntriesFast();l++){
	     SpectrometerPoint * spectro = (SpectrometerPoint*) arr2->At(l);
	     if ((spectro->GetTrackID() == j) &&  (spectro->GetDetectorID()>=111) && (spectro->GetDetectorID()<=162) ){ //stazioni PIXEL       
	       nfiglieinspettrometro[charmindex]++;
	       
               Int_t nstation = (Int_t)(spectro->GetDetectorID()- 111)/10; //numero della stazione (da 0 a 5)

	       transversedstation[nstation] = true;
	       
	       SpectrometerX[0] = spectro->GetX();
	       SpectrometerY[0] = spectro->GetY();
	       SpectrometerZ[0] = spectro->GetZ();
	       SpectrometerPX[0] = spectro->GetPx(); //uso l'angolo dell'ultimo spectrometerpoint visto, per ora
	       SpectrometerPY[0] = spectro->GetPy();
	       SpectrometerPZ[0] = spectro->GetPz();
	     }
	   }
	   for (int istation = 0; istation < 6; istation++){
	     if (transversedstation[istation]) npixel++;
	   }
	   hnpixel->Fill(npixel);
	   if (npixel==6){
	     rivelate++;
	     nfiglieinspettrometro[charmindex]++;
	   }           
	

	   //Muflux spectrometerPoint Drifttubes T3 and T4
	   for (int l = 0; l < arr3->GetEntriesFast();l++){
	     MufluxSpectrometerPoint * spectro = (MufluxSpectrometerPoint*) arr3->At(l);          
               if ((spectro->GetTrackID() == j) && (spectro->GetDetectorID()<32e+6)){ //terzo piano (con il Goliath) (riavvicinati o meno) 
	       rivelate3++;
	       SpectrometerX[2] = spectro->GetX();
	       SpectrometerY[2] = spectro->GetY();
	       SpectrometerZ[2] = spectro->GetZ();
	       SpectrometerPX[2] = spectro->GetPx();
	       SpectrometerPY[2] = spectro->GetPy();
	       SpectrometerPZ[2] = spectro->GetPz();
	       impulso3 = TMath::Sqrt(pow(SpectrometerPX[2],2) + pow(SpectrometerPY[2],2) + pow(SpectrometerPZ[2],2));

	     }
	     if ((spectro->GetTrackID() == j) && (spectro->GetDetectorID()>32e+6)){ //ultimo piano (con il Goliath)  (riavvicinati o meno)
	       //if (i % 100 == 0) cout<<"Uscita dallo spettrometro la traccia di ID: "<< spectro->GetTrackID()<< " di nome "<<spectro->PdgCode()<< " evento " << i <<endl;
	       uscite++;
	       SpectrometerX[3] = spectro->GetX();
	       SpectrometerY[3] = spectro->GetY();
	       SpectrometerZ[3] = spectro->GetZ();
	       SpectrometerPX[3] = spectro->GetPx();
	       SpectrometerPY[3] = spectro->GetPy();
               SpectrometerPZ[3] = spectro->GetPz();
	       impulso4 = TMath::Sqrt(pow(SpectrometerPX[3],2) + pow(SpectrometerPY[3],2) + pow(SpectrometerPZ[3],2));
	       nfiglieoutspettrometro[charmindex]++;
	       hp->Fill(impulso4);
	     }
          }
	 
	    
	  //calcolo le differenze fra gli angoli	  
	   if (npixel > 0 && (SpectrometerX[2] != -1000.) && (SpectrometerY[0] != -1000.) && (SpectrometerZ[2] != -1000.) && (SpectrometerX[3] != -1000.) && (SpectrometerY[3] != -1000.) && (SpectrometerZ[3] != -1000.)){
	     
	     for (int k=0; k< 4; k++){ //inserisco l'effetto di risoluzione
	       SpectrometerX[k] = SpectrometerX[k] + g->Gaus(0, resolution);
	       SpectrometerY[k] = SpectrometerY[k] + g->Gaus(0, resolution);
	     }
	     istoposizioni1->Fill(SpectrometerX[0], SpectrometerY[0]); //per lo studio della distribuzione spaziale
	     istoposizioni2->Fill(SpectrometerX[1], SpectrometerY[1]);
	     istoposizioni3->Fill(SpectrometerX[2], SpectrometerY[2]);
	     istoposizioni4->Fill(SpectrometerX[3], SpectrometerY[3]);

	     // hp->Fill(TMath::Sqrt(px**2 + py**2 + pz**2));	     

	     Theta1x = TMath::ATan(SpectrometerPX[1]/SpectrometerPZ[0]);
	     Theta2x = TMath::ATan((SpectrometerX[3]-SpectrometerX[2])/(SpectrometerZ[3]-SpectrometerZ[2]));
	     Deltathetax = Theta1x-Theta2x; //positiva per particelle cariche positivamente, negativa particelle cariche negativamente

	     //cout<<Deltathetax<<" "<<pdg->GetParticle(startpdg)->Charge()<<endl;
	     //cout<<SpectrometerX[0]<<" "<<SpectrometerY[0]<<" "<<SpectrometerX[1]<<" "<<SpectrometerY[1]<<" "<<SpectrometerX[2]<<" "<<SpectrometerY[2]<<" "<<SpectrometerX[3]<<" "<<SpectrometerY[3]<<endl;
	     if ((pdg->GetParticle(startpdg))->Charge() > 0) istodeltaxpiu ->Fill(Deltathetax);
	     if ((pdg->GetParticle(startpdg))->Charge() < 0) istodeltaxmeno ->Fill(Deltathetax);
	     istodeltax ->Fill(Deltathetax);
	     if (Deltathetax > (3 * sigmatheta)) ncaricapiu++;
	     if (Deltathetax < - (3 * sigmatheta)) ncaricameno++;
	     
	     if (((pdg->GetParticle(startpdg))->Charge() < 0) && (Deltathetax > 3*sigmatheta)) nmisid++; //misidentificazione
	     if (((pdg->GetParticle(startpdg))->Charge() > 0) && (Deltathetax <-(3*sigmatheta))) nmisid++;

	     if (TMath::Abs(Deltathetax) > (3 * sigmatheta)) cali->SetPoint(cali->GetN(), TMath::Abs(1./Deltathetax), momentum); //riempo la curva di calibrazione
	     
	     if (TMath::Abs(Deltathetax) > (3 * sigmatheta)) istodp->Fill((f2mono->Eval(TMath::Abs(1./Deltathetax))-momentum)/momentum);
	     //cout<<f2->Eval(TMath::Abs(1./Deltathetax))<<" "<<TMath::Abs(1./Deltathetax)<<endl;

	     
	     Theta1y = TMath::ATan(SpectrometerPY[0]/SpectrometerPZ[0]);
	     Theta2y = TMath::ATan((SpectrometerY[3]-SpectrometerY[2])/(SpectrometerZ[3]-SpectrometerZ[2]));
	     //Theta1y = TMath::ATan((SpectrometerY[1]-SpectrometerY[0])/21.25);
	     //Theta1y = TMath::ATan((SpectrometerY[3]-SpectrometerY[2])/22.);
	     Deltathetay = Theta1y-Theta2y;
	     // cout<<SpectrometerZ[1]-SpectrometerZ[0]<<" "<<SpectrometerZ[3]-SpectrometerZ[2]<<endl;
	     // if (TMath::Abs(Deltathetay) > 0.003) cout<<Deltathetay<<" "<<TMath::Sqrt(px**2 + py**2 + pz**2)<<endl;
	     istodeltay ->Fill(Deltathetay);
	     
	    
	   } //condizione tutti hanno parlato
	   
	   zmax = -100;
	   nemulsions = 0;
	   //parte sul target:
	   /* for (int t = 0; t <arr3->GetEntriesFast();t++){
	      BoxPoint *target = (BoxPoint*) arr3->At(t);
	      if (target->GetTrackID() == j){
	      if (target->GetZ() > zmax) zmax = target->GetZ();
	      nemulsions++;
	      }
	      }*/
	   hz -> Fill(zmax);
	     if (zmax < -2.35) hz1->Fill(zmax);
	     // cout<<" id traccia figlia di charm: "<< j << " z max: "<< zmax<<" numero evento: "<<i<<endl;
	} //condizione figlia di charm
	
	
	}//ciclo sui due charm
      
    } //condizione mumID > 0
   }//ciclo sulle tracce
   
   isto->Fill(numerofigliecariche[0]);
   isto->Fill(numerofigliecariche[1]);
   //cout<<numerofigliecariche[0]<<" "<<numerofigliecariche[1]<<" "<<nfiglieinultimaemulsione[0]<<" "<<nfiglieinultimaemulsione[1]<<endl;

   if ((nfiglieinspettrometro[0] >= 1) && (nfiglieinspettrometro[1] >=1) ) tutteinspettrometrodouble++;
   if ((nfiglieoutspettrometro[0] >= 1) && (nfiglieoutspettrometro[1] >=1) ){
     tutteoutspettrometrodouble++;
     if (isprimary) primaryok++;
   }
   if (nfiglieinspettrometro[0] >= 1) tutteinspettrometrosingle++;
   if (nfiglieoutspettrometro[0] >= 1) tutteoutspettrometrosingle++;

   if (nfiglieinspettrometro[1] >= 1) tutteinspettrometrosingle++;
   if (nfiglieoutspettrometro[1] >= 1) tutteoutspettrometrosingle++;
       
   //if ((numerofigliecariche[0] == nfiglieoutspettrometro[0]) && (numerofigliecariche[0] != 0)) tutteoutspettrometrosingle++;
   // if ((numerofigliecariche[1] == nfiglieoutspettrometro[1]) && (numerofigliecariche[1] != 0)) tutteoutspettrometrosingle++;



   //if (i%100 == 0) cout<<endl;
 }//ciclo sugli eventi
 cout<<ncut<<endl;
 cout<<"Vertici primari fuori dall'ECC: "<<nulli<<endl;
 //cout<<tutteinultimaemulsione<<endl;
 cout<<tutteinspettrometrosingle<<" "<<tutteoutspettrometrosingle<<endl;
 cout<<tutteinspettrometrodouble<<" "<<tutteoutspettrometrodouble<<endl;
 cout<<primaryok<<endl;
 cout<<"Numero di eventi: "<< nevents<<endl;
 cout<<"Numero totale di figlie di charm cariche (prodotte nel target" <<ntotalefigliecariche<<endl;
 cout<<"Rivelate:" <<rivelate<<endl;
 cout<<"Frazione di figlie di charm rivelate all'inizio dello spettrometro: (ha senso se nulli = 0) "<< (Double_t)rivelate/ntotalefigliecariche<<endl;
 //cout<<"Frazione di figlie di charm rivelate nel secondo piano dello spettrometro: "<< (Double_t)rivelate2/ntotalefigliecariche<<endl;
 // cout<<"Frazione di figlie di charm rivelate nel terzo piano dello spettrometro: "<< (Double_t)rivelate3/ntotalefigliecariche<<endl;
 cout<<"Frazione di figlie di charm rivelate alla fine dello spettrometro: (ha senso se nulli = 0) "<< (Double_t)uscite/ntotalefigliecariche<<endl;
 cout<<"Frazione di cariche misidentificate"<<(Double_t)nmisid/ntotalefigliecariche<<endl;

 cout<<"Efficienze di single: in = "<<(Double_t) tutteinspettrometrosingle/(2 * nevents - 2 * nulli)<<" out = "<<(Double_t) tutteoutspettrometrosingle/(2 * nevents - 2 * nulli)<<endl;
 cout<<"Efficienze di double: in = "<<(Double_t) tutteinspettrometrodouble/(nevents - nulli)<<" out = "<<(Double_t) tutteoutspettrometrodouble/(nevents - nulli)<<endl;
 // gStyle->SetOptStat("emr");
 TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
 
 TCanvas *c = new TCanvas();
 isto->Draw();
 isto->SetTitle("");
 isto->GetXaxis()->SetTitle("N daughters");
 TCanvas *c1 = new TCanvas();
 
 istodeltaxpiu->SetLineColor(kOrange+10);
 istodeltaxmeno->SetLineColor(kAzure+6);
 istodeltaxpiu->SetTitle("");
 istodeltaxpiu->GetXaxis()->SetTitle("#Delta #Theta_{x}");
 istodeltaxpiu->Draw();
 istodeltaxmeno->Draw("SAMES");

 leg->AddEntry("istodeltaxpiu", "Positive particles");
 leg->AddEntry("istodeltaxmeno", "Negative particles");

 leg->Draw();
 c1->Print("./temporanei/analisispettrometro/angle_xcomponent.root","root");
 TCanvas *cy = new TCanvas();
 istodeltay->Draw();
 istodeltay->SetTitle("");
 istodeltay->GetXaxis()->SetTitle("#Delta #Theta_{y}");
 istodeltay->Fit("f1","R");
 cout<<"RMS deltathetay: "<<f1->GetParameter(2)<<endl;
 cout<<"Efficienza di carica: "<<(ncaricapiu + ncaricameno)/(istodeltaxpiu->GetEntries() + istodeltaxmeno->GetEntries())<<endl;

 cy->Print("./temporanei/analisispettrometro/angle_ycomponent.root","root");
 TCanvas *c2 = new TCanvas();
 cali->GetXaxis()->SetTitle("1/DeltaThetax");
 cali->GetYaxis()->SetTitle("P");
 cali->Draw("A*");
 cali->Fit("f2mono","R");
 gStyle->SetOptFit(111111);
 c2->Print("./temporanei/analisispettrometro/calibration.root","root");
 TCanvas *cdp = new TCanvas();
 istodp->Fit("f3", "R");
 istodp->GetXaxis()->SetTitle("#Delta P/P");
 istodp->Draw();
 cdp->Print("./temporanei/analisispettrometro/momentum_resolution.root","root");

 hefflocal->Divide(hp,hpin);
 TCanvas *ceff = new TCanvas();
 hefflocal->Draw();
 TCanvas *cp = new TCanvas();
 hp->GetXaxis()->SetTitle("Gev/c");
 hp->Draw();
 cp->Print("./temporanei/analisispettrometro/momentum_daughters.root","root");
 TCanvas *cpwrong = new TCanvas();
 hpin->Draw();
 TCanvas *cangle = new TCanvas();
 hangle->Draw();
 TCanvas *canglewrong = new TCanvas();
 hanglewrong->Draw();
 TCanvas *c3 = new TCanvas();
 c3->Divide(2,2);

 c3->cd(1);
 istoposizioni1->SetTitle("First Plane");
 istoposizioni1->GetXaxis()->SetTitle("x[cm]");
 istoposizioni1->GetYaxis()->SetTitle("y[cm]");
 istoposizioni1->Draw();
 
 c3->cd(2);
 istoposizioni2->SetTitle("Second Plane");
 istoposizioni2->GetXaxis()->SetTitle("x[cm]");
 istoposizioni2->GetYaxis()->SetTitle("y[cm]");
 istoposizioni2->Draw();

 c3->cd(3);
 istoposizioni3->SetTitle("Third Plane");
 istoposizioni3->GetXaxis()->SetTitle("x[cm]");
 istoposizioni3->GetYaxis()->SetTitle("y[cm]");
 istoposizioni3->Draw();

 c3->cd(4);
 istoposizioni4->SetTitle("Fourth Plane");
 istoposizioni4->GetXaxis()->SetTitle("x[cm]");
 istoposizioni4->GetYaxis()->SetTitle("y[cm]");
 istoposizioni4->Draw();

 c3->Print("./temporanei/analisispettrometro/daughters_in_spectrometer_planes.root","root");
 

 TCanvas *c4 = new TCanvas();
 hz->Draw();
 TCanvas *c5 = new TCanvas();
 hz1->Draw();
 TCanvas *c6 = new TCanvas();
 hnpixel->Draw();
 
 return 0;
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
