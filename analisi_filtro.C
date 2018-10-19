//studio sulle prestazioni del filtro per muoni sui prodotti di decadimento dei charm (creato il 31 Maggio)
bool ischarm(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);
#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
#include "/home/utente/Scrivania/SHIPBuild/FairShip/charmdet/MuonTaggerPoint.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-13-02-ship-1/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

void analisi_filtro(TString filename = ""){
 TFile *file1;
 if (filename != "") file1 = TFile::Open(filename);
 else {
   // file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_centralbeam_10000.root");
// file1 = TFile::Open("/home/utente/SHIPBuild/sim_charm/OPERAECC/muontaggeracceptance/ship.conical.Pythia8CharmOnly-TGeant4.root");
 file1 = TFile::Open("/home/utente/Simulations/sim_charmdet/charm_events/ECC3_acceptance_morefardetectors/ship.conical.Pythia8CharmOnly-TGeant4.root");
 }
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");
 //ottenere posizione e dimensione bersaglio di neutrini tau
 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import("/home/utente/Simulations/sim_charmdet/charm_events/ECC3_acceptance_morefardetectors/geofile_full.conical.Pythia8CharmOnly-TGeant4.root"); 

 TGeoVolume *top = gGeoManager->GetTopVolume(); 
 TGeoNode *targetnode = top->GetNode("volTarget_1");
 if (targetnode == NULL){ cout<<"Error: Target volume not found"<<endl;
   return -1;
 }
 TGeoMatrix *transmatrix = targetnode->GetMatrix();
 Double_t zpostarget = transmatrix->GetTranslation()[2];
 Double_t zdimtarget = (((TGeoBBox*)targetnode->GetVolume()->GetShape())->GetDZ())*2.;
 cout<<"z position of center of target: "<<zpostarget<<"dim z: "<<zdimtarget<<endl;

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
 
 Double_t checkmin;

 const Double_t check = zendactive;
 const Double_t check1 = check;

 TH1D *hl = new TH1D("hl", "Lunghezza volo charm", 100,0,30);
 TH1D * isto = new TH1D("isto", "Numero figlie di charm cariche",5,0,5);
 // TH1D * istofiltromuoni = new TH1D("histomuons", "Numero piani attraversati", 7,0,7);
 // TH1D * istofiltroaltri = new TH1D("histohadrons", "Numero piani attraversati", 7,0,7);
 TH1D *hpmuons = new TH1D("hpmuons", "Muon momentum", 1000,0,100);
 TH1D * istofiltromuoni = new TH1D("histomuons", "Numero piani attraversati", 6,1,7); //il filtro inizia con un piano sensibile, se entrano devono aver attraversato almeno il primo
 TH1D * istofiltroaltri = new TH1D("histohadrons", "Numero piani attraversati",6,1,7);	

 TH2D *hxy = new TH2D("hxy", "posxy muoni da charm", 1000, -100, 100, 1000, -50, 50);
	    
 FairMCEventHeader *prova = NULL;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.", &prova);
 
 TClonesArray *arr0 = new TClonesArray("ShipMCTrack",10000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr0);

 TClonesArray *arr1 = new TClonesArray("SpectrometerPoint",10000);
 cbmsim->GetBranch("SpectrometerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("SpectrometerPoint",&arr1);
 
 
 TClonesArray *arr2 = new TClonesArray("MuonTaggerPoint",10000);
 cbmsim->GetBranch("MuonTaggerPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MuonTaggerPoint", &arr2);
 cbmsim->SetBranchStatus("*",1);

 Double_t Vz = 0;
 Double_t Vy = 0;
 Double_t Vx = 0;
 Double_t startz = 0;
 Double_t starty = 0;
 Double_t startx = 0;
 Double_t lunghezzavolo[2];
 bool checklunghezzavolo[2];
 
 bool off[2];
 bool selezione;

const Int_t nevents = cbmsim->GetEntries();
 Int_t ncut = 0; //numero di eventi tagliati perchè il vertice secondario è fuori dalla zona di check
 Int_t ntotalefigliecariche = 0; //numero totale di figlie cariche (non nel singolo evento)
 Int_t rivelate = 0;
 Int_t rivelate2 = 0;
 Int_t uscite = 0;
 Int_t startpdg = 0;
 Int_t intermediatepdg = 0;
 Int_t mumID = 0;
 Int_t mumpdg = 0;
 Int_t hitID = 0;
 Int_t success = 0; //numero figlie di charm in cui l'identificazione è andata a buon fine.
 Int_t numerofigliecariche[2];
 Int_t numerofiglierivelate[2];
 Int_t singleinfiltro = 0;
 Int_t doubleinfiltro = 0;
 Int_t pdgcharm[2];
 Int_t ncaricapiu = 0; //numero di figlie di charm positive in cui deltathetax supera la soglia di 3 sigma per la misura di carica
 Int_t ncaricameno = 0; //numero di figlie di charm negative come sopra

 Int_t figliecharm[2];
 Int_t nulli = 0;

 Int_t ntot = 0;
 Int_t nmuonsfromcharm = 0; //numero totale muoni figli di charm
 for (int i = 0; i < nevents; i++){
   arr0->Clear();
   arr1->Clear();
   arr2->Clear();
   cbmsim->GetEntry(i);
   selezione = true;
   Vz = prova->GetZ();
   Vy = prova->GetY();
   Vx = prova->GetX();

   if ((Vz > check) || (Vz < checkmin)){
     nulli++;
     //continue;
   }

   figliecharm[0] = 0;
   figliecharm[1] = 0;

   lunghezzavolo[0] = 0.;
   lunghezzavolo[1] = 0.;
   checklunghezzavolo[0] = false; //Fare il fill una sola volta per charm
   checklunghezzavolo[1] = false;
   
   numerofigliecariche[0] = 0;
   numerofigliecariche[1] = 0;
   numerofiglierivelate[0] = 0;
   numerofiglierivelate[1] = 0;
   
   pdgcharm[0] = 0;
   pdgcharm[1] = 0;
   off[0] = false;
   off[1] = false;

   figliecharm[0] = 0;
   figliecharm[1] = 0;
  
   for (int j = 1; j < arr0->GetEntriesFast();j++){

     ShipMCTrack *trk = (ShipMCTrack*) arr0->At(j);
     Double_t momentum = trk->GetP();
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
    if (mumID>0){
       trk = (ShipMCTrack*)arr0->At(mumID);
       mumID = trk->GetMotherId();
       mumpdg = trk->GetPdgCode();
      while ((isintermediate(mumpdg) == true) && (mumID>0)) //se la madre è uno stato intermedio controllo se lo stato intermedio è stato prodotto dal charm
	 {
	   intermediatepdg = mumpdg;
	   
	   trk = (ShipMCTrack*) arr0->At(mumID);
	   mumID = trk->GetMotherId();
	   mumpdg = trk->GetPdgCode();
	 }
      //salto il passaggio in cui lo stato intermedio è ricondotto alla madre, avendo contato i suoi prodotti
      for (int charmindex = 0; charmindex < 2; charmindex++){
      //      for (int charmindex = 0; charmindex < 1; charmindex++){//caso di un solo charm per evento
	if (mumpdg==pdgcharm[charmindex]){
	   
	  figliecharm[charmindex]++;
	  lunghezzavolo[charmindex] = TMath::Sqrt(pow((startz-Vz),2) + pow((startx-Vx),2) + pow((starty-Vy),2));
	  
	  if ((lunghezzavolo[charmindex] != 0.) &&(checklunghezzavolo[charmindex]==false)){
	    hl->Fill(lunghezzavolo[charmindex]*10);
	    checklunghezzavolo[charmindex] = true;
	  }
          
	}
	if ((mumpdg==pdgcharm[charmindex]) && (isintermediate(startpdg) == false) && (abs((pdg->GetParticle(startpdg))->Charge()) > 0.) ){

	  if (i % 1000 == 0) cout << "ID Traccia figlia di charm: " << j << "di nome: "<< pdg->GetParticle(startpdg)->GetName()<< "pdg" << startpdg<< " "<< i<<endl;

          
	 /* if (startz > check){ //controllo la posizione del vertice secondario, è all'interno della regione di controllo?
	    if (off[charmindex] == false) ncut++;
	   
	    off[charmindex] = true;
	   
       }
       if (off[charmindex] == true) continue;*/
	  int ntagged = 0; //per ogni figlia di charm conto quanti piani sono stati attraversati
/*	  for (int l = 0; l < arr1->GetEntriesFast(); l++){ //questa parte non capisco che senso abbia ora
	   SpectrometerPoint* spectro = (SpectrometerPoint*) arr1->At(l);
           if ((abs(spectro->GetY()) < 50)){
	   if ((spectro->GetTrackID() == j) && (spectro->GetZ() > 500)) ntagged++;
 }
}*/
	  for (int l = 0; l < arr2->GetEntriesFast();l++){
	    MuonTaggerPoint* filter = (MuonTaggerPoint*) arr2->At(l);
	    if (filter->GetTrackID() == j){ 
             ntagged++;	    
             // if (abs(startpdg) == 13) cout<<"evento con muone "<<i<<endl;
	     //if ((abs(startpdg) == 13) && (filter->GetZ() > 695)) hxy->Fill(filter->GetX(), filter->GetY());
	     //if ((abs(startpdg) == 13) && (filter->GetZ() < 520)) hxy->Fill(filter->GetX(), filter->GetY());
             }
	    if (filter->GetDetectorID() == 1) numerofiglierivelate[charmindex]++;            
  }
	  // cout<<pdg->GetParticle(startpdg)->GetName()<<" " << startpdg<<endl;
	  //attenzione! Il primo piano è quello sensibile, se non è stato visto neanche quello la particella non va contata!
	  if ((abs(startpdg) == 13) && (ntagged >= 1)){
            hpmuons->Fill(momentum);
	    istofiltromuoni->Fill(ntagged); //muoni
	  }
	  else if (((abs(startpdg) == 211) || (abs(startpdg) == 321)) && (ntagged >= 0)){
	    istofiltroaltri->Fill(ntagged);    //pioni o kaoni
	  }
	  if (abs(startpdg) == 13){
	      if (ntagged >= 3) success++;
	    ntot++;
	    nmuonsfromcharm++;
	    }
	  else if ((abs(startpdg) == 211) || (abs(startpdg) == 321)){
              if (ntagged < 3) success++;
	      ntot++;
	      }
	  // else if (ntagged > 0) cout<< startpdg<<" "<<ntagged<<endl;
	  numerofigliecariche[charmindex]++;
	  ntotalefigliecariche++;
	  //calcolo le differenze fra gli angoli
	}
	  
	
      }
      
    }
   }

     if (off[0] == false) isto->Fill(numerofigliecariche[0]);
     if (i % 1000 == 0) cout<<numerofigliecariche[1]<< "pdg charm 1 "<< pdgcharm[0] <<" pdg charm 2 "<< pdgcharm[1]<<endl;
     if (off[1] == false) isto->Fill(numerofigliecariche[1]);
  
     if (lunghezzavolo[0] == 0) cout <<"Ahio: non c'è lunghezza di volo per il primo charm "<<pdgcharm[0]<<" nell'evento " <<i<<endl;
     // if (lunghezzavolo[1] == 0) cout <<"Ahio: non c'è lunghezza di volo per il secondo charm "<<pdgcharm[1]<<" nell'evento " <<i<<endl;
     if (i % 1000 == 0) cout<<i<<endl;
 }
 cout<<"Numero di eventi: "<< nevents<<endl;
 cout<<"Numero totale di figlie di charm" <<ntotalefigliecariche<<endl;
 //cout<<"Single in filtro "<<singleinfiltro<<" Double in filtro: "<<doubleinfiltro<<endl;
 cout<<"Percentuale di particelle identificate correttamente"<<(double) success/ntot<<endl;
 gStyle->SetOptStat("emr");
 TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
 
 TCanvas *c = new TCanvas();
 isto->Draw();
 TCanvas *cl = new TCanvas();
 hl->Draw();
 cl->Print("./temporanei/analisifiltro/decay_length.root","root");
 hl->Fit("expo");
 TCanvas *c1 = new TCanvas();
 istofiltroaltri->SetLineColor(kBlue);
 istofiltromuoni->SetLineColor(kRed); 
 cout<<"Numero totale muoni figli di charm: "<<nmuonsfromcharm<<endl;
 cout<<"Numero altre particelle cariche: "<< istofiltroaltri->Integral()<<endl;
 istofiltromuoni->Scale(1/istofiltromuoni->Integral());
 istofiltroaltri->Scale(1/istofiltroaltri->Integral());
 istofiltroaltri->Draw("hist");
 istofiltromuoni->Draw("hist");
 istofiltromuoni->SetTitle("");
 istofiltromuoni->GetXaxis()->SetTitle("Nplanes");
 leg->AddEntry("histomuons", "Muons");
 leg->AddEntry("histohadrons", "Hadrons");

 leg->Draw();
 c1->Print("./temporanei/analisifiltro/muon_filter.root","root");
 
 cout<<endl;
 cout<<"Richiedendo almeno 2 piani:"<<endl;
 cout<<"1-Alpha = "<<istofiltromuoni->Integral(2,7)<<endl;
 cout<<"Beta = "<< istofiltroaltri->Integral(2,7)<<endl;
 cout<<endl;
 cout<<"Richiedendo almeno 3 piani:"<<endl;
 cout<<"1-Alpha = "<<istofiltromuoni->Integral(3,7)<<endl;
 cout<<"Beta = "<< istofiltroaltri->Integral(3,7)<<endl;
 cout<<endl;
 cout<<"Richiedendo almeno 4 piani:"<<endl;
 cout<<"1-Alpha = "<<istofiltromuoni->Integral(4,7)<<endl;
 cout<<"Beta = "<< istofiltroaltri->Integral(4,7)<<endl;
 cout<<"Richiedendo esattamente 6 piani:"<<endl;
 cout<<"1-Alpha = "<<istofiltromuoni->Integral(6,7)<<endl;
 cout<<"Beta = "<<istofiltroaltri->Integral(6,7)<<endl;
 
 cout<<istofiltromuoni->Integral()<< " "<< istofiltroaltri->Integral()<<endl;

 TCanvas *c2 = new TCanvas();
 hxy->Draw();

 TCanvas *cp = new TCanvas();
 hpmuons->Draw();
}
bool ischarm(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 431) || (abs(PdgCode) == 411) || (abs(PdgCode) == 4122)  || (abs(PdgCode) == 421) || (abs(PdgCode) == 4132) || (abs(PdgCode) == 4232) ||(abs(PdgCode) == 4332)|| (PdgCode == 441))   check = true; 
      return check;
      }
//lista mesoni charmati: D0(bar), D+(-), Ds0(bar), Ds+(-), Lambdac+(-), Csi+(-), Csi0(bar), Omega_co(bar), eta_c

bool isintermediate(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 223) || (abs(PdgCode) == 3224) || (abs(PdgCode) == 331) || (abs(PdgCode) == 221) || (abs(PdgCode) == 20213) || (abs(PdgCode) == 3212) ||(abs(PdgCode)==213) || (abs(PdgCode) == 113) || (abs(PdgCode) == 2224) || (abs(PdgCode) == 323)) check =  true;
  return check;
  }

//lista intermedi: omega, sigma*, eta', eta, a1,pho+-,pho0, K*+-, delta++, sigma0
