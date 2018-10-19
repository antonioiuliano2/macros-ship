//studio sulla simulazione di produzione di charm in cascata.
bool ischarm(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);
double whatmotherofcharm(Int_t PdgCode);

#include "/home/utente/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/home/utente/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
#include "/home/utente/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

void distribuzioni_charm(TString filename = "" ){ 
 TFile *file1;
 if (filename != "") file1 = TFile::Open(filename);
 else{
  //TFile *file1 = TFile::Open("./10run/ship.10.0.Pythia8CharmOnly-TGeant4_primo.root");
 file1 = TFile::Open("ship.10.0.Pythia8CharmOnly-TGeant4_centralbeam_10000.root");  
 }
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");

 const Double32_t MolTh = 0.3; //parametri da cambiare a seconda della lunghezza del bersaglio e del campionamento effettuato
 const Double32_t fFDs = 7.7/10.4;

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


 // const Double32_t check = Length - 3 *(MolTh + 2*0.0045 + 0.0205); 
 const Double32_t checkmax = 0;       //posizione in cui finisce il target.

 const int nrun = 10;
 Double32_t checkmin;
 switch(nrun){
 case 1: checkmin = -1000; //nel run 1 non ci sono check
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
 //const Double32_t checkmin = -6.0; //ultimi 60 cm
 cout<<checkmin<<endl;
 FairMCEventHeader *prova = NULL;
 cbmsim->GetBranch("MCEventHeader.");
 cbmsim->SetBranchAddress("MCEventHeader.", &prova);
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);

 TH1D *hpz = new TH1D("hpz", "componente z dell'impulso del protone madre di charm", 50, 0, 400);
 TH1D *hp = new TH1D("hp", "impulso protone madre di charm",100, 0, 400);
 TH1D *hpc1 = new TH1D("hpprimary", "impulso charm primari", 100, 0, 350);
 TH1D *hpc2 = new TH1D("hpsecondary", "impulso charm secondari", 100, 0, 350);
 TH1D *hpf1 = new TH1D("hpf1", "impulso figlie di charm primari",100, 0, 40);
 TH1D *hpf2 = new TH1D("hpf2", "impulso figlie di charm secondari",100, 0, 40);
 TH1D *hkinkangle1 = new TH1D("hkinkangle1", "angolo di kink 1 prong",100, 0, 0.6);
 TH1D *hkinkangle3 = new TH1D("hkinkangle3", "angolo di kink 3 prong",100, 0, 0.6);
 
 TH1D *hpf1_T = new TH1D("hpf1_T", "impulso trasverso figlie di charm 1 prong",100, 0, 2);
 TH1D *hpf3_T = new TH1D("hpf3_T", "impulso trasverso figlie di charm 3 prong",100, 0, 2);
 TH1D *hl = new TH1D("hl", "lunghezza volo charm",5000, 0, 5.0);
 TH1D *hvz = new TH1D("hvzprimary", "posizione vertice di produzione charm primari", 100, -50, 0);
 TH1D *hvz1 = new TH1D("hvzsecondary", "posizione vertice di produzione charm secondari", 100, -50, 0);
 TH1D *hwho = new TH1D ("hwho", "Differenti specie madri di charm secondari",10, 0, 10);	
 TH1D *hax1 = new TH1D("hax1", "Angolo x charm primari",100,-0.3,0.3);
 TH1D *hax2 = new TH1D("hax2", "Angolo x charm secondari",100,-0.3,0.3);
 TH1D *hay1 = new TH1D("hay1", "Angolo y charm primari",100,-0.3,0.3);
 TH1D *hay2 = new TH1D("hay2", "Angolo y charm secondari",100,-0.3,0.3);
 
 TF1* f = new TF1("f", "expo", 0,120);
 //TF1* f = new TF1("f", "expo", -6, 25);
 Double32_t length;
 Double32_t Vz = 0;
 Double32_t Vy = 0;
 Double32_t Vx = 0;
 Double32_t startz = 0;
 Double32_t starty = 0;
 Double32_t startx = 0;
 Double32_t endz = 0;
 Double32_t endy = 0;
 Double32_t endx = 0;
 Double32_t pz, impulso;
 Double32_t impulsofiglia;
 bool selezione;

 Int_t primari = 0;
 Int_t secondari = 0;

 Double32_t startpx; //parametri traccia iniziale
 Double32_t startpy;
 Double32_t startpz;

const Int_t nevents = cbmsim->GetEntries();
 Int_t ncut = 0; //numero di eventi tagliati perchè il vertice secondario è fuori dalla zona di check
 
 Int_t startpdg = 0;
 Int_t intermediatepdg = 0;
 Int_t mumID = 0;
 Int_t mumpdg = 0;
 Int_t charmpdg = 0;
 Int_t hitID = 0;
 
 Int_t pdgparent;
 bool parentfound;

 
 Double_t angolofiglia;
 
 Int_t pdgcharm[2];
 Double_t charmpx[2];
 Double_t charmpy[2];
 Double_t charmpz[2];
 Int_t charmlength[2];
 Double_t impulsocharm[2];
 Double_t kink[2];

 Int_t nfiglia[2]; //numero figlie dei 2 charm
 Double_t pxtotfiglia[2];
 Double_t pytotfiglia[2];
 Double_t pztotfiglia[2]; 
 Double_t impulsototfiglia[2];
 Double_t impulsotrasversofiglia[2];
 
 cout<<"Numero eventi: "<<nevents<<endl;
 for (int i = 0; i < nevents; i++){
  
  parentfound = false;

  pdgcharm[0] = 0;
  pdgcharm[1] = 0;

  nfiglia[0] = 0;
  nfiglia[1] = 0;
  
  pxtotfiglia[0] = 0.;
  pytotfiglia[0] = 0.;
  pztotfiglia[0] = 0.;

  pxtotfiglia[1] = 0.;
  pytotfiglia[1] = 0.;
  pztotfiglia[1] = 0.;

  charmlength[0] = 0;
  charmlength[1] = 0;
  
  if (i % 100 == 0) cout<<i<<endl;
   arr1->Clear();
   cbmsim->GetEntry(i);

   selezione = true;
   /*Vz = prova->GetZ();
   Vy = prova->GetY();
   Vx = prova->GetX();*/
   
   ShipMCTrack *trk = (ShipMCTrack*) arr1->At(0);
   Vx = trk->GetStartX();
   Vy = trk->GetStartY();
   Vz = trk->GetStartZ();
   cout<<"PROVA:" <<trk->GetPdgCode()<<" "<<trk->GetPz()<<endl;
  
   // if (Vz > checkmax) continue;
   //if (Vz < checkmin) continue;
   
   for (int j = 1; j < arr1->GetEntriesFast();j++){

     ShipMCTrack *trk = (ShipMCTrack*) arr1->At(j);
     
      if (ischarm(trk->GetPdgCode()) && (pdgcharm[0]==0)){  //memorizzo gli indici dei due charm
       pdgcharm[0] = trk->GetPdgCode();
       charmpx[0] = trk->GetPx();
       charmpy[0] = trk->GetPy();
       charmpz[0] = trk->GetPz();
       impulsocharm[0] = pow((charmpx[0],2) + pow(charmpy[0],2) + pow(charmpz[0],2),0.5);
     }
     if (ischarm(trk->GetPdgCode()) && (trk->GetPdgCode() != pdgcharm[0]) && (pdgcharm[1]==0)){ //il secondo lo memorizza solo se diverso dal primo
       pdgcharm[1] = trk->GetPdgCode();
       charmpx[1] = trk->GetPx();
       charmpy[1] = trk->GetPy();
       charmpz[1] = trk->GetPz();
       impulsocharm[1] = pow((charmpx[1],2) + pow(charmpy[1],2) + pow(charmpz[1],2),0.5);
     }
       mumID = trk->GetMotherId(); //prendo l'id della madre di ogni traccia (diretta, non madre della catena)
       mumpdg = trk->GetPdgCode();
       startpdg = trk->GetPdgCode();
       
       if (startpdg > 10000) continue;
       
       endx = trk->GetStartX(); //coordinate di partenza della figlia
       endy = trk->GetStartY();
       endz = trk->GetStartZ(); 
       
       startpx = trk->GetPx(); //componenti impulso figlia
       startpy = trk->GetPy();
       startpz = trk->GetPz();
    
       impulsofiglia = pow(pow(trk->GetPz(),2) + pow(trk->GetPx(),2) + pow(trk->GetPy(),2),0.5);
       angolofiglia = TMath::ATan(TMath::Sqrt(pow((startpx/startpz),2) + pow((startpy/startpz),2)));
       //if (i % 100 == 0) cout<<"Il charm è : "<<pdg->GetParticle(startpdg)->GetName()<< " "<<endl;
      
       if (mumID>=0){
	 trk = (ShipMCTrack*)arr1->At(mumID);
	 mumID = trk->GetMotherId();
	 mumpdg = trk->GetPdgCode();
	 while ((isintermediate(mumpdg) == true) && (mumID>0)) //se la madre è uno stato intermedio controllo se lo stato intermedio è stato prodotto dal charm
	   {
	     intermediatepdg = mumpdg;
	     /*if (intermediatepdg == 2224){
	       cout<<"La traccia figlia: "<<startpdg<<" inizia in "<<endx<<" "<<endy<<" "<<endz<<endl;
	       cout<<"Stato intermedio: "<<trk->GetPdgCode()<<" prodotta in "<<trk->GetStartX()<<" "<<trk->GetStartY()<<" "<<trk->GetStartZ()<<endl;
	     }*/
	     trk = (ShipMCTrack*) arr1->At(mumID);
	     mumID = trk->GetMotherId();
	     mumpdg = trk->GetPdgCode();
	   }
	 pz = trk->GetPz();
	 impulso = TMath::Sqrt(pow(trk->GetPz(),2) + pow(trk->GetPx(),2) + pow(trk->GetPy(),2));
	 //if (i % 100 == 0) cout<<"Madre di charm: "<<mumpdg<<" ha impulso P: "<<impulso<< " componente z: "<<pz<<" Vz = "<<Vz<<endl;
	 /*if (fabs(startpdg) == 431){
       Double_t rnr = gRandom->Uniform(0,1);
       if( rnr>fFDs ) cout<<" Ds da scartare?"<<mumpdg<<" "<<pz<<endl; 
       }*/
	  if((startpdg == pdgcharm[0]) || (startpdg == pdgcharm[1])){
	    //if ((pz > 390) && (mumpdg == 2212)){
            if ((pz == 400) && (mumpdg == 2212)){ //provo a restringere il campo dei primari
	      primari++;
	      hvz->Fill(Vz);
	      hpc1->Fill(impulsofiglia);
	      hpz->Fill(pz);
              hax1->Fill(TMath::ATan(startpx/startpz));
              hay1->Fill(TMath::ATan(startpy/startpz));
	    }
	    else{
	      secondari++;
	      hvz1->Fill(Vz);
	      hpc2->Fill(impulsofiglia);
	       hpz->Fill(pz);
	       hp->Fill(impulso);
	      //cout<<"Madre di charm secondario: "<<mumpdg<<" ha impulso P: "<<impulso<< " componente z: "<<pz<<" Vz = "<<Vz<<endl;
	      pdgparent = mumpdg;	
	      if (parentfound == false) {
		hwho->Fill(whatmotherofcharm(pdgparent));
		parentfound = true;
	      }
              hax2->Fill(TMath::ATan(startpx/startpz));
              hay2->Fill(TMath::ATan(startpy/startpz));
	    }
	  }

	  if((mumpdg == pdgcharm[0]) || (mumpdg == pdgcharm[1])){ //il charm qui è la madre di un'altra particella
	    startx = trk->GetStartX();
	    starty = trk->GetStartY();
	    startz = trk->GetStartZ();
	    length = TMath::Sqrt(pow((endz - startz),2) + pow((endy - starty),2) + pow((endx -startx),2));

	    mumID = trk->GetMotherId();
            charmpdg = trk->GetPdgCode();
	    
	    if(mumpdg == pdgcharm[0] && (abs(pdg->GetParticle(startpdg)->Charge()) > 0)){ //figlia del primo charm
		 nfiglia[0]++;
		 pxtotfiglia[0] += startpx;
		 pytotfiglia[0] += startpy;
		 pztotfiglia[0] += startpz;		 		 
	    }
	    if((mumpdg == pdgcharm[1]) && (abs(pdg->GetParticle(startpdg)->Charge()) > 0)){ //figlia del secondo charm
		 nfiglia[1]++;
		 pxtotfiglia[1] += startpx;
		 pytotfiglia[1] += startpy;
		 pztotfiglia[1] += startpz;
	    }
	    trk = (ShipMCTrack*)arr1->At(mumID);
	    mumpdg = trk->GetPdgCode();
	    //cout<<"Madre di charm: "<<mumpdg<<endl;
	    pz = trk->GetPz();
            if ((charmpdg != charmlength[0]) && (charmpdg != charmlength[1])){ //ricopiato da una vecchia versione
		//	cout<<i<<"PROVA"<<endl;

		hl->Fill(length);
	
		if (charmlength[0] == 0) charmlength[0] = charmpdg;
		else if (charmlength[1] == 0) charmlength[1] = charmpdg;
	      }

             if(abs(pdg->GetParticle(startpdg)->Charge()) > 0){ //solo figlie di charm cariche
	       if ((pz > 390) && (mumpdg == 2212)) hpf1->Fill(impulsofiglia);
               else hpf2->Fill(impulsofiglia);
	       // impulsotrasversofiglia1 = impulsofiglia * TMath::Sin(angolofiglia);
	       // hpf_T->Fill(impulsotrasversofiglia1);	    
	        }
	  }
	  
       }
   }//fine ciclo sulle tracce
   impulsototfiglia[0] = pow((pow(pxtotfiglia[0],2) + pow(pytotfiglia[0],2) + pow(pztotfiglia[0],2)),0.5);
   impulsototfiglia[1] = pow((pow(pxtotfiglia[1],2) + pow(pytotfiglia[1],2) + pow(pztotfiglia[1],2)),0.5);

   // impulsotrasversofiglia[0] = impulsototfiglia[0] * TMath::Sin(TMath::ATan(TMath::Sqrt(pow(pxtotfiglia[0]/pztotfiglia[0] - charmpx[0]/charmpz[0],2) + pow(pytotfiglia[0]/pztotfiglia[0] - charmpy[0]/charmpz[0],2))));
   // impulsotrasversofiglia[1] = impulsototfiglia[1] * TMath::Sin(TMath::ATan(TMath::Sqrt(pow(pxtotfiglia[1]/pztotfiglia[1] - charmpx[1]/charmpz[1],2) + pow(pytotfiglia[1]/pztotfiglia[1] - charmpy[1]/charmpz[1],2))));
   impulsotrasversofiglia[0] = impulsototfiglia[0] * TMath::Sin(TMath::ACos((charmpx[0] * pxtotfiglia[0] + charmpy[0] * pytotfiglia[0] + charmpz[0] * pztotfiglia[0])/(impulsocharm[0] * impulsototfiglia[0])));
   impulsotrasversofiglia[1] = impulsototfiglia[1] * TMath::Sin(TMath::ACos((charmpx[1] * pxtotfiglia[1] + charmpy[1] * pytotfiglia[1] + charmpz[1] * pztotfiglia[1])/(impulsocharm[1] * impulsototfiglia[1])));
   kink[0] = TMath::ACos((charmpx[0] * pxtotfiglia[0] + charmpy[0] * pytotfiglia[0] + charmpz[0] * pztotfiglia[0])/(impulsocharm[0] * impulsototfiglia[0]));
   kink[1] = TMath::ACos((charmpx[0] * pxtotfiglia[0] + charmpy[0] * pytotfiglia[0] + charmpz[0] * pztotfiglia[0])/(impulsocharm[0] * impulsototfiglia[0]));
       // cout<<impulsotrasversofiglia1<<" "<<impulsototfiglia1<<endl;
     if (nfiglia[0] == 1){
       hpf1_T->Fill(impulsotrasversofiglia[0]);
       hkinkangle1->Fill(kink[0]);
       if ((impulsotrasversofiglia[0]) == 0) cout<<i<<endl;
     }
     if (nfiglia[1] == 1){
      hpf1_T->Fill(impulsotrasversofiglia[1]);
      hkinkangle1->Fill(kink[1]);     
     }
     if (nfiglia[0] == 3){
       hpf3_T->Fill(impulsotrasversofiglia[0]);
     }
       if (nfiglia[1] == 3) hpf3_T->Fill(impulsotrasversofiglia[1]);
 }//fine ciclo sugli eventi
 //tree->Write();
 TCanvas *ckink = new TCanvas();
 hkinkangle1->Draw();
 TCanvas *c0 = new TCanvas();
 hpf1_T->Draw();
 TCanvas *c02 = new TCanvas();
 hpf3_T->Draw();
 TCanvas *c10 = new TCanvas();
 hp->SetTitle("Impulso madre di charm");
 hp->GetXaxis()->SetTitle("GeV");
 hp->Draw();
 TCanvas *c = new TCanvas();
 hpz->GetXaxis()->SetTitle("Gev");
 hpz->Draw();
 TCanvas *c1 = new TCanvas();
 hvz->Draw();
 hvz->Scale(1./hvz->Integral());
 hvz->GetXaxis()->SetTitle("cm");
 //hvz->Fit("expo");
 //TCanvas *c2 = new TCanvas();
 hvz1->SetLineColor(kRed);
 hvz1->Scale(1./hvz1->Integral());
 hvz1->Draw("SAMES");
 hvz1->GetXaxis()->SetTitle("cm");
 
 TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
 leg->AddEntry("hvzprimary", "Primaries");
 leg->AddEntry("hvzsecondary", "Secondaries");
 leg->Draw();

 TCanvas *c3 = new TCanvas();
 hpc1->GetXaxis()->SetTitle("GeV/c");
 //TCanvas *c3_2 = new TCanvas();
 hpc2->SetLineColor(kRed);
 hpc2->SetLineWidth(2);
 hpc2->SetLineStyle(2);
 hpc1->Scale(1./hpc1->Integral());
 hpc2->Draw();
 hpc2->Scale(1./hpc2->Integral());
 hpc1->Draw("SAMES");
 hpc2->GetXaxis()->SetTitle("GeV/c");

 TLegend *leg1 = new TLegend(0.1,0.7,0.48,0.9);
 leg1->AddEntry("hpprimary", "Primaries");
 leg1->AddEntry("hpsecondary", "Secondaries");
 leg1->Draw();

 TCanvas *c4 = new TCanvas();
 TAxis *Xaxis = hwho->GetXaxis();
 hwho->SetTitle("");
 Xaxis->SetBinLabel(1,"pi+"); //inserisco i nomi delle particelle sugli assi
 Xaxis->SetBinLabel(2,"pi-");
 Xaxis->SetBinLabel(3,"p");
 Xaxis->SetBinLabel(4,"pbar");
 Xaxis->SetBinLabel(5,"n");
 Xaxis->SetBinLabel(6,"nbar");
 Xaxis->SetBinLabel(7,"k+");
 Xaxis->SetBinLabel(8,"k-");
 Xaxis->SetBinLabel(9,"ks0");
 Xaxis->SetBinLabel(10,"kl0");
 hwho->Draw();

 TCanvas *c5 = new TCanvas();
 hl->Draw();

 TCanvas *c6 = new TCanvas();
 hpf1->Scale(1./hpf1->Integral());
 hpf2->Scale(1./hpf2->Integral());
 hpf2->SetLineColor(kRed);
 hpf2->SetLineWidth(2);
 hpf2->SetLineStyle(2);
 hpf2->Draw();
 hpf2->GetXaxis()->SetTitle("GeV/c");
 hpf1->Draw("SAMES");
 hpf1->GetXaxis()->SetTitle("GeV/c");

 TLegend *leg2 = new TLegend(0.1,0.7,0.48,0.9);
 leg2->AddEntry("hpf1", "Primaries");
 leg2->AddEntry("hpf2", "Secondaries");
 leg2->Draw(); 

 TCanvas *c7 = new TCanvas();
 hax1->Scale(1./hax1->Integral());
 hax1->GetXaxis()->SetTitle("rad");
 hax2->Scale(1./hax2->Integral());
 hax2->SetLineColor(kRed);
 hax2->SetLineWidth(2);
 hax2->SetLineStyle(2);
 hax1->Draw();
 hax2->Draw("SAMES");
 hax2->GetXaxis()->SetTitle("rad");
 
 TLegend *leg3 = new TLegend(0.1,0.7,0.48,0.9);
 leg3->AddEntry("hax1", "Primaries");
 leg3->AddEntry("hax2", "Secondaries");
 leg3->Draw();

 TCanvas *c8 = new TCanvas();
 hay1->Scale(1./hay1->Integral());
 hay2->Scale(1./hay2->Integral());
 hay2->SetLineColor(kRed);
 hay2->SetLineWidth(2);
 hay2->SetLineStyle(2);
 hay1->Draw();
 hay1->GetXaxis()->SetTitle("rad");
 hay2->Draw("SAMES");
 hay2->GetXaxis()->SetTitle("rad");
 
 TLegend *leg4 = new TLegend(0.1,0.7,0.48,0.9);
 leg4->AddEntry("hay1", "Primaries");
 leg4->AddEntry("hay2", "Secondaries");
 leg4->Draw(); 

 cout<<"Numero di charm primari: "<<primari<<endl;
 cout<<"Numero di charm da cascata: "<<secondari<<endl;
}
bool ischarm(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 431) || (abs(PdgCode) == 411) || (abs(PdgCode) == 4122)  || (abs(PdgCode) == 421) || (abs(PdgCode) == 4132) || (abs(PdgCode) == 4232) ||(abs(PdgCode) == 4332)|| (PdgCode == 441))   check = true; 
      return check;
      }
//lista mesoni charmati: D0(bar), D+(-), Ds0(bar), Ds+(-), Lambdac+(-), Csi+(-), Csi0(bar), Omega_co(bar), eta_c

bool isintermediate(Int_t PdgCode){
  bool check = false;
  if ((abs(PdgCode) == 223) || (abs(PdgCode) == 3222) || (abs(PdgCode) == 3224) || (abs(PdgCode) == 331) || (abs(PdgCode) == 221) || (abs(PdgCode) == 20213) || (abs(PdgCode) == 3212) ||(abs(PdgCode)==213) || (abs(PdgCode) == 113) || (abs(PdgCode) == 2224) || (abs(PdgCode) == 323)) check =  true;
  return check;
  }

double whatmotherofcharm(Int_t PdgCode){
  double codice = -1; //valore di sicurezza
  if (PdgCode == 211) codice = 0.5; //pi+
  if (PdgCode == -211) codice = 1.5; //pi-
   if (PdgCode == 2212) codice = 2.5; //p
  if (PdgCode == -2212) codice = 3.5; //pbar
  if (PdgCode == 2112) codice = 4.5; //n
  if (PdgCode == -2112) codice = 5.5; //nbar
  if (PdgCode == 321) codice = 6.5; //k+
  if (PdgCode == -321) codice = 7.5; //k-
  if (PdgCode == 310) codice = 8.5; //ks
  if (PdgCode == 130) codice = 9.5; //kl
  if (codice == -1) cout <<"Non hai considerato la particella con codice "<< PdgCode<<endl;
  return codice;
}

//lista intermedi: omega, sigma+(-), sigma*, eta', eta, a1,pho+-,pho0, K*+-, delta++, sigma0
