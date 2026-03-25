//script to check distributions after a charm simulation
//The event topology is a primary vertex followed by the decay of two charm hadrons
//to launch it

//declaring functions to be used
bool ischarm(Int_t PdgCode); //recognizing charm hadrons
bool isintermediate(Int_t PdgCode);

//include files need to be changed to refer to your file 
#include "/home/utente/Scrivania/SHIPBuild/FairShip/shipdata/ShipMCTrack.h" 
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-14-00-ship-2/include/TDatabasePDG.h"
#include "/home/utente/Scrivania/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"

void distribuzioni_charm(TString filename = "" ){ 
 TFile *file1;
 if (filename != "") file1 = TFile::Open(filename);
 else{
  //TFile *file1 = TFile::Open("./10run/ship.10.0.Pythia8CharmOnly-TGeant4_primo.root");
 file1 = TFile::Open("ship.conical.Pythia8CharmOnly-TGeant4.root");  
 }
 TTree* cbmsim = (TTree*)file1->Get("cbmsim");

 //adding some particles to pdg to avoid crashes
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
 
 TClonesArray *arr1 = new TClonesArray("ShipMCTrack",1000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr1);

 TH1D *hpz = new TH1D("hpz", "componente z dell'momentum del protone madre di charm", 50, 0, 400);
 TH1D *hp = new TH1D("hp", "momentum protone madre di charm",100, 0, 400);
 TH1D *hpc1 = new TH1D("hpprimary", "momentum charm primari", 100, 0, 350);
 TH1D *hpc2 = new TH1D("hpsecondary", "momentum charm secondari", 100, 0, 350);
 TH1D *hpf1 = new TH1D("hpf1", "momentum figlie di charm primari",100, 0, 40);
 TH1D *hpf2 = new TH1D("hpf2", "momentum figlie di charm secondari",100, 0, 40);
 TH1D *hkinkangle1 = new TH1D("hkinkangle1", "angolo di kink 1 prong",100, 0, 0.6);
 TH1D *hkinkangle3 = new TH1D("hkinkangle3", "angolo di kink 3 prong",100, 0, 0.6);
 
 TH1D *hpf1_T = new TH1D("hpf1_T", "momentum trasverso figlie di charm 1 prong",100, 0, 2);
 TH1D *hpf3_T = new TH1D("hpf3_T", "momentum trasverso figlie di charm 3 prong",100, 0, 2);
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
 Double32_t pz, momentum;
 Double32_t momentumdaughter;
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

 
 Double_t angolodaughter;
 
 Int_t pdgcharm[2];
 Double_t charmpx[2];
 Double_t charmpy[2];
 Double_t charmpz[2];
 Int_t charmlength[2];
 Double_t momentumcharm[2];
 Double_t kink[2];

 Int_t ndaughter[2]; //numero figlie dei 2 charm
 Double_t pxtotdaughter[2];
 Double_t pytotdaughter[2];
 Double_t pztotdaughter[2]; 
 Double_t momentumtotdaughter[2];
 Double_t momentumtrasversodaughter[2];
 
 cout<<"Numero eventi: "<<nevents<<endl;
 for (int i = 0; i < nevents; i++){

  //variable initialization
  parentfound = false;

  pdgcharm[0] = 0;
  pdgcharm[1] = 0;

  ndaughter[0] = 0;
  ndaughter[1] = 0;
  
  pxtotdaughter[0] = 0.;
  pytotdaughter[0] = 0.;
  pztotdaughter[0] = 0.;

  pxtotdaughter[1] = 0.;
  pytotdaughter[1] = 0.;
  pztotdaughter[1] = 0.;

  charmlength[0] = 0;
  charmlength[1] = 0;
  
  if (i % 100 == 0) cout<<i<<endl;
   arr1->Clear();
   cbmsim->GetEntry(i);
   
   ShipMCTrack *trk = (ShipMCTrack*) arr1->At(0);
   //Getting primary vertex
   Vx = trk->GetStartX();
   Vy = trk->GetStartY();
   Vz = trk->GetStartZ();
   cout<<"track pdgcode:" <<trk->GetPdgCode()<<endl;
   
   for (int j = 1; j < arr1->GetEntriesFast();j++){

     ShipMCTrack *trk = (ShipMCTrack*) arr1->At(j);

     //saving information of the charm pair produced
      if (ischarm(trk->GetPdgCode()) && (pdgcharm[0]==0)){  
       pdgcharm[0] = trk->GetPdgCode();
       charmpx[0] = trk->GetPx();
       charmpy[0] = trk->GetPy();
       charmpz[0] = trk->GetPz();
       momentumcharm[0] = pow((charmpx[0],2) + pow(charmpy[0],2) + pow(charmpz[0],2),0.5);
     }
     if (ischarm(trk->GetPdgCode()) && (trk->GetPdgCode() != pdgcharm[0]) && (pdgcharm[1]==0)){ //il secondo lo memorizza solo se diverso dal primo
       pdgcharm[1] = trk->GetPdgCode();
       charmpx[1] = trk->GetPx();
       charmpy[1] = trk->GetPy();
       charmpz[1] = trk->GetPz();
       momentumcharm[1] = pow((charmpx[1],2) + pow(charmpy[1],2) + pow(charmpz[1],2),0.5);
     }
       mumID = trk->GetMotherId(); //taking mother of each track
       mumpdg = trk->GetPdgCode();
       startpdg = trk->GetPdgCode();
       
       if (startpdg > 10000) continue;

       //start of daughter, end of mother
       endx = trk->GetStartX();
       endy = trk->GetStartY();
       endz = trk->GetStartZ(); 

       //momentum components
       startpx = trk->GetPx();
       startpy = trk->GetPy();
       startpz = trk->GetPz();
    
       momentumdaughter = pow(pow(trk->GetPz(),2) + pow(trk->GetPx(),2) + pow(trk->GetPy(),2),0.5);
       angolodaughter = TMath::ATan(TMath::Sqrt(pow((startpx/startpz),2) + pow((startpy/startpz),2)));
       //if (i % 100 == 0) cout<<"Il charm è : "<<pdg->GetParticle(startpdg)->GetName()<< " "<<endl;
      
       if (mumID>=0){
	 trk = (ShipMCTrack*)arr1->At(mumID);
	 mumID = trk->GetMotherId();
	 mumpdg = trk->GetPdgCode();
	 while ((isintermediate(mumpdg) == true) && (mumID>0)) //if mother is an intermediate state, I will check that the intermediate state is a charm
	   {
	     intermediatepdg = mumpdg;
	     trk = (ShipMCTrack*) arr1->At(mumID);
	     mumID = trk->GetMotherId();
	     mumpdg = trk->GetPdgCode();
	   }
	 pz = trk->GetPz();
	 momentum = TMath::Sqrt(pow(trk->GetPz(),2) + pow(trk->GetPx(),2) + pow(trk->GetPy(),2));

	  if((mumpdg == pdgcharm[0]) || (mumpdg == pdgcharm[1])){ 
	    startx = trk->GetStartX();
	    starty = trk->GetStartY();
	    startz = trk->GetStartZ();
	    length = TMath::Sqrt(pow((endz - startz),2) + pow((endy - starty),2) + pow((endx -startx),2));

	    mumID = trk->GetMotherId();
            charmpdg = trk->GetPdgCode();
	    
	    if(mumpdg == pdgcharm[0] && (abs(pdg->GetParticle(startpdg)->Charge()) > 0)){ //daughter of first charm
		 ndaughter[0]++;
		 pxtotdaughter[0] += startpx;
		 pytotdaughter[0] += startpy;
		 pztotdaughter[0] += startpz;		 		 
	    }
	    if((mumpdg == pdgcharm[1]) && (abs(pdg->GetParticle(startpdg)->Charge()) > 0)){ //daughter of second charm
		 ndaughter[1]++;
		 pxtotdaughter[1] += startpx;
		 pytotdaughter[1] += startpy;
		 pztotdaughter[1] += startpz;
	    }
	    trk = (ShipMCTrack*)arr1->At(mumID);
	    mumpdg = trk->GetPdgCode();
	    //cout<<"Madre di charm: "<<mumpdg<<endl;
	    pz = trk->GetPz();
            if ((charmpdg != charmlength[0]) && (charmpdg != charmlength[1])){ //lenght histogram must be filled only once for each charm
	    
		hl->Fill(length);
	
		if (charmlength[0] == 0) charmlength[0] = charmpdg;
		else if (charmlength[1] == 0) charmlength[1] = charmpdg;
	      }

             if(abs(pdg->GetParticle(startpdg)->Charge()) > 0){ //to check only for charged charm daughter
	       // momentumtrasversodaughter1 = momentumdaughter * TMath::Sin(angolodaughter);
	       // hpf_T->Fill(momentumtrasversodaughter1);	    
	        }
	  }
	  
       }
   }//end loop on tracks
   momentumtotdaughter[0] = pow((pow(pxtotdaughter[0],2) + pow(pytotdaughter[0],2) + pow(pztotdaughter[0],2)),0.5);
   momentumtotdaughter[1] = pow((pow(pxtotdaughter[1],2) + pow(pytotdaughter[1],2) + pow(pztotdaughter[1],2)),0.5);
   
   momentumtrasversodaughter[0] = momentumtotdaughter[0] * TMath::Sin(TMath::ACos((charmpx[0] * pxtotdaughter[0] + charmpy[0] * pytotdaughter[0] + charmpz[0] * pztotdaughter[0])/(momentumcharm[0] * momentumtotdaughter[0])));
   momentumtrasversodaughter[1] = momentumtotdaughter[1] * TMath::Sin(TMath::ACos((charmpx[1] * pxtotdaughter[1] + charmpy[1] * pytotdaughter[1] + charmpz[1] * pztotdaughter[1])/(momentumcharm[1] * momentumtotdaughter[1])));
   kink[0] = TMath::ACos((charmpx[0] * pxtotdaughter[0] + charmpy[0] * pytotdaughter[0] + charmpz[0] * pztotdaughter[0])/(momentumcharm[0] * momentumtotdaughter[0]));
   kink[1] = TMath::ACos((charmpx[0] * pxtotdaughter[0] + charmpy[0] * pytotdaughter[0] + charmpz[0] * pztotdaughter[0])/(momentumcharm[0] * momentumtotdaughter[0]));
       // cout<<momentumtrasversodaughter1<<" "<<momentumtotdaughter1<<endl;
     if (ndaughter[0] == 1){
       hpf1_T->Fill(momentumtrasversodaughter[0]);
       hkinkangle1->Fill(kink[0]);
       if ((momentumtrasversodaughter[0]) == 0) cout<<i<<endl;
     }
     if (ndaughter[1] == 1){
      hpf1_T->Fill(momentumtrasversodaughter[1]);
      hkinkangle1->Fill(kink[1]);     
     }
     if (ndaughter[0] == 3){
       hpf3_T->Fill(momentumtrasversodaughter[0]);
     }
       if (ndaughter[1] == 3) hpf3_T->Fill(momentumtrasversodaughter[1]);
 }//fine ciclo sugli eventi
 //tree->Write();
 TCanvas *ckink = new TCanvas();
 hkinkangle1->Draw();
 TCanvas *c0 = new TCanvas();
 hpf1_T->Draw();
 TCanvas *c02 = new TCanvas();
 hpf3_T->Draw();
 TCanvas *c10 = new TCanvas();
 hp->SetTitle("Momentum madre di charm");
 hp->GetXaxis()->SetTitle("GeV");
 hp->Draw();
 TCanvas *c = new TCanvas();
 hpz->GetXaxis()->SetTitle("Gev");
 hpz->Draw();

 TCanvas *c5 = new TCanvas();
 hl->Draw();
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

//lista intermedi: omega, sigma+(-), sigma*, eta', eta, a1,pho+-,pho0, K*+-, delta++, sigma0
