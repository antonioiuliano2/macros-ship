//studio sulla simulazione di produzione di charm in cascata (riscritto il 22 Marzo 2020 per salvare un tree invece di istogrammi)
//distribtuzioni_charm() to create a new file
//plotdistributions() to plot the distrivurions from the file

bool ischarm(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);

using namespace ROOT;
void distribuzioni_charm(TString filename = "" ){ 
 //opening file and activating reader
 TFile *file = TFile::Open("inECC_ship.conical.Pythia8CharmOnly-TGeant4_dig.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");

 //declaring variables
 TDatabasePDG *pdg = TDatabasePDG::Instance();
 Double_t pionmass = pdg->GetParticle(211)->Mass(); //used for approx invariant mass

 Double_t dx[2];
 Double_t dy[2];
 Double_t dz[2];
 Double_t gamma[2];
 Double_t startz = 0;
 Double_t starty = 0;
 Double_t startx = 0;
 Double_t endz = 0;
 Double_t endy = 0;
 Double_t endx = 0;
 const Int_t nevents = reader.GetEntries(); 
 //const Int_t nevents =100000;
  
 Int_t startpdg = 0;
 Int_t intermediatepdg = 0;
 Int_t mumID = 0;
 Int_t mumpdg = 0;
 Int_t charmpdg = 0;
 Int_t hitID = 0;

 //declaring branches for next tree
 Int_t pdgcharm[2];
 Double_t charmpx[2];
 Double_t charmpy[2];
 Double_t charmpz[2];
 Int_t charmlength[2];

 //variable for invariant mass
 Double_t Pxtot[2];
 Double_t Pytot[2];
 Double_t Pztot[2];
 Double_t Etot[2]; 
 Double_t ApproxEtot[2];
 Double_t massinv[2];
 Double_t approxmassinv[2];
// Double_t impulsocharm[2];

 Int_t nprong[2]; //numero prong dei 2 charm
 Int_t nvisible[2]; 
 Int_t ntotaldaughters[2]; //numero totale figlie dei 2 charm, anche neutri
 Int_t longdecay[2]; //did it decay in a plate after the production?
 //histogram to know plate from z position
 const Double_t zstart =  121.8880;
 const Double_t zend =  125.5620;
 const Int_t nplates = 28;
 TH1D *hplatez = new TH1D("hplatez","z Positions, bins as different plates",nplates,  zstart,  zend);
 //variables for charm daughters
 /*const int maxndaughters = 1000;
 Int_t daughpdg[maxndaughters];
 Double_t daughcharge[maxndaughters];
 Int_t daughcharmid[maxndaughters];
 Double_t daughpx[maxndaughters];
 Double_t daughpy[maxndaughters];
 Double_t daughpz[maxndaughters];

 Int_t recognizedpdg[maxndaughters];*/
 Double_t daughmomentum;

 vector<Int_t> daughpdg;
 vector<Double_t> daughcharge;
 vector<Int_t> daughcharmid;
 vector<Double_t> daughpx;
 vector<Double_t> daughpy;
 vector<Double_t> daughpz;

 vector<Int_t> recognizedpdg;

 //****************PREPARING TREE BRANCHES************************//
 TFile *distfile = new TFile("distributions_mctrue_withdaughters_longshort.root","RECREATE");
 TTree *charmlongntuple = new TTree("charmdecays","Charm Decays");//tree structure, arrays containing the information for the two charm decays
 //charm branches
 charmlongntuple->Branch("pdgcode",pdgcharm,"pdgcode[2]/I"); //charm identity
 
 charmlongntuple->Branch("px",charmpx,"px[2]/D"); //momentum information
 charmlongntuple->Branch("py",charmpy,"py[2]/D");
 charmlongntuple->Branch("pz",charmpz,"pz[2]/D");
 charmlongntuple->Branch("gamma",gamma,"gamma[2]/D");
 
 charmlongntuple->Branch("dx",dx,"dx[2]/D"); //decay length
 charmlongntuple->Branch("dy",dy,"dy[2]/D");
 charmlongntuple->Branch("dz",dz,"dz[2]/D");
 
 charmlongntuple->Branch("massinv",massinv,"massinv[2]/D");
 charmlongntuple->Branch("approxmassinv",approxmassinv,"approxmassinv[2]/D");
 charmlongntuple->Branch("longdecay",longdecay,"longdecay[2]/I"); //decay in same plate as production
 
 charmlongntuple->Branch("nprong",nprong,"nprong[2]/I");
 charmlongntuple->Branch("nvisible",nvisible,"nvisible[2]/I");
 charmlongntuple->Branch("ntotaldaughters",ntotaldaughters,"ntotaldaughters[2]/I");
 //charm daughter branches
 charmlongntuple->Branch("daugh_charmid",&daughcharmid); 
 charmlongntuple->Branch("daugh_pdg",&daughpdg); 
 charmlongntuple->Branch("daugh_charge",&daughcharge); 

 charmlongntuple->Branch("recognizedpdg",&recognizedpdg); 

 charmlongntuple->Branch("daugh_px",&daughpx); //momentum information
 charmlongntuple->Branch("daugh_py",&daughpy);
 charmlongntuple->Branch("daugh_pz",&daughpz);

 //charm daughter branches
 //charmlongntuple->Branch("daugh_px",daughpx,"daughpx[ntotaldaughters]/F");

 cout<<"Numero eventi: "<<nevents<<endl;
 for (int i = 0; i < nevents; i++){
  //resetting containers
  for (int icharm = 0; icharm <2;icharm++){

   pdgcharm[icharm] = 0;

   ntotaldaughters[icharm]=0;
   nprong[icharm] = 0; 
   nvisible[icharm] = 0;  

   charmlength[icharm] = 0;   

   Pxtot[icharm] = 0.;
   Pytot[icharm] = 0.;
   Pztot[icharm] = 0.;
   Etot[icharm] = 0.;
   ApproxEtot[icharm] = 0.;
  }
  daughpdg.clear();
  daughcharge.clear();
  daughcharmid.clear();
  daughpx.clear();
  daughpy.clear();
  daughpz.clear();

  recognizedpdg.clear();
 // ntotaldaughters = 0;
  
  if (i % 1000 == 0) cout<<i<<endl;
  reader.SetEntry(i);
   
  for (const ShipMCTrack& trk: tracks){ 

     if (ischarm(trk.GetPdgCode()) && (pdgcharm[0]==0)){  //memorizzo gli indici dei due charm
       pdgcharm[0] = trk.GetPdgCode();
       charmpx[0] = trk.GetPx();
       charmpy[0] = trk.GetPy();
       charmpz[0] = trk.GetPz();
       //impulsocharm[0] = pow((charmpx[0],2) + pow(charmpy[0],2) + pow(charmpz[0],2),0.5);
       Double_t mass = pdg->GetParticle(pdgcharm[0])->Mass();
       gamma[0] = trk.GetEnergy()/mass;
     }
     if (ischarm(trk.GetPdgCode()) && (trk.GetPdgCode() != pdgcharm[0]) && (pdgcharm[1]==0)){ //il secondo lo memorizza solo se diverso dal primo
       pdgcharm[1] = trk.GetPdgCode();
       charmpx[1] = trk.GetPx();
       charmpy[1] = trk.GetPy();
       charmpz[1] = trk.GetPz();
       //impulsocharm[1] = pow((charmpx[1],2) + pow(charmpy[1],2) + pow(charmpz[1],2),0.5);
       Double_t mass = pdg->GetParticle(pdgcharm[1])->Mass();
       gamma[1] = trk.GetEnergy()/mass;
     }
       mumID = trk.GetMotherId(); //prendo l'id della madre di ogni traccia (diretta, non madre della catena)
       mumpdg = trk.GetPdgCode();
       startpdg = trk.GetPdgCode();
       
       endx = trk.GetStartX(); //coordinate di partenza della figlia
       endy = trk.GetStartY();
       endz = trk.GetStartZ(); 
                 
   //    impulsofiglia = pow(pow(trk->GetPz(),2) + pow(trk->GetPx(),2) + pow(trk->GetPy(),2),0.5);
   //    angolofiglia = TMath::ATan(TMath::Sqrt(pow((daughpx/daughpz),2) + pow((daughpy/daughpz),2)));
       //if (i % 100 == 0) cout<<"Il charm è : "<<pdg->GetParticle(startpdg)->GetName()<< " "<<endl;
      
       if (mumID>=0){
	 ShipMCTrack mumtrk = tracks[mumID];
	 mumID = mumtrk.GetMotherId();
      	 mumpdg = mumtrk.GetPdgCode();
	 while ((isintermediate(mumpdg) == true) && (mumID>0)) //se la madre è uno stato intermedio controllo se lo stato intermedio è stato prodotto dal charm
	       {
	        intermediatepdg = mumpdg;
	        ShipMCTrack mumtrk = tracks[mumID];
	        mumID = mumtrk.GetMotherId();
	        mumpdg = mumtrk.GetPdgCode();
	       } 
     	 if((mumpdg == pdgcharm[0]) || (mumpdg == pdgcharm[1])){ //controllo se è la figlia di uno dei due charm

          int whichcharm = -1;
          if (mumpdg == pdgcharm[0]) whichcharm = 1;
          else whichcharm = 2;
          ntotaldaughters[whichcharm-1]++; //checking if it is charged

          daughcharmid.push_back(whichcharm);

          daughpx.push_back(trk.GetPx()); //componenti impulso figlia
          daughpy.push_back(trk.GetPy());
          daughpz.push_back(trk.GetPz());
          daughmomentum = TMath::Sqrt(pow(trk.GetPx(),2)+ pow(trk.GetPy(),2)+ pow(trk.GetPz(),2));


	  startx = mumtrk.GetStartX();
	  starty = mumtrk.GetStartY();
	  startz = mumtrk.GetStartZ();

          //difference with respect to charmed hadron
          dx[whichcharm-1] = endx - startx;
          dy[whichcharm-1] = endy - starty;
          dz[whichcharm-1] = endz - startz;

          //is it long or short?
          int startplate = hplatez->FindBin(startz);
          int endplate = hplatez->FindBin(endz);
          longdecay[whichcharm - 1] = (startplate == endplate)? 0 : 1;

	       // length = TMath::Sqrt(pow(dz,2) + pow(dy,2) + pow(dx,2));

          daughpdg.push_back(startpdg);
          //checking particle charge
          Double_t charge = 0.;
          Int_t recognized = 0;
          if(pdg->GetParticle(startpdg)){ //checking if it is registered
            charge = pdg->GetParticle(startpdg)->Charge()/3.; //TDatabasePDG saves charges as quark units
            recognized = 1;
            if (abs(charge)>0.){ 
	     nprong[whichcharm-1]++; //checking if it is charged
             if(daughmomentum>0.1){
                 nvisible[whichcharm-1]++; //checking if it is charged and with momentum larger than 0.1
		 Pxtot[whichcharm-1] += trk.GetPx();
		 Pytot[whichcharm-1] += trk.GetPy();
		 Pztot[whichcharm-1] += trk.GetPz();
		 Etot[whichcharm-1] += trk.GetEnergy();
         	 //approximated energy, assuming pion mass		 
      		 ApproxEtot[whichcharm-1] += TMath::Sqrt(pow(pionmass,2)+pow(trk.GetP(),2));
		}
             }
          }
          recognizedpdg.push_back(recognized);
          daughcharge.push_back(charge);
         // ntotaldaughters++; //increasing counter to store information for next daughter
       }
      }//condition about mumID
   }//fine ciclo sulle tracce
 for (int icharm = 0; icharm<2;icharm++){
  if (nvisible[icharm]>1){ //no sense otherwise
   massinv[icharm] = TMath::Sqrt(pow(Etot[icharm],2) - (pow(Pxtot[icharm],2)+pow(Pytot[icharm],2)+pow(Pztot[icharm],2)) );
   approxmassinv[icharm] = TMath::Sqrt(pow(ApproxEtot[icharm],2) - (pow(Pxtot[icharm],2)+pow(Pytot[icharm],2)+pow(Pztot[icharm],2)) );
   }
  else{
   massinv[icharm] = -1;
   approxmassinv[icharm] = -1;
   }
  }
 charmlongntuple->Fill();
 }//fine ciclo sugli eventi
 charmlongntuple->Write();
 distfile->Close();
} //fine programma
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

void plotdistributions(){
 //reading charmdecays tree produced at previous step to produce nice distributions
 TFile *charmfile = TFile::Open("alldistributions_mctrue_withdaughters.root");
 TTree *charmdecays = (TTree*) charmfile->Get("charmdecays");

 TDatabasePDG *pdgdatabase = TDatabasePDG::Instance();

 RVec<int> charmpdglist = {421,411,431,4122,4232,4132,4332};
 RVec<int> charmcolors = {kBlack, kRed, kBlue, kYellow, kMagenta, kCyan, kGreen};
 RVec<TString> charmpdgnamelist = {"D^{0}","D^{#pm}","D_{s}^{#pm}","#Lambda_{c}^{#pm}","#Sigma_{c}^{#pm}","#Sigma_{c}^{0}","#Omega_{c}^{0}"};

 if ((charmcolors.size() != charmpdglist.size()) || (charmpdgnamelist.size() != charmpdglist.size())) cout<<"Warning: sizes of RVecs not matching, please check!"<<endl;
 
 map<int,int> charmpdgbin; //for histogram of charm variety
 map<int, TH1D*> charmenergy; //for charm spectrum histogram
 map<int, TH1D*> charmdl; //decay length for each pdg
 
 TH1I *hpdg = new TH1I ("hpdg", "Produced charmed hadron pdg; Particle name", 7,0,7);
 TH1D *hdl = new TH1D ("hdl","Flight length of charmed hadrons;dl[mm]",100,0,100);

 //defining histograms
 for (int icharm = 0; icharm < charmpdglist.size(); icharm++){
  int charmpdg = charmpdglist[icharm];

  charmenergy[charmpdg] = new TH1D(TString::Format("hE%i",charmpdg), (charmpdgnamelist[icharm]+TString(";GeV")).Data(), 40,0,400);
  charmdl[charmpdg] = new TH1D(TString::Format("hdl%i",charmpdg), (charmpdgnamelist[icharm]+TString(";dl[mm]")).Data(), 100,0,100);

  charmpdgbin[charmpdg] = icharm;
  hpdg->GetXaxis()->SetBinLabel(icharm+1,charmpdgnamelist[icharm].Data());
 }
 
 const int npairs = charmdecays->GetEntries();
 const int ncharms = 2;
 //defining variables and setting addresses
 int pdg[ncharms];
 double charge[ncharms];
 double dx[ncharms], dy[ncharms], dz[ncharms];
 double px[ncharms], py[ncharms], pz[ncharms];

 int nocharged = 0, onecharged = 0, bothcharged = 0;
 
 charmdecays->SetBranchAddress("pdgcode",&pdg);
 charmdecays->SetBranchAddress("dx",&dx);
 charmdecays->SetBranchAddress("dy",&dy);
 charmdecays->SetBranchAddress("dz",&dz);
 charmdecays->SetBranchAddress("px",&px);
 charmdecays->SetBranchAddress("py",&py);
 charmdecays->SetBranchAddress("pz",&pz);


 
 for (int ientry = 0; ientry < npairs; ientry++){ //all events in the tree
  charmdecays->GetEntry(ientry);
  for (int icharm = 0; icharm < ncharms; icharm++){ //2 charms for event

   double momentum = TMath::Sqrt(pow(px[icharm],2)+ pow(py[icharm],2) + pow(pz[icharm],2));
   double mass = pdgdatabase->GetParticle(pdg[icharm])->Mass();

   charge[icharm] = pdgdatabase->GetParticle(pdg[icharm])->Charge();
   

   double energy = TMath::Sqrt(momentum * momentum + mass * mass);

   hpdg->Fill(charmpdgbin[TMath::Abs(pdg[icharm])]);
   hdl->Fill(TMath::Sqrt(pow(dx[icharm],2)+ pow(dy[icharm],2) + pow(dz[icharm],2))*10.);

   charmdl[TMath::Abs(pdg[icharm])]->Fill(TMath::Sqrt(pow(dx[icharm],2)+ pow(dy[icharm],2) + pow(dz[icharm],2))*10.);
   charmenergy[TMath::Abs(pdg[icharm])]->Fill(energy);
  }
  if (TMath::Abs(charge[0])>0 && TMath::Abs(charge[1])>0) bothcharged++;
  else if (TMath::Abs(charge[0])>0 || TMath::Abs(charge[1])>0) onecharged++;
  else nocharged++;
 }
 TCanvas *cpdg = new TCanvas();
 hpdg->Draw();
 TCanvas *cdl = new TCanvas();
 hdl->Draw();

 TCanvas *cspectrum = new TCanvas();
 for (int icharm = 0; icharm < charmpdglist.size(); icharm++){
  int charmpdg = charmpdglist[icharm];
  charmenergy[charmpdg]->SetLineColor(charmcolors[icharm]);
 
  if (icharm == 0) charmenergy[charmpdg]->Draw();
  else charmenergy[charmpdg]->Draw("SAME");
 
 }
 cspectrum->BuildLegend();
 
 TCanvas *cdlpdg = new TCanvas();
 for (int icharm = 0; icharm < charmpdglist.size(); icharm++){
  int charmpdg = charmpdglist[icharm];
  charmdl[charmpdg]->SetLineColor(charmcolors[icharm]);
 
  if (icharm == 0) charmdl[charmpdg]->Draw();
  else charmdl[charmpdg]->Draw("SAME");
 
 }
 cdlpdg->BuildLegend();

 cout<<"Over "<<hdl->GetEntries()/2.<<" events we have "<<endl;
 cout<<"Both charged "<<bothcharged<<" one charged "<<onecharged<<" no charged "<<nocharged<<endl;

}
