//studio sulla simulazione di produzione di charm in cascata (riscritto il 22 Marzo 2020 per salvare un tree invece di istogrammi)
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
 TFile *distfile = new TFile("distributions_mctrue_withdaughters.root","RECREATE");
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
 //charmlongntuple->Branch("longdecay",longdecay,'longdecay[2]/I'); //decay in same plate as production
 
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
