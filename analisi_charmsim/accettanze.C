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
 TTreeReaderArray<SciFiPoint> scifipoints(reader,"SciFiPoint");

 //declaring variables
 TDatabasePDG *pdg = TDatabasePDG::Instance();
 Int_t mytrackID;

 Int_t pdgcharm[2];
 Int_t nvisible[2];

 Int_t mumID, mumpdg, startpdg;
 Double_t momentum;

 cout<<"Numero eventi: "<<nevents<<endl;
 for (int i = 0; i < nevents; i++){
  //resetting containers
  for (int icharm = 0; icharm <2;icharm++){

   pdgcharm[icharm] = 0;
   nvisible[icharm] = 0;  
  }
 // ntotaldaughters = 0;
  
  if (i % 1000 == 0) cout<<i<<endl;
  reader.SetEntry(i);
  mytrackID=0;
  for (const ShipMCTrack& trk: tracks){ 
     mytrackID++;
     if (ischarm(trk.GetPdgCode()) && (pdgcharm[0]==0)){  //memorizzo gli indici dei due charm
       pdgcharm[0] = trk.GetPdgCode();
     }
     if (ischarm(trk.GetPdgCode()) && (trk.GetPdgCode() != pdgcharm[0]) && (pdgcharm[1]==0)){ //il secondo lo memorizza solo se diverso dal primo
       pdgcharm[1] = trk.GetPdgCode();
     }
       mumID = trk.GetMotherId(); //prendo l'id della madre di ogni traccia (diretta, non madre della catena)
       mumpdg = trk.GetPdgCode();
       startpdg = trk.GetPdgCode();      
      
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
        
          momentum = (trk.GetPx(),2 + trk.GetPy(),2 + trk.GetPz(),2);
          //checking particle charge
          Double_t charge = 0.;
          if(pdg->GetParticle(startpdg)){ //checking if it is registered
            charge = pdg->GetParticle(startpdg)->Charge()/3.; //TDatabasePDG saves charges as quark units
            if (abs(charge)>0. && daughmomentum>0.1){ 
              nvisible[whichcharm-1]++; //checking if it is charged and with momentum larger than 0.1
              hchargeddaughterP->Fill(momentum);
              //starting loop over SciFihits
              bool foundhitscifi = false;
               for (const SciFipoint& hitpoint: scifipoints){ 
                 if (hitpoint.TrackID()==mytrackID) foundhitscifi = true;
             }
             if (foundhitscifi) hfoundscifiP->Fill(momentum);
             
          }
          recognizedpdg.push_back(recognized);
          daughcharge.push_back(charge);
         // ntotaldaughters++; //increasing counter to store information for next daughter
       }
      }//condition about mumID
   }//fine ciclo sulle tracce
 }//fine ciclo sugli eventi
 hchargeddaughterP->Draw();
 hfoundscifiP->SetLineColor(kRed);
 hfoundscifiP->Draw("SAMES");
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
