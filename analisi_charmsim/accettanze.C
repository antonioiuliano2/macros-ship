//studio sulla simulazione di produzione di charm in cascata (riscritto il 22 Marzo 2020 per salvare un tree invece di istogrammi)
bool ischarm(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);

using namespace ROOT;
void accettanze(TString filename = "" ){ 
 //opening file and activating reader
 TFile *file = TFile::Open("inECC_ship.conical.Pythia8CharmOnly-TGeant4_dig.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 //const Int_t nevents = 1000;
 const Int_t nevents = reader.GetEntries();

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<PixelModulesPoint> pixelpoints(reader,"PixelModulesPoint");
 TTreeReaderArray<SciFiPoint> scifipoints(reader,"SciFiPoint");

 //primary vertex histograms
 TH2D *hvxy = new TH2D("hvxy","Primary event distribution;x[cm];y[cm]",60,-3,3,60,-3,3);
 TH1D *hvz = new TH1D("hvz","Primary event distribution;z[cm]",50,121,126);

 TH1D *hchargeddaughterP = new TH1D("hchargeddaughterP","All visible charged daughters;P[GeV/c]",100,0,100);
 TH1D *hfoundpixelP = new TH1D("hfoundpixelP","Charged daughters found in pixel;P[GeV/c]",100,0,100);
 TH1D *hfoundscifiP = new TH1D("hfoundscifiP","Charged daughters found in SciFi;P[GeV/c]",100,0,100);

 //molteplicities
 TH1I *hnvisible = new TH1I("nvisible","Number of visible daughters per decay;N",10,0,10);
 TH1I *hnfound = new TH1I("hnfound","Number of found daughters per decay;N",10,0,10);

 TH1I *hnhit_firstSciFi = new TH1I("hnhit_firstSciFi","Number of hits arrived at first SciFi per event;Nhits",40,0,40);
 TH1D *hP_firstSciFi = new TH1D("hP_firstSciFi","Momentum of hits in first SciFi;P[GeV/c]",100,0,100);
 TH2D *hxy_firstSciFi = new TH2D("hxy_firstSciFi","XY distribution of hits in first SciFi;x[cm];y[cm]",400,-20,20,400,-20,20);

 //declaring variables
 TDatabasePDG *pdg = TDatabasePDG::Instance();
 Int_t mytrackID;

 Int_t pdgcharm[2];
 Int_t nvisible[2];
 Int_t nscififound[2];
 
 Int_t nallfound = 0;
 Int_t natleastonefound = 0;
 Int_t nfirstSciFi = 0;

 Int_t mumID, mumpdg, startpdg,intermediatepdg;
 Double_t momentum;

 vector<Double_t> hitx,hity,hitz;

 cout<<"Numero eventi: "<<nevents<<endl;

 TFile *graphfile = new TFile("hitgraphs.root","RECREATE");

 for (int i = 0; i < nevents; i++){
  //resetting containers
  for (int icharm = 0; icharm <2;icharm++){

   pdgcharm[icharm] = 0;
   nvisible[icharm] = 0;
   nscififound[icharm] = 0;  
  }
  nfirstSciFi = 0;
  hitx.clear();
  hity.clear();
  hitz.clear();
  
  if (i % 100 == 0) cout<<i<<endl;
  reader.SetEntry(i);
  mytrackID=0;
  for (const ShipMCTrack& trk: tracks){ 
     //primary event (i.e. charm production) distribution
     if(mytrackID == 0){
       hvxy->Fill(trk.GetStartX(),trk.GetStartY());
       hvz->Fill(trk.GetStartZ());
     }
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
        
          momentum = trk.GetP();
          //checking particle charge
          Double_t charge = 0.;
          if(pdg->GetParticle(startpdg)){ //checking if it is registered
            charge = pdg->GetParticle(startpdg)->Charge()/3.; //TDatabasePDG saves charges as quark units
            //checking if it is charged and with momentum larger than 0.1
            if (abs(charge)>0. && momentum>0.1){              
              nvisible[whichcharm-1]++; 
              hchargeddaughterP->Fill(momentum);
              //starting loop over SciFihits
              bool foundhitscifi = false;
              bool foundhitpixel = false;

              //these loops are made for every charm daughters!

              for (const PixelModulesPoint& hitpoint: pixelpoints){ 
                 if (hitpoint.GetTrackID()==mytrackID){ 
                   foundhitpixel = true;
                   hitx.push_back(hitpoint.GetX());
                   hity.push_back(hitpoint.GetY());
                   hitz.push_back(hitpoint.GetZ());
                 }
             }
               for (const SciFiPoint& hitpoint: scifipoints){ 
                 //saving charged hits from first station
                 if (hitpoint.GetTrackID()==mytrackID){
                   if(! foundhitscifi) nscififound[whichcharm-1]++; 
                   foundhitscifi = true;
                   hitx.push_back(hitpoint.GetX());
                   hity.push_back(hitpoint.GetY());
                   hitz.push_back(hitpoint.GetZ());                    
                 }
             }
             if (foundhitscifi) hfoundscifiP->Fill(momentum);
             if (foundhitpixel) hfoundpixelP->Fill(momentum);
          }
          }
         // ntotaldaughters++; //increasing counter to store information for next daughter
       }
      }//condition about mumID
      mytrackID++;

   }//fine ciclo sulle tracce
   //ora posso fare un loop sugli hit degli SciFi per studiarne la molteplicità
   for (const SciFiPoint& hitpoint: scifipoints){ 
    if (pdg->GetParticle(hitpoint.PdgCode())&& hitpoint.GetDetectorID()==111){
     if (TMath::Abs(pdg->GetParticle(hitpoint.PdgCode())->Charge())>0){
       Double_t hitmomentum = TMath::Sqrt(pow(hitpoint.GetPx(),2)+pow(hitpoint.GetPy(),2)+pow(hitpoint.GetPz(),2));
       if (hitmomentum > 0.01){//cut to avoid to count production in my detector
       nfirstSciFi++;
       hP_firstSciFi->Fill(hitmomentum);
       hxy_firstSciFi->Fill(hitpoint.GetX(),hitpoint.GetY());
      }
     }
    }
   }
   TGraph hityz = TGraph(hity.size(),hitz.data(),hity.data());
   TGraph hitxz = TGraph(hitx.size(),hitz.data(),hitx.data());
   
   hitxz.SetTitle(Form("XZ distribution in pixel and SciFi of charm daughter for event %i",i));
   hityz.SetTitle(Form("YZ distribution in pixel and SciFi of charm daughter for event %i",i));
 
   hitxz.SetName(Form("Graphxz_%i",i));
   hityz.SetName(Form("Graphyz %i",i));

   hitxz.Write();
   hityz.Write();
   //checking fraction of found events
   for (int icharm = 0; icharm<2; icharm++){
     hnvisible->Fill(nvisible[icharm]);
     hnfound->Fill(nscififound[icharm]);

     if (nvisible[icharm] == nscififound[icharm]){ 
       nallfound++;
       if (i%1000 == 0) cout<<"Found all "<<i<<" "<<nvisible[icharm]<<endl;
      }
     else if (i%1000 == 0) cout<<"Not found all "<<i<<" "<<nvisible[icharm]<<endl;
     if (nscififound[icharm]>=1) natleastonefound++;
   }
   hnhit_firstSciFi->Fill(nfirstSciFi);
 }//fine ciclo sugli eventi

 cout<<"Over decays "<<nevents*2<<endl;
 cout<<"Found all daughters "<<nallfound<<" fraction over total "<<(Double_t) nallfound/(nevents*2)<<endl;
 cout<<"Found at least one daughter "<<natleastonefound<<" fraction over total "<<(Double_t) natleastonefound/(nevents*2)<<endl;

 TCanvas *cprim = new TCanvas();
 cprim->Divide(1,2);
 cprim->cd(1);
 hvxy->Draw("COLZ");
 cprim->cd(2);
 hvz->Draw();

 TCanvas *cmolt_firstSciFi = new TCanvas();
 hnhit_firstSciFi->Draw();
 TCanvas *chits_firstSciFi = new TCanvas();
 chits_firstSciFi->Divide(1,2);
 chits_firstSciFi->cd(1);
 hxy_firstSciFi->Draw("COLZ");
 chits_firstSciFi->cd(2);
 hP_firstSciFi->Draw();


 TCanvas *cacceptance = new TCanvas();
 hchargeddaughterP->Draw();
 hfoundpixelP->SetLineColor(kGreen);
 hfoundpixelP->Draw("SAMES");
 hfoundscifiP->SetLineColor(kRed);
 hfoundscifiP->Draw("SAMES");
 cacceptance->BuildLegend();

 TCanvas *ceff = new TCanvas();
 TEfficiency *pixeleff = new TEfficiency(*hfoundpixelP,*hchargeddaughterP);
 TEfficiency *scifieff = new TEfficiency(*hfoundscifiP,*hchargeddaughterP);

 pixeleff->SetTitle("Fraction of decay daughters seen in Pixel");
 scifieff->SetTitle("Fraction of decay daughters seen in SciFi"); 

 //pixeleff->GetHistogram()->GetXaxis()->SetRangeUser(0,1);
 pixeleff->SetLineColor(kGreen);
 scifieff->SetLineColor(kRed);
 pixeleff->Draw();
 scifieff->Draw("SAMES");
 ceff->BuildLegend();

 TCanvas *cprong = new TCanvas();
 hnfound->SetLineColor(kRed);
 hnfound->Draw();
 hnvisible->Draw("SAMES");
 cprong->BuildLegend();

 graphfile->Close();
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
