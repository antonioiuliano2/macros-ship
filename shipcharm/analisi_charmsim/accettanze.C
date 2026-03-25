//studio sulla simulazione di produzione di charm in cascata (riscritto il 22 Marzo 2020 per salvare un tree invece di istogrammi)
bool ischarm(Int_t PdgCode);
bool isintermediate(Int_t PdgCode);
vector<double> eff_formula(int found, int total);

using namespace ROOT;
void accettanze(TString filename = "" ){ 
 double minP = 0.1; //visibility daughter selection for magnetic spectrometer
 //opening file and activating reader
 TFile *file = TFile::Open("/eos/user/a/aiuliano/public/sims_FairShip/sim_charm/CH1_charmcascade_withdrifttubes_19_09_20/simulation/inECC_ship.conical.Pythia8CharmOnly-TGeant4.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 //const Int_t nevents = 1000;
 const Int_t nevents = reader.GetEntries();

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<PixelModulesPoint> pixelpoints(reader,"PixelModulesPoint");
 TTreeReaderArray<SciFiPoint> scifipoints(reader,"SciFiPoint");
 TTreeReaderArray<MufluxSpectrometerPoint> mufluxpoints(reader,"MufluxSpectrometerPoint");

 //primary vertex histograms
 TH2D *hvxy = new TH2D("hvxy","Primary event distribution;x[cm];y[cm]",60,-3,3,60,-3,3);
 TH1D *hvz = new TH1D("hvz","Primary event distribution;z[cm]",50,121,126);

 TH1D *hchargeddaughterP = new TH1D("hchargeddaughterP","All visible charged daughters;P[GeV/c]",50,0,50);
 TH1D *hfoundpixelP = new TH1D("hfoundpixelP","Charged daughters found in pixel;P[GeV/c]",50,0,50);
 TH1D *hfoundscifiP = new TH1D("hfoundscifiP","Charged daughters found in SciFi;P[GeV/c]",50,0,50);
 TH1D *hfoundmufluxP = new TH1D("hfoundmufluxP","Charged daughters found in Drift Tubes or SciFi;P[GeV/c]",50,0,50);
 
 //molteplicities
 TH1I *hnvisible = new TH1I("nvisible","Number of visible daughters per decay;N",10,0,10);
 TH1I *hnfound = new TH1I("hnfound","Number of found daughters per decay;N",10,0,10);

 TH1I *hnhit_firstSciFi = new TH1I("hnhit_firstSciFi","Number of hits arrived at first SciFi per event;Nhits",40,0,40);
 TH2D *hPstartz_firstSciFi = new TH2D("hPstartz_firstSciFi","Momentum vs startZ of hits in first SciFi; P[GeV/c];z[cm]",100,0,100,800,100,900);
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
            if (abs(charge)>0. && momentum>minP){              
              nvisible[whichcharm-1]++; 
              hchargeddaughterP->Fill(momentum);
              bool foundhitmuflux = false;
              bool foundhitscifi = false;
              bool foundhitpixel = false;

              //these loops are made for every charm daughters!
              //starting loop over pixel modules
              for (const PixelModulesPoint& hitpoint: pixelpoints){ 
                 if (hitpoint.GetTrackID()==mytrackID){ 
                   foundhitpixel = true;
                   hitx.push_back(hitpoint.GetX());
                   hity.push_back(hitpoint.GetY());
                   hitz.push_back(hitpoint.GetZ());
                 }
              }
               //starting loop over SciFihits
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
               //loop on drifttubes
              for (const MufluxSpectrometerPoint &hitpoint: mufluxpoints){
                if (hitpoint.GetTrackID()==mytrackID){
                   foundhitmuflux = true;
                }
              }//end of loop on drift tubes
             //requirements are consecutives (i.e. no sense asking for drift tubes or SciFi if pixels are lacking)
              if (foundhitpixel){
               hfoundpixelP->Fill(momentum);
               if (foundhitscifi) hfoundscifiP->Fill(momentum);
               if (foundhitmuflux || foundhitscifi) hfoundmufluxP->Fill(momentum);  
              }
          } //visibility conditions
         // ntotaldaughters++; //increasing counter to store information for next daughter
        } //check if particle pdg is recognized
       }//charm daughter condition
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

       if (hitpoint.GetTrackID() > -1){
         hPstartz_firstSciFi->Fill(hitmomentum,tracks[hitpoint.GetTrackID()].GetStartZ());
       }

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

 TCanvas *cmolt_startz = new TCanvas();
 hPstartz_firstSciFi->Draw("COLZ");

 TCanvas *cmolt_firstSciFi = new TCanvas();
 hnhit_firstSciFi->Draw();
 TCanvas *chits_firstSciFi = new TCanvas();
 chits_firstSciFi->Divide(1,2);
 chits_firstSciFi->cd(1);
 hxy_firstSciFi->Draw("COLZ");
 chits_firstSciFi->cd(2);
 hP_firstSciFi->Draw();

 int totvisibledaughters = hchargeddaughterP->GetEntries();
 int totfoundpixel = hfoundpixelP->GetEntries();
 int totfoundscifi = hfoundscifiP->GetEntries();
 int totfoundmuflux = hfoundmufluxP->GetEntries();
 vector<double> totpixeleff = eff_formula(totfoundpixel, totvisibledaughters);
 vector<double> totscifieff = eff_formula(totfoundscifi, totvisibledaughters);
 vector<double> totmufluxeff = eff_formula(totfoundmuflux, totvisibledaughters);
 
 cout<<"Magnetic spectrometer "<<endl;
 cout<<"Over a total number of visible decay daughters "<<totvisibledaughters<<endl;
 cout<<"Arriving at pixel "<<totfoundpixel<<" Ratio: "<<totpixeleff[0]<<" pm "<<totpixeleff[1]<<endl;
 cout<<"Arriving at SciFi or DT: "<<totfoundmuflux<<" Ratio: "<<totmufluxeff[0]<<" pm "<<totmufluxeff[1]<<endl;
cout<<"Arriving at SciFi or DT: "<<totfoundscifi<<" Ratio: "<<totscifieff[0]<<" pm "<<totscifieff[1]<<endl;

 TCanvas *cacceptance = new TCanvas();
 cacceptance->Divide(1,2);
 cacceptance->cd(1);
 hchargeddaughterP->Draw();
 hfoundpixelP->SetLineColor(kGreen);
 hfoundpixelP->Draw("SAMES");
 hfoundscifiP->SetLineColor(kRed);
 hfoundscifiP->Draw("SAMES");
 hfoundmufluxP->SetLineColor(kMagenta);
 hfoundmufluxP->Draw("SAMES");
 gPad->BuildLegend();
 
 cacceptance->cd(2);
 TEfficiency *pixeleff = new TEfficiency(*hfoundpixelP,*hchargeddaughterP);
 TEfficiency *scifieff = new TEfficiency(*hfoundscifiP,*hchargeddaughterP);
 TEfficiency *mufluxeff = new TEfficiency(*hfoundmufluxP,*hchargeddaughterP);
 pixeleff->SetTitle("Fraction of decay daughters seen in Pixel");
 scifieff->SetTitle("Fraction of decay daughters seen in SciFi"); 
 mufluxeff->SetTitle("Fraction of decay daughters seen in Drift Tubes or SciFi"); 
 pixeleff->SetLineColor(kGreen);
 scifieff->SetLineColor(kRed);
 pixeleff->Draw();
 scifieff->Draw("SAMES");
 mufluxeff->SetLineColor(kMagenta);
 mufluxeff->Draw("SAMES");
 gPad->BuildLegend();

 TCanvas *cprong = new TCanvas();
 hnfound->SetLineColor(kRed);
 hnfound->Draw();
 hnvisible->Draw("SAMES");
 cprong->BuildLegend();

 graphfile->Close();
} //fine programma

void projection(){
 //opening file and activating reader
 TFile *file = TFile::Open("inECC_ship.conical.Pythia8CharmOnly-TGeant4_dig.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);
 TDatabasePDG *pdg = TDatabasePDG::Instance();
 //const Int_t nevents = 1000;
 const Int_t nevents = reader.GetEntries();

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<PixelModulesPoint> pixelpoints(reader,"PixelModulesPoint");
 TTreeReaderArray<SciFiPoint> scifipoints(reader,"SciFiPoint");

 TH2D *hdR = new TH2D("hdR","Position distance between SciFi hits and corresponding pixel hit projections;dX[cm];dY[cm]",410,-20.5,20.5,110,-5.5,5.5);

 map<int,int> SciFihit_byID;
 vector<int> charged_SciFihits;

 //starting loop into the events
 for (int i = 0; i < nevents; i++){
   reader.SetEntry(i);
   SciFihit_byID.clear();
   charged_SciFihits.clear();
   //starting loop into SciFi to build the map
   for (int ihit = 0; ihit < scifipoints.GetSize();ihit++){
     if(scifipoints[ihit].GetDetectorID()==111){ 
       double charge = 0.;
       int pdgcode = scifipoints[ihit].PdgCode();
       if (pdg->GetParticle(pdgcode)){
         if (TMath::Abs(pdg->GetParticle(pdgcode)->Charge())>0) charged_SciFihits.push_back(ihit);
       }
       SciFihit_byID[scifipoints[ihit].GetTrackID()] = ihit; //saving hits from first station
     }
   }  
   if (scifipoints.GetSize() == 0) SciFihit_byID[-2] = -1; //avoid crash for empty container, dummy value
   //starting loop into pixel
   for (const PixelModulesPoint& pixelhit : pixelpoints){
     int trackID = pixelhit.GetTrackID();
     if(trackID > -1) {
       if (SciFihit_byID.count(trackID)){
        int whichhit = SciFihit_byID[trackID];
        double scifihitz = scifipoints[whichhit].GetZ();
        double scifihity = scifipoints[whichhit].GetY();
        double scifihitx = scifipoints[whichhit].GetX();
        //computing TX and TY for projections
        double pixeltx = pixelhit.GetPx()/pixelhit.GetPz();
        double pixelty = pixelhit.GetPy()/pixelhit.GetPz();

        double deltaz = scifihitz - pixelhit.GetZ();
        //computing projections
        double projectedhitx = pixeltx * deltaz + pixelhit.GetX();
        double projectedhity = pixelty * deltaz + pixelhit.GetY();
//        cout<<scifihitx- projectedhitx<<" "<<pixelhit.PdgCode()<<endl;
        hdR->Fill(scifihitx-projectedhitx,scifihity-projectedhity);


       }
     }
   }//end loop on pixel hits
 }//end loop on events

 hdR->Draw("COLZ");
}//end program

vector<double> eff_formula(int found, int total){
  vector<double> efficiency; //value and error
  
  efficiency.push_back((double) found/total);
  double efferr = TMath::Sqrt(efficiency[0] * (1- efficiency[0])/total);
  efficiency.push_back(efferr);
  
  return efficiency;
  
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
