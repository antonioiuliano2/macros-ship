//code for studing tau neutrino simulation and evaluating efficiencies for my PhD thesis (created on 22 June 2020 by A.Iuliano)
Bool_t FindBrick(Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate);
void nutau_event(){
 //getting tree and defining arrays
 TFile *file = TFile::Open("ship.conical.Genie-TGeant4.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");
 //*********DEFINITION OF HISTOGRAMS*********************//
 TH1I *hnfilmtau = new TH1I("hnfilmtau","Number of emulsion films passed by a tau lepton;Nfilms",15,0,15);
 TH2D *htauppt = new TH2D("htauppt","Transverse momentum vs momentum of tau lepton;P[GeV/c];Pt[GeV/c]",400,0,400,100,0,10);
 TH1D *htaugamma = new TH1D("htaugamma","Gamma of tau lepton;#gamma",100,0,100);

 TH2D *hvxy = new TH2D("hvxy","Transverse position of vertices",80,-40,40,80,-40,40);
 TH1D *hvz =  new TH1D("hvz","Transverse position of vertices",300,-3300,-3000);

 TH1D *hdl = new TH1D("hdl","Tau Decay length;dl[mm]",300,0,30);
 TH1D *hkink = new TH1D("hkink","Kink angle;#Theta[rad]",50,0,0.5);
 TH1D *hip = new TH1D("hip","Impact Parameter;IP[#mum]",100,0,1000);

 TH1I *hchannel = new TH1I("hchannel","Tau lepton decay channel;IChannel",4,1,5);

 TH1D *hmuonpall = new TH1D("hmuonpall","Momentum of all muons from tau lepton decay;P[GeV/c]",400,0,400);
 TH1D *hmuonangleall = new TH1D("hmuonangleall","Angle of all muons;#Theta[rad]",100,0,1);
 TH1D *hmuonpentered = new TH1D("hmuonpentered","Momentum of muons entered in first RPC;P[GeV/c]",400,0,400);
 TH1D *hmuonangleentered = new TH1D("hmuonangleentered","Angle muons entered in first RPC;#Theta[rad]",100,0,1);
 
 int trackID, ntauhits;
 double energy, mass;

 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import("geofile_full.conical.Genie-TGeant4.root");

 TDatabasePDG *pdg = TDatabasePDG::Instance();

 cout<<"Number of events"<<reader.GetEntries()<<endl;
 const int nentries = reader.GetEntries();

 const double dxtarget = 40.5;
 const double dytarget = 40.5;

 const int Ntotplates = 56; //number of passive plates in a brick
 const int Nminplates = 4;

 const double trans_mindist = 0.1;
 const double offset = 0.5; //brick starts at 0.5 and ends at 40.5

 int whichwall, whichrow, whichcolumn, whichplate; //for findbrick function
 const int ndecaychannels = 4;
 double totalweight[ndecaychannels] = {0.,0.,0.,0.};
 double geometricalweight[ndecaychannels] = {0.,0.,0.,0.};
 double localizedweight[ndecaychannels] = {0.,0.,0.,0.};
 double decaysearchweight[ndecaychannels] = {0.,0.,0.,0.};

 double muonacceptance = 0;

 bool isgeometrical, islocated, decaysearch;

 int nprimaryvisible, nprimarytracks, nvisibledaughters;

 //counters for decay channel
 int nmuons, nelectrons, nhadrons;
 int whichchannel;

 double tauendx, tauendy, tauendz;
 double taudecaylength,ip,kinkangle;

 int itrack, muonID;

 bool muoninfirstrpc;

 //***********************************START OF MAIN LOOP*************************//
 for(int ientry = 0;ientry<nentries;ientry++){
     //resetting counters
     isgeometrical = false;
     islocated = false;
     decaysearch = false;

     whichplate = 60;
     ntauhits = 0;
     nprimaryvisible = 0;
     nprimarytracks = 0;
     nvisibledaughters = 0;
 
     whichchannel = 0;
     nmuons = 0;
     nelectrons = 0;
     nhadrons = 0;

     taudecaylength = -1000.;
     ip = -1000.;
     kinkangle = -1000.;

     //muon detection variables
     muoninfirstrpc = false;
     muonID = -10;

     if (ientry%10000 == 0) cout<<"arrived at entry" <<ientry<<endl;
     reader.SetEntry(ientry);// keeps track of the number of event (from 0 to Nevents - 1)
     //event weight
     double eventweight = tracks[0].GetWeight();
     //neutrino interaction coordinates
     double vx = tracks[0].GetStartX();
     double vy = tracks[0].GetStartY();
     double vz = tracks[0].GetStartZ();
     //if neutrino interacted off of our target, go to next
     if ( (TMath::Abs(vx)>dxtarget) || (TMath::Abs(vy)>dytarget) ) continue;
     hvxy->Fill(vx,vy);
     hvz->Fill(vz);
     //where in nutautarget did the interaction happen?
     FindBrick(vx, vy, vz, whichwall, whichrow, whichcolumn, whichplate);


     //tau lepton variables
     double taup = tracks[1].GetP();
     double taupt = tracks[1].GetPt();
     double tautx = tracks[1].GetPx()/tracks[1].GetPz();
     double tauty = tracks[1].GetPy()/tracks[1].GetPz();

     //********************************FIRST CONDITION: GEOMETRICAL SELECTION********/
    if ((whichplate <= (Ntotplates - Nminplates +1)) && 
        (TMath::Abs(vx)<(dxtarget -trans_mindist)) && (TMath::Abs(vy)<(dytarget -trans_mindist)) &&
        (TMath::Abs(vx)>(offset + trans_mindist)) && (TMath::Abs(vy)>(offset + trans_mindist))
        )isgeometrical = true;

     htauppt->Fill(tracks[1].GetP(), tracks[1].GetPt());

     //computing gamma factor E/m
     energy = tracks[1].GetEnergy();
     mass = pdg->GetParticle(tracks[1].GetPdgCode())->Mass();
     htaugamma->Fill(energy/mass);

     //access the array of tracks
     itrack = 0;
     for (const ShipMCTrack& track: tracks){
         double momentum = track.GetP();
         double charge = 0.;
         double pdgcode = track.GetPdgCode();
         if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();          

         double tx = track.GetPx()/track.GetPz();
         double ty = track.GetPy()/track.GetPz();
         double tantheta = TMath::Sqrt(tx * tx + ty * ty);
         //look for charged particles from primary vertex
         if(track.GetMotherId()==0 && TMath::Abs(charge)>0){
          nprimarytracks++;
          //********************************SECOND CONDITION: VERTEX LOCATION********/
          if(tantheta<1. && momentum > 1. && pdgcode!=15 ){ 
            nprimaryvisible++;
            //cout<<ientry<<" "<<pdgcode<<" "<<momentum<<" "<<tantheta<<endl;
            }
         }
         //look for charged particles from tau decay
         if(track.GetMotherId()==1 && TMath::Abs(charge)>0){  

          //which particle is it?
          if (TMath::Abs(pdgcode)==11) nelectrons++; 
          else if (TMath::Abs(pdgcode)==13){ 
              nmuons++;
              muonID = itrack;
              hmuonpall->Fill(momentum);
              hmuonangleall->Fill(TMath::ATan(tantheta));
          } 
          else nhadrons++; 
        

          tauendx = track.GetStartX();
          tauendy = track.GetStartY();
          tauendz = track.GetStartZ();
     
          taudecaylength = TMath::Sqrt(pow(tauendx - vx,2) + pow(tauendy-vy,2)+ pow(tauendz-vz,2)); 
          kinkangle = TMath::ATan(TMath::Sqrt(pow(tx - tautx,2)+pow(ty - tauty,2)));

          //'''Impact parameter of track with respect to primary vertex (standard transverse definition)'''
 
          double dz = vz - tauendz;
          double ipx = tx * dz + tauendx - vx;
          double ipy = ty * dz + tauendy - vy;
          ip = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));

          hip->Fill(ip*1e+4);
          hkink->Fill(kinkangle);
          //*************************DECAY SEARCH**************************//
          if (taudecaylength < 0.4 && ip > 10e-4 && kinkangle > 0.02 && momentum > 0.1 && tantheta < 1.) nvisibledaughters++;

         } 
      itrack++;
     } //end of track loop
     //decay channel identification
     if (nmuons == 1) whichchannel = 1;
     else if (nelectrons == 1) whichchannel = 2;
     else if (nhadrons == 1) whichchannel = 3;
     else if (nhadrons > 1) whichchannel = 4;
     else cout<<"Unexpected channel for event "<<ientry<<" tau daughters (muons, electrons, hadrons): "<<nmuons<<" "<<nelectrons<<" "<<nhadrons<<endl;

     hchannel->Fill(whichchannel);

     if (nprimaryvisible>0) islocated = true;
     if (nvisibledaughters>0) decaysearch = true;
     //access the hits: 
     
     /*for (const TargetPoint& targetpoint: targetpoints){
        trackID = targetpoint.GetTrackID(); 
        if (trackID == 1){ //hit from tau lepton
            ntauhits++;
        }
     } //end of hit loop
     hnfilmtau->Fill(ntauhits,eventweight);*/
     hdl->Fill(taudecaylength*10);
     //somming value, when efficiency is satisfied
     if (whichchannel > 0){
      totalweight[whichchannel-1] += eventweight;
      if (isgeometrical) geometricalweight[whichchannel-1] += eventweight;
      if (isgeometrical&&islocated) localizedweight[whichchannel-1] += eventweight;
      if (isgeometrical&&islocated&&decaysearch) decaysearchweight[whichchannel-1] += eventweight;
     }
     if (whichchannel == 1){ //only for muon channel, look for muons in downstream detectors
      for (const ShipRpcPoint& rpcpoint: rpcpoints){
        trackID = rpcpoint.GetTrackID(); 
        int detID = rpcpoint.GetDetectorID();
       
        if (trackID == muonID && detID == 10000){ //hit from tau lepton
            muoninfirstrpc = true;
            hmuonpentered->Fill(tracks[trackID].GetP());
            double muonangle = 
             TMath::ATan(TMath::Sqrt(pow(tracks[trackID].GetPx()/tracks[trackID].GetPz(),2)+pow(tracks[trackID].GetPy()/tracks[trackID].GetPz(),2)));
            hmuonangleentered->Fill(muonangle);
        }

      } //end of hit loop
      if (isgeometrical&&islocated&&decaysearch&&muoninfirstrpc) muonacceptance += eventweight;
     } //end of muon channel check

 } // end of event loop
 //results
 cout<<"Analying a number of interactions: "<<hvxy->GetEntries()<<endl;
 for (int ichannel = 0; ichannel < ndecaychannels; ichannel++){
  cout<<"Fraction within fiducial volume: "<<geometricalweight[ichannel]/totalweight[ichannel]<<endl;
  cout<<"Fraction of localized vertices: "<<localizedweight[ichannel]/totalweight[ichannel]<<endl;
  cout<<"Fraction of decay search: "<<decaysearchweight[ichannel]/totalweight[ichannel]<<endl;
  cout<<endl;
 }
 cout<<"Fractions of muons from tau decays in first rpc: "<<muonacceptance/totalweight[0]<<endl;
 //***********************DRAWING HISTOGRAMS************************//
 /*TCanvas *cnfilm = new TCanvas();
 hnfilmtau->Draw();

 cout<<"Fraction of short decays"<<hnfilmtau->Integral(1,1)/hnfilmtau->Integral()<<endl;*/
 
 TCanvas *cp = new TCanvas();
 htauppt->Draw("COLZ");

 TCanvas *cgamma = new TCanvas();
 htaugamma->Draw();

 TCanvas *cv = new TCanvas();
 cv->Divide(1,2);
 cv->cd(1);
 hvxy->Draw("COLZ");
 cv->cd(2);
 hvz->Draw();
 
 TCanvas *cchannel = new TCanvas();
 hchannel->Draw();

 TCanvas *cdecay = new TCanvas();
 cdecay->Divide(2,1);
 cdecay->cd(1);
 hkink->Draw();
 cdecay->cd(2);
 hdl->Draw();

 TCanvas *cmuonRPC = new TCanvas();
 cmuonRPC->Divide(2,1);
 cmuonRPC->cd(1);
 hmuonpall->Draw();
 hmuonpentered->SetLineColor(kRed);
 hmuonpentered->Draw("SAMES");
 cmuonRPC->cd(2);
 hmuonangleall->Draw();
 hmuonangleentered->SetLineColor(kRed);
 hmuonangleentered->Draw("SAMES");
}

//find brick from position in the geometry
Bool_t FindBrick(Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate)
{
 gGeoManager->FindNode(x,y,z);
 if (gGeoManager->GetLevel() == 0) return kFALSE;//we have in the cave, no mother volume present
 const char *name = gGeoManager->FindNode(x,y,z)->GetMotherVolume()->GetName(); //go there
 if(strcmp(name, "Brick") == 0){
  NPlate = gGeoManager->GetMother(0)->GetNumber()+1; //volumes numbers start from 0, instead we are used of thinking as 1 is the first plate
  NColumn = gGeoManager->GetMother(2)->GetNumber();
  NRow = gGeoManager->GetMother(3)->GetNumber();
  NWall = gGeoManager->GetMother(4)->GetNumber();

  return kTRUE;
 }
 else return kFALSE;
}
