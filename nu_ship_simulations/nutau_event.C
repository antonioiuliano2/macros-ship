//code for studing tau neutrino simulation and evaluating efficiencies for my PhD thesis (created on 22 June 2020 by A.Iuliano)
Bool_t FindBrick(Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate);
void nutau_event(){
 //getting tree and defining arrays
 TFile *file = TFile::Open("ship.conical.Genie-TGeant4.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");
 //*********DEFINITION OF HISTOGRAMS*********************//
 TH1I *hnfilmtau = new TH1I("hnfilmtau","Number of emulsion films passed by a tau lepton;Nfilms",15,0,15);
 TH2D *htauppt = new TH2D("htauppt","Transverse momentum vs momentum of tau lepton;P[GeV/c];Pt[GeV/c]",400,0,400,100,0,10);
 TH1D *htaugamma = new TH1D("htaugamma","Gamma of tau lepton;#gamma",100,0,100);

 TH2D *hvxy = new TH2D("hvxy","Transverse position of vertices",80,-40,40,80,-40,40);
 TH1D *hvz =  new TH1D("hvz","Transverse position of vertices",300,-3300,-3000);

 TH1D *hdl = new TH1D("hdl","Tau Decay length",300,0,30);
 
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

 double totalweight = 0.;
 double geometricalweight = 0.;
 double localizedweight = 0.;

 bool isgeometrical, islocated;

 int nprimaryvisible, nprimarytracks;

 double tauendx, tauendy, tauendz;
 double taudecaylength;

 //***********************************START OF MAIN LOOP*************************//
 for(int ientry = 0;ientry<nentries;ientry++){
     //resetting counters
     isgeometrical = false;
     islocated = false;

     whichplate = 60;
     ntauhits = 0;
     nprimaryvisible = 0;
     nprimarytracks = 0;

     taudecaylength = -1000.;

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
          tauendx = track.GetStartX();
          tauendy = track.GetStartY();
          tauendz = track.GetStartZ();
     
          taudecaylength = TMath::Sqrt(pow(tauendx - vx,2) + pow(tauendy-vy,2)+ pow(tauendz-vz,2)); 
         } 
     } //end of track loop
     if (nprimaryvisible>0) islocated = true;
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
     totalweight += eventweight;
     if (isgeometrical) geometricalweight += eventweight;
     if (isgeometrical&&islocated) localizedweight += eventweight;

 } // end of event loop
 //results
 cout<<"Analying a number of interactions: "<<hvxy->GetEntries()<<endl;
 cout<<"Fraction within fiducial volume "<<geometricalweight/totalweight<<endl;
 cout<<"Fraction of localized vertices "<<localizedweight/totalweight<<endl;
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

 TCanvas *cdecay = new TCanvas();
 cdecay->Divide(1,2);
 cdecay->cd(1);
 cdecay->cd(2);
 hdl->Draw();
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
