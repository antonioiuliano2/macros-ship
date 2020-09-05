//code for studing tau neutrino simulation and evaluating efficiencies for my PhD thesis (created on 22 June 2020 by A.Iuliano)
//defined functions
TVector3 NeutrinoVertexCoordinates(const ShipMCTrack &track);
Bool_t FindBrick(Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate);

int DecayChannel(vector<int> &daughters,TTreeReaderArray<ShipMCTrack> &tracks);

bool NeutrinoVertexLocation(int trackID, const ShipMCTrack& track, vector<int> &primaryvisible);
bool GeometricalEfficiency(TVector3 Vn, double offsetxy, int Nminplates);
bool TauDecay(int trackID, const ShipMCTrack &track,vector<int> &daughters, vector<int> &visibledaughters, TVector3 Vn, double tautx, double tauty, double& taudecaylength);

//GLOBAL VARIABLES AND HISTOGRAMS
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

 const double dxtarget = 40.5;
 const double dytarget = 40.5;

 TDatabasePDG *pdg = TDatabasePDG::Instance();

//start main script
void nutau_event(){
 //getting tree and defining arrays
 TFile *file = TFile::Open("ship.conical.Genie-TGeant4.root"); 
 //if (!file) return;
 TTreeReader reader("cbmsim",file);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");
 TTreeReaderArray<HptPoint> dtpoints(reader,"HptPoint");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");
 
 int trackID, ntauhits;
 double energy, mass;

 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import("geofile_full.conical.Genie-TGeant4.root");

 cout<<"Number of events"<<reader.GetEntries()<<endl;
 const int nentries = reader.GetEntries();

 const int Nminplates = 4;

 const double trans_mindist = 0.1;

 const int ndecaychannels = 4;
 double totalweight[ndecaychannels] = {0.,0.,0.,0.};
 double geometricalweight[ndecaychannels] = {0.,0.,0.,0.};
 double localizedweight[ndecaychannels] = {0.,0.,0.,0.};
 double decaysearchweight[ndecaychannels] = {0.,0.,0.,0.};

 double muonacceptance = 0;

 bool isgeometrical, islocated, decaysearch;

 //IDs of tracks of interest (Tracks at neutrino vertex other than tau and tau daughters)
 vector<int> primaryvisible;
 vector<int> daughters;
 vector<int> visibledaughters;

 //counters for decay channel
 int whichchannel;

 double tauendx, tauendy, tauendz;
 double taudecaylength;

 int itrack;

 bool muoninfirstrpc;

 //***********************************START OF MAIN LOOP*************************//
 for(int ientry = 0;ientry<nentries;ientry++){
     //resetting counters
     isgeometrical = false;
     islocated = false;
     decaysearch = false;

     //ntauhits = 0;
     //clearing list of tracks of interest
     primaryvisible.clear();
     daughters.clear();
     visibledaughters.clear();

     taudecaylength = -1000.;

     //muon detection variables
     muoninfirstrpc = false;
     //getting entry

     if (ientry%10000 == 0) cout<<"arrived at entry" <<ientry<<endl;
     reader.SetEntry(ientry);// keeps track of the number of event (from 0 to Nevents - 1)
     //event weight
     double eventweight = tracks[0].GetWeight();
     //neutrino interaction coordinates
     TVector3 Vn = NeutrinoVertexCoordinates(tracks[0]);
     //if neutrino interacted off of our target, go to next
     if ( (TMath::Abs(Vn(0))>dxtarget) || (TMath::Abs(Vn(1))>dytarget) ) continue;

     hvxy->Fill(Vn(0),Vn(1));
     hvz->Fill(Vn(2));

     //tau lepton variables
     double taup = tracks[1].GetP();
     double taupt = tracks[1].GetPt();
     double tautx = tracks[1].GetPx()/tracks[1].GetPz();
     double tauty = tracks[1].GetPy()/tracks[1].GetPz();

     //********************************FIRST CONDITION: GEOMETRICAL SELECTION********/
     isgeometrical = GeometricalEfficiency(Vn,trans_mindist, Nminplates);

     //filling histograms about tau lepton
     htauppt->Fill(tracks[1].GetP(), tracks[1].GetPt());
     //computing gamma factor E/m
     energy = tracks[1].GetEnergy();
     mass = pdg->GetParticle(tracks[1].GetPdgCode())->Mass();
     htaugamma->Fill(energy/mass);

     //******************access the array of tracks
     itrack = 0;
     for (const ShipMCTrack& track: tracks){         
         double tx = track.GetPx()/track.GetPz();
         double ty = track.GetPy()/track.GetPz();
         double tantheta = TMath::Sqrt(tx * tx + ty * ty);
         
         //look for charged particles from primary vertex
         NeutrinoVertexLocation(itrack, track, primaryvisible);
         //look for charged tracks from tau decay lengths
         TauDecay(itrack, track, daughters, visibledaughters, Vn, tautx, tauty, taudecaylength);
         itrack++;
     } //end of track loop
     //decay channel identification
     whichchannel = DecayChannel(daughters,tracks);
     hchannel->Fill(whichchannel);
     hdl->Fill(taudecaylength*10);
     if (primaryvisible.size()>0) islocated = true;
     if (visibledaughters.size()>0) decaysearch = true;
     //access the hits: 
     
     /*for (const TargetPoint& targetpoint: targetpoints){
        trackID = targetpoint.GetTrackID(); 
        if (trackID == 1){ //hit from tau lepton
            ntauhits++;
        }
     } //end of hit loop
     hnfilmtau->Fill(ntauhits,eventweight);*/
     //somming value, when efficiency is satisfied
     if (whichchannel > 0){
      totalweight[whichchannel-1] += eventweight;
      if (isgeometrical) geometricalweight[whichchannel-1] += eventweight;
      if (isgeometrical&&islocated) localizedweight[whichchannel-1] += eventweight;
      if (isgeometrical&&islocated&&decaysearch) decaysearchweight[whichchannel-1] += eventweight;
     }

     //****************************look for tau decay daughters in downstream trackers***//
     for (const HptPoint& dtpoint: dtpoints){
         
           }

     if (whichchannel == 1){ //only for muon channel,
      //***************************look for muons in downstream detectors****************//
      for (const ShipRpcPoint& rpcpoint: rpcpoints){
        trackID = rpcpoint.GetTrackID(); 
        int detID = rpcpoint.GetDetectorID();
       
        if (trackID == visibledaughters[0] && detID == 10000){ //hit from tau lepton
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
 gPad->BuildLegend();
 cmuonRPC->cd(2);
 hmuonangleall->Draw();
 hmuonangleentered->SetLineColor(kRed);
 hmuonangleentered->Draw("SAMES");
 gPad->BuildLegend();
}


TVector3 NeutrinoVertexCoordinates(const ShipMCTrack &track){
     //building a tvector3 with neutrino vertex coordinates
     double vx = track.GetStartX();
     double vy = track.GetStartY();
     double vz = track.GetStartZ();
    
     return TVector3(vx,vy,vz);
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

//*************************************FIRST CONDITION: GEOMETRICAL EFFICIENCY*********************************//
bool GeometricalEfficiency(TVector3 Vn, double offsetxy, int Nminplates){

    double vx = Vn(0);
    double vy = Vn(1);
    double vz = Vn(2);

    const int Ntotplates = 56;
    double startbrickx = 0.5;
    double startbricky = 0.5;

    int whichwall, whichrow, whichcolumn; //for findbrick function
    int whichplate = 60;

    FindBrick(vx, vy, vz, whichwall, whichrow, whichcolumn, whichplate);
    //longitudinal edge
    if (whichplate > (Ntotplates - Nminplates +1)) return false;
    //upper transverse edges
    if (TMath::Abs(vx) > (dxtarget - offsetxy)) return false;
    if (TMath::Abs(vy) > (dytarget - offsetxy)) return false;
    //lower transverse edges
    if (TMath::Abs(vx) < (startbrickx + offsetxy)) return false;
    if (TMath::Abs(vy) < (startbricky + offsetxy)) return false;
        
    else return true;

}

bool NeutrinoVertexLocation(int trackID, const ShipMCTrack& track, vector<int> &primaryvisible){
     //********************************SECOND CONDITION: VERTEX LOCATION********/
    const double maxtantheta = 1.;
    const double minmomentum = 1.;
    //getting particle charge and momentum
    double momentum = track.GetP();
    double charge = 0.;
    double pdgcode = track.GetPdgCode();

    double tx = track.GetPx()/track.GetPz();
    double ty = track.GetPy()/track.GetPz();
    double tantheta = TMath::Sqrt(tx * tx + ty * ty);
    if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge(); 
    //is it from the primary neutrino interaction?
    if(track.GetMotherId()==0 && TMath::Abs(charge)>0){
          
          //is it visible?      
          if(tantheta<maxtantheta && momentum > minmomentum){ 
            if (pdgcode!=15) primaryvisible.push_back(trackID); //I do not add the tau neutrino since it decays too soon
            return true;
            //cout<<ientry<<" "<<pdgcode<<" "<<momentum<<" "<<tantheta<<endl;
            }
           else return false;
         }
    else return false;
}

bool TauDecay(int trackID, const ShipMCTrack &track,vector<int> &daughters, vector<int> &visibledaughters, TVector3 Vn, double tautx, double tauty, double& taudecaylength){
    //condition for a visible tau decay daughter
    const double maxdl = 0.4;
    const double minkinkangle = 0.02;
    const double minip = 10e-4;
    const double minmomentum = 0.1;
    const double maxtantheta = 1.;

    double vx = Vn(0);
    double vy = Vn(1);
    double vz = Vn(2);

    //look for charged particles from tau decay
    double momentum = track.GetP();
    double charge = 0.;
    double pdgcode = track.GetPdgCode();
    if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge(); 
    //check if it is a tau decay daughter
    if(track.GetMotherId()==1 && TMath::Abs(charge)>0){  
          daughters.push_back(trackID);
    
          double tx = track.GetPx()/track.GetPz();
          double ty = track.GetPy()/track.GetPz();
          double tantheta = TMath::Sqrt(tx * tx + ty * ty);

          if (pdgcode==13){
              hmuonpall->Fill(momentum);
              hmuonangleall->Fill(TMath::ATan(tantheta));
          }
        
          double tauendx = track.GetStartX();
          double tauendy = track.GetStartY();
          double tauendz = track.GetStartZ();
     
          taudecaylength = TMath::Sqrt(pow(tauendx - vx,2) + pow(tauendy-vy,2)+ pow(tauendz-vz,2)); 
          double kinkangle = TMath::ATan(TMath::Sqrt(pow(tx - tautx,2)+pow(ty - tauty,2)));

          //'''Impact parameter of track with respect to primary vertex (standard transverse definition)'''
 
          double dz = vz - tauendz;
          double ipx = tx * dz + tauendx - vx;
          double ipy = ty * dz + tauendy - vy;
          double ip = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));

          hip->Fill(ip*1e+4);
          hkink->Fill(kinkangle);
          //*************************DECAY SEARCH**************************//
          //cout<<"TEST "<<taudecaylength<<" "<<" "<<ip<<" "<<kinkangle<<" "<<momentum<<" "<<tantheta<<endl;
          if (taudecaylength > maxdl) return false;
          if (ip < minip) return false;
          if (kinkangle < minkinkangle) return false;
          if (momentum < minmomentum || tantheta > maxtantheta) return false;
          //if the code was not returned by all these checks, add the track to the visible list
          visibledaughters.push_back(trackID);

          return true;
         } 
         else return false;
}

int DecayChannel(vector<int> &daughters,TTreeReaderArray<ShipMCTrack> &tracks){
    //tau decay channel identification
    int nmuons = 0;
    int nhadrons = 0;
    int nelectrons = 0; 

    for (auto trackID:daughters){
       double pdgcode = tracks[trackID].GetPdgCode();
       if (TMath::Abs(pdgcode) == 11) nelectrons++;
       else if (TMath::Abs(pdgcode) == 13) nmuons++;
       else nhadrons++;
     }
     if (nmuons == 1) return 1;
     else if (nelectrons == 1) return 2;
     else if (nhadrons == 1) return 3;
     else if (nhadrons > 1) return 4;
     else{          
         cout<<"Unexpected channel for tau daughters (muons, electrons, hadrons): "<<nmuons<<" "<<nelectrons<<" "<<nhadrons<<endl;
         return -1;
     }
     return 0;

}
