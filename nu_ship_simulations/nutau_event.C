//code for studing tau neutrino simulation and evaluating efficiencies for my PhD thesis (created on 22 June 2020 by A.Iuliano)
//defined functions
TVector3 NeutrinoVertexCoordinates(const ShipMCTrack &track);
Bool_t FindBrick(Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate);

int DecayChannel(vector<int> &daughters,TTreeReaderArray<ShipMCTrack> &tracks, int nnue, int nnumu, int ievent);

bool NeutrinoVertexLocation(int trackID, const ShipMCTrack& track, vector<int> &primaryvisible, ROOT::RVec<int> signalpdgs);
bool GeometricalEfficiency(TVector3 Vn, double offsetxy, int Nminplates);
bool TauDecay(int trackID, int tauid, const ShipMCTrack &track, TVector3 Vn, double tautx, double tauty, double& taudecaylength);

void smearing (double &Xpos, double &Ypos, double spaceres);
vector<double> eff_formula(int found, int total);
//GLOBAL VARIABLES AND HISTOGRAMS
 //*********DEFINITION OF HISTOGRAMS*********************//
 TH1I *hnfilmtau = new TH1I("hnfilmtau","Number of emulsion films passed by a tau lepton;Nfilms",15,0,15);
 TH2D *htauppt = new TH2D("htauppt","Transverse momentum vs momentum of tau lepton;P[GeV/c];Pt[GeV/c]",400,0,400,100,0,10);
 TH1D *htaugamma = new TH1D("htaugamma","Gamma of tau lepton;#gamma",100,0,100);

 TH2D *hvxy = new TH2D("hvxy","Transverse position of vertices",80,-40,40,80,-40,40);
 TH1D *hvz =  new TH1D("hvz","Transverse position of vertices",300,-3300,-3000);
 
 TH1I *hdaughters = new TH1I("hdaughters","Number of daughters per decay", 5, 0,5);
 TH1I *hvisibledaughters = new TH1I("hvisibledaughters","Number of visible daughters per decay", 5, 0,5);

 TH1D *hdl = new TH1D("hdl","Tau Decay length;dl[mm]",300,0,30);
 TH1D *hkink = new TH1D("hkink","Kink angle;#Theta[rad]",50,0,0.5);
 TH1D *hip = new TH1D("hip","Impact Parameter;IP[#mum]",100,0,1000);

 TH1I *hchannel = new TH1I("hchannel","Tau lepton decay channel;IChannel",4,1,5);

 TH1D *hmuonpall = new TH1D("hmuonpall","Momentum of all muons from tau lepton decay;P[GeV/c]",400,0,400);
 TH1D *hmuonangleall = new TH1D("hmuonangleall","Angle of all muons;#Theta[rad]",100,0,1);
 TH1D *hmuonpentered = new TH1D("hmuonpentered","Momentum of muons entered in first RPC;P[GeV/c]",400,0,400);
 TH1D *hmuonangleentered = new TH1D("hmuonangleentered","Angle muons entered in first RPC;#Theta[rad]",100,0,1);
 
 TH1D *helectron_E = new TH1D("helectron_E","Energy of electrons from tau decay;E[GeV]",100,0,100);
 TH1D *helectron_NPlate = new TH1D("helectron_NPlate","Number of plate where the electron is produced;Nplate",57,1,58);
 
 TH1D *hxsagitta = new TH1D("hxsagitta","Sagitta in the xzplane;sxz[cm]",20,-0.1,0.1);
 TH1D *hysagitta_positive = new TH1D("hysagitta_positive","Sagitta for positive particles;syz[cm]",30,-15,15);
 TH1D *hysagitta_negative = new TH1D("hysagitta_negative","Sagitta for negative particles;syz[cm]",30,-15,15);
 
 //**VARIABLES***//

 const double dxtarget = 40.5;
 const double dytarget = 40.5;

 TDatabasePDG *pdg = TDatabasePDG::Instance();

//start main script
void nutau_event(){
 bool tausim = true;
 ROOT::RVec<int> signalpdgs; //particles I want to study the decay
 if (tausim) signalpdgs = {15}; //tau lepton
 else signalpdgs = {411, 431, 4122, 421, 4132, 4232, 4332, 441}; //adding 4122 makes program crash
 //getting tree and defining arrays
 TChain treechain("cbmsim"); //I add a sim of tau and antitaus
 treechain.Add("/home/utente/Simulations/tauneutrino_19June2020/ship.conical.Genie-TGeant4.root"); 
 treechain.Add("/home/utente/Simulations/tauantineutrino_22September2020/ship.conical.Genie-TGeant4.root");
 
 //if (!file) return;
 TTreeReader reader(&treechain);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");
 TTreeReaderArray<HptPoint> dtpoints(reader,"HptPoint");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");
 
 int trackID, detID, ntauhits;
 double energy, mass;

 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import("/home/utente/Simulations/tauneutrino_19June2020/geofile_full.conical.Genie-TGeant4.root");

 //const int nentries = 10000; //test with fewer entries
 const int nentries = reader.GetEntries(true);
 cout<<"Number of events"<<nentries<<endl;

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

 int nnue, nnumu;

 //counters for decay channel
 int whichchannel;

 double tauendx, tauendy, tauendz;
 double taudecaylength;

 int itrack, tauid;

 //DT spectrometer sagitta check
 const int nstations = 5;
 const int maxndaughters = 5;
 double XDT_daughter[maxndaughters][nstations],YDT_daughter[maxndaughters][nstations],ZDT_daughter[maxndaughters][nstations]; 
 const double DTspaceres = 100 * 1e-4; //50 micron in cm
 //muon rpc check
 bool muoninfirstrpc;

 //***********************************START OF MAIN LOOP*************************//
 for(int ientry = 0;ientry<nentries;ientry++){    
     //resetting counters
     isgeometrical = false;
     islocated = false;
     decaysearch = false;

     nnue = 0;
     nnumu = 0;

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


     //********************************FIRST CONDITION: GEOMETRICAL SELECTION********/
     isgeometrical = GeometricalEfficiency(Vn,trans_mindist, Nminplates);

     //******************access the array of tracks
     itrack = 0;
     tauid = -1;
     for (const ShipMCTrack& track: tracks){        
         int pdgcode = track.GetPdgCode();
         double charge = 0.;
         if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();

         double tx = track.GetPx()/track.GetPz();
         double ty = track.GetPy()/track.GetPz();
         double tantheta = TMath::Sqrt(tx * tx + ty * ty);

         double tautx, tauty, taup, taupt;
         int signalparticle = signalpdgs[signalpdgs==TMath::Abs(pdgcode)].size();
         if ( signalparticle > 0 && track.GetMotherId()==0 && tauid < 0){           
           tauid = itrack;
           taup = track.GetP();
           taupt = track.GetPt();
           tautx = track.GetPx()/track.GetPz();
           tauty = track.GetPy()/track.GetPz();
           //filling histograms about tau lepton
           htauppt->Fill(track.GetP(), track.GetPt());
           //computing gamma factor E/m
           energy = track.GetEnergy();           
           mass = pdg->GetParticle(pdgcode)->Mass();
           htaugamma->Fill(energy/mass);           
         }
         if (ientry == 975) cout<<"Start of vertex function for event "<<ientry<<" track "<<itrack<<endl;
         //look for charged particles from primary vertex
         NeutrinoVertexLocation(itrack, track, primaryvisible, signalpdgs);
         //look for charged tracks from tau decay lengths
         if (ientry == 975) cout<<"Start tau decay  function for event "<<ientry<<" track "<<itrack<<endl;

         if(track.GetMotherId()==tauid && TMath::Abs(charge)==0){ //counting neutrinos
             if (TMath::Abs(pdgcode) == 12) nnue++;
             else if (TMath::Abs(pdgcode) == 14) nnumu++;
         } 
         if(track.GetMotherId()==tauid && TMath::Abs(charge)>0){ //daughter of tau/charm particle, adding to daughter list
          daughters.push_back(itrack);
          //check if decay is visible
          if (TauDecay(itrack, tauid, track, Vn, tautx, tauty, taudecaylength)) visibledaughters.push_back(itrack);
         }
         if (ientry == 975) cout<<"End of vertex function for event "<<ientry<<" track "<<itrack<<endl;

         itrack++;        
     } //end of track loop
     if (ientry == 975) cout<<"End track loop "<<endl;
     hdaughters->Fill(daughters.size());
     hvisibledaughters->Fill(visibledaughters.size());
     //decay channel identification
     whichchannel = DecayChannel(daughters,tracks,nnumu,nnue,ientry);
     //else whichchannel = 3;
     hchannel->Fill(whichchannel);
     hdl->Fill(taudecaylength*10);
     
     //a visible electron from tau decay
     if (whichchannel == 2 && isgeometrical && visibledaughters.size() == 1){

      int electronid = visibledaughters[0];
      //in which plate the electron was produced?
      int NeWall, NeRow, NeColumn;
      int NePlate = -1;
      FindBrick(tracks[electronid].GetStartX(), tracks[electronid].GetStartY(), tracks[electronid].GetStartZ(), NeWall, NeRow, NeColumn, NePlate);
      if (NePlate > 0) helectron_NPlate->Fill(NePlate);
      //if (NePlate >0 && NePlate< 20) cout<< "Electron at event "<<ientry<<" track ID "<<electronid<<"with energy "<<tracks[electronid].GetEnergy()<<"  starting at z "<<tracks[electronid].GetStartZ()<<" vz neutrino "<<Vn(2)<<endl;
      //electron energy
      helectron_E->Fill(tracks[electronid].GetEnergy());
      
     }
     
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
     //resetting indeces to be filled with daughters positions in DT stations
     for (int idaughter = 0; idaughter < visibledaughters.size(); idaughter++){
      for (int index = 0; index < nstations; index++){
       XDT_daughter[idaughter][index] = 100000.;
       YDT_daughter[idaughter][index] = 100000.;
       ZDT_daughter[idaughter][index] = 100000.;
      }
     }
     //starting loop
     for (const HptPoint& dtpoint: dtpoints){
        trackID = dtpoint.GetTrackID();
        detID = dtpoint.GetDetectorID();
        
        int istation = detID/1000;
        //check if track is a daughter from tau decay
        for (int idaughter = 0; idaughter < visibledaughters.size(); idaughter++){
          if (trackID == visibledaughters[idaughter]){
          XDT_daughter[idaughter][istation-1] = dtpoint.GetX();
          YDT_daughter[idaughter][istation-1] = dtpoint.GetY();
          ZDT_daughter[idaughter][istation-1] = dtpoint.GetZ();
          
          smearing(XDT_daughter[idaughter][istation-1], YDT_daughter[idaughter][istation-1], DTspaceres);
          }
         }
        
        }
        
     //computing sagitta
     for (int idaughter = 0; idaughter < visibledaughters.size(); idaughter++){
      double particlecharge = pdg->GetParticle(tracks[visibledaughters[idaughter]].GetPdgCode())->Charge();
      if (ZDT_daughter[idaughter][2] < 100000. && ZDT_daughter[idaughter][1] < 100000. && ZDT_daughter[idaughter][0] < 100000.){
       double ysagitta = ((YDT_daughter[idaughter][0] + YDT_daughter[idaughter][2])/2.) - YDT_daughter[idaughter][1];
       double xsagitta = ((XDT_daughter[idaughter][0] + XDT_daughter[idaughter][2])/2.) - XDT_daughter[idaughter][1];
       hxsagitta->Fill(xsagitta);
       if (particlecharge > 0.) hysagitta_positive->Fill(ysagitta);
       else hysagitta_negative->Fill(ysagitta);      
      }
     }

     if (whichchannel == 1){ //only for muon channel,
      //***************************look for muons in downstream detectors****************//
      for (const ShipRpcPoint& rpcpoint: rpcpoints){
        trackID = rpcpoint.GetTrackID(); 
        detID = rpcpoint.GetDetectorID();
       
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
 cout<<"channel legenda: (1 = mu, 2 = e, 3 = 1h, 4 = 3h)"<<endl;
 for (int ichannel = 0; ichannel < ndecaychannels; ichannel++){
  //computing efficiencies with their values
  vector<double> geomeff = eff_formula(geometricalweight[ichannel],totalweight[ichannel]);
  vector<double> localizedeff = eff_formula(localizedweight[ichannel],totalweight[ichannel]);
  vector<double> decaysearcheff = eff_formula(decaysearchweight[ichannel],totalweight[ichannel]);
  //reporting values and errors
  cout<<"Number of interactions in channel "<<ichannel+1<<" is "<<hchannel->GetBinContent(ichannel+1)<<" Total weigth: "<<totalweight[ichannel]<<endl;
  cout<<"Fraction within fiducial volume: "<<geomeff[0]<<" with error "<<geomeff[1]<<endl;
  cout<<"Fraction of localized vertices: "<<localizedeff[0]<<" with error "<<localizedeff[1]<<endl;
  cout<<"Fraction of decay search: "<<decaysearcheff[0]<<" with error "<<decaysearcheff[1]<<endl;
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

 TCanvas *cndaughters = new TCanvas();
 hdaughters->Draw();
 hvisibledaughters->SetLineColor(kRed);
 hvisibledaughters->Draw("SAMES");
 cndaughters->BuildLegend();

 TCanvas *cdecay = new TCanvas();
 cdecay->Divide(2,1);
 cdecay->cd(1);
 hkink->Draw();
 cdecay->cd(2);
 hdl->Draw();

 TCanvas *celectron = new TCanvas();
 celectron->Divide(2,1);
 celectron->cd(1);
 helectron_E->Draw();
 celectron->cd(2);
 helectron_NPlate->Draw();

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
 
 TCanvas *csagitta = new TCanvas();
 csagitta->Divide(2,1);
 csagitta->cd(1);
 hxsagitta->Draw();
 hxsagitta->Fit("gaus");
 csagitta->cd(2);
 hysagitta_negative->Draw();
 hysagitta_positive->SetLineColor(kRed);
 hysagitta_positive->Draw("SAMES");
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

bool NeutrinoVertexLocation(int trackID, const ShipMCTrack& track, vector<int> &primaryvisible, ROOT::RVec<int> signalpdgs){
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
            if (signalpdgs[signalpdgs==TMath::Abs(pdgcode)].size() == 0) primaryvisible.push_back(trackID); //I do not add the tau lepton or charmed hadron since it decays too soon
            return true;
            //cout<<ientry<<" "<<pdgcode<<" "<<momentum<<" "<<tantheta<<endl;
            }
           else return false;
         }
    else return false;
}

bool TauDecay(int trackID, int tauid, const ShipMCTrack &track, TVector3 Vn, double tautx, double tauty, double& taudecaylength){
    if (tauid < 0) return false; //MCTrackID for tau lepton not yet found, go to next
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
    if(track.GetMotherId()==tauid && TMath::Abs(charge)>0){  
    
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
          return true;
         } 
         else return false;
}

int DecayChannel(vector<int> &daughters,TTreeReaderArray<ShipMCTrack> &tracks, int nnumu, int nnue, int ievent){
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
        cout<<"At event "<<ievent<<" Unexpected channel for tau daughters (muons, electrons, hadrons): "<<nmuons<<" "<<nelectrons<<" "<<nhadrons<<endl;
         if (nnumu > 0){          
          cout<<"At event "<<ievent<<" not stored lepton. Stored to muon channel due to found nu_mu: "<<nnumu<<endl;
          return 1;
         }
         else if (nnue > 0){          
          cout<<"At event "<<ievent<<" not stored lepton. Stored to electron channel due to found nu_e: "<<nnue<<endl;
          return 2;
         }
         else{         
          cout<<"No muon/electron neutrino found, associating it to channel 1h"<<endl;
          return 3;
        }
     }
     return 0;

}


void smearing (double &Xpos, double &Ypos, double spaceres){
 float deltaX = gRandom->Gaus(0,spaceres); //angular resolution, adding a gaussian offset to TX and TY
 float deltaY = gRandom->Gaus(0,spaceres);
 //cout<<TX<<endl;
 Xpos = Xpos + deltaX;
 Ypos = Ypos + deltaY;
}
//standard (approximate for values close to 0 and 1) efficency formula with error
vector<double> eff_formula(int found, int total){
  vector<double> efficiency; //value and error
  
  efficiency.push_back((double) found/total);
  double efferr = TMath::Sqrt(efficiency[0] * (1- efficiency[0])/total);
  efficiency.push_back(efferr);
  
  return efficiency;
  
}