#include "EfficiencyCut.h"

//reading the file, prepare TTreeReader
EfficiencyCut::EfficiencyCut(const char* inputfilename, const char* geofilename):
 simfile(),
 reader(),
 tracks(reader,"MCTrack"),
 rpcpoints(reader,"ShipRpcPoint"),
 strawtubespoints(reader,"strawtubesPoint"),
 targetpoints(reader,"TargetPoint")
{
 simfile = TFile::Open(inputfilename);
 reader.SetTree("cbmsim",simfile);    
 //importing geofile:
 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import(geofilename);
 //fucntion for mcs
 fOPERA_mcsres = new TF2("fOPERA_mcsres","([0]+[1]*x)/TMath::Sqrt(y)+ ([2]+[3]*x)+([4]+[5]*x)*TMath::Sqrt(y)");
 fOPERA_mcsres->SetParameters(0.397,0.019,0.176,0.042,-0.014,-0.003);

 hdeltaTX = new TH1D("hdeltaTX","TX difference;DeltaTX",400,-0.02,0.02);
 hdeltaTY = new TH1D("hdeltaTY","TY difference;DeltaTY",400,-0.02,0.02);
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

int EfficiencyCut::GeometricalEfficiency(double offsetxy, int Nminplates){

    double vx = tracks[0].GetStartX();
    double vy = tracks[1].GetStartY();
    double vz = tracks[2].GetStartZ();

    const int Ntotplates = 59;
    double startbrickx = 0.5;
    double startbricky = 0.5;

    int InteractionWall;
    int whichrow, whichcolumn; //for findbrick function
    int whichplate = 60;

    bool foundbrick = FindBrick(vx, vy, vz, InteractionWall, whichrow, whichcolumn, whichplate);
    if (!foundbrick){
        //cout<<"No brick found!, please check"<<endl;
        return -1;
    }
    TGeoVolume *volEmulsion = gGeoManager->GetVolume("Emulsion");
    TGeoBBox *Emulsion = (TGeoBBox*) volEmulsion->GetShape();
    double dxtarget = Emulsion->GetDX();
    double dytarget = Emulsion->GetDY();
    //double dxtarget = 20.05;
    //double dytarget = 20.05;
    //cout<<" Wall "<<InteractionWall<<" row "<< whichrow<<" column "<<whichcolumn<<" plate "<<whichplate<<endl;
    //cout<<" vx "<<vx<<" vy "<<vy<<endl;
    //longitudinal edge
    if (whichplate > (Ntotplates - Nminplates +1)) return 0;
    //upper transverse edges
    if (TMath::Abs(vx) > ((2*dxtarget) + startbrickx - offsetxy)) return 0;
    if (TMath::Abs(vy) > ((2*dytarget) + startbricky - offsetxy)) return 0;
    //lower transverse edges
    if (TMath::Abs(vx) < (startbrickx + offsetxy)) return 0;
    if (TMath::Abs(vy) < (startbricky + offsetxy)) return 0;
        
    else return 1;

}

bool EfficiencyCut::VisibleVertexLocation(double maxtantheta, double minmomentum){
     //********************************SECOND CONDITION: VERTEX LOCATION: at least one track with mom > 1 GeV and tantheta<1********/
    TDatabasePDG *pdg = TDatabasePDG::Instance();
    //getting particle charge and momentum
    int trackID = 0;
    for (auto &track:tracks){ //loop over MCTracks, find neutrino daughters
     double momentum = track.GetP();
     double charge = 0.;
     double pdgcode = track.GetPdgCode();

     double tx = track.GetPx()/track.GetPz();
     double ty = track.GetPy()/track.GetPz();
     double tantheta = TMath::Sqrt(tx * tx + ty * ty);
     if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge(); 
     //is it from the primary neutrino interaction?
     if(track.GetMotherId()==0 && TMath::Abs(charge)>0){ //charged neutrino daughter  
          if(tantheta<maxtantheta && momentum > minmomentum){ //clearly visible
            primaryvisible.push_back(trackID); //I do not add the tau lepton or charmed hadron since it decays too soon
            npl_nudaughterID[trackID] = 0;
            //cout<<ientry<<" "<<pdgcode<<" "<<momentum<<" "<<tantheta<<endl;
            } //end tantheta condition
      } //end charged daugter condition
      trackID++;
     }//end over track loop
    if (primaryvisible.size()>1) return true; //at least 1 or 2? Need one more than tau lepton
    else return false;
}

bool EfficiencyCut::DecaySearch(bool tausim, double maxdl = 0.4, double minkinkangle = 0.02, double minip = 10e-4, double mindaumomentum = 0.1, double maxdautantheta = 1.){
 //looking for decay: 
 //tausim true: we look for tau decay, 
 //tausim false:we look for charmed hadrons
 ROOT::RVec<int> signalpdgs;
 if (tausim) signalpdgs = {15}; //tau lepton
 else signalpdgs = {411, 431, 4122, 421, 4132, 4232, 4332, 441};

 int trackID = 0;
 TDatabasePDG *pdg = TDatabasePDG::Instance();
 //primary neutrino vertex coordinates
 double vx = tracks[0].GetStartX();
 double vy = tracks[0].GetStartY();
 double vz = tracks[0].GetStartZ();
 int tauid = -1; //assuming decaying particle to have trackID 1
 for (auto &track:tracks){ //loop over MCTracks, find decay daughters
    //this is the particle we are looking for (tau lepton or hadron)
    double pdgcode = track.GetPdgCode();
    bool passed_decaysearch = true;
    int signalparticle = signalpdgs[signalpdgs==TMath::Abs(pdgcode)].size();
    if ( signalparticle > 0 && track.GetMotherId()==0 && tauid < 0) tauid = trackID;
    
    if (tauid < 0){ 
        trackID++;
        continue; //MCTrackID for tau lepton not yet found, go to next
    }
    //tau lepton angles
    double tautx = tracks[tauid].GetPx()/tracks[tauid].GetPz();
    double tauty = tracks[tauid].GetPy()/tracks[tauid].GetPz();
    //condition for a visible tau decay daughter

    //look for charged particles from tau decay
    double momentum = track.GetP();
    double charge = 0.;
    if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge(); 
    //check if it is a tau decay daughter
    if(track.GetMotherId()==tauid && TMath::Abs(charge)>0){ 
        //check if the decay is "visible"
          daughters.push_back(trackID);
          double tx = track.GetPx()/track.GetPz();
          double ty = track.GetPy()/track.GetPz();
          double tantheta = TMath::Sqrt(tx * tx + ty * ty);
        
          double tauendx = track.GetStartX();
          double tauendy = track.GetStartY();
          double tauendz = track.GetStartZ();
     
          double taudecaylength = TMath::Sqrt(pow(tauendx - vx,2) + pow(tauendy-vy,2)+ pow(tauendz-vz,2)); 
          double kinkangle = TMath::ATan(TMath::Sqrt(pow(tx - tautx,2)+pow(ty - tauty,2)));

          //'''Impact parameter of track with respect to primary vertex (standard transverse definition)'''
 
          double dz = vz - tauendz;
          double ipx = tx * dz + tauendx - vx;
          double ipy = ty * dz + tauendy - vy;
          double ip = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));

          //*************************DECAY SEARCH**************************//
          //cout<<"TEST daughter ID"<<trackID<<" tauID "<<endl;
          //cout<<tauid<<" "<<taudecaylength<<" "<<ip<<" "<<kinkangle<<" "<<momentum<<" "<<tantheta<<endl;
          if (taudecaylength > maxdl) passed_decaysearch = false;
          if (ip < minip) passed_decaysearch = false;
          if (kinkangle < minkinkangle) passed_decaysearch = false;
          if (momentum < mindaumomentum || tantheta > maxdautantheta) passed_decaysearch = false;
          //if the code was not returned by all these checks, the decay daughter is accepted
          if (passed_decaysearch) visibledaughters.push_back(trackID);
         } //tau daughter check
 trackID++;
 }//end track loop
 if(visibledaughters.size()>0) return true;
 else return false;
}


int EfficiencyCut::MCSmeasurement(double maxres, int nplmin){
    //this function must be called after VisibleVertexLocation
    //looping over targetpoints
    if (targetpoints.GetSize() == 0) return -1;
    int vertexbrickID = targetpoints[0].GetDetectorID()/1e+4; //the first hit stored is the closest temporarily to nu interaction
    for (const TargetPoint& targetpoint: targetpoints){
        int detectorID = targetpoint.GetDetectorID();
        int trackID = targetpoint.GetTrackID();
        int brickID = detectorID/1e+4;         
        //selection to accept hit for counting
        if (detectorID%2 == 1.) continue; //only accept one side of emulsion film      
        if (brickID != vertexbrickID) continue; //hit does not came from the same vertex as neutrino
        if (npl_nudaughterID.find(trackID) == npl_nudaughterID.end()) continue;//hit is not from a neutrino daughter

        npl_nudaughterID[trackID]++; //all conditions are accepted, increasing counter
    }
    //checking targetpoints
    for ( const auto &myPair : npl_nudaughterID ) {
         int npl = myPair.second;
         int trackID = myPair.first;

         double momentum = tracks[trackID].GetP();

         if(npl<nplmin) continue;

         double res = fOPERA_mcsres->Eval(momentum, npl); //sigma(1/P)/(1/P)
         //cout<<"test "<<res<<" momentum "<<momentum<<" npl "<<npl<<endl;   
         if (res < maxres) measured_momentum_nudaughterID[trackID] = true;
         else measured_momentum_nudaughterID[trackID] = false;
   }
   return 0;
}

void EfficiencyCut::FillHistograms(TH1D *hnuP, TH2D *hq2_x, TH1D* hlP){
    double nuE = tracks[0].GetEnergy(); //neutrino
    double lE = tracks[1].GetEnergy(); //primary lepton
    double v = nuE - lE; //transferred energy
    double mT = 0.94; //approximate mass proton/neutron in GeV

    TVector3 nuPvec(tracks[0].GetPx(), tracks[0].GetPy(), tracks[0].GetPz());
    TVector3 lPvec(tracks[1].GetPx(), tracks[1].GetPy(), tracks[1].GetPz());
    double theta = lPvec.Angle(nuPvec); //theta angle of scattered lepton vs neutrino

    //computing our variables   
    double Q2 = 4 * nuE * lE * pow(TMath::Sin(theta/2.),2); //do not forget theta/2.
    double x = Q2/(2* mT * v); 

    //cout<<"Neutrino energy "<<nuE<<" lepton Energy "<<lE<<" angle "<<theta<<endl; 
    //cout<<"x is "<<x<<" and Q2 is "<<Q2<<endl;
    //filling histograms
    hq2_x->Fill(TMath::Log10(x),TMath::Log10(Q2),this->GetEventWeight());
    int whichbin = hq2_x->FindBin(TMath::Log10(x), TMath::Log10(Q2));
    //RVec to compute average energy at each bin
    if (energies_q2x.find(whichbin) == energies_q2x.end()) energies_q2x[whichbin] = ROOT::RVecD();
    energies_q2x[whichbin].push_back(nuE);

    //1D histogram for quick check
    hnuP->Fill(tracks[0].GetP(),this->GetEventWeight() ); //neutrino momentum
    hlP->Fill(tracks[1].GetP(),this->GetEventWeight() ); //primary lepton momentum
}

//function to compute curvature from sagitta spectrometer
ROOT::RVec<Double_t> SagittaSpectrometer(int* whichplanes, Double_t *XSpectrometer, Double_t *YSpectrometer, Double_t *ZSpectrometer){
  //0 and 4: upstream stations;
  Double_t MeanXout = (XSpectrometer[whichplanes[0]] + XSpectrometer[whichplanes[3]])/2.;
  Double_t MeanYout = (YSpectrometer[whichplanes[0]] + YSpectrometer[whichplanes[3]])/2.;

  Double_t MeanXin = (XSpectrometer[whichplanes[1]] + XSpectrometer[whichplanes[2]])/2.;
  Double_t MeanYin = (YSpectrometer[whichplanes[1]] + YSpectrometer[whichplanes[2]])/2.;

  ROOT::RVec<Double_t> Sagitta = {-99999.,-99999.};

  Sagitta[0] = MeanXin - MeanXout;
  Sagitta[1] = MeanYin - MeanYout;

  return Sagitta;
}

//in this configuration, we have two tracking stations upstream of the magnet, two downstream. Bending Angle measurement
ROOT::RVec<Double_t> BendingSpectrometer(int* whichplanes, Double_t *XSpectrometer, Double_t *YSpectrometer, Double_t * ZSpectrometer){
    //0 and 1: upstream stations;
    Double_t TXup = (XSpectrometer[whichplanes[1]] - XSpectrometer[whichplanes[0]])/(ZSpectrometer[whichplanes[1]] - ZSpectrometer[whichplanes[0]]);
    Double_t TYup = (YSpectrometer[whichplanes[1]] - YSpectrometer[whichplanes[0]])/(ZSpectrometer[whichplanes[1]] - ZSpectrometer[whichplanes[0]]);
    //2 and 3: downstream stations;
    Double_t TXdown = (XSpectrometer[whichplanes[3]] - XSpectrometer[whichplanes[2]])/(ZSpectrometer[whichplanes[3]] - ZSpectrometer[whichplanes[2]]);
    Double_t TYdown = (YSpectrometer[whichplanes[3]] - YSpectrometer[whichplanes[2]])/(ZSpectrometer[whichplanes[3]] - ZSpectrometer[whichplanes[2]]);
  
    ROOT::RVec<Double_t> DeltaT = {-99999.,-99999.};

    DeltaT[0] = TXdown - TXup;
    DeltaT[1] = TYdown - TYup;

    return DeltaT;
}

bool EfficiencyCut::SpectrometerAcceptance(double posres, double sagittares){
 //first, looping over all rpcpoints to get the lepton for ID fmuonID
 //const int nstations = 8;
 const int nstations = 21;
 Double_t XSpectrometer[nstations], YSpectrometer[nstations], ZSpectrometer[nstations];
 Bool_t InSpectrometer[nstations];
 //initialize them to false
 for (int istation = 0; istation < nstations; istation++){
  InSpectrometer[istation] = false;
 }

 for (const ShipRpcPoint& rpcpoint: rpcpoints){
  //int nstation = rpcpoint.GetDetectorID() - 1; //from 0 to 3 
  int nstation = rpcpoint.GetDetectorID();
  if (rpcpoint.GetTrackID()==fmuonID && !InSpectrometer[nstation]){ //we fill the arrays with the muon neutrino positions
  
   XSpectrometer[nstation] = rpcpoint.GetX();
   YSpectrometer[nstation] = rpcpoint.GetY();
   ZSpectrometer[nstation] = rpcpoint.GetZ();
   InSpectrometer[nstation] = true;

   //applying smearing to X and Y
   XSpectrometer[nstation] = XSpectrometer[nstation] + gRandom->Gaus(0, posres);
   YSpectrometer[nstation] = YSpectrometer[nstation] + gRandom->Gaus(0, posres);
  }//end muon track condition
 }//end rpc hit loop 

 //if the first 4 were found, checking resolution
 if (InSpectrometer[0]&&InSpectrometer[9]&&InSpectrometer[10]&&InSpectrometer[19]){
    int whichplanes[4] = {0,9,10,19}; //which planes to consider for the sagitta computation
    auto Sagitta = SagittaSpectrometer(whichplanes, XSpectrometer, YSpectrometer, ZSpectrometer);
    if (TMath::Abs(Sagitta[1]) > 3 * sagittares) return true;
    else{ 
        cout<<"In acceptance, but sagitta not enough "<<Sagitta[1]<<endl;
        return false;
    }
 }
 else if (InSpectrometer[0]&&InSpectrometer[1]&&InSpectrometer[2]&&InSpectrometer[3]){ //not arrives downstream, using upstream planes
    int whichplanes[4] = {0,1,2,3};
    auto Sagitta = SagittaSpectrometer(whichplanes, XSpectrometer, YSpectrometer, ZSpectrometer);
    if (TMath::Abs(Sagitta[1]) > 3 * sagittares) return true;
    else{ 
        cout<<"In acceptance, but sagitta not enough "<<Sagitta[1]<<endl;
        return false;
    }
 } 
 else return false;
}
//same, but using strawtubes of Decay Spectrometer this time
bool EfficiencyCut::DecaySpectrometerAcceptance(double dsposres, double dsangularres){
 //first, looping over all rpcpoints to get the lepton for ID fmuonID
 const int nstations = 4; //four stations only this time
 Double_t XSpectrometer[nstations], YSpectrometer[nstations], ZSpectrometer[nstations];
 Bool_t InSpectrometer[nstations];
 //initialize them to false
 for (int istation = 0; istation < nstations; istation++){
  InSpectrometer[istation] = false;
 }

 for (const strawtubesPoint& strawtubespoint: strawtubespoints){
  int nstation = (strawtubespoint.GetDetectorID()/1e+7)-1; //from 0 to 3 
  if (strawtubespoint.GetTrackID()==fmuonID && !InSpectrometer[nstation]){ //we fill the arrays with the muon neutrino positions
  
   XSpectrometer[nstation] = strawtubespoint.GetX();
   YSpectrometer[nstation] = strawtubespoint.GetY();
   ZSpectrometer[nstation] = strawtubespoint.GetZ();
   InSpectrometer[nstation] = true;

   //applying smearing to X and Y
   XSpectrometer[nstation] = XSpectrometer[nstation] + gRandom->Gaus(0, dsposres);
   YSpectrometer[nstation] = YSpectrometer[nstation] + gRandom->Gaus(0, dsposres);
  }//end muon track condition
 }//end strawtubes hit loop 
 if (InSpectrometer[0]){ //TEST ONLY FIRST STATION
 //if (InSpectrometer[0]&&InSpectrometer[1]&&InSpectrometer[2]&&InSpectrometer[3]){ //not arrives downstream, using upstream planes
    int whichplanes[4] = {0,1,2,3};
    auto DeltaT = BendingSpectrometer(whichplanes, XSpectrometer, YSpectrometer, ZSpectrometer);
    //studying the distributions
    hdeltaTX->Fill(DeltaT[0]);
    hdeltaTY->Fill(DeltaT[1]);
    //Filling histograms for study
   /* if (TMath::Abs(DeltaT[1]) > 3 * dsangularresres) return true;
    else{ 
        cout<<"In acceptance, but decay spectro sagitta not enough "<<DeltaT[1]<<endl;
        return false;
    }*/
    return true; //if present, always return true for now, regardless of angular deflection
 } 
 else return false;
}

int EfficiencyCut::DecayChannel(){
    //tau decay channel identification: 
    //1 -> mu
    //2 -> e
    //3 -> 1h
    //4 -> 3h
    //0 -> failed to identify
    int nmuons = 0;
    int nhadrons = 0;
    int nelectrons = 0; 

    for (auto trackID:daughters){
       double pdgcode = tracks[trackID].GetPdgCode();
       if (TMath::Abs(pdgcode) == 11) nelectrons++;
       else if (TMath::Abs(pdgcode) == 13){ 
        fmuonID = trackID; //storin to follow it later
        nmuons++;
       }
       else nhadrons++;
     }
     if (nmuons == 1) return 1;
     else if (nelectrons == 1) return 2;
     else if (nhadrons == 1) return 3;
     else if (nhadrons > 1) return 4;
     /*else{
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
     }*/
     return 0;

}