//code for studing tau neutrino simulation and evaluating efficiencies for my PhD thesis (created on 22 June 2020 by A.Iuliano)
//usage: set tausim = true and uncomment nutau and anutau CCDIS files for nutau signal sim
//.....  set tausim = false and uncomment numu CharmCCDIS file for bkg sim
//defined functions
TVector3 NeutrinoVertexCoordinates(const ShipMCTrack &track);
Bool_t FindBrick(Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate);

int DecayChannel(vector<int> &daughters,TTreeReaderArray<ShipMCTrack> &tracks, int nnue, int nnumu, int ievent);

void DecodeBrickID(Int_t detID, Int_t &NWall, Int_t &NRow, Int_t &NColumn, Int_t &NPlate, Bool_t &EmCES, Bool_t &EmBrick, Bool_t &EmTop);
bool NeutrinoVertexLocation(int trackID, const ShipMCTrack& track, vector<int> &primary, vector<int> &primaryvisible, ROOT::RVec<int> signalpdgs);
bool GeometricalEfficiency(TVector3 Vn, double offsetxy, int Nminplates, int &InteractionWall);
bool TauDecay(int trackID, int tauid, const ShipMCTrack &track, TVector3 Vn, double tautx, double tauty, double& taudecaylength);

void smearing (double &Xpos, double &Ypos, double spaceres);
vector<double> eff_formula(int foundweight, int totalweight, int Nevents_total);

double phiangle(TVector3 p1, TVector3 p2){
  //angle between the two vectors in the transverse plane (with respect to the beam direction, i.e. z)
  double transverse_dotproduct = p1.X() * p2.X() + p1.Y() * p2.Y();
  double phi = TMath::ACos(transverse_dotproduct/(p1.Pt()*p2.Pt()));
  return phi;
}
//GLOBAL VARIABLES AND HISTOGRAMS
 //*********DEFINITION OF HISTOGRAMS*********************//
 //tau kinematics
 TH1I *hnfilmtau = new TH1I("hnfilmtau","Number of emulsion films passed by a tau lepton;Nfilms",15,0,15);
 TH2D *htauppt = new TH2D("htauppt","Transverse momentum vs momentum of tau lepton;P[GeV/c];Pt[GeV/c]",400,0,400,100,0,10);
 TH1D *htaugamma = new TH1D("htaugamma","Gamma of tau lepton;#gamma",100,0,100);
 //primary vertex position
 TH2D *hvxy = new TH2D("hvxy","Transverse position of vertices",80,-40,40,80,-40,40);
 TH1D *hvz =  new TH1D("hvz","Transverse position of vertices",400,-3400,-3000);
 //decay parameters
 TH1I *hdaughters = new TH1I("hdaughters","Number of daughters per decay", 20, 0,20);
 TH1I *hvisibledaughters = new TH1I("hvisibledaughters","Number of visible daughters per decay", 20, 0,20);

 TH1D *hdl = new TH1D("hdl","Tau Decay length;dl[mm]",300,0,30);
 TH1D *hkink = new TH1D("hkink","Kink angle;#Theta[rad]",50,0,0.5);
 TH1D *hip = new TH1D("hip","Impact Parameter;IP[#mum]",100,0,1000);

 TH1I *hchannel = new TH1I("hchannel","Tau lepton decay channel;IChannel",4,1,5);
 //muon kinematics
 TH1D *hmuonpall = new TH1D("hmuonpall","Momentum of all muons from tau lepton decay;P[GeV/c]",400,0,400);
 TH1D *hmuonangleall = new TH1D("hmuonangleall","Angle of all muons;#Theta[rad]",100,0,1);
 TH1D *hmuonpentered = new TH1D("hmuonpentered","Momentum of muons entered in first RPC;P[GeV/c]",400,0,400);
 TH1D *hmuonangleentered = new TH1D("hmuonangleentered","Angle muons entered in first RPC;#Theta[rad]",100,0,1);
 //electron shower 
 TH1D *helectron_E = new TH1D("helectron_E","Energy of electrons from tau decay;E[GeV]",100,0,100);
 TH1D *helectron_NPlate = new TH1D("helectron_NPlate","Number of plate where the electron is produced;Nplate",57,1,58);
 
 //DT sagitta

 TH1D *hxsagitta_DT = new TH1D("hxsagitta_DT","Sagitta in the xzplane;sxz[cm]",20,-0.1,0.1);
 TH1D *hysagitta_DT_positive = new TH1D("hysagitta_DT_positive","Sagitta for positive particles;syz[cm]",30,-15,15);
 TH1D *hysagitta_DT_negative = new TH1D("hysagitta_DT_negative","Sagitta for negative particles;syz[cm]",30,-15,15);

 TH1D *hxsagitta_CES = new TH1D("hxsagitta_CES","Sagitta in the xzplane;sxz[cm]",20,-0.001,0.001);
 TH1D *hysagitta_CES_positive = new TH1D("hysagitta_CES_positive","Sagitta for positive particles;syz[cm]",20,-0.01,0.01);
 TH1D *hysagitta_CES_negative = new TH1D("hysagitta_CES_negative","Sagitta for negative particles;syz[cm]",20,-0.01,0.01);

 //muon occupancy
 TH2D *hoccupancy_rpc_clusters = new TH2D("hoccupancy_rpc_clusters","How many per station?;istation;occupancy",8,0,8,10,0,10);
 
 //phi angle between primary lepton and hadronic system
 TH1D *hphi = new TH1D("hphi","Phi angle between primary lepton and hadronic system;#phi[rad]",22,-1,3.4);

 //**VARIABLES***//
 //maximum dx an dy of vertex to be accepted for efficiency studies
 const double dxtarget = 40.5;
 const double dytarget = 40.5;

 TDatabasePDG *pdg = TDatabasePDG::Instance();

//start main script
int nutau_event(){
 bool tausim = true;
 bool doelectron = false; //fills histograms for electron neutrino interactions
 bool dosagitta = true; //do loop over DTs and sagitta computation
 ROOT::RVec<int> signalpdgs; //particles I want to study the decay
 if (tausim) signalpdgs = {15}; //tau lepton
 else signalpdgs = {411, 431, 4122, 421, 4132, 4232, 4332, 441}; //adding 4122 makes program crash
 //getting tree and defining arrays
 TChain treechain("cbmsim"); //I add a sim of tau and antitaus
 treechain.Add("/home/utente/Simulations/tauneutrino_19June2020/ship.conical.Genie-TGeant4.root"); 
 treechain.Add("/home/utente/Simulations/tauantineutrino_22September2020/ship.conical.Genie-TGeant4.root");
 //treechain.Add("/home/utente/Simulations/muneutrino_charm_05September2020/ship.conical.Genie-TGeant4.root"); 
 const double phimin = 2.2; //minimum phi to accept selection
 //if (!file) return;
 TTreeReader reader(&treechain);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");
 TTreeReaderArray<HptPoint> dtpoints(reader,"HptPoint");
 TTreeReaderArray<ShipRpcPoint> rpcpoints(reader,"ShipRpcPoint");
 
 int trackID, detID, pdgcode;
 int ntauhits;
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
 double largephiweight[ndecaychannels] = {0.,0.,0.,0.};
 double chargedetweight[ndecaychannels] = {0.,0.,0.,0.};

 double muonacceptance = 0., muonefficiency = 0.;
 int InteractionWall = 0;

 bool isgeometrical, islocated, decaysearch, chargeeff;
 bool largephi;

 //IDs of tracks of interest (Tracks at neutrino vertex other than tau and tau daughters)
 vector<int> primary;
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
 const int nCES = 3;
 const int nstations = 5;
 const int maxndaughters = 10;
 double XCES_daughter[maxndaughters][nCES],YCES_daughter[maxndaughters][nCES],ZCES_daughter[maxndaughters][nCES]; 
 double XDT_daughter[maxndaughters][nstations],YDT_daughter[maxndaughters][nstations],ZDT_daughter[maxndaughters][nstations]; 
 const double DTspaceres = 100 * 1e-4; //100 micron in cm
 const double CESspaceres = 1.5 * 1e-4; //5 micron in cm
 //muon rpc check
 const int nrpcstations = 8;
 double XRPC_muon[nrpcstations], YRPC_muon[nrpcstations], ZRPC_muon[nrpcstations];
 int occupancy[nrpcstations];
 bool muoninfirstrpc, muoninrpc;
 const double clusterradius = 2; //maximum cluster radius
 int nisolated_rpc_clusters;
 const double minisolated_rpc_clusters = 4; //how many hits required to be isolated

 //***********************************START OF MAIN LOOP*************************//
 for(int ientry = 0;ientry<nentries;ientry++){    
     //resetting counters and booleans
     isgeometrical = false;
     islocated = false;
     decaysearch = false;
     largephi = false;
     chargeeff = false;

     nnue = 0;
     nnumu = 0;

     InteractionWall = -1;
     //ntauhits = 0;
     //clearing list of tracks of interest
     primaryvisible.clear();
     primary.clear();
     daughters.clear();
     visibledaughters.clear();

     taudecaylength = -1000.;

     //muon detection variables
     muoninfirstrpc = false; //at start of rpc
     muoninrpc = false; //anywhere in rpc
     nisolated_rpc_clusters = 0;
     //getting entry

     if (ientry%10000 == 0) cout<<"arrived at entry" <<ientry<<endl;
     reader.SetEntry(ientry);// keeps track of the number of event (from 0 to Nevents - 1)
     //event weight
     double eventweight = tracks[0].GetWeight();
     //neutrino interaction coordinates
     TVector3 Vn = NeutrinoVertexCoordinates(tracks[0]);
     //if neutrino interacted off of our target, go to next
     if ( (TMath::Abs(Vn(0))>dxtarget) || (TMath::Abs(Vn(1))>dytarget) ) continue;

     hvxy->Fill(Vn(0),Vn(1), eventweight);
     hvz->Fill(Vn(2), eventweight);

     //tau lepton variables


     //********************************FIRST CONDITION: GEOMETRICAL SELECTION********/
     isgeometrical = GeometricalEfficiency(Vn,trans_mindist, Nminplates, InteractionWall);

     //******************access the array of tracks
     itrack = 0;
     tauid = -1;
     //*********************************START TRACK LOOP****************************//
     for (const ShipMCTrack& track: tracks){        
         pdgcode = track.GetPdgCode();
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
         //look for charged particles from primary vertex
         NeutrinoVertexLocation(itrack, track, primary, primaryvisible, signalpdgs);
         //look for charged tracks from tau decay lengths

         if(track.GetMotherId()==tauid && TMath::Abs(charge)==0){ //counting neutrinos
             if (TMath::Abs(pdgcode) == 12) nnue++;
             else if (TMath::Abs(pdgcode) == 14) nnumu++;
         } 
         if(track.GetMotherId()==tauid && TMath::Abs(charge)>0){ //daughter of tau/charm particle, adding to daughter list
          daughters.push_back(itrack);
          //check if decay is visible
          if (TauDecay(itrack, tauid, track, Vn, tautx, tauty, taudecaylength)) visibledaughters.push_back(itrack);            
         }

         itrack++;        
     } //end of track loop
     hdaughters->Fill(daughters.size());
     hvisibledaughters->Fill(visibledaughters.size());
     //decay channel identification
     if (tausim) whichchannel = DecayChannel(daughters,tracks,nnumu,nnue,ientry);
     else whichchannel = 1;
     hchannel->Fill(whichchannel);
     hdl->Fill(taudecaylength*10);
     
     //a visible electron from tau decay
     if (whichchannel == 2 && isgeometrical && visibledaughters.size() == 1 && doelectron){

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

     if (islocated && decaysearch){ //it has sense only for good DS events
       //loop over primary tracks (hadrons only), total momentum and phi
       //Nota Bene: primaryivisible contains all charged tracks with mumID 0 (tantheta < 1, P > 1 GeV), except tau/charm
       TVector3 p3tothad(0., 0., 0.);
       TVector3 p3tau(tracks[tauid].GetPx(), tracks[tauid].GetPy(), tracks[tauid].GetPz());
       double phimax = 0.;
       TVector3 p3maxhad(0., 0., 0.); //I need to find the hadron with maximum phi (to remove it later)
       for (int &primaryID:primary){
         int primarypdgcode = tracks[primaryID].GetPdgCode();      
         double pxhad = tracks[primaryID].GetPx();
         double pyhad = tracks[primaryID].GetPy();
         double pzhad = tracks[primaryID].GetPz();

         TVector3 p3had(pxhad, pyhad, pzhad);
         p3tothad = p3tothad + p3had;
         
         double phi = phiangle(p3had, p3tau);
         if (phi > phimax){ //if it has phi with tau/charm larger than current maximum, I overwrite it
           phimax = phi;
           p3maxhad = p3had;
         }
        }
      //removing maximum from total momentum and computing phi angle with tau/charm
       p3tothad = p3tothad - p3maxhad;
       double phitau; //nonsense value larger than pi, for safety
       if (primary.size()==1) phitau = -0.9;
       else phitau = phiangle(p3tothad, p3tau);
       hphi->Fill(phitau);
       if (phitau > phimin) largephi = true;
     }//end deltaphi section
     //access the hits: 
     //somming value, when efficiency is satisfied
     //****************************look for tau decay daughters in downstream trackers***//
     if (visibledaughters.size() > maxndaughters){ 
       cout<<"ERROR: number of daughters larger than maximum,exiting"<<endl;
       return 1;
     }
     if (dosagitta){
      //resetting indeces to be filled with daughters positions in DT stations
      for (int idaughter = 0; idaughter < visibledaughters.size(); idaughter++){
       for (int index = 0; index < nstations; index++){
        for (int index = 0; index < nCES; index++){
        XCES_daughter[idaughter][index] = 100000.;
        YCES_daughter[idaughter][index] = 100000.;
        ZCES_daughter[idaughter][index] = 100000.;
       }
        XDT_daughter[idaughter][index] = 100000.;
        YDT_daughter[idaughter][index] = 100000.;
        ZDT_daughter[idaughter][index] = 100000.;
       }
      }
      //loop over emulsion point to look for CES daughters
      for (const TargetPoint& targetpoint: targetpoints){
        trackID = targetpoint.GetTrackID(); 
        detID = targetpoint.GetDetectorID();
        if (trackID == 1){ //hit from tau lepton
            ntauhits++;
        }
        for (int idaughter = 0; idaughter < visibledaughters.size(); idaughter++){
          if (trackID == visibledaughters[idaughter]){//hit of a daughter
           //resetting counters
           int NWall = 0;
           int NRow = 0;
           int NColumn = 0;
           int NPlate = 0;
           bool EmCES = false;
           bool EmBrick = false;
           bool EmTop = false;
           DecodeBrickID(detID, NWall, NRow, NColumn, NPlate, EmCES, EmBrick, EmTop);
           if (EmCES && ((InteractionWall + 1) == NWall)){ //is the hit from CES detector, is it from the same wall of the neutrino interaction?
             XCES_daughter[idaughter][NPlate -1] = targetpoint.GetX();
             YCES_daughter[idaughter][NPlate -1] = targetpoint.GetY();
             ZCES_daughter[idaughter][NPlate -1] = targetpoint.GetZ();    
             smearing(XCES_daughter[idaughter][NPlate-1], YCES_daughter[idaughter][NPlate-1], CESspaceres);      
            }            
          } //end check for daughter
        } //end loop over daughters
      } //end of hit loop
      hnfilmtau->Fill(ntauhits,eventweight);
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
       double sigmasagitta = 0.0172;
       double particlecharge = pdg->GetParticle(tracks[visibledaughters[idaughter]].GetPdgCode())->Charge();
       if (ZDT_daughter[idaughter][2] < 100000. && ZDT_daughter[idaughter][1] < 100000. && ZDT_daughter[idaughter][0] < 100000.){
        double ysagitta = ((YDT_daughter[idaughter][0] + YDT_daughter[idaughter][2])/2.) - YDT_daughter[idaughter][1];
        double xsagitta = ((XDT_daughter[idaughter][0] + XDT_daughter[idaughter][2])/2.) - XDT_daughter[idaughter][1];
        hxsagitta_DT->Fill(xsagitta);
        if (particlecharge > 0.) hysagitta_DT_positive->Fill(ysagitta);
        else hysagitta_DT_negative->Fill(ysagitta);
        //I currently accept the decay as kinematically reconstructed if at least one is seen
        if (TMath::Abs(ysagitta) > (3*sigmasagitta)) chargeeff = true;      
       } //end condtion of particle passing through first three DTs

       if (ZCES_daughter[idaughter][2] < 100000. && ZCES_daughter[idaughter][1] < 100000. && ZCES_daughter[idaughter][0] < 100000.){
        double sigmasagitta_CES = 1.996e-4;

        double ysagitta = ((YCES_daughter[idaughter][0] + YCES_daughter[idaughter][2])/2.) - YCES_daughter[idaughter][1];
        double xsagitta = ((XCES_daughter[idaughter][0] + XCES_daughter[idaughter][2])/2.) - XCES_daughter[idaughter][1];
        hxsagitta_CES->Fill(xsagitta);
        if (particlecharge > 0.) hysagitta_CES_positive->Fill(ysagitta);
        else hysagitta_CES_negative->Fill(ysagitta);
        if (TMath::Abs(ysagitta) > (3*sigmasagitta_CES)) chargeeff = true;
        //I currently accept the decay as kinematically reconstructed if at least one is seen
        //if (TMath::Abs(ysagitta) > (3*sigmasagitta)) chargeeff = true;      
       } //end condtion of particle passing through first three CES
      } //end loop over visible daughters
     } //end condition if doing or not charge detection efficiency estimation
     if (whichchannel == 1){ //only for muon channel,
      int muonid;
      if (tausim) muonid = daughters[0]; //tau 1mu case, first and only charged daughter
      else muonid = 1; //numu charm case, associated lepton
      //***************************look for muons in downstream detectors****************//
      //resetting containers
      for (int index = 0; index < nrpcstations; index++){
        XRPC_muon[index] = 100000.;
        YRPC_muon[index] = 100000.;
        ZRPC_muon[index] = 100000.;
        occupancy[index] = 0;
      }
      //first loop, storing clusters
      for (const ShipRpcPoint& rpcpoint: rpcpoints){
        trackID = rpcpoint.GetTrackID(); 
        detID = rpcpoint.GetDetectorID();        
       
        if (trackID == muonid){ //hit from tau lepton
         int whichstation = detID - 10000; //0, 1, 2, 3, 4, 5, 6, 7
         if (whichstation >= 0){ //found cluster center with muon, storing position
            XRPC_muon[whichstation] = rpcpoint.GetX();
            YRPC_muon[whichstation] = rpcpoint.GetY();
            ZRPC_muon[whichstation] = rpcpoint.GetZ();
            occupancy[whichstation]++;
         }
         muoninrpc = true;
         if (detID == 10000){
            muoninfirstrpc = true;
            hmuonpentered->Fill(tracks[trackID].GetP());
            double muonangle = 
             TMath::ATan(TMath::Sqrt(pow(tracks[trackID].GetPx()/tracks[trackID].GetPz(),2)+pow(tracks[trackID].GetPy()/tracks[trackID].GetPz(),2)));
            hmuonangleentered->Fill(muonangle);
          }
        }//end condition for muonid

      } //end first rpc hit loop
      //starting second loop, to look for isolated clusters
      if (muoninrpc){
       for (const ShipRpcPoint& rpcpoint: rpcpoints){
        detID = rpcpoint.GetDetectorID();
        int whichstation = detID - 10000;   
        pdgcode = rpcpoint.PdgCode();
        double charge = 0.;
        if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
        //if charged, checking for position
        if (TMath::Abs(charge) > 0. && occupancy[whichstation]>0){
         double rpcx = rpcpoint.GetX();
         double rpcy = rpcpoint.GetY();
         double radius = TMath::Sqrt(pow(rpcx - XRPC_muon[whichstation],2) + pow(rpcy - YRPC_muon[whichstation],2));
         if(radius > clusterradius) occupancy[whichstation]++;
        }

       } //end second rpc hit loop
      }//end if muoninrpc condition
      //start loop over stations, checking how many full and outside muon cluster
      for (int istation = 0; istation < nrpcstations; istation++){
        hoccupancy_rpc_clusters->Fill(istation,occupancy[istation]);
        if (occupancy[istation]==1) nisolated_rpc_clusters++;
      }
      
      if (isgeometrical&&islocated&&decaysearch&&muoninfirstrpc) muonacceptance += eventweight; //note, we do not require the muon to transverse the DTs (chargeeff)      
     } //end of muon channel check

    if (whichchannel > 0){
      totalweight[whichchannel-1] += eventweight;
      if (isgeometrical){
        geometricalweight[whichchannel-1] += eventweight;
        if (islocated){
          localizedweight[whichchannel-1] += eventweight;
          if (decaysearch){ 
           decaysearchweight[whichchannel-1] += eventweight;
           if(largephi){
            largephiweight[whichchannel-1]+=eventweight;
            if (chargeeff){ 
              chargedetweight[whichchannel - 1 ] += eventweight;
              if (whichchannel == 1 && nisolated_rpc_clusters >= minisolated_rpc_clusters){
               if (tausim) muonefficiency += eventweight;
              } //muon efficiency (does require charge eff)
              else if(!tausim) muonefficiency += eventweight; 
            }//charge eff
          } //decay search
         }//large phi
        } //location
      }//geometrical

    }

 } // end of event loop
 //results
 int Nevents_total = hvz->GetEntries();
 cout<<endl;
 cout<<"****************PRINTING EFFICIENCY ESTIMATIONS*********************"<<endl;
 cout<<"Analying a number of interactions: "<<Nevents_total<<" total weight "<<hvz->Integral()<<endl;
 cout<<"channel legenda: (1 = mu, 2 = e, 3 = 1h, 4 = 3h)"<<endl;
 cout<<"Note: charge and muon efficiency are independently computed, but both require decay search to be true "<<endl;
 double alltotalweight = 0.; //not dependant from channel, initially set to 0
 double allgeometricalweight = 0.;
 double alllocalizedweight = 0.;
 for (int ichannel = 0; ichannel < ndecaychannels; ichannel++){
  cout<<endl;
  alltotalweight+= totalweight[ichannel];
  allgeometricalweight+=geometricalweight[ichannel];
  alllocalizedweight+= localizedweight[ichannel];
  int Nevents_ichannel = hchannel->GetBinContent(ichannel+1);
  //computing efficiencies with their values
  vector<double> vector_geomeff = eff_formula(geometricalweight[ichannel],totalweight[ichannel],Nevents_ichannel);
  vector<double> vector_localizedeff = eff_formula(localizedweight[ichannel],totalweight[ichannel],Nevents_ichannel);
  vector<double> vector_decaysearcheff = eff_formula(decaysearchweight[ichannel],totalweight[ichannel],Nevents_ichannel);
  vector<double> vector_largephieff = eff_formula(largephiweight[ichannel],totalweight[ichannel],Nevents_ichannel);
  vector<double> vector_chargedeteff = eff_formula(chargedetweight[ichannel],totalweight[ichannel],Nevents_ichannel);
  //reporting values and errors
  cout<<"Number of interactions in channel "<<ichannel+1<<" is "<<Nevents_ichannel<<" Total weigth: "<<totalweight[ichannel]<<endl;
  cout<<"Fraction within fiducial volume: "<<vector_geomeff[0]<<" with error "<<vector_geomeff[1]<<endl;
  cout<<"Fraction of localized vertices: "<<vector_localizedeff[0]<<" with error "<<vector_localizedeff[1]<<endl;
  cout<<"Fraction of decay search: "<<vector_decaysearcheff[0]<<" with error "<<vector_decaysearcheff[1]<<endl;
  cout<<"Fraction of large phi events: "<<vector_largephieff[0]<<" with error "<<vector_largephieff[1]<<endl;
  cout<<"Fraction of charge detected: "<<vector_chargedeteff[0]<<" with error "<<vector_chargedeteff[1]<<endl;
  if (ichannel == 0){ //only for tau muon decay channel or charm simulation
   vector<double> vector_muonacceptance = eff_formula(muonacceptance, totalweight[ichannel],Nevents_ichannel);
   vector<double> vector_muoneff = eff_formula(muonefficiency, totalweight[ichannel],Nevents_ichannel);
   cout<<"Fractions of muons from tau decays in first rpc: "<<vector_muonacceptance[0]<<" with error "<<vector_muonacceptance[1]<<endl;
   if (tausim) cout<<"Fractions of muons from tau decays with at least "<<minisolated_rpc_clusters<<" isolated rpc clusters: "<<vector_muoneff[0]<<" with error "<<vector_muoneff[1]<<endl;
   else cout<<"Fractions of events WITHOUT muons from tau decays with at least "<<minisolated_rpc_clusters<<" isolated rpc clusters: "<<vector_muoneff[0]<<" with error "<<vector_muoneff[1]<<endl;
  }
 }
 vector<double> vector_allgeomeff = eff_formula(allgeometricalweight,alltotalweight,Nevents_total);
 vector<double> vector_alllocalizedeff = eff_formula(alllocalizedweight,alltotalweight,Nevents_total);
 //only for ichannel == 0 (muon)
 cout<<"Fraction within fiducial volume over all: "<<vector_allgeomeff[0]<<" with error "<<vector_allgeomeff[1]<<endl;
 cout<<"Fraction of localized vertices over all: "<<vector_alllocalizedeff[0]<<" with error "<<vector_alllocalizedeff[1]<<endl;
 cout<<"*******ADVICE: IF YOU THINK YOU HAVE CHECKED ENOUGH THIS VALUES...CHECK THEM ONE MORE TIME********"<<endl;
 cout<<endl;
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
 
 TCanvas *csagitta_DT = new TCanvas();
 csagitta_DT->Divide(2,1);
 csagitta_DT->cd(1);
 hxsagitta_DT->Draw();
 hxsagitta_DT->Fit("gaus");
 csagitta_DT->cd(2);
 hysagitta_DT_negative->Draw();
 hysagitta_DT_positive->SetLineColor(kRed);
 hysagitta_DT_positive->Draw("SAMES");
 gPad->BuildLegend();

 TCanvas *csagitta_CES = new TCanvas();
 csagitta_CES->Divide(2,1);
 csagitta_CES->cd(1);
 hxsagitta_CES->Draw();
 hxsagitta_CES->Fit("gaus");
 csagitta_CES->cd(2);
 hysagitta_CES_negative->Draw();
 hysagitta_CES_positive->SetLineColor(kRed);
 hysagitta_CES_positive->Draw("SAMES");
 gPad->BuildLegend();

 TCanvas *coccupancy = new TCanvas();
 hoccupancy_rpc_clusters->Draw("COLZ");

 TCanvas *cdeltaphi = new TCanvas();
 hphi->Draw();

 return 0;
}


TVector3 NeutrinoVertexCoordinates(const ShipMCTrack &track){
     //building a tvector3 with neutrino vertex coordinates
     double vx = track.GetStartX();
     double vy = track.GetStartY();
     double vz = track.GetStartZ();
    
     return TVector3(vx,vy,vz);
}

//decode brick id function by Annarita's Target class
void DecodeBrickID(Int_t detID, Int_t &NWall, Int_t &NRow, Int_t &NColumn, Int_t &NPlate, Bool_t &EmCES, Bool_t &EmBrick, Bool_t &EmTop)
{
  Bool_t BrickorCES = 0, TopBot = 0;
  
  NWall = detID/1E7;
  NRow = (detID - NWall*1E7)/1E6;
  NColumn = (detID - NWall*1E7 -NRow*1E6)/1E4;
  Double_t b = (detID - NWall*1E7 -NRow*1E6 - NColumn*1E4)/1.E3;
  if(b < 1)
    {
      BrickorCES = 0;
      NPlate = (detID - NWall*1E7 -NRow*1E6 - NColumn*1E4 - BrickorCES*1E3)/1E1;
//      NPlate = detID - NWall*1E7 -NRow*1E6 - NColumn*1E4 - BrickorCES*1E3;
    }
  if(b >= 1)
    {
      BrickorCES = 1;
      NPlate = (detID - NWall*1E7 -NRow*1E6 - NColumn*1E4 - BrickorCES*1E3)/1E1;
//      NPlate = detID - NWall*1E7 -NRow*1E6 - NColumn*1E4 - BrickorCES*1E3;
    }
  EmTop = (detID - NWall*1E7 -NRow*1E6 - NColumn*1E4- BrickorCES*1E3- NPlate*1E1)/1E0;
  if(BrickorCES == 0)
    {
      EmCES = 1; EmBrick =0;
    }
  if(BrickorCES == 1)
    {
      EmBrick = 1; EmCES =0;
    }
  
  // cout << "NPlate = " << NPlate << ";  NColumn = " << NColumn << ";  NRow = " << NRow << "; NWall = " << NWall << endl;
  // cout << "BrickorCES = " << BrickorCES <<endl;
  // cout << "EmCES = " << EmCES << ";    EmBrick = " << EmBick << endl;
  // cout << endl;
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
bool GeometricalEfficiency(TVector3 Vn, double offsetxy, int Nminplates, int &InteractionWall){

    double vx = Vn(0);
    double vy = Vn(1);
    double vz = Vn(2);

    const int Ntotplates = 56;
    double startbrickx = 0.5;
    double startbricky = 0.5;

    int whichrow, whichcolumn; //for findbrick function
    int whichplate = 60;

    FindBrick(vx, vy, vz, InteractionWall, whichrow, whichcolumn, whichplate);
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

bool NeutrinoVertexLocation(int trackID, const ShipMCTrack& track, vector<int> &primary, vector<int> &primaryvisible, ROOT::RVec<int> signalpdgs){
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
          
          if(signalpdgs[signalpdgs==TMath::Abs(pdgcode)].size() == 0) primary.push_back(trackID);
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

          if (TMath::Abs(pdgcode)==13){
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
vector<double> eff_formula(int foundweight, int totalweight, int Nevents_total){
  vector<double> efficiency; //value and error
  
  efficiency.push_back((double) foundweight/totalweight);

  //totalweight and foundweight are weighted, I need to divide with the actual number of events simulated!
  double efferr = TMath::Sqrt(efficiency[0] * (1- efficiency[0])/Nevents_total);
  efficiency.push_back(efferr);
  
  return efficiency;
  
}