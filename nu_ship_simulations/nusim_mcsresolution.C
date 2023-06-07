//Estimate energy resolution via MCS formula from OPERA paper parametrization.
//Created by A.Iuliano in 27 May 2023

double GetParticleCharge (int pdgcode, TDatabasePDG *pdg){
  //from PDG, get charge
  double charge = 0.;
  if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
  else if (pdgcode > 1e+8) charge = 1.; //test storing heavy nuclei
  return charge;
}

void nusim_mcsresolution(){
 const int nplmin = 3; //otherwise we do not accept the track
 TF2 *fOPERA_mcsres = new TF2("fOPERA_mcsres","([0]+[1]*x)/TMath::Sqrt(y)+ ([2]+[3]*x)+([4]+[5]*x)*TMath::Sqrt(y)");
 //parameters from paper New Journal of Physics 14 (2012) 013026 
 fOPERA_mcsres->SetParameters(0.397,0.019,0.176,0.042,-0.014,-0.003);
 TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
 TString prefix("root:://eosuser.cern.ch/"); //if needed for ROOTXD
 TString simpath("/eos/user/a/aiuliano/public/sims_FairShip/sim_nutaudet/2023_05_26_numu_CCDIS_spectrosagitta/inECC_ship.conical.Genie-TGeant4.root");
 //we need a simulation with TargetPoints active
 TFile *inputfile = TFile::Open((prefix+simpath).Data());
 TTree *simtree = (TTree*) inputfile->Get("cbmsim");

 TTreeReader reader(simtree);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
 TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");

 const int nentries = reader.GetEntries(true);

 map<int,int> npl_nudaughterID;

 TH1D *hnpl = new TH1D("hnpl","npl of tracks from neutrino interactions",61,0,61);
 TH2D *hp_npl = new TProfile2D("hp_npl","Momentum and track length;P[GeV/c];npl",200,0,200,61,0,61);
 TProfile *gtotres = new TProfile("gtotres","Total resolution distribution;Enu[GeV];sigmaP/P",400,0,400,0.,10.);
 TProfile2D *gres = new TProfile2D("gres","Resolution distribution;P[GeV/c];npl",200,0,200,61,0,61,0.,10.); 

 double sigmaP, sigma1overP;

 //starting the loop
 cout<<"Start loop over events "<<nentries<<endl;
 for(int ientry = 0;ientry<nentries;ientry++){
    TVector3 totmeasuredP(0.,0.,0.);
    double totvarP = 0.;

    reader.SetEntry(ientry);
    if(targetpoints.GetSize()==0) continue;//no emulsion hits in this entry, useless for this code
    double weight = tracks[0].GetWeight();  
    int itrack = 0;
    int vertexbrickID = targetpoints[0].GetDetectorID()/1e+4; //the first hit stored is the closest temporarily to nu interaction
    //loop over mctracks
    npl_nudaughterID.clear(); //clearing map
    for (const ShipMCTrack& track: tracks){
     //select neutrino daughters and prepare a list
     int pdgcode = track.GetPdgCode();
     double charge = GetParticleCharge(pdgcode, pdg);
     //charged neutrino daughters (all, except muons)
     if (track.GetMotherId()==0 && TMath::Abs(charge) > 0 && TMath::Abs(track.GetPdgCode())!=13){ 
      npl_nudaughterID[itrack] = 0;
     }
     itrack++;
    }//ending the track loop
    //loop over target points
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
   //filling npl
   for ( const auto &myPair : npl_nudaughterID ) {
         int npl = myPair.second;
         hnpl->Fill(npl);
         int trackID = myPair.first;
         double momentum = tracks[trackID].GetP();
         TVector3 trimomentum(tracks[trackID].GetPx(),tracks[trackID].GetPy(),tracks[trackID].GetPz());
         double energy = tracks[trackID].GetEnergy();        
         if(npl<nplmin) continue;
         //sum together vectors of observed charged neutrino daughters
         totmeasuredP = totmeasuredP + trimomentum;
         double res = fOPERA_mcsres->Eval(momentum, npl); //sigma(1/P)/(1/P)
         sigma1overP = res * (1./momentum);
         sigmaP = sigma1overP * momentum * momentum; //from error propagation 

         //cout<<"evento "<<ientry<<" traccia "<<trackID<<"pdgcode: " <<tracks[trackID].GetPdgCode()<<" motherID: "<<tracks[trackID].GetMotherId()<<" energy: "<<tracks[trackID].GetEnergy()<<" npl "<<npl<<" res "<<res<<endl;         
         totvarP = totvarP + pow(sigmaP,2);
         gres->Fill(momentum,npl,res);         
         hp_npl->Fill(momentum,npl,weight);                  
   }
   gtotres->Fill(tracks[0].GetP(), TMath::Sqrt(totvarP)/totmeasuredP.Mag());
 // cout<<"Test event "<<ientry<<" MC nu energy "<<tracks[0].GetP()<<" charged measured energy "<<totmeasuredP<<" pm "<<TMath::Sqrt(totvarP)<<endl;
 }//ending the event loop
 TCanvas *cnpl = new TCanvas();
 hnpl->Draw();

 TCanvas *cgres = new TCanvas();
 gres->Draw("COLZ");

 TCanvas *cgtotres = new TCanvas();
 gtotres->Draw();

 TCanvas *cp_npl = new TCanvas();
 hp_npl->Scale(1./hp_npl->Integral());
 hp_npl->Draw("COLZ"); 
}