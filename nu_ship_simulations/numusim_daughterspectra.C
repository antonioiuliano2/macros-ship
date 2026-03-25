//reading numu simulation, after cut in acceptance, inspecting daughter distributions (March 7 2025)

double GetParticleCharge (int pdgcode, TDatabasePDG *pdg){
  //from PDG, get charge
  double charge = 0.;
  if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
  else if (pdgcode > 1e+8) charge = 1.; //test storing heavy nuclei
  return charge;
}

ROOT::RVec<Double_t> MaxAperture(ROOT::RVec<Double_t> Tracks_TX, ROOT::RVec<Double_t> Tracks_TY)
{
  //finding max aperture with a nested loop, as in FEDRA EdbVertex::MaxAperture:
  double aper=0.;
  int  ntr = Tracks_TX.size();
  ROOT::RVec<Double_t> apertures;
  double tx=0,ty=0,a=0;
  double tx_1, ty_1, tx_2, ty_2;
  for (int i=0; i<ntr-1; i++) {
    //first daughter
    tx_1 = Tracks_TX[i];
    ty_1 = Tracks_TY[i];

    for (int j=i+1; j<ntr; j++) {
      //second daughter

      tx_2 = Tracks_TX[j];
      ty_2 = Tracks_TY[j];

      tx= tx_1 - tx_2;
      ty= ty_1 - ty_2;
      a = TMath::Sqrt( tx*tx+ty*ty );
      apertures.push_back(a);
    }
  }
  return apertures; //ndr this is tan(theta)
}

void numusim_daughterspectra(){
 
 TFile *inputfile = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/sim_ecn3ship/condor_sims/2024_12_11_numu_CCDIS_mushield/inECC_ship.conical.Genie-TGeant4.root");
 TTree *simtree = (TTree*) inputfile->Get("cbmsim");
 
 TTreeReader reader(simtree);

 TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");

 TH1D * hnuE = new TH1D("hnuE","Neutrino energy;E[GeV]",400,0,400);
 

 //TFile *outputfile = new TFile("nu_and_chargeddaughters.root","RECREATE");
 //TNtuple * nuvertex_info = new TNtuple("nuvertex","Information of neutrino (track 0) and daughters from numu CCDIS","ievent:px:py:pz:x:y:z:weight");

 TH2D *htxty = new TH2D("htxty","Angular distribution of neutrino daughters;TX;TY",200,-1.,1.,200,-1.,1.);
 TH1D *hmax_aperture = new TH1D("hmax_aperture","Maximum angular aperture of neutrino daughters;max_aperture[rad]",40,0,2.);
 TH1D *hmean_aperture = new TH1D("hmean_aperture","Mean angular aperture of neutrino daughters;max_aperture[rad]",40,0,2.);
 

 const int nentries = reader.GetEntries(true);
 cout<<"Number of events"<<nentries<<endl;
 TDatabasePDG *pdg = TDatabasePDG::Instance();

 for(int ientry = 0;ientry<nentries;ientry++){   
  reader.SetEntry(ientry);
  Double_t weight = tracks[0].GetWeight(); //event weight
  ROOT::RVec<Double_t> NuChargedDaughters_TX;
  ROOT::RVec<Double_t> NuChargedDaughters_TY;

  int trackID = 0;

  for (const ShipMCTrack& track: tracks){

   int pdgcode = track.GetPdgCode();
   double charge = GetParticleCharge(pdgcode, pdg);
   int motherid = track.GetMotherId();
   //charged neutrino daughter
   if (TMath::Abs(charge) > 0. && motherid ==0){  
    
    Double_t TX = track.GetPx()/track.GetPz();
    Double_t TY = track.GetPy()/track.GetPz();

    if (trackID > 1 && track.GetP() > 0.1 && track.GetPz() > 0.){ //momentum larger than 100 MeV, no initial lepton, no backscattering

    NuChargedDaughters_TX.push_back(TX);
    NuChargedDaughters_TY.push_back(TY);
    
    htxty->Fill(TX,TY);
    }
    //storing track information
    //nuvertex_info->Fill(ientry, track.GetPx(), track.GetPy(), track.GetPz(), track.GetStartX(), track.GetStartY(), track.GetStartZ(), weight);

    
   }
   trackID++;
  }

  ROOT::RVec<Double_t> apertures = MaxAperture(NuChargedDaughters_TX, NuChargedDaughters_TY);
  //max and mean aperture
  if(apertures.size() >= 2){
   hmax_aperture->Fill(TMath::ATan(ROOT::VecOps::Max(apertures))); 
   hmean_aperture->Fill(TMath::ATan(ROOT::VecOps::Mean(apertures)));
  }
 } 

 TCanvas *c_txty = new TCanvas("ctxty","2D Angular distribution",800,800);
 //storing ntuple and closing outputfile
 htxty->Draw("COLZ");
 TCanvas *c_aperture = new TCanvas("c_aperture","Maximum aperture");
 hmax_aperture->Draw(); 
 
 TCanvas *c_meanaperture = new TCanvas("c_meanaperture","Mean aperture");
 hmean_aperture->Draw(); 

 


}
