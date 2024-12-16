//Script to evaluate the fraction of nu events contained in the target (A.I. 13 December 2024)
double GetParticleCharge (int pdgcode, TDatabasePDG *pdg){
  //from PDG, get charge
  double charge = 0.;
  if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge();
  else if (pdgcode > 1e+8) charge = 1.; //test storing heavy nuclei
  return charge;
}

void containment_nushower(){
    TDatabasePDG *pdg = TDatabasePDG::Instance(); //database of particles
    TFile *simfile = TFile::Open("inECC_ship.conical.Genie-TGeant4.root");
    TTree *simtree = (TTree*) simfile->Get("cbmsim");

    TTreeReader reader(simtree);
    TTreeReaderArray<ShipMCTrack> tracks(reader,"MCTrack");
    TTreeReaderArray<TargetPoint> targetpoints(reader,"TargetPoint");

    const int nentries = reader.GetEntries(true);

    cout<<"Start loop over number of entries "<<nentries<<endl;

    TH1D *hnu_vz = new TH1D("hnu_vz","VZ events;vz[cm]",20,-3700,-3500);
    TH1D *hnu_vz_contained = new TH1D("hnu_vz_contained","VZ events contained;vz[cm]",20,-3700,-3500);

    TH1D *hnu_p = new TH1D("hnu_p","Neutrino momenta;P[GeV/c]",40,0,400);
    TH1D *hnu_p_contained = new TH1D("hnu_p_contained","Neutrino momenta contained;P[GeV/c]",40,0,400);

    TH2D *hnu_vz_p = new TH2D("hnu_vz_p","Neutrino momenta and vz;vz[cm];P[GeV/c]",20,-3700,-3500,40,0,400);
    TH2D *hnu_vz_p_contained = new TH2D("hnu_vz_p_contained","Neutrino momenta and vz contained;vz[cm];P[GeV/c]",20,-3700,-3500,40,0,400);

    //reset indexes
    int ngoodevents = 0;
    int ncontainedevents = 0;

    for(int ientry = 0;ientry<nentries;ientry++){
        reader.SetEntry(ientry);
        //coordinates of event generations
        double nu_vx = tracks[0].GetStartX();
        double nu_vy = tracks[0].GetStartY();
        double nu_vz = tracks[0].GetStartZ();

        double nu_p = tracks[0].GetP();
        //cut of not good events
        if (nu_vx < 0.) continue;

        hnu_vz->Fill(nu_vz);
        hnu_p->Fill(nu_p);
        hnu_vz_p->Fill(nu_vz,nu_p);
        ngoodevents++;
        bool hitinlast = false;
        for (const TargetPoint& targetpoint:targetpoints){
            int pdgcode = targetpoint.PdgCode();
            int detectorID = targetpoint.GetDetectorID();

            double targetz = targetpoint.GetZ();

            double charge = GetParticleCharge(pdgcode, pdg);
            //only charged hits considered
            if (TMath::Abs(charge)>0.){
              if (targetz > -3600 && detectorID >= 10000 && detectorID <= 10002){  //does the event reach the last station?
                if (TMath::Abs(pdgcode) != 13) hitinlast=true;
              }
             }//end if charge             
        } //end for target points
        if (!hitinlast){
            ncontainedevents++; //no hit in last stations, the event is contained

            hnu_vz_contained->Fill(nu_vz);
            hnu_p_contained->Fill(nu_p);
            hnu_vz_p_contained->Fill(nu_vz,nu_p);

        } 
    } //end event loop
    cout<<"End of script "<<endl;
    cout<<"contained events "<<ncontainedevents<<" over good events "<<ngoodevents<<endl;
    cout<<"fraction contained events"<<(double) ncontainedevents/ngoodevents<<endl;
    //plotting containment studies
    TCanvas *c1 = new TCanvas();
    hnu_vz->Draw("histo");
    hnu_vz_contained->SetLineColor(kRed);
    hnu_vz_contained->Draw("histo&&SAMES");
    c1->BuildLegend();

    TCanvas *c2 = new TCanvas();
    hnu_p->Draw("histo");
    hnu_p_contained->SetLineColor(kRed);
    hnu_p_contained->Draw("histo&&SAMES");
    c2->BuildLegend();

    //efficiency plots
    TEfficiency *heff_vz_p = new TEfficiency (*hnu_vz_p_contained, *hnu_vz_p);
    TCanvas *c3 = new TCanvas();
    //hnu_vz_p_contained->Draw("COLZ");
    heff_vz_p->Draw("COLZ");

    TEfficiency *heff_vz = new TEfficiency(*hnu_vz_contained,*hnu_vz);
    TCanvas *c4 = new TCanvas();
    heff_vz->Draw();

    TEfficiency *heff_p = new TEfficiency(*hnu_p_contained,*hnu_p);
    TCanvas *c5 = new TCanvas();
    heff_p->Draw();
}