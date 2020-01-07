//before launching, load classes with
/*.L ValeDecaySearch.C
  .L FEDRAVertices.C
*/
R__LOAD_LIBRARY(/home/antonio/Scrivania/macros-ship/DecaySearchKinematics/ShipCharmDecaySearch_C.so)
//R__LOAD_LIBRARY(/home/antonio/Dottorato/CharmData/decay_search_data/FEDRAVertices_C.so)
//R__LOAD_LIBRARY(/home/antonio/Dottorato/CharmData/decay_search_data/ValeDecaySearch_C.so)

R__LOAD_LIBRARY(/home/antonio/Dottorato/CharmData/decay_search_MC/CH1R6_pot_03_12_19/secondquarter/FEDRAVertices_C.so)
R__LOAD_LIBRARY(/home/antonio/Dottorato/CharmData/decay_search_MC/CH1R6_pot_03_12_19/secondquarter/ValeDecaySearch_C.so)

//histogram definitions

TFile *histofile = new TFile("histosdecays.root","RECREATE");
TH1D *hvz = new TH1D("hvz","Z position of vertices;vz[#mum]", 50, -40000,0);
TH1D *hprob_prim = new TH1D("hprob_prim","Probability of primary vertices;Probability",50,0,1);
//true
TH1D *hnsecondaries= new TH1D("hnsecondaries","Number of secondary vertices",10,0,10);

TH1D *hprob_secondary_true = new TH1D("hprob_secondary_true","Probability of secondary vertices;Probability",100,0,1);

TH1D *hvz_secondary_true = new TH1D("hvz_secondary_true","Position of secondary vertices;vz[#mum]",50,-40000,0);
TH1D *hdlength_true = new TH1D("hdlength_true","Decay length;dl[#mum]",22,0,11000);

TH1D *hIP_true = new TH1D("hIP_true","Impact parameter with respect to primary vertex",100,0,10000);

TH1D *hmaxaperture_true = new TH1D("hmaxaperture_true","Max aperture of a vertex;rad",90,0,1.8);

TH1D *hkink_true = new TH1D("hkink_true","Kink angle;rad",100,0,1);
TH1D *haveragekink_true = new TH1D("haveragekink_true","Average Kink angle;rad",100,0,1);
TH1D *hdphi_true = new TH1D("hdphi_true","Phi angle difference;#Delta#phi[rad]",400,0,4);
//fake

TH1D *hprob_secondary_fake = new TH1D("hprob_secondary_fake","Probability of secondary vertices;Probability",50,0,1);

TH1D *hvz_secondary_fake = new TH1D("hvz_secondary_fake","Position of secondary vertices;vz[#mum]",50,-40000,0);
TH1D *hdlength_fake = new TH1D("hdlength_fake","Decay length;dl[#mum]",22,0,11000);

TH1D *hIP_fake = new TH1D("hIP_fake","Impact parameter with respect to primary vertex",100,0,10000);

TH1D *hmaxaperture_fake = new TH1D("hmaxaperture_fake","Max aperture of a vertex;rad",90,0,1.8);

TH1D *hkink_fake = new TH1D("hkink_fake","Kink angle;rad",100,0,1);
TH1D *haveragekink_fake = new TH1D("haveragekink_fake","Average Kink angle;rad",100,0,1);
TH1D *hdphi_fake = new TH1D("hdphi_fake","Phi angle difference;#Delta#phi[rad]",400,0,4);

void fillevent(int primaryvID, vector<int>* charmvID);
void filltrackevent (int primaryvID, int charmvID, int trackID);
//getting files and classes with trees

FEDRAVertices recovertices;
/*
void from_textfile(){
 //opening text file
 fstream inputdatfile;
 inputdatfile.open("/home/antonio/Dottorato/CharmData/decay_search_data/candidateslist_firstquarter.dat",std::fstream::in);

 int topology;
 int primaryvID;
 int firstcharmvID;
 int secondcharmvID; 

 int dummy;
 //for (int idummy = 0; idummy < 4; idummy++) dummy = inputdatfile.get();

 while (inputdatfile.good()){
   //getting values from line
   inputdatfile>>topology;
   inputdatfile>>primaryvID;
   inputdatfile>>firstcharmvID;
   inputdatfile>>secondcharmvID; 
   //filling the histograms
   if (topology == 21) filltrackevent(primaryvID,firstcharmvID,secondcharmvID);
   else if (topology == 22) fillevent(primaryvID,firstcharmvID,secondcharmvID);
 }
  //drawing histograms
  TCanvas *cvz = new TCanvas();
  hvz->Draw();
  TCanvas *cdecaylength = new TCanvas();
  hdlength->Draw();
  TCanvas *cmaxaperture = new TCanvas();
  hmaxaperture->Draw();
  TCanvas *ckink = new TCanvas();
  hkink->Draw();
  haveragekink->Draw();
  TCanvas *cphi = new TCanvas();
  hdphi->Draw();
  TCanvas *cIP = new TCanvas();
  hIP->Draw();
}
*/
void candidate_distributions(){
  ValeDecaySearch dsresult;
  TString *candidate_selection = new TString("1");
  dsresult.fChain->Draw(">>candidates", candidate_selection->Data());

  TEventList *candidates = (TEventList*) gDirectory->GetList()->FindObject("candidates");
  const int ncandidates = dsresult.fChain->GetEntries(candidate_selection->Data());

  cout<<"Start to analyze "<<ncandidates<<" events "<<endl;
  //starting loop on candidates
  int atleastone = 0;
  for (int icandidate = 0; icandidate < ncandidates;icandidate++){
      int entr = candidates->GetEntry(icandidate);
      dsresult.GetEntry(entr);
      int primaryvID = dsresult.vtx_fe_id;
      //int firstcharmvID = 0;
      //int secondcharmvID = 0;
      hnsecondaries->Fill(dsresult.dsvtx_vtx2_vid->size());
      if (dsresult.dsvtx_vtx2_vid->size() < 1) continue;
      atleastone++;
      int firstcharmvID = dsresult.dsvtx_vtx2_vid->at(0);
      //int secondcharmvID = dsresult.dsvtx_vtx2_vid->at(1);
      cout<<"At entry "<<entr<<"Primary vertex ID: "<<primaryvID<<" Charm Vertices IDs: "<<firstcharmvID<<" "<<dsresult.dsvtx_vtx2_vid->size()<<endl;
      fillevent(primaryvID,dsresult.dsvtx_vtx2_vid);
  }
  cout<<"Candidates with at least two secondary vertices: "<<atleastone<<" out of "<<ncandidates<<endl;
  TCanvas *csecondaries = new TCanvas();
  hnsecondaries->Draw();
  TCanvas *cvz = new TCanvas();
  cvz->cd(1);
  hvz->Draw();
  hvz->Scale(1./hvz->Integral());
  cvz->cd(2);
  hvz_secondary_fake->Draw();
  hvz_secondary_fake->Scale(1./hvz_secondary_fake->Integral());
  hvz_secondary_true->SetLineColor(kRed);
  hvz_secondary_true->Draw("SAMES");
  hvz_secondary_true->Scale(1./hvz_secondary_true->Integral());
  TCanvas *cdecaylength = new TCanvas();
  hdlength_fake->Draw();
  hdlength_fake->Scale(1./hdlength_fake->Integral());
  hdlength_true->SetLineColor(kRed);
  hdlength_true->Scale(1./hdlength_true->Integral());
  hdlength_true->Draw("SAMES");
  TCanvas *cmaxaperture = new TCanvas();
  hmaxaperture_fake->Scale(1./hdlength_fake->Integral());
  hmaxaperture_fake->Draw();
  hmaxaperture_true->Scale(1./hmaxaperture_true->Integral());
  hmaxaperture_true->SetLineColor(kRed);
  hmaxaperture_true->Draw("SAMES");
  TCanvas *ckink = new TCanvas();
  ckink->Divide(1,2);
  ckink->cd(1);
  hkink_fake->Draw();
  hkink_fake->Scale(1./hkink_fake->Integral());
  hkink_true->SetLineColor(kRed);
  hkink_true->Scale(1./hkink_true->Integral());
  hkink_true->Draw("SAMES");
  ckink->cd(2);
  haveragekink_fake->Draw();
  haveragekink_fake->Scale(1./haveragekink_fake->Integral());
  haveragekink_true->SetLineColor(kRed);
  haveragekink_true->Scale(1./haveragekink_true->Integral());
  haveragekink_true->Draw("SAMES");
  //TCanvas *cphi = new TCanvas();
  //hdphi->Draw();
  TCanvas *cIP = new TCanvas();
  hIP_fake->Draw();
  hIP_fake->Scale(1./hIP_fake->Integral());
  hIP_true->SetLineColor(kRed);
  hIP_true->Draw("SAMES");
  hIP_true->Scale(1./hIP_true->Integral());
  TCanvas *cProb = new TCanvas();
  cProb->Divide(1,2);
  cProb->cd(1);
  hprob_prim->Scale(1./hprob_prim->Integral());
  hprob_prim->Draw();
  cProb->cd(2);
  hprob_secondary_fake->Draw();
  hprob_secondary_fake->Scale(hprob_secondary_fake->Integral());
  hprob_secondary_true->SetLineColor(kRed);
  hprob_secondary_true->Scale(hprob_secondary_true->Integral());
  hprob_secondary_true->Draw("SAMES");
  histofile->Write();
}
/*
void filltrackevent (int primaryvID, int charmvID, int trackID){
    ShipCharmDecaySearch dsfunctions;
    recovertices.GetEntry(primaryvID);
    float prim_vx = recovertices.vx;
    float prim_vy = recovertices.vy;
    float prim_vz = recovertices.vz;
  
    float charm_tx,charm_ty;

    float firstsegmenttx, firstsegmentty;
    float secondsegmenttx, secondsegmentty;
    float kinkx, kinky, kinkz;
    float afterkinktx, afterkinkty;

    float deltathetamax, deltathetarms;

    TVector3 prim_vertexpos = TVector3(prim_vx,prim_vy,prim_vz);
    int ntracks = recovertices.n;
    float charm_phi[2]; //need to keep the angles for both charms, in order to compute the difference
    //look for kink point
    int nprevioussegments = 0;
    TVector3 kinkpos;
     for (int itrk = 0; itrk < ntracks; itrk++){
       if (recovertices.t__eTrack[itrk] == trackID){ //track of interest
        int nsegments = recovertices.nseg[itrk];
        float kinkangles[nsegments-1];
        for (int iseg = 0; iseg<nsegments-1;iseg++){
          firstsegmenttx = recovertices.s_eTX[(nprevioussegments + iseg)];
          firstsegmentty = recovertices.s_eTY[(nprevioussegments + iseg)];
          secondsegmenttx = recovertices.s_eTX[(nprevioussegments + iseg+1)];
          secondsegmentty = recovertices.s_eTY[(nprevioussegments + iseg+1)];
          kinkangles[iseg]=dsfunctions.KinkAngle(firstsegmenttx, firstsegmentty, secondsegmenttx, secondsegmentty);
        }
        //finding maximum and maximum location
        deltathetamax = TMath::MaxElement(nsegments-1, kinkangles);
        int locmax = TMath::LocMax(nsegments-1,kinkangles);
        deltathetarms = TMath::RMS(nsegments-1, kinkangles);
        float rmax = deltathetamax/deltathetarms;
        //angles immediately after kink
        afterkinktx = recovertices.s_eTX[(nprevioussegments + locmax+1)];
        afterkinkty = recovertices.s_eTY[(nprevioussegments + locmax+1)];
        //kink position
        kinkx = recovertices.s_eX[(nprevioussegments + locmax+1)];
        kinky = recovertices.s_eY[(nprevioussegments + locmax+1)];
        kinkz = recovertices.s_eZ[(nprevioussegments + locmax+1)];
        kinkpos = TVector3(kinkx,kinky,kinkz);
       }
       nprevioussegments = nprevioussegments + itrk * recovertices.nseg[itrk]; //increasing the counter for index;
     }
        //estimated angle
     charm_tx = (kinkx - prim_vx)/(kinkz-prim_vz);
     charm_ty = (kinky - prim_vy)/(kinkz-prim_vz);
     charm_phi[0] = TMath::ATan2(charm_ty,charm_tx);
     //loop on vertex tracks, keep only daughters
     ntracks = recovertices.n;
     hkink->Fill(deltathetamax);
     hIP->Fill(dsfunctions.IPtoVertex(prim_vertexpos,kinkpos,afterkinktx,afterkinkty));         
      //one prong, averagekink is equal to kink
     haveragekink->Fill(deltathetamax);

    //filling decay length
    hdlength->Fill(TMath::Sqrt(pow(kinkx-prim_vx,2)+pow(kinky-prim_vy,2)+pow(kinkz-prim_vz,2)));

    float charm_vx, charm_vy, charm_vz;

    int ndaughters;
    ROOT::VecOps::RVec<float> daughters_tx, daughters_ty;
    //then, position of secondary vertex
    recovertices.GetEntry(charmvID);
    charm_vx = recovertices.vx;
    charm_vy = recovertices.vy;
    charm_vz = recovertices.vz;
  
        //estimated angle
    charm_tx = (charm_vx - prim_vx)/(charm_vz-prim_vz);
    charm_ty = (charm_vy - prim_vy)/(charm_vz-prim_vz);
    charm_phi[1] = TMath::ATan2(charm_ty,charm_tx);
     //resetting counters and indices for vertex daughters
    daughters_tx.clear();
    daughters_ty.clear();
     //loop on vertex tracks, keep only daughters
    ntracks = recovertices.n;
    for (int itrk = 0; itrk < ntracks; itrk++){
        if (recovertices.incoming[itrk] == 1){         
         TVector3 startpos = TVector3(recovertices.t__eX[itrk],recovertices.t__eY[itrk],recovertices.t__eZ[itrk]); //first segment of the track

         daughters_tx.push_back(recovertices.t__eTX[itrk]);
         daughters_ty.push_back(recovertices.t__eTY[itrk]);

         hkink->Fill(dsfunctions.KinkAngle(charm_tx,charm_ty,recovertices.t__eTX[itrk],recovertices.t__eTY[itrk]));
         hIP->Fill(dsfunctions.IPtoVertex(prim_vertexpos,startpos,recovertices.t__eTX[itrk],recovertices.t__eTY[itrk]));         
        }
        //RVec::data() returns the address with the contiguous memory, read as an array
        ndaughters = daughters_tx.size();
        haveragekink->Fill(dsfunctions.AverageKinkAngle(charm_tx,charm_ty,daughters_tx.data(),daughters_ty.data(),ndaughters));
        
    }
  //filling histogram
  hmaxaperture->Fill(recovertices.maxaperture);
  hvz->Fill(prim_vz);
  hdlength->Fill(TMath::Sqrt(pow(charm_vx-prim_vx,2)+pow(charm_vy-prim_vy,2)+pow(charm_vz-prim_vz,2)));
    
 hdphi->Fill(TMath::Abs(charm_phi[1]-charm_phi[0]));
} 
*/
void fillevent(int primaryvID, vector<int> *charmvID){
    //first, position of primary vertex
//    gROOT->ProcessLine(".L /home/antonio/Scrivania/macros-ship/DecaySearchKinematics/ShipCharmDecaySearch_C.so");
    ShipCharmDecaySearch dsfunctions;
    recovertices.GetEntry(primaryvID);
    float prim_vx = recovertices.vx;
    float prim_vy = recovertices.vy;
    float prim_vz = recovertices.vz;

    TVector3 prim_vertexpos = TVector3(prim_vx,prim_vy,prim_vz);
    
    float prim_probability = recovertices.probability;
    
    float charm_vx, charm_vy, charm_vz;
    float charm_tx, charm_ty;
  
//    float charm_phi[2]; //need to keep the angles for both charms, in order to compute the difference

    int icharm = 0;
    bool truecharm = false;
    int ntracks, ndaughters;
    ROOT::VecOps::RVec<float> daughters_tx, daughters_ty;
    //loop in the two vertices
    for (int vID:*charmvID){
      truecharm = false;
    //then, position of secondary vertex
     recovertices.GetEntry(vID);
     charm_vx = recovertices.vx;
     charm_vy = recovertices.vy;
     charm_vz = recovertices.vz;


     float second_probability = recovertices.probability;
  
        //estimated angle
     charm_tx = (charm_vx - prim_vx)/(charm_vz-prim_vz);
     charm_ty = (charm_vy - prim_vy)/(charm_vz-prim_vz);
     //charm_phi[icharm] = TMath::ATan2(charm_ty,charm_tx);
     //resetting counters and indices for vertex daughters
     daughters_tx.clear();
     daughters_ty.clear();
     //loop on vertex tracks, keep only daughters
     ntracks = recovertices.n;
     for (int itrk = 0; itrk < ntracks; itrk++){
       int motherid = recovertices.MCMotherID[itrk];
       if ((motherid == 1 || motherid == 2) && recovertices.incoming[itrk] == 1) truecharm = true;
        if (recovertices.incoming[itrk] == 1){         
         TVector3 startpos = TVector3(recovertices.t__eX[itrk],recovertices.t__eY[itrk],recovertices.t__eZ[itrk]); //first segment of the track

         daughters_tx.push_back(recovertices.t__eTX[itrk]);
         daughters_ty.push_back(recovertices.t__eTY[itrk]);
         if (truecharm){
          hkink_true->Fill(dsfunctions.KinkAngle(charm_tx,charm_ty,recovertices.t__eTX[itrk],recovertices.t__eTY[itrk]));
          hIP_true->Fill(dsfunctions.IPtoVertex(prim_vertexpos,startpos,recovertices.t__eTX[itrk],recovertices.t__eTY[itrk]));         
         }
         else{
          hkink_fake->Fill(dsfunctions.KinkAngle(charm_tx,charm_ty,recovertices.t__eTX[itrk],recovertices.t__eTY[itrk]));
          hIP_fake->Fill(dsfunctions.IPtoVertex(prim_vertexpos,startpos,recovertices.t__eTX[itrk],recovertices.t__eTY[itrk])); 
         }
        }
     }
        //RVec::data() returns the address with the contiguous memory, read as an array
     
     
     ndaughters = daughters_tx.size();
    
    //filling histograms
    if (truecharm){
     haveragekink_true->Fill(dsfunctions.AverageKinkAngle(charm_tx,charm_ty,daughters_tx.data(),daughters_ty.data(),ndaughters));    
     hmaxaperture_true->Fill(recovertices.maxaperture);
     hvz_secondary_true->Fill(charm_vz);
     hprob_secondary_true->Fill(second_probability);
     float deltavzcut = 100.;
     if ((charm_vz - prim_vz) > deltavzcut ) hdlength_true->Fill(TMath::Sqrt(pow(charm_vx-prim_vx,2)+pow(charm_vy-prim_vy,2)+pow(charm_vz-prim_vz,2)));
    //passing counter to next charm
     
    }
    else{
     haveragekink_fake->Fill(dsfunctions.AverageKinkAngle(charm_tx,charm_ty,daughters_tx.data(),daughters_ty.data(),ndaughters));    
     hmaxaperture_fake->Fill(recovertices.maxaperture);
     hvz_secondary_fake->Fill(charm_vz);
     hprob_secondary_fake->Fill(second_probability);
     float deltavzcut = 100.;
     if ((charm_vz - prim_vz) > deltavzcut ) hdlength_fake->Fill(TMath::Sqrt(pow(charm_vx-prim_vx,2)+pow(charm_vy-prim_vy,2)+pow(charm_vz-prim_vz,2)));
    }
    icharm++;
 }
 //filled only once per double charm
 hvz->Fill(prim_vz);
 hprob_prim->Fill(prim_probability);
 //hdphi->Fill(TMath::Abs(charm_phi[1]-charm_phi[0]));
}