//before launching, load classes with
/*.L ValeDecaySearch.C
  .L FEDRAVertices.C
*/
R__LOAD_LIBRARY(/home/antonio/Scrivania/macros-ship/DecaySearchKinematics/ShipCharmDecaySearch_C.so)
R__LOAD_LIBRARY(/home/antonio/Dottorato/CharmData/decay_search_MC/spread_beam/firstquarter/FEDRAVertices_C.so)
R__LOAD_LIBRARY(/home/antonio/Dottorato/CharmData/decay_search_MC/spread_beam/firstquarter/ValeDecaySearch_C.so)
//histogram definitions
TH1D *hvz = new TH1D("hvz","Z position of vertices;vz[#mum]", 100, -100000,0);
TH1D *hdlength = new TH1D("hdlength","Decay length;dl[#mum]",100,0,100000);

TH1D *hIP = new TH1D("hIP","Impact parameter with respect to primary vertex",100,0,1000);

TH1D *hmaxaperture = new TH1D("hmaxaperture","Max aperture of a vertex;rad",180,0,1.8);

TH1D *hkink = new TH1D("hkink","Kink angle;rad",400,0,4);
TH1D *haveragekink = new TH1D("haveragekink","Average Kink angle;rad",400,0,4);
TH1D *hdphi = new TH1D("hdphi","Phi angle difference;#Delta#phi[rad]",400,0,4);
void fillevent(int primaryvID, int firstcharmvID, int secondcharmvID);

//getting files and classes with trees

FEDRAVertices recovertices;

void candidate_distributions(){
  ValeDecaySearch dsresult;
  TString *candidate_selection = new TString("dsvtx.nfound>=2");

  dsresult.fChain->Draw(">>candidates", candidate_selection->Data());

  TEventList *candidates = (TEventList*) gDirectory->GetList()->FindObject("candidates");
  const int ncandidates = dsresult.fChain->GetEntries(candidate_selection->Data());

  cout<<"Start to analyze "<<ncandidates<<" events "<<endl;
  //starting loop on candidates
  int atleasttwo = 0;
  for (int icandidate = 0; icandidate < ncandidates;icandidate++){
      int entr = candidates->GetEntry(icandidate);
      dsresult.GetEntry(entr);
      int primaryvID = dsresult.vtx_fe_id;
      //int firstcharmvID = 0;
      //int secondcharmvID = 0;
      if (dsresult.dsvtx_vtx2_vid->size() < 2) continue;
      atleasttwo++;
      int firstcharmvID = dsresult.dsvtx_vtx2_vid->at(0);
      int secondcharmvID = dsresult.dsvtx_vtx2_vid->at(1);
      cout<<"At entry "<<entr<<"Primary vertex ID: "<<primaryvID<<" Charm Vertices IDs: "<<firstcharmvID<<" "<<secondcharmvID<<endl;
      fillevent(primaryvID,firstcharmvID,secondcharmvID);
  }
  cout<<"Candidates with at least two secondary vertices: "<<atleasttwo<<" out of "<<ncandidates<<endl;
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

void fillevent(int primaryvID, int firstcharmvID, int secondcharmvID){
    //first, position of primary vertex
//    gROOT->ProcessLine(".L /home/antonio/Scrivania/macros-ship/DecaySearchKinematics/ShipCharmDecaySearch_C.so");
    ShipCharmDecaySearch dsfunctions;
    recovertices.GetEntry(primaryvID);
    float prim_vx = recovertices.vx;
    float prim_vy = recovertices.vy;
    float prim_vz = recovertices.vz;

    TVector3 prim_vertexpos = TVector3(prim_vx,prim_vy,prim_vz);
    
   
    
    float charm_vx, charm_vy, charm_vz;
    float charm_tx, charm_ty;
  
    float charm_phi[2]; //need to keep the angles for both charms, in order to compute the difference

    int icharm = 0;
    int ntracks, ndaughters;
    ROOT::VecOps::RVec<float> daughters_tx, daughters_ty;
    //loop in the two vertices
    for (int vID:{firstcharmvID,secondcharmvID}){
    //then, position of secondary vertex
     recovertices.GetEntry(vID);
     charm_vx = recovertices.vx;
     charm_vy = recovertices.vy;
     charm_vz = recovertices.vz;
  
        //estimated angle
     charm_tx = (charm_vx - prim_vx)/(charm_vz-prim_vz);
     charm_ty = (charm_vy - prim_vy)/(charm_vz-prim_vz);
     charm_phi[icharm] = TMath::ATan2(charm_ty,charm_tx);
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

    //filling histograms
    hmaxaperture->Fill(recovertices.maxaperture);
    hvz->Fill(prim_vz);
    hdlength->Fill(TMath::Sqrt(pow(charm_vx-prim_vx,2)+pow(charm_vy-prim_vy,2)+pow(charm_vz-prim_vz,2)));
    //passing counter to next charm
    icharm++;
 }
 hdphi->Fill(TMath::Abs(charm_phi[1]-charm_phi[0]));
}