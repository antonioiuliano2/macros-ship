using namespace ROOT;
//using namespace ROOT::RVec;
//dz selection
RVec<int> dzselection(vector<float> vtx2_dz){
  //converting to RVec, avoid doing loop later
  RVec<float> rvtx2_dz(vtx2_dz.data(),vtx2_dz.size());
  //applying condition
  RVec<float> rvtx2_goodIDs = (rvtx2_dz >= 0.); 
  return rvtx2_goodIDs;
}
//compuation of decay length
RVec<float> decaylength(RVec<float> vtx2_vx, RVec<float> vtx2_vy, RVec<float> vtx2_vz, float vx, float vy, float vz){
  //computing directly, no need for loop
  RVec<float> rvtx2_dl = pow(pow(vtx2_vx - vx,2)+ pow(vtx2_vy-vy,2)+pow(vtx2_vz-vz,2),0.5);

  return rvtx2_dl;
}

RVec<int> atleast2starting(RVec<int> vtx2_ntracks, RVec<int>vtx2_incoming){
  RVec<int> goodvertices;
  const int nvertices = vtx2_ntracks.size();
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int ntracks = vtx2_ntracks[ivtx];
    int ngoodtrks = 0;
    for (int itrk = 0; itrk < ntracks; itrk++){
      if (vtx2_incoming[itrk+nprevioustracks] == 1) ngoodtrks++;
    }
    //at least 2 good tracks
    if (ngoodtrks >= 2) goodvertices.push_back(true);
    else goodvertices.push_back(false);
    nprevioustracks += ntracks;
  }
  return goodvertices;
}
//phi medium and difference with vertex phi
RVec<float> phi_medium(RVec<int> vtx2_ntracks, RVec<float> vtx2_track_tx, RVec<float> vtx2_track_ty){
  int nprevioustracks = 0;
  const int nvertices = vtx2_ntracks.size();
  RVec<float> phi_medium_vertex;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   int ntracks = vtx2_ntracks[ivtx];
   float vtx2_track_phi = Mean(atan2(RVec<float>(&vtx2_track_ty[nprevioustracks], &vtx2_track_ty[nprevioustracks+ntracks]),
                                     RVec<float>(&vtx2_track_tx[nprevioustracks], &vtx2_track_tx[nprevioustracks+ntracks])));
   phi_medium_vertex.push_back(vtx2_track_phi);
   //increasing ntracks counter
   nprevioustracks += ntracks;
  }
  return phi_medium_vertex;
}

RVec<float> phi_difference(RVec<float> vtx2_phi_medium, float vx, float vy, float vz, RVec<float> vtx2_vx, RVec<float> vtx2_vy, RVec<float> vtx2_vz){
  RVec<float> phi_difference_vertex;
  const int nvertices = vtx2_vx.size();
  for(int ivtx = 0; ivtx < nvertices; ivtx++){
   //computing angles
   float tx = (vtx2_vx[ivtx]-vx)/(vtx2_vz[ivtx]-vz); 
   float ty = (vtx2_vy[ivtx]-vy)/(vtx2_vz[ivtx]-vz);
   float phi = TMath::ATan2(ty,tx);
   phi_difference_vertex.push_back(phi - vtx2_phi_medium[ivtx]);
  }
  return phi_difference_vertex;
}


RVec<int> atleast2starting_trk(RVec<int> vtx2_ntracks, RVec<int>vtx2_incoming){
  RVec<int> tracks_ingoodvertices;
  const int nvertices = vtx2_ntracks.size();
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int ntracks = vtx2_ntracks[ivtx];
    int ngoodtrks = 0;
    for (int itrk = 0; itrk < ntracks; itrk++){
      if (vtx2_incoming[itrk+nprevioustracks] == 1) ngoodtrks++;
    }
    //at least 2 good tracks
    if (ngoodtrks >= 2){for (int itrk = 0; itrk < ntracks; itrk++) tracks_ingoodvertices.push_back(true);}
    else {for (int itrk = 0; itrk < ntracks; itrk++) tracks_ingoodvertices.push_back(false);}
    nprevioustracks += ntracks;
  }
  return tracks_ingoodvertices;
}

RVec<int> atleast2goodtrks(RVec<int> vtx2_ntracks, RVec<int>vtx2_nseg){
  RVec<int> goodvertices;
  const int nvertices = vtx2_ntracks.size();
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int ntracks = vtx2_ntracks[ivtx];
    int ngoodtrks = 0;
    for (int itrk = 0; itrk < ntracks; itrk++){
      if (vtx2_nseg[itrk+nprevioustracks] > 2) ngoodtrks++;
    }
    //at least 2 good tracks
    if (ngoodtrks >= 2) goodvertices.push_back(true);
    else goodvertices.push_back(false);
    nprevioustracks += ntracks;
  }
  return goodvertices;
}
RVec<int> atleast2goodtrks_trk(RVec<int> vtx2_ntracks, RVec<int>vtx2_nseg){
  RVec<int> tracks_ingoodvertices;
  const int nvertices = vtx2_ntracks.size();
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int ntracks = vtx2_ntracks[ivtx];
    int ngoodtrks = 0;
    for (int itrk = 0; itrk < ntracks; itrk++){
      if (vtx2_nseg[itrk+nprevioustracks] > 2) ngoodtrks++;
    }
    //at least 2 good tracks
    if (ngoodtrks >= 2){for (int itrk = 0; itrk < ntracks; itrk++) tracks_ingoodvertices.push_back(true);}
    else {for (int itrk = 0; itrk < ntracks; itrk++) tracks_ingoodvertices.push_back(false);}
    nprevioustracks += ntracks;
  }
  return tracks_ingoodvertices;
}

RVec<float> meanlife(RVec<int>vtx2_ntracks, RVec<float> vtx2_vka, RVec<float> vtx2_dl){
  float lightspeed = 3 * 1e+10; // cm /s
  float micron2cm = 1e-4 ;
  unsigned int nprevioustracks = 0;
  const int nvertices = vtx2_ntracks.size();
  RVec<float> vtx2_tau;
  //loop on vertices;
  //RVec<float> split_vtx2_vka = vtx2_vka; //removing tracks from each vertex
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int ntracks = vtx2_ntracks[ivtx];
    //computing mean kink of tracks from that vertex
    float meankink = Mean(RVec<float>(&vtx2_vka[nprevioustracks],&vtx2_vka[nprevioustracks+ntracks])); 
    vtx2_tau.push_back(meankink*vtx2_dl[ivtx]*micron2cm/lightspeed );
    nprevioustracks += ntracks;
  }
  return vtx2_tau;
}

RVec<int> MCsameeventtrack(RVec<int>vtx2_ntracks, int vtx_mc_event, RVec<int>vtx2_mc_event,  RVec<int> vtx2_mcparentid){
  int nvertices = vtx2_ntracks.size();
  RVec<int> samevent_track;
  int nprevioustracks = 0;
  //loop into vertices
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   //loop into tracks
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      //adding track bool information
      if (vtx2_mc_event[itrk+nprevioustracks] == vtx_mc_event) samevent_track.push_back(true);    
      else samevent_track.push_back(false);
    }//close track loop
   nprevioustracks += ntracks;
  }//close vertex loop
  return samevent_track;
}
//Track selection

RVec<int> MCcharmdaughtertrack(RVec<int>vtx2_ntracks, RVec<int>vtx2_mcparentid){
  int nvertices = vtx2_ntracks.size();
  RVec<int> charmdaughter_track;
  int nprevioustracks = 0;
  //loop into vertices
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   //loop into tracks
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      //adding track information (using properties that in boolean 0 means false, and > 0 means true. We can use it as check and store the info about which parent is)
      if (vtx2_mcparentid[itrk+nprevioustracks] == 1 || vtx2_mcparentid[itrk+nprevioustracks] == 2) charmdaughter_track.push_back(vtx2_mcparentid[itrk+nprevioustracks]);    
      else charmdaughter_track.push_back(0);
    }//close track loop
   nprevioustracks += ntracks;
  }//close vertex loop
  return charmdaughter_track;
}

RVec<int> MCvertex_samevent(RVec<int>vtx2_ntracks,RVec<int>MCsameventtrack){
  int nvertices = vtx2_ntracks.size();
  RVec<int> goodvertices;
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   bool goodvertex = false;
   //loop into tracks
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      //adding track bool information
      if (MCsameventtrack[itrk+nprevioustracks]) goodvertex = true;         
    }//close track loop
   nprevioustracks += ntracks;
   goodvertices.push_back(goodvertex);
  }//close vertex loop
  return goodvertices;
}

RVec<int> MCvertex_charmdaughter(RVec<int>vtx2_ntracks,RVec<int>MCcharmdaughter){
  int nvertices = vtx2_ntracks.size();
  RVec<int> goodvertices;
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   bool goodvertex = false;
   //loop into tracks
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      //adding track bool information
      if (MCcharmdaughter[itrk+nprevioustracks]) goodvertex = true;         
    }//close track loop
   nprevioustracks += ntracks;
   goodvertices.push_back(goodvertex);
  }//close vertex loop
  return goodvertices;
}

RVec<int> MCvertex_charmdaughter_samevent(RVec<int>vtx2_ntracks,RVec<int>MCcharmdaughtertrack,RVec<int> MCsameventtrack){
  int nvertices = vtx2_ntracks.size();
  RVec<int> with_gooddaughter;
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   int goodvertex = 0;
   //loop into tracks
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      //adding track bool information
      if (MCcharmdaughtertrack[itrk+nprevioustracks] && MCsameventtrack[itrk+nprevioustracks]) goodvertex = MCcharmdaughtertrack[itrk+nprevioustracks];         
    }//close track loop
   nprevioustracks += ntracks;
   with_gooddaughter.push_back(goodvertex);
  }//close vertex loop
  return with_gooddaughter;
}

RVec<int> dzselection_trk(vector<int>vtx2_ntracks, vector<float> vtx2_dz){
  int nvertices = vtx2_ntracks.size();
  RVec<int> dz_track;
  int nprevioustracks = 0;
  //loop into vertices
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   //loop into tracks
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      //adding track bool information
      if (vtx2_dz[ivtx]>0) dz_track.push_back(true);    
      else dz_track.push_back(false);
    }//close track loop
   nprevioustracks += ntracks;
  }//close vertex loop
  return dz_track;
}

int event_topology(RVec<int> MCvertexcharmdaughter_samevent, RVec<int> positivedz){
  RVec<int> goodcharmdaughters = MCvertexcharmdaughter_samevent * positivedz; //without a positive dz it has no sense to define it as a 'good event'

  int nfirstcharm = goodcharmdaughters[goodcharmdaughters==1].size();
  int nsecondcharm = goodcharmdaughters[goodcharmdaughters==2].size();
  if (nfirstcharm > 0 && nsecondcharm > 0 ) return 2; //"golden", double charm
  else if (nfirstcharm > 0 || nsecondcharm > 0) return 1; //"silver", single charm
  else return 0; //no good event
}

//

int mostscommonidvertex(vector<int> trk_mc_id){
    RVec<int> rtrk_mc_id(trk_mc_id.data(),trk_mc_id.size());
    if (trk_mc_id.size() > 0) return Max(rtrk_mc_id);
    else return -1;
}

RVec<float> selecteddistribution(RVec<float> originaldistribution, RVec<int> selection){
  //fill a RVec only with accepted values of the RVec
    return originaldistribution[selection];
}
//selection counters
int howmanyvertices(RVec<int> selection){
  return selection[selection>0].size();
}
//start of script
void selection_decaysearch_sim(){
    //input branches
    string vtx2_mc_pid = "dsvtx.vtx2_mc_pid";

    //secondary, vertex variables
    string vtx2_ntrk = "dsvtx.vtx2_ntrk";
    string vtx2_dz = "dsvtx.vtx2_dz";
    string vtx2_mc_ev = "dsvtx.vtx2_mc_ev";

    //newbranches
    //primary, vertex variables
    string vtx_mc_ev = "vtx_mc_ev";
    string vtx_topology = "vtx_topology";
    //secondary, track variables
    string vtx2_track_charmdaugther = "dsvtx_vtx2_trk_charmdaughter";
    string vtx2_track_samevent = "dsvtx_vtx2_trk_samevent";
    string vtx2_track_positivedz = "dsvtx_vtx2_trk_positivedz";

    //secondary, vertex variables
    string vtx2_positivedz = "dsvtx_vtx2_positivedz";
    string vtx2_samevent = "dsvtx_vtx2_samevent";
    string vtx2_charmdaugther = "dsvtx_vtx2_charmdaughter";
    string vtx2_charmdaugthersamevent = "dsvtx_vtx2_charmdaughtersamevent"; 

    //condition for good charm decay vertices
    string condition = "dsvtx_vtx2_charmdaughtersamevent*dsvtx_vtx2_positivedz*dsvtx_vtx2_2goodtrks";

    //******************START OF MAIN SCRIPT*************************//

    TFile *inputfile = TFile::Open("01_17_ds_data_result.root");
    RDataFrame dsdataframe = RDataFrame("ds",inputfile);
    //computing additional variables
    auto dflength = dsdataframe.Define("dsvtx_vtx2_dl",decaylength,{"dsvtx.vtx2_vx", "dsvtx.vtx2_vy", "dsvtx.vtx2_vz","vtx.x","vtx.y","vtx.z"});
    auto dfmeanlife = dflength.Define("dsvtx_vtx2_tau",meanlife,{"dsvtx.vtx2_ntrk","dsvtx.vtx2_vka","dsvtx_vtx2_dl"});
    auto df_primmcid = dfmeanlife.Define(vtx_mc_ev,mostscommonidvertex,{"trk.mc_ev"}); //mc event of vertex

    //angular information
    auto df_meanphi = df_primmcid.Define("dsvtx_vtx2_meanphi",phi_medium,{vtx2_ntrk,"dsvtx.vtx2_tx", "dsvtx.vtx2_ty"});
    auto df_phidifference = df_meanphi.Define("dsvtx_vtx2_phidifference",phi_difference,{"dsvtx_vtx2_meanphi","vtx.x","vtx.y","vtx.z", "dsvtx.vtx2_vx", "dsvtx.vtx2_vy", "dsvtx.vtx2_vz"});

    //checking conditions for selections
    //positive dz
    auto dfcheck_posdz_trk = df_phidifference.Define(vtx2_track_positivedz,dzselection_trk,{vtx2_ntrk,vtx2_dz});
    auto dfcheck_posdz = dfcheck_posdz_trk.Define(vtx2_positivedz,dzselection,{vtx2_dz});     

    //recognize charm daughter
    auto dfcheck_charmdaughter_trk = dfcheck_posdz.Define(vtx2_track_charmdaugther,MCcharmdaughtertrack,{vtx2_ntrk,vtx2_mc_pid });
    auto dfcheck_charmdaughter = dfcheck_charmdaughter_trk.Define(vtx2_charmdaugther, MCvertex_charmdaughter,{vtx2_ntrk,vtx2_track_charmdaugther});
    //same event condition
    auto dfcheck_samevent_trk = dfcheck_charmdaughter.Define(vtx2_track_samevent,MCsameeventtrack,{vtx2_ntrk,vtx_mc_ev,vtx2_mc_ev,vtx2_mc_pid});
    auto dfcheck_samevent = dfcheck_samevent_trk.Define(vtx2_samevent, MCvertex_samevent,{vtx2_ntrk,vtx2_track_samevent});
    
    //both charmdaughter and samevent for the same track
    auto dfcheck_charmdaughtersamevent = dfcheck_samevent.Define(vtx2_charmdaugthersamevent, MCvertex_charmdaughter_samevent,{vtx2_ntrk,vtx2_track_charmdaugther,vtx2_track_samevent});
    
    //two starting tracks
   // auto dfcheck_2starting_trk = dfcheck_charmdaughtersamevent.Define("dsvtx_vtx2_2starting_trk",atleast2starting,{vtx2_ntrk,"dsvtx.vtx2_incoming"});
   // auto dfcheck_2starting = dfcheck_2starting_trk.Define("dsvtx_vtx2_2starting",atleast2starting_trk,{vtx2_ntrk,"dsvtx.vtx2_incoming"});
    //two good tracks
    auto dfcheck_twogoodtracks_trk = dfcheck_charmdaughtersamevent.Define("dsvtx_vtx2_2goodtrks_trk",atleast2goodtrks_trk,{vtx2_ntrk,"dsvtx.vtx2_tnseg"});
    auto dfcheck_twogoodtracks = dfcheck_twogoodtracks_trk.Define("dsvtx_vtx2_2goodtrks",atleast2goodtrks,{vtx2_ntrk,"dsvtx.vtx2_tnseg"});
    //event topology
    auto dfcheck_topology = dfcheck_twogoodtracks.Define(vtx_topology,event_topology,{vtx2_charmdaugthersamevent,vtx2_positivedz});

    //counting topologies
    auto ntotal = dfcheck_topology.Filter("vtx_topology>0").Count();
    auto nsilver = dfcheck_topology.Filter("vtx_topology==1").Count();
    auto ngolden = dfcheck_topology.Filter("vtx_topology==2").Count();

    cout<<"Found: "<<*ntotal<< " primaries with good MC decay vertices: "<<*nsilver<<" single charm and "<<*ngolden<<" double charm"<<endl;
    //Saving output to file
     dfcheck_topology.Snapshot("ds","annotated_ds_data_result.root");

    //applying condition and report
    dfcheck_topology = dfcheck_topology.Define("goodevent",condition);
    auto hgoodevent =   dfcheck_topology.Define("ngoodevent", howmanyvertices,{"goodevent"}).Fill<int>(TH1I("hgood", "good vertex", 10, 0, 10), {"ngoodevent"});
    TCanvas *csame = new TCanvas();
    hgoodevent->DrawClone();
    //reporting values
    int ngood_prodvertices=hgoodevent->Integral(2,10);
    cout<<"N primary vertices with at least one good production vertex: "<<ngood_prodvertices<<endl;
    int ngood_secondaries = 0;
    for (int ibin = 1; ibin<=hgoodevent->GetNbinsX();ibin++){
      ngood_secondaries+= hgoodevent->GetBinContent(ibin) * hgoodevent->GetXaxis()->GetBinLowEdge(ibin);
    }
    cout<<"N good secondary vertices: "<<ngood_secondaries<<endl;

    //show how to plot a cut distribution here
    auto hgooddl =  dfcheck_topology.Define("dsvtx_vtx2_gooddl",selecteddistribution,{"dsvtx_vtx2_dl","dsvtx_vtx2_positivedz"}).Histo1D("dsvtx_vtx2_gooddl");
    TCanvas *c = new TCanvas();
    hgooddl->DrawClone();

    //meanlife
    auto htau = dfcheck_topology.Define("dsvtx_vtx2_meanlifecharm",selecteddistribution,{"dsvtx_vtx2_tau","goodevent"}).Histo1D("dsvtx_vtx2_meanlifecharm");
    TCanvas *ctau = new TCanvas();
    htau->DrawClone();
}
