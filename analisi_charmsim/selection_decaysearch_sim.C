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

int primary_atleast2starting(RVec<int>vtx_incoming){
  int goodvertex=0;
  int ngoodtrks=0;
    for (int itrk = 0; itrk < vtx_incoming.size(); itrk++){
      if (vtx_incoming[itrk] == 1) ngoodtrks++;
    }
    //at least 2 good tracks
    if (ngoodtrks >= 2) goodvertex=1;
    else goodvertex=0;
  return goodvertex;
}

int primary_atleast2goodtrks(RVec<int>vtx_nseg){
  int goodvertex=0;
  int ngoodtrks=0;
    for (int itrk = 0; itrk < vtx_nseg.size(); itrk++){
      if (vtx_nseg[itrk] > 2) ngoodtrks++;
    }
    //at least 2 good tracks
    if (ngoodtrks >= 2) goodvertex=1;
    else goodvertex=0;
  return goodvertex;
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

RVec<float> phi_charm ( float vx, float vy, float vz, RVec<float> vtx2_vx, RVec<float> vtx2_vy, RVec<float> vtx2_vz){
  //angle between primary and secondary vertices
  RVec<float> phi_charm;
  const int nvertices = vtx2_vx.size();
  for(int ivtx = 0; ivtx < nvertices; ivtx++){
   //computing angles
   float tx = (vtx2_vx[ivtx]-vx)/(vtx2_vz[ivtx]-vz); 
   float ty = (vtx2_vy[ivtx]-vy)/(vtx2_vz[ivtx]-vz);
   phi_charm.push_back(TMath::ATan2(ty,tx));
  }
  return phi_charm;
}

RVec<float> endingtrack_angledifference(RVec<int> vtx2_ntracks, RVec<float> vtx2_phi_charm, RVec<float> vtx2_track_tx, RVec<float>vtx2_track_ty, RVec<int> vtx2_track_incoming){
  int nprevioustracks = 0;
  const int nvertices = vtx2_ntracks.size();
  RVec<float> angledifference;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int ntracks = vtx2_ntracks[ivtx];
    int nbackward = 0; //tracks ending at vertex
    float trackphi;
    for (int itrk = 0; itrk < ntracks; itrk++){
      if (vtx2_track_incoming[itrk+nprevioustracks] == 0){
       nbackward++;
       trackphi= TMath::ATan2(vtx2_track_ty[itrk+nprevioustracks],vtx2_track_tx[itrk+nprevioustracks]);
      } 
    }
    if (nbackward == 1){
      float phidifference = TMath::Abs(vtx2_phi_charm[ivtx]- trackphi);
      if (phidifference > TMath::Pi()) phidifference = 2*TMath::Pi() - phidifference;
      angledifference.push_back(phidifference);
    } 
    else if (nbackward > 1) angledifference.push_back(10); //more than one track ending at vertex
    else if (nbackward < 1) angledifference.push_back(-10); //no tracks ending at vertex
    

   //increasing ntracks counter
    nprevioustracks += ntracks;
  }//end loop on vertices
  return angledifference;
}

//phi medium and difference with vertex phi
RVec<float> phi_medium(RVec<int> vtx2_ntracks, RVec<float> vtx2_track_tx, RVec<float> vtx2_track_ty, RVec<int> vtx2_track_incoming){
  int nprevioustracks = 0;
  const int nvertices = vtx2_ntracks.size();
  RVec<float> phi_medium_vertex;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   int ntracks = vtx2_ntracks[ivtx];
   //resetting indexes
   int nforward = 0;   
   float track_phi = 0.;

   //extracting variables for single vertex
   RVec<float> vtx2_track_tx_vertex = RVec<float>(&vtx2_track_tx[nprevioustracks], &vtx2_track_tx[nprevioustracks+ntracks]);
   RVec<float> vtx2_track_ty_vertex = RVec<float>(&vtx2_track_ty[nprevioustracks], &vtx2_track_ty[nprevioustracks+ntracks]);
   RVec<int> vtx2_track_incoming_vertex = RVec<int>(&vtx2_track_incoming[nprevioustracks], &vtx2_track_incoming[nprevioustracks+ntracks]);

   for (int itrk = 0; itrk < ntracks; itrk++){
     if (vtx2_track_incoming[itrk+nprevioustracks]==0) continue;
     nforward++;
     track_phi += TMath::ATan2(vtx2_track_ty[itrk+nprevioustracks],vtx2_track_tx[itrk+nprevioustracks]);
   }
   phi_medium_vertex.push_back(track_phi/nforward);
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
   float phidifference = TMath::Abs(phi- vtx2_phi_medium[ivtx]);

   if (phidifference > TMath::Pi()) phidifference = 2*TMath::Pi() - phidifference;
   
   phi_difference_vertex.push_back(phidifference);
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

RVec<float> all_lowmomentum(RVec<int> vtx2_ntracks, RVec<float> vtx2_trk_mc_mom){
  RVec<int> all_lowmomentum;
  const int nvertices = vtx2_ntracks.size();

  float momcut = 5.;
  int nprevioustracks = 0;
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    bool haslowmomentum = true;
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      if (vtx2_trk_mc_mom[itrk+nprevioustracks] > momcut) haslowmomentum = false;
    }
    all_lowmomentum.push_back(haslowmomentum);
    nprevioustracks += ntracks;
  }
 return all_lowmomentum;
}

RVec<float> meanlife(RVec<int>vtx2_ntracks, RVec<float> vtx2_vka, RVec<float> vtx2_dl, RVec<int> vtx2_track_incoming){
  float lightspeed = 3 * 1e+10; // cm /s
  float micron2cm = 1e-4 ;
  unsigned int nprevioustracks = 0;
  const int nvertices = vtx2_ntracks.size();
  RVec<float> vtx2_tau;
  //loop on vertices;
  //RVec<float> split_vtx2_vka = vtx2_vka; //removing tracks from each vertex
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
    int ntracks = vtx2_ntracks[ivtx];

    //resetting indexes
    int nforward = 0;
    float track_kink = 0.;

    //computing mean kink of tracks from that vertex
    for (int itrk = 0; itrk < ntracks; itrk++){
     if (vtx2_track_incoming[itrk+nprevioustracks]==0) continue;
     nforward++;
     track_kink += vtx2_vka[itrk+nprevioustracks];
   }
    float meankink = track_kink/nforward;  
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
void selection_decaysearch_sim(TString inputfilename = "ds_data_result_23_01.root"){
    //DEFINE VARIABLES CONTAINING BRANCH NAMES
    //input branches
    string vtx2_mc_pid = "dsvtx.vtx2_mc_pid";

    //vertex positions
    string vtx_x = "vtx.x";
    string vtx_y = "vtx.y";
    string vtx_z = "vtx.z";

    //track variables
    string trk_mc_ev = "trk.mc_ev";
    string trk_nseg = "trk.nseg";
    string trk_incoming = "trk.incoming";

    string vtx2_dz = "dsvtx.vtx2_dz";
    string vtx2_mc_ev = "dsvtx.vtx2_mc_ev";

    //vertex positions
    string vtx2_vx = "dsvtx.vtx2_vx";
    string vtx2_vy = "dsvtx.vtx2_vy";
    string vtx2_vz = "dsvtx.vtx2_vz";
    //vertex molteplicity
    string vtx2_ntrk = "dsvtx.vtx2_ntrk";
    //secondary, track variables
    string vtx2_tx = "dsvtx.vtx2_tx";
    string vtx2_ty = "dsvtx.vtx2_ty";
    string vtx2_incoming = "dsvtx.vtx2_incoming";
    string vtx2_vka = "dsvtx.vtx2_vka";
    string vtx2_tnseg = "dsvtx.vtx2_tnseg";
    
    //*************************newbranches********************///
    //primary, vertex variables

    string vtx_mc_ev = "vtx_mc_ev";
    string vtx_topology = "vtx_topology";
    
    string vtx_2starting = "vtx_2starting";
    string vtx_2goodtrks = "vtx_2goodtrks";

    //secondary, track variables
    string vtx2_track_2starting = "dsvtx_vtx2_2starting_trk";
    string vtx2_track_2goodtrks = "dsvtx_vtx2_2goodtrks_trk";
 
    string vtx2_track_charmdaugther = "dsvtx_vtx2_trk_charmdaughter";
    string vtx2_track_samevent = "dsvtx_vtx2_trk_samevent";
    string vtx2_track_positivedz = "dsvtx_vtx2_trk_positivedz";

    //secondary, vertex variables

    string vtx2_phicharm = "dsvtx_vtx2_phicharm";
    string vtx2_endingdeltaphi = "dsvtx_vtx2_endingdeltaphi";

    string vtx2_dl = "dsvtx_vtx2_dl";
    string vtx2_tau = "dsvtx_vtx2_tau";

    string vtx2_positivedz = "dsvtx_vtx2_positivedz";
    string vtx2_samevent = "dsvtx_vtx2_samevent";
    string vtx2_charmdaugther = "dsvtx_vtx2_charmdaughter";
    string vtx2_charmdaugthersamevent = "dsvtx_vtx2_charmdaughtersamevent"; 

    string vtx2_2starting = "dsvtx_vtx2_2starting";
    string vtx2_2goodtrks = "dsvtx_vtx2_2goodtrks";

    string vtx2_meanphi = "dsvtx_vtx2_meanphi";
    string vtx2_phidifference = "dsvtx_vtx2_phidifference";
    //condition for good charm decay vertices
    string condition = "dsvtx_vtx2_charmdaughtersamevent*dsvtx_vtx2_positivedz*dsvtx_vtx2_2goodtrks";

    //******************START OF MAIN SCRIPT*************************//

    TFile *inputfile = TFile::Open(inputfilename.Data());
    RDataFrame dsdataframe = RDataFrame("ds",inputfile);
    //computing additional variables
    auto dflength = dsdataframe.Define(vtx2_dl,decaylength,{vtx2_vx, vtx2_vy, vtx2_vz,vtx_x,vtx_y,vtx_z});
    auto dfmeanlife = dflength.Define(vtx2_tau,meanlife,{vtx2_ntrk, vtx2_vka,vtx2_dl,vtx2_incoming});
    auto df_primmcid = dfmeanlife.Define(vtx_mc_ev,mostscommonidvertex,{trk_mc_ev}); //mc event of vertex

    //angular information
    auto df_meanphi = df_primmcid.Define(vtx2_meanphi,phi_medium,{vtx2_ntrk,vtx2_tx, vtx2_ty,vtx2_incoming});
    auto df_phidifference = df_meanphi.Define(vtx2_phidifference,phi_difference,{vtx2_meanphi,vtx_x,vtx_y,vtx_z, vtx2_vx, vtx2_vy, vtx2_vz});

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
    
    //vertex angle
    auto dfphi_charm = dfcheck_charmdaughtersamevent.Define(vtx2_phicharm,phi_charm,{vtx_x,vtx_y,vtx_z,vtx2_vx, vtx2_vy, vtx2_vz});

    auto dfcheck_endingtrackdeltaphi = dfphi_charm.Define(vtx2_endingdeltaphi,endingtrack_angledifference,{vtx2_ntrk,vtx2_phicharm,vtx2_tx,vtx2_ty,vtx2_incoming});

    //two starting tracks
    auto dfcheck_2starting_trk = dfcheck_endingtrackdeltaphi.Define(vtx2_track_2starting,atleast2starting,{vtx2_ntrk,vtx2_incoming});
    auto dfcheck_2starting = dfcheck_2starting_trk.Define(vtx2_2starting,atleast2starting_trk,{vtx2_ntrk,vtx2_incoming});

    auto dfcheck_primary_2starting = dfcheck_2starting.Define(vtx_2starting,primary_atleast2starting,{"trk.incoming"});
  
    //two good tracks
    auto dfcheck_twogoodtracks_trk = dfcheck_primary_2starting.Define(vtx2_track_2goodtrks,atleast2goodtrks_trk,{vtx2_ntrk,vtx2_tnseg});
    auto dfcheck_twogoodtracks = dfcheck_twogoodtracks_trk.Define(vtx2_2goodtrks,atleast2goodtrks,{vtx2_ntrk,vtx2_tnseg});
     
    auto dfcheck_primary_twogoodtracks = dfcheck_twogoodtracks.Define(vtx_2goodtrks,primary_atleast2goodtrks,{"trk.nseg"});

    //event topology
    auto dfcheck_topology = dfcheck_primary_twogoodtracks.Define(vtx_topology,event_topology,{vtx2_charmdaugthersamevent,vtx2_positivedz});

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
    auto hgooddl =  dfcheck_topology.Define("dsvtx_vtx2_gooddl",selecteddistribution,{vtx2_dl,vtx2_positivedz}).Histo1D("dsvtx_vtx2_gooddl");
    TCanvas *c = new TCanvas();
    hgooddl->DrawClone();

    //meanlife
    auto htau = dfcheck_topology.Define("dsvtx_vtx2_meanlifecharm",selecteddistribution,{vtx2_tau,"goodevent"}).Histo1D("dsvtx_vtx2_meanlifecharm");
    TCanvas *ctau = new TCanvas();
    htau->DrawClone();
}

void addlowmomentumcondition(TString inputfilename = "annotated_ds_data_result.root"){
    TFile *inputfile = TFile::Open(inputfilename.Data());
    RDataFrame dsdataframe = RDataFrame("ds",inputfile);

    auto df_newcheck = dsdataframe.Define("dsvtx_vtx2_alllowmomentum",all_lowmomentum,{"dsvtx_vtx2_ntrk","dsvtx_vtx2_trk_mc_mom"});
    df_newcheck.Snapshot("ds","annotated_ds_data_result_2.root");
}
