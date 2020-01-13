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

RVec<int> MCsamevent(vector<int>vtx2_ntracks, int vtx_mc_event, vector<int> vtx2_mc_event,  vector<int> vtx2_mcparentid){
  RVec <int> vtx2_goodIDs;
   int nvertices = vtx2_ntracks.size();
   int nprevioustracks = 0;
   //loop on secondary vertices
   for (int ivtx = 0; ivtx< nvertices;ivtx++){
     int ntracks = vtx2_ntracks[ivtx];
     //loop on tracks
     bool atleastonedaughter = false;
     for (int itrk = 0; itrk < ntracks; itrk++){
         if ((vtx2_mcparentid[itrk+nprevioustracks] == 1 || vtx2_mcparentid[itrk+nprevioustracks] == 2)&& (vtx2_mc_event[itrk+nprevioustracks] == vtx_mc_event)) atleastonedaughter = true;        
     }
     if (atleastonedaughter) vtx2_goodIDs.push_back(true);
     else vtx2_goodIDs.push_back(false);
     //increasing the counter
     nprevioustracks += ntracks;
   }
   return vtx2_goodIDs;
}

RVec<int> MCcharmdaughter(vector<int>vtx2_ntracks, vector<int> vtx2_mcparentid){
   RVec<int> vtx2_goodIDs;
   int nvertices = vtx2_ntracks.size();
   int nprevioustracks = 0;
   //loop on secondary vertices
   for (int ivtx = 0; ivtx< nvertices;ivtx++){
     int ntracks = vtx2_ntracks[ivtx];
     //loop on tracks
     bool atleastonedaughter = false;
     for (int itrk = 0; itrk < ntracks; itrk++){
         if (vtx2_mcparentid[itrk+nprevioustracks] == 1 || vtx2_mcparentid[itrk+nprevioustracks] == 2) atleastonedaughter = true;        
     }
     if (atleastonedaughter) vtx2_goodIDs.push_back(true);
     else vtx2_goodIDs.push_back(false);
     //increasing the counter
     nprevioustracks += ntracks;
   }
   return vtx2_goodIDs;
}

//Track selection

RVec<int> MCcharmdaughtertrack(vector<int>vtx2_ntracks, vector<int>vtx2_mcparentid){
  int nvertices = vtx2_ntracks.size();
  RVec<int> charmdaughter_track;
  int nprevioustracks = 0;
  //loop into vertices
  for (int ivtx = 0; ivtx < nvertices; ivtx++){
   //loop into tracks
    int ntracks = vtx2_ntracks[ivtx];
    for (int itrk = 0; itrk < ntracks; itrk++){
      //adding track bool information
      if (vtx2_mcparentid[itrk+nprevioustracks] == 1 || vtx2_mcparentid[itrk+nprevioustracks] == 2) charmdaughter_track.push_back(true);    
      else charmdaughter_track.push_back(false);
    }//close track loop
   nprevioustracks += ntracks;
  }//close vertex loop
  return charmdaughter_track;
}

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
    TFile *inputfile = TFile::Open("01_13_ds_data_result.root");
    RDataFrame dsdataframe = RDataFrame("ds",inputfile);
    //computing additional variables
    auto dflength = dsdataframe.Define("dsvtx_vtx2_dl",decaylength,{"dsvtx.vtx2_vx", "dsvtx.vtx2_vy", "dsvtx.vtx2_vz","vtx.x","vtx.y","vtx.z"});
    auto df_primmcid = dflength.Define("vtx_mc_ev",mostscommonidvertex,{"trk.mc_ev"});
    //checking conditions for selections
    auto dfcheck = df_primmcid.Define("dsvtx_vtx2_positivedz",dzselection,{"dsvtx.vtx2_dz"});
    
    auto dfcheck2 = dfcheck.Define("dsvtx_vtx2_charmdaughter", MCcharmdaughter,{"dsvtx.vtx2_ntrk","dsvtx.vtx2_mc_pid"});
    auto dfcheck2_trk = dfcheck2.Define("dsvtx_vtx2_trk_charmdaughter",MCcharmdaughtertrack,{"dsvtx.vtx2_ntrk","dsvtx.vtx2_mc_pid"});

    auto dfcheck3 = dfcheck2_trk.Define("dsvtx_vtx2_samevent", MCsamevent,{"dsvtx.vtx2_ntrk","vtx_mc_ev","dsvtx.vtx2_mc_ev","dsvtx.vtx2_mc_pid"});

    auto hsameevent =  dfcheck3.Define("compositecondition","dsvtx_vtx2_samevent*dsvtx_vtx2_positivedz").Define("nsamevent", howmanyvertices,{"compositecondition"}).Fill<int>(TH1I("hgood", "good vertex", 10, 0, 10), {"nsamevent"});
    TCanvas *csame = new TCanvas();
    hsameevent->DrawClone();
    //reporting values
    int ngood_prodvertices=hsameevent->Integral(2,10);
    cout<<"N primary vertices with at least one good production vertex: "<<ngood_prodvertices<<endl;
    int ngood_secondaries = 0;
    for (int ibin = 1; ibin<=hsameevent->GetNbinsX();ibin++){
      ngood_secondaries+= hsameevent->GetBinContent(ibin) * hsameevent->GetXaxis()->GetBinLowEdge(ibin);
    }
    cout<<"N good secondary vertices: "<<ngood_secondaries<<endl;

    //show how to plot a cut distribution here
    auto hgooddl = dfcheck3.Define("dsvtx_vtx2_gooddl",selecteddistribution,{"dsvtx_vtx2_dl","dsvtx_vtx2_positivedz"}).Histo1D("dsvtx_vtx2_gooddl");
    TCanvas *c = new TCanvas();
    hgooddl->DrawClone();
    //printout
    dfcheck3.Snapshot("ds","annotated_ds_data_result.root");

    //old size computation, for now useless
       // auto df1 = dsdataframe.Define("size",[](vector<int>myvec){return myvec.size();},{"dsvtx.vtx2_vid"});
    //auto hsize = df1.Histo1D({"hsize","size of histogram",10,0,10},"size");
    //how many surpass a selection?
   // auto entries1 = df1.Filter([](vector<int>myvec){return myvec.size()>=1;},{"dsvtx.vtx2_vid"}).Count();  
    //end of 'lazy' part, give me results
    //cout<<*entries1<<endl; 
    //TCanvas *csize = new TCanvas();
    //hsize->DrawClone();
}
