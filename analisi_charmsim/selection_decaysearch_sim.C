using namespace ROOT;
//using namespace ROOT::RVec;
//dz selection
RVec<float> dzselection(vector<float> vtx2_dz){
  //converting to RVec, avoid doing loop later
  RVec<float> rvtx2_dz(vtx2_dz.data(),vtx2_dz.size());
  //applying condition
  RVec<float> rvtx2_goodIDs = (rvtx2_dz >= 0.); 
  return rvtx2_goodIDs;
}
//compuation of decay length
RVec<float> decaylength(vector<float> vtx2_vx, vector<float> vtx2_vy, vector<float> vtx2_vz, float vx, float vy, float vz){
  //converting to RVec, avoid doing loop later
  RVec<float> rvtx2_vx(vtx2_vx.data(),vtx2_vx.size());
  RVec<float> rvtx2_vy(vtx2_vy.data(),vtx2_vy.size());
  RVec<float> rvtx2_vz(vtx2_vz.data(),vtx2_vz.size());
  //computing directly, no need for loop
  RVec<float> rvtx2_dl = pow(pow(rvtx2_vx - vx,2)+ pow(rvtx2_vy-vy,2)+pow(rvtx2_vz-vz,2),0.5);

  return rvtx2_dl;
}

vector<bool> MCsamevent(vector<int>vtx2_ntracks, int vtx_mc_event, vector<int> vtx2_mc_event,  vector<int> vtx2_mcparentid){
  vector <bool> vtx2_goodIDs;
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

vector<bool> MCcharmdaughter(vector<int>vtx2_ntracks, vector<int> vtx2_mcparentid){
   vector <bool> vtx2_goodIDs;
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

int mostscommonidvertex(vector<int> trk_mc_id){
    RVec<int> rtrk_mc_id(trk_mc_id.data(),trk_mc_id.size());
    if (trk_mc_id.size() > 0) return Max(rtrk_mc_id);
    else return -1;
}
/*
RVec<bool> selectedvertices(vector<bool> selection){
    RVec <bool> rselection = new RVec<
}
RVec<bool> selectedvertices(RVec<bool> rselection){
    return rselection;
}*/

//start of script
void selection_decaysearch_sim(){
    TFile *inputfile = TFile::Open("new_bdt2_ds_data_result.root");
    RDataFrame dsdataframe = RDataFrame("ds",inputfile);
    //computing additional variables
    auto dflength = dsdataframe.Define("dsvtx_vtx2_dl",decaylength,{"dsvtx.vtx2_vx", "dsvtx.vtx2_vy", "dsvtx.vtx2_vz","vtx.x","vtx.y","vtx.z"});
    auto df_primmcid = dflength.Define("vtx_mc_ev",mostscommonidvertex,{"trk.mc_ev"});
    //checking conditions for selections
    auto dfcheck = df_primmcid.Define("dsvtx_vtx2_positivedz",dzselection,{"dsvtx.vtx2_dz"});
    auto dfcheck2 = dfcheck.Define("dsvtx_vtx2_charmdaughter", MCcharmdaughter,{"dsvtx.vtx2_ntrk","dsvtx.vtx2_mc_pid"});
    auto dfcheck3 = dfcheck2.Define("dsvtx_vtx2_samevent", MCsamevent,{"dsvtx.vtx2_ntrk","vtx_mc_ev","dsvtx.vtx2_mc_ev","dsvtx.vtx2_mc_pid"});
    
    //printout
    dfcheck3.Snapshot("ds","annotated_ds_data_result.root");

    //drawing histograms
    TFile *newfile = TFile::Open("annotated_ds_data_result.root");
    TTree *tree = (TTree*) newfile->Get("ds");

    TString cutdefault = TString("dsvtx_vtx2_positivedz");
    TString cut00 = TString("dsvtx_vtx2_positivedz&&!dsvtx_vtx2_samevent&&!dsvtx_vtx2_charmdaughter");
    TString cut10 = TString("dsvtx_vtx2_positivedz&&!dsvtx_vtx2_samevent&&dsvtx_vtx2_charmdaughter");

    tree->Draw("dsvtx_vtx2_dz",cutdefault.Data());

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
