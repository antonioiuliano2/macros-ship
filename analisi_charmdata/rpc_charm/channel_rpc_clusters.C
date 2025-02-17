//adding a channel branch to the clusters
void channel_rpc_clusters(TString subpath){
 TString eosinputpath = TString("/eos/experiment/ship/data/rpc_charm/");
 TString afsoutputpath = TString("/afs/cern.ch/work/a/aiuliano/public/Charmdata/rpc_charm/");

 //TString subpath = TString("CHARM1/RPC_RecoTracks_run2793_s1f22afc8.root");


 TFile *RPCmap = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/Charmdata/rpc_charm/RPCchannelmapfunctions.root","read");
 TF1 *mapV = (TF1*) RPCmap->Get("fV"); //sono uguali, per ora prendo il primo
 TF1 *mapH = (TF1*) RPCmap->Get("fH");

 ROOT::RDataFrame df("RPC_RecoTracks",(eosinputpath+subpath).Data());

 auto findchannel = [mapV, mapH](ROOT::VecOps::RVec<int> &dir, ROOT::VecOps::RVec<float> &clx, ROOT::VecOps::RVec<float> &cly){
  vector<int> ichannel;
  int arraysize = clx.size();
  int ibin = 1000;
  //getting number of channel from cluster position
  for (int i = 0; i < arraysize; i++){
   if (dir[i] == 0){ 
       ibin = round(mapH->GetX(cly[i]))-1; //find the number of the channel for that value
       ichannel.push_back(ibin);
   }
   else{ 
       ibin = round(mapV->GetX(clx[i]))-1;
       ichannel.push_back(ibin);   
   }
  }
  return ichannel;
 };

 auto df1 = df.Define("cl_channel",findchannel,{"cl_dir","cl_x","cl_y"});

 // as a safety, I do not want to overwrite existing files. Better safe than sorry!
 ROOT::RDF::RSnapshotOptions opts;
 opts.fMode = "NEW";
 df1.Snapshot("RPC_RecoTracks",(afsoutputpath+subpath).Data());
 cout<<"Finished writing the tree with additional branch for channel"<<endl;
}


//creating the histograms with the maps

void createhistograms(){

 TFile * RPCmap = new TFile("/afs/cern.ch/work/a/aiuliano/public/Charmdata/rpc_charm/RPCchannelmapfunctions.root","RECREATE");

 fstream inputfile;
 inputfile.open("RPC_strip_coord_CernTB.dat");

 //variables to read from the file
 char dir;
 float channelx;
 float channely;
 float channelz;
 int rpc;
 int idchannel;
//information about number of lines to read
 const int nstations = 5;
 const int nchannelsV = 184; 
 const int nchannelsH = 116;

 TH1I *hchannelmapV [nstations];
 TH1I *hchannelmapH [nstations];

 for (int istation = 0; istation < nstations; istation++){
     hchannelmapV[istation] = new TH1I(TString::Format("hchannelmapV[%i]",istation),TString::Format("V channel map for RPC %i",istation),nchannelsV,1,nchannelsV+1);
     hchannelmapH[istation] = new TH1I(TString::Format("hchannelmapH[%i]",istation),TString::Format("H channel map for RPC %i",istation),nchannelsH,1,nchannelsH+1);
 }

 //main loop
 for (int i = 0; i < nstations; i++){
  inputfile >> rpc >> channelx >> channely >> channelz;   
  //first loop, vertical strips
  for (int jV = 1; jV <= nchannelsV; jV++){
      inputfile >> dir >> idchannel >> channelx >> channely >> channelz;
      if (jV < 6) cout<<"PROVA lettura file :"<<i<<" "<<dir<<" "<<idchannel<<" "<<channelx<<" "<<channely<<" "<<channelz<<endl;
      hchannelmapV[i]->SetBinContent(idchannel,channelx);
   }
  //second loop, horizontal strips
  for (int jH = 1; jH <= nchannelsH; jH++){
      inputfile >> dir >> idchannel >> channelx >> channely >> channelz;
      hchannelmapH[i]->SetBinContent(idchannel,channely);
   }
 

 }//end of loop over the stations

 TF1 *fV = new TF1("fV","[0] *x + [1]",0,nchannelsV);
 TF1 *fH = new TF1("fH","[0] *x + [1]",0,nchannelsH);
 //draw the histograms
 TCanvas *c[5];
 for (int istation = 0; istation < nstations; istation++){
  c[istation] = new TCanvas();
  c[istation]->Divide(1,2);

  c[istation]->cd(1);
  hchannelmapV[istation]->Draw();
  hchannelmapV[istation]->GetXaxis()->SetTitle("nchannel");
  hchannelmapV[istation]->GetYaxis()->SetTitle("x[cm]");

  c[istation]->cd(2);
  hchannelmapH[istation]->Draw();
  hchannelmapH[istation]->GetXaxis()->SetTitle("nchannel");
  hchannelmapH[istation]->GetYaxis()->SetTitle("y[cm]");
  
  hchannelmapH[istation]->Write();
  hchannelmapV[istation]->Write();

 } 

 hchannelmapV[0]->Fit(fV);
 hchannelmapH[0]->Fit(fH);

 fV->Write();
 fH->Write();

// inputfile.close();
// RPCmap->Close();
}
