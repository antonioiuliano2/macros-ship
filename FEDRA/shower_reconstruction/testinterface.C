//Small code to test my new libShower interface class
//do .L SimpleShowerRecInterface.C (online loading, for now no makefile for compilation)
//   .L testinterface.C
//Then of the following:
//buildcouples() for saving a new pvr
//testinterface() for shower reconstruction
//drawfoundshowers() for shower drawing
#define gEDBDEBUGLEVEL 2

using namespace ROOT;
void testinterface(){
    SimpleShowerRecInterface testinterface = SimpleShowerRecInterface();

    
    TFile *myfile = TFile::Open("pvrec_testdata2.root");
    testinterface.LoadPVRec(myfile);

    const int ntracks = 3;
    int tracklist[ntracks] = {1068,477799,10040};
  
    int iseglist[ntracks] = {1,0,0};
    testinterface.RecoFromTrack(ntracks,tracklist,iseglist,"b000005.0.0.0.trk_no27and20.root");
    
    
    //simulation
    //TFile *myfile = TFile::Open("pvrec_90000_40000.root");
    //  const int trackID =  133795;
    //testinterface.RecoFromTrack(ntracks,tracklist,"linked_tracks_simtest.root");
}

void drawfoundshowers(){
    SimpleShowerRecInterface testinterface = SimpleShowerRecInterface();
    int icell = 1;
    int startx = (icell / 5);
    int starty = (icell - (5*startx))* 10000;
    startx = (startx + 6)*10000;
    cout<<startx<<" "<<starty<<endl;
    TFile *pvrecfile = TFile::Open(Form("pvrecs/pvrec_%d_%d.root",startx,starty));
    testinterface.LoadPVRec(pvrecfile);
    testinterface.DrawAllShowers();
    //testinterface.DrawShower(0);

}

void buildcouples(){
    SimpleShowerRecInterface testinterface = SimpleShowerRecInterface();
    const int ibrick = 5;
    const int nplates = 29;
    //need to give list of PIDs, due to missing/hidden plates (-1 if missing)
    int PID[nplates] = {24,23,22,21,20,-1,19,-1,18,17,16,15,14,13,12,11,10,9,8,-1,7,6,5,4,3,2,-1,1,0};
    
    //selection (can be different for each plate)
    TString selections[nplates];
    
    for (int iplate = 0; iplate<nplates;iplate++) selections[iplate]=TString("TMath::Abs(s.eTX)<0.2 && TMath::Abs(s.eTY)<0.2");
    
    testinterface.BuildPVRec(ibrick, nplates, PID,selections,"pvrec_testdata2.root");
}


void prepare_simlists(int ncell = 1){

  //getting tracks from Valerio's file
  TFile *dsfile = TFile::Open("annotated_ds_data_result.root");
  TTreeReader dsreader("ds",dsfile);

  TTreeReaderArray<int> trackIDs(dsreader,"dsvtx_vtx2_tid");
  TTreeReaderArray<float> xtrack(dsreader,"dsvtx_vtx2_xt");
  TTreeReaderArray<float> ytrack(dsreader,"dsvtx_vtx2_yt");

  const int npvrecs = 30; //number of built sections
  
  RVec<int> trackscells[npvrecs];

  const int nevents = dsreader.GetEntries(); //events in tree
  //int ncell = 4; //which cell to apply reconstruction over

  for (int ievent = 0; ievent < nevents; ievent++){
      dsreader.SetEntry(ievent);
      int ntracks = trackIDs.GetSize();
      for (int itrk = 0; itrk < ntracks; itrk++){
        int xcode = (xtrack[itrk]-60000)/10000; //second quarter x from 60000 to 120000
        int ycode = (ytrack[itrk])/10000; //second quarter: y from 0 to 50000
        int code = 5*xcode + ycode; //find the right file
        if (code == ncell) trackscells[code].push_back(trackIDs[itrk]);
      }//end track loop
  }//end event loop
  //checking associations
  SimpleShowerRecInterface myinterface;
  int startx, starty;
  for (int ipvrec = 0; ipvrec < npvrecs; ipvrec++){
      if (ipvrec != ncell) continue;      
      cout<<"Found in cell "<<ipvrec<<" "<<trackscells[ipvrec].size()<<" tracks "<<endl;
      startx = (ipvrec / 5);
      starty = (ipvrec - (5*startx))* 10000;
      startx = (startx + 6)*10000;

      TFile *pvrecfile = TFile::Open(Form("pvrecs/pvrec_%d_%d.root",startx,starty));
      myinterface.LoadPVRec(pvrecfile);

      RVec<int> seglist = trackscells[ipvrec]<0; // all 0
      myinterface.RecoFromTrack(trackscells[ipvrec].size(),trackscells[ipvrec].data(),seglist.data());      
      cout<<"Predicted start and end "<<startx<<" "<<starty<<" cell number "<<ncell<<" it should be "<<ipvrec<<endl;
  }

  //renaming trees and moving them in output folder
  gSystem->Exec(Form("mv Shower.root Showreco_ds_23_01/Shower_%i_%i.root",startx,starty));
  gSystem->Exec(Form("mv Shower2.root Showreco_ds_23_01/Shower2_%i_%i.root",startx,starty));
  gSystem->Exec(Form("mv shower1.root Showreco_ds_23_01/shower1_%i_%i.root",startx,starty));



}