//Small code to test my new libShower interface class
//do .L SimpleShowerRecInterface.C (online loading, for now no makefile for compilation)
//   .L testinterface.C
//Then of the following:
//buildcouples() for saving a new pvr
//testinterface() for shower reconstruction
//drawfoundshowers() for shower drawing
#define gEDBDEBUGLEVEL 2
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
    TFile *myfile = TFile::Open("pvrec_testdata2.root");
    testinterface.LoadPVRec(myfile);
    testinterface.DrawAllShowers();
   // testinterface.DrawShower(0);

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
