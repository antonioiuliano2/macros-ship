void SimpleShowerRecInterface(){

    // Create ShowerRec Object
    EdbShowerRec * eShowerRec = new EdbShowerRec();

    // Create Initiator BT array:
    TObjArray * eInBTArray=new TObjArray();

    // Reset eShowerRec Arrays: InBTArray and RecoShowerArray....
    eShowerRec->ResetInBTArray();
    eShowerRec->ResetRecoShowerArray();

    //load first segment from track
    TFile *trackfile = TFile::Open("linked_tracks.root");
    TTree *tracks = (TTree*)trackfile->Get("tracks");

    TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
    tracks->SetBranchAddress("s",  &seg);

    tracks->GetEntry(10);
    EdbSegP *segtest = (EdbSegP*) seg->At(0);

    //set array with inBT (array of initiator base tracks)
    TFile *myfile = TFile::Open("pvrectest.root");
    EdbPVRec *pvrec_couples = (EdbPVRec*) myfile->Get("EdbPVRec");
    eInBTArray->Add(segtest);
    eShowerRec->SetInBTArray(eInBTArray);
    eShowerRec->PrintInitiatorBTs();

    //set edbpvrec
    eShowerRec->SetEdbPVRec(pvrec_couples);

    cout << " eShowerRec->SetUseAliSub(0)..." << endl;
    eShowerRec->SetUseAliSub(0);

    cout << " eShowerRec->Execute()..." << endl;

    //Start actual reconstruction
    eShowerRec->Execute();

    //Print output
    eShowerRec->PrintRecoShowerArray();
}

void drawShower(int ishower=0){
    TFile *myfile = TFile::Open("pvrectest.root");
    EdbPVRec *pvrec_couples = (EdbPVRec*) myfile->Get("EdbPVRec");
    //opening showerfile
    TFile *showerfile = TFile::Open("Shower.root");
    TTree *showertree = showerfile->Get("treebranch");

    int sizeb; 
    const int maxsize = 10000; //as in ShowerRec
    int idb[maxsize]; //IDs of basetracks
    int plateb[maxsize]; //number of plate of base track
    //setting branch addresses
    showertree->SetBranchAddress("sizeb",&sizeb);
    showertree->SetBranchAddress("idb",&idb);
    showertree->SetBranchAddress("plateb",&plateb);

    TObjArray *sarr = new TObjArray();
    //filling array with segments
    showertree->GetEntry(ishower);
    for (int iseg = 0; iseg < sizeb; iseg++){
        //getting edbseg (need to apply affine transformations)
        int vid = pvrec_couples->Vid( plateb, idb);// vid =  pid*1000000+sid, identify a segment in a unique way
        sarr->Add(pvrec_couples->GetSegment(vid));
    }

     //DISPLAY OF SEGMENTS
    const char *dsname = "Test shower reconstruction";
    ds = EdbDisplay::EdbDisplayExist(dsname);
    if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-40000., 0..);
    ds->SetVerRec(gEVR);
    ds->SetDrawTracks(4);
    ds->SetArrSegP( sarr );
    ds->Draw();

}