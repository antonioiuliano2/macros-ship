void SimpleShowerRecInterface(){

    // Create ShowerRec Object
    EdbShowerRec * eShowerRec = new EdbShowerRec();

    // Create Initiator BT array:
    TObjArray * eInBTArray=new TObjArray();

    // Reset eShowerRec Arrays: InBTArray and RecoShowerArray....
    eShowerRec->ResetInBTArray();
    eShowerRec->ResetRecoShowerArray();

    //load track
    TFile *trackfile = TFile::Open("linked_tracks.root");
    TTree *tracks = (TTree*)trackfile->Get("tracks");

    //set array with inBT (array of initiator base tracks)

    eShowerRec->SetInBTArray(arr);
    eShowerRec->PrintInitiatorBTs();

    //set edbpvrec
    TFile *myfile = TFile::Open(Form("pvrecs/pvrec_%d_%d.root",xstart,ystart));
    EdbPVRec *pvrec_couples = (EdbPVRec*) myfile->Get("EdbPVRec");
    eShowerRec->SetEdbPVRec(pvrec_couples);

    //Start actual reconstruction
    eShowerRec->Execute();

    //Print output
    eShowerRec->PrintRecoShowerArray();
}