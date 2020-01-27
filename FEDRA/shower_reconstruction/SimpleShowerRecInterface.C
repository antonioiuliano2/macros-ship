#define SimpleShowerRecInterface_C
#include "SimpleShowerRecInterface.h"


void SimpleShowerRecInterface::LoadPVRec(TFile *inputfile){
    eEdbPVRec = (EdbPVRec*) inputfile->Get("EdbPVRec");
};

void SimpleShowerRecInterface::BuildPVRec(const int ibrick, const int nplates, int* PID,  TString *selections, char * outputfilename){ //build a new PVRec according to given selection, write it to file

    EdbDataProc *dproc = new EdbDataProc();

    TFile *outputfile = new TFile(outputfilename,"RECREATE");
    eEdbPVRec = new EdbPVRec();
    //removing tracks eventually used to set patterns, container of tracks should be void;
    eEdbPVRec->eTracks = 0;
    loadcouples(ibrick, nplates, PID, selections);

    //save obtained eEdbPVRec to external files
    outputfile->cd();
    eEdbPVRec->Write();
    outputfile->Close();
    
    
    //dproc->InitVolume(100,"nseg>1 && trid<20000"); //reading some tracks to set the patterns, not done due to absence of long tracks
    // eEdbPVRec = dproc->PVR();
    

}

void SimpleShowerRecInterface::loadcouples(const int ibrick, const int nplates, int* PID, TString *selections){ //load couples from cp files, apply affine transformations    
    TString runpath = TString(Form("/ship/DESY2019/RUN%i/b00000%i/",ibrick,ibrick));
    //opening setfile (informations about transformations)
    TFile *setfile = TFile::Open((runpath+TString(Form("b00000%i.0.0.0.set.root",ibrick)).Data()));
    EdbScanSet *set = (EdbScanSet*) setfile->Get("set");

    EdbCouplesTree *ect[nplates];
    //reading couples, starting loop within plate
    for (int iplate = nplates; iplate >= 1; iplate--){
      if(PID[iplate-1]<0) continue; // id of missing/hidden plates

     EdbPlateP *p = set->GetPlate(iplate);
     float zplate = p->Z();
     EdbAffine2D *aff = new EdbAffine2D();
     set->GetAffP2P(iplate, nplates, *aff); //usually last plate is the reference one

     //setting tree
     ect[iplate-1] = new EdbCouplesTree();
     //string spaghetti code
     if (iplate <10) ect[iplate-1]->InitCouplesTree("couples",(runpath+TString(Form("p00%i/%i.%i.0.0.cp.root",iplate,ibrick,iplate))).Data(),"READ");
     else ect[iplate-1]->InitCouplesTree("couples",(runpath+TString(Form("p0%i/%i.%i.0.0.cp.root",iplate,ibrick,iplate))).Data(),"READ");

     //loop into couples (only the ones passing condition)
     cout<<"Total of couples "<<ect[iplate-1]->eTree->GetEntries() <<endl;
     ect[iplate-1]->eTree->Draw(">>goodcouples", selections[iplate-1].Data());
     TEventList *goodcouples = (TEventList*) gDirectory->GetList()->FindObject("goodcouples");
  
     const int ngoodcouples = ect[iplate-1]->eTree->GetEntries(selections[iplate-1].Data());

     //int nseg = ect[i-1]->eTree->GetEntries();
     cout<<"Reading "<<ngoodcouples<<" from plate "<<iplate<<" at z position "<<zplate<<" with selection "<<selections[iplate-1].Data()<<endl;

     //printing transformation to check that they are correct
     cout<<"zposition "<<zplate<<endl;
     cout<<"Affine transformation: "<<endl;
     aff->Print();
     

     //*************LOOP OVER SEGMENTS***************
     for (int iseg = 0; iseg< ngoodcouples; iseg++){
      //getting entry of good segment;
      int igoodsegment = goodcouples->GetEntry(iseg);
      //***Getting information about that segment***;
      ect[iplate-1]->GetEntry(igoodsegment);
      EdbSegP *seg = new EdbSegP();
      seg->Copy(*(ect[iplate-1]->eS));
      //setting z and affine transformation
      seg->SetZ(zplate);
      seg->Transform(aff);
      //setting plate information
      seg->SetPlate(iplate);
      seg->SetPID(PID[iplate-1]);
      eEdbPVRec->AddSegment(*seg);
      }
      cout<<endl;
     }
     //setting pattern IDs
     for (int ipattern = 0; ipattern < eEdbPVRec->Npatterns(); ipattern++){
      EdbSegP * onesegment = eEdbPVRec->GetPattern(ipattern)->GetSegment(0);
      eEdbPVRec->GetPattern(ipattern)->SetPID(onesegment->PID());	
     }
}

void SimpleShowerRecInterface::RecoFromTrack(int ntracks, int* tracklist, int* iseglist,char * filename){

    // Create ShowerRec Object
    EdbShowerRec * eShowerRec = new EdbShowerRec();

    
    // Print parameters
    eShowerRec->PrintParameters();
    
    // Create Initiator BT array:
    TObjArray * eInBTArray=new TObjArray();

    // Reset eShowerRec Arrays: InBTArray and RecoShowerArray....
    eShowerRec->ResetInBTArray();
    eShowerRec->ResetRecoShowerArray();

    //load first segment from track
    TFile *trackfile = TFile::Open(filename);
    TTree *tracks = (TTree*)trackfile->Get("tracks");

    TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
    tracks->SetBranchAddress("s",  &seg);


    for (int itrk = 0; itrk < ntracks; itrk++){
     tracks->GetEntry(tracklist[itrk]);
     
     if (!(seg->At(iseglist[itrk]))) cout<<"ERROR: Missing segment "<<iseglist[itrk]<<" for track "<<tracklist[itrk]<<endl;
     EdbSegP *segtest = new EdbSegP(*((EdbSegP*) seg->At(iseglist[itrk])));
     
     //set array with inBT (array of initiator base tracks)
     eInBTArray->Add(segtest);
    }
    eShowerRec->SetInBTArray(eInBTArray);
    eShowerRec->PrintInitiatorBTs();

    //set edbpvrec
    eShowerRec->SetEdbPVRec(eEdbPVRec);

    cout << " eShowerRec->SetUseAliSub(0)..." << endl;
    eShowerRec->SetUseAliSub(0);

    cout << " eShowerRec->Execute()..." << endl;

    //Start actual reconstruction
    eShowerRec->Execute();

    //Print output
    eShowerRec->PrintRecoShowerArray();

    const char *dsname = "Test shower reconstruction";

    TObjArray * segarray = new TObjArray();
    EdbTrackP *tr1 = new EdbTrackP();
   
    if (eShowerRec->GetRecoShowerArray( )->GetEntries()>0){
    for (int ishower; ishower < eShowerRec->GetRecoShowerArray( )->GetEntries();ishower++){
     EdbTrackP *mytrack = (EdbTrackP*) (eShowerRec->GetRecoShowerArray( )->At(ishower));
     ((EdbSegP*)tr1)->Copy(*mytrack);
      for(int i=0; i<mytrack->N(); i++) {
       EdbSegP *testseg = new EdbSegP();
       testseg->Copy(*((EdbSegP*)(mytrack->GetSegment(i))) );
       testseg->SetDZ(300.);
       tr1->AddSegment(new EdbSegP(*testseg ));
       tr1->AddSegmentF(new EdbSegP(*testseg ));
      }
     segarray->Add(tr1);
     }
     /*
     EdbDisplay * ds = EdbDisplay::EdbDisplayExist(dsname);
     if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-100000., 0.);

     EdbVertexRec *gEVR = new EdbVertexRec();
     ds->SetVerRec(gEVR);
     ds->SetDrawTracks(6);
     ds->SetArrTr( segarray  );
     ds->Draw();*/
     }   
}

void SimpleShowerRecInterface::DrawShower(int ishower,char * showerfilename){

    //opening showerfile
    TFile *showerfile = TFile::Open(showerfilename);
    TTree *showertree = (TTree*) showerfile->Get("treebranch");

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
    cout<<sizeb<<endl;

    TClonesArray *segments[29];
    for (int iplate = 0; iplate < 29; iplate ++){
        if (eEdbPVRec->GetPatternByPID(iplate)) segments[iplate] = eEdbPVRec->GetPatternByPID(iplate)->GetSegments();
    }

    for (int iseg = 0; iseg < sizeb; iseg++){
        EdbSegP *myseg = new EdbSegP();
        //loop on segments, find ours
        for (int i = 0; i < segments[plateb[iseg]]->GetEntriesFast();i++){
            myseg->Copy(*((EdbSegP*) segments[plateb[iseg]]->At(i)));
            myseg->SetDZ(300.);
            if (myseg->ID() == idb[iseg]){ 
		        sarr->Add(new EdbSegP(*myseg));
    	    }
        }
//        
    }
     //DISPLAY OF SEGMENTS
     
    const char *dsname = "Test shower reconstruction";
    EdbDisplay * ds = EdbDisplay::EdbDisplayExist(dsname);
    if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-100000., 0.);

    EdbVertexRec *gEVR = new EdbVertexRec();
    ds->SetVerRec(gEVR);
    ds->SetDrawTracks(6);
    ds->SetArrSegP( sarr );
    ds->Draw();    

}


void SimpleShowerRecInterface::DrawAllShowers(char * showerfilename){

    //opening showerfile
    TFile *showerfile = TFile::Open(showerfilename);
    TTree *showertree = (TTree*) showerfile->Get("treebranch");

    int sizeb; 
    const int maxsize = 10000; //as in ShowerRec
    int idb[maxsize]; //IDs of basetracks
    int plateb[maxsize]; //number of plate of base track
    //setting branch addresses
    showertree->SetBranchAddress("sizeb",&sizeb);
    showertree->SetBranchAddress("idb",&idb);
    showertree->SetBranchAddress("plateb",&plateb);


    TClonesArray *segments[29];
    for (int iplate = 0; iplate < 29; iplate ++){
        if (eEdbPVRec->GetPatternByPID(iplate)) segments[iplate] = eEdbPVRec->GetPatternByPID(iplate)->GetSegments();
    }

    TObjArray *sarr = new TObjArray();
    //filling array with segments
    cout<<"Number of showers "<<showertree->GetEntries()<<endl;
    for (int ishower; ishower < showertree->GetEntries();ishower++){
     showertree->GetEntry(ishower);
     cout<<"Number of found segments "<<sizeb<<endl;

     for (int iseg = 0; iseg < sizeb; iseg++){
        EdbSegP *myseg = new EdbSegP();
        //loop on segments, find ours
        for (int i = 0; i < segments[plateb[iseg]]->GetEntriesFast();i++){
            myseg->Copy(*((EdbSegP*) segments[plateb[iseg]]->At(i)));
            myseg->SetDZ(300.);
            if (myseg->ID() == idb[iseg]){ 
		        sarr->Add(new EdbSegP(*myseg));
    	    }
        }
//        
     }
    }//end loop on showers
     //DISPLAY OF SEGMENTS
    const char *dsname = "Test shower reconstruction";
    EdbDisplay * ds = EdbDisplay::EdbDisplayExist(dsname);
    if(!ds)  ds=new EdbDisplay(dsname,-100000.,100000.,-100000.,100000.,-100000., 0.);

    EdbVertexRec *gEVR = new EdbVertexRec();
    ds->SetVerRec(gEVR);
    ds->SetDrawTracks(6);
    ds->SetArrSegP( sarr );
    ds->Draw();    

}


ClassImp(SimpleShowerRecInterface)
