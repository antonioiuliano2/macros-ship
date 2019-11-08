/*usage:
 
 scp /Users/Giuliana/Desktop/PhD/FOOT/TEST_RECO_MIC/TrackDisplay.C scanner@nusrv9.na.infn.it:/home/scanner/foot/RECO_MC_TEST_giul/b000002
 
 */

void TrackDisplay(){
    
    // int trackID=0;
    TFile *inputfile = TFile::Open("linked_tracks.root","READ");
    TTree *tracktree = (TTree*) inputfile->Get("tracks");
    //tracktree.SetAlias("trk","t.") #points create confusion to python
    
    tracktree->BuildIndex("trid"); //to read entry with given trackid
    TObjArray * trackstodraw = new TObjArray(100);
    //setting branches
    
    int nseg=0;
    EdbSegP *trk;
    TClonesArray *segments = new TClonesArray("EdbSegP",100);
    TClonesArray *fittedsegments = new TClonesArray("EdbSegP",100);
    
    tracktree->SetBranchAddress("nseg",&nseg);
    tracktree->SetBranchAddress("t.",&trk);
    tracktree->SetBranchAddress("s",&segments);
    tracktree->SetBranchAddress("sf",&fittedsegments);

    for(int trackID=1; trackID<5; trackID++){
        segments->Clear();
        fittedsegments->Clear();
        tracktree->GetEntryWithIndex(trackID);
        //temporary object for reading the file and building EdbTrackP
        EdbTrackP * temptrack = new EdbTrackP();
        ((EdbSegP*)temptrack)->Copy(*trk);
       // temptrack->Copy(EdbTrackP(trk));
        cout<<"Track: " << trackID << " with " << nseg << " segments" <<endl;
        //start loop on segments associated to the track
        for (int i = 0; i< nseg; i++){
            temptrack->SetID(trk->ID());
            EdbSegP *seg = (EdbSegP*) segments->At(i);
            EdbSegP *segf = (EdbSegP*) fittedsegments->At(i);
            seg->SetDZ(300);
            segf->SetDZ(300);
            temptrack->AddSegmentF(new EdbSegP(*((EdbSegP*)(seg))));
	        temptrack->AddSegmentF(new EdbSegP(*((EdbSegP*)(segf))));
            temptrack->SetSegmentsTrack(temptrack->ID()); //track segments association
            temptrack->SetCounters();
        }    
	trackstodraw->Add(temptrack);
    }
    gStyle->SetPalette(1);

    const float zmin = -100000.;
    const float zmax = 0.;

    EdbDisplay * ds = new EdbDisplay("FOOT simulation FEDRA display",-50000.,50000.,-50000.,50000.,zmin,zmax);
    ds->SetDrawTracks(4); //option to draw tracks with segments
    ds->SetArrTr(trackstodraw);
    ds->Draw();
}



