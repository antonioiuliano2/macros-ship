void TrackDisplay(){
 int trackID = 1487;

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

 tracktree->GetEntryWithIndex(trackID);
 //temporary object for reading the file and building EdbTrackP
 EdbTrackP * temptrack = new EdbTrackP();
 temptrack->Copy(EdbTrackP(trk));
 cout<<nseg<<endl;
 //start loop on segments associated to the track
 for (int i = 0; i< nseg; i++){ 
     EdbSegP *seg = segments->At(i);
     EdbSegP *segf = segments->At(i);
     seg->SetDZ(300);
     segf->SetDZ(300);
     temptrack->AddSegment(seg);
     temptrack->AddSegmentF(segf);
     temptrack->SetSegmentsTrack(temptrack->ID()); //track segments association
     temptrack->SetCounters();
 }
 trackstodraw->Add(temptrack); //add to list of tracks to draw

gStyle->SetPalette(1);
EdbDisplay * ds = new EdbDisplay("Charm simulation FEDRA display",-50000.,50000.,-50000.,50000.,-4000.,80000.);
ds->SetDrawTracks(4); //option to draw tracks with segments
ds->SetArrTr(trackstodraw );
ds->Draw();


}



