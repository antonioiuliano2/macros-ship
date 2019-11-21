//testing Decoding Method for HPT and TT (Created on 15 November 2019)

void test_hpt_decoder(){
 //opening file and initializing tree
 TFile *file = TFile::Open("ship.conical.PG_13-TGeant4.root"); 
 if (!file) return;
 TTreeReader reader("cbmsim",file);   
 TTreeReaderArray<HptPoint> hpthits(reader,"HptPoint");
 TTreeReaderArray<TTPoint> tthits(reader,"TTPoint");
 //initializing objects for helper decode methods (it may be useful to use static methods in the future)
 Hpt *testhpt = new Hpt();
 TargetTracker *testTT = new TargetTracker();

 const int nentries = reader.GetEntries();

 //define lambda function with operation to do for each point
 auto dumpPoint = [testhpt, testTT](const FairMCPoint &hit){
         //get original detectorID
         int trackID = hit.GetTrackID();
         int detID = hit.GetDetectorID();
         //Reset indexes
         int nHPT = 0;
         int nplane = 0;
         bool ishor = kFALSE;
         //Ask where the hit is from, then apply corresponding decoder and print out results               
         if (strncmp(hit.GetName(),"HptPoint",8)) testhpt->DecodeVolumeID(detID,nHPT, nplane,ishor); //HPT decoder method
         else testTT->DecodeTTID(detID,nHPT, nplane,ishor);
         cout<<trackID<<"\t"<<detID<<"\t"<<nHPT<<"\t"<<nplane<<"\t"<<ishor<<endl;
 };


 //starting loop on events
 for (int i = 0; i < 10; i++){
     
     cout<<"Start dump event: "<<i<<endl;
     reader.SetEntry(i);
     cout<<"Start loop HPT"<<endl;
     cout<<"TrackID\tDetID\tnHPT\tnplane\tishor"<<endl;    
     for (const HptPoint& hit:hpthits){      
         dumpPoint(hit);
        }
    
     cout<<"Start loop TT"<<endl;
     cout<<"TrackID\tDetID\tnTT\tnplane\tishor"<<endl;
     for (const TTPoint& hit:tthits){
         dumpPoint(hit);
        }
    }
 //deleting helper objects
 delete testhpt;
 delete testTT;
}
