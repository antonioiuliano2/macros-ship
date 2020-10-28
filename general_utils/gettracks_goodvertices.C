void gettracks_goodvertices(int quarter){
 //requirements for a 'good vertex'
 //int min_ntracks = 6;
 //float min_bdtvalue = 0.15;
 const float maxbdtratio = 0.05; // (1-bdt)/ntracks < maxbdtratio 
 
 TString quarterfolders[4] = {"firstquarter","secondquarter","thirdquarter","fourthquarter"};

 //opening file and getting trees
 TFile *vertexfile = TFile::Open((quarterfolders[quarter-1]+"/vertextree_"+quarterfolders[quarter-1]+".root").Data());
 TFile *bdtfile = TFile::Open("BDT_eval/vtx_BDT_evaluated_DATA_CH2R4.root");
 
 //setting bdt tree
 TTree *bdttree = (TTree*) bdtfile->Get("bdt");

 bdttree->BuildIndex("bdt_vID"); //I need to know bdt value for a given vID
 
 float bdtvalue;
 bdttree->SetBranchAddress("bdt_value",&bdtvalue);
 //setting vtx tree
 TTreeReader vtxreader("vtx",vertexfile);
 TTreeReaderValue<int> vID(vtxreader,"vID");
 TTreeReaderValue<int> ntracks(vtxreader,"n");
 
 TTreeReaderArray<int> TrackIDs(vtxreader,"TrackID");
 TTreeReaderArray<int> nsegments(vtxreader,"nseg");
 TTreeReaderArray<int> incoming(vtxreader,"incoming"); //0 if track ends at vertex, 1 if track starts at vertex
 
 //filling a tree with ids of good tracks
 TFile *outputfile = new TFile(Form("indexestracks_goodvertices_%i.root",quarter),"NEW");
 TTree *outputtree = new TTree("goodtrks","Tracks belonging to good vertices");

 int trackID, nseg;
 int vtx_ID, vtx_ntrks;

 outputtree->Branch("trid",&trackID,"trid/I");

 outputtree->Branch("nseg",&nseg,"nseg/I");
 outputtree->Branch("quarter",&quarter,"quarter/I");

 outputtree->Branch("vtx_vID",&vtx_ID,"vtx_vID/I");
 outputtree->Branch("vtx_ntrks",&vtx_ntrks,"vtx_ntrks/I");
 outputtree->Branch("vtx_bdt_value",&bdtvalue,"vtx_bdt_value/F");

 const int nentries = vtxreader.GetEntries();
 cout<<"Starting loop over vertices "<<nentries<<endl;
 for (int ivtx = 0; ivtx < nentries; ivtx++){
     //getting entry for that vertex
     if (ivtx%10000 == 0) cout<<" arrived at vertex "<<ivtx<<endl;
     vtxreader.SetEntry(ivtx);     
     vtx_ID = *vID;
     vtx_ntrks = *ntracks;

     //getting bdt value
     int bdt_entrynumber = bdttree->GetEntryNumberWithIndex(vtx_ID);
     if (bdt_entrynumber <0) continue;
     bdttree->GetEntry(bdt_entrynumber);
     //applying condition to select good vertex
     float bdtratio = (1 - bdtvalue)/vtx_ntrks;
     if (bdtratio < maxbdtratio){
         //looping over all tracks and saving tracks starting at that vertex
         for (int itrk=0; itrk < vtx_ntrks; itrk++){
            if (incoming[itrk]==1){
                trackID = TrackIDs[itrk];
                nseg = nsegments[itrk];
                outputtree->Fill();
            }
         }//end loop over tracks
     }
 } //end loop over vertices
 outputfile->cd();
 outputtree->Write();
 outputfile->Close();
} //end main program
