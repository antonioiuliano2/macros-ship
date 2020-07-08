void gettracks_goodvertices(){
 //requirements for a 'good vertex'
 int min_ntracks = 6;
 float min_bdtvalue = 0.15;
 int quarter = 2; //which quarter files I am reading

 //opening file and getting trees
 TFile *vertexfile = TFile::Open("secondquarter/vertextree_secondquarter.root");
 TFile *bdtfile = TFile::Open("secondquarter/vtx_BDT_data_evaluated_2nd.root");
 
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
 TFile *outputfile = new TFile("indexestracks_goodvertices.root","UPDATE");
 TTree *outputtree = new TTree(Form("goodtrks_%i",quarter),Form("Tracks belonging to good vertices from quarter_%i",quarter));

 int trackID, nseg;
 int vtx_ID, vtx_ntrks;

 outputtree->Branch("trid",&trackID,"trid/I");

 outputtree->Branch("nseg",&nseg,"nseg/I");
 outputtree->Branch("quarter",&quarter,"quarter/I");

 outputtree->Branch("vtx_vID",&vtx_ID,"vtx_vID/I");
 outputtree->Branch("vtx_ntrks",&vtx_ntrks,"vtx_ntrks/I");
 outputtree->Branch("vtx_bdt_value",&bdtvalue,"vtx_bdt_value/F");

 const int nentries = vtxreader.GetEntries();
 for (int ivtx = 0; ivtx < nentries; ivtx++){
     //getting entry for that vertex
     vtxreader.SetEntry(ivtx);     
     vtx_ID = *vID;
     vtx_ntrks = *ntracks;

     //getting bdt value
     bdt_entrynumber = bdttree->GetEntryNumberWithIndex(vtx_ID);
     if (bdt_entrynumber <0) continue;
     bdttree->GetEntry(bdt_entrynumber);
     //applying condition to select good vertex
     if (vtx_ntrks >= min_ntracks && bdtvalue > min_bdtvalue){
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
