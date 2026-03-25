
namespace VERTEX_PAR
{
  float DZmax = 3000.;
  //float ProbMinV   = 0.0001;  // minimum acceptable probability for chi2-distance between tracks
  float ProbMinV   = 0.01;
  float ImpMax     = 15.;    // maximal acceptable impact parameter [microns] (for preliminary check)
  bool  UseMom     = false;  // use or not track momentum for vertex calculations
  bool  UseSegPar  = true;  // use only the nearest measured segments for vertex fit (as Neuchatel)
  int   QualityMode= 0;      // vertex quality estimation method (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)
}

void fragments_tracks(){
    
/*    TFile *inputfile = TFile::Open("vertextree.root");
    TTreeReader vtxtree("vtx",inputfile);
    //moltelplicity
    TTreeReaderValue<int> ntracks(vtxtree,"n");
    TTreeReaderArray<int> nsegments(vtxtree,"nseg");
    //MCTruth
    TTreeReaderArray<int> MCTrackID(vtxtree,"MCTrackID");
    TTreeReaderArray<EdbSegP> segments(vtxtree,"s");

    const int nvertices = vtxtree.GetEntries();

    TH1D *hDT = new TH1D("hDT","Angular difference;#Theta[rad]",20,-0.1,0.1);

    for (int ivtx = 0; ivtx < nvertices; ivtx++){
        vtxtree.SetEntry(ivtx);
        //check if same track has done vertex with itself
        if (*ntracks == 2 && MCTrackID[0]==MCTrackID[1]){
            int nsegmentsfirst = nsegments[0];
            hDT->Fill(segments[nsegmentsfirst-1].Theta()-segments[nsegmentsfirst].Theta());
        }
    }
    hDT->Draw();*/
    EdbDataProc * dproc = new EdbDataProc();
    EdbPVRec *gAli = dproc->PVR();
    EdbScanCond *scancond = new EdbScanCond();
    scancond->SetSigma0(5,5,0.003,0.003);
    gAli->SetScanCond(scancond);

    //getting list of vertices
    dproc->ReadVertexTree(*gAli,"vertextree.root","vID<50000",*scancond);
    TObjArray *vertexlist = gAli->eVTX;

    const int nvertices = vertexlist->GetEntries();

    cout<<"Starting loop over "<<nvertices<<" vertices "<<endl;
    TH2D *hdt = new TH2D("hdt","dT difference between fragments;TX;TY",20,-0.1,0.1,20,-0.1,0.1);
    TH2D *hdtcharmdaughter = new TH2D("hdtcharmdaughter","dT difference between fragments only charm daughter tracks;TX;TY",20,-0.1,0.1,20,-0.1,0.1);
    TH1I *hnaddseg = new TH1I("hnaddseg","Number of added segments;NSegments",30,0,30);
    for (int ivtx = 0; ivtx<nvertices;ivtx++){

     EdbVertex *testvertex = (EdbVertex*) vertexlist->At(ivtx);
     
     int firsttrackzpos = testvertex->GetVTa(0)->Zpos();
     EdbTrackP *firsttrack = testvertex->GetTrack(0);

     int secondtrackzpos = testvertex->GetVTa(1)->Zpos();
     EdbTrackP *secondtrack = testvertex->GetTrack(1);
     int myMCTrack = firsttrack->MCTrack();
     int addedsegs = 0;
     if (ivtx==25)     cout<<testvertex->N()<<" "<<firsttrack->MCTrack()<<" "<<secondtrack->MCTrack()<<endl;
     if(testvertex->N()==2&&firsttrack->MCTrack()==secondtrack->MCTrack()&&firsttrack->MCEvt()==secondtrack->MCEvt()){ //track doing a vertex with itself
        if (firsttrackzpos < secondtrackzpos){ //first track before second track
         for (int iseg = 0; iseg<secondtrack->N();iseg++){
             EdbSegP *segment = secondtrack->GetSegment(iseg);
             if (segment->MCTrack() == myMCTrack){ 
                 firsttrack->AddSegment(segment);
                 addedsegs++;
             }
         }         
         hnaddseg->Fill(addedsegs);
         EdbSegP * firstsegafter = secondtrack->GetSegment(0);
         EdbSegP * lastsegbefore = firsttrack->GetSegment(firsttrack->N()-1);

         hdt->Fill(firstsegafter->TX()-lastsegbefore->TX(),firstsegafter->TY()-lastsegbefore->TY());
         if(firsttrack->Aid(0)==0 || firsttrack->Aid(0)==1) hdtcharmdaughter->Fill(firstsegafter->TX()-lastsegbefore->TX(),firstsegafter->TY()-lastsegbefore->TY());
        }

        else if (firsttrackzpos>secondtrackzpos){
         for (int iseg = 0; iseg<firsttrack->N();iseg++){
             EdbSegP *segment = firsttrack->GetSegment(iseg);
             if (segment->MCTrack() == myMCTrack){ 
                 secondtrack->AddSegment(segment);
                 addedsegs++;
             }
         }
         hnaddseg->Fill(addedsegs);

         EdbSegP * firstsegafter = firsttrack->GetSegment(0);
         EdbSegP * lastsegbefore = secondtrack->GetSegment(firsttrack->N()-1);

         hdt->Fill(firstsegafter->TX()-lastsegbefore->TX(),firstsegafter->TY()-lastsegbefore->TY());
         if(firsttrack->Aid(0)==0 || firsttrack->Aid(0)==1) hdtcharmdaughter->Fill(firstsegafter->TX()-lastsegbefore->TX(),firstsegafter->TY()-lastsegbefore->TY());
        }
        testvertex->SetFlag(-500);
     }
     
    }
 //dproc->MakeVertexTree(*(gAli->eVTX),"vertextree_taggedsametrack.root");
 hnaddseg->Draw();

 TCanvas *c = new TCanvas();
 hdt->Draw("COLZ");

 TH1D * hdtx = hdt->ProjectionX();
 hdtx->SetTitle("dTX difference between fragments");
 TH1D * hdty = hdt->ProjectionY();
 hdty->SetTitle("dTY difference between fragments");
 
 TCanvas *c2d = new TCanvas();
 c2d->Divide(1,2);
 c2d->cd(1);
 hdtx->Draw();
 c2d->cd(2);
 hdty->Draw();


 TCanvas *ccharmdaughter = new TCanvas();
 hdtcharmdaughter->Draw("COLZ");

 TH1D * hdtxcharmdaughter = hdtcharmdaughter->ProjectionX();
 hdtxcharmdaughter->SetTitle("dTX difference between fragments");
 TH1D * hdtycharmdaughter = hdtcharmdaughter->ProjectionY();
 hdtycharmdaughter->SetTitle("dTY difference between fragments");
}
