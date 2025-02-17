EdbPVRec     *gAli=0;
  void draw1D_comparison(TCanvas *c, TH1D * hall, TH1D *h_matched, TH1D *h_notmatched){
      c->Divide(1,2);
      c->cd(1);
      hall->Draw();
      hall->SetLineWidth(2);
      h_matched->SetLineColor(kRed);
      h_matched->SetLineWidth(2);
      h_matched->Draw("SAMES");
      h_notmatched->SetLineColor(kGreen);
      h_notmatched->SetLineWidth(2);
      h_notmatched->Draw("SAMES");
      gPad->BuildLegend();
      c->cd(2);
      TEfficiency *heff_matched = new TEfficiency(*h_matched, *hall);
      heff_matched->SetLineColor(kRed);
      heff_matched->Draw();
      gPad->Update(); 
      heff_matched->GetPaintedGraph()->GetYaxis()->SetRangeUser(0,1.1);
      TEfficiency *heff_notmatched = new TEfficiency(*h_notmatched, *hall);
      heff_notmatched->SetLineColor(kGreen);
      heff_notmatched->Draw("SAMES");
      heff_matched->SetTitle("Ratio of matched tracks");
      heff_notmatched->SetTitle("Ratio of not matched tracks");
      gPad->BuildLegend();

  }

void matchedtracks(){
    TFile *indecesfile = TFile::Open("indexestracks_goodvertices.root");
    TTree *indecestree = (TTree*) indecesfile->Get("goodtrks");

   // int ntracks = indecestree->GetEntries();
//    int trackid;
    indecestree->BuildIndex("trid");
    //file with matched tree info for all spills
    TFile *matchedfile_8 = TFile::Open("CH1R6_Spill_8.root");
    TFile *matchedfile_9 = TFile::Open("CH1R6_Spill_9.root");
    TFile *matchedfile_10 = TFile::Open("CH1R6_Spill_10.root");
    TFile *matchedfile_11 = TFile::Open("CH1R6_Spill_11.root");
    TFile *matchedfile_12 = TFile::Open("CH1R6_Spill_12.root");
    
    const int nspills = 5;
    TTree *matchedtree[nspills];
    matchedtree[0] = (TTree*) matchedfile_8->Get("matched_trks");
    matchedtree[1] = (TTree*) matchedfile_9->Get("matched_trks");
    matchedtree[2] = (TTree*) matchedfile_10->Get("matched_trks");
    matchedtree[3] = (TTree*) matchedfile_11->Get("matched_trks");
    matchedtree[4] = (TTree*) matchedfile_12->Get("matched_trks");                
    for (int ispill = 0; ispill < nspills; ispill++)
     matchedtree[ispill]->BuildIndex("emu_trk"); //i want to know which tracks are present


    int nmatched = 0, nlast = 0;
    bool ismatched;

    //cuts
    const int maxnseg = 28;
    const float txoffset = 0.00215;
    const float tyoffset = -0.00464;
    double maxtx = 0.150;
    double maxty = 0.150;
    
    
    //reading tracks from linked_tracks.root file
    TFile *tracksfile = TFile::Open("linked_tracks.root");
    TTreeReader tracksreader("tracks",tracksfile);

    int ntracks = tracksreader.GetEntries();
    cout<<"Starting loop on "<<ntracks<<" tracks"<<endl;

    TTreeReaderValue<int> npl(tracksreader, "npl");
    TTreeReaderValue<int> nseg(tracksreader, "nseg");
    TTreeReaderArray<EdbSegP> segments(tracksreader, "s");
   
    //*********************************************list of histograms************************************************/ 
    TH2D *hxy = new TH2D("hxy","xy of last segment;x[#mu m];y[#mu m]",120,0,120000,100,0,100000);
    TH2D *hxy_matched = new TH2D("hxy_matched","xy of last segment for matched tracks; x[#mu m];y[#mu m]",120,0,120000,100,0,100000);
    TH2D *hxy_notmatched = new TH2D("hxy_notmatched","xy of last segment for not matched tracks; x[#mu m];y[#mu m]",120,0,120000,100,0,100000);
    
    TH1D *heff = new TH1D("heff","Efficiency;nseg/npl",20,0,1);
    TH1D *heff_matched = new TH1D("heff_matched","Efficency for matched tracks;nseg/npl",20,0,1);
    TH1D *heff_notmatched = new TH1D("heff_notmatched","Efficency for not matched tracks;nseg/npl",20,0,1);
    
    TH1D *hnseg = new TH1D("hnseg","Number of segments",30,0,30);
    TH1D *hnseg_matched = new TH1D("hnseg_matched","Nseg for matched tracks",30,0,30);
    TH1D *hnseg_notmatched = new TH1D("hnseg_notmatched","Nseg for not matched tracks",30,0,30);

    TH1D *htheta = new TH1D("htheta","Theta angle of last segment;#theta[rad]", 44,0,0.22);
    TH1D *htheta_matched = new TH1D("htheta_matched","Theta for matched tracks;#theta[rad]", 44,0,0.22);
    TH1D *htheta_notmatched = new TH1D("htheta_notmatched","Theta for not matched tracks;#theta[rad]", 44,0,0.22);

    TH1D *hphi = new TH1D("hphi","Phi angle of last segment;#phi[rad]", 140,-3.5,3.5);
    TH1D *hphi_matched = new TH1D("hphi_matched","Phi for matched tracks;#phi[rad]", 140,-3.5,3.5);
    TH1D *hphi_notmatched = new TH1D("hphi_notmatched","Phi for not matched tracks;#phi[rad]",140,-3.5,3.5);

    TH1D *hstartz = new TH1D("hstartz","Z position of first segment;z[#mum]", 72,-36000,0);
    TH1D *hstartz_matched = new TH1D("hstartz_matched","Z position of first segment for matched tracks;z[#mum]", 72,-36000,0);
    TH1D *hstartz_notmatched = new TH1D("hstartz_notmatched ","Z position of first segment for not matched tracks;z[#mum]",72,-36000,0);
    

    //**********************************************START LOOP*****************************************************
    for (int itrk = 0; itrk < ntracks;itrk++){
       
        //if (itrk%10000 == 0) cout<<"Arrived at track "<<itrk<<endl;
        //indecestree->GetEntry(itrk);
        tracksreader.SetEntry(itrk);

        EdbSegP lastsegment = segments[*nseg -1];
        EdbSegP firstsegment = segments[0];
    
        if (lastsegment.PID() != 0 || indecestree->GetEntryNumberWithIndex(itrk) < 0) continue; //we want tracks ending in the brick
        if (TMath::Abs(lastsegment.TX()+txoffset) > maxtx || TMath::Abs(lastsegment.TY()+tyoffset) > maxty) continue;
        if (*nseg > maxnseg) continue;

        nlast++;
        double eff= (double) *nseg / *npl;

        hxy->Fill(lastsegment.X(),lastsegment.Y());
        hnseg->Fill(*nseg);
        heff->Fill(eff);
        htheta->Fill(TMath::ATan(TMath::Sqrt(pow(lastsegment.TX()+txoffset,2)+pow(lastsegment.TY()+tyoffset,2))));
        hphi->Fill(TMath::ATan2(lastsegment.TY()+tyoffset,lastsegment.TX()+txoffset));
        hstartz->Fill(firstsegment.Z());

        //look in all files if track was matched
        ismatched = false;
        for (int ispill = 0; ispill < nspills; ispill++)
         if(matchedtree[ispill]->GetEntryNumberWithIndex(itrk) >= 0) ismatched = true;
       
        if (ismatched){
            nmatched++;
            
            hxy_matched->Fill(lastsegment.X(),lastsegment.Y());
            hnseg_matched->Fill(*nseg);
            heff_matched->Fill(eff);
            htheta_matched->Fill(TMath::ATan(TMath::Sqrt(pow(lastsegment.TX()+txoffset,2)+pow(lastsegment.TY()+tyoffset,2))));
            hphi_matched->Fill(TMath::ATan2(lastsegment.TY()+tyoffset,lastsegment.TX()+txoffset));
            hstartz_matched->Fill(firstsegment.Z());
        
        }

        else{
            
            hxy_notmatched->Fill(lastsegment.X(),lastsegment.Y());
            hnseg_notmatched->Fill(*nseg);
            heff_notmatched->Fill(eff);
            htheta_notmatched->Fill(TMath::ATan(TMath::Sqrt(pow(lastsegment.TX()+txoffset,2)+pow(lastsegment.TY()+tyoffset,2))));
            hphi_notmatched->Fill(TMath::ATan2(lastsegment.TY()+tyoffset,lastsegment.TX()+txoffset));
            hstartz_notmatched->Fill(firstsegment.Z());
        }




    }

  cout<<"Total: "<<nlast<<" Matched tracks: "<<nmatched<<endl;

  //**********************************DRAWING HISTOGRAMS**********************************//
  //1D histograms

  TCanvas *cnseg = new TCanvas();
  draw1D_comparison(cnseg, hnseg, hnseg_matched, hnseg_notmatched);

  TCanvas *ceff = new TCanvas();
  draw1D_comparison(ceff, heff, heff_matched, heff_notmatched);

  TCanvas *ctheta = new TCanvas();
  draw1D_comparison(ctheta, htheta, htheta_matched, htheta_notmatched);

  TCanvas *cphi = new TCanvas();
  draw1D_comparison(cphi, hphi, hphi_matched, hphi_notmatched);

  TCanvas *cstartz = new TCanvas();
  draw1D_comparison(cstartz, hstartz, hstartz_matched, hstartz_notmatched);

  //2D histograms

  TCanvas *cxy = new TCanvas();
  cxy->Divide(1,2);
  cxy->cd(1);
  hxy->Draw("COLZ");
  cxy->cd(2);
  hxy_matched->Draw("COLZ");
}

//checking tracks candidates for matching without the matching
void checkgoodtracks(){
    TFile *indecesfile = TFile::Open("indexestracks_goodvertices.root");
    TTree *indecestree = (TTree*) indecesfile->Get("goodtrks");
    indecestree->BuildIndex("trid");

    //cuts
    const int maxnseg = 28;
    const float txoffset = 0.00215;
    const float tyoffset = -0.00464;
    double maxtx = 0.150;
    double maxty = 0.150;
    
    int nlast = 0;
    
    //reading tracks from linked_tracks.root file
    TFile *tracksfile = TFile::Open("linked_tracks.root");
    TTreeReader tracksreader("tracks",tracksfile);

    int ntracks = tracksreader.GetEntries();
    cout<<"Starting loop on "<<ntracks<<" tracks"<<endl;

    TTreeReaderValue<int> npl(tracksreader, "npl");
    TTreeReaderValue<int> nseg(tracksreader, "nseg");
    TTreeReaderArray<EdbSegP> segments(tracksreader, "s");
   
    //*********************************************list of histograms************************************************/ 
    TH2D *hxy = new TH2D("hxy","xy of last segment;x[#mu m];y[#mu m]",120,0,120000,100,0,100000);
    TH1D *hx = new TH1D("hx","x of last segment;x[#mu m]",120,0,120000);
    TH1D *hy = new TH1D("hy","y of last segment;y[#mu m]",120,0,120000);
    TH1D *hxall = new TH1D("hxall","x of last segment all;x[#mu m]",120,0,120000);
    TH1D *hyall = new TH1D("hyall","y of last segment all;y[#mu m]",100,0,100000);
    
    TH1D *heff = new TH1D("heff","Efficiency;nseg/npl",20,0,1);
    
    TH1D *hnseg = new TH1D("hnseg","Number of segments",29,0,30);

    TH1D *htheta = new TH1D("htheta","Theta angle of last segment;#theta[rad]", 44,0,0.22);

    TH1D *hphi = new TH1D("hphi","Phi angle of last segment;#phi[rad]", 140,-3.5,3.5);
 
    //**********************************************START LOOP*****************************************************
    for (int itrk = 0; itrk < ntracks;itrk++){
       
        //if (itrk%10000 == 0) cout<<"Arrived at track "<<itrk<<endl;
        //indecestree->GetEntry(itrk);
        tracksreader.SetEntry(itrk);

        EdbSegP lastsegment = segments[*nseg -1];
        //if (lastsegment.PID()==0) cout<<lastsegment.PID()<<endl;
        if(indecestree->GetEntryNumberWithIndex(itrk) >= 0){
         hxall->Fill(lastsegment.X());
         hyall->Fill(lastsegment.Y());
        } 
        if (lastsegment.PID() > 1 || indecestree->GetEntryNumberWithIndex(itrk) < 0) continue; //we want tracks ending in the brick
        //cout<<"TEST0"<<endl;
        if (TMath::Abs(lastsegment.TX()+txoffset) > maxtx || TMath::Abs(lastsegment.TY()+tyoffset) > maxty) continue;
        //cout<<"TEST"<<endl;
        if (*nseg > maxnseg) continue;

        nlast++;
        double eff= (double) *nseg / *npl;
        hx->Fill(lastsegment.X());
        hy->Fill(lastsegment.Y());
        hxy->Fill(lastsegment.X(),lastsegment.Y());
        hnseg->Fill(*nseg);
        heff->Fill(eff);
        htheta->Fill(TMath::ATan(TMath::Sqrt(pow(lastsegment.TX()+txoffset,2)+pow(lastsegment.TY()+tyoffset,2))));
        hphi->Fill(TMath::ATan2(lastsegment.TY()+tyoffset,lastsegment.TX()+txoffset));

    }

  cout<<"Total: "<<nlast<<endl;

  //**********************************DRAWING HISTOGRAMS**********************************//
  //1D histograms

  TCanvas *cnseg = new TCanvas();
  hnseg->Draw();

  TCanvas *ceff = new TCanvas();
  heff->Draw();

  TCanvas *ctheta = new TCanvas();
  htheta->Draw();

  TCanvas *cphi = new TCanvas();
  hphi->Draw();

  TCanvas *cxyall = new TCanvas();
  cxyall->Divide(1,2);
  cxyall->cd(1);
  hx->Draw();
  cxyall->cd(2);
  hy->Draw();

  //2D histograms

  TCanvas *cxy = new TCanvas();
  hxy->Draw("COLZ");
}
