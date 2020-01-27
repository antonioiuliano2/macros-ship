//add branches with shower information to ds tree
using namespace ROOT;
void addshowerinfo(){
    //creating maps
    map<int,int> showerid;
    map<int,int> showersize;
    map<int,float> showeroutput15;
    map<int,float> showeroutput30;
    int onesizeb;
    float oneoutput15, oneoutput30;
    int startertrack;
    
    float dummy;
    int ishower;
    for (int xcode = 6; xcode < 12; xcode++ ){
     for (int ycode = 0; ycode < 5; ycode++ ){
        ifstream showerlog (Form("/home/antonio/cernbox/Synched/Charmsim_Showreco/Showreco_ds_25_01/showerlog_%i_%i.dat",xcode*10000,ycode*10000),ifstream::in);
        string dummystring;
        getline(showerlog,dummystring);
        while (showerlog.good()){
                RVec<int> showerinfo;
                showerlog>>ishower;
                showerlog >> startertrack;
                showerlog >> dummy;
                showerlog >> dummy;
                showerlog >> onesizeb;
                showerlog >> oneoutput15;
                showerlog >> oneoutput30;

                showerid[startertrack] = ishower;
                showersize[startertrack] = onesizeb;
                showeroutput15[startertrack] = oneoutput15;
                showeroutput30[startertrack] = oneoutput30;
             }
        }
    }

    //tree with true MC simulation
    TFile *simfile = TFile::Open("../inECC_ship.conical.Pythia8CharmOnly-TGeant4_dig.root");
    TTreeReader simreader("cbmsim",simfile);

    TTreeReaderArray<ShipMCTrack> tracks(simreader,"MCTrack");

    //ds tree
    TFile * inputfile = TFile::Open("annotated_ds_data_result.root","READ");
    TTree * dstree = (TTree*) inputfile->Get("ds");
    TTreeReader dsreader("ds",inputfile);

    const int nevents = dstree->GetEntries();

    TTreeReaderArray<int> trackIDs(dsreader,"dsvtx_vtx2_tid");
    TTreeReaderArray<int> MCEventID(dsreader,"dsvtx_vtx2_mc_ev");
    TTreeReaderArray<int> MCTrackID(dsreader,"dsvtx_vtx2_mc_tid");


    TFile *outputfile = new TFile("annotated_ds_data_result_testshower_v4_full.root","UPDATE");
//    vector<int> trackIDs;

    vector<int> nshower;
    vector<int> sizeb;
    vector<float> output15;
    vector<float> output30;

    vector<int> pdgcode;
    vector<float> momentum;

    dstree->SetBranchStatus("*",1);
    dstree->SetBranchStatus("dsvtx_vtx2_tid",1);
    dstree->SetBranchStatus("dsvtx_vtx2_mc_ev",1);
    dstree->SetBranchStatus("dsvtx_vtx2_mc_tid",1);
//    dstree->SetBranchAddress("dsvtx_vtx2_tid",&trackIDs);
    TTree *dstree2 = new TTree("dsshower","Result of shower reconstruction");
    TBranch * branchsizeb = dstree2->Branch("dsvtx_vtx2_trk_sizeb",&sizeb);
    dstree2->SetBranchStatus("dsvtx_vtx2_trk_sizeb",1);
    TBranch * branchnshower = dstree2->Branch("dsvtx_vtx2_trk_nshower",&nshower);
    dstree2->SetBranchStatus("dsvtx_vtx2_trk_nshower",1);
    TBranch * branchoutput15 = dstree2->Branch("dsvtx_vtx2_trk_output15",&output15);
    dstree2->SetBranchStatus("dsvtx_vtx2_trk_output15",1);
    TBranch * branchoutput30 = dstree2->Branch("dsvtx_vtx2_trk_output30",&output30);
    dstree2->SetBranchStatus("dsvtx_vtx2_trk_output30",1);

    TBranch * branchmomentum = dstree2->Branch("dsvtx_vtx2_trk_mc_momentum",&momentum);
    dstree2->SetBranchStatus("dsvtx_vtx2_trk_mc_momentum",1);
    TBranch * branchpdgcode = dstree2->Branch("dsvtx_vtx2_trk_mc_pdgcode",&pdgcode);
    dstree2->SetBranchStatus("dsvtx_vtx2_trk_mc_pdgcode",1);
 
    
    cout<<"Start loop on Valerio's events "<<endl;

    for (int ivertex = 0; ivertex < nevents; ivertex++){
        //clearing vectors
//        trackIDs.clear();
        nshower.clear();
        sizeb.clear();
        output15.clear();
        output30.clear();

        momentum.clear();
        pdgcode.clear();

        dsreader.SetEntry(ivertex);
        //start loop on tracks
        if (ivertex%100 == 0) cout<<ivertex<<endl;
        for (int itrack = 0; itrack < trackIDs.GetSize(); itrack++){
            int whichtrack = trackIDs[itrack];

            //getting true mc info for that track
            int nexteventid = MCEventID[itrack];
            int nexttrackid = MCTrackID[itrack];
           // cout<<nexteventid<<" "<<nexttrackid<<endl;
            simreader.SetEntry(nexteventid);
            ShipMCTrack MCTrack = tracks[nexttrackid];

            momentum.push_back(MCTrack.GetP());
            pdgcode.push_back(MCTrack.GetPdgCode());

          //  cout<<"RIPROVA"<<endl;
            if(showersize.count(whichtrack) > 0){
            //    cout<<"BOH "<<whichtrack<<endl;
                nshower.push_back(showerid[whichtrack]);
                sizeb.push_back(showersize[whichtrack]);
                output15.push_back(showeroutput15[whichtrack]);
                output30.push_back(showeroutput30[whichtrack]);
            //    cout<<"RIBOH"<<endl;

            }
            else{//not found track, putting -10 to all
                nshower.push_back(-10);
                sizeb.push_back(-10);
                output15.push_back(-10);
                output30.push_back(-10);
            }

        }//end tracks loop
        //filling tree with new vectors (address must be set each for STL dynamic object, see TTree reference)
/*        cout<<nshower.size()<<endl;
        branchnshower->SetAddress(&nshower);
        branchsizeb->SetAddress(&sizeb);
        branchoutput15->SetAddress(&output15);
        branchoutput30->SetAddress(&output30);
        branchpdgcode->SetAddress(&pdgcode);
        branchmomentum->SetAddress(&momentum);*/
        dstree2->Fill();
    }
    outputfile->cd();
    dstree2->Write();
    inputfile->Close();
    outputfile->Close();
}
//signal selection (electron)
RVec<int> iselectron(RVec<int> pdgcode){
    return (abs(pdgcode) == 11);
}
//signal selection (electron)
RVec<int> isnotelectron(RVec<int> pdgcode){
    return (abs(pdgcode) != 11);
}

RVec<int> outputabovecut(RVec<float> output15, float cut ){
    return (output15>cut);
}

RVec<float> selecteddistribution(RVec<float> originaldistribution, RVec<int> selection){
  //fill a RVec only with accepted values of the RVec
    return originaldistribution[selection];
}
RVec<int> selecteddistribution_int(RVec<int> originaldistribution, RVec<int> selection){
  //fill a RVec only with accepted values of the RVec
    return originaldistribution[selection];
}

void plotsignalandbackground(ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> ds, TString columnvariable, ROOT::RDF::TH1DModel histoparameters){
    auto hsignal = ds.Define((columnvariable+TString("signal")).Data(),selecteddistribution,{columnvariable.Data(),"dsvtx_vtx2_trk_signal"}).Histo1D(histoparameters,(columnvariable+TString("signal")).Data());
    auto hbackground = ds.Define((columnvariable+TString("background")).Data(),selecteddistribution,{columnvariable.Data(),"dsvtx_vtx2_trk_background"}).Histo1D(histoparameters,(columnvariable+TString("background")).Data());
    
    TCanvas *c = new TCanvas();
    hsignal->SetName((TString(hsignal->GetName())+TString(" signal")).Data());
    hbackground->SetName((TString(hbackground->GetName())+TString(" background")).Data());
    hsignal->SetTitle((TString(hsignal->GetTitle())+TString(" signal")).Data());
    hbackground->SetTitle((TString(hbackground->GetTitle())+TString(" background")).Data());

    hsignal->Scale(1./hsignal->Integral());
    hbackground->Scale(1./hbackground->Integral());

    //setting scale
    double ymaximum = hsignal->GetMaximum();
    if (hbackground->GetMaximum()>ymaximum) ymaximum = hbackground->GetMaximum();

    hsignal->GetYaxis()->SetRangeUser(0,ymaximum+(ymaximum/10.));
    hsignal->DrawClone("histo");
    hbackground->SetLineColor(kRed);
    hbackground->DrawClone("histo && SAMES");

    c->BuildLegend();
}
// if variable is an int
void plotsignalandbackground_int(ROOT::RDF::RInterface<ROOT::Detail::RDF::RLoopManager, void> ds, TString columnvariable, ROOT::RDF::TH1DModel histoparameters){

    auto hsignal = ds.Define((columnvariable+TString("signal")).Data(),selecteddistribution_int,{columnvariable.Data(),"dsvtx_vtx2_trk_signal"}).Histo1D(histoparameters,(columnvariable+TString("signal")).Data());
    auto hbackground = ds.Define((columnvariable+TString("background")).Data(),selecteddistribution_int,{columnvariable.Data(),"dsvtx_vtx2_trk_background"}).Histo1D(histoparameters,(columnvariable+TString("background")).Data());

    TCanvas *c = new TCanvas();
    //they must have different names
    hsignal->SetName((TString(hsignal->GetName())+TString(" signal")).Data());
    hbackground->SetName((TString(hbackground->GetName())+TString(" background")).Data());
    hsignal->SetTitle((TString(hsignal->GetTitle())+TString(" signal")).Data());
    hbackground->SetTitle((TString(hbackground->GetTitle())+TString(" background")).Data());

    hsignal->Scale(1./hsignal->Integral());
    hbackground->Scale(1./hbackground->Integral());
    
    //setting scale
    double ymaximum = hsignal->GetMaximum();
    if (hbackground->GetMaximum()>ymaximum) ymaximum = hbackground->GetMaximum();

    hsignal->GetYaxis()->SetRangeUser(0,ymaximum+(ymaximum/10.)); //a bit above maximum

    hsignal->DrawClone("histo");
    hbackground->SetLineColor(kRed);
    hbackground->DrawClone("histo && SAMES");

    c->BuildLegend();
}

//check how algorithm behaves for electrons and not electrons
void checkshowerefficiency(){
    TFile *inputfile = TFile::Open("annotated_ds_data_result_testshower_v3_noinputtrackduplicates.root");
    RDataFrame dsdataframe = RDataFrame("ds",inputfile);

    //defining signal branch according to selection
    auto dsdataframe_selection = dsdataframe.Define("dsvtx_vtx2_trk_signal",iselectron,{"dsvtx_vtx2_trk_pdgcode"});
    auto dsdataframe_cut = dsdataframe_selection.Define("dsvtx_vtx2_trk_background",isnotelectron,{"dsvtx_vtx2_trk_pdgcode"});
    //histograms above cut
    //auto dsdataframe_cut = dsdataframe_selection2.Define("dsvtx_vtx2_trk_outputabovecut",outputabovecut,{"dsvtx_vtx2_trk_output15"});

    plotsignalandbackground(dsdataframe_cut,"dsvtx_vtx2_trk_output15",{"houtput15","Output 15;output", 200,-15,5}); 
    plotsignalandbackground(dsdataframe_cut,"dsvtx_vtx2_trk_output30",{"houtput30", "Output 30;output", 200,-15,5}); 

    plotsignalandbackground_int(dsdataframe_cut,"dsvtx_vtx2_trk_sizeb",{"hsizeb","Dimension of shower",10,0,100});


    //prepare definitions for 2d histograms
    auto selectionsignal = dsdataframe_cut.Define("output15signal",selecteddistribution,{"dsvtx_vtx2_trk_output15","dsvtx_vtx2_trk_signal"})
                                                 .Define("momentumsignal",selecteddistribution,{"dsvtx_vtx2_trk_momentum","dsvtx_vtx2_trk_signal"});
    auto selectionbackground = dsdataframe_cut.Define("output15background",selecteddistribution,{"dsvtx_vtx2_trk_output15","dsvtx_vtx2_trk_background"})
                                                     .Define("momentumbackground",selecteddistribution,{"dsvtx_vtx2_trk_momentum","dsvtx_vtx2_trk_background"});                                                  

    auto hmomoutputsignal=selectionsignal.Profile1D({"hmomoutputsig","Momentum vs output15 signal;output15;P[GeV/c]",20,-15,5,0,400},"output15signal","momentumsignal");
    auto hmomoutputbackground=selectionbackground.Profile1D({"hmomoutputbkg","Momentum vs output15 background;output15;P[GeV/c]",20,-15,5,0,400},"output15background","momentumbackground");

    auto houtputmomsignal=selectionsignal.Profile1D({"houtputmomsig","Momentum vs output15 signal;P[GeV/c];output15",40,0,400,-15,5},"momentumsignal","output15signal");
    auto houtputmombackground=selectionbackground.Profile1D({"houtputmombkg","Momentum vs output15 background;P[GeV/c];output15",40,0,400,-15,5},"momentumbackground","output15background");

    TCanvas *c2d = new TCanvas();
    //hmomoutput->GetYaxis()->SetRangeUser(-15,5);
    c2d->Divide(2,1);
    c2d->cd(1);
    hmomoutputbackground->SetLineColor(kRed);
    hmomoutputbackground->DrawClone();
    hmomoutputsignal->DrawClone("sames");
    c2d->GetPad(1)->BuildLegend();
    c2d->cd(2);
    houtputmombackground->SetLineColor(kRed);
    houtputmombackground->DrawClone();
    houtputmomsignal->DrawClone("sames");
    c2d->GetPad(2)->BuildLegend();
}
