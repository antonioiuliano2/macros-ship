//compare values for all and only reconstructed events
/*
 Added scripts for writing information in tree friends (DO NOT COPY THEM FROM ONE FILE TO ANOTHER),
 add_dsresults() for saving the events numbers where it was reconstructed
 findquarter() to get the ones in different quarter (we will need to select only the ones in second quarter)
 comparedistributions() to produce the efficiency plots 

*/
using namespace ROOT;

RVec<int> whichcharm (RVec<int> pdgcodes){
 RVec<int> whichcharm;
 for (const int &code : pdgcodes){
  if (TMath::Abs(code)==421) whichcharm.push_back(0);
  else if (TMath::Abs(code)==411) whichcharm.push_back(1);
  else if (TMath::Abs(code)==431) whichcharm.push_back(2);
  else if (TMath::Abs(code)==4122) whichcharm.push_back(3);
  else if (TMath::Abs(code)==4132) whichcharm.push_back(4);
  else if (TMath::Abs(code)==4232) whichcharm.push_back(5);
  else if (TMath::Abs(code)==4332) whichcharm.push_back(6);
  else if (TMath::Abs(code)==441) whichcharm.push_back(7);
  else whichcharm.push_back(-1);

 }
 return whichcharm;
}

RVec<int> decayinsecondquarter (int quarter){

 if (quarter==2){ 
   RVec<int> decaysinsecond = {1,1};
   return decaysinsecond;
   }
 else{
   RVec<int> decaysinsecond = {0,0};
   return decaysinsecond;
 }

}


//show subsample of distributions according to selection

RVec<int> selectedintdistribution (RVec<int> distribution, RVec<int> selection){
  
  return distribution[selection];

}
RVec<int> selectedintdistribution_2sels (RVec<int> distribution, RVec<int> selection1, RVec<int> selection2){
  
  return distribution[selection1*selection2];

}

RVec<float> selecteddistribution (RVec<float> distribution, RVec<int> selection){
  
  return distribution[selection];

}

RVec<float> selecteddistribution_2sels (RVec<float> distribution, RVec<int> selection1, RVec<int> selection2){
  
  return distribution[selection1*selection2];

}

//rename bins according to charm topologies

void namebins(ROOT::RDF::RResultPtr<TH1D> hcharmtype){
  TDatabasePDG *pdg = TDatabasePDG::Instance();

  hcharmtype->GetXaxis()->SetBinLabel(1,pdg->GetParticle(421)->GetName());
  hcharmtype->GetXaxis()->SetBinLabel(2,pdg->GetParticle(411)->GetName());
  hcharmtype->GetXaxis()->SetBinLabel(3,pdg->GetParticle(431)->GetName());
  hcharmtype->GetXaxis()->SetBinLabel(4,pdg->GetParticle(4122)->GetName());
  hcharmtype->GetXaxis()->SetBinLabel(5,pdg->GetParticle(4132)->GetName());
  hcharmtype->GetXaxis()->SetBinLabel(6,pdg->GetParticle(4232)->GetName());
  hcharmtype->GetXaxis()->SetBinLabel(7,pdg->GetParticle(4332)->GetName());
  hcharmtype->GetXaxis()->SetBinLabel(8,pdg->GetParticle(441)->GetName());
}

void plotoriginalvsreco(ROOT::RDF::RResultPtr<TH1D> horiginal, ROOT::RDF::RResultPtr<TH1D> hreconstructed, ROOT::RDF::RResultPtr<TH1D> hds){ //to do for int variable

 TCanvas *c = new TCanvas();
 c->Divide(1,2);
 c->cd(1);
 horiginal->DrawClone();
 hreconstructed->SetLineColor(kRed);
 hreconstructed->DrawClone("sames");
 hds->SetLineColor(kBlue);
 hds->DrawClone("sames");

 c->GetPad(1)->BuildLegend();

 c->cd(2);
 TEfficiency *peff = new TEfficiency(*(hreconstructed.GetPtr()),*(horiginal.GetPtr()));
 peff->Draw();
 TEfficiency *peff2 = new TEfficiency(*(hds.GetPtr()),*(horiginal.GetPtr()));
 peff2->SetLineColor(kBlue);
 peff2->Draw("SAMES");
}

void comparedistributions(){
  TFile *file = TFile::Open("distributions_mctrue_withquarter_withds.root");
  //connecting trees through addfriend
  TTree *mcdistributions = (TTree*) file->Get("charmdecays");
  TTree *dsresults = (TTree*) file->Get("dsrecoevents");
  TTree *quarterstree = (TTree*) file->Get("quarters"); //we need to know which events belong to second quarter

  mcdistributions->AddFriend(dsresults,"ds");
  mcdistributions->AddFriend(quarterstree,"quartersinfo");

  RDataFrame distdataframe = RDataFrame(*mcdistributions);
  //only the first 10000 events are actually passed to FEDRA for reconstructions
  const int nevents = 10000;
  auto dfpassedtofedra = distdataframe.Range(0,nevents);
  auto dfdefinedselection = dfpassedtofedra.Define("decayinsecondquarter",decayinsecondquarter,{"quarter"}); 
  auto dfdecaylength = dfdefinedselection.Define("dl",[](RVec<float> dx, RVec<float> dy, RVec<float> dz){RVec<float>length=sqrt(dx*dx+dy*dy+dz*dz);return length;},{"dx","dy","dz"});
  auto dfcharmnames =  dfdecaylength.Define("charmtype",whichcharm,{"pdgcode"});

  auto hcharmtype = dfcharmnames.Define("charmtype_inquarter",selectedintdistribution,{"charmtype","decayinsecondquarter"}).Histo1D({"hcharmtype","Charmed hadron",8,0,8},"charmtype_inquarter");
  auto hcharmtype_reco = dfcharmnames.Define("charmtype_reco",selectedintdistribution_2sels,{"charmtype","reconstructed","decayinsecondquarter"}).Histo1D({"hcharmtypereco","Reco Charmed hadron",8,0,8},"charmtype_reco");
  auto hcharmtype_ds = dfcharmnames.Define("charmtype_ds",selectedintdistribution_2sels,{"charmtype","dsreco","decayinsecondquarter"}).Histo1D({"hcharmtypeds","Ds Reco Charmed hadron",8,0,8},"charmtype_ds");

  namebins(hcharmtype);
  namebins(hcharmtype_reco);
  namebins(hcharmtype_ds);

  plotoriginalvsreco(hcharmtype,hcharmtype_reco,hcharmtype_ds);

  //lambda functions, to capture directly dataframe and apply operations for drawing
  auto compareoriginalvsreco = [&dfcharmnames] (TString columnvariable,ROOT::RDF::TH1DModel histoparameters,bool isint = false){ //for all float/double variables is ok, setisint to true for int variables

    ROOT::RDF::RResultPtr<TH1D> horiginal;
    ROOT::RDF::RResultPtr<TH1D> hreco;
    ROOT::RDF::RResultPtr<TH1D> hds;
    if (!isint){ 
     horiginal = dfcharmnames.Define((columnvariable+TString("_inquarter")).Data(),selecteddistribution,{columnvariable.Data(),"decayinsecondquarter"}).
                                  Histo1D(histoparameters,(columnvariable+TString("_inquarter")).Data());

     hds = dfcharmnames.
                          Define((columnvariable+TString("_ds")).Data(),selecteddistribution_2sels,{columnvariable.Data(),"dsreco","decayinsecondquarter"}).
                          Histo1D(histoparameters,(columnvariable+TString("_ds")).Data());
     hreco = dfcharmnames.
                          Define((columnvariable+TString("_reco")).Data(),selecteddistribution_2sels,{columnvariable.Data(),"reconstructed","decayinsecondquarter"}).
                          Histo1D(histoparameters,(columnvariable+TString("_reco")).Data());
    }
    else{ 
     horiginal = dfcharmnames.Define((columnvariable+TString("_inquarter")).Data(),selectedintdistribution,{columnvariable.Data(),"decayinsecondquarter"}).
                                  Histo1D(histoparameters,(columnvariable+TString("_inquarter")).Data());

     hds = dfcharmnames.
                          Define((columnvariable+TString("_ds")).Data(),selectedintdistribution_2sels,{columnvariable.Data(),"dsreco","decayinsecondquarter"}).
                          Histo1D(histoparameters,(columnvariable+TString("_ds")).Data());     
     hreco = dfcharmnames.
                         Define((columnvariable+TString("_reco")).Data(),selectedintdistribution_2sels,{columnvariable.Data(),"reconstructed","decayinsecondquarter"}).
                         Histo1D(histoparameters,(columnvariable+TString("_reco")).Data());   
    }
    TCanvas *c = new TCanvas(); //renaming reco histogram with a different name
    c->Divide(1,2);
    c->cd(2);
    TEfficiency *peff = new TEfficiency(*(hreco.GetPtr()), *(horiginal.GetPtr()));
    TEfficiency *peff2 = new TEfficiency(*(hds.GetPtr()), *(horiginal.GetPtr()));
    peff->Draw();
    peff2->SetLineColor(kBlue);
    peff2->Draw("SAMES");
    c->cd(1);
    hreco->SetName((TString(hreco->GetName())+TString("_reco")).Data());
    hreco->SetTitle((TString(hreco->GetTitle())+TString("_reco")).Data());
    hds->SetName((TString(hreco->GetName())+TString("_ds")).Data());
    hds->SetTitle((TString(hreco->GetTitle())+TString("_ds")).Data());

    
    horiginal->DrawClone("histo");

    hreco->SetLineColor(kRed); 
    hreco->DrawClone("histo && SAMES");

    hds->SetLineColor(kBlue);
    hds->DrawClone("histo && SAMES");


    c->GetPad(1)->BuildLegend();

  };

  //drawing histograms
  compareoriginalvsreco("dl", {"hdl", "Decay length of Charmed Hadron Decay;dl[#mum]",30,0,30000});
  compareoriginalvsreco("dz", {"hdz", "DZ of Charmed Hadron Decay;z[#mum]",30,0,30000});
  compareoriginalvsreco("momentum",{"hmomentum","Momentum of charmed hadron;P[GeV/c]",40,0,400});
  compareoriginalvsreco("gamma",{"hgamma","Gamma factor",200,0,200});
  compareoriginalvsreco("nprong",{"hmolt","Number of prongs of charmed hadron;nprong",10,0,10},true);
}

void add_dsresults(){
  //opening file and getting tree (this time, remember to create a separate tree)
  TFile *inputfile = TFile::Open("annotated_ds_data_result.root");
  TTreeReader dsreader("ds",inputfile);

  TTreeReaderArray<int> mcev(dsreader, "dsvtx_vtx2_mc_ev"); 
  TTreeReaderArray<int> mcpid(dsreader, "dsvtx_vtx2_mc_pid"); 
  TTreeReaderArray<int> charmdaughter(dsreader, "dsvtx_vtx2_trk_charmdaughter"); 
  TTreeReaderArray<int> samevent(dsreader, "dsvtx_vtx2_trk_samevent"); 
  
  const int nevents = 10000;
  const int ncharm = 2;
  int dsfound[nevents][ncharm];
  for(int ievent=0;ievent < nevents; ievent++){ 
   for(int icharm=0;icharm<ncharm;icharm++){ 
	dsfound[ievent][icharm] = 0;
	}
   }
  cout<<"starting loop on vertices"<<endl;
  const int nvertices = dsreader.GetEntries();
  for (int ivertex =0; ivertex < nvertices; ivertex++){
   dsreader.SetEntry(ivertex);
   const int ntracks = mcev.GetSize();
   //start loop on tracks
   for (int itrk = 0; itrk < ntracks; itrk++){
    if (charmdaughter[itrk] && samevent[itrk]){ //found charm daughter
      dsfound[mcev[itrk]][mcpid[itrk]-1] = 1;
    }
  }
  
 }
 cout<<"End loop on vertices"<<endl;
 TFile * outputfile = new TFile("distributions_mctrue_test.root","UPDATE");
 TTree * dsrecotree = new TTree("dsrecoevents","Events reconstructed by Valerio");

 int dsreco[ncharm];
 dsrecotree->Branch("dsreco",dsreco,"dsreco[2]/I");

 for (int ievent=0; ievent<nevents;ievent++){
   for (int icharm =0; icharm < ncharm;icharm++){
    dsreco[icharm] = dsfound[ievent][icharm];
   }
    dsrecotree->Fill(); //be careful when you do the filling
  }
 outputfile->cd();
 dsrecotree->Write();
 outputfile->Close();
}



//script to find which quarter each event belongs to
void findquarter(){
 TFile *simfile = TFile::Open("inECC_ship.conical.Pythia8CharmOnly-TGeant4_dig.root");
 TTreeReader simtree("cbmsim",simfile);
  
 TTreeReaderArray<BoxPoint> emupoints(simtree,"BoxPoint");
 TTreeReaderArray<BoxPoint> spreademupoints(simtree,"EmuBaseTrks"); 

 //histograms
 TH1I* hquarter = new TH1I("hquarter","which quarter",4,1,5);
 TH2D * hxy = new TH2D("hxy","XY spread from each event;x[#mum];y[#mum]",125,0,125000,100,0,100000);

 //moving target sim parameters
 double xcut = 62500;
 double ycut = 50000;

 double spilldy = 2.0;
 double targetmoverspeed = 2.6;

 //file to update with quarter information

 int quarter = 0;

 TFile * outputfile = new TFile("distributions_mctrue_withquarter.root","UPDATE");
 TTree * outputtree = new TTree("quarters","Quarters information");
 outputtree->Branch("quarter",&quarter,"quarter/I");

 const int nevents = 10000;
 //start loop within events
 for (int ievent = 0; ievent < nevents;ievent++){
  simtree.SetEntry(ievent);
  quarter = 1;
  if (ievent%1000 == 0)
     cout<<"arrived at "<<ievent<<endl;

  double deltat = spreademupoints[0].GetTime() - emupoints[0].GetTime();

  int nspill = deltat/100;

  double pottime = deltat - nspill*100;

  double spreadx = emupoints[0].GetX() -12.5/2. + pottime * targetmoverspeed;
  double spready = emupoints[0].GetY() - 9.9/2. + nspill * spilldy;


  spreadx = spreadx* 1E+4 + 62500;
  spready = spready* 1E+4 + 49500;   

  if( spreadx > xcut) quarter = quarter + 1;
  if( spready > ycut) quarter = quarter + 2;
   
  if (spready < 0) quarter = quarter - 10;
  //saving output
   
  hquarter->Fill(quarter);
  hxy->Fill(spreadx,spready);
  outputtree->Fill();
 }
//end loop over event

TCanvas * cxy = new TCanvas();
hxy->Draw("COLZ");

for (int ibin = 0; ibin < hquarter->GetNbinsX();ibin++)
  cout<<"In quarter "<<ibin+1<<"how many events ?"<<hquarter->GetBinContent(ibin+1)<<endl;

TCanvas * cquarter = new TCanvas();
hquarter->Draw();

outputfile->cd();
outputtree->Write();
outputfile->Close();
}
