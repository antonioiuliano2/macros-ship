//compare values for all and only reconstructed events
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

//show subsample of distributions according to selection

RVec<int> selectedintdistribution (RVec<int> distribution, RVec<int> selection){
  
  return distribution[selection];

}

RVec<float> selecteddistribution (RVec<float> distribution, RVec<int> selection){
  
  return distribution[selection];

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

void plotoriginalvsreco(ROOT::RDF::RResultPtr<TH1D> horiginal, ROOT::RDF::RResultPtr<TH1D> hreconstructed){ //to do for int variable

 TCanvas *c = new TCanvas();
 c->Divide(1,2);
 c->cd(1);
 horiginal->DrawClone();
 hreconstructed->SetLineColor(kRed);
 hreconstructed->DrawClone("sames");

 c->GetPad(1)->BuildLegend();

 c->cd(2);
 TEfficiency *peff = new TEfficiency(*(hreconstructed.GetPtr()),*(horiginal.GetPtr()));
 peff->Draw();
}

void comparedistributions(){
  TFile *file = TFile::Open("distributions_mctrue.root");
  RDataFrame distdataframe = RDataFrame("charmdecays",file);
  //only the first 10000 events are actually passed to FEDRA for reconstructions
  const int nevents = 10000;
  auto dfpassedtofedra = distdataframe.Range(0,nevents);
  auto dfdecaylength = dfpassedtofedra.Define("dl",[](RVec<float> dx, RVec<float> dy, RVec<float> dz){RVec<float>length=sqrt(dx*dx+dy*dy+dz*dz);return length;},{"dx","dy","dz"});
  auto dfcharmnames =  dfdecaylength.Define("charmtype",whichcharm,{"pdgcode"});

  auto hcharmtype = dfcharmnames.Histo1D({"hcharmtype","Charmed hadron",8,0,8},"charmtype");
  auto hcharmtype_reco = dfcharmnames.Define("charmtype_reco",selectedintdistribution,{"charmtype","reconstructed"}).Histo1D({"hcharmtypereco","Reco Charmed hadron",8,0,8},"charmtype_reco");

  namebins(hcharmtype);
  namebins(hcharmtype_reco);

  plotoriginalvsreco(hcharmtype,hcharmtype_reco);

  //lambda functions, to capture directly dataframe and apply operations for drawing
  auto compareoriginalvsreco = [&dfcharmnames] (TString columnvariable,ROOT::RDF::TH1DModel histoparameters,bool isint = false){ //for all float/double variables is ok, setisint to true for int variables
    auto horiginal = dfcharmnames.Histo1D(histoparameters,(columnvariable).Data());

    ROOT::RDF::RResultPtr<TH1D> hreco;
    if (!isint) hreco = dfcharmnames.Define((columnvariable+TString("_reco")).Data(),selecteddistribution,{columnvariable.Data(),"reconstructed"}).Histo1D(histoparameters,(columnvariable+TString("_reco")).Data());
    else hreco = dfcharmnames.Define((columnvariable+TString("_reco")).Data(),selectedintdistribution,{columnvariable.Data(),"reconstructed"}).Histo1D(histoparameters,(columnvariable+TString("_reco")).Data());   
    
    TCanvas *c = new TCanvas(); //renaming reco histogram with a different name
    c->Divide(1,2);
    c->cd(2);
    TEfficiency *peff = new TEfficiency(*(hreco.GetPtr()), *(horiginal.GetPtr()));
    peff->Draw();
    c->cd(1);
    hreco->SetName((TString(hreco->GetName())+TString("_reco")).Data());
    hreco->SetTitle((TString(hreco->GetTitle())+TString("_reco")).Data());

    
    horiginal->DrawClone("histo");

    hreco->SetLineColor(kRed); 
    hreco->DrawClone("histo && SAMES");


//    c->BuildLegend();

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
  TTreeReader dsreader = TTreeReader("ds",inputfile);

  TTreeReaderArray<int> mcev = TTreeReaderArray(dsreader, "dsvtx_vtx2_mc_ev"); 
  TTreeReaderArray<int> mctrack = TTreeReaderArray(dsreader, "dsvtx_vtx2_mc_tid"); 
  TTreeReaderArray<int> charmdaughter = TTreeReaderArray(dsreader, "dsvtx_vtx2_trk_charmdaughter"); 
  TTreeReaderArray<int> sameevent = TTreeReaderArray(dsreader, "dsvtx_vtx2_trk_samevent"); 
  

}


