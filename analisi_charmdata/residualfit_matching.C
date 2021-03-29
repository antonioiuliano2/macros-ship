//Fitting with ROOFIT file provided by Chris
void residualfit_matching()
{
  TFile *inputfile = TFile::Open("scaledtree.root");
  TTree *inputtree = (TTree*) inputfile->Get("scaled_matched_tracks");

  TTree *goodtree = inputtree->CopyTree("TMath::Abs(scaled_dyr)<40");
  const double cm2micron = 1e+4;
  RooRealVar dxr("scaled_dxr", "match_dxr[#mum]", -0.025*cm2micron, 0.025*cm2micron);
  RooDataSet matcheddata("matcheddata","dataset with match_dxr",goodtree,dxr) ;

  //adding binning, try binned
  //dxr.setBins(30);
  //RooDataHist binnedData("binnedData", "matcheddata",
  //RooArgList(dxr), matcheddata);

  //define parameters and PDFs 
  RooRealVar mu("mu", "average", 0, -1000, 1000);
  RooRealVar sigma("sigma", "sigma",40, 0, 100);
  RooGaussian gauss("gauss","gaussian PDF", dxr, mu, sigma);

  RooRealVar mu2("mu2", "average of second gaussian", 0, -1000, 1000);
  RooRealVar sigma2("sigma2", "sigma of second gaussian", 20, 0, 100);
  RooGaussian gauss2("gauss2","Second gaussian PDF", dxr, mu2, sigma2);

  //polynomial, max 1 degree
  //RooRealVar slope("slope", "slope",20,-10,10);
  //RooPolynomial polynom("polynom","Pol1",dxr,RooArgList(slope));

  RooRealVar s("s", "signal yield", 1000, 0, 10000);
  RooRealVar s2("s2", "signal yield", 1000, 0, 10000);
  RooRealVar b("b", "background yield", 100, 0, 200);
  RooAddPdf sum("sum","Total distribution",RooArgList(gauss,gauss2),RooArgList(s,s2));

   RooFitResult *r = sum.fitTo(matcheddata, RooFit::Save());
  //RooFitResult *r = sum.fitTo(binnedData, RooFit::Save());
  //summary printing
  r->Print();
  //TH2* hcorr = r->correlationHist() ;
  //plot the result
  RooPlot * xFrame = dxr.frame() ;
  //binnedData.plotOn(xFrame) ;
  matcheddata.plotOn(xFrame);
  sum.plotOn(xFrame) ;
  sum.plotOn(xFrame, RooFit::Components(gauss2), RooFit::LineStyle(kDashed)) ;
  TCanvas c;
  xFrame->Draw();
  c.SaveAs("polfit.png");
  cout << "chi^2 = " << xFrame->chiSquare(6) << endl;

  //redoing AGAIN fkit on y side

  TTree *goodtree2 = inputtree->CopyTree("TMath::Abs(scaled_dxr)<40");

  RooRealVar dyr("scaled_dyr", "match_dyr[#mum]", -0.025*cm2micron, 0.025*cm2micron);
  RooDataSet matcheddata_y("matcheddata","dataset with match_dyr",goodtree2,dyr) ;

  //define parameters and PDFs 
  RooRealVar mu_y("mu_y", "average", 0, -1000, 1000);
  RooRealVar sigma_y("sigma_y", "sigma",40, 0, 100);
  RooGaussian gauss_y("gauss_y","gaussian PDF", dyr, mu, sigma);

  RooRealVar mu2_y("mu2_y", "average of second gaussian", 0, -1000, 1000);
  RooRealVar sigma2_y("sigma2_y", "sigma of second gaussian", 20, 0, 100);
  RooGaussian gauss2_y("gauss2_y","Second gaussian PDF", dyr, mu2_y, sigma2_y);

  RooRealVar s_y("s_y", "signal yield", 1000, 0, 10000);
  RooRealVar s2_y("s2_y", "signal yield", 1000, 0, 10000);
  RooRealVar b_y("b_y", "background yield", 100, 0, 200);

  RooAddPdf sum_y("sum_y","Total distribution",RooArgList(gauss_y,gauss2_y),RooArgList(s_y,s2_y));

  RooFitResult *r_y = sum_y.fitTo(matcheddata_y, RooFit::Save());
  //RooFitResult *r = sum.fitTo(binnedData, RooFit::Save());
  //summary printing
  r_y->Print();
  //TH2* hcorr = r->correlationHist() ;
  //plot the result
  RooPlot * yFrame = dyr.frame() ;
  //binnedData.plotOn(xFrame) ;
  matcheddata_y.plotOn(yFrame);
  sum_y.plotOn(yFrame) ;
  sum_y.plotOn(yFrame, RooFit::Components(gauss2_y), RooFit::LineStyle(kDashed)) ;
  TCanvas c_y;
  yFrame->Draw();
  c_y.SaveAs("polfit_y.png");
  cout << "chi^2 = " << yFrame->chiSquare(6) << endl;
}

void scaletree(){
  ROOT::RDataFrame df("matched_trks","CH1R6_allspills.root");
  auto df1 = df.Define("scaled_dxr","match_dxr*10000")
               .Define("scaled_dyr","match_dyr*10000");
  
  df1.Snapshot("scaled_matched_tracks","scaledtree.root");
}


void anglefit(){
  
  gROOT->SetStyle("Plain");
  //define variables and importing data
  TFile *inputfile = TFile::Open("CH1R6_allspills.root");
  TTree *inputtree = (TTree*) inputfile->Get("matched_trks");

  TTree *goodtree = inputtree->CopyTree("TMath::Abs(match_dtyr)<0.004");

  RooRealVar dtxr("match_dtxr", "match_dtxr", -0.015, 0.015);
  RooDataSet matcheddata("matcheddata","dataset with match_dtxr",inputtree,dtxr) ;

  //adding binning, try binned
  //dtxr.setBins(30);
  //RooDataHist binnedData("binnedData", "matcheddata",
  //RooArgList(dtxr), matcheddata);

  //define parameters and PDFs 
  RooRealVar mu("mu", "average", 0, -1, 1);
  RooRealVar sigma("sigma", "sigma", 0.0004, 0, 1);
  RooGaussian gauss("gauss","gaussian PDF", dtxr, mu, sigma);

  //polynomial, max 1 degree
  RooRealVar slope("slope", "slope",20,-10,10);
  RooPolynomial polynom("polynom","Pol1",dtxr,RooArgList(slope));

  RooRealVar s("s", "signal yield", 1000, 0, 1000);
  RooRealVar b("b", "background yield", 100, 50, 200);
  RooAddPdf sum("sum","Total distribution",RooArgList(gauss,polynom),RooArgList(s,b));

  RooFitResult *r = sum.fitTo(matcheddata, RooFit::Save());
  //RooFitResult *r = sum.fitTo(binnedData, RooFit::Save());
  //summary printing
  r->Print();
  //TH2* hcorr = r->correlationHist() ;
  //plot the result
  RooPlot * xFrame = dtxr.frame() ;
  //binnedData.plotOn(xFrame) ;
  matcheddata.plotOn(xFrame);
  sum.plotOn(xFrame) ;
  sum.plotOn(xFrame, RooFit::Components(polynom), RooFit::LineStyle(kDashed)) ;
  TCanvas c;
  xFrame->Draw();
  //c.SaveAs("polfit.pdf");
  cout << "chi^2 = " << xFrame->chiSquare(5) << endl;

  //TCanvas cfitresults;
  //hcorr->Draw("colz");
}
