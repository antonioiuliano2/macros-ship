//Exercise: after generating with an exponential, fit 
 {
  
  gROOT->SetStyle("Plain");
  //define variables, parameters and PDFs 
  RooRealVar x("x", "x", -10, 10);

  //PDFs for generation
  RooRealVar lambda("lambda", "slope", -0.1, -5., 0.);
  RooExponential expo("expo", "exponential PDF", x, lambda);

  //generate the Data
  RooDataSet * data = expo.generate(x, 20000);
  //fit the Data to the PDF
  RooFitResult *r = expo.fitTo(*data, RooFit::Save());
  //summary printing
  r->Print();
  //plot the result
  RooPlot * xFrame = x.frame() ;
  data->plotOn(xFrame) ;
  expo.plotOn(xFrame) ;
  //polsum.plotOn(xFrame, RooFit::Components(expo), RooFit::LineStyle(kDashed)) ;
  TCanvas c;
  xFrame->Draw();
  c.SaveAs("expofit.pdf")
}
