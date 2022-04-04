
void fillTH1D_textfile(TH1D *histMC, TH1D*histdata, TString inputfilename){
    ifstream in;
    in.open(inputfilename.Data());
    int ipoint;
    float x, yMC, yDT, errMC, errDT;
    while (1) {
     in >> ipoint >> x >> yMC >> errMC >> yDT >> errDT;
     cout<<"TEST "<<ipoint<<" "<<x<<" "<<yMC<<endl;
     histMC->SetBinContent(histMC->FindBin(x),yMC);
     histMC->SetBinError(histMC->FindBin(x),errMC);
     histdata->SetBinContent(histdata->FindBin(x),yDT);
     histdata->SetBinError(histdata->FindBin(x),errDT);
     if (!in.good()) break;
    }
}

void roofit_rebin6(){
 const int npoints = 30;
 float xpoints[npoints+1] = {3,9,15,21,27,38,44,50,56,62,89,101,113,125,137,164,176,188,200,212,239,251,263,275,287,314,326,338,350,362,370};
 TH1D *hMC_full = new TH1D("hMC_full","Full MC histogram;vz[mm]",npoints,xpoints); //each bin is 6 mm thick
 TH1D *hdata = new TH1D("hdata","Data histogram;vz[mm]",npoints,xpoints); //each bin is 3 mm thick
 fillTH1D_textfile(hMC_full,hdata,TString("/home/utente/Lavoro/BDT_vertices_Valerio/afterBDT_plots/errors_filips.dat"));
// hMC_full->Draw();
// hdata->Draw("SAME");
 //instancing RooFit objects
 // Declare observable x
 RooRealVar vz("vz", "vz", 0, 366);
 RooDataHist dhMC_full("dhMC_full", "dhMC_full", vz, RooFit::Import(*hMC_full));

 
 //parameters of exponential
 RooRealVar p_N("p_N", "proton yield", 800, 0, 10000);
 RooRealVar p_c("p_c", "proton slope", -0.006031, -1, 1);
 //parameters of chebychev
 RooRealVar h_N("h_N", "hadronic yield", 100, 0, 10000);
 RooRealVar h_p1("h_p1","hadronic first degree",21.6,0,50);
 RooRealVar h_p2("h_p2","hadronic second degree",-0.09531,-1,1);
 RooRealVar h_p3("h_p3","hadronic third degree",0.00012,-0.1,0.1);

 //definining distributions and summing them
 RooExponential expo("expo","expo",vz,p_c);
 RooPolynomial pol3("pol3","Pol3",vz,RooArgList(h_p1,h_p2,h_p3));
 RooAddPdf sum("sum","Total distribution",RooArgList(expo,pol3),RooArgList(p_N,h_N));

 //fitting
 RooFitResult *r = sum.fitTo(dhMC_full, RooFit::Save()); //likelihood
 //RooFitResult *r = sum.chi2FitTo(dhMC_full, RooFit::Save()); //chi squared (actually uses input errors)
 //drawing results
 RooPlot *frame = vz.frame(RooFit::Title("Imported TH1 with internal errors"));
 dhMC_full.plotOn(frame, RooFit::DataError(RooAbsData::SumW2));
 sum.plotOn(frame,RooFit::Name("total")) ;
 TCanvas c_vz;
 frame->Draw();
 
 c_vz.SaveAs("test.png");
}