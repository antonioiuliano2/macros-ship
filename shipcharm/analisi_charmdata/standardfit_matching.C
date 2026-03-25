void standardfit_matching(){
    TFile *inputfile = TFile::Open("CH1R6_allspills.root");
    TTree *inputtree = (TTree*) inputfile->Get("matched_trks");

    TCanvas *cxr = new TCanvas();
    inputtree->Draw("match_dtxr>>htxr(30,-0.015,0.015)","TMath::Abs(match_dtyr)<0.004");
    inputtree->Draw("match_dtyr>>htyr(30,-0.015,0.015)","TMath::Abs(match_dtxr)<0.003");

    TH1D *htxr = (TH1D*) gROOT->FindObject("htxr");
    htxr->SetTitle(";#Delta#theta_{xz}[rad]");

    TH1D *htyr = (TH1D*) gROOT->FindObject("htyr");
    htyr->SetTitle(";#Delta#theta_{yz}[rad]");

    TF1 *sum= new TF1("sum","gaus(0) + pol1(3)",-0.015,0.015);

    sum->SetParameter(0,10000);
    sum->SetParameter(1,0);
    sum->SetParameter(2,0.004);
    sum->SetParameter(3,100);
    sum->SetParameter(4,-2);

    sum->SetParName(0,"GausConst");
    sum->SetParName(1,"Mean");
    sum->SetParName(2,"Sigma");
    sum->SetParName(3,"p0");
    sum->SetParName(4,"p1");

    htxr->Draw("E");
    htxr->Fit(sum,"L");

    //drawing pol component
    TF1 *polynom = new TF1("polynom","pol1",-0.015,0.015);
    polynom->SetParameter(0,sum->GetParameter("p0"));
    polynom->SetParameter(1,sum->GetParameter("p1"));

    polynom->SetLineColor(kBlue);
    polynom->Draw("SAME");

    TCanvas *cyr = new TCanvas();
    TF1 *sum1= new TF1("sum1","gaus(0) + pol1(3)",-0.015,0.015);

    sum1->SetParameter(0,1000);
    sum1->SetParameter(1,0);
    sum1->SetParameter(2,0.0035);
    sum1->SetParameter(3,100);
    sum1->SetParameter(4,-3);

    sum1->SetParName(0,"GausConst");
    sum1->SetParName(1,"Mean");
    sum1->SetParName(2,"Sigma");
    sum1->SetParName(3,"p0");
    sum1->SetParName(4,"p1");

    htyr->Draw("E");
    htyr->Fit(sum1,"L");

    //drawing pol component
    TF1 *polynom1 = new TF1("polynom1","pol1",-0.015,0.015);
    polynom1->SetParameter(0,sum1->GetParameter("p0"));
    polynom1->SetParameter(1,sum1->GetParameter("p1"));

    polynom1->SetLineColor(kBlue);
    polynom1->Draw("SAME");
}

void fits_spills(){
    const int nfiles = 5;
    int nspill[nfiles] = {8,9,10,11,12};

    //for automatic plots
    TCanvas *c0 = new TCanvas();
    //input files and trees
    TFile *inputfile[nfiles];
    TTree *inputtree[nfiles];

    //starting x section
    TCanvas *cx = new TCanvas();
    cx->Divide(2,3);

    TH1D *hxr[nfiles];

    TF1 *sumx[nfiles];
    TF1 *polynomx[nfiles];
    //list of parameters initial values

    double gausnormalizationx[nfiles] = {500,1000,10000,1000,1000};
    double gausmeanx[nfiles] = {0,0,0,0,0};
    double gaussigmax[nfiles] = {50,60,40,40,40};

    double p0x[nfiles] = {10,10,10,10,10};
    double p1x[nfiles] = {-1,-1,-1,-1,-1};

    gStyle->SetOptFit(1111);
    //loop over files
    for (int ifile = 0; ifile < nfiles; ifile++){
     inputfile[ifile] = TFile::Open(Form("CH1R6_Spill_%i.root",nspill[ifile]));
     inputtree[ifile] = (TTree*) inputfile[ifile] ->Get("matched_trks");

     c0->cd();
     //drawing histogram, applying a pre-selection on the other variable
     inputtree[ifile]->Draw(Form("match_dxr*1e+4>>hxr%i(25,-250,250)",ifile),"TMath::Abs(match_dyr)*1e+4<40");
     hxr[ifile] = (TH1D*) gROOT->FindObject(Form("hxr%i",ifile));

     sumx[ifile] = new TF1(Form("sumx%i",ifile),"gaus(0) + pol1(3)",-250,250);
     sumx[ifile]->SetParameter(0,gausnormalizationx[ifile]);
     sumx[ifile]->SetParameter(1,gausmeanx[ifile]);
     sumx[ifile]->SetParameter(2,gaussigmax[ifile]);
     sumx[ifile]->SetParameter(3,p0x[ifile]);
     sumx[ifile]->SetParameter(4,p1x[ifile]);

     sumx[ifile]->SetParName(0,"GausConst");
     sumx[ifile]->SetParName(1,"Mean");
     sumx[ifile]->SetParName(2,"Sigma");
     sumx[ifile]->SetParName(3,"p0");
     sumx[ifile]->SetParName(4,"p1");

     cx->cd(ifile+1);

     hxr[ifile]->SetTitle("x matched;x[#mum]");
     hxr[ifile]->Draw("E");
     hxr[ifile]->Fit(sumx[ifile] ,"L");

     polynomx[ifile] = new TF1(Form("polynomx%i",ifile),"pol1",-250,250);
     polynomx[ifile]->SetParameter(0,sumx[ifile] ->GetParameter("p0"));
     polynomx[ifile]->SetParameter(1,sumx[ifile] ->GetParameter("p1"));

     polynomx[ifile]->SetLineColor(kBlue);
     polynomx[ifile]->Draw("SAME");
    }

    //starting y section
    TCanvas *cy = new TCanvas();
    cy->Divide(2,3);

    TH1D *hyr[nfiles];
     
    TF1 *sumy[nfiles];
    TF1 *polynomy[nfiles];
    //list of parameters initial values

    double gausnormalizationy[nfiles] = {100,1000,1000,1000,1000};
    double gausmeany[nfiles] = {0,0,0,0,0};
    double gaussigmay[nfiles] = {70,60,40,40,40};

    double p0y[nfiles] = {10,10,10,10,10};
    double p1y[nfiles] = {-1,-1,-1,-1,-1};

    gStyle->SetOptFit(1111);
    //loop over files
    for (int ifile = 0; ifile < nfiles; ifile++){
     inputfile[ifile] = TFile::Open(Form("CH1R6_Spill_%i.root",nspill[ifile]));
     inputtree[ifile] = (TTree*) inputfile[ifile] ->Get("matched_trks");

     c0->cd();
     //drawing histogram, applying a pre-selection on the other variable
     inputtree[ifile]->Draw(Form("match_dyr*1e+4>>hyr%i(25,-250,250)",ifile),"TMath::Abs(match_dxr)*1e+4<40");
     hyr[ifile] = (TH1D*) gROOT->FindObject(Form("hyr%i",ifile));

     sumy[ifile] = new TF1(Form("sumy%i",ifile),"gaus(0) + pol1(3)",-250,250);
     sumy[ifile]->SetParameter(0,gausnormalizationy[ifile]);
     sumy[ifile]->SetParameter(1,gausmeany[ifile]);
     sumy[ifile]->SetParameter(2,gaussigmay[ifile]);
     sumy[ifile]->SetParameter(3,p0y[ifile]);
     sumy[ifile]->SetParameter(4,p1y[ifile]);

     sumy[ifile]->SetParName(0,"GausConst");
     sumy[ifile]->SetParName(1,"Mean");
     sumy[ifile]->SetParName(2,"Sigma");
     sumy[ifile]->SetParName(3,"p0");
     sumy[ifile]->SetParName(4,"p1");

     cy->cd(ifile+1);

     hyr[ifile]->SetTitle("y matched;y[#mum]");
     hyr[ifile]->Draw("E");
     hyr[ifile]->Fit(sumy[ifile] ,"L");

     polynomy[ifile] = new TF1(Form("polynomy%i",ifile),"pol1",-250,250);
     polynomy[ifile]->SetParameter(0,sumy[ifile] ->GetParameter("p0"));
     polynomy[ifile]->SetParameter(1,sumy[ifile] ->GetParameter("p1"));

     polynomy[ifile]->SetLineColor(kBlue);
     polynomy[ifile]->Draw("SAME");
    }

}