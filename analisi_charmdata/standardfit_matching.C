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