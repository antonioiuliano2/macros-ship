//plot as SHiP CDS radius neutrinos
using namespace ROOT;
void drawRadiusPlot(){
    double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
    double normship = 2e+20; //reference for five years of SND DataTaking
    //listing results for different radii
    RVec<double> Radius = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};

    RVec<double> nue_ccdis = 
    {
        0.051356,
        0.0423493,
        0.0368452,
        0.0320188,
        0.0279258,
        0.0244371,
        0.0214997,
        0.0190111,
        0.0169235,
        0.0151109
    
    };
    RVec<double> numu_ccdis = 
    {
        0.219229,
        0.145314,
        0.113255,
        0.0925719,
        0.0780705,
        0.0668418,
        0.0580342,
        0.0509402,
        0.0451708,
        0.0403567
    };
    RVec<double> nutau_ccdis = 
    {
        0.00155321,
        0.0011463,
        0.0010681,
        0.000966631,
        0.000886596,
        0.000798318,
        0.000720012,
        0.00065394,
        0.000592937,
        0.00054009


    };

    RVec<double> nue_bar_charmccdis = 
    {  
        0.000498312,
        0.000431799,
        0.000383625,
        0.000348028,
        0.000313134,
        0.000281052,
        0.000252444,
        0.000226453,
        0.000203516,
        0.000183422


    };
    RVec<double> numu_bar_charmccdis = 
    {
        0.00157177,
        0.00117657,
        0.000959213,
        0.000793944,
        0.000666736,
        0.000569545,
        0.000492098,
        0.000427881,
        0.000375147

    };
 
 TGraph *gnue_ccdis = new TGraph(10, Radius.data(),nue_ccdis.data());
 TGraph *gnumu_ccdis = new TGraph(10, Radius.data(),numu_ccdis.data());
 TGraph *gnutau_ccdis = new TGraph(10, Radius.data(),nutau_ccdis.data());

 TGraph *gnue_bar_charmccdis = new TGraph(10, Radius.data(),nue_bar_charmccdis.data());
 TGraph *gnumu_bar_charmccdis = new TGraph(10, Radius.data(),numu_bar_charmccdis.data());

 gnue_ccdis->SetTitle("Electron neutrino yield;R[m]");
 gnumu_ccdis->SetTitle("Muon neutrino yield;R[m]");
 gnutau_ccdis->SetTitle("Tau neutrino yield;R[m]");

 gnue_bar_charmccdis->SetTitle("Electron neutrino yield;R[m]");
 gnumu_bar_charmccdis->SetTitle("Muon neutrino yield;R[m]");
 //normalizing to ship data taking
 gnue_ccdis->Scale(normship/normsim);
 gnumu_ccdis->Scale(normship/normsim);
 gnutau_ccdis->Scale(normship/normsim);

 gnue_bar_charmccdis->Scale(normship/normsim);
 gnumu_bar_charmccdis->Scale(normship/normsim);  
 TCanvas *cgraph = new TCanvas();
 gnumu_ccdis->SetMarkerColor(kRed);
 gnumu_ccdis->SetMarkerStyle(20);
 gnumu_ccdis->Draw("AP");
 gnumu_ccdis->GetYaxis()->SetRangeUser(1,5e+6);
 gnue_ccdis->SetMarkerColor(kBlue);
 gnue_ccdis->SetMarkerStyle(32);
 gnue_ccdis->Draw("P && SAME");
 gnutau_ccdis->SetMarkerColor(kBlack);
 gnutau_ccdis->Draw("P && SAME"); 
 gnutau_ccdis->SetMarkerStyle(21);
 cgraph->BuildLegend();
 cgraph->SetLogy();

 gnumu_ccdis->SetTitle("CCDIS Yields with mass of 1 ton square detector of different R");
}