//plot as SHiP CDS radius neutrinos
using namespace ROOT;
void drawRadiusPlot(){
    double normsim = 5e+13; //reference of simulation weights (aka. POT for one spill)
    double normship = 4e+19; //reference for one year of SND DataTaking
    //listing results for different radii
    RVec<double> Radius = {0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,1.0};
    RVec<string> FolderNames = {"0_1","0_2","0_3","0_4","0_5","0_6","0_7","0_8","0_9","1_0"};
    RVec<double> dL = Radius/TMath::Sqrt(2);

    RVec<double> nue_ccdis, numu_ccdis, nutau_ccdis, nue_bar_charmccdis, numu_bar_charmccdis;


    TString prepath("/home/utente/Simulations/nuyield_shipecn3/nuyield_AdvSND_MuonShieldTarget/tungstentarget/");
    TString posfolder = ("pos1_noSC/");
    TString resultsfolder("/plots_2/");
    //where we save our plots
    TString outputfolder(
        "/home/utente/cernbox/Synched/Archivio_cronologico/Settembre_2024/MuonShield_SND_studies/tungstentarget_radiusplots/");

    for (auto& foldername: FolderNames){
        TString fullpath = prepath + posfolder + TString("radius_")+TString(foldername) + resultsfolder;
        //opening files
        TFile *file_nue_ccdis = TFile::Open(fullpath+TString("results_nu_e_dis_cc.root"));
        TFile *file_numu_ccdis = TFile::Open(fullpath+TString("results_nu_mu_dis_cc.root"));
        TFile *file_nutau_ccdis = TFile::Open(fullpath+TString("results_nu_tau_dis_cc.root"));
        TFile *file_nue_bar_charmccdis = TFile::Open(fullpath+TString("results_nu_e_bar_dis_cc_charm.root"));
        TFile *file_numu_bar_charmccdis = TFile::Open(fullpath+TString("results_nu_mu_bar_dis_cc_charm.root"));

        //getting histograms
        TH1D * hspectrum_nu_e_intdis_cc = (TH1D*) file_nue_ccdis->Get("hspectrum_nu_e_intdis_cc");
        TH1D * hspectrum_nu_mu_intdis_cc = (TH1D*) file_numu_ccdis->Get("hspectrum_nu_mu_intdis_cc");
        TH1D * hspectrum_nu_tau_intdis_cc = (TH1D*) file_nutau_ccdis->Get("hspectrum_nu_tau_intdis_cc");
        TH1D * hspectrum_nu_e_bar_intdis_cc_charm = (TH1D*) file_nue_bar_charmccdis->Get("hspectrum_nu_e_bar_intdis_cc_charm");
        TH1D * hspectrum_nu_mu_bar_intdis_cc_charm = (TH1D*) file_numu_bar_charmccdis->Get("hspectrum_nu_mu_bar_intdis_cc_charm");
        
        //returning integrals and appending them to vectors
        nue_ccdis.push_back(hspectrum_nu_e_intdis_cc->Integral());
        numu_ccdis.push_back(hspectrum_nu_mu_intdis_cc->Integral());
        nutau_ccdis.push_back(hspectrum_nu_tau_intdis_cc->Integral());
        nue_bar_charmccdis.push_back(hspectrum_nu_e_bar_intdis_cc_charm->Integral());
        numu_bar_charmccdis.push_back(hspectrum_nu_mu_bar_intdis_cc_charm->Integral());

        //closing files
        file_nue_ccdis->Close();
        file_numu_ccdis->Close();
        file_nutau_ccdis->Close();
        file_nue_bar_charmccdis->Close();
        file_numu_bar_charmccdis->Close();
    }
    /*
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
 */
 cout<<"Checking values for safety "<<endl;
 cout<<"Nu e CCDIS"<<endl;
 cout<<nue_ccdis<<endl;
 cout<<"Nu mu CCDIS"<<endl;
 cout<<numu_ccdis<<endl;
 cout<<"Nu tau CCDIS"<<endl;
 cout<<nutau_ccdis<<endl;
 cout<<"Nu e bar CharmCCDIS"<<endl;
 cout<<nue_bar_charmccdis<<endl;
 cout<<"Nu mu bar CharmCCDIS"<<endl;
 cout<<numu_bar_charmccdis<<endl;

 TGraph *gnue_ccdis = new TGraph(10, (dL).data(),nue_ccdis.data());
 TGraph *gnumu_ccdis = new TGraph(10, (dL).data(),numu_ccdis.data());
 TGraph *gnutau_ccdis = new TGraph(10, (dL).data(),nutau_ccdis.data());

 TGraph *gnue_bar_charmccdis = new TGraph(10, (dL).data(),nue_bar_charmccdis.data());
 TGraph *gnumu_bar_charmccdis = new TGraph(10, (dL).data(),numu_bar_charmccdis.data());

 gnue_ccdis->SetTitle("Electron neutrino yield;dL[m]");
 gnumu_ccdis->SetTitle("Muon neutrino yield;dL[m]");
 gnutau_ccdis->SetTitle("Tau neutrino yield;dL[m]");

 gnue_bar_charmccdis->SetTitle("Electron neutrino yield;dL[m]");
 gnumu_bar_charmccdis->SetTitle("Muon neutrino yield;dL[m]");
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

 gnumu_ccdis->SetTitle("CCDIS Yields with same mass detector of different transverse size");

 
 cgraph->Print(outputfolder+posfolder+"nuyields_size.root");
 cgraph->Print(outputfolder+posfolder+"nuyields_size.png");
}