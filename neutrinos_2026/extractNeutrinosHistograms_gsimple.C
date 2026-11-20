
using namespace ROOT;
void extractNeutrinosHistograms_gsimple(){
    TString outfilename("neutrinos_histograms_mbiasnocharm_gsimple.root");
    gStyle->SetOptStat("nmri");
    //reading input file, getting dataframe
    TString prepath("/home/antonio/Simulations/inputs_gevgen_fnal/flux_gsimple_2026_10_22_finaltargetproduction/");
    TString filename("shuffled_Mbias2026_finaltarget_Cascade2026P8_merged_nuflux_gsimple.root");
    auto dfmeta = RDataFrame("meta", (prepath+filename).Data());
    dfmeta.Display("protons")->Print();
    auto pot_number = dfmeta.Take<double>("protons").GetValue()[0];

    double normpot = 5e+13; //reference of simulation weights (aka. POT for one spill)

    auto df = RDataFrame("flux", (prepath+filename).Data());
    auto df_mbias = df.Filter("wgt<1.001 && wgt>0.999");
    //auto df_cascade = df.Filter("wgt>1.001");
    //hard cut requiring E>1 GeV for mbias, no sense looking at lower energies
    //auto df_cut = df.Filter("E>1");
    auto df_cut = df_mbias.Filter("E>1");
    //auto df_cut = df_cascade.Filter("E>1");
    
    //computing p and pt from 3-momentum components, and log10p and log10pt
    auto df1 = df_cut.Define("p", "sqrt(px*px + py*py + pz*pz)").Define("pt", "sqrt(px*px + py*py)");
    auto df2 = df1.Define("log10p", "log10(p)").Define("log10pt", "log10(pt)");
    
    
    ROOT::RDF::TH1DModel hnu_P_model("hnu_P", "Neutrino momentum;P [GeV/c];Entries", 400, 0, 400);
    ROOT::RDF::TH2DModel hnu_log10P_log10Pt_model("hnu_log10P_log10Pt", "Neutrino log10pt vs log10p;log10p;log10pt", 
        100, -0.3, 1.7, 100, -2, 0.5);
    
    //make histograms for nue (do not forget weights)
    auto dfnue = df2.Filter("pdg==12");
    hnu_P_model.fName = "1012";
    auto hnue_P = dfnue.Histo1D(hnu_P_model, "p", "wgt");
    hnu_log10P_log10Pt_model.fName = "1212";
    auto hnue_log10P_log10Pt = dfnue.Histo2D(hnu_log10P_log10Pt_model, "log10p", "log10pt", "wgt");

    //make histograms for numu (do not forget weights)
    auto dfnumu = df2.Filter("pdg==14");
    hnu_P_model.fName = "1014";
    auto hnumu_P = dfnumu.Histo1D(hnu_P_model, "p", "wgt");  
    hnu_log10P_log10Pt_model.fName = "1214";
    auto hnumu_log10P_log10Pt = dfnumu.Histo2D(hnu_log10P_log10Pt_model, "log10p", "log10pt", "wgt");

    //make histograms for nutau (do not forget weights)
    auto dfnutau = df2.Filter("pdg==16");
    hnu_P_model.fName = "1016";
    auto hnutau_P = dfnutau.Histo1D(hnu_P_model, "p", "wgt");
    hnu_log10P_log10Pt_model.fName = "1216";
    auto hnutau_log10P_log10Pt = dfnutau.Histo2D(hnu_log10P_log10Pt_model, "log10p", "log10pt", "wgt");

    //make histograms for nue_bar (do not forget weights)
    auto dfnue_bar = df2.Filter("pdg==-12");
    hnu_P_model.fName = "2012";
    auto hnue_bar_P = dfnue_bar.Histo1D(hnu_P_model, "p", "wgt");
    hnu_log10P_log10Pt_model.fName = "2212";
    auto hnue_bar_log10P_log10Pt = dfnue_bar.Histo2D(hnu_log10P_log10Pt_model, "log10p", "log10pt", "wgt");

    //make histograms for numu (do not forget weights)
    auto dfnumu_bar = df2.Filter("pdg==-14");
    hnu_P_model.fName = "2014";
    auto hnumu_bar_P = dfnumu_bar.Histo1D(hnu_P_model, "p", "wgt");  
    hnu_log10P_log10Pt_model.fName = "2214";
    auto hnumu_bar_log10P_log10Pt = dfnumu_bar.Histo2D(hnu_log10P_log10Pt_model, "log10p", "log10pt", "wgt");

    //make histograms for nutau (do not forget weights)
    auto dfnutau_bar = df2.Filter("pdg==-16");
    hnu_P_model.fName = "2016";
    auto hnutau_bar_P = dfnutau_bar.Histo1D(hnu_P_model, "p", "wgt");
    hnu_log10P_log10Pt_model.fName = "2216";
    auto hnutau_bar_log10P_log10Pt = dfnutau_bar.Histo2D(hnu_log10P_log10Pt_model, "log10p", "log10pt", "wgt");

    //renormalize histograms to POT per spill (5e+13)
    hnue_P->Scale(normpot/pot_number);
    hnue_log10P_log10Pt->Scale(normpot/pot_number);
    
    hnumu_P->Scale(normpot/pot_number);
    hnumu_log10P_log10Pt->Scale(normpot/pot_number);
    
    hnutau_P->Scale(normpot/pot_number);
    hnutau_log10P_log10Pt->Scale(normpot/pot_number);

    hnue_bar_P->Scale(normpot/pot_number);
    hnue_bar_log10P_log10Pt->Scale(normpot/pot_number);
    
    hnumu_bar_P->Scale(normpot/pot_number);
    hnumu_bar_log10P_log10Pt->Scale(normpot/pot_number);
    
    hnutau_bar_P->Scale(normpot/pot_number);
    hnutau_bar_log10P_log10Pt->Scale(normpot/pot_number);

    //check with a TCanvas for debugging
    TCanvas *cnumu_P = new TCanvas("cnumu_P", "Muon neutrino momentum");
    hnumu_P->DrawClone();
    cnumu_P->SetLogy();

    TCanvas *cnumu_log10P_log10Pt = new TCanvas("cnumu_log10P_log10Pt", "Muon neutrino log10pt vs log10p");
    hnumu_log10P_log10Pt->DrawClone("COLZ");
    
    TFile *outfile = new TFile((prepath+outfilename).Data(), "RECREATE");
    hnue_P->Write();
    hnue_log10P_log10Pt->Write();
    hnumu_P->Write();
    hnumu_log10P_log10Pt->Write();
    hnutau_P->Write();
    hnutau_log10P_log10Pt->Write();

    hnue_bar_P->Write();
    hnue_bar_log10P_log10Pt->Write();
    hnumu_bar_P->Write();
    hnumu_bar_log10P_log10Pt->Write();
    hnutau_bar_P->Write();
    hnutau_bar_log10P_log10Pt->Write();

}