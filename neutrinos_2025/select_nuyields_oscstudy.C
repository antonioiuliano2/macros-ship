/*
Expected neutrino yields in the SHiP Scattering Detector (SD) and in the Hidden Sector Calorimeter (Calorimeter)
After applying selections for the combined acceptance for the neutrino oscillation paper.
To launch it:
 root -l
  >> .L select_nuyields_oscstudy.C
  >> SD_Selected()
  >> Calo_Selected()
Required access to SHiP EOS experiment filespace for input files.
*/
bool angular_acceptance(double vx, double vy, double vz){
    //angular threshold, obtained by 40 x 40 square centimeters at first iron section (see Section II.D)
    const double thetamax = 0.006886;
    //const double thetamax = 0.0100; //full MTC
    const double zstart = 0.;
    double dz = vz - zstart;
    if (TMath::Abs(vx) > dz * thetamax) return false;
    if (TMath::Abs(vy) > dz * thetamax) return false;
    //both checks ok, vertex in fiducial angular acceptance
    return true;

}

bool max_vz_calo(double vz){
    //require maximum z of interaction vertex for acceptance
    const double maxvz = 98.61; //m, GENIE uses SI Units, reference system SHiP geometry as gdml file
    if (vz > maxvz) return false;
    else return true;
}
bool max_vz_SD(double vz){
    //require maximum z of interaction vertex for acceptance
    const double maxvz = 28.37; //m, GENIE uses SI Units, reference system SHiP geometry as gdml file
    if (vz > maxvz) return false;
    else return true;
}

void SD_Selected(){
    //selected neutrino yields in the SHiP Scattering Detector (SD)
    //rescale from simulation normalization to total data taking expected in SHiP
    const double normsim = 4e+19; //one year of data taking simulated
    const double normship = 6e+20; //15 years of data taking
    double scale_factor = normship/normsim;
    //open GST Tree file from GENIE (converted from GHEP native format with gntpc utility)
    TString prepath = "/eos/experiment/ship/user/aiuliano/GENIE_sims/2026_08_27_4e19_allflavours_neutrino_detector/";
    ROOT::RDataFrame df("gst",prepath+TString("neutrino_detector_4e19.gst.root").Data());

    TFile *hfile = new TFile("histo_neutrinos_SD.root","RECREATE");
    //apply zcut
    auto df0 = df.Filter(max_vz_SD,{"vtxz"});
    //apply angular acceptance
    auto df1  = df0.Filter(angular_acceptance,{"vtxx","vtxy","vtxz"});

    auto hE_all = df.Filter("cc&&TMath::Abs(neu)==14").Histo1D({"hE_all","Energy of CC muon neutrinos;E[GeV]",75,0.,150.},"Ev");
    auto hE_selected = df1.Filter("cc&&TMath::Abs(neu)==14").Histo1D({"hE_selected","Energy of CC selected muon neutrinos;E[GeV]",75,0.,150.},"Ev");
    //select charged current events, and separate by neutrino flavour
    auto dfcc = df1.Filter("cc");

    auto dfnue_cc = dfcc.Filter("neu==12");
    auto dfnumu_cc = dfcc.Filter("neu==14");
    auto dfnutau_cc = dfcc.Filter("neu==16");

    auto dfnue_bar_cc = dfcc.Filter("neu==-12");
    auto dfnumu_bar_cc = dfcc.Filter("neu==-14");
    auto dfnutau_bar_cc = dfcc.Filter("neu==-16");
    //select neutral current events, we cannot identify flavour here
    auto dfnc = df1.Filter("nc");
    //we also need to select number of nc with pi^0
    auto dfnc_pi0 = dfnc.Filter("nfpi0>0");
    //retrieve the number of events for each category
    auto nue_cc = dfnue_cc.Count();
    auto numu_cc = dfnumu_cc.Count();
    auto nutau_cc = dfnutau_cc.Count();

    auto nue_bar_cc = dfnue_bar_cc.Count();
    auto numu_bar_cc = dfnumu_bar_cc.Count();
    auto nutau_bar_cc = dfnutau_bar_cc.Count();

    auto nu_nc = dfnc.Count();
    auto nu_nc_pi0 = dfnc_pi0.Count();

    //print the results
    cout<<"Scale factor: "<<scale_factor<<endl;

    cout<<"NU:"<<"\t"<< "CC"<<endl;
    cout<<"NUE"<<"\t"<<nue_cc.GetValue() * scale_factor<<endl;
    cout<<"NUMU"<<"\t"<<numu_cc.GetValue() * scale_factor<<endl;
    cout<<"NUTAU"<<"\t"<<nutau_cc.GetValue() * scale_factor<<endl;

    cout<<"ANTI_NUE"<<"\t"<<nue_bar_cc.GetValue() * scale_factor<<endl;
    cout<<"ANTI_NUMU"<<"\t"<<numu_bar_cc.GetValue() * scale_factor<<endl;
    cout<<"ANTI_NUTAU"<<"\t"<<nutau_bar_cc.GetValue() * scale_factor<<endl;

    cout<<"NUE + ANTI_NUE"<<"\t"<<(nue_cc.GetValue() + nue_bar_cc.GetValue()) * scale_factor<<endl;
    cout<<"NUMU + ANTI_NUMU"<<"\t"<<(numu_cc.GetValue() + numu_bar_cc.GetValue()) * scale_factor<<endl;
    cout<<"NUTAU + ANTI_NUTAU"<<"\t"<<(nutau_cc.GetValue() + nutau_bar_cc.GetValue()) * scale_factor<<endl;

    cout<<"ALL NC"<<"\t"<<nu_nc.GetValue() * scale_factor<<endl;
    cout<<"ALL NC with pi0"<<"\t"<<nu_nc_pi0.GetValue() * scale_factor<<endl;

}

void Calorimeter_Selected(){
    //selected neutrino yields in the SHiP Hidden Sector Calorimeter
    const double normsim = 4e+19; //one year of data taking simulated
    const double normship = 6e+20; //15 years of data taking
    double scale_factor = normship/normsim;
    //open GST Tree file from GENIE (converted from GHEP native format with gntpc utility)
    TString prepath = "/eos/experiment/ship/user/aiuliano/GENIE_sims/2026_08_27_4e19_allflavours_calorimeter/";
    ROOT::RDataFrame df("gst",prepath+TString("nuyield_4e19_calorimeter.gst.root").Data());

    TFile *hfile = new TFile("histo_neutrinos_calorimeter.root","RECREATE");
    //apply zcut
    auto df0 = df.Filter(max_vz_calo,{"vtxz"});
    //apply angular acceptance
    auto df1  = df0.Filter(angular_acceptance,{"vtxx","vtxy","vtxz"});

    auto hE_all = df.Filter("cc&&TMath::Abs(neu)==14").Histo1D({"hE_all","Energy of CC muon neutrinos;E[GeV]",75,0.,150.},"Ev");
    auto hE_selected = df1.Filter("cc&&TMath::Abs(neu)==14").Histo1D({"hE_selected","Energy of CC selected muon neutrinos;E[GeV]",75,0.,150.},"Ev");

    //select charged current events, and separate by neutrino flavour
    auto dfcc = df1.Filter("cc");

    auto dfnue_cc = dfcc.Filter("neu==12");
    auto dfnumu_cc = dfcc.Filter("neu==14");
    auto dfnutau_cc = dfcc.Filter("neu==16");

    auto dfnue_bar_cc = dfcc.Filter("neu==-12");
    auto dfnumu_bar_cc = dfcc.Filter("neu==-14");
    auto dfnutau_bar_cc = dfcc.Filter("neu==-16");

    //select neutral current events, we cannot identify flavour here
    auto dfnc = df1.Filter("nc");
    //we also need to select number of nc with pi^0
    auto dfnc_pi0 = dfnc.Filter("nfpi0>0");
    //retrieve the number of events for each category
    auto nue_cc = dfnue_cc.Count();
    auto numu_cc = dfnumu_cc.Count();
    auto nutau_cc = dfnutau_cc.Count();

    auto nue_bar_cc = dfnue_bar_cc.Count();
    auto numu_bar_cc = dfnumu_bar_cc.Count();
    auto nutau_bar_cc = dfnutau_bar_cc.Count();

    auto nu_nc = dfnc.Count();
    auto nu_nc_pi0 = dfnc_pi0.Count();

    //print the results
    cout<<"Scale factor: "<<scale_factor<<endl;

    cout<<"NU:"<<"\t"<< "CC"<<endl;
    cout<<"NUE"<<"\t"<<nue_cc.GetValue() * scale_factor<<endl;
    cout<<"NUMU"<<"\t"<<numu_cc.GetValue() * scale_factor<<endl;
    cout<<"NUTAU"<<"\t"<<nutau_cc.GetValue() * scale_factor<<endl;

    cout<<"ANTI_NUE"<<"\t"<<nue_bar_cc.GetValue() * scale_factor<<endl;
    cout<<"ANTI_NUMU"<<"\t"<<numu_bar_cc.GetValue() * scale_factor<<endl;
    cout<<"ANTI_NUTAU"<<"\t"<<nutau_bar_cc.GetValue() * scale_factor<<endl;

    cout<<"NUE + ANTI_NUE"<<"\t"<<(nue_cc.GetValue() + nue_bar_cc.GetValue()) * scale_factor<<endl;
    cout<<"NUMU + ANTI_NUMU"<<"\t"<<(numu_cc.GetValue() + numu_bar_cc.GetValue()) * scale_factor<<endl;
    cout<<"NUTAU + ANTI_NUTAU"<<"\t"<<(nutau_cc.GetValue() + nutau_bar_cc.GetValue()) * scale_factor<<endl;

    cout<<"ALL NC"<<"\t"<<nu_nc.GetValue() * scale_factor<<endl;
    cout<<"ALL NC with pi0"<<"\t"<<nu_nc_pi0.GetValue() * scale_factor<<endl;

}
