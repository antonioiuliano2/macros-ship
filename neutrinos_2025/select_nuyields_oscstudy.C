
bool angular_acceptance(double vx, double vy, double vz){
    //angular threshold, obtained by 40 x 40 square centimeters at first iron section
    const double thetamax = 0.00686;
    //const double thetamax = 0.0100; //full MTC
    const double zstart = 0.;
    double dz = vz - zstart;
    if (TMath::Abs(vx) > dz * thetamax) return false;
    if (TMath::Abs(vy) > dz * thetamax) return false;
    //both checks ok, vertex in fiducial angular acceptance
    return true;

}

bool max_vz_calo(double vz){
    //for now, only electromagnetic part
    const double maxvz = 98.61;
    if (vz > maxvz) return false;
    else return true;
}
bool max_vz_SD(double vz){
    //for now, only electromagnetic part
    const double maxvz = 28.37;
    if (vz > maxvz) return false;
    else return true;
}

void SD_Selected(){
    const double normsim = 4e+19; //one year of data taking simulated
    const double normship = 6e+20; //15 years of data taking
    double scale_factor = normship/normsim;

    TString prepath = "/eos/experiment/ship/user/aiuliano/GENIE_sims/2026_08_27_4e19_allflavours_neutrino_detector/";
    ROOT::RDataFrame df("gst",prepath+TString("neutrino_detector_4e19.gst.root").Data());

    auto df0 = df.Filter(max_vz_SD,{"vtxz"});
    auto df1  = df0.Filter(angular_acceptance,{"vtxx","vtxy","vtxz"});

    auto hzy_all = df.Histo2D({"hzy_all","zy of vertices;vz[m];vy[m]",100,26.,32.,80,-4.,4.},"vtxz","vtxy");
    auto hzy_selected = df1.Histo2D({"hzy_selected","zy of selected vertices;vz[m];vy[m]",100,26.,32.,80,-4.,4.},"vtxz","vtxy");

    //everything else remains the same
    auto dfcc = df1.Filter("cc");

    auto dfnue_cc = dfcc.Filter("neu==12");
    auto dfnumu_cc = dfcc.Filter("neu==14");
    auto dfnutau_cc = dfcc.Filter("neu==16");

    auto dfnue_bar_cc = dfcc.Filter("neu==-12");
    auto dfnumu_bar_cc = dfcc.Filter("neu==-14");
    auto dfnutau_bar_cc = dfcc.Filter("neu==-16");

    auto dfnc = df1.Filter("nc");
    //we also need nc with pi^0
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
    //cout<<"NU:"<<"\t"<<"CCDIS"<<"\t"<<"CHARM"<<"\t"<<"RES"<<" "<<"QE"<<"\t"<<"NUEEL"<<"\t"<<"NCDIS"<<endl;
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

    TCanvas *cyz = new TCanvas();
    hzy_selected->SetMarkerColor(kRed);
    hzy_all->DrawClone("SCAT");
    hzy_selected->DrawClone("SCAT&&SAME");
    cyz->Draw();


}

void Calorimeter_Selected(){
    const double normsim = 4e+19; //one year of data taking simulated
    const double normship = 6e+20; //15 years of data taking
    double scale_factor = normship/normsim;

    TString prepath = "/eos/experiment/ship/user/aiuliano/GENIE_sims/2026_08_27_4e19_allflavours_calorimeter/";
    ROOT::RDataFrame df("gst",prepath+TString("nuyield_4e19_calorimeter.gst.root").Data());

    auto df0 = df.Filter(max_vz_calo,{"vtxz"});
    auto df1  = df0.Filter(angular_acceptance,{"vtxx","vtxy","vtxz"});

    auto hzy_all = df.Histo2D({"hzy_all","zy of selected vertices",100,96.,100.,80,-4.,4.},"vtxz","vtxy");
    auto hzy_selected = df1.Histo2D({"hzy_selected","zy of selected vertices",100,96.,100.,80,-4.,4.},"vtxz","vtxy");

    //everything else remains the same
    auto dfcc = df1.Filter("cc");

    auto dfnue_cc = dfcc.Filter("neu==12");
    auto dfnumu_cc = dfcc.Filter("neu==14");
    auto dfnutau_cc = dfcc.Filter("neu==16");

    auto dfnue_bar_cc = dfcc.Filter("neu==-12");
    auto dfnumu_bar_cc = dfcc.Filter("neu==-14");
    auto dfnutau_bar_cc = dfcc.Filter("neu==-16");

    auto dfnc = df1.Filter("nc");
    //we also need nc with pi^0
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
    //cout<<"NU:"<<"\t"<<"CCDIS"<<"\t"<<"CHARM"<<"\t"<<"RES"<<" "<<"QE"<<"\t"<<"NUEEL"<<"\t"<<"NCDIS"<<endl;
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

    TCanvas *cyz = new TCanvas();
    hzy_selected->SetMarkerColor(kRed);
    hzy_all->DrawClone("SCAT");
    hzy_selected->DrawClone("SCAT&&SAME");
    cyz->Draw();


}