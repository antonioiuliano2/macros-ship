//print table of nu yield from gst file
void nuyieldtable_gst(){
    const double normsim = 4e+19; //one year of data taking simulated
    const double normship = 6e+20; //15 years of data taking
    double scale_factor = normship/normsim;

    TString prepath = "/home/antonio/Simulations/GENIE_sims/2026_08_27_4e19_allflavours_calorimeter/";
    ROOT::RDataFrame df("gst",prepath+TString("nuyield_4e19_calorimeter.gst.root").Data());
    //TString prepath = "/home/antonio/Simulations/GENIE_sims/2026_07_03_4e19_allflavours_calorimeter/";
    //ROOT::RDataFrame df("gst",prepath+TString("calorimeter_4e19.gst.root").Data());

    auto dfcc = df.Filter("cc");

    auto dfnue_cc = dfcc.Filter("neu==12");
    auto dfnumu_cc = dfcc.Filter("neu==14");
    auto dfnutau_cc = dfcc.Filter("neu==16");

    auto dfnue_bar_cc = dfcc.Filter("neu==-12");
    auto dfnumu_bar_cc = dfcc.Filter("neu==-14");
    auto dfnutau_bar_cc = dfcc.Filter("neu==-16");

    auto dfnc = df.Filter("nc");

    auto nue_cc = dfnue_cc.Count();
    auto numu_cc = dfnumu_cc.Count();
    auto nutau_cc = dfnutau_cc.Count();

    auto nue_bar_cc = dfnue_bar_cc.Count();
    auto numu_bar_cc = dfnumu_bar_cc.Count();
    auto nutau_bar_cc = dfnutau_bar_cc.Count();

    auto nu_nc = dfnc.Count();


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

    
}