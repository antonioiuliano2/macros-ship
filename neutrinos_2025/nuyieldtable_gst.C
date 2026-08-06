//print table of nu yield from gst file
void nuyieldtable_gst(){
    TString prepath = "/home/antonio/Simulations/GENIE_sims/2026_07_03_4e19_allflavours_neutrino_detector/";
    ROOT::RDataFrame df("gst",prepath+TString("neutrino_detector_4e19.gst.root").Data());

    auto dfnumu = df.Filter("neu==14");

    auto dfcc = dfnumu.Filter("cc");
    auto dfdis = dfcc.Filter("dis");

    auto nu_ccdis = dfdis.Count();


    cout<<"NU:"<<"\t"<<"CCDIS"<<"\t"<<"CHARM"<<"\t"<<"RES"<<" "<<"QE"<<"\t"<<"NUEEL"<<"\t"<<"NCDIS"<<endl;

    cout<<"NUMU"<<"\t"<<nu_ccdis.GetValue()<<endl;
}