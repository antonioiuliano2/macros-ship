//#include "Tools/Flux/GSimpleNtpFlux.h"

genie::flux::GSimpleNtpEntry gsimpleentry(float E, float pdg, float px, float py, float pz, float weight){
    genie::flux::GSimpleNtpEntry e;
    e.E      = E;
    e.pdg    = pdg;
    e.px     = px;
    e.py     = py;
    e.pz     = pz;
    e.vtxx      = px/pz*400. * 1 / 100;
    e.vtxy      = py/pz*400. * 1 / 100;
    e.vtxz      = 400. * 1 / 100;
    e.dist = 0.;
    e.wgt = weight;
    e.metakey = TString("gsimple_output.root").Hash();
    e.metakey &= 0x7FFFFFFF;
    return e;
}

void rdataframe_decay2gsimple_converter(){

    const double pot_number = 5e+13; //reference of simulation weights (aka. POT for one spill)

    // Open your flat custom input TTree
    ROOT::RDataFrame df("Decay", "/eos/experiment/ship/user/aiuliano/nuhistos_bkgproductions/bkg2026/makeCascadeTungsten_2026_06_12/all100Runs_Decay_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root");
    auto dfminE = df.Filter("E > 1."); //filter out low energy neutrinos
    auto dfneutrinos = dfminE.Filter("TMath::Abs(id) == 12 || TMath::Abs(id) == 14 || TMath::Abs(id) == 16");
    // Pack your flat variables into a single GSimpleNtpEntry object per row
    auto df_gsimple = dfneutrinos.Define("entry", gsimpleentry,{"E", "id", "px", "py", "pz", "weight"});

    // Snapshot only the 'entry' branch into a new TTree called 'flux'
    df_gsimple.Snapshot("flux", "gsimple_output.root", {"entry"});

    ROOT::RDataFrame read_df("flux","gsimple_output.root");

    // processing meta data
    auto min_weight = read_df.Min("entry.wgt");
    auto max_weight = read_df.Max("entry.wgt");
    auto max_energy = read_df.Max("entry.E");

    auto min_x = read_df.Min("entry.vtxx");
    auto max_x = read_df.Max("entry.vtxx");
    auto min_y = read_df.Min("entry.vtxy");
    auto max_y = read_df.Max("entry.vtxy");
    auto min_z = read_df.Min("entry.vtxz");
    auto max_z = read_df.Max("entry.vtxz");

    //adding a TTree for metadata
    TFile *fOut = TFile::Open("gsimple_output.root", "UPDATE");
    TTree *metaOut = new TTree("meta", "metadata for flux n-tuple");
    genie::flux::GSimpleNtpMeta *meta_entry = new genie::flux::GSimpleNtpMeta;
    metaOut->Branch("meta", &meta_entry);

    meta_entry->metakey = TString("gsimple_output.root").Hash();
    meta_entry->metakey &= 0x7FFFFFFF;

    // plane corner and plane direction of the scoring plane (The scoring plane is tilted)
    double plane_corner[] = {*min_x, *min_y, *min_z};
    double plane_dir1[] = {*max_x - *min_x, 0, *max_z - *min_z};
    double plane_dir2[] = {0, *max_y - *min_y, *max_z - *min_z};

    for (int i = 0; i < 3; i++)
      meta_entry->windowBase[i] = plane_corner[i] * 1 / 100;
    for (int i = 0; i < 3; i++)
      meta_entry->windowDir1[i] = plane_dir1[i] * 1 / 100;
    for (int i = 0; i < 3; i++)
      meta_entry->windowDir2[i] = plane_dir2[i] * 1 / 100;

    meta_entry->maxEnergy = *max_energy;
    meta_entry->maxWgt = *max_weight;
    meta_entry->minWgt = *min_weight;

    meta_entry->protons = pot_number; // Number of protons on target.
    meta_entry->infiles.push_back("/eos/experiment/ship/user/aiuliano/nuhistos_bkgproductions/bkg2026/makeCascadeTungsten_2026_06_12/all100Runs_Decay_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root");
    meta_entry->seed = 0.;

    metaOut->Fill();
    fOut->cd();
    metaOut->Write();
    fOut->Close();
    std::cout << "Done converting" << std::endl;

}
