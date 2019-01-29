void read_spill(){
    //building a RDataFrame from the inputfile
    ROOT::RDataFrame df("spill","../input/charm_spills.root");
    auto hpot = df.Histo1D("pot");
    hpot->DrawCopy();
}