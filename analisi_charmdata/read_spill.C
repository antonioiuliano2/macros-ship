void read_spill(){
    //building a RDataFrame from the inputfile
    ROOT::RDataFrame df("spill","charm_spills.root");
    
    TCanvas *c1 = new TCanvas();
    auto hemu = df.Histo1D("nemu");
    hemu->DrawCopy();
    TCanvas *c2 = new TCanvas();
    auto hpassive = df.Histo1D("npassive");
    hpassive->DrawCopy();

    //fraction of protons not interacting
    auto survivingpot = [] (int nemu, int npassive,  int name){
           //units in mm
        const float filmthickness = 0.300;
        const float emulambda = 30000; //stakeholder, need to find it 
        const float pblambda = 152.5;
        float totalemuthickness = nemu * filmthickness;
        float nlambda = totalemuthickness/emulambda + npassive/pblambda; //computing the total number of interaction lengths
        float surviving = TMath::Exp(-nlambda); //surviving protons (averaging in the N spills)
        return surviving; //percentage of surviving pot
    };
    //fraction of protons interacting in ECC
    auto signalpot = [] (int nemu, int npassive, int name){
        const float filmthickness = 0.300;
        const float emulambda = 30000; //stakeholder, need to find it 
        const float pblambda = 152.5;
        float totalemuthickness = nemu * filmthickness;
        float nlambda = totalemuthickness/emulambda + npassive/pblambda; //computing the total number of interaction lengths
        float preshowerthickness = 0.;

        if (name/10 == 0) return preshowerthickness;       
        if (name/10 == 2) preshowerthickness = 28;
        else if (name/10 > 2) preshowerthickness = 56 *(name/10-2);

        float nlambdapassive = preshowerthickness/pblambda;
        float inpreshower = (1-TMath::Exp(-nlambdapassive));

        float surviving = TMath::Exp(-nlambda); //surviving protons (averaging in the N spills)
        float inECC = 1 - surviving - inpreshower;
        return inECC; //percentage of pot interacting within the brick
    };

    auto entries = df.Filter("name/10==1").Count();
    cout<<"CHARM1 number of runs: "<<*entries<<endl;

    TCanvas *csurviving = new TCanvas();
    auto df1 = df.Define("survpercent",survivingpot,{"nemu","npassive","name"});
    auto gpot = df1.Histo2D({"h2","Surviving fraction of pot;charmrun",10,0,100,100,0,1},"name","survpercent");
    gpot->DrawCopy("COLZ");

    TCanvas *csignal = new TCanvas();
    auto df2 = df1.Define("signalpercent",signalpot,{"nemu","npassive","name"});
    auto hsignalpot = df2.Histo2D({"hsignalpot","Signal fraction of pot;charmrun",10,0,100,100,0,1},"name","signalpercent");
    hsignalpot->DrawCopy("COLZ");
}

void add_info(){ //additional information about number of emulsion and lead plates

    ROOT::RDataFrame df("spill","../input/charm_spills.root");
    
    //methods for getting the number of plates
    auto nemulsions = [](int name){
        int charmtarget = (int) name/10;
        if (charmtarget == 0) return 4;
        if (charmtarget < 3) return 29;
        else return 57;
    };

    auto npassive = [](int name){ //mm of passive slabs
        int charmtarget = (int) name/10;
        switch(charmtarget){
         case 1: return 28;
         case 2: return 56;
         case 3: return 56*2;
         case 4: return 56*3;
         case 5: return 56*4;
         case 6: return 56*5;
         case 0: return 56*4; //4 preshower blocks, no actual brick
        }
    };    
    //Defining the two new columns
    auto dnew = df.Define("nemu",nemulsions,{"name"});
    auto dfinal = dnew.Define("npassive",npassive,{"name"});

    //Drawing the histograms
    TCanvas *c1 = new TCanvas();
    auto hemu = dfinal.Histo1D("nemu");
    hemu->DrawCopy();
    TCanvas *c2 = new TCanvas();
    auto hpassive = dfinal.Histo1D("npassive");
    hpassive->DrawCopy();
    //saving the results
    dfinal.Snapshot("spill","charm_spills.root");
}