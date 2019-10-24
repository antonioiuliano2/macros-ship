void read_spill(TString runname="16"){
    //sum interacting protons and predicted charm for the spills in a selection of runs
    ROOT::RDataFrame df("spill","charm_spills_predictions.root");
    //select the run (or runs, i.e. name<20&&name>0 for first configuration)
    auto df1 = df.Filter((TString("name==")+runname).Data());

    //doing sums (not getting values yet "lazy" operations)
    auto sumpot = df1.Sum("pot");
    auto suminteractingpot = df1.Sum("signalpot");
    auto sumcharm = df1.Sum("primarycharm");

    //actually getting sum results;
    int npots = sumpot.GetValue();
    float nintpots = suminteractingpot.GetValue();
    float ncharms = sumcharm.GetValue();
    //printing the result;
    cout<<TString::Format("On a total of %d pot, estimated %.2e interactions with %.2e charm pairs", npots, nintpots, ncharms)<<endl;
    
}

void compute_spill(){
    //Actual script to compute how many protons have interacted and how many primary charm pairs have been produced

    //**************USED FUNCTIONS, [] is the C++ lambda function*************************
    //fraction of protons not interacting
    auto survivingpot = [] (int nemu, int npassive,  int name){
        //units in mm
        //thickness of each piece
        const float filmthickness = 0.300;
        float slabthickness = 1.0;        
        //interaction length of media
        const float emulambda = 600; //stakeholder, need to find it 
        float passivelambda = 175.9;
        //correction for CH1-R6 special case
        if (name == 16){ 
            passivelambda = 99.5; //CH1-R6 used W instead of PB
            slabthickness = 0.9; //less thick slabs
        }
        float totalemuthickness = nemu * filmthickness;
        float totalslabthickness = npassive * slabthickness;

        float nlambda = totalemuthickness/emulambda + totalslabthickness/passivelambda; //computing the total number of interaction lengths
        float surviving = TMath::Exp(-nlambda); //surviving protons (averaging in the N spills)
        return surviving; //percentage of surviving pot
    };
    //fraction of protons interacting in ECC
    auto signalpot = [] (int nemu, int npassive, int name){
        //units in mm
        //thickness of each piece
        const float filmthickness = 0.300;
        float slabthickness = 1.0;   
        //interaction length of media
        const float emulambda = 600; //from FairShip
        float passivelambda = 175.9;
        //correction for CH1-R6 special case
        if (name == 16){ 
            passivelambda = 99.5; //CH1-R6 used W instead of PB
            slabthickness = 0.9; //less thick slabs
        }
        float totalemuthickness = nemu * filmthickness;
        float totalslabthickness = npassive * slabthickness;

        float nlambda = totalemuthickness/emulambda + totalslabthickness/passivelambda; //computing the total number of interaction lengths
        float preshowerthickness = 0.;

        if (name/10 == 0) return preshowerthickness;       
        if (name/10 == 2) preshowerthickness = 28;
        else if (name/10 > 2) preshowerthickness = 56 *(name/10-2);

        float nlambdapassive = preshowerthickness/passivelambda;
        float inpreshower = (1-TMath::Exp(-nlambdapassive));

        float surviving = TMath::Exp(-nlambda); //surviving protons (averaging in the N spills)
        float inECC = 1 - surviving - inpreshower;
        return inECC; //percentage of pot interacting within the brick
    };
    //Charm produced in primary protoninteractions
    auto charmperpotint = [] (int name){
        float mmtocm = 0.1;
        float cmsquaretobarn = 1e+24;
        //units in CM
        float passivelambda = 175.9 * mmtocm;
        float density = 11.35;
        //values for W, CH1-R6 special case
        if (name==16){
         passivelambda = 99.5 * mmtocm;
         density = 	19.30;
        }
        const float avogadro =  6.022e+23;
        const float charmxsec = 18.1e-6; //in barn, associated charm production cross section per nucleon

        float potxsec = 1./(passivelambda*density*avogadro) * cmsquaretobarn;
        return (charmxsec/potxsec);
        

    }; 

   //****************START OF SCRIPT**************************
    //building a RDataFrame from the inputfile
    ROOT::RDataFrame df("spill","charm_spills.root");
    //drawing number of emulsion films and passive slabs
    TCanvas *c1 = new TCanvas();
    auto hemu = df.Histo1D("nemu");
    hemu->DrawCopy();
    TCanvas *c2 = new TCanvas();
    auto hpassive = df.Histo1D("npassive");
    hpassive->DrawCopy();

    auto entries = df.Filter("name/10==1").Count();
    cout<<"CHARM1 number of runs: "<<*entries<<endl;
    // getting number of surviving protons
    auto df1 = df.Define("survpercent",survivingpot,{"nemu","npassive","name"});
    auto gpot = df1.Histo2D({"h2","Surviving fraction of pot;charmrun",100,0,100,100,0,1},"name","survpercent");

    // getting number of interacting protons
    auto df2 = df1.Define("signalpercent",signalpot,{"nemu","npassive","name"});
    auto hsignalpot = df2.Histo2D({"hsignalpot","Signal fraction of pot;charmrun",100,0,100,100,0,1},"name","signalpercent");
    
    auto dfsave0 = df2.Define("survivingpot","survpercent*pot");
    auto dfsave1 = dfsave0.Define("signalpot","signalpercent*pot");
    //computing charm per pot
    auto dfsave2 = dfsave1.Define("charmfraction",charmperpotint,{"name"});
    auto dfsavefinal = dfsave2.Define("primarycharm","signalpot*charmfraction");
    //saving the result
    dfsavefinal.Snapshot("spill","charm_spills_predictions.root");

    //drawing the plots 
    TCanvas *csurviving = new TCanvas();
    gpot->GetXaxis()->SetBinLabel(10,"CH1");
    gpot->GetXaxis()->SetBinLabel(20,"CH2");
    gpot->GetXaxis()->SetBinLabel(30,"CH3");
    gpot->GetXaxis()->SetBinLabel(40,"CH4");
    gpot->GetXaxis()->SetBinLabel(50,"CH5");
    gpot->GetXaxis()->SetBinLabel(60,"CH6");
    gpot->DrawCopy("COLZ");

    TCanvas *csignal = new TCanvas();
    hsignalpot->GetXaxis()->SetBinLabel(10,"CH1");
    hsignalpot->GetXaxis()->SetBinLabel(20,"CH2");
    hsignalpot->GetXaxis()->SetBinLabel(30,"CH3");
    hsignalpot->GetXaxis()->SetBinLabel(40,"CH4");
    hsignalpot->GetXaxis()->SetBinLabel(50,"CH5");
    hsignalpot->GetXaxis()->SetBinLabel(60,"CH6");
    hsignalpot->DrawCopy("COLZ");
}

void add_info(){ //additional information about number of emulsion and lead plates

    ROOT::RDataFrame df("spill","charm_spills.root");
    
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
         default: return 0;
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
    dfinal.Snapshot("spill","charm_spills.root",{"minute","runcode","pot","spillcode","name","nemu","npassive"});
}
