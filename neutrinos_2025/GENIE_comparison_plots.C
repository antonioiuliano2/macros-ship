using namespace ROOT;

int fill_process(bool isCC, bool isRES, bool isQEL, bool isCOH, bool isDIS, bool isIMD, bool isNUEL, bool isMEC){
 //check process, then check CC and NC
    if (isCOH){ 
     if (isCC) return 1;
     else return 2;
    }
    else if (isQEL){ 
     if (isCC) return 3;
     else return 4;
    }
    else if (isRES){ 
     if (isCC) return 5;
     else return 6;
    }
    else if (isDIS){ 
     if (isCC) return 7;
     else return 8;
    }
    else if (isIMD){ 
     if (isCC) return 9;
     else return 10;
    }
    else if (isNUEL){ 
     if (isCC) return 11;
     else return 12;
    }
    else if (isMEC){ 
     if (isCC) return 13;
     else return 14;
    }
    else return 0; //no process matching event
}

void GENIE_comparison_plots(){
    EnableImplicitMT();
    //opening file, getting dataframe
    auto df = RDataFrame("gst","/home/antonio/Simulations/test_simulations/test_gevgen_pythia8decayer/G18_02a_00_000/nu_mu_gevgen_iron_5M.gst.root");

    //getting histograms

    auto h_nf = df.Histo1D<int>({"h_nf","Number of final states;n^{f};A.u./1.0 counts",70,0,70},"nf"); //Number of final state particles in hadronic system.
    auto h_El = df.Histo1D({"h_El","Final state primary lepton energy;E_{l}^{f}[GeV];A.u/3.5 GeV",100,0,350},"El"); // Final state primary lepton energy (in GeV).

    
    auto h_nfpip = df.Histo1D<int>({"h_nfpip","Number of final states pi+",40,0,40},"nfpip"); // Number of final state π +
    auto h_nfpim = df.Histo1D<int>({"h_nfpim","Number of final states pi-",40,0,40},"nfpim"); // Number of final state π -
    auto h_nfpi = df.Define("nfpi","nfpip+nfpim+nfpi0")
        .Histo1D<int>({"h_nfpi","Number of final states pi;n_{#pi}^{f};A.u./1.0 counts",40,0,40},"nfpi"); // Number of final state π
    //example one hot encoding
    auto h_process = df.Define("process",fill_process,{"cc","res","qel","coh","dis","imd","nuel","mec"})
        .Histo1D<int>({"h_process","Interaction process;process;% of all events",20,0,20},"process"); // Interaction process
    
    auto h_nfpi0 = df.Histo1D<int>({"h_nfpi0","Number of final states pi0",40,0,40},"nfpi0"); // Number of final state π 0
    
    h_nfpip->Add(h_nfpim.GetPtr());

    //drawing histograms
    TCanvas *c_nf = new TCanvas("c_nf","Number final state particles in hadronic system");
    c_nf->SetLogy();
    h_nf->Scale(1./h_nf->Integral());//normalizing to unity
    h_nf->DrawClone("hist");
    TCanvas *c_nfpip = new TCanvas("c_nfpip","Number final state pions");
    c_nfpip->SetLogy();
    h_nfpi->Scale(1./h_nfpi->Integral());//normalizing to unity
    h_nfpi->DrawClone("hist");
    TCanvas *c_El = new TCanvas("c_El","Final  state primary lepton energy");
    c_El->SetLogy();
    h_El->Scale(1./h_El->Integral());//normalizing to 3.5 GeV width (350/100 bin)
    h_El->DrawClone("hist");
    TCanvas *c_process = new TCanvas("c_process","Interaction process");
    //alphanumeric histogram
    h_process->GetXaxis()->SetBinLabel(2,"#splitline{COH}{CC}");
    h_process->GetXaxis()->SetBinLabel(3,"#splitline{COH}{NC}");
    
    h_process->GetXaxis()->SetBinLabel(4,"#splitline{QEL}{CC}");
    h_process->GetXaxis()->SetBinLabel(5,"#splitline{QEL}{NC}");
    
    h_process->GetXaxis()->SetBinLabel(6,"#splitline{RES}{CC}");
    h_process->GetXaxis()->SetBinLabel(7,"#splitline{RES}{NC}");
    
    h_process->GetXaxis()->SetBinLabel(8,"#splitline{DIS}{CC}");
    h_process->GetXaxis()->SetBinLabel(9,"#splitline{DIS}{NC}");
    
    h_process->GetXaxis()->SetBinLabel(10,"#splitline{IMD}{CC}");
    h_process->GetXaxis()->SetBinLabel(11,"#splitline{IMD}{NC}");
    
    h_process->GetXaxis()->SetBinLabel(12,"#splitline{NUEL}{CC}");
    h_process->GetXaxis()->SetBinLabel(13,"#splitline{NUEL}{NC}");
    
    h_process->GetXaxis()->SetBinLabel(14,"#splitline{MEC}{CC}");
    h_process->GetXaxis()->SetBinLabel(15,"#splitline{MEC}{NC}");
    
    gStyle->SetPaintTextFormat("4.2f");
    h_process->Scale(1./h_process->Integral()*100);//normalizing to percentage
    h_process->SetFillColor(kBlue);
    h_process->DrawClone("hist && text0");


}