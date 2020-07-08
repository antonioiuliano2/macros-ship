//Saving linking and alignment reports (first implementation in a general way on 19 December 2018)
//in the b000001 folder, there should be already a folder called 
//plots with subfolders: thicknesses, link_reports, al_reports
//TString: class which allows path concatenation. Access the char* object with Data()
TString run = "CH5-R2";

TString path = "/ship/CHARM2018/" + run +"/b000001/"; 

const int firstplate = 1;
const int lastplate = 57;

void linkreports(){
 TFile *inputfile;
 TCanvas *c;
 for (int i = firstplate; i <= lastplate; i++){
   //opening the file with the reports
   if (i < 10) inputfile = TFile::Open(Form((path+"p00%i/1.%i.0.0.cp.root").Data(),i,i));
   else inputfile = TFile::Open(Form((path+"p0%i/1.%i.0.0.cp.root").Data(),i,i));
   if (inputfile){ //check if file is present
     if (inputfile->GetListOfKeys()->Contains("report") ){ //check if histogram is present (i.e. alignment was not interrupted leading to a zombie-like file, even if not seen as zombie)
       c = (TCanvas*) inputfile->Get("report");
       
       //saving a full pdf report (very heavy)
       //if (i==firstplate) c->Print((path+"plots/link_reports/linking.pdf(").Data(),"pdf");
       //else if (i==lastplate) c->Print((path+"plots/link_reports/linking.pdf)").Data(),"pdf");
       //else c->Print((path+"plots/link_reports/linking.pdf").Data(),"pdf");
       
       //saving many png images for quick view
       c->Draw();
       c->Print(Form((path+"plots/link_reports/link_p%i.png").Data(),i),"png");  
       
       //close the file
     }//end check for key
     else cout<<"Warning: key missing for plate "<<i<<endl;
     inputfile->Close();
   }//endcheck for file existency
   else cout<<"Warning: cp file missing for plate "<<i<<endl;
 }
}

void alreports(){
  TFile *inputfile;
  TCanvas *c;
  for (int i = firstplate; i <= lastplate; i++){
    if (i == firstplate) continue;
    //opening the file with the report
    inputfile = TFile::Open(Form((path+"AFF/1.%i.0.0_1.%i.0.0.al.root").Data(),i,i-1));
    if (inputfile){ //check if file is present
      if (inputfile->GetListOfKeys()->Contains("report") ){ 
	c = (TCanvas*) inputfile->Get("report_al");
	c->Draw();
	//saving a full pdf report
	/*if (i==firstplate+1) c->Print((path+"plots/al_reports/alignment.pdf(").Data(),"pdf"); 
	  else if (i==lastplate) c->Print((path+"plots/al_reports/alignment.pdf)").Data(),"pdf");
	  else c->Print((path+"plots/al_reports/alignment.pdf").Data(),"pdf");
	*/
	//saving many png images for quick view
	c->Print(Form((path+"plots/al_reports/alignment_p%i_p%i.png").Data(),i,i-1),"png");
      }
      inputfile->Close();
      else cout<<"Warning: report canvas missing for plate "<<i<<endl;
    }
    else cout<<"Warning: file missing for plate "<<i<<endl;
  }
}

//draw all thickness
#include "thickness.C"
void thickness();

void drawallthicknesses(){
  //draw all thicknesses and plot them, with respect to number of plate
  float meanthickness_base = 0.;
  float meanthickness_emu = 0.;
  
  TGraphErrors *graphthicknesstop = new TGraphErrors();
  TGraphErrors *graphthicknessbot = new TGraphErrors();
  TGraphErrors *graphthicknessbase = new TGraphErrors();
  //dispersion instead of RMS
  TGraphErrors *graphthicknesstop_disp = new TGraphErrors();
  TGraphErrors *graphthicknessbot_disp = new TGraphErrors();
  TGraphErrors *graphthicknessbase_disp = new TGraphErrors();
  
  float thicknesstop,thicknessbot,thicknessbase;
  float thicknesserrortop, thicknesserrorbot, thicknesserrorbase;
  
  float upthmin, upthmax;
  float bothmin, bothmax;
  float basethmin, basethmax;
  
  float threshold = 50.;
  
  int ipoint = 0;

for (int i = firstplate; i <= lastplate; i++){
  /* if (i == 40){
 ipoint = ipoint +1;
 continue;
}*/
 if (i < 10) TFile *f = TFile::Open(Form((path+"p00%i/1.%i.0.0.raw.root").Data(),i,i));
 else TFile *f = TFile::Open(Form((path+"p0%i/1.%i.0.0.raw.root").Data(),i,i));
 //executing script from Valeri to draw thickness curves
 thickness();

 TCanvas *diff = gROOT->FindObject("diff");
 diff->SetName(Form("canvas%d",i));
 TCanvas *thickcanvas = (TCanvas*) diff->GetPrimitive("diff_3");
 //getting the three histograms
 TH1F *hup = thickcanvas->GetPrimitive("up");
 TH1F *hdown = thickcanvas->GetPrimitive("down");
 TH1F *hbase = thickcanvas->GetPrimitive("base");
 //removing peak at 0 (i.e. no surface acquired);
 hup->SetBinContent(1,0); 
 hbase->SetBinContent(1,0);
 hdown->SetBinContent(1,0); 
 //getting thickness values and widths
 thicknesstop = hup->GetMean();
 thicknessbot = hdown->GetMean();
 thicknessbase = hbase->GetMean();
 thicknesserrortop = hup->GetRMS();
 thicknesserrorbot = hdown->GetRMS();
 thicknesserrorbase = hbase->GetRMS();
 //use semidispersion, RMS not good due to peaks
 upthmax = hup->GetBinCenter(hup->FindLastBinAbove(threshold));
 upthmin = hup->GetBinCenter(hup->FindFirstBinAbove(threshold));
 basethmax = hbase->GetBinCenter(hbase->FindLastBinAbove(threshold));
 basethmin = hbase->GetBinCenter(hbase->FindFirstBinAbove(threshold));
 botthmax = hdown->GetBinCenter(hdown->FindLastBinAbove(threshold));
 botthmin = hdown->GetBinCenter(hdown->FindFirstBinAbove(threshold));
 
 //printout
 cout<<"Plate number: "<<i<<endl;
 cout<<"Thickness top layer: "<<thicknesstop<<" with RMS "<<thicknesserrortop<<endl;
 cout<<"Thickness plastic base: "<<thicknessbase<<" with RMS "<<thicknesserrorbase<<endl;
 cout<<"Thickness bottom layer: "<<thicknessbot<<" with RMS "<<thicknesserrorbot<<endl;

 cout<<"Used threshold for thickness range "<<threshold<<endl;
 cout<<"Top: min thickness "<<upthmin<<" max thickness "<<upthmax<<endl;
 cout<<"Base: min thickness "<<basethmin<<" max thickness "<<basethmax<<endl;
 cout<<"Bot: min thickness "<<botthmin<<" max thickness "<<botthmax<<endl;
 //plotting thicknesses 
 ipoint++;

 graphthicknesstop->SetPoint(ipoint,i,thicknesstop); 
 graphthicknessbot->SetPoint(ipoint,i,thicknessbot); 
 graphthicknessbase->SetPoint(ipoint,i,thicknessbase);
 //RMS as vertical error bars
 graphthicknesstop->SetPointError(ipoint,0.5,thicknesserrortop);
 graphthicknessbot->SetPointError(ipoint,0.5,thicknesserrorbot);
 graphthicknessbase->SetPointError(ipoint,0.5,thicknesserrorbase);

 graphthicknesstop_disp->SetPoint(ipoint,i,thicknesstop); 
 graphthicknessbot_disp->SetPoint(ipoint,i,thicknessbot); 
 graphthicknessbase_disp->SetPoint(ipoint,i,thicknessbase);
 //Semidispersion as vertical error bars
 graphthicknesstop_disp->SetPointError(ipoint,0.5,(upthmax - upthmin)/2.);
 graphthicknessbot_disp->SetPointError(ipoint,0.5,(botthmax - botthmin)/2.);
 graphthicknessbase_disp->SetPointError(ipoint,0.5,(basethmax - basethmin)/2.);

 meanthickness_emu += thicknesstop;
 meanthickness_base += thicknessbase;
 meanthickness_emu += thicknessbot;

 gROOT->GetSelectedPad()->GetCanvas()->Print((path+Form("plots/thicknesses/%i_thickness.png",i)).Data(),"png");
 /*if (i == firstplate) gROOT->GetSelectedPad()->GetCanvas()->Print((path+"plots/thicknesses/allthicknesses.pdf(").Data(),"pdf");
 else if(i == lastplate) gROOT->GetSelectedPad()->GetCanvas()->Print((path+"plots/thicknesses/allthicknesses.pdf)").Data(),"pdf");
 else gROOT->GetSelectedPad()->GetCanvas()->Print((path+"plots/thicknesses/allthicknesses.pdf").Data(),"pdf");*/
 }
 meanthickness_emu = meanthickness_emu/((lastplate - firstplate + 1)*2);
 meanthickness_base = meanthickness_base/(lastplate -firstplate +1);
 //drawing the values
 TCanvas *cgraphemu = new TCanvas();
 cgraphemu->Divide(1,2);
 cgraphemu->cd(1);
 graphthicknesstop->SetTitle("Thickness of emulsion layer top (with RMS);iplate;thickness[#mum]");
 graphthicknesstop->Draw("AP");
 cgraphemu->cd(2);
 graphthicknessbot->Draw("AP");
 graphthicknessbot->SetTitle("Thickness of emulsion layer bot (with RMS);iplate;thickness[#mum]");
 cgraphemu->Print((path+"plots/thicknesses/thicknessemu_graph.png").Data(),"png");
 cgraphemu->Close();
 TCanvas *cgraphbase = new TCanvas();
 graphthicknessbase->Draw("AP");
 graphthicknessbase->SetTitle("Thickness of emulsion layer base (with RMS);iplate;thickness[#mum]");
 cgraphbase->Print((path+"plots/thicknesses/thicknessbase_graph.png").Data(),"png");
 cgraphbase->Close();

 TCanvas *cgraphemu2 = new TCanvas();
 cgraphemu2->Divide(1,2);
 cgraphemu2->cd(1);
 graphthicknesstop_disp->SetTitle("Thickness of emulsion layer top (with semidispersion);iplate;thickness[#mum]");
 graphthicknesstop_disp->Draw("AP");
 cgraphemu2->cd(2);
 graphthicknessbot_disp->Draw("AP");
 graphthicknessbot_disp->SetTitle("Thickness of emulsion layer bot(with semidispersion);iplate;thickness[#mum]");
 cgraphemu2->Print((path+"plots/thicknesses/thicknessemu_graph_dispersion.png").Data(),"png");
 cgraphemu2->Close();
 TCanvas *cgraphbase2 = new TCanvas();
 graphthicknessbase_disp->Draw("AP");
 graphthicknessbase_disp->SetTitle("Thickness of emulsion layer base(with semidispersion);iplate;thickness[#mum]");
 cgraphbase2->Print((path+"plots/thicknesses/thicknessbase_graph_dispersion.png").Data(),"png");
 cgraphbase2->Close();

 cout<<"Average thickness of emulsion: "<<meanthickness_emu<<endl;
 cout<<"Average thickness of plastic base: "<<meanthickness_base<<endl;
}

