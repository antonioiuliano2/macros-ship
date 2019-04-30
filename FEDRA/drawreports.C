//Saving linking and alignment reports (first implementation in a general way on 19 December 2018)
//in the b000001 folder, there should be already a folder called 
//plots with subfolders: thicknesses, link_reports, al_reports
//TString: class which allows path concatenation. Access the char* object with Data()
TString run = "CH1-R6";

TString path = "/ship/CHARM2018/" + run +"/b000001/"; 

const int firstplate = 1;
const int lastplate = 29;

void linkreports(){
 TFile *inputfile;
 TCanvas *c;
 for (int i = firstplate; i <= lastplate; i++){
  //opening the file with the reports
  if (i < 10) inputfile = TFile::Open(Form((path+"p00%i/1.%i.0.0.cp.root").Data(),i,i));
  else inputfile = TFile::Open(Form((path+"p0%i/1.%i.0.0.cp.root").Data(),i,i));
  c = (TCanvas*) inputfile->Get("report");

  //saving a full pdf report (very heavy)
  if (i==firstplate) c->Print((path+"plots/link_reports/linking.pdf(").Data(),"pdf");
  else if (i==lastplate) c->Print((path+"plots/link_reports/linking.pdf)").Data(),"pdf");
  else c->Print((path+"plots/link_reports/linking.pdf").Data(),"pdf");

  //saving many png images for quick view
  //c->Print(Form((path+"plots/link_reports/link_p%i.png").Data(),i),"png");  

  //close the file
  inputfile->Close();
 }
}

void alreports(){
 TFile *inputfile;
 TCanvas *c;
 for (int i = firstplate; i <= lastplate; i++){
  if (i == firstplate) continue;
  //opening the file with the report
  inputfile = TFile::Open(Form((path+"AFF/1.%i.0.0_1.%i.0.0.al.root").Data(),i,i-1));
  c = (TCanvas*) inputfile->Get("report_al");
  c->Draw();
  //saving a full pdf report
  /*if (i==firstplate+1) c->Print((path+"plots/al_reports/alignment.pdf(").Data(),"pdf"); 
  else if (i==lastplate) c->Print((path+"plots/al_reports/alignment.pdf)").Data(),"pdf");
  else c->Print((path+"plots/al_reports/alignment.pdf").Data(),"pdf");
*/
  //saving many png images for quick view
  c->Print(Form((path+"plots/al_reports/alignment_p%i_p%i.png").Data(),i,i-1),"png");
  inputfile->Close();
 }
}

//draw all thickness
#include "thickness.C"
void thickness();

void drawallthickness(){

double meanthickness_base = 0.;
double meanthickness_emu = 0.;


for (int i = firstplate; i <= lastplate; i++){
 if (i < 10) TFile *f = TFile::Open(Form((path+"p00%i/1.%i.0.0.raw.root").Data(),i,i));
 else TFile *f = TFile::Open(Form((path+"p0%i/1.%i.0.0.raw.root").Data(),i,i));
 thickness(2);
 gROOT->GetSelectedPad()->GetCanvas()->SetName(Form("canvas%d",i));
 TCanvas *thickcanvas = (TCanvas*) gROOT->GetSelectedPad()->GetCanvas()->GetPrimitive("diff_3");
 //getting the three histograms
 TH1F *hup = thickcanvas->GetPrimitive("up");
 TH1F *hdown = thickcanvas->GetPrimitive("down");
 TH1F *hbase = thickcanvas->GetPrimitive("base");

 hup->GetXaxis()->SetRange(10,300);
 hbase->GetXaxis()->SetRange(10,300);
 hdown->GetXaxis()->SetRange(10,300);

 cout<<"Plate number: "<<i<<endl;
 cout<<"Thickness top layer: "<<hup->GetMean()<<" with RMS "<<hup->GetRMS()<<endl;
 cout<<"Thickness plastic base: "<<hbase->GetMean()<<" with RMS "<<hbase->GetRMS()<<endl;
 cout<<"Thickness bottom layer: "<<hdown->GetMean()<<" with RMS "<<hdown->GetRMS()<<endl;
 
 meanthickness_emu += hup->GetMean();
 meanthickness_base += hbase->GetMean();
 meanthickness_emu += hdown->GetMean();

 if (i == firstplate) gROOT->GetSelectedPad()->GetCanvas()->Print((path+"plots/thicknesses/allthicknesses.pdf(").Data(),"pdf");
 else if(i == lastplate) gROOT->GetSelectedPad()->GetCanvas()->Print((path+"plots/thicknesses/allthicknesses.pdf)").Data(),"pdf");
 else gROOT->GetSelectedPad()->GetCanvas()->Print((path+"plots/thicknesses/allthicknesses.pdf").Data(),"pdf");
 }
 meanthickness_emu = meanthickness_emu/((lastplate - firstplate + 1)*2);
 meanthickness_base = meanthickness_base/(lastplate -firstplate +1);
 cout<<"Average thickness of emulsion: "<<meanthickness_emu<<endl;
 cout<<"Average thickness of plastic base: "<<meanthickness_base<<endl;
}

