/*
Checking alignment performances. Modified for DESY2019 bricks format in 16/12/2019

*/
void align_check(TString runname, int lastplate=57, int firstplate=1);
//checks dz between plates after alignment
void checkdz(TString runname, int lastplate=57, int firstplate=1){
 //syntax for brickname is 002, 022, 111
 TString dir = TString("/ship/DESY2019/");
 
 TString bricknumber = TString::Format("%d",TString(runname(3,3)).Atoi());
 cout<<bricknumber<<endl;

 TString setposition = TString("/b00000"+bricknumber+"/b00000"+bricknumber+".0.0.0.set.root"); //example: /b000005/b000005.0.0.0.set.root

 TH1F *hdzA = new TH1F("hdzA",(TString("Dz obtained after alignment for brick ")+runname).Data(),300,0,3000);

 TGraph *hgraphA = new TGraph();
 hgraphA->SetName("hgraph");

 TFile *fileA = TFile::Open((dir+runname+setposition).Data());
 EdbScanSet *setA = (EdbScanSet*) fileA->Get("set");

 float dzA;
 //loop on plates
 for (int i = firstplate; i< lastplate; i++){
  dzA = setA->GetDZP2P(i,i+1);
  cout<<"from plate "<<i<<" we have dz "<<dzA<<endl;
  //filling histograms and graphs
  hdzA->Fill(dzA);
  
  hgraphA->SetPoint(i,(i+i+1)/2.,dzA);
 }
 //comparing the histograms
 TCanvas *cdz = new TCanvas("cdz");
 cdz->Divide(1,2);
 cdz->cd(1);
 hdzA->Draw();
 hdzA->GetXaxis()->SetTitle("dZ[#mum]");
 cdz->cd(2);
 hgraphA->SetTitle("Dz obtained after alignment for brick")+runname;
 hgraphA->Draw("AP*");
 hgraphA->GetXaxis()->SetTitle("nplate");
 hgraphA->GetYaxis()->SetTitle("dZ[#mum]");

}

//check alignment plots (currently xy and zphi)

void align_check(TString runname,int lastplate, int firstplate){
 TString dir = TString("/ship/DESY2019/");
 TString bricknumber = TString::Format("%d",TString(runname(3,3)).Atoi());
 TString path = "/ship/DESY2019/" + runname +"/b00000"+bricknumber+"/"; 
 TFile *inputfile;
 //when bricks are made of many plates, needed more than 1 canvas
 const int nplates = lastplate - firstplate;
 const int canvasdimx = 4;
 const int canvasdimy = 4;
 const int canvasarea = canvasdimx * canvasdimy;
 //how many canvases do I need?
 const int ncanvases=nplates/canvasarea+1; 
 TCanvas *c[ncanvases];
 TCanvas *cxy[ncanvases];
 TCanvas *cTXTY[ncanvases];
 for (int icanvas = 0; icanvas < ncanvases; icanvas++){
 c[icanvas] = new TCanvas((runname+TString::Format("_dz_phi_coarse_%d",icanvas)).Data());
 c[icanvas]->SetTitle((runname+TString::Format("_dz_phi_coarse_%d",icanvas)).Data());
 c[icanvas]->Divide(canvasdimx,canvasdimy); //NEED TO OPEN OTHER CANVASES WHEN PLATES INCREASE

 cxy[icanvas] = new TCanvas((runname+TString::Format("_dxy_final_%d",icanvas)).Data());
 cxy[icanvas]->SetTitle((runname+TString::Format("_dxy_final_%d",icanvas)).Data());
 cxy[icanvas]->Divide(canvasdimx,canvasdimy);

 cTXTY[icanvas] = new TCanvas((runname+TString::Format("_dtxty_final_%d",icanvas)).Data());
 cTXTY[icanvas]->SetTitle((runname+TString::Format("_dtxty_final_%d",icanvas)).Data());
 cTXTY[icanvas]->Divide(canvasdimx,canvasdimy);
 }

 TH2F *hzphi[nplates];
 TH2F *hxy[nplates];
 TH2F *hTXTY[nplates];

 //const int lastplate= 57;
 int iplate = 0;
 int whichcanvas = 0;
 for (int i = firstplate; i <= lastplate; i++){
  if (i == firstplate) continue;
  //computing in which canvas and subcanvas plot the histograms
  whichcanvas = iplate/canvasarea;
  c[whichcanvas]->cd(iplate-whichcanvas*canvasarea+1);
  //opening the file with the report
  inputfile = TFile::Open(Form((path+"AFF/"+bricknumber+".%i.0.0_"+bricknumber+".%i.0.0.al.root").Data(),i,i-1)); //example:AFF/2.10.0.0_2.9.0.0.al.root
  if (inputfile){ //check if file is present
  //zphi plot
   if (inputfile->GetListOfKeys()->Contains("zphi_coarse") ){ //check if histogram is present (i.e. alignment was not interrupted leading to a zombie-like file, even if not seen as zombie)
   hzphi[iplate] = (TH2F*)inputfile->Get("zphi_coarse");
   hzphi[iplate]->SetTitle(Form("plate %i to plate %i;dphi[rad];dz[#mum]",i,i-1));
   hzphi[iplate]->Draw("COLZ");   
  }
  //xy residuals
  cxy[whichcanvas]->cd(iplate-whichcanvas*canvasarea+1);
  if (inputfile->GetListOfKeys()->Contains("xy_final") ){
  hxy[iplate] = (TH2F*)inputfile->Get("xy_final");
  hxy[iplate]->SetTitle(Form("plate %i to plate %i;dx[#mum];dy[#mum]",i,i-1));
  hxy[iplate]->Draw("COLZ");
   }
  //txty residuals, not saved in AFF files, but I can get it directly from report canvas
  cTXTY[whichcanvas]->cd(iplate-whichcanvas*canvasarea+1);
  if (inputfile->GetListOfKeys()->Contains("report_al") ){
  TCanvas *creport = (TCanvas*) inputfile->Get("report_al");
  //from canvas to subcanvas, finally to histogram
  TCanvas *subcanvas = (TCanvas*)((TCanvas*)creport->GetPrimitive("c"))->GetPrimitive("c_6");
  hTXTY[iplate] = (TH2F*)subcanvas->GetPrimitive("dTXdTY");

  hTXTY[iplate]->SetTitle(Form("plate %i to plate %i;dx[#mum];dy[#mum]",i,i-1));
  cTXTY[whichcanvas]->cd(iplate-whichcanvas*canvasarea+1); //we need to go back to our canvas, otherwise we will draw in that canvas
  hTXTY[iplate]->Draw("COLZ");
   }
  }
  iplate++;
 }
}

int global_check_align(TString runname = "",int lastplate=57, int firstplate=1){

 if (runname==" "){ 
  cout<<"Script for checking alignment. Usage: check_align(\"RUN5\",lastplate,firstplate). The check for dz needs a scanset prepared with all the plates"<<endl; 
  return 1;
 }
 align_check(runname,lastplate,firstplate);
 checkdz(runname,lastplate,firstplate);
 return 0;
}
