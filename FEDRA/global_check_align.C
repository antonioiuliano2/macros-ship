//comparing the dz after alignment for a brick with tungsten and one with lead, both use Slavich emulsions
void align_check(TString runname, int lastplate=29, int firstplate=1);
void align_check(TString runname, int lastplate, int firstplate){
 //TString run = "CH1-R6";

 TString path = "/ship/CHARM2018/" + runname +"/b000001/"; 
 TFile *inputfile;
 TCanvas *c = new TCanvas((runname+TString("_dz_phi_coarse")).Data());
 c->SetTitle((runname+TString("_dz_phi_coarse")).Data());
 if (lastplate>29)c->Divide(7,8);
 else c->Divide(5,6);

 TCanvas *cxy = new TCanvas((runname+TString("_dxy_final")).Data());
 cxy->SetTitle((runname+TString("_dxy_final")).Data());
 if (lastplate>29)cxy->Divide(7,8);
 else cxy->Divide(5,6);

 TH2F *hzphi[57];
 TH2F *hxy[57];

 TGraph *peakgraph = new TGraph();
 peakgraph->SetName("peakgraph");
 TGraph *meangraph = new TGraph();
 meangraph->SetName("meangraph");
 TGraph *noisegraph = new TGraph();
 noisegraph->SetName("noisegraph");

 EdbPeak2 *combinations;

 TH1F *hdzA = new TH1F("hdzA",(TString("Dz obtained after alignment for brick ")+runname).Data(),20,1200,1400);

 TGraph *hgraphA = new TGraph();
 hgraphA->SetName("hgraph");

 float dzA;
 int xbin, ybin,zbin;

 //const int lastplate= 57;
 for (int i = firstplate; i <= lastplate; i++){
  if (i == firstplate) continue;
  c->cd(i-1);
  //opening the file with the report
  inputfile = TFile::Open(Form((path+"AFF/1.%i.0.0_1.%i.0.0.al.root").Data(),i,i-1));
  if (inputfile){ //check if file is present
   if (inputfile->GetListOfKeys()->Contains("zphi_coarse") ){ //check if histogram is present (i.e. alignment was not interrupted leading to a zombie-like file, even if not seen as zombie)
   //zphi histogram
   hzphi[i-2] = (TH2F*)inputfile->Get("zphi_coarse");
   hzphi[i-2]->SetTitle(Form("plate %i to plate %i",i,i-1));
   hzphi[i-2]->Draw("COLZ");   

   //finding the peak in z
   hzphi[i-2]->GetBinXYZ(hzphi[i-2]->GetMaximumBin(),xbin,ybin,zbin);
   dzA = TMath::Abs(hzphi[i-2]->GetYaxis()->GetBinCenter(ybin));
   hdzA->Fill(dzA);
   hgraphA->SetPoint(i,(i+i+1)/2.,dzA);
  }
  if (inputfile->GetListOfKeys()->Contains("peak2c") ){
  combinations = (EdbPeak2*) inputfile->Get("peak2c");
  cout<<"Combinations between plates: "<<i<<"and "<<i+1<<endl;
  combinations->Print();
  peakgraph->SetPoint(i,(i+i+1)/2.,combinations->Peak(0));
  meangraph->SetPoint(i,(i+i+1)/2.,combinations->Mean3(0));
  noisegraph->SetPoint(i,(i+i+1)/2.,combinations->Mean(0));
}
  //xy histogram
  cxy->cd(i-1);
  if (inputfile->GetListOfKeys()->Contains("xy_final") ){
  hxy[i-2] = (TH2F*)inputfile->Get("xy_final");
  hxy[i-2]->SetTitle(Form("plate %i to plate %i",i,i-1));
  hxy[i-2]->Draw("COLZ");
  }

  }
 }
 
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

  TCanvas *ccomb = new TCanvas("ccomb");
  peakgraph->SetTitle("Study of combinations");
  peakgraph->SetMarkerColor(kRed);
  peakgraph->Draw("AP*");
  meangraph->SetMarkerColor(kYellow);
  meangraph->Draw("P*SAME");
  noisegraph->Draw("P*SAME");
  peakgraph->GetXaxis()->SetTitle("IDPlate");
  //adding the legend
  TLegend *legend = new TLegend(0.1,0.7,0.48,0.9);
  legend->AddEntry("peakgraph","Main peak of aligned couples","lp");
  //legend->AddEntry("meangraph","Mean in a 3 X 3 region","lp");
  legend->AddEntry("noisegraph","Average of random combinations","lp");
  legend->Draw();
}

void global_check_align(){

 const int nalignedbricks = 15;

 TString alignedbricks[nalignedbricks];

 alignedbricks[0] = TString("CH1-R1");
 alignedbricks[1] = TString("CH1-R2");
 alignedbricks[2] = TString("CH1-R4");
 alignedbricks[3] = TString("CH1-R6");

 alignedbricks[4] = TString("CH2-R1");
 alignedbricks[5] = TString("CH2-R3");
 alignedbricks[6] = TString("CH2-R4");
 alignedbricks[7] = TString("CH2-R5");
 alignedbricks[8] = TString("CH2-R6");

 alignedbricks[9] = TString("CH3-R2");
 alignedbricks[10] = TString("CH4-R2");
 alignedbricks[11] = TString("CH5-R2");
 alignedbricks[12] = TString("CH5-R3");
 alignedbricks[13] = TString("CH6-R2");
 alignedbricks[14] = TString("CH6-R3");

 for (int i=0; i < nalignedbricks; i++){

  if (i < 9) check_align(alignedbricks[i],29);
  else check_align(alignedbricks[i],57);
 }

}


void checkdzold(TString runname, int nplates=29, int firstplate=1){
 //TFile *file = new TFile((runname+TString(" quicksave.root")).Data(),"RECREATE");
 TString dir = TString("/ship/CHARM2018/");

 TString brickA = runname;

 TString setposition = TString("/b000001/b000001.0.0.0.set.root");

 TH1F *hdzA = new TH1F("hdzA",(TString("Dz obtained after alignment for brick ")+brickA).Data(),20,1200,1400);

 TGraph *hgraphA = new TGraph();
 hgraphA->SetName("hgraph");


 TFile *fileA = TFile::Open((dir+brickA+setposition).Data());
 EdbScanSet *setA = (EdbScanSet*) fileA->Get("set");

 const int lastplate = nplates;
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
 hgraphA->SetTitle("Dz obtained after alignment for brick")+brickA;
 hgraphA->Draw("AP*");
 hgraphA->GetXaxis()->SetTitle("nplate");
 hgraphA->GetYaxis()->SetTitle("dZ[#mum]");

 //drawing two canvas with all the alignment checks
 align_check(brickA,nplates);
 //file->cd("");
 //hgraphA->Write();
}

