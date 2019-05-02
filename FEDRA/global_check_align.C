//comparing the dz after alignment for a brick with tungsten and one with lead, both use Slavich emulsions
void align_check(TString runname, int lastplate=29);
void comparedz(){

 TString brickA = TString("/ship/CHARM2018/CH1-R6/");
 TString brickB = TString("/ship/CHARM2018/CH1-R4/");

 TString setposition = TString("b000001/b000001.0.0.0.set.root");

 TH1F *hdzA = new TH1F("hdzA","Dz obtained after alignment for brick CH1R6",20,1200,1400);
 TH1F *hdzB = new TH1F("hdzB","Dz obtained after alignment for brick CH1R4",20,1200,1400);

 TGraph *hgraphA = new TGraph();
 TGraph *hgraphB = new TGraph();

 TFile *fileA = TFile::Open((brickA+setposition).Data());
 EdbScanSet *setA = (EdbScanSet*) fileA->Get("set");

 TFile *fileB = TFile::Open((brickB+setposition).Data());
 EdbScanSet *setB = (EdbScanSet*) fileB->Get("set");

 const int firstplate = 1;
 const int lastplate = 29;
 float dzA, dzB;
 //loop on plates
 for (int i = firstplate; i< lastplate; i++){
  dzA = setA->GetDZP2P(i,i+1);
  dzB = setB->GetDZP2P(i,i+1);
  cout<<"from plate "<<i<<" we have dz for tungsten "<<dzA<< " and for lead "<<dzB<<endl;
  //filling histograms and graphs
  hdzA->Fill(dzA);
  hdzB->Fill(dzB);
  
  hgraphA->SetPoint(i,(i+i+1)/2.,dzA);
  hgraphB->SetPoint(i,(i+i+1)/2.,dzB);
 }
 //comparing the histograms
 TCanvas *cdz = new TCanvas("cdz");
 cdz->Divide(1,2);
 cdz->cd(1);
 hdzA->Draw();
 hdzA->GetXaxis()->SetTitle("dZ[#mum]");
 cdz->cd(2);
 hdzB->Draw();
 hdzB->GetXaxis()->SetTitle("dZ[#mum]");

 //comparing the graphs
 TCanvas *cdzgraph = new TCanvas("cdzgraph");
 cdzgraph->Divide(1,2);
 cdzgraph->cd(1);
 hgraphA->SetTitle("Dz obtained after alignment for brick CH1R6");
 hgraphA->Draw("AP*");
 hgraphA->GetYaxis()->SetTitle("dZ[#mum]");
 cdzgraph->cd(2);
 hgraphB->SetTitle("Dz obtained after alignment for brick CH1R4");
 hgraphB->Draw("AP*");
 hgraphB->GetYaxis()->SetTitle("dZ[#mum]");

 //drawing two canvas with all the alignment checks
 align_check("CH1-R4");
 align_check("CH1-R6");
}

void align_check(TString runname, int lastplate){
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

 const int firstplate = 1;
 //const int lastplate= 57;
 for (int i = firstplate; i <= lastplate; i++){
  if (i == firstplate) continue;
  c->cd(i-1);
  //opening the file with the report
  inputfile = TFile::Open(Form((path+"AFF/1.%i.0.0_1.%i.0.0.al.root").Data(),i,i-1));
  if (inputfile){ //check if file is present
   if (inputfile->GetListOfKeys()->Contains("zphi_coarse") ){ //check if histogram is present (i.e. alignment was not interrupted leading to a zombie-like file, even if not seen as zombie)
   hzphi[i-2] = (TH2F*)inputfile->Get("zphi_coarse");
   hzphi[i-2]->SetTitle(Form("plate %i to plate %i",i,i-1));
   hzphi[i-2]->Draw("COLZ");   
  }

  cxy->cd(i-1);
  if (inputfile->GetListOfKeys()->Contains("xy_final") ){
  hxy[i-2] = (TH2F*)inputfile->Get("xy_final");
  hxy[i-2]->SetTitle(Form("plate %i to plate %i",i,i-1));
  hxy[i-2]->Draw("COLZ");
  }

  }
 }
/*  c->Print((path + TString("plots/al_reports/allzphi_")+runname+TString(".png")).Data());
  cxy->Print((path + TString("plots/al_reports/allxy_final_")+runname+TString(".png")).Data());*/
}

void global_check_align(){

 const int nalignedbricks = 13;

 TString alignedbricks[nalignedbricks];

 alignedbricks[0] = TString("CH1-R4");
 alignedbricks[1] = TString("CH1-R6");

 alignedbricks[2] = TString("CH2-R1");
 alignedbricks[3] = TString("CH2-R3");
 alignedbricks[4] = TString("CH2-R4");
 alignedbricks[5] = TString("CH2-R5");
 alignedbricks[6] = TString("CH2-R6");

 alignedbricks[7] = TString("CH3-R2");
 alignedbricks[8] = TString("CH4-R2");
 alignedbricks[9] = TString("CH5-R2");
 alignedbricks[10] = TString("CH5-R3");
 alignedbricks[11] = TString("CH6-R2");
 alignedbricks[12] = TString("CH6-R3");

 for (int i=0; i < nalignedbricks; i++){

  if (i < 7) align_check(alignedbricks[i],29);
  else align_check(alignedbricks[i],57);
 }

}
