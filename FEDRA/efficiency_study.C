//Script for studying detection efficiency in nuclear emulsion plates
//It reads track information from a linked_tracks.root from the same folder where it is executed
#include "TEfficiency.h"
EdbPVRec     *gAli=0;

//TCut trcut = "t.eFlag>=0 &&t.eProb>0.01&&npl >= 25 && s[0].Plate()<=3 && s[nseg-1].Plate()>=27 && TMath::Abs(s.eX-50000) < 15000 && TMath::Abs(s.eY-56000)<15000"; //identifying protons in the central region (all segments must be in this region)
//TCut trcut = "t.eFlag>=0 &&t.eProb>0.01&&npl >= 25 && TMath::Abs(s.eY-46000)<4000"; //identifying protons in the central region (all segments must be in this region)
 //identifying protons in the central region (all segments must be in this region)
//

void check(){ //quick efficiency check, as in check_tr
 TCut trcut = "t.eFlag>=0 &&t.eProb>0.01&&npl >= 25 && s[0].Plate()<=3 && s[nseg-1].Plate()>=27";
 TFile * inputfile = TFile::Open("linked_tracks.root");
 TTree * tracks = (TTree*) gFile->Get("tracks");
 //Doing selections from track tree
 TCanvas *c1 = new TCanvas();
 cout<<".........."<<tracks->GetEntries(trcut)<<endl;

 tracks->Draw("1*(nseg-2)/(npl-2):sqrt(t.eTY*t.eTY+t.eTX*t.eTX)>>heff(25)",trcut,"prof"); //dependance on theta
 
 TCanvas *c2 = new TCanvas();
 tracks->Draw("1*(nseg-2)/(npl-2):t.eY*1e-3:t.eX*1e-3>>heffmap(40,30,70,40,35,75)",trcut,"prof&&COLZ"); //efficiency map

 TCanvas *c3 = new TCanvas();
 tracks->Draw("t.eY*1e-3:t.eX*1e-3>>hxy(40,30,70,40,35,75)",trcut,"COLZ"); //simple 2D map of track positions
 
 //Getting drawn histograms
 TProfile *heff = (TProfile*) gDirectory->Get("heff");
 TProfile2D *heffmap = (TProfile2D*) gDirectory->Get("heffmap"); 
 TH2D *hxy = (TH2D*) gDirectory->Get("hxy");

 //Setting title and axis labels

 heff->SetTitle("Efficiency with tracks from npl>20");
 heff->GetXaxis()->SetTitle("tan(theta)");

 heffmap->SetTitle("Map of efficiency");
 heffmap->GetXaxis()->SetTitle("x[mm]");
 heffmap->GetYaxis()->SetTitle("y[mm]");

 hxy->SetTitle("2D track positions");
 hxy->GetXaxis()->SetTitle("x[mm]");
 hxy->GetYaxis()->SetTitle("y[mm]");
 
}

void efficiency_study(){ //efficiency estimation, as used in OPERA paper

 //TCut trcut = "t.eFlag>=0  &&t.eProb>0.01";
 TCut trcut = "t.eFlag>=0 &&t.eProb>0.01&&npl >= 25";
 const int nfilms = 29;

 TEfficiency* heff = NULL;

 TH1D *hexpected = new TH1D("hexpected", "Tracks expected to be found in each plate", 58,0,58);
 TH1D *hfound = new TH1D("hfound", "Tracks with an associated segment in each plate", 58,0,58); 

 TH2D *hxy = new TH2D("hxy", "2D position of tracks",120,0,120,100,0,100);

 //reading tracks from linked_tracks.root file
 EdbDataProc *dproc = new EdbDataProc();
 dproc->InitVolume(100, trcut); //100 is the code for track reading (as from fedra/src/libEIO/EdbDataSet.cxx)
 gAli = dproc->PVR();
 const int ntracks = gAli->eTracks->GetEntries();
 printf("ntracks after cut = %d\n", ntracks);

 //begin track loop
 EdbTrackP *trk = NULL;
 EdbSegP *seg = NULL;
 int nplate;
 for (int itrk = 0; itrk <ntracks; itrk++){
  trk = (EdbTrackP*)(gAli->eTracks->At(itrk)); //accessing track object

  for (int i = 0; i < nfilms; i++) {
    //if(i!=15)
    hexpected->Fill(i+1);
  } 
  hexpected->Fill(0);//for having y range set from 0 to 1 I add a bin with null efficiency for Plate 0
  
  for (int iseg = 0; iseg < trk->N();iseg++){ //loop on associated segments
    seg = (EdbSegP*) trk->GetSegment(iseg);
    nplate = seg->Plate();
    //if(nplate!=15)
    hfound->Fill(nplate);
  }
  
  hxy->Fill(trk->X()/1000., trk->Y()/1000.);
 } //end of track loop

 //Getting efficiency for all plates
 TCanvas *c1 = new TCanvas();
 if (TEfficiency::CheckConsistency(*hfound,*hexpected)){
  heff = new TEfficiency(*hfound, *hexpected);
  heff->Draw();
  heff->SetTitle(";npl;#epsilon");
  //estimating global efficiency and its error
  double effplate[nfilms]; //efficiency in each plate
  double mean = 0;
  double sd = 0; //standard deviation
  double sd_mean = 0; //standard deviation of the mean (standard error)
  
  for (int ibin=0; ibin < nfilms; ibin++){ 
   effplate[ibin] = heff->GetEfficiency(ibin+2);//the first bin is the dummy one
   mean+=effplate[ibin];
  }
  mean = mean/nfilms;
  for (int ibin=0; ibin < nfilms; ibin++){
   sd+=pow((effplate[ibin]-mean),2);
  }
  sd=TMath::Sqrt(sd/(nfilms-1));
  sd_mean = sd/TMath::Sqrt(nfilms);
  cout<<"Mean efficiency: "<<mean<<"with error: "<<sd_mean<<endl;
 }

 TCanvas *c2 = new TCanvas();
 hxy->Draw("COLZ");
}

void efficiency_study_alltracks(){ //efficiency estimation, dividing for number of tracks
TCut trcut = "t.eFlag>=0 &&t.eProb>0.01&&npl >= 5"; //identifying protons in the central region (all segments must be in this region)

 const int nfilms = 26;
 TEfficiency* heff = NULL;

 TH1D *hexpected = new TH1D("hexpected", "Tracks expected to be found in each plate", 26,2,28);
 TH1D *hfound = new TH1D("hfound", "Tracks with an associated segment in each plate", 26,2,28); 

 TH2D *hxy = new TH2D("hxy", "2D position of tracks", 40,30,70,40,35,75);

 //reading tracks from linked_tracks.root file
 EdbDataProc *dproc = new EdbDataProc();
 dproc->InitVolume(100, trcut); //100 is the code for track reading (as from fedra/src/libEIO/EdbDataSet.cxx)
 gAli = dproc->PVR();
 const int ntracks = gAli->eTracks->GetEntries();
 printf("ntracks after cut = %d\n", ntracks);

 //begin track loop
 EdbTrackP *trk = NULL;
 EdbSegP *seg = NULL;
 EdbSegP *firstseg = NULL;
 EdbSegP *lastseg = NULL;
 int nplate;
 for (int itrk = 0; itrk <ntracks; itrk++){
  trk = (EdbTrackP*)(gAli->eTracks->At(itrk)); //accessing track object

  firstseg = trk->GetSegmentFirst();
  lastseg = trk->GetSegmentLast();
  for (int i = firstseg->Plate(); i <= lastseg->Plate();i++){
   hexpected->Fill(i);
  }
  /*for (int i = 0; i < nfilms; i++) {
  hexpected->Fill(i+1);
  }*/ 
  hexpected->Fill(0);//for having y range set from 0 to 1 I add a bin with null efficiency for Plate 0

  for (int iseg = 0; iseg < trk->N();iseg++){ //loop on associated segments
   seg = (EdbSegP*) trk->GetSegment(iseg);
   nplate = seg->Plate();
   hfound->Fill(nplate);
  }

  hxy->Fill(trk->X()/1000., trk->Y()/1000.);
 } //end of track loop

 //Getting efficiency for all plates

 TCanvas *c1 = new TCanvas();
 if (TEfficiency::CheckConsistency(*hfound,*hexpected)){
  heff = new TEfficiency(*hfound, *hexpected);
  heff->Draw();
  heff->SetTitle("Efficiency for each plate;npl");
  //estimating global efficiency and its error
  double effplate[nfilms]; //efficiency in each plate
  double mean = 0;
  double sd = 0; //standard deviation
  double sd_mean = 0; //standard deviation of the mean (standard error)
  /*
  for (int ibin=0; ibin < heff->GetNbinsX(); ibin++){ 
   effplate[ibin] = heff->GetEfficiency(ibin+1);//the first bin is the dummy one
   mean+=effplate[ibin];
  }
  mean = mean/nfilms;
  for (int ibin=0; ibin < heff->GetNbinsX(); ibin++){
   sd+=(effplate[ibin]-mean)**2;
  }
  sd=TMath::Sqrt(sd/(nfilms-1));
  sd_mean = sd/TMath::Sqrt(nfilms);
  cout<<"Mean efficiency: "<<mean<<"with error: "<<sd_mean<<endl;*/
 }

 TCanvas *c2 = new TCanvas();
 hxy->Draw("COLZ");
}

void pos_distribution(){
 TFile * inputfile = TFile::Open("linked_tracks.root");
 TTree * tracks = (TTree*) inputfile->Get("tracks");

 TCanvas *c1 = new TCanvas();

 tracks->Draw("s[nseg-1].eY*1e-3:s[nseg-1].eX*1e-3>>hxy(125,0,125,100,0,100)","","COLZ");
 TH2D *hxy = (TH2D*) gDirectory->FindObject("hxy");

 hxy->SetTitle("xy distribution end point of tracks;x[mm];y[mm]");
 hxy->Draw("COLZ");

 gStyle->SetStatX(0.5);
 gStyle->SetStatY(0.9); 

 TCanvas *c2 = new TCanvas();

 tracks->Draw("s.ePID>>hPID");

 TH1I *hPID = (TH1I*) gDirectory->FindObject("hPID");
 hPID->SetTitle("Plate ID for each segment belonging to track;PID");
 hPID->Draw();

}
