//------------------------------------------------------
// To check the quality of raw data
// usage:
//   a) interactive: 
//         $ root -l pl56.root
//         root[0] .x check_raw.C            draw plots on the screen; use saved into rootfile canvases if any
//         root[0] .x check_raw.C(1)         the pictures will be saved as gif files
//   b) batch mode:
//         $ root -b -q check_raw.C\(2,\"pl56.root\"\)    run file will be opened for update to save canvases inside
//         $ root -b -q check_raw.C\(3,\"pl56.root\"\)    run file will be opened for update 
//                                                        to save canvases inside AND
//                                                        the pictures will be saved ALSO as gif files
//
//------------------------------------------------------
#define XSIZE 1200
#define YSIZE  800

#if ROOT_VERSION_CODE < 6 //headers should not be included in ROOT6, loaded with .pcm files
#ifndef __CINT__
#include "TCanvas.h"
#include "TFile.h"
#include "TSystem.h"
#include "TCut.h"
#include "TStyle.h"
#include "TTree.h"
#include "TH1.h"
#include "TF1.h"
#include "Riostream.h"
#include "EdbRun.h"
#include "EdbView.h"
#include "EdbCluster.h"
#include "EdbSegment.h"
using namespace std;

#endif
#endif

EdbRun* run ;
TTree* cr_tree;

//Declarations of used methods
void init();
void check_z(TCanvas *c);
void check_view(TCanvas *c);
void check_surf(TCanvas *c);
void check_raw( int output=0, char *fname=0 );
//-------------------------------------------------------
void init()
{
  cr_tree->SetAlias("puls",  "segments.GetPuls()");
  cr_tree->SetAlias("volume","segments.GetVolume()");
}

//-------------------------------------------------------
void check_z(TCanvas *c)
{
  c->cd(1);       cr_tree->Draw("eZ1:eXview:eYview");
  c->cd(2);       cr_tree->Draw("eZ1:eAreaID*121+eViewID");
  c->cd(3); {
    cr_tree->Draw("eNcl:eZframe-eZ2");
    cr_tree->SetMarkerColor(kBlue);
    cr_tree->Draw("eNcl:eZframe-eZ2","eNframesTop==0","same");
    cr_tree->SetMarkerColor(kRed);
    cr_tree->Draw("eNcl:eZframe-eZ2","eNframesBot==0","same");
  }
  cr_tree->SetLineColor(1);

  c->cd(4);       cr_tree->Draw("eZ2-eZ3");
  gStyle->SetOptStat("nemr");
  cr_tree->SetLineColor(kBlue);
  c->cd(5);       cr_tree->Draw("eNcl:eFrameID","eNframesTop==0","prof");
  gStyle->SetOptStat("ne");
  c->cd(7);       cr_tree->Draw("eZ0","eNframesTop==0");
  gStyle->SetOptStat("nemr");
  cr_tree->SetLineColor(kRed);
  c->cd(6);       cr_tree->Draw("eNcl:eFrameID","eNframesBot==0","prof");
  gStyle->SetOptStat("ne");
  c->cd(8);       cr_tree->Draw("eZ0","eNframesBot==0");
  gStyle->SetOptStat("nemr");
  cr_tree->SetLineColor(1);
}

//-------------------------------------------------------
void check_view(TCanvas *c)
{
  c->cd(1); cr_tree->Draw("eY0:eX0","eNframesTop==0&&puls>0","lego2");
  c->cd(2); cr_tree->Draw("eY0:eX0","eNframesTop==0&&puls>0&&(abs(eTy)>.006||abs(eTx)>.006)","colZ");
  c->cd(3); cr_tree->Draw("eY0:eX0","eNframesBot==0&&puls>0","lego2");
  c->cd(4); cr_tree->Draw("eY0:eX0","eNframesBot==0&&puls>0&&(abs(eTy)>.006||abs(eTx)>.006)","colZ");

  c->cd(5);  cr_tree->Draw("segments.eTx","puls>0");
  c->cd(6);  cr_tree->Draw("segments.eTy","puls>0");
  c->cd(7);  cr_tree->Draw("eTy:eTx","eNframesTop==0&&puls>0&&((abs(eTy)>.01||abs(eTx)>.01)) && abs(eTy)<.6 && abs(eTx)<.6","colZ");
  c->cd(8);  cr_tree->Draw("eTy:eTx","eNframesBot==0&&puls>0&&((abs(eTy)>.01||abs(eTx)>.01)) && abs(eTy)<.6 && abs(eTx)<.6","colZ");
}

//-------------------------------------------------------
void check_surf(TCanvas *c)
{
  c->cd(1); {
    cr_tree->Draw("eYview:eXview");
    cr_tree->Draw("eYview+eY0:eXview+eX0","eNframesTop==0&&puls>0","samecontZ");
  }
  c->cd(2);{
    cr_tree->Draw("eYview:eXview");
    cr_tree->Draw("eYview+eY0:eXview+eX0","eNframesBot==0&&puls>0","samecontZ");
  }
  c->cd(3); { 
    cr_tree->Draw("eNsegments");
    cr_tree->SetLineColor(kBlue);
    cr_tree->Draw("eNsegments>>hbot","eNframesTop==0","same");
    cr_tree->SetLineColor(kRed);
    cr_tree->Draw("eNsegments>>htop","eNframesBot==0","same");
  }
  cr_tree->SetLineColor(1);
  c->cd(4); {
    cr_tree->Draw("Sum$(eNcl)");
    cr_tree->SetLineColor(kBlue);
    cr_tree->Draw("Sum$(eNcl)","eNframesTop==0","same");
    cr_tree->SetLineColor(kRed);
    cr_tree->Draw("Sum$(eNcl)","eNframesBot==0","same");
  }
  cr_tree->SetLineColor(1);
  c->cd(5); {
    cr_tree->Draw("puls");
    cr_tree->SetLineColor(kBlue);
    cr_tree->Draw("puls","eNframesTop==0","same");
    cr_tree->SetLineColor(kRed);
    cr_tree->Draw("puls","eNframesBot==0","same");
  }
  cr_tree->SetLineColor(1);
  c->cd(6); {
    cr_tree->Draw("volume");
    cr_tree->SetLineColor(kBlue);
    cr_tree->Draw("volume","eNframesTop==0","same");
    cr_tree->SetLineColor(kRed);
    cr_tree->Draw("volume","eNframesBot==0","same");
  }
  cr_tree->SetLineColor(1);
  c->cd(7); {
    cr_tree->Draw("eSigmaX");
    cr_tree->SetLineColor(kBlue);
    cr_tree->Draw("eSigmaX","eNframesTop==0","same");
    cr_tree->SetLineColor(kRed);
    cr_tree->Draw("eSigmaX","eNframesBot==0","same");
  }
  cr_tree->SetLineColor(1);
  c->cd(8); {
    cr_tree->Draw("eSigmaY");
    cr_tree->SetLineColor(kBlue);
    cr_tree->Draw("eSigmaY","eNframesTop==0","same");
    cr_tree->SetLineColor(kRed);
    cr_tree->Draw("eSigmaY","eNframesBot==0","same");
  }
  cr_tree->SetLineColor(1);
}

//-------------------------------------------------------
void check_raw( int output=0, char *fname=0 )
{

   if(fname) cout << "fname  : " << fname  << endl;

   printf("Red  line: Top    side\n");
   printf("Blue line: Bottom side\n");

  if      (fname&&output>=2)   run = new EdbRun( fname , "UPDATE");
  else if (fname&&output<2)    run = new EdbRun( fname );
  
  if(run) cr_tree = run->GetTree() ;
  else    cr_tree = (TTree*)gDirectory->Get("Views");

  init();

  int do_z=1 , do_view=1 , do_surf=1 ;
  
  TCanvas *cz=0, *cview=0, *csurf=0;
  if( output!=2 ) {
    cz    = dynamic_cast<TCanvas*>(gDirectory->Get("raw_z"));
    cview = dynamic_cast<TCanvas*>(gDirectory->Get("raw_view"));
    csurf = dynamic_cast<TCanvas*>(gDirectory->Get("raw_surf"));
    if(cz)    { cz->Draw();    do_z=0;    }
    if(cview) { cview->Draw(); do_view=0; }
    if(csurf) { csurf->Draw(); do_surf=0; }
  }

  if(do_z) {
    cz    = new TCanvas("cz"    , "check Z", XSIZE,YSIZE);
    cz->Divide(2,4);
    check_z(cz);
  }
  if(do_view)  {
    cview = new TCanvas("view"  , "check View", XSIZE,YSIZE);
    cview->Divide(2,4);
    check_view(cview);
  }
  if(do_surf)  {
    csurf = new TCanvas("csurf" , "check surf, etc",XSIZE,YSIZE);
    csurf->Divide(2,4);
    check_surf(csurf);
  }

   if( (fname&&output==2) || (fname&&output==3) ) {
      printf("save as canvases into root file\n");
      if(cz)    cz->Write("raw_z");
      if(cview) cview->Write("raw_view");
      if(csurf) csurf->Write("raw_surf");
   } 
   if(output==1  || (fname&&output==3) ) {
      printf("save as gif pictures\n");
      gSystem->Sleep(500);
      if(cz)    cz->SaveAs("raw_z.gif");
      if(cview) cview->SaveAs("raw_view.gif");
      if(csurf) csurf->SaveAs("raw_surf.gif");
    }
   if(run) run->Close();
}

//-----------------------------------------------------------------
#ifndef __CINT__
int main( int argc, char *argv[] )
{
   int  output=3 ;
   char* fname=argv[1] ;
   
   gStyle->SetPalette(1);
   check_raw( output,  fname ) ;

   return 1 ;
}
#endif


