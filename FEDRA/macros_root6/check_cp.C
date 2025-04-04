//------------------------------------------------------
// To check the quality of couples data
// usage:
//  a) interactive: 
//      $ root -l pl56.root
//      root[0] .x check_cp.C             draw plots on the screen; use saved 
//                                        into rootfile canvases if any
//      root[0] .x check_cp.C(1)          the pictures will be saved as gif files
//      root[0] .x check_cp.C("s.eTX>0")  the pictures will be saved as gif files 
//                                        with SELECTION "s.eTX>0"
//  b) batch mode:
//    $ root -b -q check_cp.C\(2,\"pl56.root\"\)   run file will be opened for  
//                                                 update to save canvases inside
//    $ root -b -q check_cp.C\(3,\"pl56.root\"\)   run file will be opened for 
//                                                 update to save canvases inside 
//                                                 AND the pictures will be saved 
//                                                 ALSO as gif files
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
#include "EdbPVRec.h"
using namespace std;

#endif
#endif

//declaration of used functions
void init();
void check_surf( TCanvas *surf );
void check_sigma( TCanvas *cs );
void check_view( TCanvas *cs );
void check_shrinkage( TCanvas *diff );
void correct_shrinkage( TCanvas *cshr );
void check_distorsion( TCanvas *cs );

TCut csignal("1");
TCut cut1("pid2>-1&&eCHI2P<1.5");
TCut sameview("(s1.eAid[1]==s2.eAid[1])");
TCut diffview("(s1.eAid[1]!=s2.eAid[1])");
TCut diffarea("(s1.eAid[0]!=s2.eAid[0])");

TTree* couples;

#define XSIZE 1200
#define YSIZE  800

//-----------------------------------------------------------------
void init()
{
  gStyle->SetOptFit(0001);
  couples->SetAlias("tx","(s2.eX-s1.eX)/(s2.eZ-s1.eZ)");          // baseline angle tx
  couples->SetAlias("ty","(s2.eY-s1.eY)/(s2.eZ-s1.eZ)");          // baseline angle ty
  couples->SetAlias("t","sqrt(tx*tx+ty*ty)");                     // baseline theta
  couples->SetAlias("ts","sqrt(s.eTX*s.eTX + s.eTY*s.eTY)");      // basetrack theta

  couples->SetAlias("t1","sqrt(s1.eTX*s1.eTX + s1.eTY*s1.eTY)");  // theta of the s1
  couples->SetAlias("t2","sqrt(s2.eTX*s2.eTX + s2.eTY*s2.eTY)");  // theta of the s2

  couples->SetAlias("phi","atan(s.eTY/s.eTX)");                   // phi of the s
  couples->SetAlias("phi1","atan(s1.eTY/s1.eTX)");                // phi of the s1
  couples->SetAlias("phi2","atan(s2.eTY/s2.eTX)");                // phi of the s2

  couples->SetAlias("dsx1","s1.eTX-s.eTX");
  couples->SetAlias("dsy1","s1.eTY-s.eTY");
  couples->SetAlias("ds1","sqrt(dsx1*dsx1+dsy1*dsy1)");                 // absoulte angular diff    s1-s
  couples->SetAlias("dsx2","s2.eTX-s.eTX");
  couples->SetAlias("dsy2","s2.eTY-s.eTY");
  couples->SetAlias("ds2","sqrt(dsx2*dsx2+dsy2*dsy2)");                 // absoulte angular diff    s2-s

  couples->SetAlias("dst1","(s1.eTX*s.eTY-s1.eTY*s.eTX)/ts");           // transverse slope diff      s1-s
  couples->SetAlias("dsl1","sqrt(ds1*ds1-dst1*dst1)");                  // longitudinal slope diff    s1-s
  couples->SetAlias("dst2","(s2.eTX*s.eTY-s2.eTY*s.eTX)/ts");           // transverse slope diff      s2-s
  couples->SetAlias("dsl2","sqrt(ds2*ds2-dst2*dst2)");                  // longitudinal slope diff    s2-s


  couples->SetAlias("dz","s2.eZ-s1.eZ");
  couples->SetAlias("dx12","s2.eX-(s1.eX+dz*s1.eTX)");       // project s1 to s2 and calc dx
  couples->SetAlias("dy12","s2.eY-(s1.eY+dz*s1.eTY)");       // project s1 to s2 and calc dy
  couples->SetAlias("dx21","(s2.eX-dz*s2.eTX)-s1.eX");       // project s1 to s2 and calc dx
  couples->SetAlias("dy21","(s2.eY-dz*s2.eTY)-s1.eY");       // project s1 to s2 and calc dx

}

//-----------------------------------------------------------------
void check_surf( TCanvas *surf )
{
  printf("check_surf with the cut: %s \n", csignal.GetTitle() );
  
  surf->cd(1);  couples->Draw("s.eY:s.eX",csignal);
  surf->cd(2);  couples->Draw("s.eTY:s.eTX",csignal);
  surf->cd(5); {
    couples->Draw("eCHI2P", csignal );
    couples->SetLineColor(kGreen);
    couples->Draw("eCHI2P", csignal && sameview ,"same");
    couples->SetLineColor(kRed);
    couples->Draw("eCHI2P", csignal && diffview,"same");
    couples->SetLineColor(6);
    couples->Draw("eCHI2P", csignal && diffarea,"same");
  }
  couples->SetLineColor(1);
  surf->cd(3); {
    couples->SetLineColor(kBlue);
    couples->Draw("s1.eW>>hw1(15,5,20)", csignal );
    couples->SetLineColor(kRed);
    couples->Draw("s2.eW", csignal,"same");
  } 
  couples->SetLineColor(1);
  
  surf->cd(4);  couples->Draw("s.eW>>hw(25,10,35)", csignal );
  surf->cd(6);  couples->Draw("eCHI2P:s.eW>>hchiw(25,10,35,30,0,3.)", csignal,"colZ");

}

//-----------------------------------------------------------------
void check_sigma( TCanvas *cs )
{
  printf("check_sigma with the cut: %s \n", csignal.GetTitle() );

  couples->Project("htx1(20,-1.,1.,40,0.,0.1)", "abs(s1.eTX-tx):tx",csignal,"prof");
  couples->Project("htx2(20,-1.,1.,40,0.,0.1)", "abs(s2.eTX-tx):tx",csignal,"prof");
  couples->Project("hty1(20,-1.,1.,40,0.,0.1)", "abs(s1.eTY-ty):ty",csignal,"prof");
  couples->Project("hty2(20,-1.,1.,40,0.,0.1)", "abs(s2.eTY-ty):ty",csignal,"prof");
    TH1* htx1  = (TH1*) gDirectory->Get("htx1");
    TH1* htx2  = (TH1*) gDirectory->Get("htx2");
    TH1* hty1  = (TH1*) gDirectory->Get("hty1");
    TH1* hty2  = (TH1*) gDirectory->Get("hty2");

    cs->cd(1);{
    cs->GetPad(1)->SetGrid(1,1);
    htx1->SetLineColor(kRed);
    htx1->Draw();
    htx2->SetLineColor(kBlue);
    htx2->Draw("same");
  }
  cs->cd(2);{
    cs->GetPad(2)->SetGrid(1,1);
    hty1->SetLineColor(kRed);
    hty1->Draw();
    hty2->SetLineColor(kBlue);
    hty2->Draw("same");
  }
  couples->Project("htt1(10,0.,1.,40,0.,0.1)", "abs(dst1):ts",csignal,"prof");
  couples->Project("htt2(10,0.,1.,40,0.,0.1)", "abs(dst2):ts",csignal,"prof");
  couples->Project("htl1(10,0.,1.,40,0.,0.1)", "dsl1:ts",csignal,"prof");
  couples->Project("htl2(10,0.,1.,40,0.,0.1)", "dsl2:ts",csignal,"prof");
    TH1* htt1  = (TH1*) gDirectory->Get("htt1");
    TH1* htt2  = (TH1*) gDirectory->Get("htt2");
    TH1* htl1  = (TH1*) gDirectory->Get("htl1");
    TH1* htl2  = (TH1*) gDirectory->Get("htl2");

  cs->cd(3); {
    cs->GetPad(3)->SetGrid(1,1);
    htt1->SetLineColor(kRed);
    htt1->Draw();
    htt2->SetLineColor(kBlue);
    htt2->Draw("same");
  }
  cs->cd(4);{
    cs->GetPad(4)->SetGrid(1,1);
    htl1->SetLineColor(kRed);
    htl1->Draw();
    htl2->SetLineColor(kBlue);
    htl2->Draw("same");
    gStyle->SetOptStat("ne");
  }

  cs->cd(5); {
    cs->GetPad(4)->SetGrid(1,1);
    couples->Project("hs" ,"s.eW:ts" , csignal, "prof");
    couples->Project("hs1","s1.eW:ts", csignal, "prof");
    couples->Project("hs2","s2.eW:ts", csignal, "prof");
    TH1* hs  = (TH1*) gDirectory->Get("hs");
    TH1* hs1 = (TH1*) gDirectory->Get("hs1");
    TH1* hs2 = (TH1*) gDirectory->Get("hs2");
    hs->Draw();
    hs1->SetLineColor(kRed);
    hs1->Draw("same");
    hs2->SetLineColor(kBlue);
    hs2->Draw("same");
    gStyle->SetOptStat("nemr");
  }
  cs->cd(6); {
    couples->SetLineColor(kBlue);
    couples->Draw("eN1tot:ts", csignal, "prof");
    couples->SetLineColor(kRed);
    couples->Draw("eN2tot:ts", csignal, "prof same");
  } 
  couples->SetLineColor(1);

}

//-----------------------------------------------------------------
void check_view( TCanvas *cs )
{
  // check the accuracy deterioration in case when the segments are in 
  // the different views

  printf("check_view with the cut: %s \n", csignal.GetTitle() );

  cs->cd(1);
  couples->Draw("eCHI2P", csignal&&sameview);
  couples->Draw("eCHI2P", csignal&&diffview,"same");
  couples->Draw("eCHI2P", csignal&&diffarea,"same");
  cs->cd(2);
  couples->SetMarkerStyle(20);
  couples->SetAlias("same_view","(s1.eAid[1]==s2.eAid[1])");
  couples->Draw("s.eTY:same_view",csignal,"prof");
  cs->cd(3);
  couples->Draw("s.eTX:same_view",csignal,"prof");
  cs->cd(4);
  couples->Draw("s.eTX>>htxv(100)",csignal&&sameview);
  couples->Draw("s.eTX",csignal&&diffview,"same");
  couples->SetMarkerStyle(1);
  gStyle->SetOptStat("nemr");

}

//-----------------------------------------------------------------
void check_shrinkage( TCanvas *diff )
{
  printf("check_shrinkage with the cut: %s \n", csignal.GetTitle() );

  diff->cd(1);
  couples->Draw("s1.eTX-tx:tx", csignal );
  diff->cd(2);
  couples->Draw("s2.eTX-tx:tx", csignal );
  diff->cd(3);
  couples->Draw("s1.eTY-ty:ty", csignal );
  diff->cd(4);
  couples->Draw("s2.eTY-ty:ty", csignal );
  gStyle->SetOptStat("nemr");

}

//-----------------------------------------------------------------
void correct_shrinkage( TCanvas *cshr )
{
   // this function check shrinkage and/or distance between linked planes 
   // Note: do not use s.* (linked segment parameters), because them could 
   // be different from the "base angle" calculated here directly

   printf("correct_shrinkage with the cut: %s \n", csignal.GetTitle() );

   //  TCut cut1("pid2>-1&&eCHI2P<1.5");

   cshr->cd(1);   
   TH1* hsh1  ;

   couples->Draw("s1.eTX-tx:tx>>hsh1", csignal,"prof");
   hsh1 = (TH1*) gDirectory->Get("hsh1");
   hsh1->Fit("pol1","wQ","",-.4,.4);

   float p0 = hsh1->GetFunction("pol1")->GetParameter(0);
   float p1 = hsh1->GetFunction("pol1")->GetParameter(1);
   printf("side1     : p0 = %f \t p1 = %f \n",p0,p1);
   

   cshr->cd(3);   
   TH1* hsh3  ;

   char str[160]="";
   sprintf(str,"s1.eTX*(1-(%f))-(%f)-(s2.eX-s1.eX)/(s2.eZ-s1.eZ):(s2.eX-s1.eX)/(s2.eZ-s1.eZ)>>hsh3",p1,p0);
   couples->Draw(str,csignal,"prof");
   hsh3 = (TH1*) gDirectory->Get("hsh3");
   hsh3->Fit("pol1","wQ","",-.4,.4);

   p0 = hsh3->GetFunction("pol1")->GetParameter(0);
   p1 = hsh3->GetFunction("pol1")->GetParameter(1);
   printf("side1 corr: p0 = %f \t p1 = %f \n",p0,p1);

   
   cshr->cd(2);
   TH1* hsh2;

   couples->Draw("s2.eTX-tx:tx>>hsh2", csignal,"prof");
   hsh2 = (TH1*) gDirectory->Get("hsh2");
   hsh2->Fit("pol1","wQ","",-.4,.4);
  
   p0 = hsh2->GetFunction("pol1")->GetParameter(0);
   p1 = hsh2->GetFunction("pol1")->GetParameter(1);
   printf("side2     : p0 = %f \t p1 = %f \n",p0,p1);

   cshr->cd(4); 
   TH1* hsh4 ;

   sprintf(str,"s2.eTX*(1-(%f))-(%f)-(s2.eX-s1.eX)/(s2.eZ-s1.eZ):(s2.eX-s1.eX)/(s2.eZ-s1.eZ)>>hsh4",p1,p0);
   couples->Draw(str, csignal ,"prof");
   hsh4 = (TH1*) gDirectory->Get("hsh4");   
   hsh4->Fit("pol1","wQ","",-.4,.4);
  
   p0 = hsh4->GetFunction("pol1")->GetParameter(0);
   p1 = hsh4->GetFunction("pol1")->GetParameter(1);
   printf("side2 corr: p0 = %f \t p1 = %f \n",p0,p1);
  
  
   gStyle->SetOptStat("nemr");
}

//-----------------------------------------------------------------
void check_distorsion( TCanvas *cs )
{
  printf("check_distortion with the cut: %s \n", csignal.GetTitle() );

  cs->cd(1);
  couples->Draw("s1.eTX-s.eTX:s.eX", csignal );
  cs->cd(2);
  couples->Draw("s1.eTY-s.eTY:s.eY", csignal );
  cs->cd(3);
  couples->Draw("s2.eTX-s.eTX:s.eX", csignal );
  cs->cd(4);
  couples->Draw("s2.eTY-s.eTY:s.eY", csignal );
  gStyle->SetOptStat("nemr");
}
//-----------------------------------------------------------------
/*void check_all()
{
  check_surf();
  check_sigma();
  check_view();
  check_shrinkage();
  check_distorsion();
  correct_shrinkage();
}
*/
//-----------------------------------------------------------------
void check_cp(int output=0, char *fname=0, char *csignal_cut="1")
{
  csignal = csignal_cut;

  cout << "fname  : " << fname  << endl;

  printf("Check couples with the general selection as: %s\n",csignal_cut);
  printf("s1: Red  line\n");
  printf("s2: Blue line\n");
  printf("functions: check_surf, check_sigma, check_shrinkage, check_distorsion, check_view, correct_shrinkage\n");
  
  TFile *f=0;
  if      (fname&&output>=2)   f = new TFile(fname,"UPDATE");
  else if (fname&&output<2)    f = new TFile(fname);
  
  couples = (TTree*) gFile->Get("couples") ;

  init();
  //check_all();

  int do_surf=1, do_sig=1, do_view=1, do_diff=1;
  int do_shr=1;
  int do_dist=1 ;
  
  TCanvas *csurf=0, *csig=0, *cview=0, *cdiff=0, *cshr=0, *cdist=0;
  if( output!=2 ) {
    csurf   = dynamic_cast<TCanvas*>(gDirectory->Get("csurf"));
    csig    = dynamic_cast<TCanvas*>(gDirectory->Get("csig"));
    cview   = dynamic_cast<TCanvas*>(gDirectory->Get("cview"));
    cdiff   = dynamic_cast<TCanvas*>(gDirectory->Get("cdiff"));
    cshr    = dynamic_cast<TCanvas*>(gDirectory->Get("cshr"));
    cdist   = dynamic_cast<TCanvas*>(gDirectory->Get("cdist"));

    if(csurf ) { csurf ->Draw();  do_surf =0; }
    if(csig  ) { csig  ->Draw();  do_sig  =0; }
    if(cview ) { cview ->Draw();  do_view =0; }
    if(cdiff ) { cdiff ->Draw();  do_diff =0; }
    if(cshr  ) { cshr  ->Draw();  do_shr  =0; }
    if(cdist ) { cdist ->Draw();  do_dist =0; }
  }

  if(do_surf ) {
    csurf  = new TCanvas("csurf" ,"couples_surf"        ,XSIZE,YSIZE);
    csurf->Divide(2,3);
    check_surf(csurf);
  }
  if(do_sig  ) {
    csig   = new TCanvas("csig"  ,"couples_sigma"       ,XSIZE,YSIZE);
    csig ->Divide(2,3);
    check_sigma(csig);
  }
  if(do_view ) {
    cview  = new TCanvas("cview" ,"couples_view"        ,XSIZE,YSIZE);
    cview->Divide(2,2);
    check_view(cview);
  }
  if(do_diff ) {
    cdiff  = new TCanvas("cdiff" ,"couples_shrinkage"   ,XSIZE,YSIZE);
    cdiff->Divide(2,2);
    check_shrinkage(cdiff);
  }
  if(do_shr  ) {
    cshr   = new TCanvas("cshr"  ,"couples_shrinkage_corr",XSIZE,YSIZE);
    cshr ->Divide(2,2);
    correct_shrinkage(cshr);
  }
  if(do_dist ) {
    cdist  = new TCanvas("cdist" ,"couples_distortion"  ,XSIZE,YSIZE);
    cdist->Divide(2,2);
    check_distorsion(cdist) ;
  }

  if( (fname&&output==2) || (fname&&output==3) ) {
    printf("save as canvases into root file\n");
    if(csurf ) csurf ->Write("cp_surf");
    if(csig  ) csig  ->Write("cp_sigma" ); 
    if(cview ) cview ->Write("cp_view");
    if(cdiff ) cdiff ->Write("cp_diff");
    if(cshr  ) cshr  ->Write("cp_shr" ); 
    if(cdist ) cdist ->Write("cp_dist");
  } 
  if(output==1  || (fname&&output==3) ) {
    printf("save as gif pictures\n");
    gSystem->Sleep(500);
    if(csurf ) csurf ->SaveAs("cp_surf.gif" );
    if(csig  ) csig  ->SaveAs("cp_sigma.gif"); 
    if(cview ) cview ->SaveAs("cp_view.gif");
    if(cdiff ) cdiff ->SaveAs("cp_diff.gif");
    if(cshr  ) cshr  ->SaveAs("cp_shr.gif" ); 
    if(cdist ) cdist ->SaveAs("cp_dist.gif");
  }
  if(f) f->Close();
}

//-----------------------------------------------------------------
/*void check_cp(char *csignal_cut="1") //ROOT6 do not like redundant methods
{
   check_cp( 0,  0, csignal_cut ) ;
}*/

//-----------------------------------------------------------------
#ifndef __CINT__
int main( int argc, char *argv[] )
{          
   int output=3;
   char* fname=argv[1];
   char* csignal_cut="1";

   gStyle->SetPalette(1);
   check_cp( output,  fname, csignal_cut ) ;

   return 1 ;
}
#endif


