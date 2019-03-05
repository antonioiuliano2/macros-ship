//time to get cracking this (13 aprile 2018)

#include <stdio.h>
#include <TROOT.h>
#include "/afs/cern.ch/work/a/aiuliano/public/SHIPBuild/FairShip/shipdata/ShipMCTrack.h"
#include "/afs/cern.ch/work/a/aiuliano/public/SHIPBuild/FairShip/charmdet/BoxPoint.h"
#include "/afs/cern.ch/work/a/aiuliano/public/SHIPBuild/FairShip/charmdet/SpectrometerPoint.h"
//#include "/home/utente/SHIPBuild/sw/ubuntu1604_x86-64/ROOT/v6-10-06-ship-1/include/TDatabasePDG.h"
//#include "/home/utente/SHIPBuild/sw/ubuntu1604_x86-64/FairRoot/Oct17-ship-1/include/FairMCEventHeader.h"*/
//#include "/home/utente/fedra/include/EdbCouplesTree.h"
using namespace TMath;

void fromFairShip2Fedra(int nplate = 1){
 //TFile * inputfile = TFile::Open("/afs/cern.ch/work/a/aiuliano/public/sim_charm/pot/Charm1_uniform_19_11_2018/pythia8_evtgen_Geant4_1000_0.001.root");
TFile * inputfile = TFile::Open("/eos/user/a/aiuliano/sims_FairShip/sim_charm/pot/uniformonespill_onelayer_ch1_03_03_19/pythia8_Geant4_1000_0.5.root");
 TTree* cbmsim = (TTree*)inputfile->Get("cbmsim");
 gInterpreter->AddIncludePath("/afs/cern.ch/work/a/aiuliano/public/fedra/include");
 Double_t tx = 0, ty=0, xem= 0, yem = 0;
 
 TClonesArray *arr0 = new TClonesArray("ShipMCTrack",10000);
 cbmsim->GetBranch("MCTrack")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("MCTrack",&arr0);

 TClonesArray *arr1 = new TClonesArray("BoxPoint",10000);
 cbmsim->GetBranch("BoxPoint")->SetAutoDelete(kFALSE);
 cbmsim->SetBranchAddress("BoxPoint",&arr1);

 const Int_t nevents = cbmsim->GetEntries();

 // TFile * outputfile = new TFile("myfedra.out","RECREATE");
   EdbCouplesTree ect;
 //  EdbCouplesTree ectpixel;
  if (nplate <10) ect.InitCouplesTree("couples",Form("b000001/p00%i/1.%i.0.0.cp.root",nplate,nplate),"RECREATE");
  else ect.InitCouplesTree("couples",Form("b000001/p0%i/1.%i.0.0.cp.root",nplate,nplate),"RECREATE");
   // ectpixel.InitCouplesTree("couples","pixel.root","RECREATE");
   Int_t Flag = 1;
 for (int i = 0; i < nevents; i++){
   arr0->Clear();
   arr1->Clear();
   cbmsim->GetEntry(i);
   for (int j = 0; j < arr1->GetEntriesFast(); j++){
//no you dont//     if (j % 2 == 0) continue;
     BoxPoint* emupoint = (BoxPoint*)arr1->At(j);
     xem = emupoint->GetX()* 1E+4 + 62500;
     yem = emupoint->GetY()* 1E+4 + 49500;
     tx = emupoint->GetPx()/emupoint->GetPz();
     ty = emupoint->GetPy()/emupoint->GetPz();    
     if (emupoint->GetDetectorID() == nplate){
       ect.eS->Set(i,xem,yem,tx,ty,1,Flag);
       ect.eS->SetW(70.); //need a high weight to do tracking
     //if(do_invert) ect.eS->Transform(&aff_invert);
       ect.Fill();
       }
     }
   }

/*   for (int j = 0; j < arr2->GetEntriesFast(); j++){
     SpectrometerPoint* spectropoint = (SpectrometerPoint*)arr2->At(j);
     xem = spectropoint->GetX() * 1E+4; //fedra units are mum, you moron!
     yem = spectropoint->GetY() * 1E+4;
     tx = spectropoint->GetPx()/spectropoint->GetPz();
     ty = spectropoint->GetPy()/spectropoint->GetPz();
     if (spectropoint->GetDetectorID() == 101) {
     ectpixel.eS->Set(i,xem,yem,tx,ty,1,Flag);
     //if(do_invert) ect.eS->Transform(&aff_invert);
     ectpixel.Fill();
     }
   }*/ 
  ect.Close();
  //ectpixel.Close();
}
/*
void ReadGEMdata(Int_t run, Int_t nCS1, Int_t nCS2, const char *detType="nameofdetector", bool do_invert=true)
{
  const char filename[50], filename1[50], filename2[50], filename3[50]; 
  const char outname[50];
  
  char *det1 = "ECC", *det2 = "CES";
  if(strcmp(detType,det1)==0 || strcmp(detType,det2)==0)
  {
    sprintf(filename, "run_%d_%s_GEM.txt",run, detType);
    sprintf(filename1, "run_%d_%s_GEM_TRK_0.txt",run, detType);
    sprintf(filename2, "run_%d_%s_GEM_TRK_4.txt",run, detType);
    sprintf(filename3, "run_%d_%s_GEM_TRK_0_TRK_4.txt",run, detType);
    sprintf(outname, "cp_%d_%s.root",run, detType);
  }
  if(strcmp(detType,det1)!=0 && strcmp(detType,det2)!=0)
  {
    sprintf(filename, "run_%d_cs_%d_%d_GEM.txt",run, nCS1, nCS2);
    sprintf(filename1, "run_%d_cs_%d_%d_GEM_TRK_0.txt",run, nCS1, nCS2);
    sprintf(filename2, "run_%d_cs_%d_%d_GEM_TRK_4.txt",run, nCS1, nCS2);
    sprintf(filename3, "run_%d_cs_%d_%d_GEM_TRK_0_TRK_4.txt",run, nCS1, nCS2);
    sprintf(outname, "cp_%d_cs%d%d.root",run, nCS1, nCS2);
  }
  
  cout << filename << endl;
  
  FILE *f = fopen(filename,"r");
  if(!f) cout << "No file found!" << endl;
  cout << filename1 << endl;
  FILE *f1 = fopen(filename1,"r");
  if(!f1) cout << "No file found!" << endl;
  cout << filename2 << endl;
  FILE *f2 = fopen(filename2,"r");
  if(!f2) cout << "No file found!" << endl;
  cout << filename3<< endl;
  FILE *f3 = fopen(filename3,"r");
  if(!f3) cout << "No file found!" << endl;

  Int_t Flag[100000];
  Double_t posx[100000], posy[100000], posyW[100000];
  char title1[20], title2[20], title3[20], title4[20], title5[20];
  Int_t row=0, event=0;
  Double_t y =0., x=0., yWell=0.;
  
  for(Int_t i=0;i< 100000;i++)
    {
        Flag[i] =0;
        posx[i] =0.; posy[i]=0.; posyW[i]=0;       
    }       
  
  Int_t firstline=0;
  while(!feof(f))
  {
    if(firstline==0)
    {
      fscanf(f,"%s %s %s %s %s", title1, title2, title3, title4, title5);
      firstline++;
    }
    else
    {
      fscanf(f,"%d %d %lf %lf %lf", &row, &event, &posy[row], &posx[row], &posyW[row]);
      Flag[row]=1000;
    }
  }
  cout << "Here!" << endl;
  firstline =0;
  if(f1)
  {
    while(!feof(f1))
    {
      if(firstline==0)
        {
        fscanf(f1,"%s %s %s %s %s", title1, title2, title3, title4, title5);
        firstline++;
        }
      else
      {
        fscanf(f1,"%d %d %lf %lf %lf", &row,&event, &y, &x, &yWell);
        Flag[row]=1100;
      }
    }
  }
  firstline =0;
  cout << "Here!" << endl;
  
  if(f2)
  {
    while(!feof(f2))
      {
      if(firstline==0)
        {
        fscanf(f2,"%s %s %s %s %s", title1, title2, title3, title4, title5);
        firstline++;
      }
      else
      {
        fscanf(f2,"%d %d %lf %lf %lf", &row,&event, &y, &x, &yWell);
        Flag[row]=1110;
      }
    }
  }
  firstline =0;
  cout << "Here!" << endl;
  
  if(f3)
  {
    while(!feof(f3))
    {
      if(firstline==0)
      {
        fscanf(f3,"%s %s %s %s %s", title1, title2, title3, title4, title5);
        firstline++;
        }
      else
      {
        fscanf(f3,"%d %d %lf %lf %lf", &row,&event, &y, &x, &yWell);
        Flag[row]=1111;
      }
    }
  }
  firstline =0;
  cout << "Here!" << endl;
 
  
  EdbCouplesTree ect;
  ect.InitCouplesTree("couples",outname,"RECREATE");

  Double_t tx = 0, ty=0;  
  if(run==1080 || run==1085 || run==1094 || run==1101) {tx=0.; ty=0.;}
  if(run==1081 || run==1086 || run==1095 || run==1102) {tx=0.; ty=Tan(0.131);}
  if(run==1082 || run==1087 || run==1096 || run==1103) {tx=0.; ty=Tan(0.262);}
  if(run==1083 || run==1088 || run==1097 || run==1105) {tx=0.; ty=Tan(Pi()/6);}
  if(run==1084 || run==1089 || run==1098 || run==1106) {tx=0.; ty=Tan(Pi()/4);}
  
  EdbAffine2D aff_invert(-1,0,0,1,120000,0);
 
  for(Int_t i=0;i<100000;i++)
  {
    //cout <<i <<endl;
    if(posx[i]>-998&&posy[i]>-998) 
    {
      if(Flag[i]!=0)
      {       
        float xem = (posx[i]+60)*1000.;
        float yem =  (-posy[i]+50)*1000.;
        ect.eS->Set(i,xem,yem,tx,ty,1,Flag[i]);
        if(do_invert) ect.eS->Transform(&aff_invert);
        ect.Fill();
      }       
    }
  }
  
  ect.Close();
  fclose(f);
  fclose(f1);
  fclose(f2);
  fclose(f3);
}
*/
