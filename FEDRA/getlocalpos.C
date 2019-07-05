//script to get local plate positions from global ones

//global variables
TFile *setfile = 0;
EdbScanSet *set = 0;
EdbAffine2D afftoplate;

//********used functions************************+
void GetScanSetfromFile(){

 setfile = TFile::Open("/ship/CHARM2018/CH1-R6/b000001/b000001.0.0.0.set.root");
 set = (EdbScanSet*) setfile->Get("set");
}

EdbAffine2D invertaff(int nplate, int refplate){ 

 EdbAffine2D afftoplate = EdbAffine2D();
 
 set->GetAffP2P(nplate,refplate,afftoplate);
  cout<<"Print affine transformation from "<<nplate<<" to "<<refplate<<endl;

 afftoplate.Print();

 afftoplate.Invert();

 return afftoplate;

}

void getaffine(int nplate =10, int refplate = 29){
 if (!setfile) GetScanSetfromFile();
 afftoplate= invertaff(nplate,refplate);
}

void globaltolocal(float x, float y){

 float newx = afftoplate.Xtrans(x,y);
 float newy = afftoplate.Ytrans(x,y);

 cout<<"localpositions: ("<<newx<<" , "<<newy<<")"<<endl;
}

//*********automatic script which does all the steps*********

void getlocalpos(float x, float y, int nplate= 10, int refplate=29){
//from GLOBAL to LOCAL plate coordinates

 GetScanSetfromFile();
 afftoplate= getaffine(nplate,refplate);
 globaltolocal(x,y);

 setfile->Close();

}
