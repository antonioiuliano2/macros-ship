/*EdbAffine2D invertaff(int nplate, int refplate=29);

EdbAffine2D invertaff(int nplate, int refplate){ 

 TFile *setfile = TFile::Open("b000001.0.0.0.set.root");
 EdbScanSet *set = (EdbScanSet*) setfile->Get("set");

 EdbAffine2D afftoplate = EdbAffine2D();
 
 set->GetAffP2P(nplate,refplate,afftoplate);
 cout<<"Test trasformazione da 0 a 20"<<endl;

 afftoplate.Print();

 EdbAffine2D invertedaff = afftoplate.Invert();

 return invertedaff;

}*/

void getlocalpos(float x, float y, int nplate= 10, int refplate=29){
//from GLOBAL to LOCAL plate coordinates

 TFile *setfile = TFile::Open("b000001.0.0.0.set.root");
 EdbScanSet *set = (EdbScanSet*) setfile->Get("set");

 EdbAffine2D afftoplate = EdbAffine2D();
 
 set->GetAffP2P(nplate,refplate,afftoplate);
 cout<<"Test trasformazione da 0 a 20"<<endl;

 afftoplate.Print();

 afftoplate.Invert(); //it does not produce anything new
 
 float newx = afftoplate.Xtrans(x,y);
 float newy = afftoplate.Ytrans(x,y);

 cout<<"localpositions: ("<<newx<<" , "<<newy<<")"<<endl;

}
