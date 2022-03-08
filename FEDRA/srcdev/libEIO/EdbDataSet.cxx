//-- Author :  Valeri Tioukov   9.06.2003

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDataSet                                                           //
//                                                                      //
// scanned raw data set and correspondent parameters                    //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TEventList.h"
#include "TMatrix.h"
#include "EdbDataSet.h"
#include "EdbSegment.h"
#include "EdbCluster.h"
#include "EdbMath.h"
#include "EdbTraceBack.h"
#include "EdbLog.h"
#include "EdbVertex.h"
#include <TSystem.h>
#include <iostream>

using namespace std;

//ClassImp(EdbMask)
ClassImp(EdbDataPiece)
ClassImp(EdbDataSet)
ClassImp(EdbDataProc)


///==============================================================================
EdbDataPiece::EdbDataPiece()
{
  for(int i=0; i<3; i++) eLayers[i]=0;
  Set0();
}

///______________________________________________________________________________
EdbDataPiece::EdbDataPiece(int plate, int piece, char* file, int flag)
{
  for(int i=0; i<3; i++) eLayers[i]=0;
  Set0();
  ePlate=plate;
  ePiece=piece;
  AddRunFile(file);
  eFlag=flag;
}

///______________________________________________________________________________
EdbDataPiece::~EdbDataPiece()
{
   int i;
   for(i=0; i<3; i++)  if(eLayers[i])   {delete eLayers[i]; eLayers[i]=0;}
   for(i=0; i<3; i++)  if(eAreas[i])    {delete eAreas[i];  eAreas[i] =0;}
   for(i=0; i<3; i++)  if(eCond[i])     {delete eCond[i];   eCond[i]  =0;}
   for(i=0; i<3; i++)  if(eCuts[i])     {delete eCuts[i];   eCuts[i]  =0;}
   for(i=0; i<3; i++)  if(eRCuts[i])    {delete eRCuts[i];  eRCuts[i] =0;}
   if( eCouplesInd) { delete eCouplesInd; eCouplesInd=0;}
   if( eEraseMask ) { delete eEraseMask;  eEraseMask=0; }
   eRunFiles.Delete();
   CloseRun();
   CloseCPData();
}

///______________________________________________________________________________
void EdbDataPiece::Set0()
{
  eOUTPUT=0;
  ePlate=0;
  ePiece=0;
  eFlag=0;
  eAFID=0;
  eCLUST=0;
  eCutCP[0]=-1;
  int i;
  for(i=0; i<3; i++) { 
    if(eLayers[i]) delete eLayers[i];
    eLayers[i]=new EdbLayer();
  }
  for(i=0; i<3; i++) eCond[i]=0;
  for(i=0; i<3; i++) eAreas[i]=0;
  for(i=0; i<3; i++) eCuts[i]=0;
  for(i=0; i<3; i++) eRCuts[i]=0;
  eRun    = 0;
  eCouplesTree=0;
  eCouplesInd=0;
  eEraseMask=0;
}

///______________________________________________________________________________
void EdbDataPiece::CloseCPData()
{
  if (eCouplesTree) {
    TFile *f = eCouplesTree->GetDirectory()->GetFile();
    if (f) {
      if(f->IsWritable()) eCouplesTree->AutoSave();
      SafeDelete(eCouplesTree);
      SafeDelete(f);
    }
    eCouplesTree=NULL;
  }
}

///______________________________________________________________________________
void EdbDataPiece::CloseRun()
{
  if(eRun) delete eRun;
  eRun=0;
}

///______________________________________________________________________________
int EdbDataPiece::InitCouplesInd()
{
  if(eCouplesInd) delete eCouplesInd;
  eCouplesInd = new TIndexCell();

  if(!eCouplesTree) if(!InitCouplesTree("READ")) return 0;
  EdbSegP         *s1  = 0;
  EdbSegP         *s2  = 0;
  eCouplesTree->SetBranchAddress("s1."  , &s1  );
  eCouplesTree->SetBranchAddress("s2."  , &s2  );

  Long_t v[5]; // side:aid:vid:sid:entry

  int nentr = (int)(eCouplesTree->GetEntries());

  for(int i=0; i<nentr; i++) {
    eCouplesTree->GetEntry(i);
    v[0] = 1;
    v[1] = s1->Aid(0);
    v[2] = s1->Aid(1);
    v[3] = s1->Vid(1);
    v[4] = i;
    eCouplesInd->Add(5,v);
    v[0] = 2;
    v[1] = s2->Aid(0);
    v[2] = s2->Aid(1);
    v[3] = s2->Vid(1);
    v[4] = i;
    eCouplesInd->Add(5,v);
  }
  eCouplesInd->Sort();
  return eCouplesInd->N();
}

///______________________________________________________________________________
int EdbDataPiece::GetLinkedSegEntr(int side, int aid, int vid, int sid, TArrayI &entr) const
{
  if(!eCouplesInd) return 0;
  Long_t v[4];
  v[0]=side;  v[1]=aid;  v[2]=vid;  v[3]=sid;
  TIndexCell *c = eCouplesInd->Find(4,v);
  if(!c) return 0;
  int n = c->N();
  entr.Set(n);
  for(int i=0; i<n; i++) entr.AddAt( (int)(c->At(i)->Value()), i);
  return n;
}

///______________________________________________________________________________
void EdbDataPiece::Print()
{
  printf("Piece: %s\n",GetName());
  printf("%d %d \n", ePlate,ePiece);
  int i;
  for(i=0; i<Nruns(); i++)
    printf("%s\n",GetRunFile(i));
  for(i=0; i<3; i++)  if(eLayers[i])  eLayers[i]->Print();
  for(i=0; i<3; i++)  if(eCond[i])    eCond[i]->Print();
  for(i=0; i<3; i++)
    if(eCuts[i])
      for(int j=0; j<NCuts(i); j++)  GetCut(i,j)->Print();
}

///______________________________________________________________________________
void EdbDataPiece::AddRunFile( const char *name )
{
  TObjString *str = new TObjString(name);
  eRunFiles.Add(str);
}

///______________________________________________________________________________
const char *EdbDataPiece::GetRunFile( int i ) const
{
  if(Nruns()<i+1) return 0;
  return ((TObjString *)eRunFiles.At(i))->GetName();
}

///______________________________________________________________________________
const char *EdbDataPiece::MakeName()
{
  char name[8];
  sprintf(name,"%2.2d_%3.3d", ePlate,ePiece);
  SetName(name);
  return GetName();
}

///______________________________________________________________________________
const char *EdbDataPiece::MakeNameCP(const char *dir)
{
  eFileNameCP=dir;
  eFileNameCP+=MakeName();
  eFileNameCP+=".cp.root";
  return eFileNameCP.Data();
}

///______________________________________________________________________________
EdbLayer *EdbDataPiece::GetMakeLayer(int id)
{
  if(id<0) return 0;
  if(id>2) return 0;
  if(!GetLayer(id))  eLayers[id] = new EdbLayer();
  return GetLayer(id);
}

///______________________________________________________________________________
void EdbDataPiece::AddSegmentCut(int layer, int xi, float var[10])
{
  if(!eCuts[layer])  eCuts[layer] = new TObjArray();
  eCuts[layer]->Add( new EdbSegmentCut(xi,var) );
}

///______________________________________________________________________________
void EdbDataPiece::AddSegmentCut(int layer, int xi, float min[5], float max[5])
{
  if(!eCuts[layer])  eCuts[layer] = new TObjArray();
  EdbSegmentCut *cut=new EdbSegmentCut();
  cut->SetXI(xi);
  cut->SetMin(min);
  cut->SetMax(max);
  eCuts[layer]->Add( cut );
}

///______________________________________________________________________________
void EdbDataPiece::AddRCut(int layer, TCut &cut)
{
  if(layer<0) return;
  if(layer>2) return;
  
  if(!eRCuts[layer]) eRCuts[layer] = new TCut(cut);
  else (*(eRCuts[layer]))+=cut;
}

///______________________________________________________________________________
void EdbDataPiece::AddCutCP(float var[6])
{
  for(int i=0; i<6; i++) eCutCP[i]=var[i];
}

///______________________________________________________________________________
int EdbDataPiece::NCuts(int layer)
{
  if(!eCuts[layer])  return 0;
  return eCuts[layer]->GetEntriesFast();
}

///______________________________________________________________________________
EdbScanCond *EdbDataPiece::GetMakeCond(int id)
{
  if(id<0) return 0;
  if(id>2) return 0;
  if(!GetCond(id))  eCond[id] = new EdbScanCond();
  return GetCond(id);
}


///______________________________________________________________________________
void EdbDataPiece::MakeNamePar(const char *dir)
{
  eFileNamePar=dir;
  eFileNamePar += MakeName();
  eFileNamePar += ".par";
}

///______________________________________________________________________________
int EdbDataPiece::TakePiecePar()
{
  return ReadPiecePar( eFileNamePar.Data() );
}

///______________________________________________________________________________
int EdbDataPiece::ReadPiecePar(const char *file)
{
  // read parameters from par-file
  // return: 0 if ok
  //        -1 if file access failed

  char buf[512];
  char key[256];
  char name[512];

  FILE *fp = fopen(file,"r");
  if (!fp) {
    Log(2,"ReadPiecePar","ERROR open file: %s", file);
    return -1;
  }
  else Log(3,"ReadPiecePar","Read piece parameters from file: %s", file );

  int id,mode;
  float z,zmin,zmax,shr;
  float a11,a12,a21,a22,b1,b2;
  float x1,x2,x3,x4;
  float var[10];

  while (fgets(buf, sizeof(buf), fp)) {
    for (Int_t i = 0; i < (Int_t)strlen(buf); i++) 
      if (buf[i]=='#')  {
	buf[i]='\0';                       // cut out comments starting from #
	break;
      }
    
    if( sscanf(buf,"%s",key)!=1 )                             continue;

    if      ( !strcmp(key,"INCLUDE")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	ReadPiecePar(name);
      }
    else if ( !strcmp(key,"ZLAYER")   )
      {
	sscanf(buf+strlen(key),"%d %f %f %f",&id,&z,&zmin,&zmax);
	GetMakeLayer(id)->SetZlayer(z,zmin,zmax);
      }
    else if ( !strcmp(key,"SHRINK")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&shr);
	GetMakeLayer(id)->SetShrinkage(shr);
      }
    else if ( !strcmp(key,"AFFXY")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f",&id,&a11,&a12,&a21,&a22,&b1,&b2);
	GetMakeLayer(id)->SetAffXY(a11,a12,a21,a22,b1,b2);
      }
    else if ( !strcmp(key,"AFFTXTY")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f",&id,&a11,&a12,&a21,&a22,&b1,&b2);
	GetMakeLayer(id)->SetAffTXTY(a11,a12,a21,a22,b1,b2);
      }
    else if ( !strcmp(key,"SIGMA0")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f",&id,&x1,&x2,&x3,&x4);
	GetMakeCond(id)->SetSigma0(x1,x2,x3,x4);
      }
    else if ( !strcmp(key,"DEGRAD")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&x1);
	GetMakeCond(id)->SetDegrad(x1);
      }
    else if ( !strcmp(key,"BINS")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f",&id,&x1,&x2,&x3,&x4);
	GetMakeCond(id)->SetBins(x1,x2,x3,x4);
      }
    else if ( !strcmp(key,"RAMP0")  )
      {
	sscanf(buf+strlen(key),"%d %f %f",&id,&x1,&x2);
	GetMakeCond(id)->SetPulsRamp0(x1,x2);
      }
    else if ( !strcmp(key,"RAMP04")  )
      {
	sscanf(buf+strlen(key),"%d %f %f",&id,&x1,&x2);
	GetMakeCond(id)->SetPulsRamp04(x1,x2);
      }
    else if ( !strcmp(key,"CHI2MAX")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&x1);
	GetMakeCond(id)->SetChi2Max(x1);
      }
    else if ( !strcmp(key,"CHI2PMAX")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&x1);
	GetMakeCond(id)->SetChi2PMax(x1);
      }
    else if ( !strcmp(key,"CHI2MODE")  )
      {
	sscanf(buf+strlen(key),"%d %d",&id,&mode);
	GetMakeCond(id)->SetChi2Mode(mode);
      }
    else if ( !strcmp(key,"OFFSET")  )
      {
	sscanf(buf+strlen(key),"%d %f %f",&id,&x1,&x2);
	GetMakeCond(id)->SetOffset(x1,x2);
      }
    else if ( !strcmp(key,"SIGMAGR")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f",&id,&x1,&x2,&x3,&x4);
	GetMakeCond(id)->SetSigmaGR(x1,x2,x3);
	SetCutGR(x4);
      }
    else if ( !strcmp(key,"RADX0")  )
      {
	sscanf(buf+strlen(key),"%d %f",&id,&x1);
	GetMakeCond(id)->SetRadX0(x1);
      }
    else if ( !strcmp(key,"XCUT")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f %f %f %f %f",&id,
	       var,var+1,var+2,var+3,var+4,var+5,var+6,var+7,var+8,var+9);
	AddSegmentCut(id,0,var);
      }
    else if ( !strcmp(key,"ICUT")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f %f %f %f %f",&id,
	       var,var+1,var+2,var+3,var+4,var+5,var+6,var+7,var+8,var+9);
	AddSegmentCut(id,1,var);
      }
    else if ( !strcmp(key,"RCUT")  )
      {
	char rcut[256];
	sscanf(buf+strlen(key),"%d %s",&id, rcut );
	TCut cut(rcut);
	AddRCut(id,cut);
      }
    else if ( !strcmp(key,"CUTCP")  )
      {
	sscanf(buf+strlen(key),"%d %f %f %f %f %f %f",&id,
	       var,var+1,var+2,var+3,var+4,var+5);
	AddCutCP(var);
      }
    else if ( !strcmp(key,"AFID")  )
      {
	sscanf(buf+strlen(key),"%d",&id);
	eAFID=id;
      }
    else if ( !strcmp(key,"OUTPUT")  )
      {
	sscanf(buf+strlen(key),"%d",&id);
	eOUTPUT=id;
      }
    else if ( !strcmp(key,"VOLUME0")  )
      {
	float x0=0,y0=0,z0=0,tx=0,ty=0;
	sscanf(buf+strlen(key),"%f %f %f %f %f",&x0,&y0,&z0,&tx,&ty);
	SetVolume0(x0,y0,z0,tx,ty);
      }
    else if ( !strcmp(key,"VOLUMEA")  )
      {
	float dx=0,dy=0;
	sscanf(buf+strlen(key),"%f %f",&dx,&dy);
	SetVolumeA(dx,dy);
      }
    else if ( !strcmp(key,"CLUST")  )
      {
	int icl=0;
	sscanf(buf+strlen(key),"%d",&icl);
	eCLUST=icl;
      }
  }
  fclose(fp);

  return 0;
}

///______________________________________________________________________________
void EdbDataPiece::SetVolume0(float x0, float y0, float z0, float tx, float ty)
{
  float z = GetLayer(0)->Z();
  GetLayer(0)->SetXY( x0+(z-z0)*tx, y0+(z-z0)*ty );
  GetLayer(0)->SetTXTY(tx,ty);
}

///______________________________________________________________________________
void EdbDataPiece::CorrectShrinkage(int layer, float shr)
{
  GetLayer(layer)->SetShrinkage( shr*GetLayer(layer)->Shr() );
}

///______________________________________________________________________________
int EdbDataPiece::UpdateShrPar(int layer)
{
  const char *file=eFileNamePar.Data();

  FILE *fp=fopen(file,"a");
  if (fp==NULL)   {
    Log(1,"EdbDataPiece::UpdateShrPar","ERROR open file: %s", file);
    return -1;
  }else Log(2,"EdbDataPiece::UpdateShrPar","Update parameters file with SHRINK %d: %s", layer, file );

  char str[64];
  sprintf(str,"SHRINK \t %d \t %f \n",layer, GetLayer(layer)->Shr() );
  fprintf(fp,"\n%s",str);

  fclose(fp);
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::UpdateAffPar(int layer, EdbAffine2D &aff)
{
  const char *file=eFileNamePar.Data();

  FILE *fp=fopen(file,"a");
  if (!fp) {
    Log(1,"EdbDataPiece::UpdateAffPar","ERROR open file: %s", file);
    return -1;
  } else Log(2,"EdbDataPiece::UpdateAffPar","Update parameters file with AFFXY: %s", file );

  char str[124];
  sprintf(str,"AFFXY \t %d \t %f %f %f %f %f %f\n",layer,
	 aff.A11(),aff.A12(),aff.A21(),aff.A22(),aff.B1(),aff.B2() );
  fprintf(fp,"\n%s",str);

  fclose(fp);
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::UpdateZPar(int layer, float z)
{
  const char *file=eFileNamePar.Data();

  FILE *fp=fopen(file,"a");
  if (!fp) {
    Log(1,"EdbDataPiece::UpdateZPar", "ERROR open file: %s", file);
    return -1;
  }
  else Log(2,"EdbDataPiece::UpdateZPar", "Update parameters file with ZLAYER: %s", file );

  char str[124];
  sprintf(str,"ZLAYER \t %d \t %f %f %f\n",layer,
	 z,0.,0. );
  fprintf(fp,"\n%s",str);

  fclose(fp);
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::UpdateAffTPar(int layer, EdbAffine2D &aff)
{

  EdbAffine2D *a = GetLayer(layer)->GetAffineTXTY();
  a->Transform(&aff);

  const char *file=eFileNamePar.Data();

  FILE *fp=fopen(file,"a");
  if (!fp) {
    Log(1,"EdbDataPiece::UpdateAffTPar","ERROR open file: %s", file);
    return -1;
  }
  else Log(2,"EdbDataPiece::UpdateAffTPar","\nUpdate parameters file with AFFTXTY: %s\n\n", file );

  char str[124];
  sprintf(str,"AFFTXTY \t %d \t %f %f %f %f %f %f\n",layer,
	 a->A11(),a->A12(),a->A21(),a->A22(),a->B1(),a->B2() );
  fprintf(fp,"\n%s",str);

  fclose(fp);
  return 1;
}

///______________________________________________________________________________
void EdbDataPiece::WriteCuts()
{
  TString file = eFileNamePar+".C";

  FILE *fp=fopen(file,"w");
  if (!fp) {
    Log(1," EdbDataPiece::WriteCuts","ERROR open file: %s", file.Data());
    return;
  }
  else Log(2,"EdbDataPiece::WriteCuts","Put Cuts to file: %s", file.Data() );

  fprintf(fp,"{\n");

  char str[256];
  for(int i=0; i<3; i++)
    if(eCuts[i])
      for(int j=0; j<NCuts(i); j++)  {
	GetCut(i,j)->CutLine(str,i,j);
	fprintf(fp,"%s",str);
      }

  fprintf(fp,"}\n");
  fclose(fp);
}

///______________________________________________________________________________
int EdbDataPiece::PassCutCP(float var[6])
{
  if(eCutCP[0]<0) return 1;
  for(int i=0; i<6; i++)    if(var[i]>eCutCP[i])  return 0;
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::PassCuts(int id, float var[5])
{
  int nc = NCuts(id);
  for(int i=0; i<nc; i++)
    if( !(GetCut(id,i)->PassCut(var)) )  return 0;
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::TakeRawSegment(EdbView *view, int id, EdbSegP &segP, int side)
{
  EdbSegment *seg = view->GetSegment(id);

  float var[5];
  var[0] = seg->GetX0();
  var[1] = seg->GetY0();
  var[2] = seg->GetTx();
  var[3] = seg->GetTy();
  var[4] = seg->GetPuls();

  if( !PassCuts(side,var) )     return 0;

  float pix, chi2;
  if(eCLUST) {
    pix = GetRawSegmentPix(seg);
    segP.SetVolume( pix );
    chi2 = CalculateSegmentChi2( seg,
				 GetCond(1)->SigmaXgr(),  //TODO: side logic
				 GetCond(1)->SigmaYgr(), 
				 GetCond(1)->SigmaZgr());

    if(chi2>GetCutGR())  return 0;
    segP.SetChi2( chi2 );    
  }

  EdbLayer  *layer=GetLayer(side);
  seg->Transform(layer->GetAffineXY());                          //internal view transformation
  if(eAFID) seg->Transform( view->GetHeader()->GetAffine() );

  float x,y,z,tx,ty,puls;
  tx   = seg->GetTx()/layer->Shr();
  ty   = seg->GetTy()/layer->Shr();
  x    = seg->GetX0() + layer->Zmin()*tx;
  y    = seg->GetY0() + layer->Zmin()*ty;
  z    = layer->Z() + layer->Zmin();
  if(eAFID==0) {
    x+=view->GetXview();
    y+=view->GetYview();
  }
  puls = seg->GetPuls();

  EdbAffine2D *aff = layer->GetAffineTXTY();
  float txx = aff->A11()*tx+aff->A12()*ty+aff->B1();
  float tyy = aff->A21()*tx+aff->A22()*ty+aff->B2();
  segP.Set( seg->GetID(),x,y,txx,tyy,puls,0);
  segP.SetZ( z );
  segP.SetDZ( seg->GetDz()*layer->Shr() );
  segP.SetW( puls );
  segP.SetVolume( seg->GetVolume() );
  segP.SetChi2( seg->GetSigmaY() );   // make sence in case of fedra tracking 
  segP.SetMC( view->GetHeader()->GetEvent(),view->GetHeader()->GetTrack() );
  return 1;
}

///______________________________________________________________________________
float EdbDataPiece::CalculateSegmentChi2( EdbSegment *seg, float sx, float sy, float sz )
{
	//TODO: remove this function from here (already in EdbViewRec)

	//assumed that clusters are attached to segments
  double chi2=0;
  EdbCluster *cl=0;
  TObjArray *clusters = seg->GetElements();
  if(!clusters) return 0;
  int ncl = clusters->GetLast()+1;
  if(ncl<=0)     return 0;

  float xyz1[3], xyz2[3];             // segment line parametrized as 2 points
  float xyz[3];
  bool inside=true;

  xyz1[0] = seg->GetX0() /sx;
  xyz1[1] = seg->GetY0() /sy;
  xyz1[2] = seg->GetZ0() /sz;
  xyz2[0] = (seg->GetX0() + seg->GetDz()*seg->GetTx()) /sx;
  xyz2[1] = (seg->GetY0() + seg->GetDz()*seg->GetTy()) /sy;
  xyz2[2] = (seg->GetZ0() + seg->GetDz())              /sz;

  double d;
  for(int i=0; i<ncl; i++ ) {
    cl = (EdbCluster*)clusters->At(i);
    xyz[0] = cl->GetX()/sx;
    xyz[1] = cl->GetY()/sy;
    xyz[2] = cl->GetZ()/sz;
    d = EdbMath::DistancePointLine3(xyz,xyz1,xyz2, inside);
    chi2 += d*d;
  }

  return TMath::Sqrt(chi2/ncl);
}

///______________________________________________________________________________
float EdbDataPiece::GetRawSegmentPix( EdbSegment *seg )
{
  //assumed that clusters are attached to segments
  float pix=0;
  EdbCluster *cl=0;
  TObjArray *clusters = seg->GetElements();
  if(!clusters) return 0;
  int ncl = clusters->GetLast()+1;
  for(int i=0; i<ncl; i++ ) {
    cl = (EdbCluster*)clusters->At(i);
    pix += cl->GetArea();
  }
  return pix;
}

///______________________________________________________________________________
int EdbDataPiece::TakeCPSegment( EdbSegCouple &cp, EdbSegP &seg)
{
  float var[6];

  var[0] = cp.N1();
  var[1] = cp.N1tot();
  var[2] = cp.N2();
  var[3] = cp.N2tot();
  var[4] = cp.CHI2();
  var[5] = cp.CHI2P();

  if( !PassCutCP(var) )     return 0;

  var[0] = seg.X();
  var[1] = seg.Y();
  var[2] = seg.TX();
  var[3] = seg.TY();
  var[4] = seg.W();

  if( !PassCuts(0,var) )     return 0;

  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::GetCPData( EdbPattern *pat, EdbPattern *p1, EdbPattern *p2)
{
  Log(2,"EdbDataPiece::GetCPData","z = %f", GetLayer(0)->Z());
  pat->SetID(0);
  EdbSegP    segP;

  TTree *tree=eCouplesTree;
  EdbSegCouple    *cp = 0;
  EdbSegP         *s1 = 0;
  EdbSegP         *s2 = 0;
  EdbSegP         *s  = 0;

  TBranch *b_cp=0, *b_s=0, *b_s1=0, *b_s2=0;
  b_cp = tree->GetBranch("cp");
  b_s  = tree->GetBranch("s.");
  b_s1 = tree->GetBranch("s1.");
  b_s2 = tree->GetBranch("s2.");

  b_cp->SetAddress( &cp  );
  b_s->SetAddress(  &s   );
  b_s1->SetAddress( &s1  );
  b_s2->SetAddress( &s2  );


  int nseg = 0;
  int nentr = (int)(tree->GetEntries());
  for(int i=0; i<nentr; i++ ) {
    tree->GetEntry(i);
    b_cp->GetEntry(i);
    b_s->GetEntry(i);
    if( !TakeCPSegment(*cp,*s) )      continue;
    if(pat) {
      s->SetZ( s->Z() + pat->Z() );   /// TO CHECK !!!
      //s->SetPID( ePlate*10 );       /// TO CHECK !!!
      EdbTraceBack::SetBaseTrackVid( *s, ePlate, ePiece, i );
      s->SetChi2(cp->CHI2P());
      pat->AddSegment( *s  );
      nseg++;
    }
    if(p1)  { 
      b_s1->GetEntry(i);                             // !!!
      s1->SetZ( s1->Z() + pat->Z() );
      //s1->SetPID( ePlate*10 +1 );    /// TO CHECK !!!
      p1->AddSegment(  *s1 ); 
      nseg++; 
    }
    if(p2)  { 
      b_s2->GetEntry(i);                             // !!!
      s2->SetZ( s2->Z() + pat->Z() );
      //s2->SetPID( ePlate*10 + 2 );   /// TO CHECK !!!
      p2->AddSegment(  *s2 ); 
      nseg++; 
    }
  }

  Log(2,"EdbDataPiece::GetCPData","%d (of %d) segments are readed", nseg,nentr );

  return nseg;
}

///______________________________________________________________________________
int EdbDataPiece::GetCPData_new( EdbPattern *pat, EdbPattern *p1, EdbPattern *p2, TIndex2 *trseg )
{
  Log(3,"EdbDataPiece::GetCPData_new","z = %f ", GetLayer(0)->Z());
  if(!eCouplesTree)  return  0;
  pat->SetID(0);
  EdbSegP    segP;

  TTree *tree=eCouplesTree;
  EdbSegCouple    *cp = 0;
  EdbSegP         *s1 = 0;
  EdbSegP         *s2 = 0;
  EdbSegP         *s  = 0;

  TBranch *b_cp=0, *b_s=0, *b_s1=0, *b_s2=0;                        // !!!
  b_cp = tree->GetBranch("cp");                   // !!!
  b_s  = tree->GetBranch("s.");                   // !!!
  b_s1 = tree->GetBranch("s1.");                  // !!!
  b_s2 = tree->GetBranch("s2.");                  // !!!

  b_cp->SetAddress( &cp  );
  b_s->SetAddress(  &s   );
  b_s1->SetAddress( &s1  );
  b_s2->SetAddress( &s2  );

  int nseg = 0;
  int nentr = (int)(tree->GetEntries());  if(nentr<1) return 0;

  TCut *cut = GetRCut(0);
  TEventList *lst =0;
  if(cut)       tree->Draw(">>lst", *cut );
  else          tree->Draw(">>lst", "" );
  lst = (TEventList*)(gDirectory->GetList()->FindObject("lst"));

  if(!lst) {Log(1,"EdbDataPiece::GetCPData_new","ERROR!: lst do not found! empty couples tree??"); return 0;}

  int nlst =lst->GetN();

  int entr=0;
  for(int i=0; i<nlst; i++ ) {
    entr = lst->GetEntry(i);

    if(trseg) {           //exclude segments participating in tracks
      if( (trseg->Find(ePlate*1000+ePiece,entr) >= 0) )  continue;
    }

    if(eEraseMask) if(eEraseMask->At(entr)) continue;

    b_cp->GetEntry(entr);                            // !!!
    b_s->GetEntry(entr);                             // !!!
    if( !TakeCPSegment(*cp,*s) )      continue;
    if(pat) {
      s->SetZ( s->Z() + pat->Z() );   /// TO CHECK !!!
      //s->SetPID( ePlate*10 );       /// TO CHECK !!!
      s->SetVid(ePlate*1000+ePiece,entr);
      s->SetChi2(cp->CHI2P());
      pat->AddSegment( *s  );
      nseg++;
    }
    if(p1)  { 
      b_s1->GetEntry(entr);                             // !!!
      s1->SetZ( s1->Z() + pat->Z() );
      //s1->SetPID( ePlate*10 +1 );    /// TO CHECK !!!
      p1->AddSegment(  *s1 ); 
      nseg++; 
    }
    if(p2)  { 
      b_s2->GetEntry(entr);                             // !!!
      s2->SetZ( s2->Z() + pat->Z() );
      //s2->SetPID( ePlate*10 + 2 );   /// TO CHECK !!!
      p2->AddSegment(  *s2 ); 
      nseg++; 
    }
  }

  SafeDelete(lst);

  if(cut) Log(2,"EdbDataPiece::GetCPData_new","select %d of %d segments by cut %s",nlst, nentr, cut->GetTitle() );
  else    Log(2,"EdbDataPiece::GetCPData_new","%d (of %d) segments accepted", nseg,nentr );

  return nseg;
}

///______________________________________________________________________________
int EdbDataPiece::CorrectAngles()
{
  int n=0;
  TTree *cptree=0;
  cptree=EdbDataPiece::InitCouplesTree(GetNameCP(),"READ");
  if( cptree )  n = CorrectAngles(cptree);
  Log(2,"CorrectAngles","in piece: %s  using %d basetracks",GetNameCP(),n);
  return n;
}

///______________________________________________________________________________
int EdbDataPiece::CorrectAngles(TTree *tree)
{
  EdbSegCouple    *cp = 0;
  EdbSegP         *s1 = 0;
  EdbSegP         *s2 = 0;
  EdbSegP         *s  = 0;
  tree->SetBranchAddress("cp"  , &cp  );
  tree->SetBranchAddress("s1." , &s1  );
  tree->SetBranchAddress("s2." , &s2  );
  tree->SetBranchAddress("s."  , &s   );

  EdbAffine2D *aff = new EdbAffine2D();

  int nentr = (int)(tree->GetEntries());
  TArrayF x(nentr);
  TArrayF y(nentr);
  TArrayF x1(nentr);
  TArrayF y1(nentr);
  TArrayF x2(nentr);
  TArrayF y2(nentr);

  int nseg = 0;
  Log(2,"EdbDataPiece::CorrectAngles","nentr = %d",nentr);
  for(int i=0; i<nentr; i++ ) {
    tree->GetEntry(i);

    if(cp->N1tot()>1)  continue;
    if(cp->N2tot()>1)  continue;
    if(cp->CHI2()>1.5) continue;

    x1[nseg] = s1->TX();
    y1[nseg] = s1->TY();
    x2[nseg] = s2->TX();
    y2[nseg] = s2->TY();
    x[nseg]  = (s2->X()-s1->X()) / (s2->Z()-s1->Z());
    y[nseg]  = (s2->Y()-s1->Y()) / (s2->Z()-s1->Z());
    nseg++;
  }

  aff->CalculateTurn( nseg,x1.fArray,y1.fArray,x.fArray,y.fArray );
  UpdateAffTPar(1,*aff);
  aff->CalculateTurn( nseg,x2.fArray,y2.fArray,x.fArray,y.fArray );
  UpdateAffTPar(2,*aff);

  delete aff;

  // Closing tree and file. A.C.

  TFile *f = tree->GetDirectory()->GetFile();
  SafeDelete(tree);
  if (f) SafeDelete(f);

  return nseg;
}

///______________________________________________________________________________
int EdbDataPiece::CheckCCD(int maxentr)
{
  if (eRun ) delete eRun;
  eRun =  new EdbRun( GetRunFile(0),"READ" );
  if(!eRun) { Log(1,"EdbDataPiece::CheckCCD","ERROR open file: %s",GetRunFile(0)); return -1; }
  EdbView    *view = eRun->GetView();
  EdbSegment *seg;

  int npeak=0;
  
  TMatrix matr(1000,1000);
  int ix,iy;

  int i,j;
  for(i=0; i<1000; i++) 
    for(j=0; j<1000; j++)
      matr[i][j]=0;

  int ncheck=0;
  int nentr = TMath::Min(maxentr,eRun->GetEntries());
  Log(2,"EdbDataPiece::CheckCCD","nentr=%d",nentr);
  for (i=0; i<nentr; i++) {
    view = eRun->GetEntry(i);
    int nseg=view->Nsegments();
    for (j=0; j<nseg; j++) {
      seg = view->GetSegment(j);
      if( seg->GetTx() >  .05 )   continue;
      if( seg->GetTx() < -.05 )   continue;
      if( seg->GetTy() >  .05 )   continue;
      if( seg->GetTy() < -.05 )   continue;
      ix = (Int_t)(seg->GetX0()+500.);
      iy = (Int_t)(seg->GetY0()+500.);
      (matr[ix][iy])++;
      ncheck++;
    }
  }

  Log(2,"EdbDataPiece::CheckCCD","ncheck=%d",ncheck);
  for(i=0; i<200; i++ ) {           //eliminate upto 200 CCD defects
    if(!RemoveCCDPeak(matr)) break;
    npeak++;
  }
  return npeak;
}

///______________________________________________________________________________
int EdbDataPiece::RemoveCCDPeak(TMatrix &matr)
{
  double mean=0;
  int filled=0;
  float max=0;
  int ix=0, iy=0;
  int nc = matr.GetNcols();
  int nr = matr.GetNrows();
  for(int i=0; i<nc; i++) 
    for(int j=0; j<nr; j++) {
      if(matr[i][j]<=0) continue;
      filled++;
      mean+=matr[i][j];
      if(max<matr[i][j]) {
	max = matr[i][j];
	ix = i;
	iy = j;
      }
    }
  mean/=filled;
  Log(2,"EdbDataPiece::RemoveCCDPeak","mean = %f \t max[%d,%d]=%d",mean,ix,iy,(int)max);

  float vmin[5],vmax[5];
  EdbSegmentCut cut;
  if( max > 10.*mean ) {
    matr[ix][iy]=0;
    vmin[0] = ix-500.;    vmax[0] = vmin[0]+1;
    vmin[1] = iy-500.;    vmax[1] = vmin[1]+1;
    vmin[2] = -.6;        vmax[2] = .6;
    vmin[3] = -.6;        vmax[3] = .6;
    vmin[4] =   0;        vmax[4] = 17;
    cut.SetXI(0);
    cut.SetMin(vmin);
    cut.SetMax(vmax);
    UpdateSegmentCut(cut);
    return 1;
  }
  return 0;
}

///______________________________________________________________________________
int EdbDataPiece::UpdateSegmentCut(EdbSegmentCut cut)
{
  const char *file=eFileNamePar.Data();

  FILE *fp=fopen(file,"a");
  if (!fp) {
    Log(1,"EdbDataPiece::UpdateSegmentCut","ERROR open file: %s", file);
    return -1;
  }

  char str[124];
  if(cut.XI()==0) {
    Log(2,"EdbDataPiece::UpdateSegmentCut","Update parameters file with XCUT: %s", file );
    sprintf(str,"XCUT \t %d \t %f %f %f %f %f %f %f %f %f %f\n",1,
	    cut.Min(0),cut.Max(0), 
	    cut.Min(1),cut.Max(1), 
	    cut.Min(2),cut.Max(2), 
	    cut.Min(3),cut.Max(3), 
	    cut.Min(4),cut.Max(4)
	    );
    fprintf(fp,"\n%s",str);
    sprintf(str,"XCUT \t %d \t %f %f %f %f %f %f %f %f %f %f\n",2,
	    cut.Min(0),cut.Max(0), 
	    cut.Min(1),cut.Max(1), 
	    cut.Min(2),cut.Max(2), 
	    cut.Min(3),cut.Max(3), 
	    cut.Min(4),cut.Max(4)
	    );
    fprintf(fp,"\n%s",str);
  } else  if(cut.XI()==1) {
    Log(2,"EdbDataPiece::UpdateSegmentCut","Update parameters file with ICUT: %s", file );
    sprintf(str,"ICUT \t %d \t %f %f %f %f %f %f %f %f %f %f\n",1,
	    cut.Min(0),cut.Max(0), 
	    cut.Min(1),cut.Max(1), 
	    cut.Min(2),cut.Max(2), 
	    cut.Min(3),cut.Max(3), 
	    cut.Min(4),cut.Max(4)
	    );
    fprintf(fp,"\n%s",str);
    sprintf(str,"XCUT \t %d \t %f %f %f %f %f %f %f %f %f %f\n",2,
	    cut.Min(0),cut.Max(0), 
	    cut.Min(1),cut.Max(1), 
	    cut.Min(2),cut.Max(2), 
	    cut.Min(3),cut.Max(3), 
	    cut.Min(4),cut.Max(4)
	    );
    fprintf(fp,"\n%s",str);
  }

  fclose(fp);
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::GetRawData(EdbPVRec *ali)
{
  //TODO: irun logic

  CloseRun();
  eRun =  new EdbRun( GetRunFile(0),"READ" );
  if(!eRun) { Log(1,"EdbDataPiece::GetRawData","ERROR open file: %s",GetRunFile(0)); return -1; }

  EdbPattern *pat1 = new EdbPattern( 0.,0., GetLayer(1)->Z() + GetLayer(0)->Z() );
  EdbPattern *pat2 = new EdbPattern( 0.,0., GetLayer(2)->Z() + GetLayer(0)->Z() );

  EdbViewHeader *head=0;
  EdbSegP        segP;
  EdbView       *view = eRun->GetView();
  int      nseg=0, nrej=0;
  int      nsegV=0;
  int      side=0;

  int nentr = eRun->GetEntries();
  for(int iv=0; iv<nentr; iv++ ) {
    head = eRun->GetEntryHeader(iv);
    if(      head->GetNframesTop()==0 ) side=2;
    else if( head->GetNframesBot()==0 ) side=1;

    if( !AcceptViewHeader(head) )  continue;


    if(eCLUST)       {
      view = eRun->GetEntry(iv,1,1,1);
      view->AttachClustersToSegments();
    }
    else             view = eRun->GetEntry(iv);

    nsegV = view->Nsegments();
    for(int j=0;j<nsegV;j++) {
      if(!TakeRawSegment(view,j,segP,side)) {
	nrej++;
	continue;
      }
      nseg++;
      segP.SetVid(iv,j);
      segP.SetAid(view->GetAreaID(),view->GetViewID(),side);

      if(side==1) pat1->AddSegment( segP );
      else if(side==2) pat2->AddSegment( segP );
    }

  }
  ali->AddPattern(pat1);
  ali->AddPattern(pat2);
  return nseg;
}

///______________________________________________________________________________
int EdbDataPiece::GetAreaData(EdbPVRec *ali, int aid, int side)
{
  TIndexCell *elist   = eAreas[side]->At(aid);
  if(!elist)  return 0;
  EdbPattern *pat = new EdbPattern( 0.,0., GetLayer(side)->Z() );
  pat->SetID(side);


  EdbSegP    segP;
  EdbView    *view = eRun->GetView();
  int   nseg=0, nrej=0;
  int   entry;
  int   nsegV;

  int niu=elist->N();
  for(int iu=0; iu<niu; iu++) {
    entry = elist->At(iu)->Value();
    if(eCLUST)       {
      view = eRun->GetEntry(entry,1,1,1);
      view->AttachClustersToSegments();
    }
    else             view = eRun->GetEntry(entry);

    nsegV = view->Nsegments();

    for(int j=0;j<nsegV;j++) {
      if(!TakeRawSegment(view,j,segP,side)) {
	nrej++;
	continue;
      }
      nseg++;
      segP.SetVid(entry,j);
      segP.SetAid(view->GetAreaID(),view->GetViewID(),side);
      pat->AddSegment( segP);
    }
  }

  Log(2,"EdbDataPiece::GetAreaData","Area: %d ( %d%%)  %d \t views: %d \t nseg: %d \t rejected: %d", 
	   aid,100*aid/eAreas[side]->N(1),side,niu,nseg, nrej );
  pat->SetSegmentsZ();
  ali->AddPattern(pat);
  return nseg;
}

///______________________________________________________________________________
int EdbDataPiece::MakeLinkListArea(int irun)
{
  if (eRun ) delete eRun;
  eRun =  new EdbRun( GetRunFile(irun),"READ" );
  if(!eRun) { Log(1," EdbDataPiece::MakeLinkListArea","ERROR open file: %s",GetRunFile(irun)); return -1; }

  for(int i=0; i<3; i++) {
    if(eAreas[i])    delete eAreas[i];
    eAreas[i]= new TIndexCell();
  }

  int nentr = eRun->GetEntries();
  Log(2,"EdbDataPiece::MakeLinkListArea","Make views entry map,  nentr = %d",nentr);

  Long_t v[2];   // areaID,entry
  EdbViewHeader *head=0;

  for(int iv=0; iv<nentr; iv++ ) {
    head = eRun->GetEntryHeader(iv);
    v[0]=head->GetAreaID();
    v[1]=iv;
    if(head->GetNframesTop()==0) {         // fill down views
      eAreas[2]->Add(2,v);
    }
    else if(head->GetNframesBot()==0) {    // fill up views
       eAreas[1]->Add(2,v);
    }
  }

  eAreas[1]->Sort();
  eAreas[2]->Sort();
  return eAreas[1]->N(1);
}

///______________________________________________________________________________
int EdbDataPiece::AcceptViewHeader(const EdbViewHeader *head)
{
  EdbSegP p;
  if(eAFID==0) {
    p.SetX( head->GetXview() );
    p.SetY( head->GetYview() );
  } else {
    p.SetX( 0. );
    p.SetY( 0. );
    p.Transform( head->GetAffine() );
  }
  p.SetX( p.X() - GetLayer(0)->X() );
  p.SetY( p.Y() - GetLayer(0)->Y() );
  p.Transform( GetLayer(0)->GetAffineXY() );
  if( p.X() < -GetLayer(0)->DX() )    return 0;
  if( p.X() >  GetLayer(0)->DX() )    return 0;
  if( p.Y() < -GetLayer(0)->DY() )    return 0;
  if( p.Y() >  GetLayer(0)->DY() )    return 0;
  return 1;
}

///______________________________________________________________________________
int EdbDataPiece::InitCouplesTree(const char *mode)
{
  if ((eCouplesTree=InitCouplesTree(eFileNameCP,mode))) return 1;
  return 0;
}

//______________________________________________________________________________
TTree *EdbDataPiece::InitCouplesTree(const char *file_name, const char *mode)
{
  static TString tree_name("couples");
  TTree *tree = NULL;

  // checking for the existing directory

  FileStat_t buf;

  if (gSystem->GetPathInfo(gSystem->DirName(file_name), buf)) {
    if (gEDBDEBUGLEVEL > 0) {
      cout << "ERROR! Directory " << gSystem->DirName(file_name) 
	   << " does not exist.\n";
      cout << "Please run scanning procedure first for the plate " 
	   << gSystem->BaseName(gSystem->DirName(file_name)) << endl;
    }
    return 0;

  }


  TFile *f = new TFile(file_name, mode);

  if (f->IsOpen()) tree = (TTree*)f->Get(tree_name);
  else return 0;

  if (!tree) {
    f->cd();
    tree = new TTree(tree_name, tree_name);
    tree->SetMaxTreeSize(15000000000LL);   //set 15 Gb file size limit)
    //tree->SetMaxVirtualSize( 512 * 1024 * 1024 ); // default is 64000000

    int pid1=0,pid2=0;
    float xv=0,yv=0;
    EdbSegCouple *cp=0;
    EdbSegP      *s1=0;
    EdbSegP      *s2=0;
    EdbSegP      *s=0;
      
    tree->Branch("pid1",&pid1,"pid1/I");
    tree->Branch("pid2",&pid2,"pid2/I");
    tree->Branch("xv",&xv,"xv/F");
    tree->Branch("yv",&yv,"yv/F");
    tree->Branch("cp","EdbSegCouple",&cp,32000,99);
    tree->Branch("s1.","EdbSegP",&s1,32000,99);
    tree->Branch("s2.","EdbSegP",&s2,32000,99);
    tree->Branch("s." ,"EdbSegP",&s,32000,99);
    tree->Write();
    // tree->SetAutoSave(2000000);
  }

  if(!tree) Log(1,"EdbDataPiece::InitCouplesTree","ERROR!!! InitCouplesTree: can't initialize tree at %s as %s",file_name,mode);

  return tree;
}


///______________________________________________________________________________
int EdbDataPiece::MakeLinkListCoord(int irun)
{
  if(eRun) delete eRun;
  eRun =  new EdbRun( GetRunFile(irun),"READ" );
  if(!eRun) { Log(1,"EdbDataPiece::MakeLinkListCoord","ERROR open file: %s",GetRunFile(0)); return -1; }

  
  for(int i=0; i<3; i++) {
    if(eAreas[i])    delete eAreas[i];
    eAreas[i]= new TIndexCell();
  }

  int nentr = eRun->GetEntries();
  Log(2,"EdbDataPiece::MakeLinkListCoord","Make views coordinate map,  nentr = %d",nentr);

  TIndexCell upc;
  TIndexCell downc;
  Long_t v[3];   // x,y,entry
  EdbViewHeader *head=0;

  Long_t xx=0, yy=0;
  float cx = 2000., cy = 2000.;  // 2x2 mm cells
  float dx = 400. , dy = 400.;   // 400 microns margins
  float xv=0,yv=0;
  int mx[9] = {0, 0, 0,-1,  1, -1,-1, 1, 1};
  int my[9] = {0,-1, 1, 0,  0, -1, 1,-1, 1};

  for(int iv=0; iv<nentr; iv++ ) {

    head = eRun->GetEntryHeader(iv);
    if( !AcceptViewHeader(head) )  continue;

    yv = head->GetXview();
    xv = head->GetYview();
    xx = (Long_t)(xv/cx);
    yy = (Long_t)(yv/cx);
    v[0] = xx;
    v[1] = yy;
    v[2] = iv;

    if(head->GetNframesBot()==0) {              // fill up views
      upc.Add(3,v);
    }
    else if(head->GetNframesTop()==0) {         // fill down views
      downc.Add(3,v);                           // add center
      int im;
      for( im=1; im<5; im++ ) {             // add sides margins
        v[0] = (Long_t)(( xv+dx*mx[im] ) / cx);
        v[1] = (Long_t)(( yv+dy*my[im] ) / cy);
        if( (v[0] != xx) || (v[1] != yy) )
            downc.Add(3,v);
      }
      for( im=5; im<9; im++ ) {             // add angles margins
        v[0] = (Long_t)(( xv+dx*mx[im] ) / cx);
        v[1] = (Long_t)(( yv+dy*my[im] ) / cy);
        if( (v[0] != xx) && (v[1] != yy) )
          //      if(!downc.Find(3,v))
            downc.Add(3,v);
      }

    }

  }

  upc.Sort();
  downc.Sort();

  TIndexCell *clx=0;
  TIndexCell *cly=0;

  int areac=0;
  int nix,niy,nie;
  nix=upc.N(1);
  for(int ix=0; ix<nix; ix++) {
    clx = upc.At(ix);
    xx = clx->Value();
    niy=clx->N(1);
    for(int iy=0; iy<niy; iy++) {
      cly  = clx->At(iy);
      yy = cly->Value();
      areac++;

      nie = cly->N(1);
      int ie;
      for(ie=0; ie<nie; ie++) {
        v[0]=areac;
        v[1]  = cly->At(ie)->Value();
        eAreas[1]->Add(2,v);
      }

      cly = downc.Find(xx)->Find(yy);
      if(!cly) continue;
      nie=cly->N(1);
      for(ie=0; ie<nie; ie++) {
        v[0] = areac;
        v[1] = cly->At(ie)->Value();
        eAreas[2]->Add(2,v);
      }
    }
  }

  eAreas[1]->Sort();
  eAreas[2]->Sort();

  return eAreas[1]->N(1);
}

///==============================================================================
EdbDataSet::EdbDataSet()
{
  Set0();
}

///------------------------------------------------------------------------------
EdbDataSet::EdbDataSet(const char *file)
{
  Set0();
  if(ReadDataSetDef(file)<0);
  Set0();
}

///______________________________________________________________________________
EdbDataSet::~EdbDataSet()
{
  if(ePieces.GetEntries()) ePieces.Delete();
}

///______________________________________________________________________________
void EdbDataSet::Set0()
{
  eInputList  = "runs.lst";
  eAnaDir     = "./";
  eParDir     = "./";
  eDBFileName = "pieces_DB.root";
  eDBFile     = 0;
  ePieces.SetOwner();
}

///______________________________________________________________________________
int EdbDataSet::ReadDataSetDef(const char *file)
{
  char            buf[256];
  char            key[256];
  char            name[256];

  FILE *fp=fopen(file,"r");
  if (fp==NULL)   {
    Log(1,"EdbDataSet::ReadDataSetDef","ERROR open file: %s", file);
    return(-1);
  }else
    Log(2,"EdbDataSet::ReadDataSetDef", "Read Data Set Definitions from: %s", file );

  while( fgets(buf,256,fp)!=NULL ) {
    for( int i=0; i<(int)strlen(buf); i++ ) 
      if( buf[i]=='#' )  {
	buf[i]='\0';                       // cut out comments starting from #
	break;
      }

   if( sscanf(buf,"%s",key)!=1 )                             continue;

   if      ( !strcmp(key,"OUTPUT_DATA_DIR")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eAnaDir = name;
      }
   else if ( !strcmp(key,"INPUT_RUNS_LIST")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eInputList=name;
      }
   else if ( !strcmp(key,"PARAMETERS_DIR")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eParDir=name;
      }
   else if ( !strcmp(key,"DBFILNAME")   )
      {
	sscanf(buf+strlen(key),"%s",name);
	eDBFileName=name;
      }
  }
  fclose(fp);
  return GetRunList(eInputList.Data());
}

///______________________________________________________________________________
void EdbDataSet::PrintRunList()
{
  for(int i=0; i<ePieces.GetEntriesFast(); i++){
    ((EdbDataPiece*)ePieces.At(i))->Print();
  }
}

///______________________________________________________________________________
void EdbDataSet::Print()
{
  EdbDataPiece *p=0;
  printf("EdbDataSet with %d pieces:",N());
  for(int i=0; i<N(); i++){
    p = ((EdbDataPiece*)ePieces.At(i));
    printf("%s %s\n",p->GetName(), p->GetRunFile(0));
  }
}


///______________________________________________________________________________
void EdbDataSet::WriteRunList()
{
  if(!eDBFile) eDBFile = new TFile(eDBFileName.Data(),"UPDATE");
  //else eDBFile->Cd();

  for(int i=0; i<ePieces.GetEntriesFast(); i++){
    ((EdbDataPiece*)ePieces.At(i))->MakeName();
  }

  ePieces.Write("pieces",1);
}

///______________________________________________________________________________
EdbDataPiece *EdbDataSet::FindPiece(const char *name)
{
  const char *nn;
  EdbDataPiece *piece=0;
  for(int i=0; i<ePieces.GetEntriesFast(); i++){
    piece = (EdbDataPiece*)ePieces.At(i);
    nn = piece->GetName();
    if(!strcmp(nn,name)) return piece;
  }
  return 0;
}

///______________________________________________________________________________
int EdbDataSet::GetRunList(const char *file)
{
  char            buf[256];
  char            filename[256];
  int             nrun=0;

  FILE *fp=fopen(file,"r");
  if (!fp)   {
    Log(1,"EdbDataSet::GetRunList","ERROR open file: %s", file);
    return(-1);
  }
  else Log(2,"EdbDataSet::GetRunList", "Read runs list from file: %s", file );

  EdbDataPiece *piece=0;
  EdbDataPiece *pp=0;
  int plateID,pieceID,flag;

  int ntok=0;
  while( fgets(buf,256,fp)!=NULL ) {
    ntok = sscanf(buf,"%d %d %s %d",&plateID,&pieceID,filename,&flag);
    if(ntok!=4) break;
    if(flag<=0) continue;
    nrun++;
    piece = new EdbDataPiece(plateID,pieceID,filename,flag);
    piece->MakeName();
    if( (pp=FindPiece(piece->GetName())) ) {
      pp->AddRunFile( filename );
      delete piece;
    } else {
      piece->MakeNameCP(GetAnaDir());
      piece->MakeNamePar(GetParDir());
      if(piece->TakePiecePar()>=0)
	ePieces.Add(piece);
      else {
	Log(1,"EdbDataSet::GetRunList","Missing par file for piece!!!");
	if(fp) fclose(fp);
	return -1;
      }
    }
  }
  fclose(fp);

  return nrun;
}

///==============================================================================
EdbDataProc::EdbDataProc()
{
  eDataSet=0;
  eNoUpdate=1;
  ePVR     = new EdbPVRec();
}

///==============================================================================
EdbDataProc::EdbDataProc(const char *file)
{
  eDataSet = new EdbDataSet(file);
  ePVR     = 0;
  eNoUpdate=0;
}

///------------------------------------------------------------------------------
EdbDataProc::EdbDataProc(int npl, TArrayI &ids, TArrayF &zs)
{
	// Special constructor (requested by Nicolai) set everything as default as possible

  eDataSet  = new EdbDataSet();
  ePVR      = 0;
  eNoUpdate = 0;

  if( npl>ids.GetSize() )  return;
  if( npl>zs.GetSize() )   return;

  ePVR = new EdbPVRec();
  EdbPattern *pat =0;
  for(int i=0; i<npl; i++) {
    pat = new EdbPattern( 0,0, zs.At(i) );
    pat->SetPID(ids.At(i));
    ePVR->AddPattern(pat);
  }
  ePVR->SetPatternsID();

  ePVR->SetScanCond( new EdbScanCond() );
  ePVR->SetCouplesAll();
}

///------------------------------------------------------------------------------
EdbDataProc::~EdbDataProc()
{
  if(eDataSet) delete eDataSet;
  if(ePVR)     delete ePVR;
}

///------------------------------------------------------------------------------
int EdbDataProc::CheckCCD()
{
  if(!eDataSet) return 0;
  EdbDataPiece *piece;
  int ndef=0;
  int np=eDataSet->N();
  for(int i=0; i<np; i++) {
    piece = eDataSet->GetPiece(i);
    ndef = piece->CheckCCD();
    if(ndef<0) { Log(2,"EdbDataProc::CheckCCD","skip piece"); continue; }
    Log(1,"EdbDataProc::CheckCCD","piece %s: eliminated defects: %d",piece->GetName(),ndef); 
    piece->WriteCuts();
    piece->CloseRun();
  }
  return np;
}

//----------------------------------------------------------------------------
int EdbDataProc::Link()
{
  if (!eDataSet) return 0;
  EdbDataPiece *piece;
  Int_t np = eDataSet->N();
  for (Int_t i = 0; i < np; i++) {
    piece = eDataSet->GetPiece(i);
    if (piece->TakePiecePar() >= 0) {
      if (gEDBDEBUGLEVEL > 2) piece->Print();
      piece->WriteCuts();
      if (piece->Flag() == 1) Link(*piece);
    }
    piece->CloseRun();
  }
  return np;
}

//______________________________________________________________________________
int EdbDataProc::Link(EdbDataPiece &piece)
{
  EdbPVRec  *ali;

  //const char *file_name=piece.GetNameCP();
  //TTree *cptree=EdbDataPiece::InitCouplesTree(file_name,"RECREATE");
  //if(!cptree) return 0;

  if (!piece.InitCouplesTree("RECREATE")) return 0;

  EdbScanCond *cond = piece.GetCond(1);
  int    ntot=0, nareas=0;
  float  shr1=1,shr2=1;
  double shrtot1=0., shrtot2=0.;
  int    nshr=0,nshrtot=0;

  for( int irun=0; irun<piece.Nruns(); irun++ ) {
    nareas = piece.MakeLinkListCoord(irun);
    //nareas = piece.MakeLinkListArea(irun);
    Log(2,"EdbDataProc::Link","%d areas mapped",nareas);
    if (nareas <= 0) continue;

    for (Int_t i = 0; i < nareas; i++ ) {

      ali      = new EdbPVRec();
      ali->SetScanCond( cond );

      piece.GetAreaData(ali,i,1);
      piece.GetAreaData(ali,i,2);

      //ali->SetSegmentsErrors();

      ali->SetCouplesAll();
      ali->SetChi2Max(cond->Chi2PMax());

      for(int ic=0; ic<ali->Ncouples(); ic++) 
	ali->GetCouple(ic)->SetCHI2mode(cond->Chi2Mode());

      ali->Link();

      if( ShrinkCorr() ) {
	shr1 = 1;
	shr2 = 1;
	nshr = CheckShrinkage( ali,0, shr1, shr2 );
	if(nshr) {
	  nshrtot += nshr;
	  shrtot1 += nshr*shr1;
	  shrtot2 += nshr*shr2;
	}
      }

      FillCouplesTree( piece.eCouplesTree, ali,piece.GetOUTPUT());  //!!! 0

      delete ali;
    }
    ntot+=nareas;
    if(nshrtot>3) {
      shrtot1 = shrtot1/nshrtot;
      shrtot2 = shrtot2/nshrtot;
      if(nshrtot<20) Log(1,"EdbDataProc::Link","WARNING: unreliable shrinkage correction - low statistics");
      Log(2,"EdbDataProc::Link","Shrinkage correction(%d): %f %f", nshrtot, (float)shrtot1,(float)shrtot2);
      piece.CorrectShrinkage( 1, (float)shrtot1 );
      piece.CorrectShrinkage( 2, (float)shrtot2 );
      if(!NoUpdate())   piece.UpdateShrPar(1);
      if(!NoUpdate())   piece.UpdateShrPar(2);
    }
  }
  int ncp = piece.eCouplesTree->GetEntries();
  piece.CloseCPData();
  return ncp;
}

///______________________________________________________________________________
void EdbDataProc::FillCouplesTree( TTree *tree, EdbPVRec *al, int fillraw )
{
  tree->GetDirectory()->cd();

  Log(2,"EdbDataProc::FillCouplesTree","fill couples tree...");

  EdbPatCouple *patc=0;
  float xv = al->X();
  float yv = al->Y();

  int pid1,pid2;
  EdbSegCouple *cp=0;
  EdbSegP      *s1=0;
  EdbSegP      *s2=0;
  EdbSegP      *s=0;

  tree->SetBranchAddress("pid1",&pid1);
  tree->SetBranchAddress("pid2",&pid2);
  tree->SetBranchAddress("xv"  ,&xv);
  tree->SetBranchAddress("yv"  ,&yv);
  tree->SetBranchAddress("cp"  ,&cp);
  tree->SetBranchAddress("s1." ,&s1);
  tree->SetBranchAddress("s2." ,&s2);
  tree->SetBranchAddress("s."  ,&s );

  if(fillraw) {
    // **** fill tree with raw segments ****
    EdbPattern *pat=0;
    int nic;
    int nip=al->Npatterns();
    for( int ip=0; ip<nip; ip++ ) {
      pat  = al->GetPattern(ip);
      pid1 = pat->ID();
      pid2 = -1;
      nic=pat->N();
      for( int ic=0; ic<nic; ic++ ) {
        s1 = pat->GetSegment(ic);
        tree->Fill();
      }
    }
  }

  // **** fill tree with found couples ****

  s = new EdbSegP();

  int nip=al->Ncouples();
  for( int ip=0; ip<nip; ip++ ) {
    patc = al->GetCouple(ip);
    pid1 = patc->Pat1()->ID();
    pid2 = patc->Pat2()->ID();

    int nic=patc->Ncouples();
    for( int ic=0; ic<nic; ic++ ) {
      cp = patc->GetSegCouple(ic);
      s1 = patc->Pat1()->GetSegment(cp->ID1());
      s2 = patc->Pat2()->GetSegment(cp->ID2());
      s = cp->eS;
      tree->SetBranchAddress("s."  ,&s );

      s->SetID(tree->GetEntries());             // basetrack id will be the tree entry number
      EdbTraceBack::SetBaseTrackVid( *s, 0, 0, tree->GetEntries() );   //TODO: plate, piece if available
      tree->Fill();
    }
  }

}

///______________________________________________________________________________
void EdbDataProc::CloseCouplesTree(TTree *tree)
{
  tree->AutoSave();
  TFile *f=0;
  f = tree->GetCurrentFile();
  if (f) SafeDelete(f);
  tree=0;
}

///______________________________________________________________________________
int EdbDataProc::CheckShrinkage(EdbPVRec *ali, int couple, float &shr1, float &shr2)
{
  EdbPatCouple  *patc = ali->GetCouple(couple);
  if(!patc) return 0;

  EdbSegCouple *sc=0;
  EdbSegP *s1=0, *s2=0;
  double dz,tx,ty,t,t1,t2,sumt1=0,sumt2=0;
  int    nsum=0;

  int  nc=patc->Ncouples();
  for( int ic=0; ic<nc; ic++ ) {
    sc = patc->GetSegCouple(ic);

    if(sc->CHI2()>1.5)  continue;

    s1 = patc->Pat1()->GetSegment(sc->ID1());
    s2 = patc->Pat2()->GetSegment(sc->ID2());

    dz = s2->Z() - s1->Z();
    tx = (s2->X() - s1->X())/dz;
    ty = (s2->Y() - s1->Y())/dz;
 
    t  = TMath::Sqrt( tx*tx + ty*ty );
    if(t<.1) continue;
    if(t>.45) continue;
    t1 = TMath::Sqrt( s1->TX()*s1->TX() + s1->TY()*s1->TY() );
    t2 = TMath::Sqrt( s2->TX()*s2->TX() + s2->TY()*s2->TY() );

    nsum++;
    sumt1 += t1/t;
    sumt2 += t2/t;
  }

  if(nsum<1)  return 0;

  shr1 = sumt1/nsum;
  shr2 = sumt2/nsum;

  return nsum;
}

///______________________________________________________________________________
void EdbDataProc::CorrectAngles()
{
  EdbDataPiece *piece;
  int npieces = eDataSet->N();
  Log(2,"EdbDataProc::CorrectAngles","npieces = %d",npieces);
  if(!npieces) return;

  for(int i=0; i<npieces; i++ ) {
    piece = eDataSet->GetPiece(i);
    if(piece) piece->CorrectAngles();
  }
}

///______________________________________________________________________________
int EdbDataProc::InitVolumeRaw(EdbPVRec    *ali)
{
  EdbScanCond *cond = eDataSet->GetPiece(0)->GetCond(1);
  ali->SetScanCond( cond );
  
  EdbDataPiece *piece;
  int npieces = eDataSet->N();
  Log(2,"EdbDataProc::InitVolumeRaw","npieces = %d",npieces);
  if(!npieces) return 0;

  for(int i=0; i<npieces; i++ ) {
    piece = eDataSet->GetPiece(i);
    piece->GetRawData(ali);
  }

  float x0 = eDataSet->GetPiece(npieces-1)->GetLayer(0)->X();
  float y0 = eDataSet->GetPiece(npieces-1)->GetLayer(0)->Y();
  ali->Centralize(x0,y0);

  for(int j=0; j<npieces; j++ ) {
    for(int ip=0; ip<2; ip++ ) {
      int i = 2*j+ip;
      ali->GetPattern(i)->SetSegmentsZ();
      ali->GetPattern(i)->Transform(    eDataSet->GetPiece(j)->GetLayer(0)->GetAffineXY()   );
      ali->GetPattern(i)->TransformA(   eDataSet->GetPiece(j)->GetLayer(0)->GetAffineTXTY() );
      ali->GetPattern(i)->TransformShr( eDataSet->GetPiece(j)->GetLayer(0)->Shr()  );
    }
  }
  ali->SetPatternsID();
  ali->SetSegmentsErrors();

  ali->SetCouplesAll();
  ali->SetChi2Max(cond->Chi2PMax());
  ali->SetOffsetsMax(cond->OffX(),cond->OffY());
  return npieces;
}

///______________________________________________________________________________
EdbPVRec *EdbDataProc::ExtractDataVolume( EdbSegP &seg, int plmin, int plmax,
					  float acc[4],
					  int datatype   )
{
  // datatype: 0   - base track
  //           1   - up
  //           2   - down
  //           100 - raw?

  if(!ePVR)               return 0;
  int npat = ePVR->Npatterns();

  EdbPVRec *ali = new EdbPVRec();

  Log(2,"EdbDataProc::ExtractDataVolume","Select Data %d in EdbPVRec with %d patterns:",
	 datatype, npat);

  float dz, min[5],  max[5];
  EdbPattern *pat=0;

  for(int i=0; i<npat; i++) {

    pat = ePVR->GetPattern(i);
    if(!pat)                   continue;

    dz = pat->Z() - seg.Z();

    min[0] = seg.X() + dz*seg.TX() - acc[0];
    min[1] = seg.Y() + dz*seg.TY() - acc[1];
    max[0] = seg.X() + dz*seg.TX() + acc[0];
    max[1] = seg.Y() + dz*seg.TY() + acc[1];

    min[2] = -.6;    max[2] = .6;
    min[3] = -.6;    max[3] = .6;
    min[4] =  0.;    max[4] = 100.;

    ali->AddPattern( pat->ExtractSubPattern(min,max) );
  }
  
  return ali;
}

///______________________________________________________________________________
EdbPVRec *EdbDataProc::ExtractDataVolumeF( EdbTrackP &tr, float binx, float bint,
					   int datatype   )
{
  if(!ePVR)               return 0;
  int npat = ePVR->Npatterns();
  EdbPVRec *ali = new EdbPVRec();
  Log(2,"EdbDataProc::ExtractDataVolumeF","Select Data %d in EdbPVRec with %d patterns:",
	 datatype, npat);

  EdbSegP ss; // the "selector" segment 
  ss.SetCOV( tr.GetSegment(0)->COV() );

  float dz  = (tr.GetSegment(tr.N()-1)->Z() - tr.GetSegment(0)->Z());
  float tx  = (tr.GetSegment(tr.N()-1)->X() - tr.GetSegment(0)->X())/dz;
  float ty  = (tr.GetSegment(tr.N()-1)->Y() - tr.GetSegment(0)->Y())/dz;

  ss.SetTX(tx);
  ss.SetTY(ty);
  ss.SetX(tr.X());
  ss.SetY(tr.Y());
  ss.SetZ(tr.Z());

  if(gEDBDEBUGLEVEL>2) ss.Print();

  EdbPattern *pat  = 0;
  EdbPattern *spat = 0;
  int nseg =0;

  TObjArray arr;

  for(int i=0; i<npat; i++) {
    pat = ePVR->GetPattern(i);
    if(!pat)                   continue;
    ss.PropagateTo(pat->Z());
    nseg += pat->FindCompliments(ss,arr,binx,bint);
  }

  spat = new EdbPattern(0,0,0);
  for(int i=0; i<arr.GetEntriesFast(); i++)
    spat->AddSegment( *((EdbSegP*)(arr.At(i))) );

  ali->AddPattern( spat );

  Log(2,"EdbDataProc::ExtractDataVolumeF","%d segments are selected",nseg);
  return ali;
}

///______________________________________________________________________________
EdbPVRec *EdbDataProc::ExtractDataVolume( EdbTrackP &tr, float binx, float bint,
					  int datatype   )
{
  // datatype: 0   - base track
  //           1   - up
  //           2   - down
  //           100 - raw?

  if(!ePVR)               return 0;
  int npat = ePVR->Npatterns();

  EdbPVRec *ali = new EdbPVRec();

  Log(2,"EdbDataProc::ExtractDataVolume","Select Data %d in EdbPVRec with %d patterns",
	 datatype, npat);

  float min[5],  max[5];

  float dx  = binx*TMath::Sqrt(tr.SX()*tr.N());
  float dy  = binx*TMath::Sqrt(tr.SY()*tr.N());
  float dtx = bint*TMath::Sqrt(tr.STX()*tr.N());
  float dty = bint*TMath::Sqrt(tr.STY()*tr.N());

  float dz  = (tr.GetSegment(tr.N()-1)->Z() - tr.GetSegment(0)->Z());
  float tx  = (tr.GetSegment(tr.N()-1)->X() - tr.GetSegment(0)->X())/dz;
  float ty  = (tr.GetSegment(tr.N()-1)->Y() - tr.GetSegment(0)->Y())/dz;

  Log(2,"EdbDataProc::ExtractDataVolume","select segments with acceptance: %f %f [microns]   %f %f [mrad]",
	 dx,dy,dtx*1000,dty*1000);

  EdbPattern *pat  = 0;
  EdbPattern *spat = 0;
  int nseg =0;

  for(int i=0; i<npat; i++) {

    pat = ePVR->GetPattern(i);
    if(!pat)                   continue;

    dz = pat->Z() - tr.Z() + 107.;        //TODO

    min[0] = tr.X() + dz*tx - dx;
    max[0] = tr.X() + dz*tx + dx;
    min[1] = tr.Y() + dz*ty - dy;
    max[1] = tr.Y() + dz*ty + dy;
    min[2] = tx - dtx;
    max[2] = tx + dtx;
    min[3] = ty - dty;
    max[3] = ty + dty;
    min[4] =  0.;    
    max[4] = 100.;

    spat = pat->ExtractSubPattern(min,max);
    nseg += spat->N();
    ali->AddPattern( spat );
  }
  Log(2,"EdbDataProc::ExtractDataVolume","%d segments are selected",nseg);
  return ali;
}

///______________________________________________________________________________
int EdbDataProc::InitVolume(int datatype, const char *rcut)
{
  // datatype: 0    - couples basetrack data
  //           10   - couples full data
  //           100  - tracks only
  //           1000 - tracks and basetracks

  if(!ePVR) ePVR = new EdbPVRec();
  if(datatype==100)  return InitVolumeTracks(ePVR, rcut);

  if(datatype==1000) {
    InitVolumeTracks(ePVR, rcut);
    TIndex2 *trseg=MakeTracksSegmentsList( *ePVR );
    delete ePVR;
    ePVR = new EdbPVRec();
    int n = InitVolume(ePVR, datatype, trseg);
    n += InitVolumeTracks(ePVR, rcut);
    return n;
  }

  return InitVolume(ePVR, datatype);
}


///______________________________________________________________________________
int EdbDataProc::InitVolumeTracks(EdbPVRec    *ali, const char *rcut)
{
  EdbScanCond *cond=0;
  if(eDataSet) {
    cond= eDataSet->GetPiece(0)->GetCond(0);
    ali->SetScanCond( cond );
  
    EdbDataPiece *piece;
    int npieces = eDataSet->N();
    Log(2,"EdbDataProc::InitVolumeTracks","npieces = %d",npieces);
    if(!npieces) return 0;

    if( ali->Npatterns() == 0) {                    // case of segments already readed
      EdbAffine2D *a=0;
      EdbPattern  *pat=0;
      for(int i=0; i<npieces; i++ ) {
	piece = eDataSet->GetPiece(i);
	pat = new EdbPattern( 0.,0., piece->GetLayer(0)->Z(),100 );
	a = piece->GetLayer(0)->GetAffineXY();
	pat->SetKeep( a->A11(),a->A12(),a->A21(),a->A22(), a->B1(),a->B2() );
	pat->SetPID(i);
	ali->AddPattern( pat );
      }
      ali->SetPatternsID();
    }
  }
  else {
    cond = new EdbScanCond();
    ali->SetScanCond( cond );
  }

  ReadTracksTree( *ali, "linked_tracks.root", rcut );

  ali->SetSegmentsTracks();  // reset the id of tracks????? to check that it is really necessary!!
  ali->SetSegmentsErrors();
  ali->SetCouplesAll();
  ali->SetChi2Max(cond->Chi2PMax());
  ali->SetOffsetsMax(cond->OffX(),cond->OffY());
  return ali->Npatterns();
}

///______________________________________________________________________________
void EdbDataProc::FineAlignmentTracks()
{
  InitVolume(100);
  EdbPVRec *ali = PVR();

  EdbAffine2D aff,afft;

  int fctr=0,fctr0=10000, fcMin=50;

  for( int i=0; i<ali->Npatterns(); i++ ) {
    fctr = ali->FineCorrF(i,aff,afft);
    if(fctr<fctr0) fctr0=fctr;
    Log(2,"EdbDataProc::FineAlignmentTracks","fctr = %d",fctr);
  }

  if(fctr0>fcMin) {
    for( int i=0; i<ali->Npatterns(); i++) {
      ali->GetPattern(i)->GetKeep(aff);
      if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateAffPar(0,aff);

      //if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateAffTPar(0,afft);
    }
  }
}

///______________________________________________________________________________
int EdbDataProc::InitVolume(EdbPVRec    *ali, int datatype, TIndex2 *trseg)
{
  if (!eDataSet) Log(2,"EdbDataProc::InitVolume","Warning. eDataSet is NULL"); 
  
  EdbScanCond *cond = eDataSet->GetPiece(0)->GetCond(0);
  ali->SetScanCond( cond );
  
  EdbDataPiece *piece;
  int npieces = eDataSet->N();
  Log(2,"EdbDataProc::InitVolume","npieces = %d",npieces);
  if(!npieces) return 0;

  //TTree *cptree=0;
  int i;

  EdbPattern *pat=0;
  EdbPattern *p1=0;
  EdbPattern *p2=0;

  for(i=0; i<npieces; i++ ) {
    piece = eDataSet->GetPiece(i);
    if(!piece->InitCouplesTree("READ"))  Log(1," EdbDataProc::InitVolume","no tree %d",i);

    pat = new EdbPattern( 0.,0., piece->GetLayer(0)->Z(),100 );
    pat->SetPID(i);
    //pat->SetSegmentsPID();  //TO CHECK!!
    if(datatype==10) {
      p1 = new EdbPattern( 0.,0., piece->GetLayer(0)->Z() + piece->GetLayer(1)->Z() );
      p2 = new EdbPattern( 0.,0., piece->GetLayer(0)->Z() + piece->GetLayer(2)->Z() );
      p1->SetPID(i);
      p2->SetPID(i);
    }

    piece->GetCPData_new( pat,p1,p2,trseg );

    ali->AddPattern( pat );
    if(datatype==10) {
      ali->AddPattern( p1 );
      ali->AddPattern( p2 );
    }
    piece->CloseCPData();
  }

  float x0 = eDataSet->GetPiece(npieces-1)->GetLayer(0)->X();
  float y0 = eDataSet->GetPiece(npieces-1)->GetLayer(0)->Y();
  //ali->Centralize();
  ali->Centralize(x0,y0);

  int npat = ali->Npatterns();
  for(i=0; i<npat; i++ ) {
    pat = ali->GetPattern(i);
    pat->SetSegmentsZ();   /// TO CHECK!!!
    pat->Transform(    eDataSet->GetPiece(pat->PID())->GetLayer(0)->GetAffineXY()   );
    pat->TransformA(   eDataSet->GetPiece(pat->PID())->GetLayer(0)->GetAffineTXTY() );
    pat->TransformShr( eDataSet->GetPiece(pat->PID())->GetLayer(0)->Shr()  );
  }
  ali->SetPatternsID();

  ali->SetSegmentsErrors();

  ali->SetCouplesAll();
  ali->SetChi2Max(cond->Chi2PMax());
  ali->SetOffsetsMax(cond->OffX(),cond->OffY());
  return npieces;
}

///______________________________________________________________________________
void EdbDataProc::Align(int doAlign)
{
  EdbPVRec    *ali  = new EdbPVRec();
  InitVolume(ali);

  ali->Align(doAlign);

  if(gEDBDEBUGLEVEL>1) ali->PrintAff();
  EdbAffine2D  aff;
  for(int i=0; i<ali->Npatterns(); i++) {
    ali->GetPattern(i)->GetKeep(aff);
    if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateAffPar(0,aff);
  }
}

///______________________________________________________________________________
void EdbDataProc::AjustZ(int doZ)
{
  EdbPVRec    *ali  = new EdbPVRec();
  InitVolume(ali);
  ali->Link();
  ali->FillTracksCell();
  ali->SelectLongTracks(doZ);
  ali->MakeSummaryTracks();
  ali->FineCorrZnew();
  for(int i=0; i<ali->Npatterns(); i++) 
    if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateZPar(0,ali->GetPattern(i)->Z());
}

///______________________________________________________________________________
int EdbDataProc::LinkTracksWithFlag(EdbPVRec *ali, float p, float probmin, int nsegmin, int maxgap, int flag, float mass)
{
  // p       - momentum expected for tracks
  // probmin - minimum probability to accept tracks
  // nsegmin - min number of segments/track
  // maxgap  - max gap permitted for propagation
  // flag    - to be assigned for the tracks found in this pass

  int ntr0 = ali->Ntracks();
  ali->Link();
  ali->FillTracksCell();
  ali->MakeTracks(nsegmin,flag);

  int       noProp=0;
  if(p<0)   noProp=1;
  //if(p<0.01) p=4.;  //TODO: this protection should be out of this function?

  int ntr = ali->Ntracks();
  float X0 =  ali->GetScanCond()->RadX0();
  EdbTrackP *tr=0;
  for(int itr=ntr0; itr<ntr; itr++) {
    tr = ali->GetTrack(itr);
    tr->ClearF();
    if(p<0) tr->SetP(1.);
    else    tr->SetP(p);
    tr->SetM(mass);
    tr->FitTrackKFS(false,X0);
  }

  if(!noProp)
    if( ali->Npatterns()>2 ) {
      ali->FillCell(50,50,0.015,0.015);
      for(int i=0; i<10; i++) 
	if( ali->PropagateTracks(ali->Npatterns()-1,2, probmin, maxgap ) <1) break;
    }

  /*
  ntr = ali->Ntracks();
  for(int i=0; i<ntr; i++) {
    tr = ali->GetTrack(i);
    if(tr->Flag()<0) continue;
    if(tr->N()<5)  tr->SetP(0.5); 
    else           tr->SetP(tr->P_MS());
  }

  for(int i=0; i<10; i++) 
    if( ali->PropagateTracks(ali->Npatterns()-1,2, probmin, maxgap ) <1) break;
  */

  return ntr-ntr0;
}

///______________________________________________________________________________
void EdbDataProc::LinkTracks( int alg, float p )
{
  EdbPVRec    *ali  = new EdbPVRec();
  InitVolume(ali);

  Log(2,"EdbDataProc::LinkTracks","tracking:  alg = %d \t p= %f", alg, p);

  ali->SetCouplesPeriodic(0,1);
  int ntr=0;

  ntr = LinkTracksWithFlag( ali, p, 0.05, 2, 3, 0 );
  //ntr = LinkTracksWithFlag( ali, p/4., 0.01, 2, 3, 1 );
  //if(merge>0) ali->MergeTracks(merge);

  if(alg==-1) { // this is a test option produce couples tree from the first 2 linked patterns
    TTree *cptree=EdbDataPiece::InitCouplesTree("linked_couples.root","RECREATE");
    FillCouplesTree(cptree, ali,0);
    CloseCouplesTree(cptree);
  }

  float mass=0.139;                    //TODO!
  ali->FitTracks( p, mass );         // is important to call it before MakeTracksTree!

  MakeTracksTree(ali);
  
  // After the EdbPVRec is created and filled, set it also as the EdbDataProc own object.
  SetPVR(ali);
  if (gEDBDEBUGLEVEL>2) {
    ali->Print();
  }
  Log(2,"EdbDataProc::LinkTracks","EdbDataProc::LinkTracks...Done.");
  return;
}

///______________________________________________________________________________
void EdbDataProc::LinkTracksC( int alg, float p )
{
  // special tracking procedure for carbonium dataset

  EdbPVRec    *ali  = new EdbPVRec();
  InitVolume(ali);

  Log(2,"EdbDataProc::LinkTracksC","carbonium tracking:  alg = %d \t p= %f", alg, p);

  ali->SetCouplesPeriodic(0,1);
  LinkTracksWithFlag( ali, p, 0.01, 2, 3, 0 );

  ali->SetCouplesPeriodic(0,3);
  LinkTracksWithFlag( ali, p, 0.01, 2, 6, 1 );

  ali->SetCouplesPeriodic(1,3);
  LinkTracksWithFlag( ali, p, 0.01, 2, 6, 2 );

  ali->SetCouplesPeriodic(2,3);
  LinkTracksWithFlag( ali, p, 0.01, 2, 6, 3 );

  float mass=0.139;                    //TODO!
  ali->FitTracks( p, mass );         // is important to call it before MakeTracksTree!

  MakeTracksTree(ali);
}

///______________________________________________________________________________
void EdbDataProc::LinkRawTracks( int alg )
{
  EdbPVRec    *ali  = new EdbPVRec();
  InitVolumeRaw(ali);
  ali->Link();
  Log(2,"EdbDataProc::LinkRawTracks","link ok");

//    if(alg==1) {
//      ali->MakeHoles();
//      ali->Link();
//      printf("link ok\n");
//    }

  ali->FillTracksCell();

  TTree *cptree=EdbDataPiece::InitCouplesTree("linked_couples.root","RECREATE");
  FillCouplesTree(cptree, ali,0);
  CloseCouplesTree(cptree);

  MakeTracksTree(ali);
}

///______________________________________________________________________________
void EdbDataProc::FineAlignment(int doFine)
{
  EdbPVRec    ali;
  InitVolume(&ali);
  ali.Link();
  Log(2,"EdbDataProc::FineAlignment","link ok");
  ali.FillTracksCell();

  EdbAffine2D aff;

  int fctr=0;
  int fcMin = 49;

  int i;
  ali.SelectLongTracks(ali.Npatterns());
  ali.MakeSummaryTracks();
  for( i=0; i<ali.Npatterns(); i++ ) {
    fctr = ali.FineCorrXY(i,aff,doFine);
  }

  ali.SelectLongTracks(ali.Npatterns());
  ali.MakeSummaryTracks();
  for( i=0; i<ali.Npatterns(); i++ ) {
    fctr = ali.FineCorrXY(i,aff,doFine);
  }

  if(fctr>fcMin) {
    for( i=0; i<ali.Npatterns(); i++) {
      ali.GetPattern(i)->GetKeep(aff);
      if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateAffPar(0,aff);
    }
  }

  ali.SelectLongTracks(ali.Npatterns());
  ali.MakeSummaryTracks();
  float dz=0;
  float z = ali.GetPattern(ali.Npatterns()-1)->Z();
  for( i=ali.Npatterns()-2; i>=0; i-- ) {
    fctr = ali.FineCorrZ(i,dz);
    if(fctr<=fcMin) break;
    z -= dz;
    Log(2,"EdbDataProc::FineAlignment","dz = %f  z = %f",dz,z);
    if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateZPar(0,z);
  }

  ali.SelectLongTracks(ali.Npatterns());
  ali.MakeSummaryTracks();
  for( i=0; i<ali.Npatterns(); i++ ) {
    fctr = ali.FineCorrTXTY(i,aff);
    if(fctr<=fcMin) break;
    if(gEDBDEBUGLEVEL>2) aff.Print();
    if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateAffTPar(0,aff);
  }

  ali.SelectLongTracks(ali.Npatterns());
  ali.MakeSummaryTracks();
  dz=0;
  z = ali.GetPattern(ali.Npatterns()-1)->Z();
  for( i=ali.Npatterns()-2; i>=0; i-- ) {
    fctr = ali.FineCorrZ(i,dz);
    if(fctr<=fcMin) break;
    z -= dz;
    Log(2,"EdbDataProc::FineAlignment","dz = %f  z = %f",dz,z);
    if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateZPar(0,z);
  }

  ali.SelectLongTracks(ali.Npatterns());
  ali.MakeSummaryTracks();
  float shr=0;
  for( i=0; i<ali.Npatterns(); i++ ) {
    fctr = ali.FineCorrShr(i,shr);
    if(fctr<=fcMin) break;
    eDataSet->GetPiece(i)->CorrectShrinkage(0,shr);
    if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateShrPar(0);
  }

}

///______________________________________________________________________________
void EdbDataProc::AlignLinkTracks(int alg, int doAlign)
{
  EdbPVRec    *ali  = new EdbPVRec();
  InitVolume(ali);

  ali->Align(doAlign);
  if(gEDBDEBUGLEVEL>1) ali->PrintAff();
  EdbAffine2D  aff;
  for(int i=0; i<ali->Npatterns(); i++) {
    ali->GetPattern(i)->GetKeep(aff);
    if(!NoUpdate())   eDataSet->GetPiece(i)->UpdateAffPar(0,aff);
  }

  ali->Link();
  Log(2,"EdbDataProc::AlignLinkTracks","link ok");

//    if(alg==1) {
//      ali->MakeHoles();
//      ali->Link();
//      printf("link ok\n");
//    } 

  ali->FillTracksCell();

  TTree *cptree=EdbDataPiece::InitCouplesTree("linked_couples.root","RECREATE");
  FillCouplesTree(cptree, ali,0);
  CloseCouplesTree(cptree);

  MakeTracksTree(ali);
}

//______________________________________________________________________________
int EdbDataProc::MakeTracksTree(EdbPVRec *ali, const char *file)
{
  if(!ali) return 0;
  TObjArray *trarr = ali->eTracks;
  if(!trarr) return 0;
  float xv=ali->X();
  float yv=ali->Y();
  
  return MakeTracksTree(*trarr,xv,yv, file);
}

int EdbDataProc::MakeVertexTree(TObjArray &vtxarr, const char *file)
{
  Log(2,"EdbDataProc::MakeVertexTree","write vertices into %s ... ",file);
  TFile fil(file,"RECREATE");
  TTree *vtx= new TTree("vtx","Reconstructed vertices in emulion");
  
  EdbSegP      *tr = new EdbSegP();
  EdbTrackP *track = NULL;

  TClonesArray *tracks = new TClonesArray("EdbSegP"); //tracks are saved as EdbSegP, like in MakeTracksTree
  TClonesArray *segments  = new TClonesArray("EdbSegP");
  TClonesArray *segmentsf = new TClonesArray("EdbSegP");

  //list of variables to be stored

  Float_t vx, vy, vz; //true reconstructed vertex position
  TMatrixD vCOV(3,3); //covariance matrix vertex
  Float_t meanvx, meanvy, meanvz; //track connection point
  Int_t vID = 0;
  Float_t maxaperture;
  Float_t probability;
  Int_t n;
  Int_t flag; //we need also the flag to check vertex type
  const Int_t maxdim = 1000; //maximum number of tracks associated to a vertex that can be saved in the tree
  //big arrays for containers of track variables
  Int_t TrackID[maxdim];
  Int_t nholes[maxdim];
  Int_t maxgap[maxdim];
  Int_t nseg[maxdim];
  Int_t npl[maxdim];
  Float_t TX[maxdim];
  Float_t TY[maxdim];
  Float_t impactparameter[maxdim];
  Int_t incoming[maxdim]; //coming from the vertex
  //MC information
  Int_t MCEventID[maxdim];
  Int_t MCTrackID[maxdim];
  Int_t MCTrackPdgCode[maxdim];
  Int_t MCMotherID[maxdim];

  //list of branches
  vtx->Branch("vID",&vID,"vID/I");
  vtx->Branch("flag",&flag,"flag/I");
  vtx->Branch("vx",&vx,"vx/F");
  vtx->Branch("vy",&vy,"vy/F");
  vtx->Branch("vz",&vz,"vz/F");
  vtx->Branch("vCOV",&vCOV);
  vtx->Branch("meanvx",&meanvx,"meanvx/F");
  vtx->Branch("meanvy",&meanvy,"meanvy/F");
  vtx->Branch("meanvz",&meanvz,"meanvz/F");
  vtx->Branch("maxaperture",&maxaperture,"maxaperture/F");
  //vtx->Branch("maxrmsthetaspace",&maxrmsthetaspace,"maxrmsthetaspace/F");
  vtx->Branch("probability",&probability,"probability/F");
  vtx->Branch("n",&n,"n/I");
  vtx->Branch("t.",&tracks); // track data is now filled after using Copy method
  //tracks->Branch("t.","EdbSegP",&track,32000,99);
  vtx->Branch("s", &segments);
  vtx->Branch("sf",&segmentsf);
  //track variables (they are array with the number of tracks as size)
  vtx->Branch("TrackID",&TrackID,"TrackID[n]/I");
  vtx->Branch("nseg",&nseg,"nseg[n]/I");
  vtx->Branch("npl",&npl,"npl[n]/I");
  vtx->Branch("nholes",&nholes,"nholes[n]/I"); //even more variables in this tree
  vtx->Branch("maxgap",&maxgap,"maxgap[n]/I"); 
  vtx->Branch("incoming",&incoming,"incoming[n]/I");
  vtx->Branch("impactparameter",&impactparameter,"impactparameter[n]/F");
  //inserting MCtrue information
  vtx->Branch("MCEventID", &MCEventID, "MCEventID[n]/I");
  vtx->Branch("MCTrackID",&MCTrackID,"MCTrackID[n]/I");
  vtx->Branch("MCTrackPdgCode",&MCTrackPdgCode,"MCTrackPdgCode[n]/I");
  vtx->Branch("MCMotherID",&MCMotherID,"MCMotherID[n]/I");

  int nvtx = vtxarr.GetEntriesFast();
  //START TREE FILLING
  for(int ivtx=0; ivtx<nvtx; ivtx++) {
    int itotalseg = 0;
    tracks->Clear("C");
    segments->Clear("C");
    segmentsf->Clear("C");

    EdbVertex* vertex = (EdbVertex*) (vtxarr.At(ivtx));
    if(vertex->Flag()<0) continue; //saving only 'true' vertices in the tree file and in the object
    vx=vertex->VX();
    vy=vertex->VY();
    vz=vertex->VZ();
    meanvx=vertex->X();
    meanvy=vertex->Y();
    meanvz=vertex->Z();
    
    //covariance matrix, need to do it by hand due to no known conversion command from SMatrix to TMatrix
    for (int irow = 0; irow < 3; irow++){
     for (int icolumn = 0; icolumn < 3; icolumn++){
      vCOV(irow,icolumn) = vertex->V()->VCOV()(irow,icolumn);
     }
    }
    
    n=vertex->N();
    maxaperture = vertex->MaxAperture();
    probability = vertex->V()->prob();  
    flag = vertex->Flag();
    //adding vertex to list to be saved
    //loop on tracks //now it can be done offline (again)
    for (int itrk = 0; itrk < n; itrk++){
     //getting track and track variables to fill branches
     track = vertex->GetTrack(itrk);
     tr->Copy(*track);     
//     tr->ForceCOV(track->COV());
     if(tr) new((*tracks)[itrk])  EdbSegP( *tr ); //adding track to trackclonesarray
     ((EdbSegP*) tracks->At(itrk))->ForceCOV(track->COV()); //default copy does NOT save COV corrrectly
     TrackID[itrk] = track->Track(); //the eTrack attribute of EdbSegP now allows association to original root tree
     nseg[itrk] = track->N();
     npl[itrk]  = track->Npl();
     Int_t zpos = vertex->GetVTa(itrk)->Zpos();
     incoming[itrk] = zpos;
     nholes[itrk] = track->N0();
     maxgap[itrk] = track->CheckMaxGap();
     impactparameter[itrk] = vertex->GetVTa(itrk)->Imp();
     //Storing MCTrue information (of course for real data these values have no sense)
     MCEventID[itrk] = track->MCEvt();
     MCTrackID[itrk] = track->MCTrack();
     if(MCEventID[itrk]>-1){ //vertices from MC simulation, used Vid[0] and Aid[0] to store these information
      MCTrackPdgCode[itrk] = track->GetSegment(0)->Vid(0); //due to MakeTracksTree, Vid[0] is always 0
      MCMotherID[itrk] = track->Aid(0); //used to store MotherID information
     }
     else { //vertices from data, set default values 
      MCTrackPdgCode[itrk] = -999;
      MCMotherID[itrk] = -999;
     }
    //w    = track->Wgrains();
     EdbSegP *s=0,*sf=0;
    //loop on segments
     for(int is=0; is<nseg[itrk]; is++) {
      s = track->GetSegment(is);     
      if(s) new((*segments)[itotalseg])  EdbSegP( *s );
      ((EdbSegP*) segments->At(itotalseg))->ForceCOV(s->COV()); //default copy does NOT save COV corrrectly
      sf = track->GetSegmentF(is);
      if(sf) new((*segmentsf)[itotalseg])  EdbSegP( *sf );
      ((EdbSegP*) segmentsf->At(itotalseg))->ForceCOV(sf->COV()); //default copy does NOT save COV corrrectly
      itotalseg++;
      }

   // track->SetVid( 0, tracks->GetEntries() );  // put track counter in t.eVid[1]
     }
    vtx->Fill();
    vID++; //now the number of elements in the tree and vertices is the same. vID starts from 0
    }
  vtx->Write();
  fil.Close();
  Log(2,"EdbDataProc::MakeVertexTree","%d vertices are written",nvtx);
  return nvtx; 
}
int EdbDataProc::ReadVertexTree( EdbVertexRec &vertexrec, const char     *fname, const char *rcut, map<int,EdbTrackP*> trackID_map)
{
  TFile f(fname);
  if(f.IsZombie()) { Log(1,"EdbDataProc::ReadVertexTree","Error open file %s", fname);  return 0; }
  
  TTree *vtx = (TTree*)f.Get("vtx");

  //vertex variables
  Float_t vx, vy, vz;
  Int_t vID = 0;
  Int_t n = 0;
  Int_t flag = 0;
  //track variables

  const Int_t maxdim = 1000; //maximum number of tracks associated to a vertex that can be saved in the tree
  //big arrays for containers of track variables
  Int_t TrackID[maxdim];
  Int_t nholes[maxdim];
  Int_t maxgap[maxdim];
  Int_t nseg[maxdim];
  Float_t TX[maxdim];
  Float_t TY[maxdim];
  Float_t impactparameter[maxdim];
  Int_t incoming[maxdim]; //coming from the vertex
  //MC information
  Int_t MCEventID[maxdim];
  Int_t MCTrackID[maxdim];
  Int_t MCMotherID[maxdim];
/*
  Int_t   trid=0;
  Int_t   nseg=0;
  Int_t   npl=0;
  Int_t   n0=0;
  Float_t xv=0.;
  Float_t yv=0.;*/

  int nentr = (int)(vtx->GetEntries());
  TCut cut = rcut;
  vtx->Draw(">>lst", cut );
  TEventList *lst = (TEventList*)gDirectory->GetList()->FindObject("lst");
  int nlst =lst->GetN();
  if(cut) Log(2,"EdbDataProc::ReadVtxTree","select %d of %d vertices by cut %s",nlst, nentr, cut.GetTitle() );

  TClonesArray *tracks = new TClonesArray("EdbSegP");
  TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
  TClonesArray *segf = new TClonesArray("EdbSegP", 60);
  
  EdbSegP *trk=0;
  EdbSegP *s1=0; 
  EdbSegP *s1f=0;
  //vertex variables
  vtx->SetBranchAddress("vID",&vID);
  vtx->SetBranchAddress("flag",&flag);
  vtx->SetBranchAddress("vx",&vx);
  vtx->SetBranchAddress("vy",&vy);
  vtx->SetBranchAddress("vz",&vz);
  vtx->SetBranchAddress("n",&n);
  //arrays
  vtx->SetBranchAddress("nseg",&nseg);
  vtx->SetBranchAddress("TrackID",&TrackID);
  vtx->SetBranchAddress("incoming",&incoming);
  //TClones arrays
  vtx->SetBranchAddress("sf", &segf);
  vtx->SetBranchAddress("s",  &seg);
  vtx->SetBranchAddress("t.", &tracks);

  EdbPattern *pat=0;
  int entr=0;
  //now each tree entry is a vertex

  TObjArray * vertextracks = new TObjArray(maxdim);
  float expz = 0.;

  EdbPVRec *ali = vertexrec.ePVR; //I get EdbPVR from EdbVertexRec from now on
  for (int j=0; j<nlst; j++){
    //RESET COUNTERS!
    int itotalseg = 0; //counter for all the segments
    vertextracks->Clear("C");

    entr = lst->GetEntry(j);
    vtx->GetEntry(entr);
    EdbVertex *v1 = new EdbVertex();
    //loop on all tracks associated to the vertex
    for (int itrk = 0; itrk < n; itrk++){
     EdbTrackP *tr1;
     if(!trackID_map[TrackID[itrk]]){
      tr1= new EdbTrackP();
      trk = (EdbSegP*) tracks->At(itrk); //getting track itrk
      ((EdbSegP*)tr1)->Copy(*trk);
      tr1->ForceCOV(trk->COV());
      tr1->SetM(0.139);                 //TODO

      for(int i=0; i<nseg[itrk]; i++) {
       s1 = (EdbSegP*)(seg->At(itotalseg));
       s1f = (EdbSegP*)(segf->At(itotalseg));
       pat = ali->GetPattern( s1->PID() );
       if(!pat) { 
	      Log(1,"EdbDataProc::ReadVertexTree","WARNING: no pattern with pid %d: creating new one!",s1->PID()); 
	      pat = new EdbPattern( 0., 0., s1->Z() );
	      pat->SetID(s1->PID());
	      pat->SetScanID(s1->ScanID());
	      ali->AddPatternAt(pat,s1->PID());
      }

        //This way it should be more clear which segment instance is added to the track
        EdbSegP *segtobeadded = pat->AddSegment(*s1);
        EdbSegP *segftobeadded = new EdbSegP(*s1f);
        segtobeadded->ForceCOV(s1->COV());
        segftobeadded->ForceCOV(s1f->COV());

        tr1->AddSegment(segtobeadded);
        tr1->AddSegmentF(segftobeadded);
        itotalseg++; //increasing counter of number of segments
      }
      tr1->SetSegmentsTrack(tr1->ID());
      tr1->SetCounters();
     //tr1->FitTrackKFS(true);
      tr1->SetTrack(TrackID[itrk]); //providing trackid to eTrack so it will not be lost when trackID resets
      ali->AddTrack(tr1);
      trackID_map[TrackID[itrk]] = tr1;
     }
     else{ //track already added to ali, add it only to vertex
      tr1 = trackID_map[TrackID[itrk]];
      itotalseg+= nseg[itrk];
     }
     vertextracks->Add(tr1);     
     //v1 = vertexrec->AddTrackToVertex(v1, tr1, incoming[itrk]);
    }
    //setting vertex parameters and saving vertex
    //v1 = vertexrec.Make1Vertex(*vertextracks,vz) (copied here instead, by putting zpos directly as for in vertexing);
    v1->SetXYZ( 0,0, vz );
    int ntr = vertextracks->GetEntries();
    for(int i=0; i<ntr; i++) {
      EdbTrackP *t = (EdbTrackP*)vertextracks->At(i);
      EdbVTA *vta = new EdbVTA(t,v1);
      vta->SetFlag(2);
      v1->AddVTA(vta);
      vta->SetZpos(incoming[i]);
      t->AddVTA(vta);
    }
    if( vertexrec.MakeV(*v1) )  vertexrec.AddVertex(v1);
    //else { SafeDelete(v); return 0; }                                // vertex is not valid
    //return v;
    v1->SetID(vID);
    v1->SetFlag(flag);
    /* make vertex manually, old version
    int ntr = vertextracks->GetEntries();
    v1->SetXYZ(vx,vy,vz);  
    for(int i=0; i<ntr; i++) {
     EdbTrackP *t = (EdbTrackP*)vertextracks->At(i);
     EdbVTA *vta = new EdbVTA(t,v1);
     vta->SetFlag(2);
     v1->AddVTA(vta);
     (t->Z() >= v1->Z())? vta->SetZpos(1) : vta->SetZpos(0);
     t->AddVTA(vta);
    } */
    
    ali->AddVertex(v1);
  }

  Log(2,"EdbDataProc::ReadVtxTree","%d vertices are read",nlst);
  return nlst;
}


 int EdbDataProc::ReadVertexTree( EdbVertexRec &vertexrec, const char     *fname, const char *rcut,TObjArray *builttracks){

  map<int,EdbTrackP*>emptymap;
  if (builttracks){
   int ntracks = builttracks->GetEntries();
  
   for (int itrk = 0; itrk < ntracks; itrk++){
        EdbTrackP * track = (EdbTrackP*) builttracks->At(itrk);
        emptymap[track->Track()] = track; //adding this trackID to map
    }
  }

  int nlst = ReadVertexTree( vertexrec, fname, rcut,emptymap);
  return nlst;
 }

 EdbVertex*  EdbDataProc::GetVertexFromTree( EdbVertexRec &vertexrec, const char     *fname, const int vertexID )
{
  //reading vertex tree, getting only the vertex I need
  map<int,EdbTrackP*>emptymap;
  ReadVertexTree(vertexrec, fname, Form("vID==%i",vertexID),emptymap); //vertex is added to EdbPVRec

  EdbVertex *myvertex;
  EdbPVRec *ali = vertexrec.ePVR;
  //looking for myvertex (when called more than once, ali keeps adding vertices)
  for (int ivtx = 0; ivtx<ali->eVTX->GetEntries();ivtx++){
    EdbVertex *v1 = (EdbVertex*) ali->eVTX->At(ivtx);
    if (v1->ID()==vertexID) myvertex = v1;
  }
  
  return myvertex;

 }

//______________________________________________________________________________
int EdbDataProc::MakeTracksTree(TObjArray &trarr, float xv, float yv, const char *file)
{
  Log(2,"EdbDataProc::MakeTracksTree","write tracks into %s ... ",file);
  TFile fil(file,"RECREATE");
  TTree *tracks= new TTree("tracks","tracks");

  EdbTrackP    *track = new EdbTrackP(8);   // why is this track initialised with 8 segments by default ???
  EdbSegP      *tr = new EdbSegP();
  TClonesArray *segments  = new TClonesArray("EdbSegP");
  TClonesArray *segmentsf = new TClonesArray("EdbSegP");

  int   nseg,trid,npl,n0;
  float w=0.;

  tracks->Branch("trid",&trid,"trid/I");
  tracks->Branch("nseg",&nseg,"nseg/I");
  tracks->Branch("npl",&npl,"npl/I");
  tracks->Branch("n0",&n0,"n0/I");
  tracks->Branch("xv",&xv,"xv/F");
  tracks->Branch("yv",&yv,"yv/F");
  tracks->Branch("w",&w,"w/F");
  tracks->Branch("t.","EdbSegP",&tr,32000,99); // track data is now filled after using Copy method
  //tracks->Branch("t.","EdbSegP",&track,32000,99);
  tracks->Branch("s", &segments);
  tracks->Branch("sf",&segmentsf);

  int ntr = trarr.GetEntriesFast();
  
  for(int itr=0; itr<ntr; itr++) {
    
    track = (EdbTrackP*)(trarr.At(itr));
    tr->Copy(*track); //Using Copy method from EdbSegP avoids the cast
    tr->ForceCOV(track->COV());
    trid = track->ID();
    nseg = track->N();
    npl  = track->Npl();
    n0   = track->N0();
    
    segments->Clear("C");
    segmentsf->Clear("C");
    nseg = track->N();
    w    = track->Wgrains();
    EdbSegP *s=0,*sf=0;
    for(int is=0; is<nseg; is++) {
      s = track->GetSegment(is);
      if(s) new((*segments)[is])  EdbSegP( *s );
      sf = track->GetSegmentF(is);
      if(sf) new((*segmentsf)[is])  EdbSegP( *sf );
    }

    tr->SetVid( 0, tracks->GetEntries() );  // put track counter in t.eVid[1]
    tracks->Fill();
    // track->Clear(); // if this Clear() is written at this stage, the corresponding track in
    // the EdbPattern Volume (ali, eTracks) will also be cleared! This behavior is not
    // wanted, is it? So we do comment this out here ... (Frank, 09 26 2016)
  }

  tracks->Write();
  fil.Close();
  Log(2,"EdbDataProc::MakeTracksTree","%d tracks are written",ntr);
  return ntr; 
}

//---------------------------------------------------------------------------
int EdbDataProc::ReadTracksTree( EdbPVRec &ali,
				 const char     *fname,
				 //				 int      nsegMin,
				 //				 float    probMin,
				 const char *rcut )
{
  TFile f(fname);
  if(f.IsZombie()) { Log(1,"EdbDataProc::ReadTracksTree","Error open file %s", fname);  return 0; }
  
  TTree *tracks = (TTree*)f.Get("tracks");
  Int_t   trid=0;
  Int_t   nseg=0;
  Int_t   npl=0;
  Int_t   n0=0;
  Float_t xv=0.;
  Float_t yv=0.;

  int nentr = (int)(tracks->GetEntries());
  TCut cut = rcut;
  tracks->Draw(">>lst", cut );
  TEventList *lst = (TEventList*)gDirectory->GetList()->FindObject("lst");
  int nlst =lst->GetN();
  if(cut) Log(2,"EdbDataProc::ReadTracksTree","select %d of %d tracks by cut %s",nlst, nentr, cut.GetTitle() );

  TClonesArray *seg  = new TClonesArray("EdbSegP", 60);
  TClonesArray *segf = new TClonesArray("EdbSegP", 60);
  EdbSegP *trk=0;
  EdbSegP *s1=0;
  EdbSegP *s1f=0;

  tracks->SetBranchAddress("trid", &trid);
  tracks->SetBranchAddress("nseg", &nseg);
  tracks->SetBranchAddress("npl", &npl);
  tracks->SetBranchAddress("n0", &n0);
  tracks->SetBranchAddress("xv", &xv);
  tracks->SetBranchAddress("yv", &yv);
  tracks->SetBranchAddress("sf", &segf);
  tracks->SetBranchAddress("s",  &seg);
  tracks->SetBranchAddress("t.", &trk);

  EdbPattern *pat=0;
  int entr=0;
  for (int j=0; j<nlst; j++){
    entr = lst->GetEntry(j);
    tracks->GetEntry(entr);

    EdbTrackP *tr1 = new EdbTrackP();
    ((EdbSegP*)tr1)->Copy(*trk);
    tr1->ForceCOV(trk->COV());
    tr1->SetM(0.139);                 //TODO

    for(int i=0; i<nseg; i++) {
      s1 = (EdbSegP*)(seg->At(i));
      s1f = (EdbSegP*)(segf->At(i));
      pat = ali.GetPattern( s1->PID() );
      if(!pat) { 
	Log(1,"EdbDataProc::ReadTracksTree","WARNING: no pattern with pid %d: creating new one!",s1->PID()); 
	pat = new EdbPattern( 0., 0., s1->Z() );
	pat->SetID(s1->PID());
	pat->SetScanID(s1->ScanID());
	ali.AddPatternAt(pat,s1->PID());
      }
      //This way it should be more clear which segment instance is added to the track
      EdbSegP *segtobeadded = pat->AddSegment(*s1);
      EdbSegP *segftobeadded = new EdbSegP(*s1f);
      segtobeadded->ForceCOV(s1->COV());
      segftobeadded->ForceCOV(s1f->COV());

      tr1->AddSegment( segtobeadded );
      tr1->AddSegmentF( segftobeadded );
    }
    tr1->SetSegmentsTrack(tr1->ID());
    tr1->SetCounters();
    //tr1->FitTrackKFS(true);
    tr1->SetTrack(trid); //providing trackid to eTrack so it will not be lost when trackID resets
    ali.AddTrack(tr1);
  }

  Log(2,"EdbDataProc::ReadTracksTree","%d tracks are read",nlst);
  return nlst;
}

//---------------------------------------------------------------------------
TIndex2 *EdbDataProc::MakeTracksSegmentsList( EdbPVRec &ali)
{
  TIndex2 *itracks = 0;
  int nsegtot=0;
  EdbTrackP *tr=0;
  int ntr = ali.eTracks->GetEntries();

  if(!ntr) return 0;
  Double_t *w    = new Double_t[500000];
  itracks = new TIndex2();

  for(int i=0; i<ntr; i++) {
    tr = (EdbTrackP*)(ali.eTracks->At(i));
    if(!tr) continue;
    EdbSegP *s=0;
    int nseg=tr->N();
    for(int j=0; j<nseg; j++) {
      s = tr->GetSegment(j);
      if(!s) continue;
      w[nsegtot] = TIndex2::BuildValue( s->Vid(0), s->Vid(1) );
      nsegtot++;
    }
  }
  itracks->Set(0);
  itracks->BuildIndex(nsegtot,w);
  delete[] w;
  return itracks;
}
