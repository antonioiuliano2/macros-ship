//-- Author :  Valeri Tioukov   22/12/2006
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbScanProc                                                          //
//                                                                      //
// scanned data processing library                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <stdio.h>
#include <time.h>
#include <sys/stat.h>

#include "Varargs.h"
#include "TROOT.h"
#include "TSystem.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TEnv.h"
#include "TChain.h"
#include "EdbLog.h"
#include "EdbTrackFitter.h"
#include "EdbScanProc.h"
#include "EdbTestAl.h"
#include "EdbPlateAlignment.h"
#include "EdbAlignmentMap.h"
#include "EdbRunTracking.h"
#include "EdbCouplesTree.h"
#include "EdbLinking.h"

using namespace std;
using namespace TMath;

ClassImp(EdbScanProc)
//----------------------------------------------------------------
EdbScanProc::EdbScanProc()
{
}

//----------------------------------------------------------------
void EdbScanProc::Print()
{
  printf("ProcDirClient: %s\n",eProcDirClient.Data());
}

//----------------------------------------------------------------
bool EdbScanProc:: AddAFFtoScanSet(EdbScanSet &sc, int b1, int p1, int s1, int e1,int b2, int p2, int s2, int e2)
{
  int id1[4]={b1,p1,s1,e1};
  int id2[4]={b2,p2,s2,e2};
  return AddAFFtoScanSet(sc,id1,id2);
}

//----------------------------------------------------------------
bool EdbScanProc:: AddAFFtoScanSet(EdbScanSet &sc, EdbID eid1, EdbID eid2)
{
  int id1[4]; eid1.Get(id1);
  int id2[4]; eid2.Get(id2);
  return AddAFFtoScanSet(sc,id1,id2);
}

//----------------------------------------------------------------
bool EdbScanProc::ReadAffToLayer(EdbLayer &la, EdbID id1, EdbID id2)
{
  // read affine tranformation from file, into layer
  EdbDataPiece piece;
  TString parfile;
  MakeAffName(parfile,id1,id2);
  piece.eFileNamePar = parfile;
  if (piece.TakePiecePar() < 0) { Log(2,"ReadAFFtoLayer","file %s do not found...", parfile.Data() );    return 0; }
  la.Copy( *(piece.GetLayer(0)) );
  la.SetZcorr(piece.GetLayer(0)->Z());
  la.SetZlayer(0.,0.,0.);
  if(gEDBDEBUGLEVEL>2) la.Print();
  return true;
}

//----------------------------------------------------------------
bool EdbScanProc:: AddAFFtoScanSet(EdbScanSet &sc, int id1[4], int id2[4])
{
  // read affine tranformation from file, form "plate", add "plate" to EdbScanSet:ePC

  EdbDataPiece piece;  
  TString parfile;
  MakeAffName(parfile,id1,id2);
  piece.eFileNamePar = parfile;
  float dz=0;
  EdbAffine2D *a=0;
  if (piece.TakePiecePar() < 0) {
	dz = 1300.*( id2[1]-id1[1] );
	a = new EdbAffine2D();
	Log(2,"AddAFFtoScanSet","add default transform with dz = %f",dz);
  } else {
	dz = piece.GetLayer(0)->Z();
	a = piece.GetLayer(0)->GetAffineXY();
	Log(2,"AddAFFtoScanSet","add read transform with dz = %f",dz);
  }
  if(gEDBDEBUGLEVEL>2) a->Print();
  EdbPlateP *p = new EdbPlateP();   // in this function  1-one plate, 2-next plate
  p->GetLayer(1)->SetID(id1[1]);
  p->GetLayer(1)->SetZlayer(-dz,-dz,-dz);
  p->GetLayer(1)->SetAffXY(a->A11(),a->A12(),a->A21(),a->A22(),a->B1(),a->B2());
  p->GetLayer(2)->SetID(id2[1]);
  p->GetLayer(2)->SetZlayer(0,0,0);
  p->GetLayer(2)->SetAffXY(1,0,0,1,0,0);

  sc.ePC.Add(p);

  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::ReadPiecePar(EdbID id, EdbPlateP &plate)
{
    EdbDataPiece piece;
    TString str;
    MakeFileName(str,id,"par");
    piece.ReadPiecePar(str);
    for(int i=0; i<3; i++) {
      EdbLayer *l=piece.GetLayer(i);
      if(l) {
        EdbAffine2D *a = l->GetAffineXY();
        plate.GetLayer(i)->SetAffXY( a->A11(),a->A12(),a->A21(),a->A22(),a->B1(),a->B2() );
        a = l->GetAffineTXTY();
        plate.GetLayer(i)->SetAffTXTY( a->A11(),a->A12(),a->A21(),a->A22(),a->B1(),a->B2() );
        plate.GetLayer(i)->SetShrinkage( l->Shr() );
      }
    }
    return true;
}

//----------------------------------------------------------------
int EdbScanProc::AssembleScanSet(EdbScanSet &sc)
{
  if(sc.eIDS.GetSize() < 1)   return 0;
  else if (sc.eIDS.GetSize() == 1) {  // add nominal plate
    EdbPlateP *plate = new EdbPlateP();
    plate->SetID(sc.GetID(0)->ePlate);
    Float_t z=0, dz0=214,dz1=45,dz2=45;  //TODO!
    plate->SetZlayer(z, z - dz0/2 + dz1, z+dz0/2+dz2);                
    plate->GetLayer(0)->SetZlayer(0,-dz0/2,dz0/2);       // internal plate coord
    plate->GetLayer(2)->SetZlayer(-dz0/2,-dz0/2-dz2,-dz0/2);
    plate->GetLayer(1)->SetZlayer( dz0/2, dz0/2, dz0/2+dz1);
    sc.eB.AddPlate(plate);
    sc.MakePIDList();
    return 1;
  }

  sc.ePC.Delete();  // clear plate array

  EdbID *id1=0, *id2=0;
  
  for (Int_t i = 1; i < sc.eIDS.GetSize(); i++) {
    id1 = (EdbID *)(sc.eIDS.At(i-1));
    id2 = (EdbID *)(sc.eIDS.At(i));
    if (!AddAFFtoScanSet(sc, *id1, *id2)) return -1;
  }

  return sc.AssembleBrickFromPC();
}

//----------------------------------------------------------------
void EdbScanProc::CopyParSet(EdbID idset1,  EdbID idset2)
{
  // assuming that exist the scan sets for idset1 and idset2 
  // copy the plate par files with shrinkage corrections from 1 to 2
  EdbScanSet *ss1 = ReadScanSet(idset1);   if(!ss1) return;
  EdbScanSet *ss2 = ReadScanSet(idset2);   if(!ss2) return;
  int n = ss1->eIDS.GetSize();
  for(int i=0; i<n; i++) {
    EdbID *id1  = ss1->GetID(i);
    EdbID *id2  = ss2->FindPlateID(id1->ePlate);
    if(id2) CopyPar(*id1,*id2);
  }
  
}

//----------------------------------------------------------------
int EdbScanProc::ReadFoundSegment(EdbID id,  EdbSegP &s, int flag)
{
  EdbPattern pat;
  if (!ReadFound(pat, id, flag)) return 0;
  EdbSegP *ss = pat.FindSegment(s.ID());
  if (!ss) return 0;
  s.Copy(*ss);
  return 1;
}

//----------------------------------------------------------------
int EdbScanProc::ReadFoundTrack(EdbScanSet &sc,  EdbTrackP &track, int flag)
{
  // read track.ID() from pieces in sc.IDS .found.root and apply transformations from sc

  EdbPlateP  *plate;
  int count=0;
  int n = sc.eIDS.GetSize();
  EdbPattern pat;
  pat.AddSegment(track.ID(), 0,0,0,0);
  EdbSegP *s = pat.GetSegment(0);

  for(int i=0; i<n; i++) {
    EdbID *id    = (EdbID *)(sc.eIDS.At(i));
    if (ReadFoundSegment(*id,*s,flag) <= 0) continue;
    count++;
    plate = sc.GetPlate(id->ePlate);
    pat.Transform(    plate->GetAffineXY()   );
    pat.TransformA(   plate->GetAffineTXTY() );
    pat.TransformShr( plate->Shr() );
    s->SetZ(plate->Z());
    s->SetDZ(-214);                            //TODO!!!
    s->SetPID(id->ePlate);
    track.AddSegment( new EdbSegP(*s) );
  }
  track.SetCounters();
  return count;
}

//----------------------------------------------------------------
int EdbScanProc::ReadFoundTracks(EdbScanSet &sc,  EdbPVRec &ali, int flag)
{
  // read all tracks  from found.root listed in sc.IDS and apply transformations from sc.eB
  // return the total number of segments added

  EdbPlateP  *plate;
  int n = sc.eIDS.GetSize();
  int count=0;

  for(int i=0; i<n; i++) {
    EdbID *id  = sc.GetID(i);
    EdbPattern pat;
    ReadFound(pat, *id, flag);
    plate = sc.GetPlate(id->ePlate);
    pat.SetPID(id->ePlate);
    pat.SetID(id->ePlate);
    pat.SetSegmentsPID();
    pat.Transform(    plate->GetAffineXY()   );
    pat.TransformA(   plate->GetAffineTXTY() );
    pat.TransformShr( plate->Shr() );
    pat.SetZ(plate->Z());
    pat.SetSegmentsZ();
    pat.SetSegmentsDZ(plate->DZ());  // to check the direction!
    
    for(int j=0; j<pat.N(); j++) {
      EdbSegP *s = pat.GetSegment(j);
      EdbTrackP *track = ali.FindTrack(s->ID());
      if(track)  track->AddSegment( new EdbSegP(*s) );
      else       ali.AddTrack( new EdbTrackP(new EdbSegP(*s)) );
      count++;
    }
  }

  int ntr = ali.Ntracks();
  for(int i=0; i<ntr; i++) ali.GetTrack(i)->SetCounters();

  LogPrint(sc.eB.ID(),2,"ReadFoundTracks","%d  segments read\n",count );

  return count;
}

//----------------------------------------------------------------
int EdbScanProc::ReadManFoundTracks(EdbScanSet &sc,  EdbPVRec &ali, int flag)
{
  // read all tracks  from *man.found.txt listed in sc.IDS and apply transformations from sc.eB
  // return the total number of segments added

  EdbPlateP  *plate;
  int n = sc.eIDS.GetSize();
  int count=0;

  for(int i=0; i<n; i++) {
    EdbPattern pat;
    EdbID *id  = sc.GetID(i);                 if(!id) continue;
    if( !ReadPatTXT(pat, *id, "man.found.txt",flag) ) continue;
    plate = sc.GetPlate(id->ePlate);       if(!plate) continue;
    pat.SetPID(id->ePlate);
    pat.SetID(id->ePlate);
    pat.SetSegmentsPID();
    pat.Transform(    plate->GetAffineXY()   );
    pat.TransformA(   plate->GetAffineTXTY() );
    pat.TransformShr( plate->Shr() );
    pat.SetZ(plate->Z());
    //pat.SetSegmentsZ();
    //pat.SetSegmentsDZ(plate->DZ());  // to check the direction!
    
    for(int j=0; j<pat.N(); j++) {
      EdbSegP *s = pat.GetSegment(j);
      s->SetVid(id->ePlate,0);
      int side = s->Flag()%10;
      if(side==0) s->SetZ(plate->Z());
      if(side==1) s->SetZ(plate->Z()+107);
      if(side==2) s->SetZ(plate->Z()-107);
      EdbTrackP *track = ali.FindTrack(s->ID());
      if(track)  track->AddSegment( new EdbSegP(*s) );
      else       ali.AddTrack( new EdbTrackP(new EdbSegP(*s)) );
      count++;
    }
  }
  
  n = ali.Ntracks();
  for(int i=0; i<n; i++) {
    ali.GetTrack(i)->SetCounters();
  }

  LogPrint(sc.eB.ID(),2,"ReadManFoundTracks","%d  segments read\n",count );

  return count;
}

//----------------------------------------------------------------
void EdbScanProc::MakeEraseFile(EdbID id,  EdbPattern &pat)
{
  // input: pat with the segments to be erased - assumed that s.eVid[1] is the entry number in the couples tree
  // ouput: file id.er.root with the EdbMask object "mask" inside
  int entrmax=0;
  for(int i=0; i<pat.N(); i++) if(pat.GetSegment(i)->Vid(1)>entrmax) entrmax=pat.GetSegment(i)->Vid(1);
  if(!entrmax)  return;
  EdbMask mask(entrmax+1);
  for(int i=0; i<pat.N(); i++) mask.SetAt(pat.GetSegment(i)->Vid(1),1);
  TString str;
  MakeFileName(str,id,"er.root");
  TFile f(str.Data(),"RECREATE");
  mask.Write("mask");
  f.Close();
}

//----------------------------------------------------------------
int EdbScanProc::ReadScanSetCP(EdbID id,  EdbPVRec &ali, TCut c, bool do_erase, bool do_assemble, int minplate, int maxplate)
{
  // read data for scanset defined with id  apply cut c and fill ali
  // if(do_erase)    - exclude segmets from the erase mask if it is exist (default is true) 
  // if(do_assemble) - read affine parfiles (default is true) 
  
    EdbScanSet *ss = ReadScanSet(id);
    ss->Brick().SetID(id.eBrick);
    if(do_assemble) AssembleScanSet(*ss);
    return ReadScanSetCP(*ss, ali, c, do_erase, minplate, maxplate);
}

//----------------------------------------------------------------
int EdbScanProc::ReadScanSetCP(EdbScanSet &sc,  EdbPVRec &ali, TCut c, bool do_erase, int minplate, int maxplate)
{
  // read data from scanset sc with cut c and fill ali
  // sc.eIDS is used as an id list
  // sc.eB   is used to get the brick geometry and affine transformations - must be filled before
  // if(do_erase) - exclude segmets from the erase mask if it is exist (default is true) 

  if(minplate>maxplate) { 
    Log(1,"EdbScanProc::ReadScanSetCP","ERROR: minplate (%d) greater then maxplate(%d)!",minplate, maxplate);
    return 0;
  }
  int cnt=0;
  int n = sc.eIDS.GetSize();    // number of pieces to get
  EdbID      *id;
  EdbPlateP  *plate;
  EdbPattern *pat;
  for(int i=0; i<n; i++) {
    id    = (EdbID *)(sc.eIDS.At(i));
    if( !(minplate==-1000&&maxplate==-1000) ) {
      if( id->ePlate < minplate) continue;
      if( id->ePlate > maxplate) continue;
    }
    plate = sc.GetPlate(id->ePlate);
    pat   = new EdbPattern();
    cnt += ReadPatCPnopar( *pat, *id, c, do_erase);

    pat->SetScanID(*id);
    pat->SetPID(id->ePlate);
    for(int j=0; j<pat->N(); j++) pat->GetSegment(j)->SetVid( id->ePlate, pat->GetSegment(j)->Vid(1) );
    pat->SetZ(plate->Z());
    pat->SetSegmentsZ();
    pat->Transform(    plate->GetAffineXY()   );
    pat->TransformA(   plate->GetAffineTXTY() );
    pat->TransformShr( plate->Shr() );
    pat->SetSegmentsPlate(id->ePlate);
    ali.AddPattern( pat );
  }

  ali.SetPatternsID();
  int np = ali.Npatterns();
  for(int i=0; i<np; i++)  ali.GetPattern(i)->SetSegmentsPID();  // PID of the segment must be ID of the pattern!

  // to be moved into the processing part???
  ali.SetSegmentsErrors();
  ali.SetCouplesAll();
  ali.SetChi2Max(ali.GetScanCond()->Chi2PMax());
  ali.SetOffsetsMax(ali.GetScanCond()->OffX(),ali.GetScanCond()->OffY());
  
  return cnt;
}

//----------------------------------------------------------------
bool EdbScanProc::MakeAFFSet(EdbScanSet &sc)
{
  // create AFF dir if do not exist 
  // put all affine transformations inside for each plates couple

  if(!CheckAFFDir(sc.eID.eBrick)) return false;
  if(sc.eIDS.GetSize()<2) return 0;
  EdbID *id1,*id2;
  for(int i=0; i<sc.eIDS.GetSize()-1; i++) {
    id1 = (EdbID *)(sc.eIDS.At(i));
    id2 = (EdbID *)(sc.eIDS.At(i+1));
    EdbAffine2D aff;
    sc.GetAffP2P(id1->ePlate,id2->ePlate, aff);
    float dz = sc.GetDZP2P(id1->ePlate,id2->ePlate);
    TString str;
    MakeAffName(str,*id1,*id2);
    char card[128];
    sprintf(card,"ZLAYER 0 %f 0 0",dz);
    LogPrint(id1->eBrick,2,"MakeAFFSet","%s as %s", str.Data(),card);
    AddParLine(str.Data(),card);
    sprintf(card,"AFFXY 0 %f %f %f %f %f %f", aff.A11(), aff.A12(), aff.A21(), aff.A22(), aff.B1(), aff.B2() );
    LogPrint(id1->eBrick,2,"MakeAFFSet","%s as %s", str.Data(),card);
    AddParLine(str.Data(),card);
  }
  return 1;
}

//----------------------------------------------------------------
bool  EdbScanProc::MakeParSet(EdbScanSet &sc)
{
  if(sc.eIDS.GetSize()<1)   return false;
  EdbID *id=0;
  for(int i=0; i<sc.eIDS.GetSize(); i++) {
    id = (EdbID *)(sc.eIDS.At(i));
    CheckProcDir(*id);
    TString name;
    MakeFileName(name,*id,"par");
    char card[128];
    sprintf(card,"SHRINK 1 %f",sc.GetPlate(id->ePlate)->GetLayer(1)->Shr());
    AddParLine(name,card,1);
    sprintf(card,"SHRINK 2 %f",sc.GetPlate(id->ePlate)->GetLayer(2)->Shr());
    AddParLine(name,card);
  }
  return 1;
}

//----------------------------------------------------------------
bool  EdbScanProc::PrepareSetStructure(EdbScanSet &sc)
{
  // input: sc with brick and all id's defined
  // function -check or create the directory structure for this brick,
  //          -create par-files for all plates
  //          -create aff.par files for all couples in the defined order
  MakeParSet(sc);
  MakeAFFSet(sc);
  return true;
}

//----------------------------------------------------------------
void EdbScanProc::AlignSetNewNopar(EdbID id, TEnv &cenv )
{
  EdbScanSet *ss = ReadScanSet(id);
  if(!ss) { Log(1,"AlignSetNewNopar","Error! set for %s not found",id.AsString()); return; }
  AlignSetNewNopar(*ss, cenv);
  ss->eB.ResetAff();
  WriteScanSet(id,*ss);
  UpdateSetWithAff(id,id);
}

//----------------------------------------------------------------
int EdbScanProc::AlignSetNewNopar(EdbScanSet &sc, TEnv &cenv )
{
  if(sc.eIDS.GetSize()<2) return 0;
  int n=0;
  int minPlate = cenv.GetValue("fedra.align.minPlate"  ,-999);
  int maxPlate = cenv.GetValue("fedra.align.maxPlate"  , 999);

  bool do_map = cenv.GetValue("fedra.align.do_map"  , 0);
  for(int i=0; i<sc.eIDS.GetSize()-1; i++) {
    EdbID *id1 = sc.GetID(i);
    EdbID *id2 = sc.GetID(i+1);
    EdbPlateP *plate2 = sc.GetPlate(id2->ePlate);
    if(id1->ePlate<minPlate||id1->ePlate>maxPlate) continue;
    if(id2->ePlate<minPlate||id2->ePlate>maxPlate) continue;
    EdbAffine2D aff;
    float dz = -1300;
    if(sc.GetAffP2P(id1->ePlate, id2->ePlate, aff))
      dz = sc.GetDZP2P(id1->ePlate, id2->ePlate);
    if (!do_map) n += AlignNewNopar(*id1, *id2, cenv, &aff, dz);
    else n+= AlignNewLayernopar(*id1, *id2, cenv, *plate2, &aff, dz);
  }
  return n;
}

//-------------------------------------------------------------------
int EdbScanProc::AlignNewNopar(EdbID id1, EdbID id2, TEnv &cenv, EdbAffine2D *aff, float dz)
{
  // Align 2 patterns. All necessary information should be in the envfile
  // Convension about Z(setted while process): the z of id2 is 0, the z of id1 is (-deltaZ) where
  // deltaZ readed from aff.par file in a way that pattern of id1 projected 
  // to deltaZ correspond to pattern of id2

  int npat=0;
  
  EdbPlateAlignment av;
  av.eOffsetMax =   cenv.GetValue("fedra.align.OffsetMax"   , 500. );
  av.SetSigma(      cenv.GetValue("fedra.align.SigmaR"      , 13.  ), 
	            cenv.GetValue("fedra.align.SigmaT"      , 0.008) );
  av.eDoFine      = cenv.GetValue("fedra.align.DoFine"      , 1);
  av.eDZ          = cenv.GetValue("fedra.align.DZ"          , 120.);
  av.eDPHI        = cenv.GetValue("fedra.align.DPHI"        , 0.008 );
  bool  do_erase  = cenv.GetValue("fedra.align.erase"    , false );
  const char *cut = cenv.GetValue("fedra.readCPcut"         , "eCHI2P<2.5&&s.eW>18&&eN1==1&&eN2==1&&s.Theta()>0.05&&s.Theta()<0.5");
  av.eSaveCouples = cenv.GetValue("fedra.align.SaveCouples" , 1);

  EdbPattern p1,p2;
  ReadPatCPnopar( p1, id1, cut, do_erase);
  ReadPatCPnopar( p2, id2, cut, do_erase);
  if(aff) { aff->Print(); p1.Transform(aff);}
 
  TString dataout;  MakeAffName(dataout,id1,id2,"al.root");
  av.InitOutputFile( dataout );
  av.Align( p1, p2 , dz);
  av.CloseOutputFile();
  UpdateAFFPar( id1, id2, av.eCorrL[0], aff );
  return npat;
}

//-------------------------------------------------------------------
int EdbScanProc::AlignNewLayernopar(EdbID id1, EdbID id2, TEnv &cenv, EdbPlateP &plate2, EdbAffine2D *aff, float dz)
{
  // Align 2 patterns. All necessary information should be in the envfile
  // Convension about Z(setted while process): the z of id2 is 0, the z of id1 is (-deltaZ) where
  // deltaZ readed from aff.par file in a way that pattern of id1 projected 
  // to deltaZ correspond to pattern of id2

  int npat=0;

  int nx = cenv.GetValue("fedra.align.NX", 1);
  int ny = cenv.GetValue("fedra.align.NY", 1);
  int nlayers = nx * ny;
  
  EdbPlateAlignment av[nlayers];
  for (int ialign = 0; ialign < nlayers; ialign ++){
   av[ialign].eOffsetMax =   cenv.GetValue("fedra.align.OffsetMax"   , 500. );
   av[ialign].SetSigma(      cenv.GetValue("fedra.align.SigmaR"      , 13.  ), 
       	            cenv.GetValue("fedra.align.SigmaT"      , 0.008) );
   av[ialign].eDoFine      = cenv.GetValue("fedra.align.DoFine"      , 1);
   av[ialign].eDZ          = cenv.GetValue("fedra.align.DZ"          , 120.);
   av[ialign].eDPHI        = cenv.GetValue("fedra.align.DPHI"        , 0.008 );
   av[ialign].eSaveCouples = cenv.GetValue("fedra.align.SaveCouples" , 1);
  }
  bool  do_erase  = cenv.GetValue("fedra.align.erase"    , false );
  const char *cut = cenv.GetValue("fedra.readCPcut"         , "eCHI2P<2.5&&s.eW>18&&eN1==1&&eN2==1&&s.Theta()>0.05&&s.Theta()<0.5");

  EdbPattern p1,p2;
  ReadPatCPnopar( p1, id1, cut, do_erase);
  ReadPatCPnopar( p2, id2, cut, do_erase);

  TString mapout;  MakeAffName(mapout,id1,id2,"map.root");
  TFile *corrmapfile = new TFile(mapout.Data(),"RECREATE");

  EdbCorrectionMap *corrmap = new EdbCorrectionMap();
  EdbCorrectionMap *corrmap_diff = new EdbCorrectionMap();
  //corrmap->Init(nx,p1.Xmin(),p1.Xmax(),ny,p1.Ymin(),p1.Ymax()); //If i need to get limits from patterns
  corrmap->Init(nx,0,120000,ny,0,100000); 
  corrmap_diff->Init(nx,0,120000,ny,0,100000); 

  if(aff) { aff->Print(); p1.Transform(aff);}
  //applying current local correction to plate 2
  if(plate2.Map().Ncell()>0) {
    int nseg = p2.N();
    Log(2,"AlignSetNewNoLayerPar","Apply correction to map with cells %i",plate2.Map().Ncell());
    for(int j=0; j<nseg; j++) {
            EdbSegP *s = p2.GetSegment(j);
            plate2.CorrectSegLocal(*s);
          }
    corrmap->ApplyCorrections(plate2.Map());
  }
  //instead of doing only once, we do it in a loop
  for (int ilayer = 0; ilayer < nlayers; ilayer++){
    Log(2,"AlignSetNewNoLayerpar","Starting alignment for layer number %i",ilayer);
    //extracting patterns form layer, X(ilayer) gives the position of middle point
    float min[5] = {-1000,-1000,-1000,-1000,-1000}; //X,Y,TX,TY,W
    float max[5] = {1000,1000,1000,1000,1000};

    min[0] = corrmap->X(corrmap->IX(ilayer)) - (corrmap->Xbin()/2.);
    max[0] = corrmap->X(corrmap->IX(ilayer))+ (corrmap->Xbin()/2.);
    min[1] = corrmap->Y(corrmap->IY(ilayer))- (corrmap->Ybin()/2.);
    max[1] = corrmap->Y(corrmap->IY(ilayer)) + (corrmap->Ybin()/2.);

    //corrmap->GetLayer(ilayer)->SetAffXY(aff->A11(),aff->A12(),aff->A21(),aff->A22(),aff->B1(),aff->B2()); //we do not include the global affine transformations in the map

    EdbPattern *p1_layer = p1.ExtractSubPattern(min,max);
    EdbPattern *p2_layer = p2.ExtractSubPattern(min,max);

    TString dataout;  MakeAffName(dataout,id1,id2,Form("al%i.root",ilayer));
    av[ilayer].InitOutputFile( dataout );
    av[ilayer].Align( *p1_layer, *p2_layer , dz);
    av[ilayer].CloseOutputFile();
    //saving number of coincidences and final transformations to map 
    corrmap->GetLayer(ilayer)->ApplyCorrections(1.,av[ilayer].eCorrL[0].Zcorr() - dz, *(av[ilayer].eCorrL[0].GetAffineXY()), *(av[ilayer].eCorrL[0].GetAffineTXTY()));
    corrmap_diff->GetLayer(ilayer)->ApplyCorrections(1.,av[ilayer].eCorrL[0].Zcorr()-dz, *(av[ilayer].eCorrL[0].GetAffineXY()), *(av[ilayer].eCorrL[0].GetAffineTXTY()));
    corrmap->GetLayer(ilayer)->SetNcp(av[ilayer].eNcoins);
    corrmap_diff->GetLayer(ilayer)->SetNcp(av[ilayer].eNcoins);
   }
  //UpdateAFFPar( id1, id2, av[0].eCorrL[0], aff );
  corrmapfile->cd();
  corrmap->Write("corrmap");
  corrmap_diff->Write("corrmap_diff");
  corrmapfile->Close();

  return npat;
}

//----------------------------------------------------------------
void EdbScanProc::CheckSetQualityRaw( EdbID idss )
{
  EdbScanSet *ss = ReadScanSet(idss);
  int n = ss->eIDS.GetSize();  if(n<1) return;
  for(int i=0; i<n; i++) {
    EdbID *id = (EdbID *)(ss->eIDS.At(i));
    TString name;
    MakeFileName(name,*id,"raw.root");
    EdbRunAccess a(name.Data()); 
    a.InitRun();
    a.CheckRunLine();
  }
  delete ss;
}

//----------------------------------------------------------------
void EdbScanProc::AlignSet(EdbID id, int npre, int nfull, const char *opt  )
{
  EdbScanSet *ss = ReadScanSet(id);
  if(!ss) { Log(1,"AlignSet","Error! set for %s do not found",id.AsString()); return; }
  AlignSet(*ss, npre,nfull, opt);
  ss->eB.ResetAff();
  WriteScanSet(id,*ss);
  UpdateSetWithAff(id,id);
}

//----------------------------------------------------------------
int EdbScanProc::AlignSet(EdbScanSet &sc, int npre, int nfull, const char *opt )
{
  if(sc.eIDS.GetSize()<2) return 0;
  int n=0;
  EdbID *id1,*id2;
  for(int i=0; i<sc.eIDS.GetSize()-1; i++) {
    id1 = (EdbID *)(sc.eIDS.At(i));
    id2 = (EdbID *)(sc.eIDS.At(i+1));
    int id14[4]; id1->Get(id14);
    int id24[4]; id2->Get(id24);
    n += AlignAll(id14, id24, npre, nfull, opt);
  }
  return n;
}

//----------------------------------------------------------------
int EdbScanProc::TrackSetBT(EdbScanSet &sc, TEnv &cenv)
{
  EdbScanCond cond;
  TCut  c        =    cenv.GetValue("fedra.readCPcut"     , "eCHI2P<2.5&&s.eW>13&&eN1==1&&eN2==1&&s1.eFlag>=0&&s2.eFlag>=0");
  int   nsegmin  =    cenv.GetValue("fedra.track.nsegmin"  , 2 );
  int   ngapmax  =    cenv.GetValue("fedra.track.ngapmax"  , 4 );
  float probmin  =    cenv.GetValue("fedra.track.probmin"  , 0.01 );
  float momentum =    cenv.GetValue("fedra.track.momentum" , 2 );
  float mass     =    cenv.GetValue("fedra.track.mass"     , 0.14 );
  bool  do_erase =    cenv.GetValue("fedra.track.erase"    , false );
  cond.SetSigma0(     cenv.GetValue("fedra.track.Sigma0"         , "3 3 0.005 0.005") );
  cond.SetPulsRamp0(  cenv.GetValue("fedra.track.PulsRamp0"      , "15 20") );
  cond.SetPulsRamp04( cenv.GetValue("fedra.track.PulsRamp04"     , "15 20") );
  cond.SetDegrad(     cenv.GetValue("fedra.track.Degrad"         , 4) );
  cond.SetRadX0(      cenv.GetValue("fedra.track.RadX0"          , 5810.) );

  if(sc.eIDS.GetSize()<2) return 0;
  int n=0;
  EdbID *id=0;
  for(int i=0; i<sc.eIDS.GetSize(); i++) {
    id = (EdbID *)(sc.eIDS.At(i));
    MakeInPar(*id,"tracking");
  }

  EdbDataProc dproc;
  EdbPVRec *ali = dproc.PVR();
  ali->SetScanCond( &cond );
  ReadScanSetCP(sc,  *ali, c, do_erase);
  ali->Print();

  n = dproc.LinkTracksWithFlag( ali, momentum, probmin, nsegmin, ngapmax, 0 );
  ali->FitTracks( momentum, mass );

  EdbID ido( *((EdbID *)(sc.eIDS.At(0))) );
  ido.ePlate=0;
  TString name; 
  MakeFileName(name,ido,"trk.root",false);
  dproc.MakeTracksTree(ali,name);

  return n;
}

//----------------------------------------------------------------
int EdbScanProc::ReadTracksTree(EdbID id, EdbPVRec &ali, TCut cut)
{
  TString name; 
  MakeFileName(name,id,"trk.root",false);
  return ReadTracksTree(name.Data(), ali, cut);
}

//----------------------------------------------------------------
int EdbScanProc::ReadTracksTree(const char *name, EdbPVRec &ali, TCut cut)
{
  int n=0;
  EdbDataProc dproc;
  n = dproc.ReadTracksTree(ali, name, cut);
  return n;
}


//----------------------------------------------------------------
int EdbScanProc::LinkSet(EdbScanSet &sc, int npre, int nfull, int correct_ang)
{
  int n=0;
  EdbID *id;
  for(int i=0; i<sc.eIDS.GetSize(); i++) {
    id = (EdbID *)(sc.eIDS.At(i));
    n += LinkRunAll(*id,npre,nfull,correct_ang);
  }
  return n;
}

//----------------------------------------------------------------
void EdbScanProc::CheckFiles( EdbScanSet &sc, const char *suffix )
{
  EdbID *id;
  for(int i=0; i<sc.eIDS.GetSize(); i++)  {
    id = (EdbID*)sc.eIDS.At(i);
    TString str;
    MakeFileName(str,*id,suffix);
    gSystem->Exec(Form("ls -l %s",str.Data()));
  }
}

//----------------------------------------------------------------
void EdbScanProc::MakeScannedIDList( EdbID id0, EdbScanSet &sc, int pl_from, int pl_to, const char *suffix )
{
  // check in all directories from pl_from to pl_to for a files with versions defined in id0 and the given suffix
  // if found - add EdbID the scan set
  sc.eIDS.Delete();
  EdbID id(id0);
  int step = pl_from>pl_to? -1:1;
  for(int i=pl_from; i != pl_to+step; i+=step)  {
    id.ePlate=i;
    TString str;
    MakeFileName(str,id,suffix);
    if( gSystem->AccessPathName(str, kReadPermission) )  
      {   Log(3,"MakeScannedIDList","not found %s",str.Data()); continue; }
    else  Log(3,"MakeScannedIDList","use file %s",str.Data());
    sc.eIDS.Add(new EdbID(id));
  }
}

//----------------------------------------------------------------
EdbMask *EdbScanProc::ReadEraseMask(EdbID id)
{
  TString str;
  MakeFileName(str,id,"er.root");
  Log(3,"EdbScanProc::ReadEraseMask"," %s",str.Data());
  if( gSystem->AccessPathName(str, kReadPermission) ) return 0;  //can not access file!
  TFile f(str.Data());
  if(!f.IsOpen()) return 0;

  EdbMask *m = 0;
  f.GetObject("mask",m);
  f.Close();
  return m;
}

//----------------------------------------------------------------
int EdbScanProc::ReadPatCPnopar(EdbPattern &pat, EdbID id, TCut cut, bool do_erase, bool read_mt)
{
  TString cpfile;
  MakeFileName(cpfile,id,"cp.root");
  EdbMask* mask = 0;
  if(do_erase) mask = ReadEraseMask(id);
  return ReadPatCPnopar(pat,cpfile.Data(),cut, mask, read_mt );
}

//----------------------------------------------------------------
int EdbScanProc::ReadPatCPnopar(EdbPattern &pat, const char *cpfile, TCut cut, EdbMask *mask, bool read_mt)
{
  EdbCouplesTree ect;
  if(!ect.InitCouplesTree("couples",cpfile,"READ")) return 0;
  ect.eCut=cut;
  ect.eEraseMask = mask;

  int nread=0;
  if(!read_mt) nread = ect.GetCPData( &pat );
  else {             // add microtracks as "digits" to the basetrack
    EdbPattern p1, p2;
    nread = ect.GetCPData( &pat,&p1,&p2,0 );
    for(int i=0; i<pat.N(); i++) {
      EdbSegP *s = pat.GetSegment(i);
      s->addEMULDigit( new EdbSegP(*(p1.GetSegment(i))) );
      s->addEMULDigit( new EdbSegP(*(p2.GetSegment(i))) );
      s->EMULDigitArray()->SetOwner();
    }
  }
  ect.Close();
  return nread;
}

//----------------------------------------------------------------
void EdbScanProc::PrepareVolumesPred( int idIN[4], EdbPattern &points,
				      int before, int after, int pmin, int pmax, EdbScanSet *sc )
{
  // Create prediction patterns for the stopping points in the segments of the pattern "points".
  // For each point should be correctly defined  (eID ePID eX eY eSX eSY).
  // Input: idIN   - identifier for predictions to be created
  //        points - stopping points in form of the segments (last found segment)
  //        before,after - number of plates before and after stopping point
  //        pmin,pmax    - brick limits (normally 1-57)
  //        sc           - EdbScanSet with affine transformations
  // Output: x.x.x.x.pred.root for all necessary plates

  TIndexCell cell;
  Long_t  v[2];  // pl,id
  EdbSegP *s=0;
  for(int i=0; i<points.N(); i++) {
    s = points.GetSegment(i);
    for(int ip=Max(s->PID()-before,pmin); ip<=Min(s->PID()+after,pmax); ip++) {
      v[0] = (Long_t)(ip);
      v[1] = (Long_t)(i);
      cell.Add(2,v);
    }
  }

  cell.Sort();
  int count=0;
  int plate;
  for(int ip=0; ip<cell.GetEntriesFast(); ip++) {
    plate = cell.At(ip)->Value();                             // current prediction plate
    EdbPattern pat;
    for(int iv=0; iv<cell.At(ip)->GetEntriesFast(); iv++) {
      s = points.GetSegment(cell.At(ip)->At(iv)->Value());    // s.PID() - stopping plate
      pat.AddSegment(*s);
      pat.GetSegment(iv)->SetPID(plate);
      if(sc) {
	EdbAffine2D aff;
	if(sc->GetAffP2P(s->PID(), plate, aff) ) pat.GetSegment(iv)->Transform(&aff);
      }
    }
    int id[4]={idIN[0],plate,idIN[2],idIN[3]};
    WritePred(pat,id);
    count += pat.N();
  }
  
  LogPrint(idIN[0],2,"PrepareVolumesPred","%d.%d.%d.%d: for %d  volumes generated %d predictions with settings: before=%d after=%d pmin=%d pmax=%d\n",
	   idIN[0],idIN[1],idIN[2],idIN[3],points.N(), count, before,after,pmin,pmax );
}

//----------------------------------------------------------------
bool EdbScanProc::FlashRawDir(EdbScanClient &scan, int id[4])
{
  // move all rwc and rwd files from the raw scanning directory into the new subdir

  char str[256];
  TDatime dt;
  sprintf(str,"%s/rw_%u",scan.GetRawDirClient(),dt.Get());
  LogPrint(id[0],2,"FlashRawDir","%d.%d.%d.%d: move all into %s", id[0],id[1],id[2],id[3],str);
  if(!gSystem->OpenDirectory(str))   
    if( gSystem->MakeDirectory(str) == -1) 
      {
	LogPrint(id[0],2,"FlashRawDir","WARNING! %d.%d.%d.%d: FAILED creating directory %s", 
		 id[0],id[1],id[2],id[3],str);
	return false;
      }
  char str2[256];

#ifdef WIN32
  sprintf(str2,"ren %s/raw.* %s",scan.GetRawDirClient(),str);
#else
  sprintf(str2,"mv %s/raw.* %s",scan.GetRawDirClient(),str);
#endif

  gSystem->Exec(str2);
  return true;
}

//----------------------------------------------------------------
int EdbScanProc::LoadPlate(EdbScanClient &scan, int id[4], int attempts)
{
  int status=0;
  FlashRawDir(scan,id);
  TString map;
  if(!GetMap(id[0],map)) {
    LogPrint(id[0],1,"LoadPlate","ERROR: map file does not exist! Stop here. *** %d.%d.%d  status = %d ***",
	     id[0],id[1],id[2],status);
    return status;
  }
  status= scan.LoadPlate(id[0],id[1],map.Data(),attempts);
  LogPrint(id[0],1,"LoadPlate","******************** %d.%d.%d  status = %d ************************************",
	   id[0],id[1],id[2],status);
  return status;
}

//------------------------------------------------------------------------------------------
int EdbScanProc::RemoveDublets( EdbPattern &pin, EdbPattern &pout, int brick )
{
  // input:  pin  - predictions pattern
  // output: pout - predictions pattern with dublets removed
  if(pin.N()<2)     return 0;
  float r,rt;
  float rmin = 0.1;    // [micron] TODO - pass as parameter?
  float rtmin= 0.0001; // [rad]    TODO - pass as parameter?
  OptimizeScanPath(pin, pout,brick);
  EdbSegP *s=0,*s1=0;
  int n= pout.N();
  for(int i=n-1; i>0; i--) {
    s  = pout.GetSegment(i);
    s1 = pout.GetSegment(i-1);
    r  = Sqrt((s->X()-s1->X())*(s->X()-s1->X()) + (s->Y()-s1->Y())*(s->Y()-s1->Y()));
    if( r>rmin   )   continue;
    rt = Sqrt((s->TX()-s1->TX())*(s->TX()-s1->TX()) + (s->TY()-s1->TY())*(s->TY()-s1->TY()));
    if( rt>rtmin )   continue;
    pout.GetSegments()->RemoveAt(i);
  }
  pout.GetSegments()->Compress();
  LogPrint(brick,2,"RemoveDublets","%d segments before -> %d segments after: %d dublets removed", n, pout.N(), n-pout.N());
  return n-pout.N();
}

//------------------------------------------------------------------------------------------
void EdbScanProc::OptimizeScanPath(EdbPattern &pin, EdbPattern &pout, int brick)
{
  // input:  pin  - predictions pattern
  // output: pout - predictions pattern with optimized path (should be empty at the beginning)

  int n = pin.N();
  if(pout.N()) pout.GetSegments()->Delete();
  if(n>2) {
    EdbSegP *s;
    TIndexCell cell;
    float xmin = pin.Xmin()-0.000001, xmax = pin.Xmax()+0.000001;
    float ymin = pin.Ymin()-0.000001, ymax = pin.Ymax()+0.000001;
    float eps  = Sqrt(3.*(xmax-xmin)*(ymax-ymin)/n);
    float binx = (xmax-xmin)/((int)((xmax-xmin)/eps));
    Long_t  v[3];  // x,y,i
    for(int i=0; i<n; i++ ) {
      s = pin.GetSegment(i);
      v[0] = (Long_t)((s->X()-xmin)/binx);
      v[1] = (Long_t)((s->Y()-ymin));
      if(v[0]%2 != 0) v[1] = -v[1];         // serpentina along y
      v[2] = i;
      cell.Add(3,v);
    }
    cell.Sort();
    for(int ix=0; ix<cell.GetEntriesFast(); ix++)
      for(int iy=0; iy<cell.At(ix)->GetEntriesFast(); iy++)
	for(int ii=0; ii<cell.At(ix)->At(iy)->GetEntriesFast(); ii++) {
	  s = pin.GetSegment(cell.At(ix)->At(iy)->At(ii)->Value());
	  pout.AddSegment(*s);
	}
  }
  if( pout.N()>0 )
    if( pin.SummaryPath() > pout.SummaryPath() ) {          // good optimization
      LogPrint(brick,2,"OptimizeScanPath","with %d  predictions: gain in path[mm] before/after = %.1f/%.1f = %.1f",
	       n, pin.SummaryPath()/1000,pout.SummaryPath()/1000, pin.SummaryPath()/pout.SummaryPath());
      return;
    }
  if(pout.N()) pout.GetSegments()->Delete();
  for(int i=0; i<n; i++) pout.AddSegment(*(pin.GetSegment(i)));
}

//----------------------------------------------------------------
int EdbScanProc::ConvertAreas(EdbScanClient &scan, int id[4], int flag, const char *opt)
{
  // can be called separately in case of missed conversion
  EdbPattern pred;
  ReadPred(pred, id, flag);
  EdbRun *run = InitRun(id);
  if(!run) return 0;
  int scanned = scan.ConvertAreas(id,pred,*run,opt);
  LogPrint(id[0],1,"ConvertAreas","%d.%d.%d.%d  with %d predictions with flag %d; %d views stored", 
	   id[0],id[1],id[2],id[3],pred.N(),flag,run->GetEntries());
  run->Close();
  delete run;
  return scanned;
}

//----------------------------------------------------------------
int EdbScanProc::ScanAreas(EdbScanClient::ScanType st, EdbScanClient &scan, int id[4], int flag, const char *opt)
{
  EdbPattern pred;
  ReadPred(pred, id, flag);
  return ScanAreas(st, scan, pred, id, opt);
}

//----------------------------------------------------------------
int EdbScanProc::ScanAreas(EdbScanClient::ScanType st, EdbScanClient &scan, EdbPattern &pred, int id[4], const char *opt) // NEW!!!
{
  LogPrint(id[0],1,"ScanAreas","%d.%d.%d.%d  with %d predictions", id[0],id[1],id[2],id[3],pred.N());
  EdbPattern predopt;
  OptimizeScanPath(pred,predopt,id[0]);

  bool createRun = !scan.ServerCreatesRootFile(); // if is created by server side - no need to create it
  char runname[512];
  char runNameSrv[512];
  
  EdbRun *run = InitRun(id, runname, runNameSrv, createRun);
  if(!run && createRun) return 0;

  if(!createRun && eProcDirServer.Length())//if ProcDir not set - server cant create proper files
    scan.SetProcTgtServer(runNameSrv);
  scan.SetProcPthServer(eProcDirServer.Data());
    
  int scanned = scan.ScanAreas(st, id,predopt,run,opt);
  LogPrint(id[0],1,"ScanAreas","%d.%d.%d.%d  %d predictions scanned; run with %d views stored", id[0],id[1],id[2],id[3],scanned, (createRun)? run->GetEntries(): (-1) );
	
  if(run){
    run->Close();
    delete run;
  }

  if(!createRun){
    if(!scan.ServerCreatesTarget()){//server creates *.root in tmp dir move server-side created file in target loaction (<*>/brick/plate/*.*.*.*.raw.root)
      TString serverCreatedRunName = scan.GetServerCreatedRunName();
      if(serverCreatedRunName.Length()==0){
	printf("Server had to create run but we don't know where it is.\n");
	return 0;
      }
      WaitFileReady(serverCreatedRunName.Data());
      char str[1024];
  #ifdef WIN32
      sprintf(str,"move /F %s %s", serverCreatedRunName.Data(), runname);
  #else
      sprintf(str,"mv -f %s %s", serverCreatedRunName.Data(), runname);
  #endif
      gSystem->Exec(str);
    }else{//Server creates *.root on place
      WaitFileReady(runname);
    }
  }
  return scanned;
}

//----------------------------------------------------------------
int EdbScanProc::CopyFile(int id1[4], int id2[4], const char *suffix, bool overwrite)
{
  // copy piece file from id1 to id2
  TString name1, name2;
  MakeFileName(name1,id1,suffix);
  MakeFileName(name2,id2,suffix);
  int status = gSystem->CopyFile(name1.Data(), name2.Data(), overwrite);
  LogPrint(id1[0],2,"CopyFile","status=%d from %s to %s", status, name1.Data(),name2.Data() );
  //sleep(10);
  if( gSystem->AccessPathName(name2, kReadPermission) ) return 0; //can not access file!
  return 1;
}

//----------------------------------------------------------------
int EdbScanProc::RemoveFile(EdbID id, const char *suffix)
{
  // remove file
  TString name;
  MakeFileName(name,id,suffix);
  char str[256];
  sprintf(str,"%s/file_to_remove",eProcDirClient.Data());
  int status = gSystem->Rename(name.Data(), str);
  LogPrint(id.eBrick,2,"RemoveFile","status=%d from %s to %s", status, name.Data(), str );
  return status;
}

//----------------------------------------------------------------
int EdbScanProc::CopyAFFPar(int id1c[4], int id2c[4], int id1p[4], int id2p[4], bool overwrite)
{
  // copy AFF/xc_yc.par to AFF/xp_yp.par
  TString name1, name2;
  MakeAffName(name1,id1c,id2c);
  MakeAffName(name2,id1p,id2p);
  LogPrint(id1c[0],2,"CopyAFFPar","from %s to %s", name1.Data(),name2.Data());
  return gSystem->CopyFile(name1.Data(), name2.Data(), overwrite);
}

//----------------------------------------------------------------
bool EdbScanProc::UpdateAFFPar( EdbID id1, EdbID id2, EdbLayer &l, EdbAffine2D *aff0 )
{
  // update the aff par file with the values defined in the EdbLayer

  TString parout;  MakeAffName(parout,id1,id2);
  LogPrint(id1.eBrick,2,"UpdateAFFPar","%s ", parout.Data());
  char card[80];

  TDatime t;
  sprintf(card,"\n## %s", t.AsSQLString());
  if(!AddParLine(parout.Data(),card)) return false;

  sprintf(card,"ZLAYER \t %d \t %f %f %f",l.ID(), l.Zcorr() ,0.,0. );
  if(!AddParLine(parout.Data(),card)) return false;

  EdbAffine2D &aff = *(l.GetAffineXY());
  if(aff0) aff.Transform(aff0);

  //sprintf(card,"AFFXY \t %d \t %f %f %f %f %f %f", l.ID(), 
   sprintf(card,"AFFXY \t %d \t %g %g %g %g %g %g", l.ID(), 
                  aff.A11(),aff.A12(),aff.A21(),aff.A22(),aff.B1(),aff.B2() );
  if(!AddParLine(parout.Data(),card)) return false;

  EdbAffine2D &afftxy = *(l.GetAffineTXTY());
  sprintf(card,"AFFTXTY \t %d \t %g %g %g %g %g %g", l.ID(), 
	  afftxy.A11(),afftxy.A12(),afftxy.A21(),afftxy.A22(),afftxy.B1(),afftxy.B2() );
  if(!AddParLine(parout.Data(),card)) return false;
  
  sprintf(card,"SHRINK \t %d \t %f",l.ID(), l.Shr() );  if(!AddParLine(parout.Data(),card)) return false;

  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::UpdatePlatePar( EdbID id, EdbLayer &l)
{
  // update the plate par file with the values defined in the EdbLayer

  TString parout;  MakeFileName(parout,id,"par");
  LogPrint(id.eBrick,2,"UpdatePar","%s ", parout.Data());
  char card[80];

  TDatime t;
  sprintf(card,"\n## %s", t.AsSQLString());             if(!AddParLine(parout.Data(),card)) return false;
  sprintf(card,"SHRINK \t %d \t %f",l.ID(), l.Shr() );  if(!AddParLine(parout.Data(),card)) return false;

  EdbAffine2D &aff = *(l.GetAffineXY());
  sprintf( card,"AFFXY \t %d %s", l.ID(), aff.AsString() );  if(!AddParLine(parout.Data(),card)) return false;

  aff = *(l.GetAffineTXTY());
  sprintf( card,"AFFTXTY \t %d %s", l.ID(), aff.AsString() );  if(!AddParLine(parout.Data(),card)) return false;

  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::SetAFF0(int id1[4], int id2[4])
{
  TString str;
  MakeAffName(str,id1,id2);
  char card[64];
  sprintf(card,"AFFXY 0 1. 0. 0. 1. 0. 0.");
  LogPrint(id1[0],2,"SetAFF0","%s as %s", str.Data(),card);
  return AddParLine(str.Data(),card);
}

//----------------------------------------------------------------
bool EdbScanProc::SetAFFDZ(int id1[4], int id2[4], float dz)
{
  TString str;
  MakeAffName(str,id1,id2);
  char card[64];
  sprintf(card,"ZLAYER 0 %f 0 0",dz);
  LogPrint(id1[0],2,"SetAFFDZ","%s as %s", str.Data(),card);
  return AddParLine(str.Data(),card);
}

//----------------------------------------------------------------
int EdbScanProc::LinkRunAll(int id[4], int npre, int nfull, int correct_ang )
{
  LogPrint(id[0],1,"LinkRunAll","%d.%d.%d.%d  %d prelinking + %d fullinking", 
	   id[0],id[1],id[2],id[3],npre,nfull);
  Int_t nc = 0;
  if (npre > 0) {
    MakeInPar(id, "prelinking");    // make input par file (x.x.x.x.in.par) for the current ID including the prelinking par file
    for (Int_t i = 0; i < npre; i++) {
      nc = LinkRun(id, 0);          // will be done (pre)linking and updated x.x.x.x.par file
      if (correct_ang) CorrectAngles(id);
    }
  }
  if (nfull > 0) {
    MakeInPar(id, "fulllinking");   // make input par file including the fulllinking par file
    for (Int_t i = 0; i < nfull; i++) {
      //SafeDelete(gDIFF);
      //TFile f("diff.root","RECREATE");
      //gDIFF = new TNtuple("diff","diff","x1:y1:tx1:ty1:w1:x2:y2:tx2:ty2:w2:z:aid10:aid11:aid20:aid21");
      
      nc = LinkRun(id,1);         // will be done (full)linking and DO NOT updated x.x.x.x.par file

      //gDIFF->AutoSave();
      //f.Close();
      //SafeDelete(gDIFF);
    }
  }
  LogPrint(id[0],1,"LinkRunAll","%d couples stored", nc);
  return nc;
}

//----------------------------------------------------------------
bool EdbScanProc::ProjectFound(int id1[4],int id2[4])
{
  //take xp.xp.xp.xp.found.root and produce yp.yp.yp.yp.pred.root using the AFF/xp_yp.par
  LogPrint(id1[0],2,"ProjectFound","from %d.%d.%d.%d to %d.%d.%d.%d", 
	   id1[0],id1[1],id1[2],id1[3],id2[0],id2[1],id2[2],id2[3]);
  EdbPattern pat;
  ReadFound(pat,id1);
  //pat.SetZ(0);
  //pat.SetSegmentsZ();
  for (Int_t i = 0; i < pat.N(); i++) {
    pat.GetSegment(i)->SetErrors(50., 50., 0., .1, .1);  // area to be scanned for this prediction
  }
  ApplyAffZ(pat,id1,id2);
  if (!WritePred(pat,id2)) return false;
  WritePatTXT(pat,id2,"man.pred.txt");
  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::CorrectPredWithFound(int id1[4],int id2[4], const char *opt, int patmin)
{
  // take p1.found.root, apply AFF/p1_p2.par, read p2.found.root, align and update AFF/p1_p2.par

  EdbPattern p1;   ReadFound(p1,id1);
  EdbPattern p2;   ReadFound(p2,id2);

  if(p1.N()<2||p2.N()<2) {
    LogPrint(id1[0],1,"CorrectPredWithFound","ERROR: correction is impossible - too small pattern: %d", p1.N(), p2.N());
    return false;
  }
  else   if(p1.N()<patmin || p2.N()<patmin) 
    LogPrint(id1[0],1,"CorrectPredWithFound","WARNING: unreliable correction - pattern is too small: %d < %d", p1.N(), p2.N(),patmin);
  
  EdbAffine2D aff;  float dz; 
  if(!GetAffZ(aff, dz, id1,id2))  return false;
  p1.Transform(&aff);
  p1.SetZ(-dz);   p1.SetSegmentsZ();
  p2.SetZ(0);     p2.SetSegmentsZ();

  MakeInPar(id2,"fullalignment");
  EdbDataPiece piece2;
  InitPiece(piece2, id2);

  EdbPVRec ali;
  ali.AddPattern(&p1);
  ali.AddPattern(&p2);
  EdbScanCond *cond = piece2.GetCond(0);
  cond->SetChi2Mode(3);
  ali.SetScanCond( cond );
  ali.SetPatternsID();
  ali.SetSegmentsErrors();
  ali.SetCouplesAll();
  ali.SetChi2Max(cond->Chi2PMax());
  ali.SetOffsetsMax(cond->OffX(),cond->OffY());

  ali.Align(2);
  //ali.Align(0);
  int nal = ali.GetCouple(0)->Ncouples();

  if(nal<patmin) {
    LogPrint(id1[0],1,"CorrectAffWithPred","WARNING: pattern is too small: %d < %d: do not update par file!", nal, patmin);
    return false;
  }

  TString parfileOUT;
  MakeAffName(parfileOUT,id1,id2);
  piece2.eFileNamePar = parfileOUT;
  ali.GetPattern(0)->GetKeep(aff);
  piece2.UpdateAffPar(0,aff);

  if( strstr(opt,"-z")) {
    EdbDataProc proc;
    proc.LinkTracksWithFlag( &ali, 10., 0.05, 2, 3, 0 );
    printf("befire corr z: z1=%f z2=%f \n",ali.GetPattern(0)->Z(), ali.GetPattern(1)->Z());
    ali.FineCorrZnew();
    printf("after corr z: z1=%f z2=%f \n",ali.GetPattern(0)->Z(), ali.GetPattern(1)->Z());
    piece2.UpdateZPar(0,-ali.GetPattern(0)->Z());
  }
 
  LogPrint(id1[0],1,"CorrectPredWithFound","from %d.%d.%d.%d to %d.%d.%d.%d: used %d (out of %d predictions and %d found) for correction", 
	   id1[0],id1[1],id1[2],id1[3],id2[0],id2[1],id2[2],id2[3],  nal, p1.N(), p2.N() );
  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::CorrectAffWithPred(int id1[4],int id2[4], const char *opt, int patmin, const char *parfile)
{
  // take p1.found.root, apply AFF/p1_p2.par, align to p2 and update AFF/p1_p2.par
  EdbPattern pat;
  ReadFound(pat,id1);
  if(pat.N()<patmin)
    LogPrint(id1[0],1,"CorrectAffWithPred","WARNING: unreliable correction - pattern is too small: %d < %d", pat.N(),patmin);
  else if(pat.N()<2) {
    LogPrint(id1[0],1,"CorrectAffWithPred","ERROR: correction is impossible - too small pattern: %d", pat.N());
    return false;
  }

  EdbAffine2D aff;
  float dz; 
  if(!GetAffZ(aff, dz, id1,id2))  return false;
  pat.Transform(&aff);
  pat.SetZ(-dz);
  pat.SetSegmentsZ();

  MakeInPar(id2,parfile);
  EdbPattern p2;
  EdbDataPiece piece2;
  InitPiece(piece2, id2);
  ReadPiece(piece2, p2);
  p2.SetZ(0);
  p2.SetSegmentsZ();

  EdbPVRec ali;
  ali.AddPattern(&pat);
  ali.AddPattern(&p2);
  EdbScanCond *cond = piece2.GetCond(0);
  cond->SetChi2Mode(3);
  ali.SetScanCond( cond );
  ali.SetPatternsID();
  ali.SetSegmentsErrors();
  ali.SetCouplesAll();
  ali.SetChi2Max(cond->Chi2PMax());
  ali.SetOffsetsMax(cond->OffX(),cond->OffY());

  ali.Align(2);
  //ali.Align(0);
  int nal = ali.GetCouple(0)->Ncouples();

  if(nal<patmin) {
    LogPrint(id1[0],1,"CorrectAffWithPred","WARNING: pattern is too small: %d < %d: do not update par file!", nal, patmin);
    return false;
  }

  TString parfileOUT;
  MakeAffName(parfileOUT,id1,id2);
  piece2.eFileNamePar = parfileOUT;
  ali.GetPattern(0)->GetKeep(aff);
  piece2.UpdateAffPar(0,aff);
  if( strstr(opt,"-z")) {
    EdbDataProc proc;
    proc.LinkTracksWithFlag( &ali, 10., 0.05, 2, 3, 0 );
    ali.FineCorrZnew();
    piece2.UpdateZPar(0,-ali.GetPattern(0)->Z());
  }
 
  LogPrint(id1[0],1,"CorrectAffWithPred","from %d.%d.%d.%d to %d.%d.%d.%d: used %d (out of %d predictions) for correction", 
	   id1[0],id1[1],id1[2],id1[3],id2[0],id2[1],id2[2],id2[3],  nal, pat.N() );
  return true;
}

//----------------------------------------------------------------
int EdbScanProc::FindPredictions(int id[4], int flag, int maxholes)
{
  //find predictions of yp.yp.yp.yp.pred.root in yp.yp.yp.yp.cp.root and produce yp.yp.yp.yp.found.root
  // flag: -1 - all predictions
  //        0 - only found in the previous plate (no holes)
  //        1 - 1 holes before
  //        2 - 2 holes before
  //        3 - 3 holes before

  int nfound=0;
  EdbPattern pred, found;
  ReadPred(pred,id, flag);
  nfound =  FindPredictions(pred, id,found, maxholes);
  WriteFound(found,id);

  //EdbTestAl ta;
  //ta.HDistance(pred,found);
  //ta.FillTree(0);
  //EdbAffine2D aff;
  //ta.MakeTrans(aff,0,"(dx*dx+dy*dy)>0.00001&&abs(dx)<100&&abs(dy)<100");

  return nfound;
}

//----------------------------------------------------------------
int EdbScanProc::TestAl(int id1[4], int id2[4])
{
  MakeInPar(id1,"prealignment");
  MakeInPar(id2,"prealignment");
  TString name;
  MakeFileName(name,id1,"in.par");
  TString parfileOUT;
  MakeAffName(parfileOUT,id1,id2);
  parfileOUT.Prepend("INCLUDE ");
  AddParLine(name.Data(),	parfileOUT.Data());
  
  EdbPattern p1,p2;

  ReadPatCP(p2,id2);
  p2.SetZ(0);
  p2.SetSegmentsZ();

  ReadPatCP(p1,id1);
  EdbAffine2D aff1;
  float dz;
  GetAffZ(aff1,dz,id1,id2);
  p1.SetZ(-dz);
  p1.SetSegmentsZ();

  return TestAl(p1,p2);
}

//----------------------------------------------------------------
int EdbScanProc::TestAl(const char *cpfile1, const char *cpfile2, TCut &cut, float dz, EdbAffine2D *aff)
{
  EdbPattern p1,p2;
  ReadPatCPnopar(p2, cpfile2, cut);
  p2.SetZ(0);  p2.SetSegmentsZ();
  ReadPatCPnopar(p1, cpfile1, cut);
  if(aff) p1.Transform(aff);
  p1.SetZ(-dz);  p1.SetSegmentsZ();
  return TestAl(p1,p2);
}

//----------------------------------------------------------------
int EdbScanProc::TestAl(EdbPattern &p1, EdbPattern &p2)
{
  Log(2,"EdbScanProc::TestAl","align patterns %d and %d", p1.N(), p2.N());
  EdbTestAl ta;
  TFile ftree("testal.root","RECREATE");
  ta.eBinTree = new TNtuple("bintree","bin tree","dz:phi:dx:dy:bin");
  ta.eITMAX=50;  // default value
  ta.HDistance(p1,p2);

  float bin[4]={250,250, 250,0.001};                  // default values for normal alignment (expected dz=1300)
  ta.eDmin[0]=-25000; ta.eDmin[1]=-25000; ta.eDmin[2]=  -25; ta.eDmin[3]=-0.015;
  ta.eDmax[0]= 25000; ta.eDmax[1]= 25000; ta.eDmax[2]=   25; ta.eDmax[3]= 0.015;

  FILE *f = fopen("testal.par","r");
  if(f) {
    for(int i=0; i<4; i++) {
      float min=0,max=0,b=0;
      if(!(fscanf(f,"%f %f %f",&min,&max,&b )==3))  {Log(1,"EdbScanProc::TestAl","ERROR: read from testal.par"); return 0;}
      else { 
	ta.eDmin[i]=min; ta.eDmax[i]=max;  bin[i]=b; 
      }
    }
    fclose(f);
  }
  for(int i=0; i<4; i++) ta.eN[i] = (Int_t)((ta.eDmax[i]-ta.eDmin[i]-bin[i]/2.)/bin[i])+1;
  for(int i=0; i<4; i++) ta.eDmax[i] = ta.eDmin[i]+bin[i]*ta.eN[i];
  for(int i=0; i<4; i++) printf("%d \t%f %f %f %d\n",i, ta.eDmin[i],ta.eDmax[i],bin[i],ta.eN[i]);

  ta.CheckMaxBin();
  EdbAffine2D aff;
  aff.Rotate(-ta.eD0[3]);
  aff.ShiftX(ta.eD0[0]);
  aff.ShiftY(ta.eD0[1]);
  //if(gEDBDEBUGLEVEL>1) {
    //aff.Print();
    //ta.FillTree(-ta.eD0[2]);
  //}

  ta.eBinTree->Write("bintree");
  ftree.Close();

  return 0;
}
  
//----------------------------------------------------------------
void EdbScanProc::MergeSetSBT(EdbID id)
{
  EdbScanSet  *ss = ReadScanSet(id);
  if(ss) MergeSetSBT(id, *ss);
  delete ss;
}
  
//----------------------------------------------------------------
void EdbScanProc::MergeSetSBT(EdbID id, EdbScanSet &ss)
{
  // merge all sbt files for a given scan set into bxxxxx.x.x.x.sbt.root
  if(!CheckBrickDir(id)) return;
  int npl = ss.eIDS.GetSize();
  if(npl<1) return;
  TString filename;
  TChain allsbt("sbt");

  for( int i=0; i<npl; i++ ) {
    EdbID *idp  = ss.GetID(i);
    Log(2,"MergeSetSBT", "Adding plate %d...",idp->ePlate);
    MakeFileName(filename,*idp,"sbt.root");
    allsbt.Add(filename.Data());
  }
  MakeFileName(filename,id,"sbt.root",false);
  Log(2,"MergeSetSBT","Writing %s file...",filename.Data());
  allsbt.Merge(filename.Data());
}

//----------------------------------------------------------------
int EdbScanProc::WriteSBTrack(EdbTrackP &sb, int path, EdbID id)
{
  if(!CheckBrickDir(id)) return 0;
  TString name;
  MakeFileName(name,id,"sb.root",false);
  TFile f(name.Data(),"UPDATE");
  if(!f.IsOpen()) return 0;
  sb.Write(Form("sb_%d",path));
  f.Close();
  return 1;
}

//----------------------------------------------------------------
int EdbScanProc::WriteSBTracks(TObjArray &tracks, EdbID id)
{
  if(!CheckBrickDir(id)) return 0;
  TString name;
  MakeFileName(name,id,"sb.root",false);
  //TFile f(name.Data(),"UPDATE");
  TFile f(name.Data(),"RECREATE");
  if(!f.IsOpen()) return 0;
  tracks.Write("sbtracks",1);
  f.Close();
  LogPrint(id.eBrick,2,"WriteSBTracks","%d tracks are written into %s", tracks.GetEntriesFast(), name.Data());
  return 1;
}

//----------------------------------------------------------------
TObjArray *EdbScanProc::ReadSBTracks(EdbID id)
{
  if(!CheckBrickDir(id)) return 0;
  TString name;
  MakeFileName(name,id,"sb.root",false);
  TFile f(name.Data(),"READ");
  if(!f.IsOpen()) return 0;
  TObjArray *tracks = (TObjArray*)f.Get("sbtracks");
  f.Close();
  LogPrint(id.eBrick,2,"ReadSBTracks","%d tracks read from %s", tracks->GetEntriesFast(), name.Data());
  return tracks;
}

//----------------------------------------------------------------
int EdbScanProc::WriteScanSet(EdbID id, EdbScanSet &ss)
{
  if(!CheckBrickDir(id)) return 0;
  TString name;
  MakeFileName(name,id,"set.root",false);
  TFile f(name.Data(),"UPDATE");
  if(!f.IsOpen()) return 0;
  ss.Write("set");
  f.Purge(3);  f.Close();
  return 1;
}

//----------------------------------------------------------------
EdbScanSet *EdbScanProc::ReadScanSetGlobal( EdbID id, bool x_marks )
{
  // read scan set and use map.OL and map.XG files to obtain transformations in global OPERA RS
  EdbScanSet *ss = ReadScanSet(id);
  if(ss) {
    if(gEDBDEBUGLEVEL>2) ss->Print();
    EdbMarksSet msOL;
    if (x_marks) ReadMarksSet(msOL,id.eBrick,"map.LL",'_','L');
    else ReadMarksSet(msOL,id.eBrick,"map.OL",'_','S');
    if(gEDBDEBUGLEVEL>2) msOL.Print();
    float MINX = msOL.eXmin;
    float MAXX = msOL.eXmax;
    float MINY = msOL.eYmin;
    float MAXY = msOL.eYmax;
    EdbMarksSet msXG;
    ReadMarksSet(msXG,id.eBrick,"map.XG",'_');
    float ZEROX = msXG.eXmin;
    float ZEROY = msXG.eYmin;
    if(gEDBDEBUGLEVEL>2) msXG.Print();
   
    int n = ss->eIDS.GetSize();
    for(int i=0; i<n; i++) {
      EdbID *pid  = ss->GetID(i);
      EdbPlateP *plate = ss->GetPlate(pid->ePlate);
      float RX = 0.5*(MINX+MAXX) - ZEROX;
      float RY = 0.5*(MINY+MAXY) - ZEROY;
      EdbAffine2D *aff = plate->GetAffineXY();
      float MAPXX = aff->A11();
      float MAPXY = aff->A12();
      float MAPYX = aff->A21();
      float MAPYY = aff->A22();
      float MAPDX = aff->B1() - RX + aff->A11() * RX + aff->A12() * RY;
      float MAPDY = aff->B2() - RY + aff->A21() * RX + aff->A22() * RY;
      aff->Set(MAPXX,MAPXY,MAPYX,MAPYY,MAPDX,MAPDY);
    }
    if(gEDBDEBUGLEVEL>2) ss->Print();
  }
  return ss;
}

//----------------------------------------------------------------
EdbScanSet *EdbScanProc::ReadScanSet(EdbID id)
{
  TString name;
  MakeFileName(name,id,"set.root",false);
  TFile f(name.Data());
  Log(2,"ReadScanSet","from %s",name.Data());
  if(!f.IsOpen()) return 0;
  TObject *ss = f.Get("set");
  if(ss) return (EdbScanSet *)ss;
  else   return 0;
}

//----------------------------------------------------------------
int EdbScanProc::WritePatTXT(EdbPattern &pred, int id[4], const char *suffix, int flag)
{
  // write ascii predictions file as .../bXXXXXX/pYYY/a.a.a.a.suffix
  //        man      - for manual check by sysal
  //                 - overvise complete format
  EdbSegP *s=0;
  TString str;
  MakeFileName(str,id,suffix);
  FILE *f = fopen(str.Data(),"w");
  if(!f) { LogPrint(id[0],1,"WritePatTXT","ERROR! can not open file %s", str.Data()); return 0; }
  for(int i=0; i<pred.N(); i++) {
    s = pred.GetSegment(i);
    if (flag > -1 && flag != s->Flag())  continue;
    if( strcmp(suffix,"man") >=0 )
      fprintf(f,"%8d %11.2f %11.2f %8.4f %8.4f %d\n",
	      s->ID(),s->X(),s->Y(),s->TX(),s->TY(),s->Flag());
    else 
      fprintf(f,"%8d %11.2f %11.2f %8.4f %8.4f %9.2f %9.2f %8.4f %8.4f %d\n",
	      s->ID(),s->X(),s->Y(),s->TX(),s->TY(),s->SX(),s->SY(),s->STX(),s->STY(),s->Flag());
  }
  fclose(f);
  LogPrint(id[0],2,"WritePatTXT","%s with %d predictions with flag: %d", str.Data(),pred.N(), flag);
  return pred.N();
}

//----------------------------------------------------------------
int EdbScanProc::WriteSBcndTXT(int id[4], const char *suffix)
{
  // write ascii found candidates file as .../bXXXXXX/pYYY/a.a.a.a.suffix
  //                 - new id = 1000*id + 100* [BT:0; MT1:1; MT2:2]+ #
  TString sbtname;
  MakeFileName(sbtname,id,"sbt.root");
  
  TString outname;
  MakeFileName(outname,id,suffix);
  
  
  FILE *f = fopen(outname.Data(),"w");
  if(!f) { LogPrint(id[0],1,"WriteSBcndTXT","ERROR! can not open file %s", outname.Data()); return 0; }
  int num=0;
  
  
  EdbRunTracking t;
  TTree *st = t.InitSBtree(sbtname.Data(), "OPEN");
  int n = st->GetEntries();
  EdbSegP *s=0;
  for(int i =0; i<n; ++i){
    t.GetSBtreeEntry(i, *st);
    int pid = t.eNext.ID();
    for(int j = 0; j<t.eScnd.N(); ++j){
      s=t.eScnd.GetSegment(j);
      fprintf(f, "%d\t%.3lf\t%.3lf\t%.5lf\t%.5lf\t0\n", pid*1000+j, s->X(), s->Y(), s->TX(), s->TY());
      num++;
    }
    for(int j = 0; j<t.eS1cnd.N(); ++j){
      s=t.eS1cnd.GetSegment(j);
      fprintf(f, "%d\t%.3lf\t%.3lf\t%.5lf\t%.5lf\t10001\n", pid*1000+100+j, s->X(), s->Y(), s->TX(), s->TY());
      num++;
    }
    for(int j = 0; j<t.eS2cnd.N(); ++j){
      s=t.eS2cnd.GetSegment(j);
      fprintf(f, "%d\t%.3lf\t%.3lf\t%.5lf\t%.5lf\t10002\n", pid*1000+200+j, s->X(), s->Y(), s->TX(), s->TY());
      num++;
    }
  }
  t.CloseSBtree(st);
  
  
  fclose(f);
  LogPrint(id[0],2,"WriteSBcndTXT","%s with %d predictions", outname.Data(),num);
  return num;
}

//----------------------------------------------------------------
int EdbScanProc::ReadPatTXT(EdbPattern &pred, int id[4], const char *suffix, int flag)
{
  TString str;
  MakeFileName(str,id,suffix);

  ReadPatTXT(str.Data(),pred,flag);
  EdbID sid(id[0],id[1],id[2],id[3]);
  pred.SetSegmentsScanID(sid);

  LogPrint(id[0], 2, "ReadPatTXT","%s with %d predictions with flag: %d", 
	   str.Data(), pred.N(), flag);

  return pred.N();
}


//----------------------------------------------------------------
int EdbScanProc::ReadPatTXT(const char *file, EdbPattern &pred, int flag0)
{
  // read ascii predictions file as .../bXXXXXX/pYYY/a.a.a.a.suffix
  //        man      - for manual check by sysal
  //                 - overvise complete format
  // flag convention:
  //     -1 read all: basetracks, microtracks and holes (extrapolated predictions)
  //    >=0 read only segments with this flag
  //    -11 read basetracks or microtracks, skip holes (expected behaviour for display users)

  EdbSegP s;
  Int_t   ids = 0, flag = 0, ic = 0;
  Float_t x=0,y=0,tx=0,ty=0,sx=0,sy=0,stx=0,sty=0;

  char    buffer[256];

  FILE *f = fopen(file, "r");

  if (!f) {
    Log(1,"ReadPatTXT","ERROR! can not open file %s", file); 
    return 0; 
  }
  
  fgets (buffer, sizeof(buffer), f);
  int ncolumns = sscanf(buffer,"%d %f %f %f %f %f %f %f %f %d", 
			&ids,&x,&y,&tx,&ty,&sx,&sy,&stx,&sty,&flag);

  if (ncolumns!=6 && ncolumns!=10) {
    Log(1,"ReadPatTXT","ERROR! cannot recognize the format of file %s", file); 
    return 0; 
  }

  rewind(f);

  while (fgets (buffer, sizeof(buffer), f)) {
    if (sscanf(buffer,"%d %f %f %f %f %d", 
	       &ids,&x,&y,&tx,&ty,&flag) !=ncolumns )      break;
    if (flag0 > -1 && flag0 != flag)                       continue;
    if (flag0 == -11 && EdbRunTracking::GetMTHoles(flag)>0) continue;
    int side = flag%10;
    float ww=50; 
    s.SetDZ(-214);
    if(side!=0)       ww=25;
    if(side==1)  s.SetDZ(45);
    if(side==2)  s.SetDZ(-45);

    s.Set(ids,x,y,tx,ty,ww,flag);
    s.SetSide( side );
    if (ncolumns==10) s.SetErrors(sx,sy,0.,stx,sty);
    else s.SetErrors(50,50,0.,0.6,0.6);
    pred.AddSegment(s);
    ic++;
  }

  fclose(f);

  return ic;
}

//----------------------------------------------------------------
int EdbScanProc::WritePatRoot(EdbPattern &pred, int id[4], const char *suffix, int flag)
{
  // write root predictions file as .../bXXXXXX/pYYY/a.a.a.a.suffix
  int n = pred.N();
  TString str;
  MakeFileName(str,id,suffix);

  // checking for the existing directory

  FileStat_t buf;

  if (gSystem->GetPathInfo(gSystem->DirName(str), buf)) {
    if (gEDBDEBUGLEVEL > 0) {
      cout << "ERROR! Directory " << gSystem->DirName(str) 
	   << " does not exist.\n";
      cout << "Please run scanning procedure for the plate number " << id[1] 
	   << endl;
    }
    return 0;
  }

  TFile f(str.Data(), "RECREATE");

  if (flag < 0) pred.Write("pat");
  else {
    EdbPattern pat;
    pat.SetID(pred.ID());
    pat.SetPID(pred.PID());
    pat.SetX(pred.X());
    pat.SetY(pred.Y());
    pat.SetZ(pred.Z());
    for (Int_t i = 0; i < pred.N(); i++) 
      if (pred.GetSegment(i)->Flag()==flag) pat.AddSegment(*(pred.GetSegment(i)));
    n = pat.N();
    pat.Write("pat");
  }
  f.Close();
  LogPrint(id[0],2,"WritePatRoot","%s with %d predictions with flag: %d", 
	   str.Data(),n,flag);
  return n;
}

//----------------------------------------------------------------
int EdbScanProc::ReadPatRoot(EdbPattern &pred, int id[4], const char *suffix, int flag)
{
  // read root predictions file as .../bXXXXXX/pYYY/a.a.a.a.suffix
  // flag convention:
  //     -1 read all: basetracks, microtracks and holes (extrapolated predictions)
  //    >=0 read only segments with this flag
  //    -11 read basetracks or microtracks, skip holes (expected behaviour for display users)

  int n=0;
  TString str;
  MakeFileName(str,id,suffix);

  // checking for the existing directory

  FileStat_t buf;

  if (gSystem->GetPathInfo(gSystem->DirName(str), buf)) {
    if (gEDBDEBUGLEVEL > 0) {
      cout << "ERROR! Directory " << gSystem->DirName(str) 
	   << " does not exist.\n";
      cout << "Please run scanning procedure for the plate number " << id[1] 
	   << endl;
    }
    return 0;
  }

  TFile f(str.Data());
  EdbPattern *p=0;
  f.GetObject("pat",p);
  if (!p) {f.Close(); return 0;}
  for(int i=0; i<p->N(); i++) {
    EdbSegP *s = p->GetSegment(i);
    if (flag > -1  && flag != s->Flag()) continue;
    if (flag == -11 && EdbRunTracking::GetMTHoles(s->Flag())>0) continue;
    pred.AddSegment(*s)->SetSide( s->Flag()%10 );
    n++;
  }
  EdbID sid(id[0],id[1],id[2],id[3]);
  pred.SetSegmentsScanID(sid);
  f.Close();
  p->Delete();
  LogPrint(id[0],2,"ReadPatRoot","%s with %d predictions with flag: %d", str.Data(),n,flag);
  return n;
}

//----------------------------------------------------------------
bool EdbScanProc::CheckDirWritable(const char *dir)
{
  if( gSystem->AccessPathName(dir, kWritePermission) )   //can not access file!
    {
      Log(1,"CheckDirWritable","ERROR: can not open output directory: %s !!!",dir);
      return 0;
    }
  return 1;
}

//----------------------------------------------------------------
bool EdbScanProc::CheckDir(const char *dir, bool create)
{
  // check the existance of the directory dir
  // if do not exist and create==true (default) - create it 
  // return true if the directory exists or succesfully created

  void *dirp=0; // pointer to the directory
  dirp = gSystem->OpenDirectory(dir);
  if(!dirp) {
    if(create) {
      if(gSystem->MakeDirectory(dir)==0) Log(2,"EdbScanProc::CheckDir","create directory %s", dir);
      else                               Log(1,"EdbScanProc::CheckDir","ERROR! can not create directory %s", dir);
    }
    else return false;
  } else  { gSystem->FreeDirectory(dirp); dirp=0;}
  dirp = gSystem->OpenDirectory(dir);
  if(!dirp) {
     Log(1,"EdbScanProc::CheckDir","ERROR! directory %s is not created!", dir);
    return false;
  }
  gSystem->FreeDirectory(dirp);
  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::CheckAFFDir(int brick, bool create)
{
  //return true if dir ../bXXXXXX/AFF is exists, if create==true (default) create it if necessary 
  char str[256];
  sprintf(str,"%s/b%6.6d/AFF", eProcDirClient.Data(),brick);
  return CheckDir(str,create);
}

//----------------------------------------------------------------
bool EdbScanProc::CheckBrickDir( EdbID id, bool create )
{
  //return true if dir ../bXXXXXX exist, if create==true (default) create it if necessary 
  char str[256];
  sprintf(str,"%s/b%6.6d", eProcDirClient.Data(),id.eBrick);
  if(!CheckDir(str,create)) return false;
  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::CheckPlateDir( EdbID id, bool create )
{
  //return true if dir ../bXXXXXX/pXXX exist, if create==true (default) create it if necessary 
  char str[256];
  sprintf(str,"%s/b%6.6d/p%3.3d", eProcDirClient.Data(),id.eBrick,id.ePlate);
  if(!CheckDir(str,create)) return false;
  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::CheckProcDir( int id[4], bool create)
{
  //return true if dir ../bXXXXXX/pXXX exist, if create==true (default) create it if necessary 
  EdbID idd( id[0], id[1], id[2], id[3] );
  if(!CheckBrickDir(idd,create)) return false;
  if(!CheckPlateDir(idd,create)) return false;
  return true;
}

//----------------------------------------------------------------
char *EdbScanProc::BrickDir( int brick )
{
  char *str = new char[256];
  sprintf(str,"%s/b%6.6d",
          eProcDirClient.Data(),brick);
  return str;
}

//----------------------------------------------------------------
void EdbScanProc::MakeFileName(TString &s, int ID[4], const char *suffix, bool inplate)
{
  //make file pathname as .../bXXXXXX/pYYY/a.a.a.a.suffix if inplate==true
  //otherwise as .../bXXXXXX/a.a.a.a.suffix
  char str[256];
  if (inplate)
    sprintf(str,"%s/b%6.6d/p%3.3d/%d.%d.%d.%d.%s",
            eProcDirClient.Data(),ID[0], ID[1], ID[0], ID[1], ID[2], ID[3],suffix);
  else
    sprintf(str,"%s/b%6.6d/b%6.6d.%d.%d.%d.%s",
            eProcDirClient.Data(),ID[0], ID[0], ID[1], ID[2], ID[3],suffix);
    
  s=str;
}

//----------------------------------------------------------------
void EdbScanProc::MakeFileNameSrv(TString &s, int ID[4], const char *suffix, bool inplate)
{
  //make file pathname as .../bXXXXXX/pYYY/a.a.a.a.suffix if inplate==true
  //otherwise as .../bXXXXXX/a.a.a.a.suffix
  char str[256];
  if (inplate)
    sprintf(str,"%s/b%6.6d/p%3.3d/%d.%d.%d.%d.%s",
            eProcDirServer.Data(),ID[0], ID[1], ID[0], ID[1], ID[2], ID[3],suffix);
  else
    sprintf(str,"%s/b%6.6d/b%6.6d.%d.%d.%d.%s",
            eProcDirServer.Data(),ID[0], ID[0], ID[1], ID[2], ID[3],suffix);
    
  s=str;
}

//----------------------------------------------------------------
void EdbScanProc::MakeAffName(TString &s, int id1[4], int id2[4], const char *suffix)
{
  //make affine file name as ../bXXXXXX/AFF/a.a.a.a_b.b.b.b.aff.par
  char str[256];
  sprintf(str,"%s/b%6.6d/AFF/%d.%d.%d.%d_%d.%d.%d.%d.%s",
	  eProcDirClient.Data(), id1[0],
	  id1[0], id1[1], id1[2], id1[3],
	  id2[0], id2[1], id2[2], id2[3],
	  suffix);
  s=str;
}

//----------------------------------------------------------------
bool EdbScanProc::GetMap(int brick, TString &map)
{
  // get map string from the map file of this brick : .../bXXXXXX/bXXXXXX.map
  char str[256];
  sprintf(str,"%s/b%6.6d/b%6.6d.map", eProcDirClient.Data(),brick,brick);
  LogPrint(brick,1,"GetMap"," from file: %s\n",str);
  FILE *f = fopen(str,"r");
  if(!f)     {  LogPrint(brick,1,"GetMap","no map file: %s !!!\n",str); return false; }
  else if (fgets (str, 256, f) == NULL) 
    {  LogPrint(brick,1,"GetMap","error reading map file: %s !!!\n",str); return false; }
  LogPrint(brick,1,"GetMap","%s\n",str);
  map=str;
  fclose(f);
  return true;
}

//----------------------------------------------------------------
bool EdbScanProc::WaitFileReady(const char* fname_){//waits file copied/moved, in ready state
time_t t0 = time(0);
  const int dt=5; //wait max 5 sec till file ready
  const int dtupd=500; //wait 500 msec to check if file is updating
  const int dtslp=100; //wait 100 msec before recheck
  
  bool ready = false;
  LogPrint(0,2,"WaitFileReady","source file \"%s\"\n",fname_);
  while(t0+dt > time(0) && !ready){
    FILE *f=fopen(fname_, "r");
    if(f != NULL){
      fclose(f);
      //printf("file there\n");
      struct stat statbuf1, statbuf2;
      int v=stat(fname_, &statbuf1);
      if (v != -1) {
        if(statbuf1.st_size != 0){
          gSystem->Sleep(dtupd);
          if (stat(fname_, &statbuf2) != -1) {
            if(statbuf2.st_size == statbuf1.st_size){
              ready=true;
            }else{
              //printf("file updating\n");
            }
          }
        }
      }
      //printf("%d\n", v);
    }else{
      //printf("file not there\n");
      gSystem->Sleep(dtslp);

    }
  }
  LogPrint(0,2,"WaitFileReady","%s\n", ready ? "file is ready." : "file waiting timeout!");
  return ready;
}

EdbRun *EdbScanProc::InitRun(int id[4], char* runname_, char* runnamesrv_, bool createrun_)
{
  // create new run file as eProcDirClient/bXXXXXX/pYYY/x.y.s.p.raw.root
  if(!CheckProcDir(id)) return 0;
  TString str;
  TString strSrv;

  MakeFileName(str,id,"raw.root");   // the file will have the requested name 
  MakeFileNameSrv(strSrv,id,"raw.root");   // the file will have the requested name 

  if( !gSystem->AccessPathName(str.Data(), kFileExists) ) {   // if the file with the same name exists it will be saved as *root.xxx.save
    TString str2;
    for(int ic=0; ic<1000; ic++) {
      str2 = str; str2+="."; str2+=ic; str2+=".save";
      if( !gSystem->AccessPathName(str2.Data(), kFileExists) ) continue;
      else                                                     break;
    }
    //gSystem->CopyFile(str.Data(), str2.Data());
    char strbuf[1024];
#ifdef WIN32
    sprintf(strbuf,"move /F %s %s", str.Data(), str2.Data());
#else
    sprintf(strbuf,"mv -f %s %s", str.Data(), str2.Data());
#endif
    gSystem->Exec(strbuf);
    LogPrint(id[0], 3,"EdbScanProc::InitRun"," %s\n",str2.Data());
  }
  LogPrint(id[0], 3,"EdbScanProc::InitRun"," %s\n",str.Data());
  if(runname_ != NULL){
    strcpy(runname_, str.Data());
  }
  if(runnamesrv_ != NULL){
    strcpy(runnamesrv_, strSrv.Data());
  }
  EdbRun* run = NULL;
  if(createrun_)
    run = new EdbRun(str.Data(),"RECREATE");
  return run;
}
//-------------------------------------------------------------------
bool EdbScanProc::AddParLine(const char *file, const char *line, bool recreate)
{
  // add string to par file
  FILE *f=0;
  if(recreate) f = fopen(file,"w");
  else         f = fopen(file,"a+");
  if (!f) return false;
  fprintf(f,"%s\n",line);
  fclose(f);
  return true;
}

//-------------------------------------------------------------------
void EdbScanProc::MakeInParSet(EdbID id, const char *option)
{
  EdbScanSet *ss = ReadScanSet(id);
  int n = ss->eIDS.GetSize();
  for(int i=0; i<n; i++) {
    EdbID *id  = ss->GetID(i);
	if(id) MakeInPar(*id,option);
  }
}

//-------------------------------------------------------------------
bool EdbScanProc::MakeInPar(int id[4], const char *option)
{
  // prepare x.x.x.x.in.par file for the process defined in the option string
  TString name;
  MakeFileName(name,id,"in.par");
  FILE *f = fopen(name.Data(),"w");
  if(!f) {
    LogPrint(id[0],2,"MakeInPar","ERROR! can't open file: %s",name.Data() );
    return false;
  }
  if(eParDir.IsNull()) { 
    eParDir=eProcDirClient; eParDir+="/parset";  // for backword scripts compatibility
    Log(3,"EdbScanProc::MakeInPar","take default parameters from: %s",eParDir.Data());
  }
  fprintf(f,"INCLUDE %s/opera_emulsion.par\n",eParDir.Data());
  fprintf(f,"INCLUDE %s/%s.par\n",eParDir.Data(),option);
  TString nm;
  MakeFileName(nm,id,"par");
  fprintf(f,"INCLUDE %s\n",nm.Data());
  fclose(f);
  return true;
}

//-------------------------------------------------------------------
int EdbScanProc::CorrectAngles(int id[4])
{
  EdbDataPiece piece;
  InitPiece(piece,id);
  TString parfileOUT;
  MakeFileName(parfileOUT,id,"par");
  piece.eFileNamePar = parfileOUT;
  return piece.CorrectAngles();
}

//-------------------------------------------------------------------
int EdbScanProc::LinkRun(int id[4], int noUpdate)
{
  // the x.x.x.x.in.par file must be prepared before
  EdbDataPiece piece;
  InitPiece(piece,id);
  EdbDataProc  proc;
  proc.SetNoUpdate(noUpdate);
  TString parfileOUT;
  MakeFileName(parfileOUT,id,"par");
  piece.eFileNamePar = parfileOUT;
  return proc.Link(piece);
}

//----------------------------------------------------------------------------------------
int EdbScanProc::FindPredictionsRawSet(EdbID idp, EdbScanSet &ss, int last_plate)
{
  bool direction = idp.ePlate<last_plate ? true : false;
  int count=0;
  EdbID *idthis = ss.FindPlateID(idp.ePlate);
  if(!idthis) return 0;
  EdbID *idnext = 0;
  while(1) {
    FindPredictionsRaw(*idthis, *idthis);
    count++;
    idnext = ss.FindNextPlateID( idthis->ePlate, direction );  // find next plate
    if(!idnext) break;
    if(  direction && idnext->ePlate>last_plate ) break;
    if( !direction && idnext->ePlate<last_plate ) break;
    ProjectFound(*idthis,*idnext);
    idthis=idnext;
  }
  return count;
}

//----------------------------------------------------------------------------------------
int EdbScanProc::FindPredictionsRaw(EdbID idp, EdbID idr)
{
  // find raw microtracks for the predictions of idp in raw data of idr
  // Input:  idp.pred.root, idr.raw.root
  // Output: idp.found.root, idp.found.raw.txt

  EdbPattern pred;
  EdbPattern found;
  ReadPred(pred,idp);
  EdbRunAccess ra;
  InitRunAccess(ra,idr);

  EdbScanCond condBT;
  EdbScanCond condMT;
  SetDefaultCondBT(condBT);
  SetDefaultCondMT(condMT);
  float delta_theta = 0.1;
  float puls_min    = 7;
  float chi2max     = 1.6;
  int nfound = FindPredictionsRaw(pred,found,ra, condBT,condMT, delta_theta, puls_min, chi2max);
  WriteFound(found, idp);

  return nfound;
}

//----------------------------------------------------------------------------------------
int EdbScanProc::FindPredictionsRaw( EdbPattern &pred, EdbPattern &fnd, EdbRunAccess &ra, 
				     EdbScanCond &condBT, EdbScanCond &condMT, 
				     float delta_theta, float puls_min, float puls_mt, float chi2max, FILE *out )
{
  // find raw microtracks for the predictions "pred" in run "ra"
  // Input:  pred,ra
  // Output: fnd - basetracks with position taken from the best microtrack (if any) and the angle of the predicted basetrack

  Log(2,"FindPredictionsRaw"," search for %d predictions \n", pred.N());

  TFile ftree("micro.root","RECREATE");
  EdbSegP      *s_b    = new EdbSegP();
  EdbSegP      *sf_b  = new EdbSegP();
  TClonesArray *pat1_b = new TClonesArray("EdbSegP");
  TClonesArray *pat2_b = new TClonesArray("EdbSegP");
  TTree micro("micro","micro");
  micro.Branch("s.","EdbSegP",&s_b,32000,99);
  micro.Branch("sf.","EdbSegP",&sf_b,32000,99);
  micro.Branch("s1.",&pat1_b,32000,99);
  micro.Branch("s2.",&pat2_b,32000,99);

  if(!out) out = fopen("micro.txt","w");
  if(!out) return 0;
  if(gEDBDEBUGLEVEL>2) ra.Print();
  
  for(int i=0; i<pred.N(); i++) {
    ra.ClearCuts();
    EdbSegP s;
    s.Copy( *(pred.GetSegment(i)) );
    s.SetZ(107.);                                      // TODO!
    s.SetErrors();
    float xmin[5]={-500, -500, s.TX()-delta_theta, s.TY()-delta_theta,  puls_min };         //TODO!!
    float xmax[5]={ 500,  500, s.TX()+delta_theta, s.TY()+delta_theta,  50       };
    condBT.FillErrorsCov( s.TX(), s.TY(), s.COV() );
    fprintf(out,"\n%8.8d   %11.2f %11.2f %7.4f %7.4f\n",s.ID(),s.X(),s.Y(),s.TX(),s.TY());

    EdbPVRec aview; //with 2 patterns of preselected microtracks
    aview.AddPattern( new EdbPattern(0,0,214));  // TODO! sequence??
    aview.AddPattern( new EdbPattern(0,0,0)  );

    for(int side=1; side<=2; side++) {
      Log(3,"FindPredictionsRaw","side = %d\n",side);
      EdbPattern pat;
      ra.AddSegmentCut(side,1,xmin,xmax);
      ra.GetPatternXY( s, side,  pat );

      for(int i=0; i<pat.N(); i++) {
	EdbSegP *ss = pat.GetSegment(i);
	ss->SetErrors();
	condMT.FillErrorsCov( s.TX(), s.TY(), ss->COV() );
      }
      pat.FillCell(10,10,0.01,0.01);  //divide view on this cells
      
      TArrayF   chi2arr(1000);  //TODO!
      TObjArray found;
      int nf= FindCompliments(s,pat,found, chi2max, chi2arr);
      for(int j=0; j<found.GetEntriesFast(); j++) {
	EdbSegP *s2 = (EdbSegP *)(found.At(j));
	s2->SetChi2(chi2arr[j]);                        // TODO -???
	aview.GetPattern(side-1)->AddSegment(*s2);
	float offx=s2->X()-(s.X()+s.TX()*(s2->Z()-s.Z()));
	float offy=s2->Y()-(s.Y()+s.TY()*(s2->Z()-s.Z()));
	fprintf(out,"s%1d(%2d)%4d %11.2f %11.2f %7.4f %7.4f %6.1f %7.2f %7.2f %7.4f %7.4f %6.3f %3.0f\n",
		side,nf,s2->ID(),s2->X(),s2->Y(),s2->TX(),s2->TY(),s2->Z(),
		offx,offy,
		s2->TX()-s.TX(),s2->TY()-s.TY(),chi2arr[j],s2->W());
      }
    }

    EdbSegP sfmt;      // container for the found best microtrack passed the cuts
    sfmt.SetChi2(10.*chi2max);

    float rlim   = 20;        // TODO
    float chi, chimin = 1.5;  // TODO  for bt selection
    EdbSegP *s1b=0, *s2b=0;   // the best bt
    EdbSegP s3;

    for(int is1=0; is1<aview.GetPattern(0)->N(); is1++) {
      EdbSegP *s1 = aview.GetPattern(0)->GetSegment(is1);
      if( sfmt.Chi2() > s1->Chi2() && s1->W() >= puls_mt ) {sfmt.Copy(*s1); sfmt.SetFlag(1);}
    }
    for(int is2=0; is2<aview.GetPattern(1)->N(); is2++) {
      EdbSegP *s2 = aview.GetPattern(1)->GetSegment(is2);
      if( sfmt.Chi2() > s2->Chi2() && s2->W() >= puls_mt ) {sfmt.Copy(*s2); sfmt.SetFlag(2);}
    }
    
    for(int is1=0; is1<aview.GetPattern(0)->N(); is1++) {
      for(int is2=0; is2<aview.GetPattern(1)->N(); is2++) {
	EdbSegP *s1 = aview.GetPattern(0)->GetSegment(is1);
	EdbSegP *s2 = aview.GetPattern(1)->GetSegment(is2);

	float dx1=s1->X()-(s.X()+s.TX()*(s1->Z()-s.Z()));
	float dy1=s1->Y()-(s.Y()+s.TY()*(s1->Z()-s.Z()));
	float dx2=s2->X()-(s.X()+s.TX()*(s2->Z()-s.Z()));
	float dy2=s2->Y()-(s.Y()+s.TY()*(s2->Z()-s.Z()));
	float r = Sqrt( (dx1-dx2)*(dx1-dx2) + (dy1-dy2)*(dy1-dy2) ); 
	fprintf(out,"r = %7.2f  ",r);
	if(r<rlim) {  // has good BT
	  s3.Copy(s);
	  s3.SetX( 0.5*(s1->X() + s2->X()) );
	  s3.SetY( 0.5*(s1->Y() + s2->Y()) );
	  s3.SetZ( 0.5*(s1->Z() + s2->Z()) );
	  s3.SetTX( (s2->X() - s1->X()) / (s2->Z() - s1->Z()) );
	  s3.SetTY( (s2->Y() - s1->Y()) / (s2->Z() - s1->Z()) );
	  s3.SetFlag(0);

	  //s3.Print();
	  //s.Print();
	  chi = EdbTrackFitter::Chi2Seg(&s3, &s);
	  //printf("chi = %7.4f\n",chi);
	  fprintf(out,"chi = %7.4f\n",chi);
	  if(chi<chimin) {                           //select the best basetrack
	    chimin = chi;
	    s1b = s1;
	    s2b = s2;
	  }
	}
      }
    }

    EdbSegP sf;                // container for the found track

    int bth, mth, tb; // basetrack holes, mt-holes, top/bottom;

    sf.Copy(s);
    if(s1b&&s2b) {
      sf.SetX( 0.5*(s1b->X() + s2b->X()) );
      sf.SetY( 0.5*(s1b->Y() + s2b->Y()) );
      sf.SetZ( 0.5*(s1b->Z() + s2b->Z()) );
      sf.SetTX( (s2b->X() - s1b->X()) / (s2b->Z() - s1b->Z()) );
      sf.SetTY( (s2b->Y() - s1b->Y()) / (s2b->Z() - s1b->Z()) );
      sf.SetFlag(0);                                              // if bt found : bth=0, mth=0, tb=0
      sf.SetChi2(chimin);
    } else if(sfmt.Chi2()<chi2max) {  // found good microtrack
      //float zmean = ra.GetLayer(0)->Z();
      float zmean = s.Z();
      //printf("Zmt = %f  zmean = %f   Zs = %f \n",sfmt.Z(),zmean,s.Z());
      sf.SetX( sfmt.X() + s.TX()*(zmean-sfmt.Z()) );
      sf.SetY( sfmt.Y() + s.TY()*(zmean-sfmt.Z()) );
      sf.SetZ(zmean);
      bth = s.Flag()/10000;  bth++;
      mth = 0;
      tb = sfmt.Flag();
      sf.SetFlag(bth*10000+mth*100+tb);                   // if mt found : bth++, mth=0, tb=1/2
    } else {
      bth = s.Flag()/10000;      bth++;
      mth = (s.Flag()/100)%100;  mth++;
      tb  =  s.Flag()%10;
      sf.SetFlag(bth*10000+mth*100+tb);         // hole: if not found: bth++, mth++, tb= keep last value
    }

    s_b->Copy(s);
    sf_b->Copy(sf);
    pat1_b = aview.GetPattern(0)->GetSegments();
    pat2_b = aview.GetPattern(1)->GetSegments();
    micro.SetBranchAddress("s1."  , &pat1_b );
    micro.SetBranchAddress("s2."  , &pat2_b );
    micro.Fill();

    fnd.AddSegment(sf);
  }
  
  //ftree.cd();
  micro.Write();
  ftree.Close();
  fclose(out);
  return 1; //TODO!
}

//----------------------------------------------------------------------------------------
int EdbScanProc::FindCompliments( EdbSegP &s, EdbPattern &pat, TObjArray &found, float chi2max, TArrayF &chiarr )
{
  // return found sorted by increasing chi2

  int nfound=0;
  int maxcand=chiarr.GetSize();
  TArrayF   chi2arr(maxcand);
  TObjArray arr(maxcand);
  TArrayI   ind(maxcand);
  
  int nseg = pat.FindCompliments(s,arr,30,200);  // acceptance (prelim): s.SX()*30; s.STX*200
  //printf("\nnseg = %d\n",nseg);
  if(nseg>maxcand)  {
    Log(1,"FindCompliments","Warning!: Too many segments %d, accept only the first %d", nseg, maxcand);
    nseg = maxcand;
  }
  if(nseg<=0) return 0;

  EdbSegP *s2=0;
  for(int j=0; j<nseg; j++) {
    s2 = (EdbSegP *)arr.At(j);
    EdbSegP s3;
    s3.Copy(s);
    chi2arr[j] = EdbTrackFitter::Chi2Seg(&s3, s2);
  }
  TMath::Sort(nseg,chi2arr.GetArray(),ind.GetArray(),0);
  for(int j=0; j<nseg; j++) {
    s2 = (EdbSegP *)arr.At(ind[j]);
    if(chi2arr[ind[j]] > chi2max ) break;
    chiarr[j] = chi2arr[ind[j]];
    s2->SetMC(s.MCEvt(),s.MCTrack());
    found.Add(s2);
    nfound++;
  }

  //printf("nfound = %d\n",nfound);
  return nfound;
}

//----------------------------------------------------------------------------------------
void EdbScanProc::SetDefaultCondBT(EdbScanCond &cond)
{
  cond.SetSigma0( 10., 10., 0.007, 0.007 );   // sigma0 "x, y, tx, ty" at zero angle
  cond.SetDegrad( 5. );                       // sigma(tx) = sigma0*(1+degrad*tx)
  cond.SetBins(0, 0, 0, 0); //???                  // bins in [sigma] for checks
  cond.SetPulsRamp0(  5., 5. );               // in range (Pmin:Pmax) Signal/All is nearly linear
  cond.SetPulsRamp04( 5., 5. );
  cond.SetChi2Max( 6.5 );
  cond.SetChi2PMax( 6.5 );
  cond.SetRadX0( 5810. );
  cond.SetName("OPERA_basetrack");
}

//----------------------------------------------------------------------------------------
void EdbScanProc::SetDefaultCondMT(EdbScanCond &cond)
{
  cond.SetSigma0( 1., 1., 0.025, 0.025 );   // sigma0 "x, y, tx, ty" at zero angle
  cond.SetDegrad( 5. );                       // sigma(tx) = sigma0*(1+degrad*tx)
  cond.SetBins(0, 0, 0, 0);  //???                 // bins in [sigma] for checks
  cond.SetPulsRamp0(  5., 5. );               // in range (Pmin:Pmax) Signal/All is nearly linear
  cond.SetPulsRamp04( 5., 5. );
  cond.SetChi2Max( 6.5 );
  cond.SetChi2PMax( 6.5 );
  cond.SetRadX0( 5810. );
  cond.SetName("OPERA_microtrack");
}

//-------------------------------------------------------------------
bool EdbScanProc::InitRunAccess(EdbRunAccess &ra, int id[4], bool do_update)
{
  // initialize the EdbRunAccess object useful for the raw data handling
  EdbDataPiece p;
  if(!InitPiece(p,id)) return false;
  if(gEDBDEBUGLEVEL>2) p.Print();
  ra.eAFID = p.eAFID;
  if( !ra.InitRun(p.GetRunFile(0), do_update) ) {
    LogPrint(id[0],1,"InitRunAccess","ERROR open file %s !!!",p.GetRunFile(0));
    return false;
  } else
    LogPrint(id[0],2,"InitRunAccess"," %s with %d views",p.GetRunFile(0), ra.GetRun()->GetEntries() );

  ra.GetLayer(1)->SetZlayer( p.GetLayer(1)->Z(),p.GetLayer(1)->Zmin(),p.GetLayer(1)->Zmax());
  ra.GetLayer(2)->SetZlayer( p.GetLayer(2)->Z(),p.GetLayer(2)->Zmin(),p.GetLayer(2)->Zmax());
  ra.GetLayer(1)->SetShrinkage( p.GetLayer(1)->Shr());
  ra.GetLayer(2)->SetShrinkage( p.GetLayer(2)->Shr());

  EdbAffine2D *a1 =  p.GetLayer(1)->GetAffineTXTY();
  EdbAffine2D *a2 =  p.GetLayer(2)->GetAffineTXTY();
  ra.GetLayer(1)->SetAffTXTY( a1->A11(),a1->A12(),a1->A21(),a1->A22(), a1->B1(), a1->B2() );
  ra.GetLayer(2)->SetAffTXTY( a2->A11(),a2->A12(),a2->A21(),a2->A22(), a2->B1(), a2->B2() );

  return true;
}

//-------------------------------------------------------------------
bool EdbScanProc::InitPiece(EdbDataPiece &piece, int id[4])
{
  // set raw, cp and par for the piece according to id
  TString runfile, cpfile, parfile;
  MakeFileName(runfile,id,"raw.root");
  MakeFileName(cpfile,id,"cp.root");
  MakeFileName(parfile,id,"in.par");
  piece.AddRunFile(runfile);
  piece.eFileNameCP  = cpfile;
  piece.eFileNamePar = parfile;
  if(piece.TakePiecePar()<0) Log(1,"InitPiece","Warning: file %s does not exist!",parfile.Data());
  return true;
}

//-------------------------------------------------------------------
int EdbScanProc::ReadPatCP(EdbPattern &pat, int id[4], TCut cut)
{
  // read CP file ("base" segments) applying all cuts and transformations from x.x.x.x.in.par
  // the Z of the pat and all segments will be z of layer 0 defined in the par file(s)
  EdbDataPiece piece;
  InitPiece(piece, id);
  piece.GetLayer(0)->SetZlayer(piece.GetLayer(0)->Z(), 0,0);
  piece.AddRCut(0,cut);
  int n = ReadPiece(piece, pat);
  pat.SetZ(piece.GetLayer(0)->Z());
  pat.SetSegmentsZ();

  EdbID sid(id[0],id[1],id[2],id[3]);
  pat.SetSegmentsScanID(sid);  

  return n;
}

//-------------------------------------------------------------------
int EdbScanProc::ReadPiece(EdbDataPiece &piece, EdbPattern &pat)
{
  //assuming that piece was initialised correctly
  if(!piece.InitCouplesTree("READ")) return 0;
  piece.GetCPData_new( &pat,0,0,0 );
  pat.SetSegmentsZ();
  pat.Transform(    piece.GetLayer(0)->GetAffineXY()   );
  pat.TransformA(   piece.GetLayer(0)->GetAffineTXTY() );
  pat.TransformShr( piece.GetLayer(0)->Shr()  );
  return 1;
}

//-------------------------------------------------------------------
int EdbScanProc::FindPredictions(EdbPattern &pred, int id[4], EdbPattern &found, int maxholes)
{
  // find predictions pred in couples tree of id and prepare pattern "found"
  // assumed that pred are transformed and projected into the coord system of id
  // Input:   pred - pattern with predictions
  //            id - the data piece to be processed
  //      maxholes - the maximum number of holes (missed segments) for doing extrapolation
  // Output: found - pattern with found tracks
  //         x.x.x.x.found.txt summary file with all candidats
  //

  //   Probably obsolete function - to investigate if it in use now - in most 
  //   of cases it can be substituted by EdbRunTracking::FindPredictions            (VT, AC)
  //

  EdbPVRec ali;
  EdbPattern *pat=0;
  EdbDataPiece piece;

  // predicted:
  pat = new EdbPattern(pred.X(),pred.Y(),0,100);
  for(int i=0; i<pred.N(); i++) pat->AddSegment(*(pred.GetSegment(i)));
  pat->SetPID(0);
  pat->SetSegmentsZ();      // z=0  (the same)
  ali.AddPattern(pat);

  // scanned:
  InitPiece(piece, id);
  EdbPattern *patbt = new EdbPattern(0.,0., 0,100 );
  EdbPattern *pat1  = new EdbPattern(0.,0., 0,100 );
  EdbPattern *pat2  = new EdbPattern(0.,0., 0,100 );

  if(!piece.InitCouplesTree("READ")) return 0;
  piece.GetCPData_new( patbt,pat1,pat2,0 );
  patbt->SetSegmentsZ();
  patbt->Transform(    piece.GetLayer(0)->GetAffineXY()   );
  patbt->TransformA(   piece.GetLayer(0)->GetAffineTXTY() );
  patbt->TransformShr( piece.GetLayer(0)->Shr()  );
  pat1->SetSegmentsZ();
  pat1->Transform(    piece.GetLayer(0)->GetAffineXY()   );
  pat1->TransformA(   piece.GetLayer(0)->GetAffineTXTY() );
  pat2->SetSegmentsZ();
  pat2->Transform(    piece.GetLayer(0)->GetAffineXY()   );
  pat2->TransformA(   piece.GetLayer(0)->GetAffineTXTY() );

  //ReadPiece(piece, *pat);
  patbt->SetSegmentsZ();      // z=0 (the same)
  patbt->SetPID(1);
  ali.AddPattern(patbt);

  EdbScanCond *cond = piece.GetCond(0);
  cond->SetChi2Mode(3);
  ali.SetScanCond( cond );
  ali.SetPatternsID();
  //ali.SetSegmentsErrors();
  ali.SetCouplesAll();
  ali.SetChi2Max(cond->Chi2PMax());

  EdbSegP *s=0;
  for(int ip=0; ip<ali.Npatterns(); ip++)
    for(int i=0; i<ali.GetPattern(ip)->N(); i++) {
      s = ali.GetPattern(ip)->GetSegment(i);
      s->SetErrors();
      cond->FillErrorsCov( s->TX(), s->TY(), s->COV() );
    }

  TString str;
  MakeFileName(str,id,"found.txt");
  FILE *f = fopen(str.Data(),"w");

  TString strmt;
  MakeFileName(strmt,id,"found.mt.txt");
  FILE *fmt = fopen(strmt.Data(),"w");

  int maxcand=100;
  TArrayF chiarr(maxcand);
  TArrayI ind(maxcand);
  TArrayI count(maxcand);
  TArrayI cnsel(maxcand);

  pat = ali.GetPattern(1);
  pat->FillCell(20,20,0.01,0.01);
  int nseg=0;
  TObjArray arr;
  EdbSegP *s2=0;
  EdbSegP s3;
  for(int i=0; i<pred.N(); i++) {
    s = ali.GetPattern(0)->GetSegment(i);
    arr.Clear();
    nseg = pat->FindCompliments(*s,arr,cond->BinX(),cond->BinTX());
    if(nseg>maxcand)       continue;
    count[nseg]++;
    int nsel=0;
    if(nseg>=0) {
      for(int j=0; j<nseg; j++) {
	s2 = (EdbSegP *)arr.At(j);
	s3.Copy(*s2);
	chiarr[j] = EdbTrackFitter::Chi2Seg(&s3, s);
      }
      TMath::Sort(nseg,chiarr.GetArray(),ind.GetArray(),0);
      for(int j=0; j<nseg; j++) {
	s2 = (EdbSegP *)arr.At(ind[j]);
	if(chiarr[ind[j]] > cond->Chi2PMax() ) break;
	nsel=j+1;
	s2->SetMC(s->MCEvt(),s->MCTrack());
      }

      fprintf(f,"\n%8.8d %11.2f %11.2f %8.4f %8.4f %d\n",
	      s->ID(),s->X(),s->Y(),s->TX(),s->TY(), nsel);
      fprintf(fmt,"\n%8.8d %11.2f %11.2f %8.4f %8.4f %d\n",
	      s->ID(),s->X(),s->Y(),s->TX(),s->TY(), nsel);
      for(int j=0; j<nsel; j++) {
	s2 = (EdbSegP *)arr.At(ind[j]);
	fprintf(f,"%8d %11.2f %11.2f %8.4f %8.4f %6.2f %3.0f\n",
		j+1,s2->X(),s2->Y(),s2->TX(),s2->TY(),chiarr[ind[j]],s2->W());
	fprintf(fmt,"%8d %11.2f %11.2f %8.4f %8.4f %6.2f %3.0f\n",
		j+1,s2->X(),s2->Y(),s2->TX(),s2->TY(),chiarr[ind[j]],s2->W());
	int imt=-1;
	if(patbt->GetSegments()->FindObject(s2)) imt = patbt->GetSegments()->IndexOf(s2);
	if(imt>0) {
	  EdbSegP *smt = pat1->GetSegment(imt);
	  fprintf(fmt,"s1:%5d %11.2f %11.2f %8.4f %8.4f %3.0f\n",
		  j+1,smt->X(),smt->Y(),smt->TX(),smt->TY(),smt->W());
	  smt = pat2->GetSegment(imt);
	  fprintf(fmt,"s2:%5d %11.2f %11.2f %8.4f %8.4f %3.0f\n",
		  j+1,smt->X(),smt->Y(),smt->TX(),smt->TY(),smt->W());
	}
      }

    }
    cnsel[nsel]++;
    if(nsel>0) {
      found.AddSegment(*((EdbSegP *)arr.At(ind[0])));   // add the best segment
      found.GetSegmentLast()->SetFlag(0);               // reset flag if found good candidate
      found.GetSegmentLast()->SetID(s->ID());
    }
    else if(s->Flag()<maxholes) {
      found.AddSegment(*(s));                       // add itself in case of hole
      found.GetSegmentLast()->SetFlag(s->Flag()+1); // flag is the number of missed plates // OLD FLAG DEFINITION
      found.GetSegmentLast()->SetID(s->ID());
    }
  }
  fclose(f);
  fclose(fmt);
  delete pat1; pat1=0;
  delete pat2; pat2=0;

  printf("Total: %d predictions, %d basetracks in scanned pattern\n",pred.N(), pat->N() );
  int sum=0;
  printf("Before chi2 cut: \n" );
  for(int i=0; i<maxcand; i++) 
    if(count[i]>0) {
      printf("count(%5d)= %5d\n",i, count[i] );
      sum+=count[i];
    }
  printf("sum = %d\n",sum );
  sum=0;
  printf("After chi2 cut: \n" );
  for(int i=0; i<maxcand; i++) 
    if(cnsel[i]>0) {
      printf("cnsel(%5d)= %5d\n",i, cnsel[i] );
      sum+=cnsel[i];
    }
  printf("sum = %d\n",sum );

  LogPrint(id[0],1,"FindPredictions","%d.%d.%d.%d:  %d out of %d predictions are found (%d-zero, %d-single, %d-multy),  maxholes=%d", 
	   id[0],id[1],id[2],id[3],sum-cnsel[0],pred.N(),cnsel[0],cnsel[1], sum-cnsel[0]-cnsel[1], maxholes);

  return sum-cnsel[0];
}

//-------------------------------------------------------------------
bool  EdbScanProc::GetAffZ(EdbAffine2D &aff, float &dz, int id1[4],int id2[4])
{
  // read affine transformations and deltaZ from x.x.x.x_y.y.y.y.aff.par
  TString parfile;
  MakeAffName(parfile,id1,id2);
  EdbDataPiece piece;
  piece.eFileNamePar = parfile;
  if (piece.TakePiecePar() < 0) return false;
  EdbAffine2D *a = piece.GetLayer(0)->GetAffineXY();
  if(!a) return false;
  aff.Set( a->A11(),a->A12(),a->A21(),a->A22(),a->B1(),a->B2() );
  dz = piece.GetLayer(0)->Z();
  return true;
}

//-------------------------------------------------------------------
bool  EdbScanProc::ApplyAffZ(EdbPattern &pat,int id1[4],int id2[4])
{
  // read affine transformations and deltaZ from x.x.x.x_y.y.y.y.aff.par and apply it to pat
  EdbAffine2D aff;
  float       dz;
  if( !GetAffZ( aff, dz, id1, id2) ) return false;
  if(gEDBDEBUGLEVEL>2) {
    printf("ApplyAffZ: dz = %f\n",dz);
    aff.Print();
  }
  pat.ProjectTo(dz);
  pat.Transform(&aff);
  return true;
}

//----------------------------------------------------------------
int EdbScanProc::AlignAll(int id1[4], int id2[4], int npre, int nfull, const char *opt)
{
  int nal=0;
  if (npre > 0) {
    MakeInPar(id1,"prealignment");
    MakeInPar(id2,"prealignment");
    for (Int_t i = 0; i < npre; i++) {
      // find affine transformation from id1 to id2 and update par of id1
      nal = Align(id1, id2, "");	
      if (nal == -1) return -1;
    }
  }
  if (nfull > 0) {
    MakeInPar(id1,"fullalignment");
    MakeInPar(id2,"fullalignment");
    for (Int_t i = 0; i < nfull; i++) {
      // find affine transformation from id1 to id2 and update par of id1
      nal = Align(id1, id2, opt);	
      if (nal == -1) return -1;
    }
  }
  LogPrint(id1[0],1,"AlignAll","%d.%d.%d.%d to %d.%d.%d.%d with %d pre and %d full.  The final pattern is: %d", 
	   id1[0],id1[1],id1[2],id1[3],id2[0],id2[1],id2[2],id2[3],npre,nfull,nal);
  return nal;
}

//-------------------------------------------------------------------
int EdbScanProc::Align(int id1[4], int id2[4], const char *option)
{
  // Align 2 patterns, assumed that already exists file x.x.x.x_y.y.y.y.aff.par with deltaZ inside.
  // Convension about Z(setted while process): the z of id2 is 0, the z of id1 is (-deltaZ) where
  // deltaZ readed from aff.par file in a way that pattern of id1 projected 
  // to deltaZ correspond to pattern of id2

  int npat=0;
  TString name;
  MakeFileName(name,id1,"in.par");
  TString parfileOUT;
  MakeAffName(parfileOUT,id1,id2);
  parfileOUT.Prepend("INCLUDE ");
  AddParLine(name.Data(), parfileOUT.Data());
  
  EdbPVRec ali;
  EdbPattern *pat=0;
  EdbDataPiece piece1,piece2;
  
  InitPiece(piece1, id1);
  piece1.GetLayer(0)->SetZlayer(-1*piece1.GetLayer(0)->Z(), 0,0);
  pat = new EdbPattern(0.,0., piece1.GetLayer(0)->Z(),100 );
  if (!ReadPiece(piece1, *pat)) {delete pat; return -1;}
  pat->SetPID(0);
  pat->SetSegmentsZ();
  ali.AddPattern(pat);
  
  InitPiece(piece2, id2);
  piece2.GetLayer(0)->SetZlayer(0, 0,0);
  pat = new EdbPattern(0.,0., piece2.GetLayer(0)->Z(),100 );
  if (!ReadPiece(piece2, *pat)) {delete pat; return -1;}
  pat->SetPID(1);
  pat->SetSegmentsZ();
  ali.AddPattern(pat);

  Log(2,"EdbScanProc::Align","Z1 = %f z2 = %f with option: %s",
      piece1.GetLayer(0)->Z(),piece2.GetLayer(0)->Z(), option );
  
  EdbScanCond *cond = piece1.GetCond(0);
  cond->SetChi2Mode(3);
  ali.SetScanCond( cond );
  ali.SetPatternsID();
  ali.SetSegmentsErrors();
  ali.SetCouplesAll();
  ali.SetChi2Max(cond->Chi2PMax());
  ali.SetOffsetsMax(cond->OffX(),cond->OffY());

  if( strstr( option,"-a") && strstr( option,"-a2")==0 )    ali.Align(0);
  else                                                      ali.Align(2); 

  ali.PrintAff();
  npat = ali.GetCouple(0)->Ncouples();

  EdbAffine2D  aff;
  MakeAffName(parfileOUT,id1,id2);
  piece1.eFileNamePar = parfileOUT;
  ali.GetPattern(0)->GetKeep(aff);
  piece1.UpdateAffPar(0,aff);
  
  TString  cpfile;
  MakeFileName(cpfile,id1,"al.cp.root");
  EdbDataProc proc;
  proc.LinkTracksWithFlag( &ali, 10., 0.05, 2, 3, 0 );

  TTree *cptree=EdbDataPiece::InitCouplesTree(cpfile.Data(),"RECREATE");

  proc.FillCouplesTree(cptree, &ali,0);
  proc.CloseCouplesTree(cptree);
  
  //ali.FitTracks( 10, 0.139 );         // is important to call it before MakeTracksTree!
  //MakeFileName(cpfile,id1,"al.tr.root");
  //proc.MakeTracksTree(&ali,cpfile.Data());
  if( strstr(option,"-z") ) {
    ali.FineCorrZnew();
    piece1.UpdateZPar(0,-ali.GetPattern(0)->Z());
  }
  
  return npat;
 }

//______________________________________________________________________________
void EdbScanProc::UpdateSetWithAff(EdbID idset, EdbID idset1, EdbID idset2)
{
  // in this function the geometry of the brick of idset is updated with the affine transformations
  // found for the couple idset1 idset2 (same plate rescanning transformations)
  EdbScanSet *ss  = ReadScanSet(idset);    if(!ss)  return;
  EdbScanSet *ss1 = ReadScanSet(idset1);   if(!ss1) return;
  EdbScanSet *ss2 = ReadScanSet(idset2);   if(!ss2) return;
  int n = ss->eIDS.GetSize();
  for(int i=0; i<n; i++) {
    EdbID *id   = ss->GetID(i);
    EdbID *id1  = ss1->FindPlateID(id->ePlate);    if(!id1) continue;
    EdbID *id2  = ss2->FindPlateID(id->ePlate);    if(!id2) continue;

    EdbLayer la;    ReadAffToLayer( la, *id1, *id2);
    EdbPlateP *p = ss->GetPlate(id->ePlate);
    if(p) {
      p->SetZcorr(la.Zcorr());
      p->GetAffineXY()->Transform(la.GetAffineXY());
      p->GetAffineTXTY()->Transform(la.GetAffineTXTY());
      p->SetShrinkage(la.Shr());
    }
  }
  WriteScanSet( idset ,*ss );
  SafeDelete(ss);   SafeDelete(ss1);   SafeDelete(ss2);
}

//______________________________________________________________________________
void EdbScanProc::UpdateSetWithAff(EdbID idset, EdbAffine2D aff)
{
  // in this function the geometry of the brick of idset is updated with the given affine transformations
  EdbScanSet *ss  = ReadScanSet(idset);    if(!ss)  return;
  int n = ss->eIDS.GetSize();
  for(int i=0; i<n; i++) {
    EdbID *id   = ss->GetID(i);
    EdbPlateP *p = ss->GetPlate(id->ePlate);
    if(p) {
      p->GetAffineXY()->Transform(&aff);
    }
  }
  WriteScanSet( idset ,*ss );
  SafeDelete(ss);
}

//______________________________________________________________________________
void EdbScanProc::UpdateSetWithAff(EdbID idset, EdbID idsetu )
{
  // in this function the geometry of the brick of idset is updated with the affine transformations
  // found for plate-to-plate alignment of idset
  EdbScanSet  *ss  = ReadScanSet(idset);     if(!ss)  return;
  EdbScanSet *ssu  = ReadScanSet(idsetu);    if(!ssu)  return;
  int n = ssu->eIDS.GetSize();            if(n<2)  return;
  for(int i=0; i<n-1; i++) {
    EdbID *id1   = ssu->GetID(i);
    EdbID *id2   = ssu->GetID(i+1);
    EdbLayer la;    if(!ReadAffToLayer( la, *id1, *id2)) continue;
    
    TString mapfilename; MakeAffName(mapfilename,*id1,*id2,"map.root");
    EdbCorrectionMap *map; 
    TFile *mapfile = TFile::Open(mapfilename.Data());
    if(mapfile){ 
     map = (EdbCorrectionMap*) mapfile->Get("corrmap");
     std::cout<<"Adding correction map"<<std::endl;
     if(map) la.ApplyCorrectionsLocal(*map);
    }
    
    ss->UpdateBrickWithP2P(la,id1->ePlate,id2->ePlate);
  }
  WriteScanSet( idset ,*ss );
  SafeDelete(ss);  SafeDelete(ssu);
}

//----------------------------------------------------------------
void EdbScanProc::UpdateSetWithPlatePar(EdbID idset)
{
 EdbScanSet  *ss  = ReadScanSet(idset);     if(!ss)  return;
 UpdateSetWithPlatePar(*ss);
 WriteScanSet( idset ,*ss );
 SafeDelete(ss);
}

//----------------------------------------------------------------
void EdbScanProc::UpdateSetWithPlatePar(EdbScanSet &ss)
{
  int n = ss.eIDS.GetSize();
  Log(3,"EdbScanProc::UpdateSetWithPlatePar","set with %d plates",n);
  for(int i=0; i<n; i++) {
    EdbID *id  = ss.GetID(i);
    EdbPlateP  *plate = ss.GetPlate(id->ePlate);
    if(!plate)   Log(1,"EdbScanProc::UpdateSetWithPlatePar","ERROR! plate %d do not found!",id->ePlate);
    //printf("plate %d before\n",i); plate->Print();
    if(plate) ReadPiecePar( *id, *plate);
    //printf("after\n"); plate->Print();
  }
}

//______________________________________________________________________________
void EdbScanProc::MakeAlignSetSummary(EdbID idset1, EdbID idset2, const char *file, const char *mode )
{
  // assuming that exist the scan sets for idset1 and idset2 
  // read id1_id2.aff.par and make a summary tree

  Log(2,"MakeAlignSetSummary","open file %s for %s",file,mode);
  TFile *f = file?  new TFile(file,mode): 0;
  if(!f) return;
  TTree *tree = (TTree*)f->Get("alsum");
  if(!tree) {
    tree  = new TTree("alsum","Alignment Summary");
    Int_t   peak=0;
    //Float_t x0,y0,dx,dy,dz1,dz2;
    EdbID *id1=0, *id2=0;
    EdbLayer *corr = 0;
    EdbPeak2 *peak2c = 0;
    Float_t xoffset=0,yoffset=0;
    tree->Branch("id1", "EdbID", &id1);
    tree->Branch("id2", "EdbID", &id2);
    tree->Branch("corr", "EdbLayer", &corr);
    tree->Branch("peak2c","EdbPeak2", &peak2c);
    tree->Branch("peak",  &peak,"peak/I");
    tree->Branch("xoffset",  &xoffset,"xoffset/F");
    tree->Branch("yoffset",  &yoffset,"yoffset/F");
  }
  EdbScanSet *ss1 = ReadScanSet(idset1);   if(!ss1) return;
  EdbScanSet *ss2 = ReadScanSet(idset2);   if(!ss2) return;
  gStyle->SetPalette(1);
  int n = ss1->eIDS.GetSize();
  for(int i=0; i<n; i++) {
    EdbID *id1  = ss1->GetID(i);
    EdbID *id2  = ss2->FindPlateID(id1->ePlate);
    if(id2) UpdateAlignSummaryTree(*id1,*id2,*tree);
  }
  tree->AutoSave();
  f->Close();
}

//______________________________________________________________________________
void EdbScanProc::MakeAlignSetSummary(EdbID idset)
{
  // assuming that exist the scan sets for idset  
  // read id.n_id.n+1.aff.par and make a summary tree

  TString name;
  MakeFileName(name,idset,"align.pdf",false);
  Log(2,"MakeAlignSetSummary","%s",name.Data());
  gStyle->SetOptDate(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1001111);
  EdbScanSet *ss = ReadScanSet(idset);    if(!ss) return;
  int n = ss->eIDS.GetSize();          if(n<2)  return;
  for(int i=0; i<n-1; i++) {
    EdbID *id1  = ss->GetID(i);
    EdbID *id2  = ss->GetID(i+1);

    TString dataout;  MakeAffName(dataout,*id1,*id2,"al.root");
    TFile *f = new TFile(dataout,"READ");
    if(!f) continue;

    EdbLayer *corr   = (EdbLayer*)f->Get("corr_layer1");
    EdbPeak2 *peak2c = (EdbPeak2*)f->Get("peak2c");
    float xcenter1 = corr->X();
    float ycenter1 = corr->Y();
    EdbAffine2D *aXY=corr->GetAffineXY();
    float xoffset= aXY->A11()*xcenter1 + aXY->A12()*ycenter1 + aXY->B1() - xcenter1;
    float yoffset= aXY->A21()*xcenter1 + aXY->A22()*ycenter1 + aXY->B2() - ycenter1;

    Log(1,"UpdateAlignSummaryTree","peak: %7.0f/%7.2f/%7.3f   dx,dy,dz: %7.3f %7.3f %7.3f  for %s_%s", 
       peak2c->Peak(),peak2c->Mean3(),peak2c->Mean(), xoffset,yoffset,corr->Zcorr(), id1->AsString(),id2->AsString() );

    TCanvas *c = (TCanvas*)f->Get("report_al");
    if(c) {
      c->SetName(Form("%s_%s",id1->AsString(), id2->AsString()));
    //c->Draw();
      if(i==0&&n>2)         c->Print(Form("%s(",name.Data()),"pdf");
      else if(i==n-2&&n>2)  c->Print(Form("%s)",name.Data()),"pdf");
      else             c->Print(name,"pdf");
    }
    f->Close();
  }
}

//______________________________________________________________________________
void EdbScanProc::UpdateAlignSummaryTree(EdbID id1s, EdbID id2s, TTree &tree)
{
  TString dataout;  MakeAffName(dataout,id1s,id2s,"al.root");
  TFile *f = new TFile(dataout,"READ");
  if(!f) return;
  EdbLayer *corr = (EdbLayer*)f->Get("corr_layer1");
  EdbPeak2 *peak2c = (EdbPeak2*)f->Get("peak2c");
  if(!corr)   {Log(1,"UpdateAlignSummaryTree","Warning no corr for %s_%s",id1s.AsString(), id2s.AsString()); return;}
  if(!peak2c) {Log(1,"UpdateAlignSummaryTree","Warning no peak2c for %s_%s",id1s.AsString(), id2s.AsString()); return;}
  EdbAffine2D *aXY=corr->GetAffineXY();
  int peak = (int)(peak2c->Peak(0));
  float xcenter1 = corr->X();
  float ycenter1 = corr->Y();
  float xoffset= aXY->A11()*xcenter1 + aXY->A12()*ycenter1 + aXY->B1() - xcenter1;
  float yoffset= aXY->A21()*xcenter1 + aXY->A22()*ycenter1 + aXY->B2() - ycenter1;

  EdbID *id1=&id1s, *id2=&id2s;
  tree.SetBranchAddress("id1" ,  &id1);
  tree.SetBranchAddress("id2" ,  &id2);
  tree.SetBranchAddress("corr" ,  &corr);
  tree.SetBranchAddress("peak2c",  &peak2c);
  tree.SetBranchAddress("peak",  &peak);
  tree.SetBranchAddress("xoffset",  &xoffset);
  tree.SetBranchAddress("yoffset",  &yoffset);
  tree.Fill();
  Log(1,"UpdateAlignSummaryTree","peak: %7.0f/%7.2f/%7.3f   dx,dy,dz: %7.3f %7.3f %7.3f  for %s_%s", 
      peak2c->Peak(),peak2c->Mean3(),peak2c->Mean(), xoffset,yoffset,corr->Zcorr(), id1s.AsString(),id2s.AsString() );

  TCanvas *c = (TCanvas*)f->Get("report_al");
  if(c) {
    c->SetName(Form("%s_%s",id1s.AsString(), id2s.AsString()));
    c->Draw();
  }
  SafeDelete(f);
}

///______________________________________________________________________________
void EdbScanProc::MakeLinkSetSummary(EdbID idset)
{
  // assuming that exist the scan sets for idset  
  // read id.n_id.n+1.aff.par and make a summary tree
 
  TString name;
  MakeFileName(name,idset,"link.pdf",false);
 
  EdbScanSet *ss = ReadScanSet(idset);    if(!ss) return;
  int n = ss->eIDS.GetSize();          if(n<1)  return;
  gStyle->SetOptDate(1);
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1001111);
  for(int i=0; i<n; i++) {
    EdbID *id  = ss->GetID(i);
    TString dataout;  MakeFileName(dataout, *id, "cp.root");
    TFile *f = new TFile(dataout,"READ");
    if(!f) return;
    TCanvas *c = (TCanvas*)f->Get("report");
    if(c) {
      c->SetName(Form("%s",id->AsString()));
      if(i==0&&n>1)         c->Print(Form("%s(",name.Data()),"pdf");
      else if(i==n-1&&n>1)  c->Print(Form("%s)",name.Data()),"pdf");
      else                  c->Print(name,"pdf");
   }
    SafeDelete(f);
  }
  //tree->AutoSave();
  //f->Close();
}

//______________________________________________________________________________
void EdbScanProc::AlignRawSet(EdbID idset1, EdbID idset2, TEnv &cenv)
{
  // assuming that exist the scan sets for idset1 and idset2 
  // for each plate do the alignment raw from 1 to 2 and write id1_id2.aff.par
  EdbScanSet *ss1 = ReadScanSet(idset1);   if(!ss1) return;
  EdbScanSet *ss2 = ReadScanSet(idset2);   if(!ss2) return;
  int n = ss1->eIDS.GetSize();
  for(int i=0; i<n; i++) {
    EdbID *id1  = ss1->GetID(i);
    EdbID *id2  = ss2->FindPlateID(id1->ePlate);
    if(id2) AlignRaw(*id1,*id2,cenv);
  }
}

//______________________________________________________________________________
int EdbScanProc::AlignRaw(EdbID id1, EdbID id2, TEnv &cenv, EdbAffine2D *applyAff)
{
  // Align raw segments patterns. Typical application: align 2 different scannings of the same emulsion plate.

  int npat=0;
  Log(1,"\nAlignRaw","%s_%s", id1.AsString(), id2.AsString() );
  int       side1      = cenv.GetValue("fedra.alignRaw.side1"     , 1);
  int       side2      = cenv.GetValue("fedra.alignRaw.side2"     , 1);
  float        z1      = cenv.GetValue("fedra.alignRaw.Z1"        , 0);
  float        z2      = cenv.GetValue("fedra.alignRaw.Z2"        , 0);
  float     wmin1      = cenv.GetValue("fedra.alignRaw.Wmin1"     , 8);
  float     wmin2      = cenv.GetValue("fedra.alignRaw.Wmin2"     , 8);
  float      sigmaR    = cenv.GetValue("fedra.alignRaw.sigmaR"    , 2.5 );
  float      sigmaT    = cenv.GetValue("fedra.alignRaw.sigmaT"    , 0.005 );
  float      offsetMax = cenv.GetValue("fedra.alignRaw.offsetMax" , 500. );
  float      thetaMax  = cenv.GetValue("fedra.alignRaw.thetaMax"  ,   1. );
  int        path1     = cenv.GetValue("fedra.alignRaw.path1"      , -1 );
  int        path2     = cenv.GetValue("fedra.alignRaw.path2"      , -1 );

  float      xmin1  = cenv.GetValue("fedra.alignRaw.xmin1"  ,   -500. );
  float      xmax1  = cenv.GetValue("fedra.alignRaw.xmax1"  ,    500. );
  float      xmin2  = cenv.GetValue("fedra.alignRaw.xmin2"  ,   -500. );
  float      xmax2  = cenv.GetValue("fedra.alignRaw.xmax2"  ,    500. );
  float      ymin1  = cenv.GetValue("fedra.alignRaw.ymin1"  ,   -500. );
  float      ymax1  = cenv.GetValue("fedra.alignRaw.ymax1"  ,    500. );
  float      ymin2  = cenv.GetValue("fedra.alignRaw.ymin2"  ,   -500. );
  float      ymax2  = cenv.GetValue("fedra.alignRaw.ymax2"  ,    500. );

  float      DZ     = cenv.GetValue("fedra.alignRaw.DZ"     ,     25.   );
  float      DPHI   = cenv.GetValue("fedra.alignRaw.DPHI"   ,     0.003 );


  EdbRunAccess r1;  if(!InitRunAccess(r1,id1)) return 0;
  EdbRunAccess r2;  if(!InitRunAccess(r2,id2)) return 0;
  r1.eAFID = r2.eAFID = cenv.GetValue("fedra.alignRaw.AFID"   , 1);

  float min1[5] = {xmin1,ymin1,-thetaMax,-thetaMax, wmin1 };
  float max1[5] = {xmax1,ymax1, thetaMax, thetaMax, 100   };
  r1.AddSegmentCut(side1, 1, min1, max1);
  float min2[5] = {xmin2,ymin2,-thetaMax,-thetaMax, wmin2 };
  float max2[5] = {xmax2,ymax2, thetaMax, thetaMax, 100   };
  r2.AddSegmentCut(side2, 1, min2, max2);

  EdbPattern p1, p2;
  r1.GetPatternDataForPrediction( path1, side1, p1 );

  if(applyAff) p1.Transform(applyAff);

  if(path2>=0)  r2.GetPatternDataForPrediction( path2, side2, p2 );
  else  {
    float xmin=p1.Xmin(),  xmax=p1.Xmax();
    float ymin=p1.Ymin(),  ymax=p1.Ymax();
    EdbSegP s(0, (xmin+xmax)/2., (ymin+ymax)/2., 0,0);
    float dx=(xmax-xmin)/2., dy=(ymax-ymin)/2.;
    float size = Sqrt( dx*dx+dy*dy ) + offsetMax+200.;
    r2.GetPatternXY(s, side2, p2, size);
  }

  EdbPlateAlignment av;
  av.SetSigma(sigmaR,sigmaT);
  av.eOffsetMax = offsetMax;
  av.eDZ        = DZ;
  av.eDPHI      = DPHI;
  av.eDoFine = 1;
  TString dataout;  MakeAffName(dataout,id1,id2,"al.root");
  av.InitOutputFile( dataout );

  av.Align( p1, p2, z2-z1);
  av.CloseOutputFile();

  //av.eCorrL[0].Print();
  UpdateAFFPar( id1, id2, av.eCorrL[0] );
  return npat;
 }

//______________________________________________________________________________
int EdbScanProc::ReadMarksSet(EdbMarksSet &ms, int brick, const char *filename, char spacer, char shape)
{
  // Reads map file and copy its informations into an EdbMarksSet object
  // With "spacer" you can choose the character which is assumed to be
  // between the words in the map file. Default is '_'
  char str[256];
  sprintf(str,"%s/b%6.6d/b%6.6d.%s",
	  eProcDirClient.Data(),brick,brick,filename);

  ms.ReadMap(str,spacer,shape);

  return(1);
}

//______________________________________________________________________________
int EdbScanProc::WriteMarksSet(EdbMarksSet &ms, int brick, const char *filename, char spacer, char shape, int plate)
{
  // Reads an EdbMarksSet object and uses its content to write a map file
  // With "spacer" you can choose the character that you want to insert
  // between the words in the map file. Default is '_'
  char str[256];
  sprintf(str,"%s/b%6.6d/b%6.6d.%s",
	  eProcDirClient.Data(),brick,brick,filename);

  ms.WriteMap(str,spacer,shape,plate);

  return(1);
}

//______________________________________________________________________________
void EdbScanProc::LogPrint(int brick, int level, const char *location, const char *va_(fmt), ...)
{
// Print message to the logfile and to stdout.
  if(gEDBLOGFILE) {
    printf("WARNING in LogPrint! logfile seems to be opened. Trying to close it...\n");
    fclose(gEDBLOGFILE);
    gEDBLOGFILE=0;
  }
  char str[512];
  sprintf(str,"%s/b%6.6d/b%6.6d.log", eProcDirClient.Data(), brick,brick);
  gEDBLOGFILE = fopen(str,"a");
  if(!gEDBLOGFILE) printf("ERROR in LogPrint! can not open logfile: %s\n",str);
  
  va_list ap;
  va_start(ap,va_(fmt));
  Log0(level, location, va_(fmt), ap);
  va_end(ap);

  if(gEDBLOGFILE) fclose(gEDBLOGFILE);
  gEDBLOGFILE=0;
}

//--------------------------------------------------------------------
void EdbScanProc::GetPatternSide( EdbID id, int side, EdbLayer &la, const char *segcut, int afid, EdbPattern &p)
{
  int runside = 3-side;
  Log(2,"EdbScanProc::GetPatternSide","for id %s with cut %s", id.AsString(), segcut);
  EdbRunAccess r;
  TString runfile;
  MakeFileName(runfile,id,"raw.root");
  if( !r.InitRun(runfile) ) return;
  r.eAFID = afid;
  *(r.GetLayer(runside)) = la;
  r.AddSegmentCut( 1, segcut );
  r.GetLayer(runside)->Print();
  r.GetPatternDataForPrediction( -1, runside, p );
}

//--------------------------------------------------------------------
bool EdbScanProc::InitRunAccessNew(EdbRunAccess &r, EdbID idset, int idplate, bool do_update)
{
  EdbScanSet  *set = ReadScanSet(idset);
  EdbPlateP *plate = set->GetPlate(idplate);
  if(!plate) return 0;
  EdbID id = idset; id.ePlate = idplate;
  return InitRunAccessNew(r, id, *plate, do_update);
}

//--------------------------------------------------------------------
bool EdbScanProc::InitRunAccessNew(EdbRunAccess &r, EdbID id, EdbPlateP &plate, bool do_update)
{
  // use only scanset file (no *.par) 
  // before this function one should define variables:
  //             r.eInvertSides (default is 0 - no invert)
  // after:
  //             r.eAFID        (default is 1 - use view aff)
  //             r.AddSegmentCut(...) , etc
  TString runfile;
  MakeFileName(runfile,id,"raw.root");
  if( !r.InitRun(runfile, do_update) )
    {     
      LogPrint(id.eBrick,1,"InitRunAccess","ERROR open file %s !!!",runfile.Data());
      return false;
    } 
    else
      LogPrint(id.eBrick,2,"InitRunAccess"," %s with %d views",runfile.Data(), r.GetRun()->GetEntries() );
  r.GetLayer(2)->Copy( *(plate.GetLayer(1)) );
  r.GetLayer(1)->Copy( *(plate.GetLayer(2)) );
  r.GetLayer(0)->Copy( *((EdbLayer*)(&plate)) );
  //r.GetLayer(2)->Print();
  //r.GetLayer(1)->Print();
  //r.GetLayer(0)->Print();
  return true;
}

//--------------------------------------------------------------------
void EdbScanProc::LinkRunTest( EdbID id, EdbPlateP &plate, TEnv &cenv)
{
  EdbRunAccess r;
  r.eInvertSides=cenv.GetValue("fedra.link.read.InvertSides"      , 0);
  r.eHeaderCut = cenv.GetValue("fedra.link.read.HeaderCut"      , "1");
  r.eHeaderCut.Print();
  r.eAFID           =  cenv.GetValue("fedra.link.AFID"      , 1);
  printf("EdbScanProc::LinkRunTest ** AFID=%d\n", r.eAFID);
  InitRunAccessNew(r,id,plate);
  r.eWeightAlg  =  cenv.GetValue("fedra.link.read.WeightAlg"      , 0  );
  r.AddSegmentCut(1,cenv.GetValue("fedra.link.read.ICUT"      , "-1") );
  r.SetPixelCorrection( cenv.GetValue("fedra.link.PixelCorr"      , "0 1. 1.") );
  r.eTracking =  cenv.GetValue("fedra.link.Tracking"      , -1);

  EdbPattern p1, p2;
  p1.SetScanID(id); p1.SetSide(2);
  p2.SetScanID(id); p2.SetSide(1);
  r.GetPatternDataForPrediction( -1, 2, p1 );
  r.GetPatternDataForPrediction( -1, 1, p2 );

  EdbLinking link;
  TString cpfile;
  MakeFileName(cpfile,id,"cp.root");
  link.InitOutputFile( cpfile );

  if( cenv.GetValue("fedra.link.CheckUpDownOffset"      , 1) ) r.CheckUpDownOffsets()->Write();
  if(r.eDoViewAnalysis)   {
     r.eHViewXY[1].DrawH2("ViewXY1","XY segments distribution in a view coord side 1")->Write();
     r.eHViewXY[2].DrawH2("ViewXY2","XY segments distribution in a view coord side 2")->Write();
     r.CheckViewSize();
     //r.CheckStepSize();
  }

  link.Link( p2, p1, *(plate.GetLayer(2)), *(plate.GetLayer(1)), cenv );
  link.CloseOutputFile();
  if(link.eDoCorrectShrinkage || link.eDoCorrectAngles) {
     UpdatePlatePar( id, link.eL1 );  //TODO: check up/down id
     UpdatePlatePar( id, link.eL2 );
   }
}

//--------------------------------------------------------------------
void EdbScanProc::AlignOverlaps(EdbID id, EdbPattern &p1,EdbPattern &p2, TEnv &cenv, const char *suff)
{
  float      sigmaR    = cenv.GetValue("fedra.alignRaw.sigmaR"    , 0.7 );
  float      sigmaT    = cenv.GetValue("fedra.alignRaw.sigmaT"    , 0.019 );
  float      offsetMax = cenv.GetValue("fedra.alignRaw.offsetMax" , 10. );
  float      DZ        = cenv.GetValue("fedra.alignRaw.DZ"        , 10.   );
  float      DPHI      = cenv.GetValue("fedra.alignRaw.DPHI"      , 0.0015 );
  EdbPlateAlignment av;
  av.SetSigma(sigmaR,sigmaT);
  av.eOffsetMax = offsetMax;
  av.eDZ        = DZ;
  av.eDPHI      = DPHI;
  av.eDoFine = 1;

  TString dataout;  MakeFileName(dataout,id,suff);
  av.InitOutputFile( dataout );
  av.eSaveCouples=true;
  av.Align( p1, p2, 0);
  av.CloseOutputFile();
}

//--------------------------------------------------------------------
void EdbScanProc::CheckViewOverlaps(EdbID id, TEnv &cenv)
{
  float pulsMin1 = 6;
  float pulsMin2 = 6;
  EdbRunAccess r;  if(!InitRunAccess(r,id)) return;
  r.CheckViewStep();
  EdbPatternsVolume vol; vol.AddPattern(new EdbPattern()); vol.AddPattern(new EdbPattern());
  r.GetVolumeArea(vol, 1);
  r.CheckViewSize();

  EdbPattern p1l,p1r, p2l,p2r;
  r.SetCutLeft(1, pulsMin1);  r.GetPatternDataForPrediction( -1, 1, p1l );
  r.SetCutLeft(2, pulsMin1);  r.GetPatternDataForPrediction( -1, 2, p2l );
  r.SetCutRight(1, pulsMin1); r.GetPatternDataForPrediction( -1, 1, p1r );
  r.SetCutRight(2, pulsMin1); r.GetPatternDataForPrediction( -1, 2, p2r );

  EdbPattern p1t,p1b, p2t,p2b;
  r.SetCutTop(1, pulsMin2);     r.GetPatternDataForPrediction( -1, 1, p1t );
  r.SetCutTop(2, pulsMin2);     r.GetPatternDataForPrediction( -1, 2, p2t );
  r.SetCutBottom(1, pulsMin2);  r.GetPatternDataForPrediction( -1, 1, p1b );
  r.SetCutBottom(2, pulsMin2);  r.GetPatternDataForPrediction( -1, 2, p2b );

  AlignOverlaps(id, p1l,p1r, cenv, "al.vlr1.root");
  AlignOverlaps(id, p2l,p2r, cenv, "al.vlr2.root");
  AlignOverlaps(id, p1t,p1b, cenv, "al.vtb1.root");
  AlignOverlaps(id, p2t,p2b, cenv, "al.vtb2.root");
}

//--------------------------------------------------------------------
void EdbScanProc::LinkRunNew( EdbID id, EdbPlateP &plate, TEnv &cenv)
{
  TString rawfile, cpfile;
  MakeFileName(rawfile,id,"raw.root");
  MakeFileName(cpfile,id,"cp.root");
  EdbAlignmentMap amap( cpfile.Data(), "RECREATE");
  amap.eEnv = &cenv;
  amap.Link( rawfile.Data(), plate );
}

//----------------------------------------------------------------
void EdbScanProc::LinkSetNewTest(EdbScanSet &sc, TEnv &cenv )
{
  if(sc.eIDS.GetSize()<1) return;
  for(int i=0; i<sc.eIDS.GetSize(); i++) {
    EdbID *id = (EdbID *)(sc.eIDS.At(i));        if(!id)    continue;
    EdbPlateP *plate = sc.GetPlate(id->ePlate);  if(!plate) continue;
    LinkRunTest(*id, *plate, cenv);
  }
}

//----------------------------------------------------------------
void EdbScanProc::LinkSetNew(EdbScanSet &sc, TEnv &cenv )
{
  if(sc.eIDS.GetSize()<1) return;
  for(int i=0; i<sc.eIDS.GetSize(); i++) {
    EdbID *id = (EdbID *)(sc.eIDS.At(i));
    EdbPlateP *plate = sc.GetPlate(id->ePlate);
    LinkRunNew(*id, *plate, cenv);
  }
}

//----------------------------------------------------------------
int EdbScanProc::MakeTracksPred(TObjArray &tracks, EdbID id, EdbLayer &layer)
{
  // Extrapolate each track of array tracks to layer.Z(), apply layer.Aff() and prepare 
  // id.pred.root 

  int ntr = tracks.GetEntriesFast();
  Log(2,"MakeTracksPred","%d tracks",ntr );
  id.Print(); layer.Print();

  EdbPattern pred(0,0, layer.Z(), ntr);
  EdbSegP s;
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = (EdbTrackP*)tracks.At(i);
    t->MakePredictionTo( layer.Z(), s );
    pred.AddSegment(s);
  }
  EdbAffine2D aff(*(layer.GetAffineXY()));
  aff.Invert();
  pred.Transform( &aff );
  WritePred(pred,id);
  return pred.N();
}

//----------------------------------------------------------------
void EdbScanProc::ExtractRawVolume(EdbID id, EdbID idnew, EdbSegP pred, int plateid, TEnv &cenv)
{
  // from EdbScanSet id extract volume around pred and save into idnew
  // Note: the z of pred is assumed to be the z of plate in this function
  
  int   nplBefore =   cenv.GetValue("fedra.ExtractVol.NplBefore"  , 4    );
  int   nplAfter  =   cenv.GetValue("fedra.ExtractVol.NplAfter"   , 5    );
  float dR        =   cenv.GetValue("fedra.ExtractVol.DR"         , 1000.);
  Log(2,"ExtractRawVolume","%s --> %s  ref plate: %d (%d %d) at (%f %f)",
      id.AsString(),idnew.AsString(), plateid, plateid-nplBefore, plateid+nplAfter, pred.X(),pred.Y() );
  if(gEDBDEBUGLEVEL>2) cenv.Print();

  EdbScanSet *ss = ReadScanSet(id);   if(!ss) return;
  EdbScanSet ssnew; ssnew.Copy(*ss);
  ssnew.eB.SetID(idnew.eBrick);
  ssnew.eIDS.Clear();
  ssnew.SetID(idnew);

  EdbPlateP *pl = ss->GetPlate(plateid);
  if(!pl) Log(1,"ExtractRawVolume","ERROR: the plate %d is missing in scan set",plateid);
  Log(3,"ExtractRawVolume","assign to the prediction z = %f of the plate %d",pl->Z(),plateid );
  pred.SetZ(pl->Z());

  int n = ss->eIDS.GetSize();
  for(int i=0; i<n; i++ ) {
    EdbID *id = ss->GetID(i);
    if( (id->ePlate >= plateid-nplBefore) && (id->ePlate <= plateid+nplAfter) ) {
      EdbID *newid = new EdbID(*id);
      newid->eBrick = idnew.eBrick;
      newid->eMajor = idnew.eMajor;
      newid->eMinor = idnew.eMinor;
      ssnew.eIDS.Add( newid );
    }
  }
  ssnew.MakePIDList();

  ExtractRawVolume(*ss, ssnew, pred, dR);
}

//----------------------------------------------------------------
int EdbScanProc::FindRawTrack( EdbTrackP &pred,  EdbTrackP &found, EdbID idset, int plate)
{
  // found segments will be added to track tr
  TEnv env;    //env.SaveLevel(kEnvLocal);
  return FindRawTrack( pred,found,idset,plate,env );
}

//----------------------------------------------------------------
int EdbScanProc::FindRawTrack( EdbTrackP &pred,  EdbTrackP &found, EdbID idset, int plate, TEnv &env )
{
  // found segments will be added to track tr

  EdbRunTracking rt;
  rt.eDeltaRview           = env.GetValue( "fedra.RawTrack.DeltaRview",          700.       );
  rt.eDeltaTheta           = env.GetValue( "fedra.RawTrack.DeltaTheta",            0.15     );
  rt.eDeltaR               = env.GetValue( "fedra.RawTrack.DeltaR",                10.       );
  rt.ePreliminaryPulsMinMT = env.GetValue( "fedra.RawTrack.PreliminaryPulsMinMT",  4.       );
  rt.ePreliminaryChi2MaxMT = env.GetValue( "fedra.RawTrack.PreliminaryChi2MaxMT",  5.       );
  rt.ePulsMinMT            = env.GetValue( "fedra.RawTrack.PulsMinMT",             10.       );
  rt.eChi2MaxMT            = env.GetValue( "fedra.RawTrack.Chi2MaxMT",             2.6      );
  rt.ePulsMinBT            = env.GetValue( "fedra.RawTrack.PulsMinBT",            15.       );
  rt.eChi2MaxBT            = env.GetValue( "fedra.RawTrack.Chi2MaxBT",            2.5      );
  rt.eDegradPos            = env.GetValue( "fedra.RawTrack.DegradPos",             3.       );
  rt.eDegradSlope          = env.GetValue( "fedra.RawTrack.DegradSlope",           0.001    );
  rt.eAFID                 = env.GetValue( "fedra.RawTrack.AFID"       ,           1        );
  SetDefaultCondBT(rt.eCondBT);
  SetDefaultCondMT(rt.eCondMT);

  EdbScanSet *ss = ReadScanSet(idset);          if(!ss) return 0;
  EdbID      *id = ss->FindPlateID(plate);      if(!id) return 0;
  EdbPlateP  *pl = ss->GetPlate(plate);         if(!pl) return 0;
  
  TString runfile;
  MakeFileName(runfile,*id,"raw.root");
  if( !rt.InitRun(runfile) ) return 0;
  
  rt.GetLayer(2)->Copy( *(pl->GetLayer(2)) );
  rt.GetLayer(1)->Copy( *(pl->GetLayer(1)) );
  
  int status = rt.FindTrack(pred,found, *pl);
  //found.AddSegment( new EdbSegP( rt.eS1 ) );
  //found.AddSegment( new EdbSegP( rt.eS2 ) );
  return status;
}

//----------------------------------------------------------------
void EdbScanProc::ExtractRawVolume(EdbScanSet &ss, EdbScanSet &ssnew, EdbSegP &pred, float dR)
{
  // from EdbScanSet ss extract volume for plates defined in ssnew around pred with dR
  int n = ssnew.eIDS.GetSize();    if(n<1) return;
  if(!CheckBrickDir(ssnew.eID,true))          return;
  WriteScanSet(ssnew.eID,ssnew);

  for(int i=0; i<n; i++ ) {
    EdbID *newid = ssnew.GetID(i);
    newid->Print();
    EdbPlateP *plate = ssnew.GetPlate(newid->ePlate);
    EdbSegP newpred(pred);
    newpred.PropagateTo(plate->Z());

    EdbID *id = ss.FindPlateID(newid->ePlate);
    EdbRunAccess ra; if(!InitRunAccess(ra,*id)) return;   // TODO: check transformations used to select views!

    if(CheckProcDir(*newid,true)) {
      TString runfile;
      MakeFileName(runfile,*newid,"raw.root");
      ra.CopyRawDataXY( newpred.X(), newpred.Y(), dR,runfile.Data() );
    }
  }

}

//----------------------------------------------------------------
void EdbScanProc::SetServerRunName(const char* fname_){
	eServerCreatedRunName = fname_;
};
const char* EdbScanProc::GetServerRunName()const{
	return eServerCreatedRunName.Data();
};

//----------------------------------------------------------------
void EdbScanProc::ReadUncorrectedBTforFoundTracks( EdbPVRec &ali )
{
  // get from couples trees and set in the ali the uncorrected basetracks
  ali.SetSegmentsTracks();
  int npat = ali.Npatterns();
  for(int i=0; i<npat; i++)
  {
    EdbPattern *pat = ali.GetPattern(i);
    int nseg = pat->N();    if(nseg<1) return;
    EdbID idpl((pat->GetSegment(0)->ScanID()));
    TString cpfile;
    MakeFileName(cpfile,idpl,"cp.root");
    EdbCouplesTree ect;
    if(!ect.InitCouplesTree("couples",cpfile,"READ")) continue;
    ect.eAcceptMask = new EdbMask(ect.eTree->GetEntries());
    EdbSegP **segments = new EdbSegP*[ect.eTree->GetEntries()];
    for(int j=0; j<nseg; j++)       ect.eAcceptMask->SetAt( pat->GetSegment(j)->Vid(1), 1);
    EdbPattern newpat;
    int nnew= ect.GetCPDataAcceptedMask(&newpat);
    for(int j=0; j<nnew; j++)  segments[newpat.GetSegment(j)->Vid(1)] = newpat.GetSegment(j);
    for(int j=0; j<nseg; j++) {
      EdbSegP *sc = pat->GetSegment(j);
      EdbSegP *sn = segments[sc->Vid(1)];
      //sc->PrintNice();
      //sn->PrintNice();
      sc->SetX(sn->X());
      sc->SetY(sn->Y());
      sc->SetTX(sn->TX());
      sc->SetTY(sn->TY());
    }
    Log(2,"EdbScanProc::ReadUncorrectedBTforFoundTracks","%d -> %d segments for %s", 
        pat->N(), newpat.N(), idpl.AsString() );
  }

}
