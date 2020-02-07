//-- Author :  Valeri Tioukov & Yury Petukhov  10.02.2004
 
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbVertex                                                            //
//                                                                      //
// Class for vertex reconstruction                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TFile.h"
#include "TRandom.h"
#include "TMap.h"
#include <TArrayF.h>

#include "vt++/CMatrix.hh"
#include "vt++/VtVector.hh"
#include "vt++/VtDistance.hh"
#include "smatrix/Functions.hh"
#include "smatrix/SVector.hh"

#include "EdbVertex.h"
#include "EdbPVRec.h"
#include "EdbPhys.h"
#include "EdbMath.h"
#include "TIndexCell.h"
#include "EdbLog.h"

ClassImp(EdbVTA);
ClassImp(EdbVertexPar);
ClassImp(EdbVertex);
ClassImp(EdbVertexRec);

using namespace TMath;
using namespace MATRIX;
using namespace VERTEX;

//______________________________________________________________________________
EdbVTA::EdbVTA()
{
  Set0();
}

//______________________________________________________________________________
EdbVTA::EdbVTA( EdbVTA& vta )
{
  eTrack=vta.GetTrack();
  eVertex=vta.GetVertex();
  eZpos=vta.Zpos();
  eFlag=vta.Flag();
  eImp=vta.Imp();
  eDist=vta.Dist();
}

//______________________________________________________________________________
EdbVTA::EdbVTA( EdbTrackP *tr, EdbVertex *v )
{
  Set0();
  eTrack=tr;
  eVertex=v;
}

//______________________________________________________________________________
void EdbVTA::Set0()
{
  eTrack=0;
  eVertex=0;
  eZpos=0;
  eFlag=0;
  eImp=0.;
  eDist=0.;
}
//________________________________________________________________________
EdbVTA::~EdbVTA()
{
  if(eTrack && eFlag == 2) { eTrack->ClearVTA(this); }
  if(eVertex && (eFlag == 2)) { (eVertex->VTa())->Remove(this); }
  if(eVertex && (eFlag < 2)) { (eVertex->VTn())->Remove(this); }
}

//______________________________________________________________________________
void EdbVTA::Print()
{
  printf("VTA with: track=%5d  vertex= %5d  zpos =%2d  flag =%2d  imp =%8.4f  dist = %f\n", 
	 eTrack->ID(), eVertex->ID(), eZpos, eFlag, eImp, eDist );
}

//______________________________________________________________________________
void EdbVTA::AddVandT()
{
    GetTrack()->AddVTA(this);
    GetVertex()->AddVTA(this);
}
//______________________________________________________________________________
EdbVertex::EdbVertex()
{
  eID= 0;
  eV = 0;
  eX = 0.;
  eY = 0.;
  eZ = 0.;
  eFlag = 0;
  eQuality=0.;
  eMCEvt=-999;
}

//________________________________________________________________________
EdbVertex::~EdbVertex()
{
  Clear();
}

//________________________________________________________________________
void EdbVertex::ClearV()
{
  //clear VtVertex object
  if (eV) {
    eV->clear();     // should delete also tracks ownered by vertex
    SafeDelete(eV);
  }
}

//________________________________________________________________________
void EdbVertex::Clear()
{
  eVTa.Delete("slow");
  eVTn.Delete("slow");
  ClearV();
  eX = 0.;
  eY = 0.;
  eZ = 0.;
  eID = 0;
  eFlag = 0;
  eQuality = 0.;
  eMCEvt=-999;
}

//________________________________________________________________________
EdbVTA   *EdbVertex::GetMaxImpVTA()
{
// return vta with a biggest impact par
  int ntr = N();
  EdbVTA *maxvta=0;
  if(ntr<1)       maxvta=0;
  else if(ntr==1) maxvta=GetVTa(0);
  else {
    float maximp=0;
    for(int i=0; i<ntr; i++) {
      EdbVTA *vta = GetVTa(i);
      float   imp = vta->Imp();
      if( imp > maximp )  { maximp=imp; maxvta=vta; }
    }
  }
  return maxvta;
}

//________________________________________________________________________
EdbSegP   *EdbVertex::GetTrackV(int i, bool usesegpar)
{
	EdbTrackP *t = GetTrack(i); if(!t) return 0;
	EdbSegP *s = t->TrackExtremity(Zpos(i), usesegpar);
  if( s->P()>=0 && (s->P() != t->P()) ) Log(1,"GetTrackV","Warning! segment momentum=%f is not equal to the track momentum=%f",s->P(), t->P() );
	return s;
}

//______________________________________________________________________________
Int_t EdbVertex::CheckDiscardedTracks()
{
  int ndsc = 0;
  EdbVTA *vta = 0;
  for (int i=0; i<N(); i++)
  {
    vta = GetVTa(i);
    if (vta) {
      if( GetTrack(i)->Vertex(vta->Zpos()) != this ) ndsc++;
    }
    else ndsc++;
  }
  return ndsc;
}

//------------------------------------------------------------------------------
Int_t EdbVertex::EstimateVertexFlag()
{
  int flag=-1;
  int n0 = 0, n1=0, nn=0;
  EdbVTA *vta = 0;
  for (int i=0; i<N(); i++)
  {
    vta = GetVTa(i);
    if     ( vta->Zpos()==0 ) n0++;
    else if( vta->Zpos()==1 ) n1++;
    else                      nn++;
  }

  if      ( n0>0 && n1 >0)    flag=1;         // end   & start
  else if ( n0==0 && n1>0)    flag=0;         // start & start
  else if ( n0>0 && n1==0)    flag=2;         // end   & end
  else   flag= -1;
  if(nn) flag= -1;              // flag -1 should newer happened: if so - should be debugged

  SetFlag(flag);
  Log(3,"EdbVertex::EstimateVertexFlag","%d   with n0=%d n1=%d, nn=%d", flag, n0,n1,nn );
  return flag;
}

//______________________________________________________________________________
void EdbVertex::ResetTracks()
{
  // Assign the eVTAS or eVTAE of the tracks to the current vertex
  EdbVTA *vta = 0;
  for (int i=0; i<N(); i++)
  {
    vta = GetVTa(i);
    if (vta) GetTrack(i)->AddVTA(vta);
  }
}

//______________________________________________________________________________
Bool_t EdbVertex::TrackInVertex( EdbTrackP *t )
{
  int   ntr = N();
  if(!ntr) return 0;
  for (int i=0; i<ntr; i++)  if(t == GetTrack(i)) return 1;
  return 0;
}

// //______________________________________________________________________________
// float EdbVertex::MaxImpact()
// {
//   int   ntr = N();  if(ntr<2)                      return 0;
//   float imp=0.;
//   for (int i=0; i<ntr; i++)    if(Impact(i)>imp) imp=Impact(i);
//   return imp;
// }

//______________________________________________________________________________
Float_t EdbVertex::MinDist()
{
  float mind=999999999999.;
  int   ntr = N();
  if(ntr<2)                      return mind;
  for (int i=0; i<ntr; i++) {
    float d = GetVTa(i)->Dist();
    if(Abs(d)<mind) mind = Abs(d);
  }
  return mind;
}

//______________________________________________________________________________
float EdbVertex::MaxAperture()
{
  float aper=0.;
  int   ntr = N();
  if(ntr<2)                      return aper;
  EdbTrackP *t1=0;
  EdbTrackP *t2=0;

  float tx=0,ty=0,a=0;
  for (int i=0; i<ntr-1; i++) {
    t1 = GetTrack(i);
    for (int j=i+1; j<ntr; j++) {
      t2 = GetTrack(j);

      tx= t1->TX() - t2->TX();
      ty= t1->TY() - t2->TY();
      a = TMath::Sqrt( tx*tx+ty*ty );
      if( a>aper) aper=a;
    }
  }
  return aper;
}
//________________________________________________________________________
void EdbVertex::ClearNeighborhood()
{
  int nn = Nn();
  if (nn>0) eVTn.Clear("nodelete");
  for(int i=0; i<nn; i++) delete GetVTn(i);
}
//________________________________________________________________________
int EdbVertex::Compare(const TObject *o) const
{
    /*printf("Inside compare\n");*/
    if      ( eQuality >  ((EdbVertex *)o)->eQuality )  return -1;
    else if ( eQuality == ((EdbVertex *)o)->eQuality )  return  0;
    else	      	       				return  1;
}
//________________________________________________________________________
bool EdbVertex::IsEqual(const TObject *o) const
{
    /*printf("Inside isequal\n");*/
    if   ( eID      != ((EdbVertex *)o)->eID )      return false;
    if   ( eQuality != ((EdbVertex *)o)->eQuality ) return false;
    else      	             			    return true;
}

//________________________________________________________________________
EdbTrackP *EdbVertex::MeanTrack()
{
  // calculate mean track trajectory waited with momentum
  int ntr = N();
  EdbTrackP *mean = new EdbTrackP();
  float psum=0;
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = GetTrack(i);
    mean->Set( 0 , mean->X()  + t->X()*t->P()
                 , mean->Y()  + t->Y()*t->P()
                 , mean->TX() + t->TX()*t->P()
                 , mean->TY() + t->TY()*t->P()
                 , mean->W()  + t->W()*t->P()
                 , 0 );
    mean->SetZ( mean->Z()  + t->Z()*t->P() );
    psum += t->P();
  }
  mean->Set( 0 , mean->X()  / psum
               , mean->Y()  / psum
               , mean->TX() / psum
               , mean->TY() / psum
               , mean->W()  / psum
               , 0 );
  mean->SetZ( mean->Z() / psum );
  mean->SetP(psum);
  return mean;
}

//________________________________________________________________________
void EdbVertex::Print()
{
	int ntr = N();
	printf( "\n********************** Vertex %d  with flag %d   and %d tracks **************************\n", ID(), Flag(), ntr );
	printf( "Fit quality     : probability = %f    Chi2 = %f\n", V()->prob(), V()->chi2() );
	printf( "Vertex Position : %12.3f %12.3f %12.3f \n", VX(), VY(), VZ() );
	printf( "Position Errors : %12.3f %12.3f %12.3f\n", V()->vxerr(), V()->vyerr(), V()->vzerr() );
	//printf("---------------------------------------------------------\n");
	printf( "Track   ID   Nseg   Mass       P       Chi2/ndf    Prob     Chi2Contrib   Impact         Ztr\n");
	for(int i=0; i<ntr; i++) {
		EdbTrackP *tr = GetTrack(i);
		float ztr =  tr->TrackExtremity( Zpos(i) )->Z();
		printf("%4d  %4d  %4d   %7.4f  %7.2f    %5.2f     %7.4f      %7.3f   %8.2f   %10.2f\n",
		       i, tr->ID(), tr->N(), tr->M(), tr->P(),
		       tr->Chi2()/tr->N(), tr->Prob(), V()->track_chi2(i), Impact(i), ztr );
	}
	/*
	printf("------- mean \"jet\" direction: ---------\n");
	EdbTrackP *mean = MeanTrack();
	printf("         X            Y            Z       TX       TY         P\n");
	printf("%12.1f %12.1f %12.1f %8.4f %8.4f %10.2f\n", 
	        mean->X(), mean->Y(), mean->Z(), mean->TX(), mean->TY(), mean->P() );
	SafeDelete(mean);
	*/
	printf("****************************************************************************************\n");
}
//________________________________________________________________________
void EdbVertex::AddVTA(EdbVTA *vta)
{
  if(vta->Flag()!=2) eVTn.Add(vta);
  else               eVTa.Add(vta);
}
//________________________________________________________________________
void EdbVertex::RemoveVTA(EdbVTA *vta)
{
  if(vta->Flag()!=2) eVTn.Remove(vta);
  else               eVTa.Remove(vta);
}
//________________________________________________________________________
int EdbVertex::Nv()
{
  // return the number of linked vertex: tracks attached by the other edge to the different vertex
  EdbTrackP *tr = 0;
  EdbVertex *vc = 0;
  int nv = 0;
  for (int i=0; i<N(); i++)
    {
      if ((tr = GetTrack(i))) {
	if      ((Zpos(i) == 1) && (vc = tr->VertexE())) nv++;
	else if ((Zpos(i) == 0) && (vc = tr->VertexS())) nv++;
      }
    }
  return nv;
}

//________________________________________________________________________
EdbVertex *EdbVertex::GetConnectedVertexForTrack(int it)
{
    EdbTrackP *tr = 0;
    EdbVertex *vc = 0;
    if (it<N())
    {
	if ((tr = GetTrack(it)))
	{
	    if      ((Zpos(it) == 1) && (vc = tr->VertexE()))
	    {
		return vc;
	    }
	    else if ((Zpos(it) == 0) && (vc = tr->VertexS()))
	    {
		return vc;
	    }
	}
    }
    return 0;
}
//________________________________________________________________________
EdbVertex *EdbVertex::GetConnectedVertex(int nv)
{
    EdbTrackP *tr = 0;
    EdbVertex *vc = 0;
    int n = 0;
    for (int i=0; i<N(); i++)
    {
	if ((tr = GetTrack(i)))
	{
	    if      ((Zpos(i) == 1) && (vc = tr->VertexE()))
	    {
		if (n == nv) return vc;
		n++;
	    }
	    else if ((Zpos(i) == 0) && (vc = tr->VertexS()))
	    {
		if (n == nv) return vc;
		n++;
	    }
	}
    }
    return 0;
}

//________________________________________________________________________
float EdbVertex::CheckImpGeom( const EdbTrackP *tr )
{
    float pv[3] = {VX(), VY(), VZ()};
    float p1[3] = { tr->X(), tr->Y(), tr->Z() };
    float p2[3] = { tr->X() + tr->TX()*1000., tr->Y() + tr->TY()*1000., tr->Z()+1000. };
    bool inside=0;
    float imp = EdbMath::DistancePointLine3(pv, p1,p2, &inside);
    return imp;
}

//________________________________________________________________________
float EdbVertex::CheckImp( const EdbTrackP *tr )
{
    Track *t = new Track();
    Vertex *v = this->V();
    Edb2Vt(*tr, *t);
    return distance(*t,*v);
}

//________________________________________________________________________
EdbVTA *EdbVertex::CheckImp( const EdbTrackP *tr , float ImpMax, int zpos, float dist)
{
    EdbVTA *vta = 0;
    EdbTrackP *tr1 = (EdbTrackP *)tr;
    if (!tr) return vta;
    Track *t = new Track();
    Vertex *v = this->V();
    Edb2Vt(*tr, *t);
    float imp = distance(*t,*v);
    if (imp > ImpMax) { delete t; return vta;}
    vta = new EdbVTA(tr1, this);
    vta->SetZpos(zpos);
    vta->SetFlag(0);
    vta->SetImp(imp);
    vta->SetDist(dist);
    AddVTA(vta);
    delete t;
    return vta;
}

//________________________________________________________________________
bool EdbVertex::EstimateVertexMath( float& xv, float& yv, float& zv, float& d )
{
  int nt = N(); 
  if(!nt) return false;

  double tx_sum  = 0.;
  double x_sum   = 0.;
  double xw_sum  = 0.;
  double xtx_sum = 0.;
  double tx2_sum = 0.;
  double ty_sum  = 0.;
  double y_sum   = 0.;
  double yw_sum  = 0.;
  double yty_sum = 0.;
  double ty2_sum = 0.;
  const EdbTrackP  *tr = 0;
  const EdbSegP    *seg = 0;

  double x,y,tx,ty,xweight,yweight,xweight2,yweight2;

  // fill cumulants
  for( int i = 0; i < nt; i++ ) {
    
    tr = GetTrack(i);

    seg = tr->TrackExtremity( Zpos(i)); // usesegpar??

    x        = seg->X();
    y        = seg->Y();
    tx       = seg->TX();
    ty       = seg->TY();
    xweight  = 1./(seg->COV())(0,0);
    yweight  = 1./(seg->COV())(1,1);
    xweight2 = xweight*xweight;
    yweight2 = yweight*yweight;

    tx_sum  += tx * xweight;
    x_sum   += x  * xweight;
    xw_sum  += xweight;
    xtx_sum += x  * tx * xweight2;
    tx2_sum += tx * tx * xweight2;

    ty_sum  += ty * yweight;
    y_sum   += y  * yweight;
    yw_sum  += yweight;
    yty_sum += y  * ty * yweight2;
    ty2_sum += ty * ty * yweight2;

  } // for track

  double det = -tx2_sum - ty2_sum + tx_sum*tx_sum/xw_sum + ty_sum*ty_sum/yw_sum;

  if(det == 0.) {
    return false;
  }

  zv = ( xtx_sum + yty_sum - tx_sum*x_sum/xw_sum - ty_sum*y_sum/yw_sum ) / det;
  xv = ( x_sum + tx_sum * zv ) / xw_sum;
  yv = ( y_sum + ty_sum * zv ) / yw_sum;


//   float zTolerance=300.;

//   for( int i = 0; i < nt; i++ ) {
//       tr = GetTrack(i);
//       if (Zpos(i)) seg = tr->TrackZmin();
//       else	   seg = tr->TrackZmax();
//       if( zv > (seg->Z() + zTolerance) )  return false;
//       if( zv < (seg->Z() - zTolerance) )  return false;
//   }

  double drx;
  double dry;
  double drz;
  double drt;
  double drms = 0.;

  for( int i = 0; i < nt; i++ ) {

    tr = GetTrack(i);

     seg = tr->TrackExtremity( Zpos(i)); // usesegpar??
    
    drx = seg->X() - xv;
    dry = seg->Y() - yv;
    drz = seg->Z() - zv;
    drt = (drx*seg->TX() + dry*seg->TY() + drz);
    drms += drx*drx + dry*dry + drz*drz -
      (drt*drt)/(1.+seg->TX()*seg->TX()+seg->TY()*seg->TY());
  }

  d = TMath::Sqrt(drms/nt);

  return true;
}

//________________________________________________________________________
void EdbVertex::Edb2Vt(const EdbTrackP& tr, Track& t, float X0, float m)
{
  Edb2Vt( *((EdbSegP*)&tr), t, X0, m);
}

//________________________________________________________________________
void EdbVertex::Edb2Vt(const EdbSegP& tr, Track& t, float X0, float m)
{
  // Input: EdbSegP tr - track parameters near vertex
  //        X0         - rad length for ms estimation
  //        m          - mass of the particle
  //          if X0 or m are negative - ignore multiple scattering
  // Output: VERTEX:Track t - object propagated to the vertex position with the estimated errors matrix
  // Used: eX,eY,eZ - the reference point of this vertex

  double dz = eZ - tr.Z();
  double tx = (double)tr.TX();
  double ty = (double)tr.TY();
  double x  = (double)tr.X() + tx*dz - eX;
  double y  = (double)tr.Y() + ty*dz - eY;
  double z  = 0.;
  float  p  = tr.P();

  VtSymMatrix dms(4);   // multiple scattering matrix
  dms.clear();
  if ( X0 > 0. && m > 0.)
    {
      double dPb = dz*TMath::Sqrt(1.+tx*tx+ty*ty); // thickness of the Pb+emulsion cell in microns
      double theta0sq = EdbPhysics::ThetaMS2( p, m, dPb, X0 );
      //printf("( p, m, dPb, X0 dz) = %f %f %f %f %f theta0 = %g\n", p, m, dPb, X0, dz, TMath::Sqrt(theta0sq));
      dms(0,0) = theta0sq*dz*dz/3.;
      dms(1,1) = dms(0,0);
      dms(2,2) = theta0sq;
      dms(3,3) = dms(2,2);
      dms(2,0) = theta0sq*dz/2.;
      dms(3,1) = dms(2,0);
      dms(0,2) = dms(2,0);
      dms(1,3) = dms(2,0);
    }

  VtSqMatrix pred(4);        //propagation matrix for track parameters (x,y,tx,ty)
  pred.clear();
  pred(0,0) = 1.;
  pred(1,1) = 1.;
  pred(2,2) = 1.;
  pred(3,3) = 1.;
  pred(0,2) = dz;
  pred(1,3) = dz;

  VtSymMatrix cov(4);             // covariance matrix for seg0
  for(int k=0; k<4; k++) 
    for(int l=0; l<4; l++) cov(k,l) = (tr.COV())(k,l);

  VtSymMatrix covpred(4);         // covariance matrix for prediction
  covpred = pred*(cov*(pred.T()))+dms;

  CMatrix covp;             // covariance matrix for the track
  covp.clear();
  for(int k=0; k<4; k++) 
    for(int l=0; l<4; l++) covp(k,l) = covpred(k,l);
  covp(4,4) = (tr.COV())(4,4);

  t.set(x,  y,  z,  tx,  ty,  (double)p, covp);
  t.rm((double)m);
}

//________________________________________________________________________
float EdbVertex::Chi2Track(EdbTrackP *track, int zpos, float X0)
{
  // Chi2-distance from track to already existing vertex

  if (!track) return 0.;
  double distchi2 = -2.;
  EdbSegP *seg = track->TrackExtremity( zpos); // usesegpar??
  if (eV)
    {
      if (track->NF() <= 0) return -1.;
      if (track->P() <= 0.) track->SetP(1.);
      if (track->M() <= 0.) track->SetM(.1395);
      if (track->SP() <= 0.) track->SetErrorP(1.);
      Track *t=new Track();
      Edb2Vt( *seg, *t, X0, track->M() );
      distchi2 = -3.;
      if (eV->valid())
	{
	  t->rm(track->M());
	  distchi2 = eV->distance(*t);
	} 
      SafeDelete(t);
    }
  return (float)distchi2;
}

//________________________________________________________________________
float EdbVertex::DistTrack(EdbTrackP *track, int zpos, float X0)
{
  // distance from track to already fitted vertex
  if (!track) return 0.;
  EdbSegP *seg = track->TrackExtremity( zpos); // usesegpar??
  if (!seg)  return 0.;
  return DistSeg(seg,X0);
}

//________________________________________________________________________
float EdbVertex::DistSeg(EdbSegP *seg, float X0)
{
  // distance from segment to already fitted vertex
  if(seg)   
    if(eV)  
      if(eV->valid()) {
	Track t;
	Edb2Vt( *seg, t, X0, 0. );
	return (float)distance(t, *eV);
      }
  return 100000.;
}

//========================================================================
EdbVertexPar::EdbVertexPar()
{
  eZbin       = 100.;      // microns
  eAbin       = 0.01;      // rad
  eDZmax      = 3000.;     // microns
  eProbMin    = 0.01;      // i.e 1%
  eImpMax     = 50.;       // microns
  eImpMaxV    = 10.;       // microns
  eUseMom     = false;     // do not use it
  eUseSegPar  = false;     // use fitted track parameters
  eUseKalman  = true;      // use or not Kalman for the vertex fit. Default is true
  eUseLimits  = false;     // if true - look for the vertex only inside limits defined by eVmin:eVmax, default is false
}

//========================================================================
void EdbVertexRec::Set0()
{
  eEdbTracks = 0;
  eVTX       = 0;
  ePVR       = 0;
  eVertex    = 0;
  eWorking   = 0;
  (gROOT->GetListOfSpecials())->Add(this);    // To check if this can cause memory leak!
}

//________________________________________________________________________
EdbVertexRec::~EdbVertexRec()
{
  if ((gROOT->GetListOfSpecials())->FindObject(this))
  {
    (gROOT->GetListOfSpecials())->Remove(this);
  }
  Reset();
}

//________________________________________________________________________
void EdbVertexRec::Reset()
{
  SafeDelete(eVTX);
  eVTA.Clear("nodelete");
}

//________________________________________________________________________
Bool_t EdbVertexRec::EstimateVertexQuality( EdbVertex &vtx )
{
  // TODO! razobratsia s etimi qualitiami!!

  float quality=0.;
  vtx.SetQuality(quality);
  Vertex *v = vtx.V();
  if(!v)           return 0;
  if(!v->valid())  return 0;

  if (eQualityMode == 0)  quality = v->prob()/(v->vtx_cov_x()+v->vtx_cov_y());
  else if (eQualityMode == 1)    
    {
      double rms=v->rmsDistAngle();
      if (rms != 0.)  quality= (float)(1./rms);
      else quality =  10.e+35;
    }
  else quality= 1.;
  
  vtx.SetQuality(quality);
  return 1;
}

//________________________________________________________________________
Bool_t EdbVertexRec::EstimateVertexPosition( EdbVertex &v )
{
  // make approximate (without matrix) estimation of the vertex and set the reference point of 
  // the vertex XYZ in the estimated position

  int nt = v.N(); 
  if(nt<2) return false;

  EdbSegP *s1=0,*s2=0;
  bool zpos1,zpos2;
  float vxyz[3], vsum[3];
  for(int i=0; i<3; i++) vsum[i]=0;
  int count=0;
  for(int i1=0; i1<nt-1; i1++) {
    s1    = v.GetTrackV(i1,eUseSegPar);
    zpos1 = v.Zpos(i1);
    for(int i2=1; i2<nt; i2++) {
      s2    = v.GetTrackV(i2,eUseSegPar);
      zpos2 = v.Zpos(i2);
      if( CheckImpact(s1,s2,zpos1,zpos2, vxyz) > 2.*eImpMax) continue;
      for(int i=0; i<3; i++) vsum[i]+=vxyz[i];
      count++;
    }
  }
  if(count)  for(int i=0; i<3; i++) vsum[i]/=count;
  else {                    // take just mean tracks position
    for(int i=0; i<nt; i++) {
      EdbSegP *s  = v.GetTrackV(i,eUseSegPar);
      vsum[0] += s->X(); vsum[1] += s->Y(); vsum[2] += s->Z();
    }
    for(int i=0; i<3; i++) vsum[i]/=nt;
  }
  v.SetXYZ(vsum[0],vsum[1],vsum[2]);
  Log(3,"EdbVertexRec::EstimateVertexPosition","%f %f %f",vsum[0],vsum[1],vsum[2]);
  return true;
}

//________________________________________________________________________
int EdbVertexRec::MakeV( EdbVertex &edbv, bool isRefit )
{
  // create new VtVertex and add tracks to this one
  // if isRefit - use input vertex position to improve the fit (default is false)

  if (ePVR) if (ePVR->IsA() != EdbPVRec::Class()) ePVR = 0;
  if(!ePVR) {Log(1,"EdbVertexRec::MakeV","ERROR: EdbPVRec not defined, use SetPVRec(...)"); return 0;}
  
  int n = edbv.N();
  if (n<2)         return 0;

  if(isRefit) edbv.SetXYZ(edbv.VX(),edbv.VY(),edbv.VZ());
  else        EstimateVertexPosition(edbv);
  edbv.ClearV();
  Vertex *v = new Vertex();
  v->use_kalman(eUseKalman);
  v->use_momentum(eUseMom);

  float X0 =  ePVR->GetScanCond()->RadX0();
  if(!eUseMom) X0 = -1.;  // ignore multiple scattering contribution if eUseMom is false

  EdbSegP *seg=0;
  for (int i=0; i<n; i++)
    {
      seg = edbv.GetTrackV(i,eUseSegPar);
      Track *t = new Track();
      //seg->PrintNice();
      //printf("%f\n",edbv.Z());
      edbv.Edb2Vt(*seg, *t, X0, edbv.GetTrack(i)->M());
      v->add_track(*t);
    }
  if (!v->findVertexVt())       { Log(3,"MakeV","can not find VtVertex" ); return 0; }
  if (!(v->valid()))            { Log(3,"MakeV","VtVertex is not valid" ); return 0; }
  edbv.SetV(v);
  for(int i=0; i<n; i++)
    {
      edbv.GetVTa(i)->SetDist( edbv.VZ() -  edbv.GetTrackV(i,eUseSegPar)->Z() );
      edbv.GetVTa(i)->SetImp( distance(v->get_track(i),*v) );
//      printf("dist = %f imp = %f\n", edbv.GetVTa(i)->Dist(), edbv.GetVTa(i)->Imp() );
//      edbv.GetTrack(i)->PrintNice();
    }
  EstimateVertexQuality(edbv);
  Log( 3,"EdbVertexRec::MakeV","impmax = %f",edbv.MaxImpact() );
  return 1;
}

//________________________________________________________________________
EdbVertex *EdbVertexRec::StripBadTracks( EdbVertex &vtx, float impMax, int ntrMin )
{
  EdbVertex *v      = &vtx;
  int        ntr0    = v->N();
  int        ntr    = v->N();
  while( v && ntr>ntrMin && v->MaxImpact()>impMax ) 
  {
    Log(2,"EdbVertexRec::StripBadTracks"," ntr = %d     maximp = %f", ntr,v->MaxImpact() );
    EdbVTA *vta = v->GetMaxImpVTA();
    v = RemoveVTAFromVertex( *v, *vta );
    vta->SetFlag(0);  v->AddVTA(vta);    // add it as VTn
    ntr = v->N();
  }  
  
  Log(2,"EdbVertexRec::StripBadTracks","  %d -> %d + %d     maximp = %f", ntr0, v->N(), v->Nn() ,v->MaxImpact() );
  return v;
}

//________________________________________________________________________
EdbVertex *EdbVertexRec::Make1Vertex(TObjArray &tracks, float zexpected)
{
 // make a single vertex using tracks array
  
  EdbVertex *v = new EdbVertex();
  v->SetXYZ( 0,0, zexpected );
  int ntr = tracks.GetEntries();
  for(int i=0; i<ntr; i++) {
    EdbTrackP *t = (EdbTrackP*)tracks.At(i);
    EdbVTA *vta = new EdbVTA(t,v);
    vta->SetFlag(2);
    v->AddVTA(vta);
    (t->Z() >= v->VZ())? vta->SetZpos(1) : vta->SetZpos(0);
    t->AddVTA(vta);
  }
  if( MakeV(*v) )  AddVertex(v);
  else { SafeDelete(v); return 0; }                                // vertex is not valid
  return v;
}

//________________________________________________________________________
EdbVTA *EdbVertexRec::AddTrack( EdbVertex &edbv, EdbTrackP *track, int zpos )
{
  // add track to already existing vertex if prob > eProbMin
  // if vertex do not exist yet - calculate medium x,y,z

  if (!track)        return 0;
  EdbSegP *seg = track->TrackExtremity( zpos, eUseSegPar);
  if (!seg)          return 0;
  int ntr = edbv.N();
  if(edbv.TrackInVertex(track)) return 0;
  //for (int i=0; i<ntr; i++) if(track == edbv.GetTrack(i)) return 0;

  EdbVTA *vta = new EdbVTA(track, &edbv);
  vta->SetZpos(zpos);
  vta->SetFlag(2);
  edbv.AddVTA(vta);

  if(ntr<1)               // this is the first track in the vertex
    {
      edbv.SetXYZ( seg->X(), seg->Y(), seg->Z());
      return vta;
    }
  else
    {
      if( MakeV(edbv) ) 
	if( (edbv.V()->prob() > eProbMin) ) 	  return vta;
    }
  edbv.RemoveVTA(vta);
  if(ntr>1) 
    if( !MakeV(edbv) ) 	Log(1,"EdbVertexRec::AddTrack","vertex lost! ntr = %d",ntr);
  SafeDelete(vta);
  return 0;
}

//______________________________________________________________________________
Int_t EdbVertexRec::FindSimilarTracks(EdbTrackP &track, TObjArray &found, int nsegmin, float dMax, float dTheta, float dZmax)
{
  // find all tracks close to the "track" and return them in "found"
  // nsegmin - min number of segments for the interesting tracks
  // gMax    - max 3-d distance between track lines
  // dTheta  - max spatial angle between track lines
  // dZmax   - max distance in z between track lines

  using namespace TMath;
  int ntr  = eEdbTracks->GetEntries();
  if(ntr<1) return 0;

  EdbTrackP *t=0;
  EdbSegP   *s=track.TrackZmin(eUseSegPar);
  EdbSegP   *e=track.TrackZmax(eUseSegPar);
  EdbSegP   *ts=0, *te=0;  // start and end for the other tracks

  track.FitTrack();
  EdbSegP *t1=(EdbSegP*)(&track);
  EdbSegP *t2;
  int   nfound=0;
  float pv[3], imp, dtheta;
  for(int itr=0; itr<ntr; itr++)   {
    t = (EdbTrackP*)(eEdbTracks->At(itr));
    if (t == &track)                                          continue;
    if (t->Flag() < 0)                                        continue;
    ts = t->TrackZmin(eUseSegPar);
    te = t->TrackZmax(eUseSegPar);
    if( Min(ts->X(),te->X()) - Max(s->X(),e->X())   > dMax )  continue;
    if( Min(s->X(),e->X())   - Max(ts->X(),te->X()) > dMax )  continue;
    if( Min(ts->Y(),te->Y()) - Max(s->Y(),e->Y())   > dMax )  continue;
    if( Min(s->Y(),e->Y())   - Max(ts->Y(),te->Y()) > dMax )  continue;
    if( Min(ts->Z(),te->Z()) - Max(s->Z(),e->Z())   > dZmax ) continue;
    if( Min(s->Z(),e->Z())   - Max(ts->Z(),te->Z()) > dZmax ) continue;
    
    t->FitTrack();
    t2 = (EdbSegP*)t;
    dtheta = Sqrt( (t2->TX()-t1->TX())*(t2->TX()-t1->TX()) + (t2->TY()-t1->TY())*(t2->TY()-t1->TY()) );
    if(dtheta>dTheta)                                         continue;
    imp = CheckImpact( t1,t2,1,1, pv);
    if(imp>dMax)                                              continue;

    found.Add(t);
    nfound++;
  }

  Log(2,"EdbVertexRec::FindSimilarTracks","%d tracks found",nfound);
  if(gEDBDEBUGLEVEL>1) {
    for(int i=0; i<found.GetEntries(); i++) {
      t = (EdbTrackP*)found.At(i);
      t2 = (EdbSegP*)t;
      dtheta = Sqrt( (t2->TX()-t1->TX())*(t2->TX()-t1->TX()) + (t2->TY()-t1->TY())*(t2->TY()-t1->TY()) );
      imp = CheckImpact( t1,t2,1,1, pv);
      Log(3,"EdbVertexRec::FindSimilarTracks",
          "id =%6d  imp = %7.3f  dtheta = %7.3f  nseg =%3d\n", t->ID(), imp, dtheta,t->N());
   }
  }

  return nfound;
}


//______________________________________________________________________________
bool EdbVertexRec::CompatibleSegments( EdbSegP &pred, EdbSegP &stest, 
				       float impact, float dthetaMax, float dxy, 
				       float zminT, float zmaxT, float zminV, float zmaxV)
{
  float impMin   = 5;      // limits for the tracks separability
  float thetaMin = 0.005; //

  if( TMath::Abs(stest.X()-pred.X())   > dxy )     return 0;
  if( TMath::Abs(stest.Y()-pred.Y())   > dxy )     return 0;
  if( stest.Z()                        < zminT )   return 0;
  if( stest.Z()                        > zmaxT )   return 0;

  float dtheta = TMath::Sqrt( (stest.TX()-pred.TX())*(stest.TX()-pred.TX()) + (stest.TY()-pred.TY())*(stest.TY()-pred.TY()) );
  if(dtheta>dthetaMax)                             return 0;
  float pv[3];
  bool parallel;
  float imp = CheckImpactN( &pred,&stest,pv, parallel, eDZmax);
  if(imp>impact)                                   return 0;
  
  //printf("v: %f %f %f\n", pv[0],pv[1],pv[2]);
  //pred.PrintNice();
  
  if( dtheta<thetaMin && imp > 3*impMin)     return 0;   // parallel not close  tracks ( can be interesting for decay search??)
  
  if( dtheta>thetaMin)  {  // check z of the vertex
    float marg = impMin/dtheta;
    if( pv[2]    <  zminV-marg )            return 0;
    if( pv[2]    >  zmaxV+marg )            return 0;
  }

  Log(2,"EdbVertexRec::CompatibleSegments","ids = %6d and %6d    PH: %3d  imp = %7.3f  dtheta = %7.3f dzT = %7.3f dzV = %7.3f", 
      pred.ID(), stest.ID(), int(stest.W()), imp, dtheta, stest.Z()-pred.Z(), pv[2]-pred.Z());
  return 1;
}

//______________________________________________________________________________
Int_t EdbVertexRec::FindSimilarSegments( EdbSegP &spred, TObjArray &found, EdbPattern &pat,
					float impact, float dthetaMax, float dxy, 
					float zminT, float zmaxT, float zminV, float zmaxV)
{
  // Find all segments from path compatible with the pred segment and add them to found array

  int  nseg  = pat.N();  if(nseg<1) return 0;
  int  nfound=0;
  for(int i=0; i<nseg; i++)   {
    EdbSegP *s = pat.GetSegment(i);
    if (s->Flag() < 0)                                        continue;  // sure?
    if( CompatibleSegments( spred,*s, impact, dthetaMax, dxy, zminT, zmaxT, zminV, zmaxV) )  found.Add(s);
    else continue;
    nfound++;
  }
  Log(2,"EdbVertexRec::FindSimilarSegments","%d segments are found",nfound);
  return nfound;
}

//______________________________________________________________________________
Int_t EdbVertexRec::FindSimilarTracksE( EdbSegP &spred, TObjArray &found, bool startend,
					float impact, float dthetaMax, float dxy, 
					float zminT, float zmaxT, float zminV, float zmaxV)
{
  // Find all tracks with the requested extremity close to the pred segment and add them to found array
  //
  //   startend     - selecting tracks extremity to be used: 0-tracks starts (zmin)  1-tracks ends (zmax)
  //   impact       - max 3D distance between segment lines
  //   dthetaMax    - max spatial angle between track lines
  //   dxy          - preliminary distance cut between segments
  //   zminT,zmaxT  - limits in Z for the tracks extremity
  //   zminV,zmaxV  - limits in Z for the estimated vertex position

  int ntr  = eEdbTracks->GetEntries();
  if(ntr<1) return 0;
  int   nfound=0;

  EdbSegP *t1=&spred;
  for(int itr=0; itr<ntr; itr++)   {
    EdbTrackP *t = (EdbTrackP*)(eEdbTracks->At(itr));
    if (t->Flag() < 0)                                        continue;
    EdbSegP *te = t->TrackExtremity(startend,eUseSegPar);  // select track extremity

    if( CompatibleSegments( *t1,*te, impact, dthetaMax, dxy, zminT, zmaxT, zminV, zmaxV) )  found.Add(t);
    else continue;

    nfound++;
  }
  Log(2,"EdbVertexRec::FindSimilarTracksE","%d tracks found",nfound);
  return nfound;
}

//______________________________________________________________________________
int EdbVertexRec::FindVertex()
{
  // Note: in this function is assumed that all tracks selections are already done
  // ProbMin - minimal probability for chi2-distance between tracks

  //if(!ePVR) ePVR = ((EdbPVRec *)(gROOT->GetListOfSpecials()->FindObject("EdbPVRec")));
  if (ePVR) if (ePVR->IsA() != EdbPVRec::Class()) ePVR = 0;
  if(!ePVR) {Log(1,"EdbVertexRec::FindVertex","Error! EdbPVRec not defined, use SetPVRec(...)"); return 0;}

  EdbVertex *edbv = 0;
  TIndexCell starts,ends;              // "ist:entry"   "iend:entry"
  FillTracksStartEnd( starts, ends );

  if(gEDBDEBUGLEVEL>1) printf("-----Search 2-track vertexes----------------------------\n");

  int nvtx   = 0;

  if(gEDBDEBUGLEVEL>1) printf(" End-Begin tracks combinations:\n");
  nvtx += LoopVertex(ends  , starts,  0, 1 );

  if(gEDBDEBUGLEVEL>1) printf(" Begin-Begin tracks combinations:\n");
  nvtx += LoopVertex(starts, starts,  1, 1 );

  if(gEDBDEBUGLEVEL>1) printf(" End-End tracks combinations:\n");
  nvtx += LoopVertex(ends  , ends,    0, 0 );

  int nvtxt = 0;
  if (eVTX) nvtxt = eVTX->GetEntries();

  if(nvtx!=nvtxt) printf("ERROR: EdbVertexRec::FindVertex():  nvtx =%d nvtxt =%d\n",nvtx,nvtxt);

  //for (Int_t i = 0; i < nvtxt; i++) GetVertex(i)->SetID(i);

  if (nvtxt) eVTX->Sort(nvtxt-1);

  for (Int_t i = nvtx-1; i >= 0; i--) {
    edbv = GetVertex(i);
    if (!edbv) continue;
    edbv->SetID(i);
    edbv->ResetTracks();
  }

  if(gEDBDEBUGLEVEL>1) printf("----------------- %d vtx ---------------------------------------\n", nvtx);

  return nvtx;
}

//______________________________________________________________________________
void EdbVertexRec::FillTracksStartEnd(TIndexCell &starts, TIndexCell &ends )
{
  // fill tracks starts and ends lookup tables "z:entry"
  // inside sorted tracks: starts - minimal Z; ends - maximal Z

  EdbTrackP *tr=0;
  int ntr  = eEdbTracks->GetEntries();
  Long_t  v[2];

  EdbSegP *s=0;
  for(int itr=0; itr<ntr; itr++)   {
    tr = (EdbTrackP*)(eEdbTracks->At(itr));
    if (tr->Flag() < 0)                     continue;
    s = tr->TrackZmin(eUseSegPar);
    if(eUseLimits)  if( !IsInsideLimits(*s) ) continue;
    v[0] = (Long_t)(s->Z()/eZbin);
    v[1] = itr;
    starts.Add(2,v);
    s = tr->TrackZmax(eUseSegPar);
    if(eUseLimits)  if( !IsInsideLimits(*s) ) continue;
    v[0] = (Long_t)(s->Z()/eZbin);
    v[1] = itr;
    ends.Add(2,v);
  }
  starts.Sort();
  ends.Sort();
}

//______________________________________________________________________________
int EdbVertexRec::LoopVertex( TIndexCell &list1, TIndexCell &list2, 
			      int zpos1, int zpos2 )
{

  // zpos1  - the direction flag for the first  track 1-start, 0-end
  // zpos2  - the direction flag for the second track
  // in cycles is assumed that members of list1 has z <= members of list2

  Log(3,"EdbVertexRec::LoopVertex","  Selection: dZmax=%.0f Abin=%.3f ProbMin=%f zBin=%.0f usemom=%d",
         eDZmax, eAbin, eProbMin, eZbin, eUseMom);

  int nvtx     = 0; 
  int ncombin  = 0;
  int ncount   = 0;

  TIndexCell *c1=0,  *c2=0;
  EdbTrackP  *tr1=0, *tr2=0;
  int         itr1,  itr2;

  int   nz1 = list1.GetEntries();
  int   nz2 = list2.GetEntries();
  float z1, z2;

  //int ntot = nz1*nz2;
  //printf("  2-track vertexes search in progress... %3d%%", 0);

  for(int iz1=0; iz1<nz1; iz1++)   {           // first z-group
    c1 = list1.At(iz1);
    z1 = c1->Value()*eZbin;
    int nc1=c1->GetEntries();

    for(int iz2=0; iz2<nz2; iz2++)   {         // second z-group
      c2 = list2.At(iz2);
      z2 = c2->Value()*eZbin;

      if( z2 < z1 )           continue;
      if( z2-z1 > eDZmax )    break;

      ncount++;
      //printf("\b\b\b\b%3d%%",(int)((double)ncount/double(ntot)*100.));
      fflush(stdout);
      
      int nc2=c2->GetEntries();
      for(int ic1=0; ic1<nc1; ic1++) {      // first z-group entries

	itr1 = c1->At(ic1)->Value();
	tr1  = (EdbTrackP*)((*eEdbTracks)[itr1]);
	if(!tr1)             continue;

	int ic2start=0;
	if(c1==c2) ic2start=ic1+1;

	for(int ic2=ic2start; ic2<nc2; ic2++) {    // second z-group entries
      	  ncombin++;

	  itr2 = c2->At(ic2)->Value();
	  if(itr2==itr1)       continue; 
	  tr2  = (EdbTrackP*)((*eEdbTracks)[itr2]);
	  if(!tr2)             continue;

	  EdbVertex *vtx = ProbVertex2( tr1, tr2, zpos1, zpos2 );
	  if(vtx) {
	    AddVertex(vtx);
	    nvtx++;
	  }

	}
      }
    }
  }

  //printf("\b\b\b\b%3d%%\n",100);

  Log(3,"EdbVertexRec::LoopVertex","  %6d pairs -> %d vertices accepted\n", ncombin, nvtx);

  return nvtx;
}

//______________________________________________________________________________
float EdbVertexRec::CheckImpact( EdbSegP *s1,   EdbSegP *s2,
				 int zpos1,     int zpos2,    float pv[3] )
{
  // Return: the distance between 2 lines defined by s1 and s2
  // Input: s1,s2,zpos1,zpos2, where
  //        zpos1, zpos2: 0 - end   of the track (to be propagated forward in z) 
  //        zpos1, zpos2: 1 - start of the track (to be propagated backward in z) 
  // Output: pv - "vertex position" - the nearest point to the both lines

  float p1[3], p2[3];  // first line
  float p3[3], p4[3];  // second line
  float pa[3], pb[3];  // shortest line
  double mua, mub;

  float dz;

  p1[0] = s1->X();
  p1[1] = s1->Y();
  p1[2] = s1->Z();
  dz = zpos1? -eDZmax : eDZmax;
  p2[0] = s1->X() + dz*s1->TX();
  p2[1] = s1->Y() + dz*s1->TY();
  p2[2] = s1->Z() + dz;
  
  p3[0] = s2->X();
  p3[1] = s2->Y();
  p3[2] = s2->Z();
  dz = zpos2? -eDZmax : eDZmax;
  p4[0] = s2->X() + dz*s2->TX();
  p4[1] = s2->Y() + dz*s2->TY();
  p4[2] = s2->Z() + dz;

  if( !EdbMath::LineLineIntersect( p1, p2, p3, p4, pa, pb, mua, mub ) )   return 10e+10;
  if (pv) for (int i = 0; i<3; i++) pv[i] = 0.5*(pa[i] + pb[i]);
  return EdbMath::Magnitude3( pa, pb );
}


//______________________________________________________________________________
float EdbVertexRec::CheckImpactN( EdbSegP *s1,   EdbSegP *s2,
			          float pv[3], bool &parallel, float dzMax )
{
  // Return: the distance between 2 lines defined by s1 and s2
  // Output: pv - "vertex position" - the nearest point to the both lines

  float imp = 1e+10;
  float p1[3], p2[3];  // first line
  float p3[3], p4[3];  // second line
  float pa[3], pb[3];  // shortest line
  double mua, mub;

  p1[0] = s1->X() - dzMax*s1->TX();
  p1[1] = s1->Y() - dzMax*s1->TY();
  p1[2] = s1->Z() - dzMax;
  p2[0] = s1->X() + dzMax*s1->TX();
  p2[1] = s1->Y() + dzMax*s1->TY();
  p2[2] = s1->Z() + dzMax;
  
  p3[0] = s2->X() - dzMax*s2->TX();
  p3[1] = s2->Y() - dzMax*s2->TY();
  p3[2] = s2->Z() - dzMax;
  p4[0] = s2->X() + dzMax*s2->TX();
  p4[1] = s2->Y() + dzMax*s2->TY();
  p4[2] = s2->Z() + dzMax;
  
  if( EdbMath::LineLineIntersect( p1, p2, p3, p4, pa, pb, mua, mub ) ) {
    parallel = false;
    for (int i = 0; i<3; i++) pv[i] = 0.5*(pa[i] + pb[i]);
    imp =  EdbMath::Magnitude3( pa, pb );
  } else  { // the lines are parallel: take as vertex the mean point, calc. the distance to line
    parallel = true;
    pv[0] = 0.5*(s1->X()+s2->X());
    pv[1] = 0.5*(s1->Y()+s2->Y());
    pv[2] = 0.5*(s1->Z()+s2->Z());
    bool inside=0;
    imp = 2. * EdbMath::DistancePointLine3(pv, p1,p2, &inside);
  }
  return imp;
}

//______________________________________________________________________________
EdbVertex *EdbVertexRec::ProbVertex2( EdbTrackP *tr1, EdbTrackP *tr2,
				      int zpos1,      int zpos2 )
{
  // Check if 2 tracks can form the vertex. If yes - return the pointer to the new EdbVertex object
  Log(3,"EdbVertexRec::ProbVertex2", "try tracks %d and %d", tr1->ID(), tr2->ID() );
  
  if(!tr1)     return 0;
  if(!tr2)     return 0;
  if(tr1==tr2) return 0;
  EdbSegP *s1 = tr1->TrackExtremity(zpos1,eUseSegPar);
  EdbSegP *s2 = tr2->TrackExtremity(zpos2,eUseSegPar);
//zpos2? tr2->TrackZmin(eUseSegPar):tr2->TrackZmax(eUseSegPar);
  if(!s1)      return 0;
  if(!s2)      return 0;

  // check the dZ position agreement
  float dz     = TMath::Abs(s2->Z()-s1->Z());
  if( dz > eDZmax )                                       return 0;

  // check the dX dY position agreement
  int   isign;
  float dtx,dty,deltaZ=0;
  if(zpos1!=zpos2) {              // start-end,   end-start
    deltaZ = (dz+eZbin);
    isign = -1;
  } else {                        // start-start, end-end
    deltaZ = (eDZmax-dz/2.);
    isign = +1;
  }
  dtx = TMath::Abs(s2->TX() - isign*s1->TX())+eAbin;
  if( TMath::Abs(s2->X()-s1->X()) > dtx*deltaZ )          return 0;
  dty = TMath::Abs(s2->TY() - isign*s1->TY())+eAbin;
  if( TMath::Abs(s2->Y()-s1->Y()) > dty*deltaZ )          return 0;

  // check the impact
  float vestim[3];
  //float imp = CheckImpact( s1, s2, zpos1, zpos2, vestim );
  bool parallel;
  float imp = CheckImpactN( s1, s2, vestim, parallel, eDZmax );
  Log(3,"EdbVertexRec::ProbVertex2", "impact = %f", imp );
  if(imp>eImpMax)                                         return 0;

  // create the new EdbVertex
  EdbVertex *vtx =  new EdbVertex();
  vtx->SetXYZ(vestim[0],vestim[1],vestim[2]);

  EdbVTA    *vta1 =  new EdbVTA(tr1,vtx);
  vta1->SetZpos(zpos1);
  vta1->SetFlag(2);
  vtx->AddVTA(vta1);

  EdbVTA    *vta2 =  new EdbVTA(tr2,vtx);
  vta2->SetZpos(zpos2);
  vta2->SetFlag(2);
  vtx->AddVTA(vta2);
  
  if(MakeV(*vtx)) {
    vtx->SetFlag( EstimateVertexFlag(zpos1,zpos2) );
    //if(imp<eImpMaxV)     // accept the vertex even if prob is small (parallel tracks)
    //  return vtx;
    //else                 // do additional checks
    if( vtx->V()->prob() >= eProbMin ) {
       if( Log(3,"EdbVertexRec::ProbVertex2", "prob = %f",  vtx->V()->prob()) )	      vtx->Print();
       if( CheckDZ2( s1->Z(), s2->Z(), zpos1,zpos2, vtx->VZ() ) )
	     return vtx;
    }
  }
  SafeDelete(vta1);
  SafeDelete(vta2);
  SafeDelete(vtx); 
  return 0;
}

//______________________________________________________________________________
int EdbVertexRec::EstimateVertexFlag(int zpos1, int zpos2)
{
    if      (zpos1 == 0 && zpos2 == 1)   return 1;         // end   & start
    else if (zpos1 == 1 && zpos2 == 0)   return 1;         // end   & start
    else if (zpos1 == 1 && zpos2 == 1)   return 0;         // start & start
    else if (zpos1 == 0 && zpos2 == 0)   return 2;         // end   & end
    return -1;
}

//______________________________________________________________________________
Bool_t EdbVertexRec::IsInsideLimits(EdbSegP &s)
{
  // return 1 if the segment position (x,y,z) is inside the limits defined by eVmin,eVmax
  if(s.X()<eVmin.X()) return 0;
  if(s.Y()<eVmin.Y()) return 0;
  if(s.Z()<eVmin.Z()) return 0;
  if(s.X()>eVmax.X()) return 0;
  if(s.Y()>eVmax.Y()) return 0;
  if(s.Z()>eVmax.Z()) return 0;
  return 1;
}

//______________________________________________________________________________
Bool_t EdbVertexRec::CheckDZ2(float z1, float z2, int zpos1, int zpos2, float z )
{
  // return 1 if the vertex position (z) is in agreement with limits defined by  eZbin and eDZmax

  float zvmin,zvmax;
  if (zpos1 == 0 && zpos2 == 1)            // ends & starts
    {
      zvmin = TMath::Min(z1,z2) - 0.5*eZbin;
      zvmax = TMath::Max(z1,z2) + 0.5*eZbin;
      if (z < zvmin || z > zvmax)  return 0;
    }
  else if (zpos1 == 1 && zpos2 == 1)        // starts & starts
    {
      zvmax = TMath::Min(z1,z2) + eZbin;
      zvmin = TMath::Max(z1,z2) - eDZmax;
      if (z > zvmax || z < zvmin)	  return 0;
    }
  else if (zpos1 == 0 && zpos2 == 0)        // ends & ends
    {
      zvmin = TMath::Max(z1 ,z2);
      zvmax = TMath::Min(z1 ,z2) + eDZmax + eZbin;
      if (z < zvmin || z > zvmax )	  return 0;
    }
  return 1;
}

//______________________________________________________________________________
int EdbVertexRec::ProbVertexN()
{
  if(gEDBDEBUGLEVEL>1) printf("*** on entry to ProbVertexN: nv = %d\n", eVTX->GetEntries());
  ProbVertexNpos(0);  // find n-track vtx attached to left edges
  if(gEDBDEBUGLEVEL>1) printf("*** npoz     0  ProbVertexN: nv = %d\n", eVTX->GetEntries());
  ProbVertexNpos(1);  // to right edges
  if(gEDBDEBUGLEVEL>1) printf("*** npoz     1  ProbVertexN: nv = %d\n", eVTX->GetEntries());

  CheckVTX();         // rank the vertices and reassign tracks according to the major rank
  if(gEDBDEBUGLEVEL>0) StatVertexN();      // print vertex statistics
  return 0;
}

//______________________________________________________________________________
int EdbVertexRec::ProbVertexNpos(int zpos)
{
  // cycle by all vertices, check if it is possible to join some of them by common track, do it

  if (!eVTX) return 0;
  int nv2 = eVTX->GetEntries();    // number of vertices
  if(nv2<2)  return 0;

  Log(2,"EdbVertexRec::ProbVertexNpos","%d vertices as input",nv2);

  // first group vertices with a common track

  TMap maptr;  //key is track, value is TObjArray of vertices
  TObjArray *arr=0; 
  EdbVertex *vtx=0;
  EdbTrackP *t=0;
  for( int i=0; i<nv2; i++) {
    vtx = GetVertex(i);
    if(!vtx) Log(1,"EdbVertexRec::ProbVertexNpos","vertex not found!");
    for( int j=0; j<vtx->N(); j++) {
      t = vtx->GetTrack(j);
      if(!t) Log(1,"EdbVertexRec::ProbVertexNpos","track not found!");
      if(!(vtx->Zpos(j)==zpos)) continue;
      arr = (TObjArray*)maptr.GetValue(t);
      if(!arr) {
	arr = new TObjArray();
	maptr.Add(t,arr);
      }
      arr->Add(vtx);
    }
  }

  TIter next(maptr.GetTable());
  TPair *a;
  while ((a = (TPair*) next())) {   // cycle by all keys (tracks)
    t  = (EdbTrackP*)a->Key();
    arr = (TObjArray*)a->Value();
    int nv = arr->GetEntries();
    //printf( "nv = %d\n", nv );
    if(nv<2)    continue;
    TObjArray arrvta;        // group of vta of tracks attached to t
    for(int iv=0; iv<nv; iv++) {
      vtx = (EdbVertex*)arr->At(iv);
      for(int it=0; it<vtx->N(); it++) {	arrvta.Add(vtx->GetVTa(it));  }
    }
    EdbVertex *newvtx=TestVTAGroup(arrvta);
    if(newvtx) {
      Log(3,"EdbVertexRec::ProbVertexNpos","add new vtx with %d tracks and flag = %d\n",newvtx->N(),newvtx->Flag() );
      eVTX->Add(newvtx);
    }
  }

  Log(2,"EdbVertexRec::ProbVertexNpos","%d entries in the map for zpos = %d",maptr.GetEntries(), zpos);

  return 0;
}

//______________________________________________________________________________
void  EdbVertexRec::CheckVTX()
{
  // rank the vertices and reassign tracks according to the major vertex weight

  if(!eVTX)   return;
  int nvtx = eVTX->GetEntries();
  Log(2,"EdbVertexRec::CheckVTX","%d vertices",nvtx);
  if(nvtx<1)  return;

  TArrayF   weight(nvtx);  // the vertex "weight" = 10*ntr+prob
  TArrayI   ind(nvtx);
  EdbVertex *vtx=0;
  for(int i=0; i<nvtx; i++)  {
    vtx = GetVertex(i);
    //vtx->SetFlag(0);       // to check!!!
    weight[i] = 10*vtx->N() + vtx->V()->prob();
  }
  TMath::Sort(nvtx,weight.GetArray(),ind.GetArray(),0); // sort in ascending order
  for(int i=0; i<nvtx; i++)  {
    vtx = GetVertex(ind[i]);
    vtx->ResetTracks();
  }

  // discard vertices with the detached tracks
  for(int i=0; i<nvtx; i++)  {
    vtx = GetVertex(ind[i]);
    int ndisc = vtx->CheckDiscardedTracks();
    if(ndisc>0) {
      Log(2,"EdbVertexRec::CheckVTX","discard vtx i=%d ntr=%d flag before: %d  disc tracks: %d",i, vtx->N(), vtx->Flag(), ndisc);
      vtx->SetFlag(-10);
    } else {
      vtx->EstimateVertexFlag();
    }
  }

  // reassign the vertex id's
  for (Int_t i = 0; i < nvtx; i++) GetVertex(i)->SetID(i);

  //eVTA.Clear();  //TODO?
}

//______________________________________________________________________________
EdbVertex *EdbVertexRec::TestVTAGroup(TObjArray &arrvta)
{
  // Try to create  N-prong vertex from vta's group
  // Input: array of preselected vta's
  // return the new vertex if successful

  int nvta=arrvta.GetEntries();
  Log(3,"EdbVertexRec::TestVTAGroup","nvta = %d\n",nvta);

  EdbVTA *vta=0, *newvta=0;
  EdbVertex *newvtx = new EdbVertex();
  for(int i=0; i<nvta; i++) {
    vta = (EdbVTA*)arrvta.At(i);
    if( newvtx->TrackInVertex(vta->GetTrack()) )  continue;
    newvta =  new EdbVTA(vta->GetTrack(),newvtx);
    newvta->SetZpos(vta->Zpos());
    newvta->SetFlag(2);
    newvtx->AddVTA(newvta);
  }

  if(MakeV(*newvtx))
    if( newvtx->V() )
      if( newvtx->V()->valid() ) {
	//printf("prob = %f\n", newvtx->V()->prob() );
	if( newvtx->V()->prob() >= eProbMin )   // accept new N-tracks vertex
	  {
	    for(int i=0; i<nvta; i++)	{
	      if( ((EdbVTA*)arrvta.At(i))->GetVertex() != newvtx )
					  ((EdbVTA*)arrvta.At(i))->GetVertex()->SetFlag(-10);
	    }
	    //printf("accepted!\n");
	    return newvtx;
	  }
      }
  SafeDelete(newvtx);   // discard new vertex
  return 0;
}

//______________________________________________________________________________
int EdbVertexRec::ProbVertexN_old()
{
  // deprecated function (VT: 7/05/2008. Keeped for back-compatibility tests. 
  // After the complete testing of the new ProbVertexN this function can be removed

  EdbVTA *vta = NULL, *vta1 = NULL, *vta2 = NULL;
  EdbVertex *edbv1 = NULL;
  EdbVertex *edbv2 = NULL;
  Vertex *v = 0;
  EdbTrackP *tr = 0;
  EdbTrackP *tr2 = 0;
  Int_t zpos = 0;
  int nvtx = 0;
  int nadd = 0;
  int ncombin = 0;
  int ncombinv = 0;
  bool wasadded = false;
  float dz = 0.;
  
  if (eVTX) {
    nvtx = eVTX->GetEntries();
    for (Int_t i = 0; i < nvtx; i++) {
      edbv1 = GetVertex(i);
      if (edbv1) {
	if (edbv1->N() > 2) {
	  for (Int_t j = 0; j<edbv1->N(); j++)  eVTA.Remove(edbv1->GetVTa(j));
	  for (Int_t j = 0; j<edbv1->Nn(); j++) eVTA.Remove(edbv1->GetVTn(j));
	  tr  = edbv1->GetTrack(0);
	  tr2 = edbv1->GetTrack(1);
	  edbv1->Clear();
	  vta1 = AddTrack(*edbv1, tr, edbv1->Zpos(0));
	  vta2 = AddTrack(*edbv1, tr2, edbv1->Zpos(1));
	  MakeV(*edbv1);
	  v = edbv1->V();
	  v->findVertexVt();

	  if (!eQualityMode)
	    edbv1->SetQuality(v->prob()/(v->vtx_cov_x()+v->vtx_cov_y()));
	  else if (eQualityMode == 1) {
	    Double_t rms = v->rmsDistAngle();
	    if (rms) edbv1->SetQuality((Float_t)(1./rms));
	    else edbv1->SetQuality(10.e+35);
	  }
	  else edbv1->SetQuality(1.);

	  tr->AddVTA(vta1);
	  tr2->AddVTA(vta2);
	  AddVTA(vta1);
	  AddVTA(vta2);
	}
	else {
	  if (edbv1->Flag() < 0) {
	    zpos = edbv1->Zpos(0) + edbv1->Zpos(1);
	    if (!zpos) edbv1->SetFlag(2);
	    else if (zpos == 1) edbv1->SetFlag(1);
	    else if (zpos == 2) edbv1->SetFlag(0);
	  }
	}
      }
    }
    edbv1 = 0;
  }
  else return 0;

  zpos = 0;

  nvtx = eVTX->GetEntries();
  printf("-----Merge 2-track vertex pairs to N-track vertexes-----\n");
  printf("N-track vertexes search in progress... %3d%%", 0);

  int nprint = (int)(0.05*(double)nvtx);
  if (nprint <= 0) nprint = 1;

  for (Int_t i1 = 0; i1 < nvtx; i1++) {
    wasadded = false;
    edbv1 = GetVertex(i1);
    if (!(i1%nprint)) {
      printf("\b\b\b\b%3d%%",(int)((double)i1/double(nvtx)*100.));
      fflush(stdout);
    }
    if (!edbv1) continue;
    if (edbv1->Flag() == -10) continue;
    Int_t nt1 = edbv1->N();
    bool exist = false;
    if (nt1 == 2) {
      for (Int_t ic1 = 0; ic1 < nt1; ic1++) {
	tr = edbv1->GetTrack(ic1);
	zpos = edbv1->Zpos(ic1);
	if (zpos) {
	  if (tr->VertexS()) {
	    if (nt1 < tr->VertexS()->N()) {
	      exist = true;
	      break;
	    }
	  }
	}
	else {
	  if (tr->VertexE()) {
	    if (nt1 <  tr->VertexE()->N()) {
	      exist = true;
	      break;
	    }
	  }
	}
      }

      if (exist) {
	edbv1->SetFlag(-10);
	continue;
      }
    }
    for (Int_t i2 = i1+1; i2<nvtx; i2++) {
      edbv2 = GetVertex(i2);
      if (!edbv2) continue;
      if (edbv2->Flag() == -10) continue;
      if (edbv2->N() == 2) {
	// printf(" v1 id %d, v2 id %d\n", edbv1->ID(), edbv2->ID()); 
	nt1 = edbv1->N();
	int nt2 = edbv2->N();
	int it1=0;
	int nomatch = 1;
	while (it1 < nt1 && nomatch) {
	  int it2=0;
	  tr = edbv1->GetTrack(it1);
	  while ( (it2<nt2) && nomatch) {
	    if (edbv2->GetTrack(it2) == tr && 
		edbv1->Zpos(it1) == edbv2->Zpos(it2)) {
	      ncombin++;
	      if (!it2) {
		tr2 = edbv2->GetTrack(1);
		zpos = edbv2->Zpos(1);
	      }
	      else if (it2 == 1) {
		tr2 = edbv2->GetTrack(0);
		zpos = edbv2->Zpos(0);
	      }

	      exist = false;
	      for (int ic1=0; ic1<edbv1->N(); ic1++)
		if (tr2 == edbv1->GetTrack(ic1)) exist = true;

	      if (zpos) {
		if (tr2->VertexS()) {
		  if (tr2->VertexS()->N() > edbv1->N()) exist = true;
		}
	      }
	      else {
		if (tr2->VertexE()) {
		  if (tr2->VertexE()->N() > edbv1->N()) {
		    exist = true;
		  }
		}
	      }
	      if (!exist) {
		ncombinv++;
		if (zpos) dz = edbv1->VZ() - tr2->TrackZmin(eUseSegPar)->Z();
		else      dz = tr2->TrackZmax(eUseSegPar)->Z() - edbv1->VZ();
		if(dz <= eZbin)
		  if ((vta = AddTrack(*edbv1, tr2, zpos))) {
		    nomatch = 0;
		    wasadded = true;
		    edbv2->SetFlag(-10);
		    tr2->AddVTA(vta);
		    AddVTA(vta);
		    int vfl=edbv1->Flag();
		    if      (vfl==0&&zpos==0) edbv1->SetFlag(1);
		    else if (vfl==2&&zpos==1) edbv1->SetFlag(1);
		    // printf("Add track ID %d from vertex %d to vertex %d\n",
		    //				    tr2->ID(), i2, i1);
		  }
	      }
	      else {
		nomatch = 0;
	      }
	      edbv2->SetFlag(-10);
	    } // if one of tracks vertex 2 equal any track in vertex 1
	    it2++;
	  } // tracks in vertex 2
	  it1++;
	} // tracks in vertex 1
      } // if vertex 2 has rank 2
    } // second vertex loop
    if (wasadded) nadd++;
  }  // first vertex loop

  printf("\b\b\b\b%3d%%\n",100);

  printf("  %6d 2-track vertex pairs with common track\n", ncombin);
  printf("  %6d pairs when common track not yet attached\n  %6d N-track vertexes with Prob > %f\n",
	 ncombinv, nadd, eProbMin);
  printf("--------------------------------------------------------\n");

  for (int i1=0; (i1<nvtx); i1++) {
    edbv1 = GetVertex(i1);
    if (!edbv1) continue;
    if (edbv1->Flag() == -10) continue;
    edbv1->ResetTracks();
  }

  StatVertexN();
  return nadd;
}


//---------------------------------------------------------
void EdbVertexRec::StatVertexN()
{
  Int_t nvt = Nvtx();
  if (!nvt) return;
  TArrayI navtx(10);
  EdbVertex *v=0;
  Int_t ntv = 0;
  for (Int_t i = 0; i < nvt; i++) {
    v = GetVertex(i);
    if (!v || v->Flag() < 0) continue;
    ntv = v->N();
    if (ntv > 11) ntv = 11;
    navtx[ntv-2]++;
  }
  for (ntv = 0; ntv < 10; ntv++) {
    if (ntv < 9)
      printf("%5d vertexes with number of tracks  = %2d was found\n",
	     navtx[ntv], ntv+2);
    else
      printf("%5d vertexes with number of tracks >= %2d was found\n",
	     navtx[ntv], ntv+2);
  }
}

//---------------------------------------------------------
int EdbVertexRec::LinkedVertexes()
{
  // calculate the number of linked vertices (with Nv()>0) and set flag+3 for them
  int nvt = Nvtx();
  if (!nvt) return 0;
  EdbVertex *v = 0;
  int nvl = 0;
  for (int iv=0; iv<nvt; iv++) {
    v = GetVertex(iv);
    if (v)
      if (v->Flag() != -10)
	if (v->Nv() != 0) {
	  if (v->Flag() < 3) v->SetFlag(v->Flag()+3);
	  nvl++;
	}
  }
  return nvl;
}

//---------------------------------------------------------
int EdbVertexRec::SelVertNeighbor( EdbVertex *v, int seltype, float RadMax, int Dpat, TObjArray *ao)
{
  EdbSegP ss; // the virtual "vertex" segment

  if (!v) return 0;

  float x = v->VX();
  float y = v->VY();
  float z = v->VZ();

  ss.SetX(x);
  ss.SetY(y);
  ss.SetZ(z);
  ss.SetTX(0.);
  ss.SetTY(0.);
  ss.SetErrors(RadMax*RadMax, RadMax*RadMax, 0., 0., 0., 0.);

  TObjArray arr(20);
  ePVR->FindComplimentsVol(ss,arr,1,1,Dpat);
  
  int nseg = arr.GetEntries();
  EdbTrackP *tr = 0;
  EdbSegP   *s  = 0;
  int trflg  = 0;
  int trind  = 0;
  int ntr = eEdbTracks->GetEntries();

  int nadd = 0;
  for (int i=0; i<nseg; i++)
  {
    s = (EdbSegP *)(arr.At(i));
    if (!s) continue;
    tr = 0;
    trind = s->Track();
    trflg = 0;
    if ( trind >= 0 && trind < ntr)
    {
	if ((tr = (EdbTrackP *)eEdbTracks->At(trind))) 
	{
	    trflg = tr->Flag();
	    if (trflg != -10 && seltype == 0 && tr->MCEvt() >= -999)
	    {
//		if (!(ao->FindObject(tr)))
//		{
		    if (ao)
		    {
			tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
			ao->Add(tr);
		    }
		    nadd++;
//		}
		continue;
	    }
	}
    } 
    if((trind < 0 || trflg == -10)&&(seltype == 1))
    {
	if (ao) ao->Add(s);
	nadd++;
    }
  }
 
  if (seltype == 0)
  {
    nseg = ao->GetEntries();
    for (int i=0; i<nseg; i++)
    {
	tr = (EdbTrackP *)(ao->At(i)); 
	if (tr && (tr->MCEvt() < -999)) tr->SetMC( -tr->MCEvt()-2000, tr->MCTrack());
    }
  }
 
  return nadd;
}

//---------------------------------------------------------
int EdbVertexRec::SelSegNeighbor( EdbSegP *sin, int seltype, float RadMax, int Dpat, TObjArray *ao)
{
  EdbSegP ss; // the virtual "vertex" segment

  if (!sin) return 0;

  ss.SetX(sin->X());
  ss.SetY(sin->Y());
  ss.SetZ(sin->Z());
  ss.SetTX(sin->TX());
  ss.SetTY(sin->TY());
  ss.SetErrors(RadMax*RadMax, RadMax*RadMax, 0., 0., 0., 0.);

  TObjArray arr(1000);
  
  ePVR->FindComplimentsVol(ss,arr,1,1,Dpat);

  int nseg = arr.GetEntries();
  EdbTrackP *tr = 0;
  EdbSegP   *s  = 0;
  int trflg  = 0;
  int trind  = 0;
  int ntr = 0;
  if (eEdbTracks) ntr = eEdbTracks->GetEntries();

  int nadd = 0;
  for (int i=0; i<nseg; i++)
  {
    s = (EdbSegP *)(arr.At(i));
    if (!s) continue;
    tr = 0;
    trind = s->Track();
    trflg = 0;
    if ( trind >= 0 && trind < ntr)
    {
	if ((tr = (EdbTrackP *)eEdbTracks->At(trind))) 
	{
	    trflg = tr->Flag();
	    if (trflg != -10 && seltype == 0 && tr->MCEvt() >= -999)
	    {
//		if (!(ao->FindObject(tr)))
//		{
		    if (ao)
		    {
			tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
			ao->Add(tr);
		    }
		    nadd++;
//		}
		continue;
	    }
	}
    } 
    if((trind < 0 || trflg == -10)&&(seltype == 1))
    {
      if (ao) ao->Add(s);
	nadd++;
    }
  }
 
  if (seltype == 0)
  {
    nseg = ao->GetEntries();
    for (int i=0; i<nseg; i++)
    {
	tr = (EdbTrackP *)(ao->At(i)); 
	if (tr && (tr->MCEvt() < -999)) tr->SetMC( -tr->MCEvt()-2000, tr->MCTrack());
    }
  }

  return nadd;
}
//______________________________________________________________________________
int EdbVertexRec::AddSegmentToVertex(EdbSegP *s, float ImpMax, float ProbMin, float Mom)
{
    EdbVTA *vta = 0;
    EdbVertex *ePrevious = 0;

    if (eWorking == 0)
    {
	eWorking = new EdbVertex();
	int ntr = eVertex->N();
	int i = 0, n = 0;
	for(i=0; i<ntr; i++)
	{
	    if ((vta = AddTrack(*(eWorking), eVertex->GetTrack(i), eVertex->Zpos(i))))
	    {
		eVertex->GetTrack(i)->AddVTA(vta);
		n++;
	    }
	}
	if (n < 2)
	{
	  SafeDelete(eWorking);
	  eVertex->ResetTracks();
	  printf("Can't create working copy of the vertex!\n");
	  fflush(stdout);
	  return 0;
	}

	if (!MakeV(*(eWorking)))
	{
	  SafeDelete(eWorking);
	  eVertex->ResetTracks();
	  printf("Can't create working copy of the vertex!\n");
	  fflush(stdout);
	  return 0;
	}
    }
    else
    {
	ePrevious = eWorking;
	eWorking = new EdbVertex();
	int ntr = ePrevious->N();
	int i = 0, n = 0;
	for(i=0; i<ntr; i++)
	{
	    if ((vta = AddTrack(*(eWorking),(ePrevious)->GetTrack(i), (ePrevious)->Zpos(i))))
	    {
		(ePrevious->GetTrack(i))->AddVTA(vta);
		n++;
	    }
	}
	if (n < 2)
	{
	  SafeDelete(eWorking);
	  if (ePrevious)
	    {
	      eWorking = ePrevious;
	      eWorking->ResetTracks();
	    }
	    else
	      {
		eVertex->ResetTracks();
	      }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return 0;
	}

	if (!MakeV(*(eWorking)))
	{
	  SafeDelete(eWorking);

	    if (ePrevious)
	    {
		eWorking = ePrevious;
		eWorking->ResetTracks();
	    }
	    else
	    {
		eVertex->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return 0;
	}
    }
    float mass = 0.139; //pion
    EdbTrackP *Tr = new EdbTrackP(s, mass);
    Tr->SetP(Mom);
    Tr->FitTrackKFS();
    float ImpMaxSave = eImpMax;
    eImpMax = ImpMax;
    float ProbMinSave = eProbMin;
    eProbMin = ProbMin;
    if ((vta = AddTrack(*(eWorking), Tr, 1)))
    {
	if ( Tr->Z() >= eWorking->VZ() ) vta->SetZpos(1);
	else vta->SetZpos(0);
	Tr->AddVTA(vta);
	EdbVertex *eW = eWorking;
	eW->SetID(eVertex->ID());
	eW->V()->rmsDistAngle();
	int trind = eEdbTracks->GetEntries();
	Tr->SetID(trind);
	eEdbTracks->Add(Tr);
	Tr->SetSegmentsTrack();
	eImpMax = ImpMaxSave;
	eProbMin = ProbMinSave;
    }
    else
    {
	printf("Track not added! May be Prob < ProbMin=%f. Change ProbMin with 'TrackParams' button!\n",ProbMin);
	fflush(stdout);
	delete Tr;

	SafeDelete(eWorking);
	if (ePrevious)
	{
	    eWorking = ePrevious;
	    eWorking->ResetTracks();
	}
	else
	{
	    eVertex->ResetTracks();
	}
	eImpMax = ImpMaxSave;
	eProbMin = ProbMinSave;
	return 0;
    }
    return 1;
}
//______________________________________________________________________________
int EdbVertexRec::VertexPolish(EdbVertex *v, int refill, float RadMax, int Dpat, float ImpMax, float ProbMin, float Mom)
{
    if (refill) if (!VertexNeighbor(v, RadMax, Dpat, ImpMax)) return 0;

    EdbVTA *vta = 0;
    EdbTrackP *tn = 0;
    EdbSegP *sn = 0;
    EdbVertex *w = 0;

    int nt = v->N();
    int nn = v->Nn();
    int naddtot = 0, nmod = 0;
    int news = 0;
    if (nn)
    {
	// first of all try to propagate existing tracks with new momentum
	double p = Mom;
	int nadd = 0;
	for(int i=0; i<nt; i++)
	{
	    nadd = 0;
	    tn = v->GetTrack(i);
	    for (int ip=0; ip<2; ip++)
	    {
	     p = Mom/(ip+1);
	     tn->SetErrorP(0.2*0.2*p*p);
	     tn->SetP(p);
	     if (v->Zpos(i)) nadd += ePVR->PropagateTrack( *tn, true,  0.01, 3, 0 );
	     else            nadd += ePVR->PropagateTrack( *tn, false, 0.01, 3, 0 );
	    }
	    if (nadd) nmod++;
	    naddtot += nadd;
	}

	// then attach single segments with small impact to the vertex
	eVertex = v;
	for(int i=0; i<nn; i++)
	{
	    vta = eVertex->GetVTn(i);
            if (vta->Flag()  == 1) // neighbour segment
            {
        	sn = (EdbSegP *)vta->GetTrack();
		news += AddSegmentToVertex(sn, ImpMax, ProbMin, Mom);
            }
	}
	w = eVertex;
	if (eWorking != 0) w = eWorking;
	nt = w->N();
	p = Mom;
	nadd = 0;
	// then try to propagate one-segment tracks with new momentum
	for(int i=0; i<nt; i++)
	{
	    nadd = 0;
	    tn = w->GetTrack(i);
	    if (tn->N() > 1) continue;
	    for (int ip=0; ip<2; ip++)
	    {
	     p = Mom/(ip+1);
	     tn->SetErrorP(0.2*0.2*p*p);
	     tn->SetP(p);
	     if (w->Zpos(i)) nadd += ePVR->PropagateTrack( *tn, true,  0.01, 3, 0 );
	     else            nadd += ePVR->PropagateTrack( *tn, false, 0.01, 3, 0 );
	    }
	    if (nadd) nmod++;
	    naddtot += nadd;
	}
    }
    printf("%d single segments are attached, %d tracks are propagated (total %d segments are added).\n",
           news, nmod, naddtot);
    fflush(stdout);
    return (news+nmod);
}
//______________________________________________________________________________
void EdbVertexRec::AcceptPolish()
{
    AcceptModifiedVTX(eVertex, eWorking);
}
//______________________________________________________________________________
void EdbVertexRec::RejectPolish()
{
    CancelModifiedVTX(eVertex, eWorking);
}
//______________________________________________________________________________
int EdbVertexRec::VertexTuning(int seltype)
{
  //if(!ePVR) ePVR = ((EdbPVRec *)(gROOT->GetListOfSpecials()->FindObject("EdbPVRec")));
  if (ePVR) if (ePVR->IsA() != EdbPVRec::Class()) ePVR = 0;
  if(!ePVR) {printf("Warning: EdbVertexRec::VertexNeighbor: EdbPVRec not defined, use SetPVRec(...)\n"); return 0;}

  int iv = 0, nntr = 0, nn = 0;
  int nvt = 0;
  int ntr = 0;
  if (eEdbTracks) ntr = eEdbTracks->GetEntries();
  if (!ntr) return 0;
  if (eVTX) nvt = eVTX->GetEntries();
  if (!nvt) return 0;

  EdbVertex *v1 = 0, *v2 = 0;
  EdbVertex *v1n = 0, *v2n = 0;
  EdbVTA *vta = 0;
  EdbTrackP *tr1 = 0, *tr2 = 0;
  int ntr1 = 0, ntr2 = 0, was = 0;
  double impa1[50] = {0.}, chia1[50] = {0.}, imp1max = 0., chi1max = 0.;
  double impa2[50] = {0.}, chia2[50] = {0.}, imp2max = 0., chi2max = 0.;
  double cri1max = 0., cri2max = 0.;
  int it1impmax = -1, it1chimax = -1, it1max = -1;
  int it2impmax = -1, it2chimax = -1, it2max = -1;
  double vchisumorig = 0., vdistsumorig = 0.;
  double v1chiorig = 0., v1distorig = 0.;
  double v2chiorig = 0., v2distorig = 0., crit = 0., critorig = 0.;

  for (iv=0; iv<nvt; iv++) {  // loop on all vertexes
    v1 = GetVertex(iv);
    if (v1)
    {
	    if (v1->Flag() < 0) continue;
	    if ((!(v1->V()))) continue;
	    ntr1 = v1->N();
	    imp1max = -1.;
	    chi1max = -1.;
	    it1max  = -1;
	    for(int it1=0; it1<v1->N() && it1<50; it1++)
	    {
		tr1 = v1->GetTrack(it1);
		impa1[it1] = v1->ImpTrack(it1);
		if (impa1[it1] > imp1max)
		{
		    imp1max = impa1[it1];
		    it1impmax = it1;
		}
		chia1[it1] = v1->Chi2Track(tr1, v1->Zpos(it1), 0.);
		if (chia1[it1] > chi1max)
		{
		    chi1max = chia1[it1];
		    it1chimax = it1;
		}
	    }
	    if (seltype == 0)
	    {
		it1max = it1impmax;
	    }
	    else
	    {
		it1max = it1chimax;
	    }
	    nntr = v1->Nn();
	    v1chiorig = v1->V()->chi2();
	    v1distorig = v1->V()->rmsDistAngle();
	    if (nntr)
	    {
	      for(int i=0; i<nntr; i++)
	      {
		vta = v1->GetVTn(i);
                if (vta->Flag() == 0)    // neighbour track
                {
//                       EdbTrackP *tn = vta->GetTrack();
                }
                else if (vta->Flag()  == 1) // neighbour segment
                {
//                       EdbSegP *sn = (EdbSegP *)vta->GetTrack();
                }
                else if (vta->Flag()  == 3) // neighbour vertex 
                {
                       v2 = (EdbVertex *)vta->GetTrack();
		       if ((!v2) || v2->Flag() == -10) continue;
		       if ((!(v2->V()))) continue;
		       ntr2 = v2->N();
		       v2chiorig = v2->V()->chi2();
		       v2distorig = v2->V()->rmsDistAngle();
		       vchisumorig = v1chiorig + v2chiorig;
		       vdistsumorig = v1distorig + v2distorig;
		       if (seltype == 0)
		       {
			    critorig = vdistsumorig;
		       }
		       else
		       {
			    critorig = vchisumorig;
		       }
		       was = 0;
		       if (ntr1==2 && ntr2>2)
		       {
			    imp2max = 0.;
			    chi2max = 0.;
			    it2max  = -1;
			    for(int it2=0; it2<ntr2 && it2<50; it2++)
			    {
				tr2 = v2->GetTrack(it2);
				impa2[it2] = v2->ImpTrack(it2);
				if (impa2[it2] > imp2max)
				{
				    imp2max = impa2[it2];
				    it2impmax = it2;
				}
				chia2[it2] = v2->Chi2Track(tr2, v2->Zpos(it2), 0.);
				if (chia2[it2] > chi2max)
				{
				    chi2max = chia2[it2];
				    it2chimax = it2;
				}
			    }
			    if (seltype == 0)
			    {
				it2max = it2impmax;
			    }
			    else
			    {
				it2max = it2chimax;
			    }
			    if (it2max < 0) continue;
			    if (v2->GetConnectedVertexForTrack(it2max)==v1) continue;
			    crit = MoveTrackToOtherVertex(v2, it2max, v1, seltype, &v2n, &v1n);
			    was = 1;
		       }
		       else if (ntr1>2 && ntr2==2)
		       {
			    imp1max = 0.;
			    chi1max = 0.;
			    it1max  = -1;
			    for(int it1=0; it1<ntr1 && it1<50; it1++)
			    {
				tr1 = v1->GetTrack(it1);
				impa1[it1] = v1->ImpTrack(it1);
				if (impa1[it1] > imp1max)
				{
				    imp1max = impa1[it1];
				    it1impmax = it1;
				}
				chia1[it1] = v1->Chi2Track(tr1, v1->Zpos(it1), 0.);
				if (chia1[it1] > chi1max)
				{
				    chi1max = chia1[it1];
				    it1chimax = it1;
				}
			    }
			    if (seltype == 0)
			    {
				it1max = it1impmax;
			    }
			    else
			    {
				it1max = it1chimax;
			    }
			    if (it1max < 0) continue;
			    if (v1->GetConnectedVertexForTrack(it1max)==v2) continue;
			    crit = MoveTrackToOtherVertex(v1, it1max, v2, seltype, &v2n, &v1n);
			    was = 1;
		       }
		       else if (ntr1>2 && ntr2>2)
		       {
			    imp2max = 0.;
			    chi2max = 0.;
			    it2max  = -1;
			    for(int it2=0; it2<ntr2 && it2<50; it2++)
			    {
				tr2 = v2->GetTrack(it2);
				impa2[it2] = v2->ImpTrack(it2);
				if (impa2[it2] > imp2max)
				{
				    imp2max = impa2[it2];
				    it2impmax = it2;
				}
				chia2[it2] = v2->Chi2Track(tr2, v2->Zpos(it2), 0.);
				if (chia2[it2] > chi2max)
				{
				    chi2max = chia2[it2];
				    it2chimax = it2;
				}
			    }
			    if (seltype == 0)
			    {
				it2max = it2impmax;
				cri1max = imp1max;
				cri2max = imp2max;
			    }
			    else
			    {
				it2max = it2chimax;
				cri1max = chi1max;
				cri2max = chi2max;
			    }
			    if (it1max < 0) continue;
			    if (it2max < 0) continue;
			    if (cri2max > cri1max)
			    {
				if (v2->GetConnectedVertexForTrack(it2max)==v1) continue;
				crit = MoveTrackToOtherVertex(v2, it2max, v1, seltype, &v2n, &v1n);
				was = 1;
			    }
			    else
			    {
				if (v1->GetConnectedVertexForTrack(it1max)==v2) continue;
				crit = MoveTrackToOtherVertex(v1, it1max, v2, seltype, &v2n, &v1n);
				was = 1;
			    } // variants of changing 
		       } // multiplicities
		       if (was)
		       {
		        if ((crit < critorig) && v1n && v2n)
		        {
			    AcceptModifiedVTX(v1,v1n);
			    AcceptModifiedVTX(v2,v2n);
			    v1 = v1n;
			    v2 = v2n;
		    	    if (v2->Flag() > -11) v2->SetFlag(-v2->Flag()-11);
			    nn++;
			    //break;
		        }
		        else
		        {
			    CancelModifiedVTX(v1,v1n);
			    CancelModifiedVTX(v2,v2n);
		        }
		       }
                } // neighbor vertex
	      } // loop on neighbor
	    } // neighbor exist
    } // good v1
  } // loop on vertexes

  nvt = eVTX->GetEntries();

  for (iv=0; iv<nvt; iv++) {
    v1 = GetVertex(iv);
    if (v1)
    {
	    if (v1->Flag()<-10)
	    {
		v1->SetFlag(-v1->Flag()-11);
	    }
    }
  }

  return nn;
} 
//______________________________________________________________________________
double EdbVertexRec::MoveTrackToOtherVertex(EdbVertex *v2, int it2max, EdbVertex *v1, int seltype,
                                            EdbVertex **v2no, EdbVertex **v1no)
{
//    printf("Rearrange vertexies %d and %d\n",v1->ID(),v2->ID());
    *v1no = 0;
    *v2no = 0;
    if (!v1 || !v2) return 0.;
    if (v2->N() < 3) return 0.;
    if (it2max >= v2->N()) return 0.;
    EdbVertex *v1n = 0, *v2n = 0;
    EdbTrackP *tr2 = 0;
    int zpos = 0;
    double dx, dy, dz, dist1, dist2, imp, vchisum, vdistsum;
    tr2 = v2->GetTrack(it2max);
    dx = v1->VX() - tr2->TrackZmin()->X();
    dy = v1->VY() - tr2->TrackZmin()->Y();
    dz = v1->VZ() - tr2->TrackZmin()->Z();
    dist1 = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
    dx = v1->VX() - tr2->TrackZmax()->X();
    dy = v1->VY() - tr2->TrackZmax()->Y();
    dz = v1->VZ() - tr2->TrackZmax()->Z();
    dist2 = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
    if (dist2 < dist1)
    {
	zpos = 0;
    }
    else
    {
	zpos = 1;
    }
    imp = v1->DistTrack(tr2, zpos);
    if (imp > 1.2*eImpMax) return 0.; 
    if ((v2n = RemoveTrackFromVertex(v2, it2max)))
    {
	*v2no = v2n;
	v1n = AddTrackToVertex(v1, tr2, zpos);
	if (!v1n)
	{
//	    CancelModifiedVTX(v2,v2n);
	    return 0.;
	}
	*v1no = v1n;
	vchisum = v1n->V()->chi2() + v2n->V()->chi2();
	vdistsum = v1n->V()->rmsDistAngle() + v2n->V()->rmsDistAngle();
	if (seltype == 0)
	{
	    return vdistsum;
	}
	else
	{
	    return vchisum;
	}
    }
    return 0.;
}

//______________________________________________________________________________
int EdbVertexRec::RefitAll()
{
  // use already found vertex position to improve the fit
  int nv = eVTX->GetEntries();
  EdbVertex *vtx = 0;
  int cnt=0;
  for (int i=0; i<nv; i++)   {
    vtx = GetVertex(i);
    if(vtx)
      if(vtx->Flag()>-1) {
	MakeV(*vtx,true);
	cnt++;
      }
  }
  return cnt;
}

//______________________________________________________________________________
EdbVertex *EdbVertexRec::AddTrackToVertex(EdbVertex *eVertex, EdbTrackP *eTr, int zpos)
{
    EdbVTA *vta = 0;
    EdbVertex *old = 0;
    EdbVertex *eWorking = 0;
    if (!eVertex)
    {
    
	    printf("No working vertex selected!\n");
	    fflush(stdout);
	    return 0;
    }
    if (!eTr)
    {
    
	    printf("No working track selected!\n");
	    fflush(stdout);
	    return 0;
    }
    if ((old = eTr->VertexS()) && (zpos == 1))
    {
    
	    printf("Track alredy connected to a vertex by this edge!\n");
	    fflush(stdout);
	    return 0;
    }
    if ((old = eTr->VertexE()) && (zpos == 0))
    {
    
	    printf("Track alredy connected to a vertex by this edge!\n");
	    fflush(stdout);
	    return 0;
    }
    if (eWorking == 0)
    {
	eWorking = new EdbVertex();
	int ntr = eVertex->N();
	int i = 0, n = 0;
	for(i=0; i<ntr; i++)
	{
	    if ((vta = AddTrack(*(eWorking), (eVertex)->GetTrack(i), (eVertex)->Zpos(i))))
	    {
		(eVertex)->GetTrack(i)->AddVTA(vta);
		n++;
	    }
	}
	if (n < 2)
	{
	  SafeDelete(eWorking);
	  (eVertex)->ResetTracks();
	  printf("Can't create working copy of the vertex!\n");
	  fflush(stdout);
	  return 0;
	}

	if (!MakeV(*eWorking))
	{
	  SafeDelete(eWorking);
	  eVertex->ResetTracks();
	  printf("Can't create working copy of the vertex!\n");
	  fflush(stdout);
	  return 0;
	}
    }
    if ((vta = AddTrack(*eWorking, eTr, zpos)))
    {
	eTr->AddVTA(vta);
	eWorking->SetID(eVertex->ID());
    }
    else
    {
//	printf("Track not added! May be Prob < ProbMin.\n");
//	fflush(stdout);
      SafeDelete(eWorking);
      eVertex->ResetTracks();
      return 0;
    }
    return eWorking;
}

//_____________________________________________________________________________
EdbVertex *EdbVertexRec::RemoveVTAFromVertex(EdbVertex &v, EdbVTA &vta)
{
  int ntr = v.N();
  TObjArray tracks;
  for(int i=0; i<ntr; i++) {
    if( v.GetVTa(i) != &vta ) tracks.Add( v.GetVTa(i)->GetTrack() );
  }
  EdbVertex *vnew = Make1Vertex( tracks, v.Z() );
  
  int nn = v.Nn();
  for(int i=0; i<nn; i++) {
    vnew->AddVTA( v.GetVTn(i) );
  }
  return vnew;
}

//_____________________________________________________________________________
EdbVertex *EdbVertexRec::RemoveTrackFromVertex(EdbVertex *eVertex, int itr)
{
    if (!eVertex)
    {
    
	    printf("No working vertex selected!\n");
	    fflush(stdout);
	    return 0;
    }
    EdbVTA *vta = 0;
    EdbVertex *eWorking = 0;     // is a bug?? redeclaration of a member; noted VT:25/08/2011
    int n = 0;
    int ntr = 0;
    if (eWorking == 0)
    {
	ntr = eVertex->N();
	if (ntr < 3)
	{
    
	    printf("Working vertex has 2 prongs only!\n");
	    fflush(stdout);
	    return 0;
	}
	eWorking = new EdbVertex();
	int i = 0;
	for(i=0; i<ntr; i++)
	{
	    if (i == itr)
	    {
		(eVertex->GetTrack(i))->ClearVTA(eVertex->GetVTa(i));
		continue;
	    }
	    if ((vta = AddTrack( *eWorking, eVertex->GetTrack(i), eVertex->Zpos(i))))
	    {
		(eWorking->GetTrack(n))->AddVTA(vta);
		n++;
	    }
	}
    }
    if ((n < 2)||(n == ntr))
    {
	printf("Can't create working copy of the vertex!\n");
	fflush(stdout);
	SafeDelete(eWorking);
	eVertex->ResetTracks();
	return 0;
    }

    if (MakeV(*eWorking))
    {
	EdbVertex *eW = eWorking;
	eW->ResetTracks();
	eW->SetID(eVertex->ID());
    }
    else
    {
	printf("Can't create working copy of the vertex!\n");
	fflush(stdout);
	SafeDelete(eWorking);
	eVertex->ResetTracks();
	return 0;
    }
    return eWorking;
}
//_____________________________________________________________________________
void EdbVertexRec::AcceptModifiedVTX(EdbVertex *eVertex, EdbVertex *eWorking)
{
    if (eWorking && eVertex)
    {
        EdbVertex *eW = eWorking;
	int ind = eVertex->ID();
	eW->SetID(ind);
	eW->SetQuality(eW->V()->prob()/
			   (eW->V()->vtx_cov_x()+eW->V()->vtx_cov_y()));
	int ntr = eVertex->N();
	for(int i=0; i<ntr; i++)
	{
	    eVTA.Remove(eVertex->GetVTa(i));
	}

	eW->ResetTracks();
	ntr = eW->N();
	int ifl = 0;
	for(int i=0; i<ntr; i++)
	{
		if (eW->Zpos(i)) ifl = ifl | 1; 
		else		 ifl = ifl | 2; 
	}
	ifl = 4 - ifl;
	if (ifl == 3) ifl = 0;
	if (eW->Nv()) ifl += 3;
	eW->SetFlag(ifl);
	if (eVTX) eVTX->AddAt(eW, ind);
	for(int i=0; i<ntr; i++)
	{
	    AddVTA(eW->GetVTa(i));
	}
	ntr = eVertex->Nn();
	for(int i=0; i<ntr; i++)
	{
	    eW->AddVTA(eVertex->GetVTn(i));
	}
//	if (eVertex) delete eVertex;
	eW->ResetTracks();
    }
}
//_____________________________________________________________________________
void EdbVertexRec::CancelModifiedVTX(EdbVertex *eVertex, EdbVertex *eWorking)
{
  SafeDelete(eWorking);
  if (eVertex) eVertex->ResetTracks();
}
//______________________________________________________________________________
int EdbVertexRec::VertexNeighbor(float RadMax, int Dpat, float ImpMax)
{
  //if(!ePVR) ePVR = ((EdbPVRec *)(gROOT->GetListOfSpecials()->FindObject("EdbPVRec")));
  if (ePVR) if (ePVR->IsA() != EdbPVRec::Class()) ePVR = 0;
  if(!ePVR) {printf("Warning: EdbVertexRec::VertexNeighbor: EdbPVRec not defined, use SetPVRec(...)\n"); return 0;}

  int nn   = 0, iv = 0;
//  int i = 0, nntr = 0;
  int nvt = 0;
  int ntr = 0;
  if (eVTX) nvt = eVTX->GetEntries();
  if (eEdbTracks) ntr = eEdbTracks->GetEntries();
  if (!nvt) return 0;
  if (!ntr) return 0;

  EdbVertex *v   = 0;

  for (iv=0; iv<nvt; iv++) {
    v = GetVertex(iv);
    if (v)
    {
	    nn += VertexNeighbor(v, RadMax, Dpat, ImpMax);
//	    nntr = v->Nn();
//	    for(i=0; i<nntr; i++) AddVTA(v->GetVTn(i));
    }
  }

  return nn;
} 

//______________________________________________________________________________
int EdbVertexRec::VertexNeighbor(EdbVertex *v, float RadMax, int Dpat, float ImpMax)
{
  //if(!ePVR) ePVR = ((EdbPVRec *)(gROOT->GetListOfSpecials()->FindObject("EdbPVRec")));
  if (ePVR) if (ePVR->IsA() != EdbPVRec::Class()) ePVR = 0;
  if(!ePVR) {printf("Warning: EdbVertexRec::VertexNeighbor: EdbPVRec not defined, use SetPVRec(...)\n"); return 0;}

  EdbVTA    *vta = 0;
  EdbTrackP *tr  = 0;
  const EdbSegP   *ss  = 0;
  const EdbSegP   *se  = 0;
  EdbTrackP *trv = 0;
  EdbVertex *ve = 0;
  float dx = 0., dy = 0.;
  float dz = 0;
  int 	    zpos = 0;
  int 	    nn   = 0, nntr = 0;
  float     distxs, distys, distzs1, distzs, dists, dist = 0.;
  float     distxe = 0., distye = 0., distze1= 0., distze = 0., diste = 0.;
  float     xvert = 0, yvert = 0, zvert = 0;
  Float_t   Zbin = TMath::Abs((ePVR->GetPattern(1))->Z() - (ePVR->GetPattern(0))->Z());
  TObjArray an(20);

  if (v->Flag() != -10)
	{
	    v->ClearNeighborhood();
	    xvert = v->VX();
	    yvert = v->VY();
	    zvert = v->VZ();
	    // Select tracks neigborhood
	    an.Clear();
	    int ntr = v->N();
	    for(int i=0; i<ntr; i++)
	    {
		tr = v->GetTrack(i);
		tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
		if (v->Zpos(i) == 1)
		    ve = tr->VertexE(); 
		else
		    ve = tr->VertexS(); 
		if (ve && ve->Flag() >= -99 && ve->Flag() != -10)
		{
		    dx = ve->VX() - xvert;
		    dy = ve->VY() - yvert;
		    if (TMath::Sqrt(dx*dx + dy*dy) > RadMax) continue;
		    dz = ve->VZ() - zvert;
		    if (TMath::Abs(dz) > Dpat*Zbin) continue;
		    vta = new EdbVTA((EdbTrackP *)ve, v);
		    vta->SetFlag(3);
		    vta->SetDist(TMath::Sqrt(dx*dx+dy*dy+dz*dz));
		    v->AddVTA(vta);
		    ve->SetFlag(-ve->Flag()-200);
		    for(int j=0; j<ve->N(); j++)
		    {
			trv = ve->GetTrack(j);
			if (trv->MCEvt() < -999) continue;
			if (ve->Zpos(j) == 0)
			    ss = trv->TrackZmax();
			else
			    ss = trv->TrackZmin();
			dx = ss->X() - xvert;
			dy = ss->Y() - yvert;
			dz = ss->Z() - zvert;
			if (TMath::Sqrt(dx*dx+dy*dy) > RadMax)    continue;
			if (TMath::Abs(dz)        > Dpat*Zbin) continue;
			dist = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
			vta = v->CheckImp(trv, 10.e+10, ve->Zpos(j), dist);
			if (vta) trv->SetMC(-trv->MCEvt()-2000, trv->MCTrack());
		    }
		}
	    }
	    int nvn = SelVertNeighbor(v, 0, RadMax, Dpat, &an); 
	    for (int it=0; it<nvn; it++) {
		    tr = (EdbTrackP*)(an.At(it));
		    if (tr)
		    {
			    if (tr->MCEvt() < -999) continue;
			    ss = tr->TrackZmin();
			    distxs = (xvert - ss->X());
			    distxs *= distxs;
			    distys = (yvert - ss->Y());
			    distys *= distys;
			    distzs1 = (zvert - ss->Z());
			    distzs = distzs1*distzs1;
			    dists  =  distxs + distys + distzs;
			    dists  =  TMath::Sqrt(dists); 
			    se = tr->TrackZmax();
			    distxe = (xvert - se->X());
			    distxe *= distxe;
			    distye = (yvert - se->Y());
			    distye *= distye;
			    distze1 = (zvert - se->Z());
			    distze = distze1*distze1;
			    diste  =  distxe + distye + distze;
			    diste  =  TMath::Sqrt(diste);
			    if (diste < dists)
			    {
				if (TMath::Sqrt(distxe+distye) > RadMax)    continue;
				if (TMath::Abs(distze1)        > Dpat*Zbin) continue;
				dist = diste;
				zpos = 0;
			    }
			    else
			    {
				if (TMath::Sqrt(distxs+distys) > RadMax)    continue;
				if (TMath::Abs(distzs1)        > Dpat*Zbin) continue;
				dist = dists;
				zpos = 1;
			    }
			    vta = v->CheckImp(tr, ImpMax, zpos, dist);
			    if (vta)
			    {
				nn++;
				tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
				ve = tr->VertexS(); 
				if (ve && ve->Flag() >= -99 && ve->Flag() != -10 && ve != v)
				{
				    dx = ve->VX() - xvert;
				    dy = ve->VY() - yvert;
				    if (TMath::Sqrt(dx*dx + dy*dy) <= RadMax)
				    {
					dz = TMath::Abs(ve->VZ() - zvert);
					if (dz <= Dpat*Zbin)
					{
					    ve->SetFlag(-ve->Flag()-200);
					    vta = new EdbVTA((EdbTrackP *)ve, v);
					    vta->SetFlag(3);
					    vta->SetDist(TMath::Sqrt(dx*dx+dy*dy+dz*dz));
					    v->AddVTA(vta);
					    for(int j=0; j<ve->N(); j++)
					    {
						trv = ve->GetTrack(j);
						if (trv->MCEvt() < -999) continue;
						if (ve->Zpos(j) == 0)
						    ss = trv->TrackZmax();
						else
						    ss = trv->TrackZmin();
						dx = ss->X() - xvert;
						dy = ss->Y() - yvert;
						dz = ss->Z() - zvert;
						if (TMath::Sqrt(dx*dx+dy*dy) > RadMax)    continue;
						if (TMath::Abs(dz)        > Dpat*Zbin) continue;
						dist = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
						vta = v->CheckImp(trv, 10.e+10, ve->Zpos(j), dist);
    						if (vta) trv->SetMC(-trv->MCEvt()-2000, trv->MCTrack());
					    }
					}
				    }
				}
				ve = tr->VertexE(); 
				if (ve && ve->Flag() >= -99 && ve->Flag() != -10 && ve != v)
				{
				    dx = ve->VX() - xvert;
				    dy = ve->VY() - yvert;
				    if (TMath::Sqrt(dx*dx + dy*dy) <= RadMax)
				    {
					dz = TMath::Abs(ve->VZ() - zvert);
					if (dz <= Dpat*Zbin)
					{
					    ve->SetFlag(-ve->Flag()-200);
					    vta = new EdbVTA((EdbTrackP *)ve, v);
					    vta->SetFlag(3);
					    vta->SetDist(TMath::Sqrt(dx*dx+dy*dy+dz*dz));
					    v->AddVTA(vta);
					    for(int j=0; j<ve->N(); j++)
					    {
						trv = ve->GetTrack(j);
						if (trv->MCEvt() < -999) continue;
						if (ve->Zpos(j) == 0)
						    ss = trv->TrackZmax();
						else
						    ss = trv->TrackZmin();
						dx = ss->X() - xvert;
						dy = ss->Y() - yvert;
						dz = ss->Z() - zvert;
						if (TMath::Sqrt(dx*dx+dy*dy) > RadMax)    continue;
						if (TMath::Abs(dz)        > Dpat*Zbin) continue;
						dist = TMath::Sqrt(dx*dx+dy*dy+dz*dz);
						vta = v->CheckImp(trv, 10.e+10, ve->Zpos(j), dist);
    						if (vta) trv->SetMC(-trv->MCEvt()-2000, trv->MCTrack());
					    }
					}
				    }
				}
			    }
		    }
	    }
	    ntr = v->N();
	    for(int i=0; i<ntr; i++)
	    {
		tr = v->GetTrack(i);
		tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
	    }
	    nntr = v->Nn();
	    for(int i=0; i<nntr; i++)
	    {
		if ((vta = v->GetVTn(i)))
		{
		    if (vta->Flag() == 0) //track
		    {
			tr = vta->GetTrack();
    			if (tr->MCEvt() < -999 ) tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
		    }
		    else if (vta->Flag() == 3) //vertex
		    {
    			ve = (EdbVertex *)vta->GetTrack();
    			if (ve->Flag() < -99 ) ve->SetFlag(-ve->Flag()-200);
		    }
		}
	    }
	    // Select segments neigborhood
	    an.Clear();
	    nvn = SelVertNeighbor(v, 1, RadMax, Dpat, &an); 
	    for (int it=0; it<nvn; it++) {
		    ss = (EdbSegP*)(an.At(it));
		    if (ss)
		    {
			    distxs = (xvert - ss->X());
			    distxs *= distxs;
			    distys = (yvert - ss->Y());
			    distys *= distys;
			    distzs1 = (zvert - ss->Z());
			    distzs = distzs1*distzs1;
			    dists  =  distxs + distys + distzs;
			    dists  =  TMath::Sqrt(dists); 
//			    if (TMath::Sqrt(distxs+distys) > RadMax)    continue;
//			    if (TMath::Abs(distzs1)        > Dpat*Zbin) continue;
//			    vta = v->CheckImp((EdbTrackP *)ss, ImpMax, zpos, dists);
			    if (v->DistSeg((EdbSegP *)ss) > ImpMax) continue;
			    vta = new EdbVTA((EdbTrackP *)ss, v);
			    vta->SetZpos(1);
			    vta->SetFlag(1);
			    vta->SetImp(0.);
			    vta->SetDist(dists);
			    v->AddVTA(vta);
			    nn++;
		    }
	    }
  }
  return nn;
} 

//______________________________________________________________________________
int EdbVertexRec::SegmentNeighbor(EdbSegP *s, float RadMax, int Dpat, float ImpMax, float SegWmin, TObjArray *arrs, TObjArray *arrt, TObjArray *arrv)
{
  //if(!ePVR) ePVR = ((EdbPVRec *)(gROOT->GetListOfSpecials()->FindObject("EdbPVRec")));
  if (ePVR) if (ePVR->IsA() != EdbPVRec::Class()) ePVR = 0;
  if(!ePVR) {printf("Warning: EdbVertexRec::SegmentNeighbor: EdbPVRec not defined, use SetPVRec(...)\n"); return 0;}

  EdbTrackP *tr  = 0, *trown = 0;
  const EdbSegP   *ss  = 0;
  const EdbSegP   *se  = 0;
  int 	    nn   = 0;
  float     distxs, distys, distzs1, distzs, dists;
  float     distxe = 0., distye = 0., distze1= 0., distze = 0., diste = 0.;
  float     xseg = 0, yseg = 0, zseg = 0;
  EdbTrackP *trv = 0;
  EdbVertex *ve = 0;
  float dx = 0., dy = 0.;
  float dz = 0;
  Float_t   Zbin = TMath::Abs(ePVR->GetPattern(1)->Z() - ePVR->GetPattern(0)->Z());
  TObjArray an(200);

	    xseg = s->X();
	    yseg = s->Y();
	    zseg = s->Z();
	    an.Clear();
	    if (s->Track() >= 0)
	    {
		trown = (EdbTrackP *)ePVR->eTracks->At(s->Track()); 
	    }
	    // Select tracks neigborhood
	    int nvn = SelSegNeighbor(s, 0, RadMax, Dpat, &an); 
	    if (trown) an.Add(trown);
	    nvn = an.GetEntries();
	    for (int it=0; it<nvn; it++) {
		    tr = (EdbTrackP*)(an.At(it));
		    if (tr)
		    {
			    if (tr != trown)
			    {
				ss = tr->TrackZmin();
				distxs = (xseg - ss->X());
				distxs *= distxs;
				distys = (yseg - ss->Y());
				distys *= distys;
				distzs1 = (zseg - ss->Z());
				distzs = distzs1*distzs1;
				dists  =  distxs + distys + distzs;
				dists  =  TMath::Sqrt(dists); 
				se = tr->TrackZmax();
				distxe = (xseg - se->X());
				distxe *= distxe;
				distye = (yseg - se->Y());
				distye *= distye;
				distze1 = (zseg - se->Z());
				distze = distze1*distze1;
				diste  =  distxe + distye + distze;
				diste  =  TMath::Sqrt(diste);
				if (diste < dists)
				{
				    if (TMath::Sqrt(distxe+distye) > RadMax)    continue;
				    if (TMath::Abs(distze1)        > Dpat*Zbin) continue;
				    if (Tdistance(*(const EdbSegP *)s, *se)           > ImpMax)    continue;
				}
				else
				{
				    if (TMath::Sqrt(distxs+distys) > RadMax)    continue;
				    if (TMath::Abs(distzs1)        > Dpat*Zbin) continue;
				    if (Tdistance(*(const EdbSegP *)s, *ss)           > ImpMax)    continue;
				}
			    }
			    if (arrt) arrt->Add(tr);
			    nn++;
//--------------------------Add vertex neighborhood
    			    tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
			    ve = tr->VertexS(); 
			    if (ve && ve->Flag() >= -99 && ve->Flag() != -10)
			    {
				dx = ve->VX() - xseg;
				dy = ve->VY() - yseg;
				if (TMath::Sqrt(dx*dx + dy*dy) <= RadMax)
				{
				    dz = ve->VZ() - zseg;
				    if (TMath::Abs(dz) <= Dpat*Zbin)
				    {
					ve->SetFlag(-ve->Flag()-200);
					if (arrv) arrv->Add(ve);
					for(int j=0; j<ve->N(); j++)
					{
					    trv = ve->GetTrack(j);
					    if (trv->MCEvt() < -999) continue;
					    if (ve->Zpos(j) == 0)
					        ss = trv->TrackZmax();
					    else
					        ss = trv->TrackZmin();
					    dx = ss->X() - xseg;
					    dy = ss->Y() - yseg;
					    dz = ss->Z() - zseg;
					    if (TMath::Sqrt(dx*dx+dy*dy) > RadMax)    continue;
					    if (TMath::Abs(dz)        > Dpat*Zbin) continue;
					    //if (TDistance(*(const EdbSegP *)s, *ss)      > ImpMax)    continue;
					    if (arrt) arrt->Add(trv);
    					    trv->SetMC(-trv->MCEvt()-2000, trv->MCTrack());
					}
				    }
				}
			    }
			    ve = tr->VertexE(); 
			    if (ve && ve->Flag() >= -99 && ve->Flag() != -10)
			    {
				dx = ve->VX() - xseg;
				dy = ve->VY() - yseg;
				if (TMath::Sqrt(dx*dx + dy*dy) <= RadMax)
				{
				    dy = ve->VZ() - zseg;
				    if (TMath::Abs(dz) <= Dpat*Zbin)
				    {
					ve->SetFlag(-ve->Flag()-200);
					if (arrv) arrv->Add(ve);
					for(int j=0; j<ve->N(); j++)
					{
					    trv = ve->GetTrack(j);
					    if (trv->MCEvt() < -999) continue;
					    if (ve->Zpos(j) == 0)
					        ss = trv->TrackZmax();
					    else
					        ss = trv->TrackZmin();
					    dx = ss->X() - xseg;
					    dy = ss->Y() - yseg;
					    dz = ss->Z() - zseg;
					    if (TMath::Sqrt(dx*dx+dy*dy) > RadMax) continue;
					    if (TMath::Abs(dz)        > Dpat*Zbin) continue;
					    //if (TDistance(*(const EdbSegP *)s, *ss)      > ImpMax)    continue;
					    if (arrt) arrt->Add(trv);
    					    trv->SetMC(-trv->MCEvt()-2000, trv->MCTrack());
					}
				    }
				}
			    }
		    }
	    }
//----------Restore MCEvt
	    int nv = arrv->GetEntries();
	    for (int i=0; i<nv; i++)
	    {
		ve = (EdbVertex *)arrv->At(i);
		if (ve->Flag() < -99 ) ve->SetFlag(-ve->Flag()-200);
	    }
	    int ntr = arrt->GetEntries();
	    for (int i=0; i<ntr; i++)
	    {
		tr = (EdbTrackP *)arrt->At(i);
		if (tr->MCEvt() < -999 ) tr->SetMC(-tr->MCEvt()-2000, tr->MCTrack());
	    }
	    // Select segments neigborhood
	    //nvn = SelSegNeighbor(s, 1, RadMax, Dpat, arrs); 
	    an.Clear();
	    nvn = SelSegNeighbor(s, 1, RadMax, Dpat, &an); 
	    nvn = an.GetEntries();
	    for (int is=0; is<nvn; is++) {
		    ss = (EdbSegP*)(an.At(is));
		    if (ss)
		    {
			    if (ss != s)
			    {
            if ( ss->W() < SegWmin) continue;
        if (Tdistance(*(const EdbSegP *)s, *ss) > ImpMax) continue;
        arrs->Add((TObject *)ss);
				nn++;
			    }
		    }
	    }
//	    printf("Selected %d segments\n", nvn);
//	    nn += nvn;
  return nn;
} 

//________________________________________________________________________
double EdbVertexRec::Tdistance(const Track& t1, const Track& t2) {
  //
  // geometrical distance between 2 track lines in the space XYZ
  //
  const double a  = t1.tx();
  const double b  = t1.ty();
  const double c  = 1.;
  const double a1 = t2.tx();
  const double b1 = t2.ty();
  const double c1 = 1.;

  const double det = square(a*b1 - b*a1) + square(b*c1 - c*b1) + square(c*a1 - a*c1);
  const SVector<double,3> dx = t2.xvec() - t1.xvec();
  // are tracks parallel?
  if(det==0) return mag(cross(dx,t1.evec()));

  const double det2 = dx[0]*(b*c1 - c*b1) + dx[1]*(c*a1 - a*c1) + dx[2]*(a*b1 - b*a1);
  
  return fabs(det2/sqrt(det));
}

//________________________________________________________________________
double EdbVertexRec::Tdistance(const EdbSegP &s1, const EdbSegP &s2) 
{
  EdbVertex edbv;
  Track t1,t2;
  edbv.Edb2Vt(s1,t1);
  edbv.Edb2Vt(s2,t2);
  return Tdistance(t1,t2);
}

//________________________________________________________________________
double EdbVertexRec::TdistanceChi2(const EdbTrackP &tr1, const EdbTrackP &tr2) 
{
  EdbVertex edbv;
  Track t1,t2;
  float X0 =  ePVR->GetScanCond()->RadX0();
  edbv.Edb2Vt(tr1,t1,X0,tr1.M());
  edbv.Edb2Vt(tr2,t2,X0,tr2.M());
  return distanceChi2(t1,t2);
}

//________________________________________________________________________
double EdbVertexRec::TdistanceChi2(const EdbSegP &s1, const EdbSegP &s2, float m)
{
  EdbVertex edbv;
  Track t1,t2;
  float X0 =  ePVR->GetScanCond()->RadX0();
  edbv.Edb2Vt(s1,t1,X0,m);
  edbv.Edb2Vt(s2,t2,X0,m);
  return distanceChi2(t1,t2);
}

//________________________________________________________________________
EdbTrackP *EdbVertexRec::GetEdbTrack(const int index)
{
  if (eEdbTracks) return (EdbTrackP *)eEdbTracks->At(index);
  else            return 0;
}

//________________________________________________________________________
int EdbVertexRec::CheckTrack(EdbTrackP &track, int zpos)
{
  // loop for all tracks and check if some of them can form the vertex with tr
  int ntr  = eEdbTracks->GetEntries();
  if(ntr<1) return 0;
  int nvtx=0; 

  EdbTrackP *t=0;
  for(int itr=0; itr<ntr; itr++)   {
    t = (EdbTrackP*)(eEdbTracks->At(itr));
    if (t == &track)                                          continue;
    if (t->Flag() < 0)                                        continue;

    EdbVertex *vtx0 = ProbVertex2( &track, t, zpos, 0 );
    if(vtx0) {
      AddVertex(vtx0);
      nvtx++;
    }
    EdbVertex *vtx1 = ProbVertex2( &track, t, zpos, 1 );
    if(vtx1) {
      AddVertex(vtx1);
      nvtx++;
    }
    
  }
  Log(2,"EdbVertexRec::CheckTrack","%d couples found",nvtx);
  return nvtx;
}
