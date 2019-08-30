//-- Author of initial drawing part : Igor Kreslo     27.11.2003
//-- Author                         : Yury Petukhov   25.03.2005
//   Based on AliDisplay class (AliRoot framework - ALICE CERN)
//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDisplay                                                           //
//                                                                      //
// Class to display pattern volume in 3D                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TGButton.h"
#include "TGLabel.h"
#include "TGFont.h"
#include "TGFrame.h"
#include "TGResourcePool.h"
#include "TGClient.h"
#include "TRootEmbeddedCanvas.h"
#include "EdbLog.h"
#include "EdbMomentumEstimator.h"
#include "EdbDisplay.h"

ClassImp(EdbDisplay);
ClassImp(EdbSegG);
ClassImp(EdbTrackG);
ClassImp(EdbVertexG);

//_____________________________________________________________________________
void EdbVertexG::DumpVertex()
{
  if (eV) eV->Print();
}

//_____________________________________________________________________________
void EdbVertexG::InspectVertex()
{
  if (eV) eV->Inspect();
}

//_____________________________________________________________________________
const char *EdbVertexG::GetTitle() const
{
    static char title[80];
    if (eV->V())
	sprintf(title, "Vertex ID %d, Prongs %d, Prob %f", eV->ID(), eV->N(), eV->V()->prob());
    else
	sprintf(title, "Vertex ID %d, Prongs %d, Prob %f", eV->ID(), eV->N(), 1.);
    return title;
}

//_____________________________________________________________________________
const char *EdbVertexG::GetName() const
{
    static char name[] = "Vertex";
    return name;
}

//_____________________________________________________________________________
char *EdbVertexG::GetObjectInfo(int px, int py) const
{
    static char coordinates[80];
    if (eV->V())
	sprintf(coordinates, "X = %.1f, Y = %.1f, Z = %.1f", eV->VX(), eV->VY(), eV->VZ());
    else
	sprintf(coordinates, "X = %.1f, Y = %.1f, Z = %.1f", eV->X(), eV->Y(), eV->Z());
    return coordinates;
}

//_____________________________________________________________________________
void EdbTrackG::DumpTrack()
{
  if (eTr) eTr->Print();
}

//_____________________________________________________________________________
void EdbTrackG::InspectTrack()
{
  if (eTr) eTr->Inspect();
}

//_____________________________________________________________________________
void EdbTrackG::EstimateMomentum()
{
  if (!eTr) return;
  EdbMomentumEstimator me;
  me.PMSang(*eTr);
  me.DrawPlots();
}

//_____________________________________________________________________________
const char *EdbTrackG::GetTitle() const
{
    static char title[80];
    if (eTr) sprintf(title, "Track ID %d, OrigTrackID %d, Nseg %d, Prob %.4f, P %.3f", eTr->ID(), eTr->Track(), eTr->N(), eTr->Prob(), eTr->P());
    else strcpy(title, "Track address not defined");
    return title;
}

//_____________________________________________________________________________
const char *EdbTrackG::GetName() const
{
    static char name[] = "Track";
    return name;
}

//_____________________________________________________________________________
char *EdbTrackG::GetObjectInfo(int px, int py) const
{
    static char coordinates[80];
    float tx = 0., ty = 0., z = 0.;
    int zpos = 1;
    if (GetMarkerColor() == kRed) zpos = 0;
    if( zpos == 0 )
      {
             tx = (eTr->TrackZmax())->TX();
             ty = (eTr->TrackZmax())->TY();
             z  = (eTr->TrackZmax())->Z();
      }
    else
      {
             tx = (eTr->TrackZmin())->TX();
             ty = (eTr->TrackZmin())->TY();
             z  = (eTr->TrackZmin())->Z();
      }
    sprintf(coordinates, "X = %.1f, Y = %.1f, Z = %.1f, TX = %.3f, TY = %.3f", eTr->X(), eTr->Y(), z, tx, ty);
    return coordinates;
}

//_____________________________________________________________________________
EdbSegG::EdbSegG( EdbSegP &s ) : TPolyLine3D(2)
{
  eSeg=0;
  eD = 0;
  float dz = TMath::Abs(s.DZ())/2.; if(dz<0.001) dz=10;
  SetPoint(0,
	   s.X() - s.TX()*dz,
	   s.Y() - s.TY()*dz,
	   s.Z() -           dz);
  SetPoint(1,
	   s.X() + s.TX()*dz,
	   s.Y() + s.TY()*dz,
	   s.Z() +           dz);

  SetLineColor(kRed);
  SetLineWidth(1);
}
//_____________________________________________________________________________
void EdbSegG::DumpSegment()
{
  if (eSeg) eSeg->Print();
}

//_____________________________________________________________________________
void EdbSegG::InspectSegment()
{
  if (eSeg) eSeg->Inspect();
}

//_____________________________________________________________________________
const char *EdbSegG::GetTitle() const
{
    static char title[80];
    if (eSeg) sprintf(title, "Segment ID %d, PID %d, Track %d, MCTrack %d, MCEvent %d", eSeg->ID(), eSeg->PID(), eSeg->Track(), eSeg->MCTrack(), eSeg->MCEvt());
    else strcpy(title, "Segment address not defined");
    return title;
}

//_____________________________________________________________________________
const char *EdbSegG::GetName() const
{
    static char name[] = "Segment";
    return name;
}

//_____________________________________________________________________________
char *EdbSegG::GetObjectInfo(int px, int py) const
{
    static char coordinates[80];
    if (eSeg) sprintf(coordinates, "X = %.1f, Y = %.1f, Z = %.1f, TX = %.3f, TY = %.3f PH = %d", 
		      eSeg->X(), eSeg->Y(), eSeg->Z(), eSeg->TX(), eSeg->TY(), (int)eSeg->W() );
    else strcpy(coordinates, "Segment address not defined");
    return coordinates;
}

//________________________________________________________________________
void EdbDisplay::Set0()
{
  //eVerRec = ((EdbVertexRec *)(gROOT->GetListOfSpecials()->FindObject("EdbVertexRec")));
  //if (eVerRec) if (eVerRec->IsA() != EdbVertexRec::Class()) eVerRec = 0;
  //if (!eVerRec) {printf("Warning: EdbDisplay:Set0: EdbVertexRec not defined, use SetVerRec(...)\n");}
  eArrSegG = 0;
  eArrSegP = 0;
  eArrTr   = 0;
  eArrV    = 0;
  eArrSegPSave = 0;
  eArrTrSave   = 0;
  eArrVSave    = 0;
  eDrawTracks = 1;
  eDrawVertex = 0;
  eColors  = 0;
  eDZs     = 0;
  eVertex  = 0;
  eWorking = 0;
  ePrevious= 0;
  eSegment = 0;
  eTrack   = 0;
  eTrack1   = 0;
  eTrack2   = 0;
  eSegPM = 0;
  eWait_Answer = false;
  eIndVert = -1;
  eIndVertSave = -1;
  eRadMax = 10000.;
  eDpat = 2;
  eImpMax = 1000.;
  eSegWmin = 9;
  eP = 0.;
  eM = 0.;
  eTImpMax = 10000.;
  eTProbMin = 0.0;
  fNumericEntries[0] = 0;
  fNumericEntries[1] = 0;
  fNumericEntries[2] = 0;
  eFromPlate=1;
  eToPlate=57;
  if (fCanvas) fCanvas->Connect("Closed()", "EdbDisplay", this, "Delete()");
}
//________________________________________________________________________
EdbDisplay *EdbDisplay::EdbDisplayExist(const char *title)
{
  TSeqCollection *l = gROOT->GetListOfSpecials();
  if (!l) return 0;
  EdbDisplay *ds = (EdbDisplay*)(l->FindObject(title));
  return ds;
}
//________________________________________________________________________
EdbDisplay::~EdbDisplay()
{
  DeleteModifiedVTX();

  if (gROOT->GetListOfSpecials()->FindObject(this))
    gROOT->GetListOfSpecials()->Remove(this);

  SafeDelete(fCanvasVTX);
  SafeDelete(fCanvasTRK);
  SafeDelete(eWorking);
  SafeDelete(ePrevious);

  if (eArrSegPSave) SafeDelete(eArrSegP);
  if (eArrVSave)    SafeDelete(eArrV);
  if (eArrTrSave)   SafeDelete(eArrTr);

  //SafeDelete(eArrSegG);
}

//________________________________________________________________________
void EdbDisplay::GuessRange(float margZmin,float margZmax,float margR )
{
  // guess minimum display range to fit all objects

  float xmin=0,xmax=1,ymin=0,ymax=1,zmin=0,zmax=1;

  if (eArrSegG) {
    EdbSegG *s=0;
    for(int i=0; i<eArrSegG->GetEntries(); i++) {
      s = (EdbSegG *)(eArrSegG->At(i));
      if(!s) continue;
      if(xmax-xmin<margR) {
        xmax=s->X()+margR;
        xmin=s->X()-margR;
        ymax=s->Y()+margR;
        ymin=s->Y()-margR;
        zmax=s->Z()+margZmax;
        zmin=s->Z()-margZmin;
      }
      if(xmax<s->X()+margR) xmax=s->X()+margR;
      if(xmin>s->X()-margR) xmin=s->X()-margR;
      if(ymax<s->Y()+margR) ymax=s->Y()+margR;
      if(ymin>s->Y()-margR) ymin=s->Y()-margR;
      if(zmax<s->Z()+margZmax) zmax=s->Z()+margZmax;
      if(zmin>s->Z()-margZmin) zmin=s->Z()-margZmin;
    }
  }
  
  if (eArrSegP) {
    EdbSegP *s=0;
    for(int i=0; i<eArrSegP->GetEntries(); i++) {
      s = (EdbSegP *)(eArrSegP->At(i));
      if(!s) continue;
      if(xmax-xmin<margR) {
        xmax=s->X()+margR;
        xmin=s->X()-margR;
        ymax=s->Y()+margR;
        ymin=s->Y()-margR;
        zmax=s->Z()+margZmax;
        zmin=s->Z()-margZmin;
      }
      if(xmax<s->X()+margR) xmax=s->X()+margR;
      if(xmin>s->X()-margR) xmin=s->X()-margR;
      if(ymax<s->Y()+margR) ymax=s->Y()+margR;
      if(ymin>s->Y()-margR) ymin=s->Y()-margR;
      if(zmax<s->Z()+margZmax) zmax=s->Z()+margZmax;
      if(zmin>s->Z()-margZmin) zmin=s->Z()-margZmin;
    }
  }

  if (eArrV) {
    for(int i=0; i<eArrV->GetEntries(); i++) {
      EdbVertex *v = (EdbVertex *)(eArrV->At(i));
      if(!v) continue;
      if(!v->V()) continue;
      if(!v->V()->valid()) continue;
      if(xmax-xmin<margR) {
	xmax=v->VX()+margR;
	xmin=v->VX()-margR;
	ymax=v->VY()+margR;
	ymin=v->VY()-margR;
	zmax=v->VZ()+margZmax;
	zmin=v->VZ()-margZmin;
      }
      if(xmax<v->VX()+margR) xmax=v->VX()+margR;
      if(xmin>v->VX()-margR) xmin=v->VX()-margR;
      if(ymax<v->VY()+margR) ymax=v->VY()+margR;
      if(ymin>v->VY()-margR) ymin=v->VY()-margR;
      if(zmax<v->VZ()+margZmax) zmax=v->VZ()+margZmax;
      if(zmin>v->VZ()-margZmin) zmin=v->VZ()-margZmin;
    }
  }
  
  if (eArrTr) {
    EdbTrackP *t;
    EdbSegP   *s;
    for(int i=0; i<eArrTr->GetEntries(); i++) {
      t = (EdbTrackP *)eArrTr->At(i);
      if(!t) continue;
      for(int j=0; j<t->N(); j++) {
	s = t->GetSegment(j);
	if(!s) continue;
	if(xmax-xmin<margR) {
	  xmax=s->X()+margR;
	  xmin=s->X()-margR;
	  ymax=s->Y()+margR;
	  ymin=s->Y()-margR;
	  zmax=s->Z()+margZmax;
	  zmin=s->Z()-margZmin;
	}
	if(xmax<s->X()+margR) xmax=s->X()+margR;
	if(xmin>s->X()-margR) xmin=s->X()-margR;
	if(ymax<s->Y()+margR) ymax=s->Y()+margR;
	if(ymin>s->Y()-margR) ymin=s->Y()-margR;
	if(zmax<s->Z()+margZmax) zmax=s->Z()+margZmax;
	if(zmin>s->Z()-margZmin) zmin=s->Z()-margZmin;
      }
    }
  }
  if(gEDBDEBUGLEVEL>1) {
    printf("Guess Range:\n");
    printf("X: %10.1f %10.1f\n",xmin,xmax);
    printf("Y: %10.1f %10.1f\n",ymin,ymax);
    printf("Z: %10.1f %10.1f\n",zmin,zmax);
  }
  SetRange(xmin,xmax,ymin,ymax,zmin,zmax);
  fZooms=0;
}

//________________________________________________________________________
void EdbDisplay::Delete()
{
  DeleteModifiedVTX();

  if (gROOT->GetListOfSpecials()->FindObject(this))
    gROOT->GetListOfSpecials()->Remove(this);

  SafeDelete(fCanvasVTX);
  SafeDelete(fCanvasTRK);
  SafeDelete(eWorking);
  SafeDelete(ePrevious);

  if (eArrSegPSave) SafeDelete(eArrSegP);
  if (eArrVSave)    SafeDelete(eArrV);
  if (eArrTrSave)   SafeDelete(eArrTr);
}


//________________________________________________________________________
void EdbDisplay::SetArrSegP(TObjArray *arr)
{
  eArrSegP=arr;

  float x = 0., y = 0., z = 0;
  float wx0 = fVx0, wy0 = fVy0, wz0 = fVz0;
  float wx1 = fVx1, wy1 = fVy1, wz1 = fVz1;

  EdbSegP *seg=0;
  if( eArrSegP ) {
    int nseg = eArrSegP->GetEntries();
    for(int j=0;j<nseg;j++) {
      seg = (EdbSegP*)(eArrSegP->At(j));
      if (seg)
      {
	x = seg->X();
	y = seg->Y();
	z = seg->Z();
	if (fView == 0)
	{
	    wx0 = x < wx0 ? x : wx0;
	    wy0 = y < wy0 ? y : wy0;
	    wz0 = z < wz0 ? z : wz0;
	    wx1 = x > wx1 ? x : wx1;
	    wy1 = y > wy1 ? y : wy1;
	    wz1 = z > wz1 ? z : wz1;
	}
      }
    }
  }
  if (fView == 0)
  {
    if (wx0 < fVx0) fVx0 = wx0 - 500.;
    if (wy0 < fVy0) fVy0 = wy0 - 500.;
    if (wz0 < fVz0) fVz0 = wz0 - 500.;
    if (wx1 > fVx1) fVx1 = wx1 + 500.;
    if (wy1 > fVy1) fVy1 = wy1 + 500.;
    if (wz1 > fVz1) fVz1 = wz1 + 500.;
  }
}
//________________________________________________________________________
void EdbDisplay::SetArrTr(TObjArray *arr)
{
  eArrTr=arr;

  float x = 0., y = 0., z = 0;
  float wx0 = fVx0, wy0 = fVy0, wz0 = fVz0;
  float wx1 = fVx1, wy1 = fVy1, wz1 = fVz1;

  EdbTrackP *tr=0;
  if( eArrTr ) {
    int ntr = eArrTr->GetEntries();
    for(int j=0;j<ntr;j++) {
      tr = (EdbTrackP*)(eArrTr->At(j));
      if(tr)
      {
	if (fView == 0)
	{
	    if (tr->NF())
	    {
		x = tr->TrackZmin()->X();
		y = tr->TrackZmin()->Y();
		z = tr->TrackZmin()->Z();
	    }
	    else
	    {
		eRadMax = 0;
		x = tr->GetSegmentFirst()->X();
		y = tr->GetSegmentFirst()->Y();
		z = tr->GetSegmentFirst()->Z();
	    }
	    wx0 = x < wx0 ? x : wx0;
	    wy0 = y < wy0 ? y : wy0;
	    wz0 = z < wz0 ? z : wz0;
	    wx1 = x > wx1 ? x : wx1;
	    wy1 = y > wy1 ? y : wy1;
	    wz1 = z > wz1 ? z : wz1;
	    if (tr->NF())
	    {
		x = tr->TrackZmax()->X();
		y = tr->TrackZmax()->Y();
		z = tr->TrackZmax()->Z();
	    }
	    else
	    {
		eRadMax = 0;
		x = tr->GetSegmentLast()->X();
		y = tr->GetSegmentLast()->Y();
		z = tr->GetSegmentLast()->Z();
	    }
	    wx0 = x < wx0 ? x : wx0;
	    wy0 = y < wy0 ? y : wy0;
	    wz0 = z < wz0 ? z : wz0;
	    wx1 = x > wx1 ? x : wx1;
	    wy1 = y > wy1 ? y : wy1;
	    wz1 = z > wz1 ? z : wz1;
	}
      }
    }
  }
  if (fView == 0)
  {
    if (wx0 < fVx0) fVx0 = wx0 - 500.;
    if (wy0 < fVy0) fVy0 = wy0 - 500.;
    if (wz0 < fVz0) fVz0 = wz0 - 500.;
    if (wx1 > fVx1) fVx1 = wx1 + 500.;
    if (wy1 > fVy1) fVy1 = wy1 + 500.;
    if (wz1 > fVz1) fVz1 = wz1 + 500.;
  }
}
//________________________________________________________________________
void EdbDisplay::SetArrV(TObjArray *arrv)
{
  eArrV=arrv;

  float x = 0., y = 0., z = 0;
  float wx0 = fVx0, wy0 = fVy0, wz0 = fVz0;
  float wx1 = fVx1, wy1 = fVy1, wz1 = fVz1;

  EdbVertex *v=0;
  if( eArrV ) {
    int nv = eArrV->GetEntries();
    for(int j=0;j<nv;j++) {
      v = (EdbVertex*)(eArrV->At(j));
      if (v)
      {
	if (fView == 0)
	{
	    if (v->V())
	    {
		x = v->VX();
		y = v->VY();
		z = v->VZ();
	    }
	    else
	    {
		eRadMax = 0;
		x = v->X();
		y = v->Y();
		z = v->Z();
	    }
	    wx0 = x < wx0 ? x : wx0;
	    wy0 = y < wy0 ? y : wy0;
	    wz0 = z < wz0 ? z : wz0;
	    wx1 = x > wx1 ? x : wx1;
	    wy1 = y > wy1 ? y : wy1;
	    wz1 = z > wz1 ? z : wz1;
	}
      }
    }
  }
  if (fView == 0)
  {
    if (wx0 < fVx0) fVx0 = wx0 - 500.;
    if (wy0 < fVy0) fVy0 = wy0 - 500.;
    if (wz0 < fVz0) fVz0 = wz0 - 500.;
    if (wx1 > fVx1) fVx1 = wx1 + 500.;
    if (wy1 > fVy1) fVy1 = wy1 + 500.;
    if (wz1 > fVz1) fVz1 = wz1 + 500.;
  }
}
//________________________________________________________________________
void EdbDisplay::Refresh()
{
  EdbSegP *seg=0;
  if( eArrSegP ) {
    int nseg = eArrSegP->GetEntries();
    printf("%d segments to draw...\n",nseg);
    for(int j=0;j<nseg;j++) {
      seg = (EdbSegP*)(eArrSegP->At(j));
      if (seg)
      {
        SegLine(seg)->Draw();
      }
    }
  }


  EdbTrackP *tr=0;
  if( eArrTr ) {
    int ntr = eArrTr->GetEntries();
    printf("%d tracks to draw...\n",ntr);
    for(int j=0;j<ntr;j++) {
      tr = (EdbTrackP*)(eArrTr->At(j));
      if(tr)
      {
        TrackDraw(tr);
      }
    }
  }

  EdbVertex *v=0;
  if( eArrV ) {
    int nv = eArrV->GetEntries();
    printf("%d vertices to draw...\n",nv);
    for(int j=0;j<nv;j++) {
      v = (EdbVertex*)(eArrV->At(j));
      if (v)
      {
        VertexDraw(v);
      }
    }
  }

  if( eArrSegG ) {
    int nseg = eArrSegG->GetEntries();
    printf("%d Graph segments to draw...\n",nseg);
    for(int j=0;j<nseg;j++) {
      EdbSegG *sg = (EdbSegG*)(eArrSegG->At(j));
      if (sg)        sg->Draw();
    }
  }

  if (eSegment)
  {
    eSegPM = new TPolyMarker3D(1);
    eSegPM->SetMarkerStyle(kFullCircle);
    float dz = TMath::Abs(eSegment->DZ()/2.);
    eSegPM->SetPoint(0, 
		        eSegment->X() + eSegment->TX()*dz,
			eSegment->Y() + eSegment->TY()*dz,
			eSegment->Z() +                dz);
    eSegPM->SetMarkerColor(kGreen);
    eSegPM->SetMarkerSize(1.2);
    eSegPM->SetBit(kCannotPick);
    eSegPM->Draw();
  }
}

//________________________________________________________________________
void EdbDisplay::DrawRef(float start[3], float end[3])
{
  // reference vector to indicate the direction (ie scanback track) can be drawn

  TPolyLine3D *line = new TPolyLine3D(2);
  line->SetPoint(0, start[0],start[1], start[2] );
  line->SetPoint(1,   end[0],  end[1],   end[2] );
  line->SetLineColor(kGreen);
  line->SetLineWidth(fLineWidth);
  line->SetBit(kCannotPick);
  line->Draw();

  TPolyMarker3D *mark = new TPolyMarker3D(1);
  mark->SetPoint(0, start[0],start[1], start[2] );
  mark->SetMarkerColor(kGreen);
  mark->SetMarkerStyle(kMultiply);
  mark->SetMarkerSize(1.2);
  mark->Draw();
}

//________________________________________________________________________
void EdbDisplay::VertexDraw(EdbVertex *vv)
{
  float xv,yv,zv;
  if (!vv) return;
  //if (vv->Flag() == -10) {
  Log(2,"EdbDisplay::VertexDraw","id=%d  %d-prong  flag = %d",  vv->ID(), vv->N(), vv->Flag() );
    //return;
    //}
  EdbVertexG *v = new EdbVertexG(this);
  v->SetVertex( vv );

  xv = vv->X();
  yv = vv->Y();
  zv = vv->Z();

  v->SetPoint(0, xv, yv, zv);
  v->SetMarkerStyle(kFullCircle);
  if (fStyle/2 == 1) v->SetMarkerColor(kBlack);
  else               v->SetMarkerColor(kWhite);
//  v->Draw();

  if (vv->V())
  {
//    v = new EdbVertexG(this);
//    v->SetVertex( vv );

    xv = vv->VX();
    yv = vv->VY();
    zv = vv->VZ();
    v->SetPoint(0, xv, yv, zv);
  }
  if (eWorking == vv)
  {
	v->SetMarkerColor(kGreen);
  }
  else if (eVertex == vv)
  {
	if (!eWorking) v->SetMarkerColor(kGreen);
	else
	{
	    if (fStyle/2 == 1) v->SetMarkerColor(kBlack);
	    else               v->SetMarkerColor(kWhite);
	}
  }
  else
  {
	if (fStyle/2 == 1) v->SetMarkerColor(kBlack);
	else               v->SetMarkerColor(kWhite);
  }

  if(vv->Flag()==-10)  v->SetMarkerStyle(kOpenStar);
  else                 v->SetMarkerStyle(kFullStar);
  v->SetMarkerSize(1.2);
  v->Draw();


  if(eDrawVertex>0) {
    TPolyLine3D *line=0;
    const EdbSegP *seg=0;
    EdbTrackP *tr=0;
    bool measured_segment = false;
    if (eVerRec) measured_segment = eVerRec->eUseSegPar;
    for(int i=0; i<vv->N(); i++ ) {
      tr = vv->GetTrack(i);
      if (!(tr->NF())) measured_segment = true;

      seg = vv->GetTrackV(i,measured_segment);         if(!seg) continue;
      line = new TPolyLine3D(2);
      line->SetPoint(0, xv,yv,zv );
      line->SetPoint(1, seg->X(), seg->Y(), seg->Z());
      if (fStyle/2 == 1) line->SetLineColor(kBlack);
      else               line->SetLineColor(kYellow);
      line->SetLineWidth(2);
      line->SetBit(kCannotPick);
      line->Draw();
      if(eDrawTracks>10)  DrawSegmentExtrapolationLine( *seg, seg->Z(), zv );
    }
    for(int i=0; i<vv->Nn(); i++ ) {    //draw auxillary tracks if any
      EdbTrackP *t = vv->GetTrackN(i);
      TrackDraw(t);
      if(eDrawTracks>10)  DrawSegmentExtrapolationLine( *((EdbSegP *)t), t->Z(), zv );
    }
  }
}

//________________________________________________________________________
void EdbDisplay::DrawSegmentExtrapolationLine(const EdbSegP &s, float zmin, float zmax)
{
   EdbSegP start(s);
   EdbSegP end(s);
   start.PropagateTo(zmin);
   start.PropagateTo(zmax);
   TPolyLine3D *trline = new TPolyLine3D(2);
   trline->SetPoint(0, start.X(),start.Y(),start.Z() );
   trline->SetPoint( 1,   end.X(),  end.Y(),  end.Z() );
   if (fStyle/2 == 1) trline->SetLineColor(kBlack);
   else               trline->SetLineColor(kWhite);
   trline->SetLineWidth(1);
   trline->SetBit(kCannotPick);
   trline->Draw();
   Log(3,"EdbDisplay::DrawSegmentExtrapolationLine", "(%f %f %f) -> (%f %f %f)", start.X(),start.Y(),start.Z(), end.X(),  end.Y(),  end.Z());
}

//________________________________________________________________________
void EdbDisplay::TrackDraw(EdbTrackP *tr, Color_t kColor)
{
  // eDrawTracks: 1 - draw fitted track dotted line only
  //              2 - draw also white marker at zmin
  //              3 - draw also red marker at zmax
  //              4 - draw also measured segments
  //              5 - draw fitted line and  measured segments, no markers 
  //              6 - draw measured segments only
  //              7 - draw measured segments only
  //              8 - draw only solid white track line
  //             14 - as "4" plus track extrapolation line

  if (!tr) return;
  TPolyLine3D *line=0;
  const EdbSegP *seg=0;
  float dz = 0.;

  if(eDrawTracks%10>0 && eDrawTracks%10<6) {   // only dotted track line
    line = new TPolyLine3D(tr->N());
    if (tr->NF())    {
      for(int is=0; is<tr->NF(); is++) {
	seg = tr->GetSegmentF(is);
	if(seg) 
	{
	    if (is == 0)
	    {
		dz = TMath::Abs(seg->DZ());
		line->SetPoint(is, seg->X()+seg->TX()*dz,
				   seg->Y()+seg->TY()*dz,
				   seg->Z()+          dz );
	    }
	    else
	    {
		line->SetPoint(is, seg->X(), seg->Y(), seg->Z() );
	    }
	}
      }
    }
    else      {
      for(int is=0; is<tr->N(); is++) {
	seg = tr->GetSegment(is);
	if(seg) 
	{
	    if (is == 0)
	    {
		dz = TMath::Abs(seg->DZ());
		line->SetPoint(is, seg->X()+seg->TX()*dz,
				   seg->Y()+seg->TY()*dz,
				   seg->Z()+          dz);
	    }
	    else
	    {
		line->SetPoint(is, seg->X(), seg->Y(), seg->Z() );
	    }
	}
      }
    }
    if (fStyle/2 == 1) line->SetLineColor(kBlack);
    else               line->SetLineColor(kColor);
    line->SetLineWidth(fLineWidth);
    line->SetBit(kCannotPick);
    if (tr->Flag() != -10)
	line->SetLineStyle(3);
    else
	line->SetLineStyle(4);
    line->Draw();
  }

  if(eDrawTracks%10>3 && eDrawTracks%10<8 && tr->N()>0 ) {
    for(int is=0; is<tr->N(); is++) {
      seg = tr->GetSegment(is);
      if (seg) SegLine(seg)->Draw();
    }
  }

  if(eDrawTracks%10>7) {
    line = new TPolyLine3D(tr->N());
    if (tr->NF()){
      for(int is=0; is<tr->NF(); is++) {
	seg = tr->GetSegmentF(is);
	if(seg) line->SetPoint(is, seg->X(), seg->Y(), seg->Z() );
      }
    }    else    {
      for(int is=0; is<tr->N(); is++) {
	seg = tr->GetSegment(is);
	if(seg) line->SetPoint(is, seg->X(), seg->Y(), seg->Z() );
      }
    }
    if (fStyle/2 == 1) line->SetLineColor(kBlack);
    else               line->SetLineColor(kColor);
    line->SetLineWidth(fLineWidth);
    line->SetLineStyle(1);
    line->Draw();
  }

  if(eDrawTracks%10>1 && eDrawTracks%10<5) {
    EdbTrackG *pms = new EdbTrackG(1, this);
    pms->SetTrack( tr );
    pms->SetMarkerStyle(kOpenCircle);

    if (tr->NF())      seg = tr->TrackZmin();
    else               seg = tr->GetSegmentFirst();
    if(seg) {
      pms->SetPoint(0, seg->X(), seg->Y(), seg->Z() );
      if (fStyle/2 == 1) pms->SetMarkerColor(kBlack);
      else               pms->SetMarkerColor(kColor);
      pms->SetMarkerSize(1.2);
      pms->Draw();
      //printf("White track edge, Id %d, Nseg %d.\n", tr->ID(), tr->N());
    }
  }

  if(eDrawTracks%10>2 && eDrawTracks%10<5) {
    EdbTrackG *pme = new EdbTrackG(1, this);
    pme->SetTrack( tr );
    pme->SetMarkerStyle(kOpenCircle);

    if (tr->NF())      seg = tr->TrackZmax();
    else               seg = tr->GetSegmentLast();
    if(seg) {
      pme->SetPoint(0, seg->X(), seg->Y(), seg->Z() );
//      dz = TMath::Abs(seg->DZ());
//      pme->SetPoint(0, seg->X()+seg->TX()*dz,
//    		       seg->Y()+seg->TY()*dz,
//		       seg->Z()+          dz);
      pme->SetMarkerColor(kRed);
      pme->SetMarkerSize(1.2);
      pme->Draw();
    }
  }

}

///------------------------------------------------------------
EdbSegG *EdbDisplay::SegLine(const EdbSegP *seg)
{
  if (!seg) return 0;
  float dz = TMath::Abs(seg->DZ())/2.;
  EdbSegG *line = new EdbSegG(2, this);
  //line->SetPoint(0, seg->X(), seg->Y(), seg->Z()-dz );
  line->SetPoint(0,
                 seg->X() - seg->TX()*dz,
                 seg->Y() - seg->TY()*dz,
                 seg->Z() -           dz);
  line->SetPoint(1,
                 seg->X() + seg->TX()*dz,
                 seg->Y() + seg->TY()*dz,
                 seg->Z() +           dz);

  int npieces=eToPlate-eFromPlate;
  line->SetLineColor(gStyle->GetColorPalette(150+int(46.*(1.-1.*(seg->PID()-eFromPlate)/npieces))));  
  //line->SetLineColor( 2+seg->PID()%7 );
  Width_t lwf = int(seg->W()/10.);
  if (lwf > 3) lwf = 2;
  line->SetLineWidth(lwf*fLineWidth);
  if (seg->Track() <= -2) line->SetLineWidth(3*lwf*fLineWidth);
  line->SetSeg(seg);
  return line;
}

//_____________________________________________________________________________
void EdbSegG::AddAsTrackToVertex()
{
  char text[512];
  EdbVTA *vta = 0;
//  if (eD) if (!(eD->eVerRec)) eD->eVerRec = ((EdbVertexRec *)(gROOT->GetListOfSpecials()->FindObject("EdbVertexRec")));
  if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
  if (!eD->eVerRec) {printf("Error: EdbDisplay:AddAsTrackToVertex: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
  if (eSeg && eD)
  {
    if (eD->eWait_Answer) return;
    if (!(eD->eVertex))
    {
    
	    printf("No working vertex selected!\n");
	    fflush(stdout);
	    return;
    }
    if (eSeg->Track() >= 0)
    {
    
	    printf("Segment already belong to a track!\n");
	    fflush(stdout);
	    return;
    }
    EdbVertex *ePreviousSaved = eD->ePrevious;
    if (eD->eWorking == 0)
    {
	eD->eWorking = new EdbVertex();
	int ntr = eD->eVertex->N();
	int i = 0, n = 0;
	for(i=0; i<ntr; i++)
	{
	    if ((vta = (eD->eVerRec)->AddTrack(*(eD->eWorking), (eD->eVertex)->GetTrack(i), (eD->eVertex)->Zpos(i))))
	    {
		(eD->eVertex)->GetTrack(i)->AddVTA(vta);
		n++;
	    }
	}
	if (n < 2)
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}

	if (!((eD->eVerRec)->MakeV(*(eD->eWorking))))
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}
    }
    else
    {
//	if (ePrevious) delete ePrevious;
	eD->ePrevious = eD->eWorking;
	eD->eWorking = new EdbVertex();
	int ntr = eD->ePrevious->N();
	int i = 0, n = 0;
	for(i=0; i<ntr; i++)
	{
	    if ((vta = (eD->eVerRec)->AddTrack(*(eD->eWorking),(eD->ePrevious)->GetTrack(i), (eD->ePrevious)->Zpos(i))))
	    {
		((eD->ePrevious)->GetTrack(i))->AddVTA(vta);
		n++;
	    }
	}
	if (n < 2)
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
		eD->ePrevious = ePreviousSaved;
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}

	if (!((eD->eVerRec)->MakeV(*(eD->eWorking))))
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
		eD->ePrevious = ePreviousSaved;
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}
    }
    float mass = 0.139; //pion
    float momentum = 0.1; // 100 MeV
    if (eD->eP > 0.) momentum = eD->eP;
    if (eD->eM > 0.) mass = eD->eM;
    EdbTrackP *Tr = new EdbTrackP((EdbSegP *)eSeg, mass);
    Tr->SetP(momentum);
    Tr->FitTrackKFS();
    float ImpMaxSave = (eD->eVerRec)->eImpMax;
    (eD->eVerRec)->eImpMax = eD->eTImpMax;
    float ProbMinSave = (eD->eVerRec)->eProbMin;
    (eD->eVerRec)->eProbMin = eD->eTProbMin;
//    printf(" seg id %d x %f y %f z %f\n", eSeg->ID(), eSeg->X(), eSeg->Y(), eSeg->Z());
    if ((vta = (eD->eVerRec)->AddTrack(*(eD->eWorking), Tr, 1)))
    {
	if ( Tr->Z() >= (eD->eWorking)->VZ() ) vta->SetZpos(1);
	else vta->SetZpos(0);
//	printf(" tr  id %d x %f y %f z %f\n", Tr->ID(), Tr->X(), Tr->Y(), Tr->Z());
	Tr->AddVTA(vta);
//	(eD->eWorking)->ResetTracks();
	(eD->eArrTr)->Add(Tr);
	EdbVertex *eW = eD->eWorking;
	eW->SetID(eD->eVertex->ID());
	eW->V()->rmsDistAngle();
	sprintf(text,"New     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	eD->DrawOldBut("Original");
	if (eD->ePrevious)
	{
	    eD->DrawNewVTX(text);
	    eD->DrawNewBut("Modified");
	    eW = eD->ePrevious;
	    eW->V()->rmsDistAngle();
	    sprintf(text,"Pre     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	    eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	    eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	    eD->DrawPreVTX(text);
	    eD->DrawPreBut("Previous");
	}
	else
	{
	    eD->DrawPreVTX(text);
	    eD->DrawPreBut("Modified");
	    if (eD->eVertex->ID() >= 0)
	    {
		eD->DrawAcc();
		eD->DrawCan();
		eD->DrawUnd();
	    }
	}
//	eD->DrawOldBut("Original");
	eD->DrawVTXTracks("Modified", eD->eWorking);
	if (eD->eIndVert >= 0) 
	{
//	    eD->eArrV->RemoveAt(eD->eIndVert);
	    eD->eArrV->AddAt(eD->eWorking, eD->eIndVert);
	}	
	(eD->eCreatedTracks).Add(Tr);
	eD->Draw();
	if (ePreviousSaved) delete ePreviousSaved;
	ePreviousSaved = 0;
	(eD->eVerRec)->eImpMax = ImpMaxSave;
	(eD->eVerRec)->eProbMin = ProbMinSave;
    }
    else
    {
	printf("Track not added! May be Prob < ProbMin=%f. Change ProbMin with 'TrackParams' button!\n", eD->eTProbMin);
	fflush(stdout);
	delete Tr;
	delete eD->eWorking;
	if (eD->ePrevious)
	{
	    eD->eWorking = eD->ePrevious;
	    (eD->eWorking)->ResetTracks();
	    eD->ePrevious = ePreviousSaved;
	}
	else
	{
	    eD->eWorking = 0;
	    (eD->eVertex)->ResetTracks();
	}
	(eD->eVerRec)->eImpMax = ImpMaxSave;
	(eD->eVerRec)->eProbMin = ProbMinSave;
	return;
    }
  }
}
//_____________________________________________________________________________
void EdbSegG::AddToNewTrack()
{
  if (eSeg && eD)
  {
	if (eSeg->Track() >= 0)
	{
	    printf("This segment alredy belong to a track!\n");
	    fflush(stdout);
	    return;
	}
	if (eD->eTrack)
	{
	    if ((eD->eTrack)->NF())
	    {
		delete eD->eTrack;
		eD->eTrack = 0;
	    }
	}
	if (!eD->eTrack) eD->eTrack = new EdbTrackP();
	(eD->eTrack)->AddSegment((EdbSegP *)eSeg);
  }
}
//_____________________________________________________________________________
void EdbSegG::AddToNewTrackAndFit()
{
  if (eSeg && eD)
  {
	if (eSeg->Track() >= 0)
	{
	    printf("This segment already belong to a track!\n");
	    fflush(stdout);
	    return;
	}
	if (!eD->eTrack) eD->eTrack = new EdbTrackP();
	(eD->eTrack)->AddSegment((EdbSegP *)eSeg);
	float mass = 0.139; //pion
	float momentum = 1.; // 1000 MeV
	if (eD->eP > 0.) momentum = eD->eP;
	if (eD->eM > 0.) mass = eD->eM;
	eD->eTrack->SetM(mass);
	eD->eTrack->SetP(momentum);
	float X0 = 0.;
	if (eD->eVerRec) if ((eD->eVerRec)->ePVR) X0 = (((eD->eVerRec)->ePVR)->GetScanCond())->RadX0();
	(eD->eTrack)->FitTrackKFS(true, X0, 0);
	if (!eD->eArrTr) eD->eArrTr = new TObjArray();
	(eD->eArrTr)->Add(eD->eTrack);
	if (eD->eArrTrSave) (eD->eArrTrSave)->Add(eD->eTrack);
	//if (eD->eArrSegP) (eD->eArrSegP)->Remove((TObject *)eSeg);
	//if (eD->eArrSegPSave) (eD->eArrSegPSave)->Remove((TObject *)eSeg);
	eD->Draw();
  }
}
//_____________________________________________________________________________
void EdbSegG::RemoveFromTrack()
{
  if (eSeg && eD)
  {
	int ind = 0;
	if ((ind = eSeg->Track()) < 0)
	{
	    printf("This segment not included in a track!\n");
	    fflush(stdout);
	    return;
	}
	if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
	if (!eD->eVerRec) {printf("Error: EdbDisplay:RemoveFromTrack: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
	TObjArray *etr = 0;
	if (eD->eVerRec) etr = (eD->eVerRec)->eEdbTracks;
	if (!etr)
	{
	    printf("No information about tracks array!\n");
	    fflush(stdout);
	    return;
	}
	int trind = etr->GetEntries();
	if (ind >= trind)
	{
	    printf("Wrong track index in segment!\n");
	    fflush(stdout);
	    return;
	}
	EdbTrackP *tr = (EdbTrackP *)(etr->At(ind));
	if (!tr)
	{
	    printf("Wrong track address in array!\n");
	    fflush(stdout);
	    return;
	}
	if (tr->N() < 2)
	{
	    printf("Only one segment in track - delete track insted!\n");
	    fflush(stdout);
	    return;
	}
	if (tr->VTAS() || tr->VTAE())
	{
	    printf("Track belong to a vertex - impossible operate with it!\n");
	    fflush(stdout);
	    return;
	}
	tr->RemoveSegment((EdbSegP *)eSeg);
	float X0 = 0.;
	if (eD->eVerRec) if ((eD->eVerRec)->ePVR) X0 = (((eD->eVerRec)->ePVR)->GetScanCond())->RadX0();
	tr->FitTrackKFS(true, X0, 0);
	if (!(eD->eArrSegP)) eD->eArrSegP = new TObjArray();
	if(!((eD->eArrSegP)->FindObject(eSeg))) (eD->eArrSegP)->Add((EdbSegP *)eSeg);
	if (eD->eArrTrSave)
	{
	    if (!(eD->eArrSegPSave)) eD->eArrSegPSave = new TObjArray();
	    if(!((eD->eArrSegPSave)->FindObject(eSeg))) eD->eArrSegPSave->Add((EdbSegP *)eSeg);
	}
	eD->Draw();
  }
}

//_____________________________________________________________________________
void EdbSegG::SplitTrack()
{
  if (eSeg && eD)
  {
	int ind = 0;
	if ((ind = eSeg->Track()) < 0)
	{
	    printf("This segment not included in a track!\n");
	    fflush(stdout);
	    return;
	}
	if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
	if (!eD->eVerRec) {printf("Error: EdbDisplay:SplitTrack: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
	TObjArray *etr = 0;
	if (eD->eVerRec) etr = (eD->eVerRec)->eEdbTracks;
	if (!etr)
	{
	    printf("No information about tracks array!\n");
	    fflush(stdout);
	    return;
	}
	int trind = etr->GetEntries();
	if (ind >= trind)
	{
	    printf("Wrong track index in segment!\n");
	    fflush(stdout);
	    return;
	}
	EdbTrackP *tr = (EdbTrackP *)(etr->At(ind));
	if (!tr)
	{
	    printf("Wrong track address in array!\n");
	    fflush(stdout);
	    return;
	}
	if (tr->N() < 2)
	{
	    printf("Only one segment in track - delete track insted!\n");
	    fflush(stdout);
	    return;
	}
	if (tr->VTAS() || tr->VTAE())
	{
	    printf("Track belong to a vertex - impossible operate with it!\n");
	    fflush(stdout);
	    return;
	}
	if (eD->eTrack1)
	{
	    printf("Intermediate track already exist - fix it!\n");
	    fflush(stdout);
	    return;
	}
	if (eD->eTrack2)
	{
	    printf("Intermediate track already exist - fix it!\n");
	    fflush(stdout);
	    return;
	}
 	EdbSegP *seg = 0;
	eD->eTrack1 = new EdbTrackP();
	eD->eTrack2 = new EdbTrackP();
    	for(int is=0; is<tr->N(); is++) {
	    seg = tr->GetSegment(is);
	    if (seg->Z() <= eSeg->Z()) (eD->eTrack1)->AddSegment((EdbSegP *)seg);
	    else                       (eD->eTrack2)->AddSegment((EdbSegP *)seg);
	}
	tr->SetFlag(-10);
	tr->SetSegmentsTrack(-2-(tr->ID()+1));
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(tr))
	{
	    eD->eArrTr->Remove(tr);
	    eD->eArrTr->Compress();
	    if (!(eD->eArrSegP)) eD->eArrSegP = new TObjArray();
    	    for(int is=0; is<tr->N(); is++) {
		seg = tr->GetSegment(is);
		if(eD->eArrSegP) if(!((eD->eArrSegP)->FindObject(seg))) eD->eArrSegP->Add(seg);
	    }
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(tr))
	{
	    eD->eArrTrSave->Remove(tr);
	    eD->eArrTrSave->Compress();
	    if (!(eD->eArrSegPSave)) eD->eArrSegPSave = new TObjArray();
    	    for(int is=0; is<tr->N(); is++) {
		seg = tr->GetSegment(is);
		if(eD->eArrSegPSave) if(!((eD->eArrSegPSave)->FindObject(seg))) eD->eArrSegPSave->Add(seg);
	    }
	}
	eD->Draw();
	float mass = tr->M();
	float momentum = tr->P();
	if (eD->eP > 0.) momentum = eD->eP;
	if (eD->eM > 0.) mass = eD->eM;

	eD->eTrack1->SetM(mass);
	eD->eTrack1->SetP(momentum);
	float X0 = 0.;
	if (eD->eVerRec) if ((eD->eVerRec)->ePVR) X0 = (((eD->eVerRec)->ePVR)->GetScanCond())->RadX0();
	eD->eTrack1->FitTrackKFS(true, X0, 0);
	if (!eD->eArrTr) eD->eArrTr = new TObjArray();
	(eD->eArrTr)->Add(eD->eTrack1);
	if (eD->eArrTrSave) (eD->eArrTrSave)->Add(eD->eTrack1);

	eD->eTrack2->SetM(mass);
	eD->eTrack2->SetP(momentum);
	eD->eTrack2->FitTrackKFS(true, X0, 0);
	(eD->eArrTr)->Add(eD->eTrack2);
	if (eD->eArrTrSave) (eD->eArrTrSave)->Add(eD->eTrack2);

	eD->Draw();
  }
}
//_____________________________________________________________________________
void EdbSegG::SetAsWorking()
{
  EdbDisplay *eDs = 0;
  EdbSegP    *eSegs = 0;
  eDs = eD;
  eSegs = (EdbSegP *)eSeg;
  if (eDs && eSegs)
  {
    if (eDs->eSegment == eSegs) return;
    if (eDs->eWait_Answer) return;
    if (eDs->eVertex)
    {
	if (eDs->eWorking)
	{
	    eDs->DialogModifiedVTX();
	    return;
	}
	else 
	{
	    eDs->CancelModifiedVTX();
	}
    }
    if (eDs->eSegPM)
    { 
	if (eDs->fPad->GetListOfPrimitives()->FindObject(eDs->eSegPM))
	{
	    eDs->fPad->GetListOfPrimitives()->Remove(eDs->eSegPM);
	}
	delete eDs->eSegPM;
	eDs->eSegPM = 0;
    }
    if (!(eDs->eArrSegP)) 
    {
	eDs->eArrSegP = new TObjArray(20);
	eDs->eArrSegP->Add((TObject *)eSegs);
	eDs->Draw();
    }
    else
    { 

	if (!(eDs->eArrSegP->FindObject(eSegs)))
	{
	    eDs->eArrSegP->Add((TObject *)eSegs);
	    eDs->Draw();
	}
    }
    eDs->eSegment = (EdbSegP *)eSegs;
    eDs->DrawEnv();
    eDs->eSegPM = new TPolyMarker3D(1);
    eDs->eSegPM->SetMarkerStyle(kFullCircle);
    if(eSegs) {
      float dz = TMath::Abs(eSegs->DZ()/2.);
      eDs->eSegPM->SetPoint(0, 
    		    eSegs->X() + eSegs->TX()*dz,
		    eSegs->Y() + eSegs->TY()*dz,
		    eSegs->Z() +            dz);
      eDs->eSegPM->SetMarkerColor(kGreen);
      eDs->eSegPM->SetMarkerSize(1.2);
      eDs->eSegPM->AppendPad();
      eDs->eSegPM->Draw();
    }
  }
}
//_____________________________________________________________________________
void EdbSegG::InfoSegVert()
{
  if (!(eD->eVertex))
  {
    
	    printf("No working vertex selected!\n");
	    fflush(stdout);
	    return;
  }
  int zpos = 1;
  EdbVertex *v = eD->eVertex;
  if (eD->eWorking) v = eD->eWorking;
  char CanvasTRKName[140];
  strcpy(CanvasTRKName, "TRK-");
  strcat(CanvasTRKName, (eD->fCanvas)->GetName());
  if ((eD->fCanvasTRK = (TCanvas *)(gROOT->GetListOfCanvases()->FindObject(CanvasTRKName))))
  {
    (eD->fCanvasTRK)->SetTitle("Segment - Vertex relation parameters");
    (eD->fCanvasTRK)->Clear();
    (eD->fCanvasTRK)->Modified();
    (eD->fCanvasTRK)->Update();
  }
  else
  {
    int xpos = (eD->fCanvas)->GetWindowTopX()+(eD->fCanvas)->GetWw();
    int ypos = (eD->fCanvas)->GetWindowTopY();
    eD->fCanvasTRK = new TCanvas(CanvasTRKName, "Segment - Vertex relation parameters",
			     -xpos, ypos, 640, 330);
    (eD->fCanvasTRK)->ToggleEventStatus();
  }
  if (eD->fVTXTRKInfo)
  {
    (eD->fVTXTRKInfo)->Clear();
  }
  else
  {
    eD->fVTXTRKInfo = new TPaveText(0.05, 0.05, 0.95, 0.95);
    (eD->fVTXTRKInfo)->ResetBit(kCanDelete);
  }
  char line[128];
  EdbSegP *s = (EdbSegP *)eSeg;  
  TText *t = 0;

  strcpy(line, "  Segment   ID          X          Y          Z          TX         TY        PH");
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  sprintf(line,"            %-4d        %-8.1f   %-8.1f   %-8.1f   %-7.4f    %-7.4f    %-4d",
	  s->ID(), s->X(), s->Y(), s->Z(), s->TX(), s->TY(), (int)s->W() );
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

//  t = (eD->fVTXTRKInfo)->AddText("");
//  t->SetTextColor(kBlack);
//  t->SetTextSize(0.03);
//  t->SetTextAlign(12);
//  t->SetTextFont(102);

  strcpy(line, "  Vertex    ID    Mult  X          Y          Z          Dist   Chi2     Prob");
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  sprintf(line,"            %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
    v->ID(), v->N(), v->VX(), v->VY(), v->VZ(), v->V()->dist(),
    v->V()->chi2()/v->V()->ndf(), v->V()->prob());
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

//  t = (eD->fVTXTRKInfo)->AddText("");
//  t->SetTextColor(kBlack);
//  t->SetTextSize(0.03);
//  t->SetTextAlign(12);
//  t->SetTextFont(102);

  float dx   = v->VX() - s->X();
  float dy   = v->VY() - s->Y();
  float dz   = v->VZ() - s->Z();
  float dist = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

  float mass = 0.139; //pion
  float momentum = 0.1; // 100 MeV
  if (eD->eP > 0.) momentum = eD->eP;
  if (eD->eM > 0.) mass = eD->eM;
  EdbTrackP *tr = new EdbTrackP((EdbSegP *)eSeg, mass);
  tr->SetP(momentum);
  tr->FitTrackKFS();

  float X0 = 0.;
  if (eD->eVerRec) if ((eD->eVerRec)->ePVR)  X0 = (((eD->eVerRec)->ePVR)->GetScanCond())->RadX0();
  float chi2 = v->Chi2Track(tr, zpos, X0);
  float impa = v->DistTrack(tr, zpos, X0);
  delete tr;
  sprintf(line, "  Segment - Vertex  impact = %-6.1f, chi2 = %-7.1f, distance = %-8.1f", impa, chi2, dist);
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kRed);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  t = (eD->fVTXTRKInfo)->AddText("");
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  (eD->fVTXTRKInfo)->Draw();
  (eD->fCanvasTRK)->Modified();
  (eD->fCanvasTRK)->Update();
}

//_____________________________________________________________________________
void EdbSegG::InfoSegSeg()
{
  if (!(eD->eSegment))
  {
    
	    printf("No working segment selected!\n");
	    fflush(stdout);
	    return;
  }
  EdbSegP *s = eD->eSegment;
  char CanvasTRKName[140];
  strcpy(CanvasTRKName, "SEG-");
  strcat(CanvasTRKName, (eD->fCanvas)->GetName());
  if ((eD->fCanvasTRK = (TCanvas *)(gROOT->GetListOfCanvases()->FindObject(CanvasTRKName))))
  {
    (eD->fCanvasTRK)->SetTitle("Segment - Segment relation parameters");
    (eD->fCanvasTRK)->Clear();
    (eD->fCanvasTRK)->Modified();
    (eD->fCanvasTRK)->Update();
  }
  else
  {
    int xpos = (eD->fCanvas)->GetWindowTopX()+(eD->fCanvas)->GetWw();
    int ypos = (eD->fCanvas)->GetWindowTopY();
    eD->fCanvasTRK = new TCanvas(CanvasTRKName, "Segment - Segment relation parameters",
			     -xpos, ypos, 640, 330);
    (eD->fCanvasTRK)->ToggleEventStatus();
  }
  if (eD->fVTXTRKInfo)
  {
    (eD->fVTXTRKInfo)->Clear();
  }
  else
  {
    eD->fVTXTRKInfo = new TPaveText(0.05, 0.05, 0.95, 0.95);
    (eD->fVTXTRKInfo)->ResetBit(kCanDelete);
  }
  char line[128];
  TText *t = 0;

  strcpy(line, "  Segment   ID          X          Y          Z          TX         TY");
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  sprintf(line,"            %-4d        %-8.1f   %-8.1f   %-8.1f   %-7.4f    %-7.4f",
		      s->ID(), s->X(), s->Y(), s->Z(),
		      s->TX(), s->TY());
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

//  t = (eD->fVTXTRKInfo)->AddText("");
//  t->SetTextColor(kBlack);
//  t->SetTextSize(0.03);
//  t->SetTextAlign(12);
//  t->SetTextFont(102);

  strcpy(line, "  Segment   ID          X          Y          Z          TX         TY");
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  sprintf(line,"            %-4d        %-8.1f   %-8.1f   %-8.1f   %-7.4f    %-7.4f",
		      eSeg->ID(), eSeg->X(), eSeg->Y(), eSeg->Z(),
		      eSeg->TX(), eSeg->TY());
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);


//  t = (eD->fVTXTRKInfo)->AddText("");
//  t->SetTextColor(kBlack);
//  t->SetTextSize(0.03);
//  t->SetTextAlign(12);
//  t->SetTextFont(102);

  float dx   = eSeg->X() - s->X();
  float dy   = eSeg->Y() - s->Y();
  float dz   = eSeg->Z() - s->Z();
  float dist = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

  float mass = 0.139; //pion
  float momentum = .1; // 100 MeV
  if (eD->eP > 0.) momentum = eD->eP;
  if (eD->eM > 0.) mass = eD->eM;
  EdbTrackP *tr1 = new EdbTrackP((EdbSegP *)s, mass);
  tr1->SetP(momentum);
  tr1->FitTrackKFS();
  EdbTrackP *tr2 = new EdbTrackP((EdbSegP *)eSeg, mass);
  tr2->SetP(momentum);
  tr2->FitTrackKFS();
  EdbSegP *seg1 = (EdbSegP *)(tr1->TrackZmax());
  seg1->SetP(momentum);
  EdbSegP *seg2 = (EdbSegP *)(tr2->TrackZmax());
  seg2->SetP(momentum);
  if (eD->eVerRec)
  {
    float chi2 = (eD->eVerRec)->TdistanceChi2(*seg1, *seg2, mass);
    float impa = (eD->eVerRec)->Tdistance(*s, *eSeg);
    sprintf(line, "  Segment - Segment  impact = %-6.1f, chi2 = %-7.1f, distance = %-8.1f", impa, chi2, dist);
  }
  else
  {
    sprintf(line, "  Impossible to calculate impact and chi2 - No eVerRec defined!");
    printf("Error: EdbDisplay:AddAsTrackToVertex: EdbVertexRec not defined, use SetVerRec(...)\n");
    fflush(stdout);
  }
  delete tr1;
  delete tr2;
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kRed);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  t = (eD->fVTXTRKInfo)->AddText("");
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  (eD->fVTXTRKInfo)->Draw();
  (eD->fCanvasTRK)->Modified();
  (eD->fCanvasTRK)->Update();
}
//_____________________________________________________________________________
void EdbVertexG::TestVertex()
{
  printf("\n");
  eV->Print();
  EdbTrackP *t=0;
  EdbMomentumEstimator me;
  for(int i=0; i<eV->N(); i++) {
    t = eV->GetTrack(i);
    //t->PrintNice();
    printf("nseg = %d p = %f\n",t->N(), (float)(me.PMSang(*t)) );
  }
}
//_____________________________________________________________________________
void EdbVertexG::SetAsWorking()
{
  char text[512];
  EdbTrackP *tr = 0;
  EdbDisplay *eDs = 0;
  EdbVertex  *eVs = 0;
  eDs = eD;
  eVs = eV;
  if (eDs && eVs)
  {
    if (eDs->eVertex == eVs) return;
    if (eDs->eWait_Answer) return;
    if (eDs->eWorking)
    {
	eDs->DialogModifiedVTX();
	return;
    }
    if (eDs->eVertex)
    {
    	eDs->CancelModifiedVTX();
    }
    if (eDs->eSegPM)
    { 
	eDs->ClearSegmentEnv();
    }
    SetMarkerColor(kGreen);
    if (!(eDs->eArrV)) 
    {
	eDs->eArrV = new TObjArray(20);
	eDs->eArrV->Add((TObject *)eVs);
	if (!(eDs->eArrTr))  eDs->eArrTr = new TObjArray(20);
	for (int i=0; i<eVs->N(); i++)
	{
	    tr = eVs->GetTrack(i);
	    if(!(eDs->eArrTr->FindObject(tr))) eDs->eArrTr->Add(tr);
	}
	eDs->Draw();
    }
    else
    {
	if (!((eDs->eArrV)->FindObject(eVs)))
	{
	    eDs->eArrV->Add(eVs);
	    if (!(eDs->eArrTr))  eDs->eArrTr = new TObjArray(20);
	    for (int i=0; i<eVs->N(); i++)
	    {
		tr = eVs->GetTrack(i);
		if(!(eDs->eArrTr->FindObject(tr))) eDs->eArrTr->Add(tr);
	    }
	    eDs->eVertex = 0;
	    eDs->Draw();
	}
    }
    eDs->eSegment = 0;
    eDs->eVertex = eVs;
    eDs->ePrevious = 0;
    eDs->CreateCanvasVTX();
    eVs->V()->rmsDistAngle();
    sprintf(text,"Orig    %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
    eVs->ID(), eVs->N(), eVs->VX(), eVs->VY(), eVs->VZ(), eVs->V()->dist(),
    eVs->V()->chi2()/eVs->V()->ndf(), eVs->V()->prob());
    eDs->DrawOldVTX(text);
    eDs->DrawVTXTracks("Original", eDs->eVertex);
    eDs->DrawEnv();
    if (eDs->eArrV) eDs->eIndVert = (eDs->eArrV)->IndexOf(eVs);
    else            eDs->eIndVert = -1;
    (eDs->eCreatedTracks).Clear();
    eDs->Draw();
  }
}
//_____________________________________________________________________________
void EdbTrackG::SetAsWorkingVertex()
{
  printf("1\n");
  EdbTrackP *tr = 0;
  EdbDisplay *eDs = 0;
  EdbVTA *vta = 0;
  EdbVertex  *eVs = 0, *old = 0;
  int zpos = 1;
  eDs = eD;
  if (eDs && eTr)
  {
    if (eDs->eWait_Answer) return;
    if (eDs->eWorking)
    {
      eDs->DialogModifiedVTX();
      return;
    }
    if (eDs->eVertex)
    {
      if (eDs->eVertex->ID()<0 && eDs->eVertex->N() >= 2)
      {
        eDs->DialogModifiedVTX();
        return;
      }
    }
    if (eDs->eVertex)
    {
      eDs->CancelModifiedVTX();
    }
    if (eDs->eSegPM)
    { 
      eDs->ClearSegmentEnv();
    }
    if (GetMarkerColor() == kRed) zpos = 0;
    if ((old = eTr->VertexS()) && (zpos == 1))
    {
    
      printf("Track alredy connected to a vertex by this edge!\n");
      fflush(stdout);
      return;
    }
    if ((old = eTr->VertexE()) && (zpos == 0))
    {
    
      printf("Track alredy connected to a vertex by this edge!\n");
      fflush(stdout);
      return;
    }
    eVs = new EdbVertex();
    if ((vta = (eD->eVerRec)->AddTrack(*eVs, eTr, zpos)))
    {
      eTr->AddVTA(vta);
    }
    else
    {
      delete eVs;
      printf("Can't connect track to a vertex!\n");
      fflush(stdout);
      return;
    }
    double dz = 0.;
    EdbSegP *seg = 0;
    if( zpos == 0 )
    {
      seg = eTr->TrackZmax();
      dz = TMath::Abs(seg->DZ());
    }
    else
    {
      seg = eTr->TrackZmin();
    }
    if(!seg)
    {
      eVs->SetXYZ(eTr->X(), eTr->Y(), eTr->Z());
    }
    else
    {
      eVs->SetXYZ(seg->X()+seg->TX()*dz, seg->Y()+seg->TY()*dz, seg->Z()+dz);
    }
    printf("2\n");
    eVs->SetID(-1);
    SetMarkerColor(kGreen);
    eDs->eSegment = 0;
    eDs->eVertex = eVs;
    eDs->eWorking = 0;
    eDs->ePrevious = 0;
    eDs->eIndVert = -1;
    (eDs->eCreatedTracks).Clear();
    if (!(eDs->eArrV)) 
    {
      eDs->eArrV = new TObjArray(20);
      eDs->eArrV->Add((TObject *)eVs);
      if (!(eDs->eArrTr))  eDs->eArrTr = new TObjArray(20);
      for (int i=0; i<eVs->N(); i++)
      {
        tr = eVs->GetTrack(i);
        if(!(eDs->eArrTr->FindObject(tr))) eDs->eArrTr->Add(tr);
      }
      eDs->Draw();
    }
    else
    {
      if (!((eDs->eArrV)->FindObject(eVs)))
      {
        eDs->eArrV->Add(eVs);
        if (!(eDs->eArrTr))  eDs->eArrTr = new TObjArray(20);
        for (int i=0; i<eVs->N(); i++)
        {
          tr = eVs->GetTrack(i);
          if(!(eDs->eArrTr->FindObject(tr))) eDs->eArrTr->Add(tr);
        }
        eDs->Draw();
      }
    }
    //eD->DrawEnv();
    eD->DrawAcc();
    eD->DrawCan();
    eD->DrawUnd();
    eD->Draw();
    printf("3\n");
  }
}

//_____________________________________________________________________________
void EdbVertexG::DeleteVertex()
{
  EdbDisplay *eDs = 0;
  EdbVertex  *eVs = 0;
  eDs = eD;
  eVs = eV;
  if (eDs && eVs)
  {
    if (eDs->eWait_Answer) return;
    if (eDs->eWorking == eVs)
    {
	eDs->DialogModifiedVTX();
	return;
    }
    if (eDs->eVertex == eVs)
    {
    	eDs->CancelModifiedVTX();
    }
    if (eDs->eArrV) 
    {
	eDs->eArrV->Remove((TObject *)eVs);
	eDs->eArrV->Compress();
	eDs->Draw();
    }
    else
    {
	if ((eDs->eArrV)->FindObject(eVs))
	{
	    eDs->eArrV->Remove(eVs);
	    eDs->eArrV->Compress();
	    eDs->eVertex = 0;
	    eDs->Draw();
	}
    }
    delete eVs;
    eDs->Draw();
  }
}

//_____________________________________________________________________________
void EdbVertexG::RemoveKink()
{
  EdbDisplay *eDs = 0;
  EdbVertex  *eVs = 0;
  eDs = eD;
  eVs = eV;
  if (eDs && eVs)
  {
    if (eDs->eWait_Answer) return;
    if (eVs->N() != 2)
    {
	printf("Wrong vertex type - not a kink-like!\n");
	fflush(stdout);
	return;
    }
    EdbTrackP *etr = 0;
    for (int it=0; it<eVs->N(); it++)
    {
	etr = eVs->GetTrack(it);
	if (etr->VTAS() && etr->VTAE())
	{
	    if (etr->VertexS() == eVs && etr->VertexE()->Flag() >= 0)
	    {
		printf("Vertex track (at zmax) belong to another vertex (ID=%d) too - impossible delete it!\n", etr->VertexE()->ID());
		fflush(stdout);
		return;
	    }
	    if (etr->VertexE() == eVs && etr->VertexS()->Flag() >= 0)
	    {
		printf("Vertex track (at zmin) belong to another vertex (ID=%d) too - impossible delete it!\n", etr->VertexS()->ID());
		fflush(stdout);
		return;
	    }
	}
    }
    if (eDs->eWorking == eVs)
    {
	eDs->DialogModifiedVTX();
	return;
    }
    if (eDs->eVertex == eVs)
    {
    	eDs->CancelModifiedVTX();
    }
    if (eDs->eArrV) 
    {
	eDs->eArrV->Remove((TObject *)eVs);
	eDs->eArrV->Compress();
	eDs->Draw();
    }
    else
    {
	if ((eDs->eArrV)->FindObject(eVs))
	{
	    eDs->eArrV->Remove(eVs);
	    eDs->eArrV->Compress();
	    eDs->eVertex = 0;
	    eDs->Draw();
	}
    }
    EdbSegP *seg = 0;
    if (!eD->eTrack) eD->eTrack = new EdbTrackP();
    for (int it=0; it<eVs->N(); it++)
    {
	etr = eVs->GetTrack(it);
    	for(int is=0; is<etr->N(); is++) {
	    seg = etr->GetSegment(is);
	    (eD->eTrack)->AddSegment((EdbSegP *)seg);
	}
	etr->SetFlag(-10);
	etr->SetSegmentsTrack(-2-(etr->ID()+1));
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(etr))
	{
	    eD->eArrTr->Remove(etr);
	    eD->eArrTr->Compress();
	    if (!(eD->eArrSegP)) eD->eArrSegP = new TObjArray();
    	    for(int is=0; is<etr->N(); is++) {
		seg = etr->GetSegment(is);
		if(eD->eArrSegP) if(!((eD->eArrSegP)->FindObject(seg))) eD->eArrSegP->Add(seg);
	    }
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(etr))
	{
	    eD->eArrTrSave->Remove(etr);
	    eD->eArrTrSave->Compress();
	    if (!(eD->eArrSegPSave)) eD->eArrSegPSave = new TObjArray();
    	    for(int is=0; is<etr->N(); is++) {
		seg = etr->GetSegment(is);
		if(eD->eArrSegPSave) if(!((eD->eArrSegPSave)->FindObject(seg))) eD->eArrSegPSave->Add(seg);
	    }
	}
    }
    if(eD->eArrTr) eD->eArrTr->Compress();
    if(eD->eArrTrSave) eD->eArrTrSave->Compress();
    eVs->SetFlag(-11);
    eDs->Draw();
    float mass = etr->M();
    float momentum = etr->P();
    if (eD->eP > 0.) momentum = eD->eP;
    if (eD->eM > 0.) mass = eD->eM;
    (eD->eTrack)->SetM(mass);
    (eD->eTrack)->SetP(momentum);
    float X0 = 0.;
    if (eD->eVerRec) if ((eD->eVerRec)->ePVR) X0 = (((eD->eVerRec)->ePVR)->GetScanCond())->RadX0();
    (eD->eTrack)->FitTrackKFS(true, X0, 0);
    if (!eD->eArrTr) eD->eArrTr = new TObjArray();
    (eD->eArrTr)->Add(eD->eTrack);
    if (eD->eArrTrSave) (eD->eArrTrSave)->Add(eD->eTrack);
    eD->Draw();
  }
}
//_____________________________________________________________________________
void EdbTrackG::UndoNewTrack()
{
    if (eD->eTrack && eD->eTrack == eTr)
    {
	if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
	if (!eD->eVerRec) {printf("Error: EdbDisplay:UndoNewTrack: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
	TObjArray *etr = 0;
	if (eD->eVerRec) etr = eD->eVerRec->eEdbTracks;

    	for(int is=0; is<eD->eTrack->N(); is++) {
	    EdbSegP *seg = eD->eTrack->GetSegment(is);
	    if (seg->Track() == -2) seg->SetTrack(-1);
	    else if (seg->Track() < -2)
	    {
		int ind = -(seg->Track()+2) - 1;
		if (!etr)
		{
		    printf("No information about tracks array!\n");
		    fflush(stdout);
		    return;
		}
		int trind = etr->GetEntries();
		if (ind >= trind)
		{
		    printf("Wrong track index in segment!\n");
		    fflush(stdout);
		    return;
		}
		EdbTrackP *tr = (EdbTrackP *)(etr->At(ind));
		if (!tr)
		{
		    printf("Wrong track address in array!\n");
		    fflush(stdout);
		    return;
		}
		tr->SetFlag(0);
		tr->SetSegmentsTrack();
		if(eD->eArrTr) if(!((eD->eArrTr)->FindObject(tr)))
		{
		    eD->eArrTr->Add(tr);
		}
		if(eD->eArrTrSave) if(!((eD->eArrTrSave)->FindObject(tr)))
		{
		    eD->eArrTrSave->Add(tr);
		}
	    }
	}
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(eTr))
	{
	    eD->eArrTr->Remove(eTr);
	    eD->eArrTr->Compress();
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(eTr))
	{
	    eD->eArrTrSave->Remove(eTr);
	    eD->eArrTr->Compress();
	}
	eD->Draw();
	delete eD->eTrack;
	eD->eTrack = 0;
    }
    else
    {
	printf("No intermediate track or this one fixed already!\n");
	fflush(stdout);
	return;
    }
}
//_____________________________________________________________________________
void EdbTrackG::UndoSplit()
{
    if (eD->eTrack1 && eD->eTrack2 && ((eD->eTrack1 == eTr) || (eD->eTrack2 == eTr)))
    {
	if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
	if (!eD->eVerRec) {printf("Error: EdbDisplay:UndoSplit: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
	TObjArray *etr = 0;
	if (eD->eVerRec) etr = eD->eVerRec->eEdbTracks;

    	for(int is=0; is<eD->eTrack1->N(); is++) {
	    EdbSegP *seg = eD->eTrack1->GetSegment(is);
	    if (seg->Track() == -2) seg->SetTrack(-1);
	    else if (seg->Track() < -2)
	    {
		int ind = -(seg->Track()+2) - 1;
		if (!etr)
		{
		    printf("No information about tracks array!\n");
		    fflush(stdout);
		    return;
		}
		int trind = etr->GetEntries();
		if (ind >= trind)
		{
		    printf("Wrong track index in segment!\n");
		    fflush(stdout);
		    return;
		}
		EdbTrackP *tr = (EdbTrackP *)(etr->At(ind));
		if (!tr)
		{
		    printf("Wrong track address in array!\n");
		    fflush(stdout);
		    return;
		}
		tr->SetFlag(0);
		tr->SetSegmentsTrack();
		if(eD->eArrTr) if(!((eD->eArrTr)->FindObject(tr)))
		{
		    eD->eArrTr->Add(tr);
		}
		if(eD->eArrTrSave) if(!((eD->eArrTrSave)->FindObject(tr)))
		{
		    eD->eArrTrSave->Add(tr);
		}
	    }
	}
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(eD->eTrack1))
	{
	    eD->eArrTr->Remove(eD->eTrack1);
	    eD->eArrTr->Compress();
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(eD->eTrack1))
	{
	    eD->eArrTrSave->Remove(eD->eTrack1);
	    eD->eArrTr->Compress();
	}
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(eD->eTrack2))
	{
	    eD->eArrTr->Remove(eD->eTrack2);
	    eD->eArrTr->Compress();
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(eD->eTrack2))
	{
	    eD->eArrTrSave->Remove(eD->eTrack2);
	    eD->eArrTr->Compress();
	}
	eD->Draw();
	delete eD->eTrack1;
	eD->eTrack1 = 0;
	delete eD->eTrack2;
	eD->eTrack2 = 0;
    }
    else
    {
	printf("No parts of splitted track or these ones fixed already!\n");
	fflush(stdout);
	return;
    }
}
//_____________________________________________________________________________
void EdbTrackG::FixNewTrack()
{
  if (eTr && eD)
  {
	if (eD->eTrack)
	{
		if (eTr == eD->eTrack)
		{
		    TObjArray *etr = 0;
		    if (eD->eVerRec) etr = (eD->eVerRec)->eEdbTracks;
		    int trind = 0;
		    if (etr) trind = etr->GetEntries();
		    (eD->eTrack)->SetID(trind);
		    if (etr) etr->Add(eD->eTrack);
		    (eD->eTrack)->SetSegmentsTrack();
		    eD->eTrack = 0;
		    eD->Draw();
		}
	}
	if (eD->eTrack1)
	{
		if (eTr == eD->eTrack1)
		{
		    TObjArray *etr = 0;
		    if (eD->eVerRec) etr = (eD->eVerRec)->eEdbTracks;
		    int trind = 0;
		    if (etr) trind = etr->GetEntries();
		    (eD->eTrack1)->SetID(trind);
		    if (etr) etr->Add(eD->eTrack1);
		    (eD->eTrack1)->SetSegmentsTrack();
		    eD->eTrack1 = 0;
		    eD->Draw();
		}
	}

	if (eD->eTrack2)
	{
		if (eTr == eD->eTrack2)
		{
		    TObjArray *etr = 0;
		    if (eD->eVerRec) etr = (eD->eVerRec)->eEdbTracks;
		    int trind = 0;
		    if (etr) trind = etr->GetEntries();
		    (eD->eTrack2)->SetID(trind);
		    if (etr) etr->Add(eD->eTrack2);
		    (eD->eTrack2)->SetSegmentsTrack();
		    eD->eTrack2 = 0;
		    eD->Draw();
		}
	}
  }
}
//_____________________________________________________________________________
void EdbTrackG::UndoRemoveKink()
{
    if (eD && eTr)
    {
	if (!(eD->eTrack))
	{
	    printf("No joined track or this one fixed already!\n");
	    fflush(stdout);
	    return;
	}
	if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
	if (!eD->eVerRec) {printf("Error: EdbDisplay:UndoRemoveKink: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
	TObjArray *etr = 0;
	if (eD->eVerRec) etr = eD->eVerRec->eEdbTracks;

    	for(int is=0; is<eD->eTrack->N(); is++) {
	    EdbSegP *seg = eD->eTrack->GetSegment(is);
	    if (seg->Track() == -2) seg->SetTrack(-1);
	    else if (seg->Track() < -2)
	    {
		int ind = -(seg->Track()+2) - 1;
		if (!etr)
		{
		    printf("No information about tracks array!\n");
		    fflush(stdout);
		    return;
		}
		int trind = etr->GetEntries();
		if (ind >= trind)
		{
		    printf("Wrong track index in segment!\n");
		    fflush(stdout);
		    return;
		}
		EdbTrackP *tr = (EdbTrackP *)(etr->At(ind));
		if (!tr)
		{
		    printf("Wrong track address in array!\n");
		    fflush(stdout);
		    return;
		}
		tr->SetFlag(0);
		tr->SetSegmentsTrack();
		if(eD->eArrTr) if(!((eD->eArrTr)->FindObject(tr)))
		{
		    eD->eArrTr->Add(tr);
		}
		if(eD->eArrTrSave) if(!((eD->eArrTrSave)->FindObject(tr)))
		{
		    eD->eArrTrSave->Add(tr);
		}
		EdbVertex *v = 0;
		if (tr->VTAS())
		{
		    v = tr->VertexS();
		    if (v)
		    {
			if (v->Flag() == -11)
			{
			    v->SetFlag(0);
			    if(eD->eArrV) if(!((eD->eArrV)->FindObject(v)))
			    {
				eD->eArrV->Add(v);
			    }
			    if(eD->eArrVSave) if(!((eD->eArrVSave)->FindObject(v)))
			    {
				eD->eArrVSave->Add(v);
			    }
			}
		    }
		}
		if (tr->VTAE())
		{
		    v = tr->VertexE();
		    if (v)
		    {
			if (v->Flag() == -11)
			{
			    v->SetFlag(0);
			    if(eD->eArrV) if(!((eD->eArrV)->FindObject(v)))
			    {
				eD->eArrV->Add(v);
			    }
			    if(eD->eArrVSave) if(!((eD->eArrVSave)->FindObject(v)))
			    {
				eD->eArrVSave->Add(v);
			    }
			}
		    }
		}
	    }
	}
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(eD->eTrack))
	{
	    eD->eArrTr->Remove(eD->eTrack);
	    eD->eArrTr->Compress();
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(eD->eTrack))
	{
	    eD->eArrTrSave->Remove(eD->eTrack);
	    eD->eArrTr->Compress();
	}
	eD->Draw();
	delete eD->eTrack;
	eD->eTrack = 0;
    }
    else
    {
	printf("No parts of splitted track or these ones fixed already!\n");
	fflush(stdout);
	return;
    }
}
//_____________________________________________________________________________
void EdbTrackG::DeleteTrack()
{
  if (eTr && eD)
  {
	if (eTr->VTAS() || eTr->VTAE())
	{
	    printf("Track belong to a vertex - impossible delete it!\n");
	    fflush(stdout);
	    return;
	}
	eTr->SetFlag(-10);
	eTr->SetSegmentsTrack(-1);
	EdbSegP *seg = 0;
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(eTr))
	{
	    eD->eArrTr->Remove(eTr);
	    eD->eArrTr->Compress();
	    if (!(eD->eArrSegP)) eD->eArrSegP = new TObjArray();
    	    for(int is=0; is<eTr->N(); is++) {
		seg = eTr->GetSegment(is);
		if(eD->eArrSegP) if(!((eD->eArrSegP)->FindObject(seg))) eD->eArrSegP->Add(seg);
	    }
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(eTr))
	{
	    eD->eArrTrSave->Remove(eTr);
	    eD->eArrTr->Compress();
	    if (!(eD->eArrSegPSave)) eD->eArrSegPSave = new TObjArray();
    	    for(int is=0; is<eTr->N(); is++) {
		seg = eTr->GetSegment(is);
		if(eD->eArrSegPSave) if(!((eD->eArrSegPSave)->FindObject(seg))) eD->eArrSegPSave->Add(seg);
	    }
	}
	eD->Draw();
  }
}
//_____________________________________________________________________________
void EdbTrackG::AddToNewTrack()
{
  if (eTr && eD)
  {
	if (eTr->VTAS() || eTr->VTAE())
	{
	    printf("Track belong to a vertex - impossible destroy it!\n");
	    fflush(stdout);
	    return;
	}
	if (eD->eTrack)
	{
	    if ((eD->eTrack)->NF())
	    {
		delete eD->eTrack;
		eD->eTrack = 0;
	    }
	}
	EdbSegP *seg = 0;
	if (!eD->eTrack) eD->eTrack = new EdbTrackP();
    	for(int is=0; is<eTr->N(); is++) {
	    seg = eTr->GetSegment(is);
	    (eD->eTrack)->AddSegment((EdbSegP *)seg);
	}
	eTr->SetFlag(-10);
	eTr->SetSegmentsTrack(-2-(eTr->ID()+1));
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(eTr))
	{
	    eD->eArrTr->Remove(eTr);
	    eD->eArrTr->Compress();
	    if (!(eD->eArrSegP)) eD->eArrSegP = new TObjArray();
    	    for(int is=0; is<eTr->N(); is++) {
		seg = eTr->GetSegment(is);
		if(eD->eArrSegP) if(!((eD->eArrSegP)->FindObject(seg))) eD->eArrSegP->Add(seg);
	    }
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(eTr))
	{
	    eD->eArrTrSave->Remove(eTr);
	    eD->eArrTrSave->Compress();
	    if (!(eD->eArrSegPSave)) eD->eArrSegPSave = new TObjArray();
    	    for(int is=0; is<eTr->N(); is++) {
		seg = eTr->GetSegment(is);
		if(eD->eArrSegPSave) if(!((eD->eArrSegPSave)->FindObject(seg))) eD->eArrSegPSave->Add(seg);
	    }
	}
	eD->Draw();
  }
}
//_____________________________________________________________________________
void EdbTrackG::AddToNewTrackAndFit()
{
  if (eTr && eD)
  {
	if (eTr->VTAS() || eTr->VTAE())
	{
	    printf("Track belong to a vertex - impossible destroy it!\n");
	    fflush(stdout);
	    return;
	}
	EdbSegP *seg = 0;
	if (!eD->eTrack) eD->eTrack = new EdbTrackP();
    	for(int is=0; is<eTr->N(); is++) {
	    seg = eTr->GetSegment(is);
	    (eD->eTrack)->AddSegment((EdbSegP *)seg);
	}
	eTr->SetFlag(-10);
	eTr->SetSegmentsTrack(-2-(eTr->ID()+1));
	if(eD->eArrTr) if((eD->eArrTr)->FindObject(eTr))
	{
	    eD->eArrTr->Remove(eTr);
	    eD->eArrTr->Compress();
	    if (!(eD->eArrSegP)) eD->eArrSegP = new TObjArray();
    	    for(int is=0; is<eTr->N(); is++) {
		seg = eTr->GetSegment(is);
		if(eD->eArrSegP) if(!((eD->eArrSegP)->FindObject(seg))) eD->eArrSegP->Add(seg);
	    }
	}
	if(eD->eArrTrSave) if((eD->eArrTrSave)->FindObject(eTr))
	{
	    eD->eArrTrSave->Remove(eTr);
	    eD->eArrTrSave->Compress();
	    if (!(eD->eArrSegPSave)) eD->eArrSegPSave = new TObjArray();
    	    for(int is=0; is<eTr->N(); is++) {
		seg = eTr->GetSegment(is);
		if(eD->eArrSegPSave) if(!((eD->eArrSegPSave)->FindObject(seg))) eD->eArrSegPSave->Add(seg);
	    }
	}
	eD->Draw();
	float mass = eTr->M();
	float momentum = eTr->P();
	if (eD->eP > 0.) momentum = eD->eP;
	if (eD->eM > 0.) mass = eD->eM;
	eD->eTrack->SetM(mass);
	eD->eTrack->SetP(momentum);
	float X0 = 0.;
	if (eD->eVerRec) if ((eD->eVerRec)->ePVR) X0 = (((eD->eVerRec)->ePVR)->GetScanCond())->RadX0();
	eD->eTrack->FitTrackKFS(true, X0, 0);
	if (!eD->eArrTr) eD->eArrTr = new TObjArray();
	(eD->eArrTr)->Add(eD->eTrack);
	if (eD->eArrTrSave) (eD->eArrTrSave)->Add(eD->eTrack);
	eD->Draw();
  }
}
//_____________________________________________________________________________
void EdbTrackG::RemoveTrackFromVertex()
{
  char text[512];
  if (eTr && eD)
  {
    if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
    if (!eD->eVerRec) {printf("Error: EdbDisplay:RemoveTrackFromVertex: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
    if (eD->eWait_Answer) return;
    if (!(eD->eVertex))
    {
    
	    printf("No working vertex selected!\n");
	    fflush(stdout);
	    return;
    }
    EdbVTA *vta = 0;
    EdbVertex *ePreviousSaved = eD->ePrevious;
    int n = 0;
    int ntr = 0;
    if (eD->eWorking == 0)
    {
	ntr = eD->eVertex->N();
	if (ntr < 3)
	{
    
	    printf("Working vertex has 2 prongs only!\n");
	    fflush(stdout);
	    return;
	}
	eD->eWorking = new EdbVertex();
	int i = 0;
	for(i=0; i<ntr; i++)
	{
	    if (eD->eVertex->GetTrack(i) == eTr)
	    {
		eTr->ClearVTA(eD->eVertex->GetVTa(i));
		continue;
	    }
	    if ((vta = (eD->eVerRec)->AddTrack( *(eD->eWorking), (eD->eVertex)->GetTrack(i), (eD->eVertex)->Zpos(i))))
	    {
		((eD->eVertex)->GetTrack(i))->AddVTA(vta);
		n++;
	    }
	}
    }
    else
    {
	ntr = eD->eWorking->N();
	if (ntr < 3)
	{
    
	    printf("Working vertex has 2 prongs only!\n");
	    fflush(stdout);
	    return;
	}
	eD->ePrevious = eD->eWorking;
	eD->eWorking = new EdbVertex();
	int i = 0;
	for(i=0; i<ntr; i++)
	{
	    if (eD->ePrevious->GetTrack(i) == eTr)
	    {
		eTr->ClearVTA(eD->ePrevious->GetVTa(i));
		continue;
	    }
	    if ((vta = (eD->eVerRec)->AddTrack(*(eD->eWorking),(eD->ePrevious)->GetTrack(i), (eD->ePrevious)->Zpos(i))))
	    {
		((eD->ePrevious)->GetTrack(i))->AddVTA(vta);
		n++;
	    }
	}
    }
    if ((n < 2)||(n == ntr))
    {
	delete eD->eWorking;
	if (eD->ePrevious)
	{
	    eD->eWorking = eD->ePrevious;
	    (eD->eWorking)->ResetTracks();
	    eD->ePrevious = ePreviousSaved;
	}
	else
	{
	    eD->eWorking = 0;
	    (eD->eVertex)->ResetTracks();
	}
	printf("Can't create working copy of the vertex!\n");
	fflush(stdout);
	return;
    }

    if ((eD->eVerRec)->MakeV(*(eD->eWorking)))
    {
	EdbVertex *eW = eD->eWorking;
	eW->ResetTracks();
	eW->SetID(eD->eVertex->ID());
	eW->V()->rmsDistAngle();
	sprintf(text,"New     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	eD->DrawOldBut("Original");
	if (eD->ePrevious)
	{
	    eD->DrawNewVTX(text);
	    eD->DrawNewBut("Modified");
	    eW = eD->ePrevious;
	    sprintf(text,"Pre     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	    eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	    eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	    eD->DrawPreVTX(text);
	    eD->DrawPreBut("Previous");
	}
	else
	{
	    eD->DrawPreVTX(text);
	    eD->DrawPreBut("Modified");
	    if (eD->eVertex->ID() >= 0)
	    {
		eD->DrawAcc();
		eD->DrawCan();
		eD->DrawUnd();
	    }
	}
//	eD->DrawOldBut("Original");
	eD->DrawVTXTracks("Modified", eD->eWorking);
	if (eD->eArrV && (eD->eIndVert >= 0)) 
	{
//	    eD->eArrV->RemoveAt(eD->eIndVert);
	    eD->eArrV->AddAt(eD->eWorking, eD->eIndVert);
	}	
//	if ((eD->eCreatedTracks).FindObject(eTr))
//	{
//	    (eD->eCreatedTracks).Remove(eTr);
//	    if ( eD->eArrTr)
//		if ((eD->eArrTr)->FindObject(eTr)) (eD->eArrTr)->Remove(eTr);
//	    if ( eD->eArrTrSave)
//		if ((eD->eArrTrSave)->FindObject(eTr)) (eD->eArrTrSave)->Remove(eTr);
//	    delete eTr;
//	}
	eD->Draw();
	if (ePreviousSaved) delete ePreviousSaved;
	ePreviousSaved = 0;
    }
    else
    {
	delete eD->eWorking;
	if (eD->ePrevious)
	{
	    eD->eWorking = eD->ePrevious;
	    (eD->eWorking)->ResetTracks();
	    eD->ePrevious = ePreviousSaved;
	}
	else
	{
	    eD->eWorking = 0;
	    (eD->eVertex)->ResetTracks();
	}
	printf("Can't create working copy of the vertex!\n");
	fflush(stdout);
    }
  }
}


//_____________________________________________________________________________
void EdbDisplay::RemoveTrackFromTable(int ivt)
{
  char text[512];
  EdbTrackP *etr = 0;
  if (eWait_Answer) return;
  if (!(eVertex)) return;
  Log(3,"EdbDisplay::RemoveTrackFromTable","%d tracks in vertex before",eVertex->N());
  
  EdbVTA *vta = 0;
  EdbVertex *ePreviousSaved = ePrevious;
  int n = 0;
  int ntr = 0;
  if (eWorking == 0)
    {
      ntr = eVertex->N();
      if (ntr < 3)
	{
	  
	  printf("Working vertex has 2 prongs only!\n");
	  fflush(stdout);
	  return;
	}
      if (eVerRec) if (eVerRec->IsA() != EdbVertexRec::Class()) eVerRec = 0;
      if (!eVerRec) {
	printf("Error: EdbDisplay:RemoveTrackFromTable: EdbVertexRec not defined, use SetVerRec(...)\n"); 
	fflush(stdout); return;
      }
      eWorking = new EdbVertex();
      int i = 0;
      etr = eVertex->GetTrack(ivt);
      for(i=0; i<ntr; i++)
	{
	  if (i == ivt)
	    {
	      etr->ClearVTA(eVertex->GetVTa(i));
	      continue;
	    }
	  if ((vta = (eVerRec)->EdbVertexRec::AddTrack( *(eWorking), (eVertex)->GetTrack(i), (eVertex)->Zpos(i))))
	    {
	      ((eVertex)->GetTrack(i))->AddVTA(vta);
	      n++;
	    }
	}
    }
  else
    {
      ntr = eWorking->N();
      if (ntr < 3)
	{
	  
	  printf("Working vertex has 2 prongs only!\n");
	  fflush(stdout);
	  return;
	}
      if (eVerRec) if (eVerRec->IsA() != EdbVertexRec::Class()) eVerRec = 0;
      if (!eVerRec) {printf("Error: EdbDisplay:RemoveTrackFromTable: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
      etr = eWorking->GetTrack(ivt);
      ePrevious = eWorking;
      eWorking = new EdbVertex();
      int i = 0;
      for(i=0; i<ntr; i++)
	{
	  if (i == ivt)
	    {
	      etr->ClearVTA(ePrevious->GetVTa(i));
	      continue;
	    }
	  if ((vta = (eVerRec)->EdbVertexRec::AddTrack(*(eWorking),(ePrevious)->GetTrack(i), (ePrevious)->Zpos(i))))
	    {
	      ((ePrevious)->GetTrack(i))->AddVTA(vta);
	      n++;
	    }
	}
    }
  if ((n < 2)||(n == ntr))
    {
      delete eWorking;
      if (ePrevious)
	{
	  eWorking = ePrevious;
	  (eWorking)->ResetTracks();
	  ePrevious = ePreviousSaved;
	}
      else
	{
	  eWorking = 0;
	  (eVertex)->ResetTracks();
	}
      printf("Can't create working copy of the vertex!\n");
      fflush(stdout);
      return;
    }

  if ((eVerRec)->MakeV(*(eWorking)))
    {
      EdbVertex *eW = eWorking;
      eW->ResetTracks();
      eW->SetID(eVertex->ID());
      eW->V()->rmsDistAngle();
      sprintf(text,"New     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	      eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	      eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
      DrawOldBut("Original");
      if (ePrevious)
	{
	  DrawNewVTX(text);
	  DrawNewBut("Modified");
	  eW = ePrevious;
	  sprintf(text,"Pre     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
		  eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
		  eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	  DrawPreVTX(text);
	  DrawPreBut("Previous");
	}
      else
	{
	  DrawPreVTX(text);
	  DrawPreBut("Modified");
	  if (eVertex->ID() >= 0)
	    {
	      DrawAcc();
	      DrawCan();
	      DrawUnd();
	    }
	}
      //	DrawOldBut("Original");
      TButton *rm = fRemBut[ivt];
      fRemBut[ivt] = 0;
      DrawVTXTracks("Modified", eWorking);
      if (eArrV && (eIndVert >= 0)) 
	{
	  //	    eArrV->RemoveAt(eIndVert);
	  eArrV->AddAt(eWorking, eIndVert);
	}	
      //	if (eCreatedTracks.FindObject(etr))
      //	{
      //	    eCreatedTracks.Remove(etr);
      //	    if ( eArrTr)
      //		if (eArrTr->FindObject(etr)) eArrTr->Remove(etr);
      //	    if ( eArrTrSave)
      //		if (eArrTrSave->FindObject(etr)) eArrTrSave->Remove(etr);
      //	    delete etr;
      //	} 
      fCanvas->cd();
      Draw();
      fPad->Update(); 
      if (ePreviousSaved) delete ePreviousSaved;
      ePreviousSaved = 0;
      fCanvasVTX->cd();
      if (rm) delete rm;
    }
  else
    {
      delete eWorking;
      if (ePrevious)
	{
	  eWorking = ePrevious;
	  (eWorking)->ResetTracks();
	  ePrevious = ePreviousSaved;
	}
      else
	{
	  eWorking = 0;
	  (eVertex)->ResetTracks();
	}
      printf("Can't create new vertex after track removing!\n");
      fflush(stdout);
    }
}

//_____________________________________________________________________________
void EdbTrackG::AddTrackToVertex()
{
  char text[512];
  int zpos = 1, zpos2 = 1;
  EdbVTA *vta = 0;
  EdbVertex *old = 0, *eVs = 0;
  EdbTrackP *eTr2 = 0;
  if (eTr && eD)
  {
    if (eD->eWait_Answer) return;
    if (!(eD->eVertex))
    {
    
	    printf("No working vertex selected!\n");
	    fflush(stdout);
	    return;
    }
    if (GetMarkerColor() == kRed) zpos = 0;
    if ((old = eTr->VertexS()) && (zpos == 1))
    {
    
	    printf("Track alredy connected to a vertex by this edge!\n");
	    fflush(stdout);
	    return;
    }
//    {
//	if (old != eD->eVertex && old != eD->ePrevious && old != eD->eWorking) return;
//    } 
    if ((old = eTr->VertexE()) && (zpos == 0))
    {
    
	    printf("Track alredy connected to a vertex by this edge!\n");
	    fflush(stdout);
	    return;
    }
//    {
//	if (old != eD->eVertex && old != eD->ePrevious && old != eD->eWorking) return;
//    } 
    if (eD->eVerRec) if (eD->eVerRec->IsA() != EdbVertexRec::Class()) eD->eVerRec = 0;
    if (!eD->eVerRec) {printf("Error: EdbDisplay:AddTrackToVertex: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
    double ImpMaxSave = 0.;
    double ProbMinSave = 0.;
    if (eD->eVertex->N() == 1 && eD->eVertex->ID() < 0)
    {
	ImpMaxSave = (eD->eVerRec)->eImpMax;
	(eD->eVerRec)->eImpMax = eD->eTImpMax;
	ProbMinSave = (eD->eVerRec)->eProbMin;
	(eD->eVerRec)->eProbMin = eD->eTProbMin;
	eTr2 = (eD->eVertex)->GetTrack(0);
	zpos2 = (eD->eVertex)->Zpos(0);
	if((eVs = eD->eVerRec->ProbVertex2(eTr2, eTr, zpos2, zpos)))
	{
	    if (eD->eArrV) 
	    {
		if (eD->eArrV->FindObject(eD->eVertex))
		{
		    eD->eArrV->Remove(eD->eVertex);
		    eD->eArrV->Add(eVs);
		}
	    }
	    delete eD->eVertex;
	    eVs->SetID(-1);
	    eD->eVertex = eVs;
	    eD->CreateCanvasVTX();
	    eVs->V()->rmsDistAngle();
	    sprintf(text,"Creat   %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	    eVs->ID(), eVs->N(), eVs->VX(), eVs->VY(), eVs->VZ(), eVs->V()->dist(),
	    eVs->V()->chi2()/eVs->V()->ndf(), eVs->V()->prob());
	    eD->DrawOldVTX(text);
	    eD->DrawVTXTracks("Created", eD->eVertex);
	    eD->DrawEnv();
	    eD->Draw();
	}
	else
	{
	    printf("EdbTrackG::AddTrackToVertex: New vertex not created! May be Prob < ProbMin=%f. Change ProbMin with 'TrackParams' button!\n",eD->eTProbMin);
	}
	(eD->eVerRec)->eImpMax = ImpMaxSave;
	(eD->eVerRec)->eProbMin = ProbMinSave;
	return;
    }
    EdbVertex *ePreviousSaved = eD->ePrevious;
    if (eD->eWorking == 0)
    {
	eD->eWorking = new EdbVertex();
	int ntr = eD->eVertex->N();
	int i = 0, n = 0;
	for(i=0; i<ntr; i++)
	{
	    if ((vta = (eD->eVerRec)->AddTrack(*(eD->eWorking), (eD->eVertex)->GetTrack(i), (eD->eVertex)->Zpos(i))))
	    {
		(eD->eVertex)->GetTrack(i)->AddVTA(vta);
		n++;
	    }
	}
	if (n < 2)
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}

	if (!((eD->eVerRec)->MakeV(*(eD->eWorking))))
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}
    }
    else
    {
//	if (ePrevious) delete ePrevious;
	eD->ePrevious = eD->eWorking;
	eD->eWorking = new EdbVertex();
	int ntr = eD->ePrevious->N();
	int i = 0, n = 0;
	for(i=0; i<ntr; i++)
	{
	    if ((vta = (eD->eVerRec)->AddTrack(*(eD->eWorking),(eD->ePrevious)->GetTrack(i), (eD->ePrevious)->Zpos(i))))
	    {
		((eD->ePrevious)->GetTrack(i))->AddVTA(vta);
		n++;
	    }
	}
	if (n < 2)
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
		eD->ePrevious = ePreviousSaved;
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}

	if (!((eD->eVerRec)->MakeV(*(eD->eWorking))))
	{
	    delete eD->eWorking;
	    if (eD->ePrevious)
	    {
		eD->eWorking = eD->ePrevious;
		(eD->eWorking)->ResetTracks();
		eD->ePrevious = ePreviousSaved;
	    }
	    else
	    {
		eD->eWorking = 0;
		(eD->eVertex)->ResetTracks();
	    }
	    printf("Can't create working copy of the vertex!\n");
	    fflush(stdout);
	    return;
	}
    }
    ImpMaxSave = (eD->eVerRec)->eImpMax;
    (eD->eVerRec)->eImpMax = eD->eTImpMax;
    ProbMinSave = (eD->eVerRec)->eProbMin;
    (eD->eVerRec)->eProbMin = eD->eTProbMin;
    if ((vta = (eD->eVerRec)->AddTrack(*(eD->eWorking), eTr, zpos)))
    {
	eTr->AddVTA(vta);
	EdbVertex *eW = eD->eWorking;
	eW->SetID(eD->eVertex->ID());
	eW->V()->rmsDistAngle();
	sprintf(text,"New     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	eD->DrawOldBut("Original");
	if (eD->ePrevious)
	{
	    eD->DrawNewVTX(text);
	    eD->DrawNewBut("Modified");
	    eW = eD->ePrevious;
	    eW->V()->rmsDistAngle();
	    sprintf(text,"Pre     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	    eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	    eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	    eD->DrawPreVTX(text);
	    eD->DrawPreBut("Previous");
	}
	else
	{
	    eD->DrawPreVTX(text);
	    eD->DrawPreBut("Modified");
	    if (eD->eVertex->ID() >= 0)
	    {
		eD->DrawAcc();
		eD->DrawCan();
		eD->DrawUnd();
	    }
	}
//	eD->DrawOldBut("Original");
	eD->DrawVTXTracks("Modified", eD->eWorking);
	if (eD->eArrV && (eD->eIndVert >= 0)) 
	{
//	    eD->eArrV->RemoveAt(eD->eIndVert);
	    eD->eArrV->AddAt(eD->eWorking, eD->eIndVert);
	}	
	eD->Draw();
	if (ePreviousSaved) delete ePreviousSaved;
	ePreviousSaved = 0;
	(eD->eVerRec)->eImpMax = ImpMaxSave;
	(eD->eVerRec)->eProbMin = ProbMinSave;
    }
    else
    {
	printf("Track not added! May be Prob < ProbMin. Change ProbMin with 'TrackParams' button!\n");
	fflush(stdout);
	delete eD->eWorking;
	if (eD->ePrevious)
	{
	    eD->eWorking = eD->ePrevious;
	    (eD->eWorking)->ResetTracks();
	    eD->ePrevious = ePreviousSaved;
	}
	else
	{
	    eD->eWorking = 0;
	    (eD->eVertex)->ResetTracks();
	}
	(eD->eVerRec)->eImpMax = ImpMaxSave;
	(eD->eVerRec)->eProbMin = ProbMinSave;
	return;
    }
  }
}
//_____________________________________________________________________________
void EdbTrackG::InfoTrackVert()
{
  if (!(eD->eVertex))
  {
    
	    printf("No working vertex selected!\n");
	    fflush(stdout);
	    return;
  }
  int zpos = 1;
  if (GetMarkerColor() == kRed) zpos = 0;
  EdbVertex *v = eD->eVertex;
  if (eD->eWorking) v = eD->eWorking;
  char CanvasTRKName[140];
  strcpy(CanvasTRKName, "TRK-");
  strcat(CanvasTRKName, (eD->fCanvas)->GetName());
  if ((eD->fCanvasTRK = (TCanvas *)(gROOT->GetListOfCanvases()->FindObject(CanvasTRKName))))
  {
    (eD->fCanvasTRK)->SetTitle("Track - Vertex relation parameters");
    (eD->fCanvasTRK)->Clear();
    (eD->fCanvasTRK)->Modified();
    (eD->fCanvasTRK)->Update();
  }
  else
  {
    int xpos = (eD->fCanvas)->GetWindowTopX()+(eD->fCanvas)->GetWw();
    int ypos = (eD->fCanvas)->GetWindowTopY();
    eD->fCanvasTRK = new TCanvas(CanvasTRKName, "Track - Vertex relation parameters",
			     -xpos, ypos, 640, 330);
//    (eD->fCanvasTRK)->ToggleEventStatus();
  }
  if (eD->fVTXTRKInfo)
  {
    (eD->fVTXTRKInfo)->Clear();
  }
  else
  {
    eD->fVTXTRKInfo = new TPaveText(0.05, 0.05, 0.95, 0.95);
    (eD->fVTXTRKInfo)->ResetBit(kCanDelete);
  }
  char line[128];
  EdbTrackP *tr = eTr;  
  TText *t = 0;

  strcpy(line, " Track     ID   Nseg   Mass       P       Chi2/ndf    Prob");
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  sprintf(line,"         %4d  %4d   %7.4f  %7.2f    %5.2f     %7.4f",
		      tr->ID(), tr->N(), tr->M(), tr->P(),
		      tr->Chi2()/tr->N(), tr->Prob());
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

//  t = (eD->fVTXTRKInfo)->AddText("");
//  t->SetTextColor(kBlack);
//  t->SetTextSize(0.03);
//  t->SetTextAlign(12);
//  t->SetTextFont(102);

  strcpy(line, "  Vertex    ID    Mult  X          Y          Z          Dist   Chi2     Prob");
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  sprintf(line,"            %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
    v->ID(), v->N(), v->VX(), v->VY(), v->VZ(), v->V()->dist(),
    v->V()->chi2()/v->V()->ndf(), v->V()->prob());
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

//  t = (eD->fVTXTRKInfo)->AddText("");
//  t->SetTextColor(kBlack);
//  t->SetTextSize(0.03);
//  t->SetTextAlign(12);
//  t->SetTextFont(102);

  EdbSegP *seg = tr->TrackExtremity( zpos );//, eUseSegPar);
  float dx   = v->VX() - seg->X();
  float dy   = v->VY() - seg->Y();
  float dz   = v->VZ() - seg->Z();
  float dist = TMath::Sqrt(dx*dx + dy*dy + dz*dz);

  float mass = tr->M();
  float momentum = tr->P();
  if (eD->eP > 0.) seg->SetP(eD->eP);
  if (eD->eM > 0.) tr->SetM(eD->eM);

  float X0 = 0.;
  if (eD->eVerRec) if ((eD->eVerRec)->ePVR)  X0 = (((eD->eVerRec)->ePVR)->GetScanCond())->RadX0();

  float chi2 = v->Chi2Track(tr, zpos, X0);
  float impa = v->DistTrack(tr, zpos, X0);

  seg->SetP(momentum);
  tr->SetM(mass);
 
  sprintf(line, "  Track - Vertex  impact = %-6.1f, chi2 = %-7.1f, distance = %-8.1f", impa, chi2, dist);
  t = (eD->fVTXTRKInfo)->AddText(line);
  t->SetTextColor(kRed);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  t = (eD->fVTXTRKInfo)->AddText("");
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);

  (eD->fVTXTRKInfo)->Draw();
  (eD->fCanvasTRK)->Modified();
  (eD->fCanvasTRK)->Update();
}
//_____________________________________________________________________________
void EdbDisplay::CancelModifiedVTX()
{
    if (eArrVSave || eArrTrSave || eArrSegPSave) DrawAllObjects();
    if (eWorking)
    {
	delete eWorking;
	eWorking = 0;
    }
    if (ePrevious)
    {
	delete ePrevious;
	ePrevious = 0;
    }
    if (eVertex)
    {
	if (eVertex->ID() < 0)
	{
	    if (eArrV) 
	    {
	    	eArrV->Remove(eVertex);
		eArrV->Compress();
	    }	
	    delete eVertex;
	    eVertex = 0;
	}
	else
	{
	    (eVertex)->ResetTracks();
	}
    }
    TList *li = fTrigPad->GetListOfPrimitives();
    if (li)
    {
	li->Remove(fUndButton);
	li->Remove(fAccButton);
	li->Remove(fCanButton);
	li->Remove(fEnvButton);
    }

    EdbTrackP *tr = 0;

    for (int i = 0; i < eCreatedTracks.GetSize(); i++)
    {
	tr = (EdbTrackP *)eCreatedTracks.At(i);
	if (tr)
	{
	    if (eArrTr) if (eArrTr->FindObject(tr))
	    {
		eArrTr->Remove(tr);
		eArrTr->Compress();
	    }
	    if (eArrTrSave) if (eArrTrSave->FindObject(tr))
	    {
		eArrTrSave->Remove(tr);
		eArrTrSave->Compress();
	    }
	    delete tr;
	} 
    }
    eCreatedTracks.Clear();

    if (eArrV && (eIndVert >= 0) && eVertex) 
    {
//	    eArrV->RemoveAt(eIndVert);
	    eArrV->AddAt(eVertex, eIndVert);
    }	
    fTrigPad->cd();
    fTrigPad->Update();
    fTrigPad->Draw();
    fPad->cd();
    if (eWait_Answer) CloseDialogModifiedVTX();
    if (fCanvasVTX) fCanvasVTX->Close();
    fCanvasVTX = 0;
    eIndVert = -1;
    eVertex = 0;
    Draw();
}

//_____________________________________________________________________________
void EdbDisplay::DeleteModifiedVTX()
{
    if (eWorking)
    {
	delete eWorking;
	eWorking = 0;
    }
    if (ePrevious)
    {
	delete ePrevious;
	ePrevious = 0;
    }
    if (eVertex)
    {
	if (eVertex->ID() < 0)
	{
	    if (eArrV) 
	    {
	    	eArrV->Remove(eVertex);
		eArrV->Compress();
	    }	
	    delete eVertex;
	    eVertex = 0;
	}
	else
	{
	    (eVertex)->ResetTracks();
	}
    }

    EdbTrackP *tr = 0;

    for (int i = 0; i < eCreatedTracks.GetSize(); i++)
    {
	tr = (EdbTrackP *)eCreatedTracks.At(i);
	if (tr)
	{
	    if (eArrTr) if (eArrTr->FindObject(tr))
	    {
		eArrTr->Remove(tr);
		eArrTr->Compress();
	    }
	    if (eArrTrSave) if (eArrTrSave->FindObject(tr))
	    {
		eArrTrSave->Remove(tr);
		eArrTrSave->Compress();
	    }
	    delete tr;
	} 
    }
    eCreatedTracks.Clear();

    if (eArrVSave || eArrTrSave || eArrSegPSave)
    {
	delete eArrV;
	eArrV = eArrVSave;
	eArrVSave = 0;
	delete eArrTr;
	eArrTr = eArrTrSave;
	eArrTrSave = 0;
	delete eArrSegP;
	eArrSegP = eArrSegPSave;
	eArrSegPSave = 0;
	eIndVert = eIndVertSave;
    }
    if (eArrV && (eIndVert >= 0) && eVertex) 
    {
//	    eArrV->RemoveAt(eIndVert);
	    eArrV->AddAt(eVertex, eIndVert);
	    eIndVert = -1;
    }	
}
//_____________________________________________________________________________
void EdbDisplay::UndoModifiedVTX()
{
    char text[512];
    EdbTrackP *LastCreated = 0;
    int CreatedInd = eCreatedTracks.GetEntries();
    if (CreatedInd > 0) LastCreated = (EdbTrackP *)eCreatedTracks.At(CreatedInd-1); 
    if (ePrevious)
    {
	int InWork = 0;
	if (LastCreated)
	{
	    for (int i=0; i<eWorking->N(); i++)
		if (eWorking->GetTrack(i) == LastCreated)
		    InWork = 1;
	}
	int InPrev = 0;
	if (LastCreated)
	{
	    for (int i=0; i<ePrevious->N(); i++)
		if (ePrevious->GetTrack(i) == LastCreated)
		    InPrev = 1;
	}
	delete eWorking;
	if (InWork && !InPrev)
	{
	    eCreatedTracks.Remove(LastCreated);
	    if (eArrTr) if (eArrTr->FindObject(LastCreated))
	    {
		eArrTr->Remove(LastCreated);
		eArrTr->Compress();
	    }
	    if (eArrTrSave) if (eArrTrSave->FindObject(LastCreated))
	    {
		eArrTrSave->Remove(LastCreated);
		eArrTrSave->Compress();
	    }
	    delete LastCreated;
	}
	eWorking = ePrevious;
	(eWorking)->ResetTracks();
	EdbVertex *eW = eWorking;
	eW->V()->rmsDistAngle();
	sprintf(text,"New     %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	DrawPreVTX(text);
	DrawPreBut("Modified");
	ClearNewVTX();
	if (fNewBut)
	{
	    delete fNewBut;
	    fNewBut = 0;
	}
	ePrevious = 0;
	DrawVTXTracks("Modified", eWorking);
	if (eIndVert >= 0) 
	{
//	    eArrV->RemoveAt(eIndVert);
	    eArrV->AddAt(eWorking, eIndVert);
	}
	Draw();	
    }
    else if (eWorking)
    {
	EdbTrackP *tr = eWorking->GetTrack(eWorking->N() - 1);
	delete eWorking;
	if (tr && eCreatedTracks.FindObject(tr))
	{
	    eCreatedTracks.Remove(tr);
	    if (eArrTr) if (eArrTr->FindObject(tr))
	    {
		eArrTr->Remove(tr);
		eArrTr->Compress();
	    }
	    if (eArrTrSave) if (eArrTrSave->FindObject(tr))
	    {
		eArrTrSave->Remove(tr);
		eArrTrSave->Compress();
	    }
	    delete tr;
	}
	(eVertex)->ResetTracks();
	EdbVertex *eW = eVertex;
	eW->V()->rmsDistAngle();
	sprintf(text,"Orig    %-4d  %-4d  %-8.1f   %-8.1f   %-8.1f   %-6.1f %-7.1f  %-7.5f",
	eW->ID(), eW->N(), eW->VX(), eW->VY(), eW->VZ(), eW->V()->dist(),
	eW->V()->chi2()/eW->V()->ndf(), eW->V()->prob());
	DrawOldVTX(text);
	ClearPreVTX();
	if (fPreBut)
	{
	    delete fPreBut;
	    fPreBut = 0;
	}
	if (fOldBut)
	{
	    delete fOldBut;
	    fOldBut = 0;
	}
	eWorking = 0;
	DrawVTXTracks("Original", eVertex);
	if (eArrV && (eIndVert >= 0)) 
	{
//	    eArrV->RemoveAt(eIndVert);
	    eArrV->AddAt(eVertex, eIndVert);
	}
	TList *li = fTrigPad->GetListOfPrimitives();
	if (li)
	{
	    li->Remove(fUndButton);
	    li->Remove(fAccButton);
	    li->Remove(fCanButton);
	}
	Draw();	
    }
}
//_____________________________________________________________________________
void EdbDisplay::AcceptModifiedVTX()
{
    EdbVertex *eVe = 0;
    if (eArrVSave || eArrTrSave || eArrSegPSave) DrawAllObjects();
    if (ePrevious && (ePrevious != eVertex) && (ePrevious != eWorking))
    {
	delete ePrevious;
	ePrevious = 0;
    }
    if (eWorking && eVertex)
    {
	//if (!eVerRec) eVerRec = ((EdbVertexRec*)(gROOT->GetListOfSpecials()->FindObject("EdbVertexRec")));
	if (eVerRec) if (eVerRec->IsA() != EdbVertexRec::Class()) eVerRec = 0;
	if (!eVerRec) {printf("Error: EdbDisplay:AcceptModifiedVTX: EdbVertexRec not defined, use SetVerRec(...)\n"); return;}
        EdbVertex *eW = eWorking;
	int ind = eVertex->ID();
	if (ind < 0)
	{
	    if (eVerRec->eVTX)
	    {
		ind = eVerRec->eVTX->GetEntries();
	    }
	    else
	    {
	        eVerRec->eVTX = new TObjArray(10);
	    	ind = 0;
	    }
	}
	eW->SetID(ind);
	eW->SetQuality(eW->V()->prob()/
			   (eW->V()->vtx_cov_x()+eW->V()->vtx_cov_y()));
	if (eVerRec)
	{
	    int ntr = eVertex->N();
//	    printf("ind %d ntr %d\n", ind, ntr);
//	    fflush(stdout);
	    for(int i=0; i<ntr; i++)
	    {
//		printf("    %d\n", i);
//		fflush(stdout);
		(eVerRec->eVTA).Remove(eVertex->GetVTa(i));
	    }
	}
	int indd = -1;
//	if (!eArrV) eArrV = new TObjArray(10);
	if (eArrV) indd = eArrV->IndexOf(eVertex);
	if (indd >= 0)
	{
//	    eArrV->RemoveAt(indd);
	    eArrV->AddAt(eW, indd);
	}
	eVe = eVertex;
	EdbTrackP *tr = 0;

	TObjArray *etr = 0;
	if (eVerRec) etr = eVerRec->eEdbTracks;
	int trind = 0;
	if (etr) trind = etr->GetEntries();
	for (int i = 0; i < eCreatedTracks.GetSize(); i++)
	{
	    tr = (EdbTrackP *)(eCreatedTracks.At(i));
	    if (tr)
	    {
//		    printf("i %d id %d\n", i, tr->ID());
//		    fflush(stdout);
//		if (tr->VTAS() || tr->VTAE())
//		{
		    tr->SetID(trind++);
		    if (etr) etr->Add(tr);
		    if (eArrTr) eArrTr->Add(tr);
//		    printf("     id %d\n", tr->ID());
//		    fflush(stdout);
//		}
//		else delete tr;
	    } 
	}
	eCreatedTracks.Clear();

	delete eVertex; //?????????
	eVertex = 0;

	eW->ResetTracks();

	int ntr = eW->N();
	int ifl = 0;
	for(int i=0; i<ntr; i++)
	{
		tr = eW->GetTrack(i);
		if (eTrack)
		{
		  if (tr == eTrack)
		  {
		    TObjArray *etr = 0;
		    if (eVerRec) etr = eVerRec->eEdbTracks;
		    int trind = 0;
		    if (etr) trind = etr->GetEntries();
		    eTrack->SetID(trind);
		    if (etr) etr->Add(eTrack);
		    eTrack->SetSegmentsTrack();
		    eTrack = 0;
		  }
		}
		if (eW->Zpos(i)) ifl = ifl | 1; 
		else		 ifl = ifl | 2; 
		if (tr && eArrTr)
		{
		    indd = eArrTr->IndexOf(tr);
		    if (indd < 0)
		    {
			eArrTr->Add(tr);
		    }
		}
	}
	ifl = 4 - ifl;
	if (ifl == 3) ifl = 0;
	if (eW->Nv()) ifl += 3;
	eW->SetFlag(ifl);
	if (eVerRec)
	{
//	    eVerRec->eVTX->RemoveAt(ind);
	    if (eVerRec->eVTX) eVerRec->eVTX->AddAt(eW, ind);
	    for(int i=0; i<ntr; i++)
	    {
		eVerRec->AddVTA(eW->GetVTa(i));
	    }
	}
    }
    if (eVertex && eVertex->ID() < 0 && !eWorking)
    {
	int ind = 0;
	if (eVerRec->eVTX)
	    {
		ind = eVerRec->eVTX->GetEntries();
	    }
	else
	    {
	        eVerRec->eVTX = new TObjArray(10);
	    	ind = 0;
	    }
	eVertex->SetID(ind);
//	if (!eArrV) eArrV = new TObjArray(10);
	if (eArrV) eArrV->Add(eVertex);
	EdbTrackP *tr = 0;
	TObjArray *etr = 0;
	if (eVerRec) etr = eVerRec->eEdbTracks;
	int trind = 0;
	if (etr) trind = etr->GetEntries();
	for (int i = 0; i < eCreatedTracks.GetSize(); i++)
	{
	    tr = (EdbTrackP *)(eCreatedTracks.At(i));
	    if (tr)
	    {
//		if (tr->VTAS() || tr->VTAE())
//		{
		    tr->SetID(trind++);
		    if (etr) etr->Add(tr);
		    if (eArrTr) eArrTr->Add(tr);
//		}
//		else delete tr;
	    } 
	}
	eCreatedTracks.Clear();

	eVertex->ResetTracks();
	int ntr = eVertex->N();
	int ifl = 0;
	int indd = 0;
	for(int i=0; i<ntr; i++)
	{
		tr = eVertex->GetTrack(i);
		if (eTrack)
		{
		  if (tr == eTrack)
		  {
		    TObjArray *etr = 0;
		    if (eVerRec) etr = eVerRec->eEdbTracks;
		    int trind = 0;
		    if (etr) trind = etr->GetEntries();
		    eTrack->SetID(trind);
		    if (etr) etr->Add(eTrack);
		    eTrack->SetSegmentsTrack();
		    eTrack = 0;
		  }
		}
		if (eVertex->Zpos(i)) ifl = ifl | 1; 
		else		      ifl = ifl | 2; 
		if (tr && eArrTr)
		{
		    indd = eArrTr->IndexOf(tr);
		    if (indd < 0)
		    {
			eArrTr->Add(tr);
		    }
		}
	}
	ifl = 4 - ifl;
	if (ifl == 3) ifl = 0;
	if (eVertex->Nv()) ifl += 3;
	eVertex->SetFlag(ifl);
	if (eVerRec)
	{
//	    eVerRec->eVTX->RemoveAt(ind);
	    if (!(eVerRec->eVTX)) eVerRec->eVTX = new TObjArray(10);
	    eVerRec->eVTX->Add(eVertex);
	    for(int i=0; i<ntr; i++)
	    {
		eVerRec->AddVTA(eVertex->GetVTa(i));
	    }
	}
    }
    eWorking = 0;
    eVertex = 0;
    TList *li = fTrigPad->GetListOfPrimitives();
    if (li)
    {
	li->Remove(fUndButton);
	li->Remove(fAccButton);
	li->Remove(fCanButton);
	li->Remove(fEnvButton);
    }
    fTrigPad->cd();
    fTrigPad->Update();
    fTrigPad->Draw();

    fPad->cd();
    if (eWait_Answer) CloseDialogModifiedVTX();
    if (fCanvasVTX) fCanvasVTX->Close();
    fCanvasVTX = 0;
    Draw();
    //if (eVe) delete eVe;
    eIndVert = -1;
}
//_____________________________________________________________________________
void EdbDisplay::DrawVertexEnvironment()
{
    if (eRadMax == 0)
    {
	printf("No neighborhood in EdbPVRec object!\n");
	fflush(stdout);
	return;
    }
    //if (!eVerRec) eVerRec = ((EdbVertexRec *)(gROOT->GetListOfSpecials()->FindObject("EdbVertexRec")));
    if (eVerRec) if((eVerRec->IsA()) != EdbVertexRec::Class()) eVerRec = 0;
    if (!eVerRec)
       {printf("Error: EdbDisplay:DrawVertexEnvironment: EdbVertexRec not defined, use SetVerRec(...)\n"); fflush(stdout); return;}
    if (!eVertex && !eSegment)
       {printf("Error: EdbDisplay:DrawVertexEnvironment: Working vertex (or segment) not defined\n"); fflush(stdout); return;}
    
    fTrigPad->cd();
    fTrigPad->GetListOfPrimitives()->Remove(fEnvButton);
    if (fStyle/2 == 1)
	fAllButton->SetFillColor(33);
    else
	fAllButton->SetFillColor(38);
    fTrigPad->GetListOfPrimitives()->Add(fAllButton);
    fAllButton->SetPad(0.05,0.47,0.95,0.56);
    fAllButton->Draw();
    fTrigPad->Modified(kTRUE);
    fTrigPad->Update();
    fTrigPad->Draw();
    fPad->cd();

    EdbVertex *eW = eVertex;
    EdbTrackP *tr = 0;
    if (eWorking) eW = eWorking;
    float Rmax = eRadMax, ImpMax = eImpMax;
    int Dpat = (int)eDpat;

    eArrTrSave  = eArrTr; 
    eArrSegPSave  = eArrSegP; 
    eArrVSave  = eArrV;
    eIndVertSave = eIndVert; 
    eArrV = new TObjArray(20);
    eArrTr = new TObjArray(20);
    eArrSegP = new TObjArray(20);

    if (eW)
    {
	eVerRec->VertexNeighbor(eW, Rmax, Dpat, ImpMax);
	eArrV->Add(eW);
	int ntr = eW->N();
	for(int i=0; i<ntr; i++)
	{
	    eArrTr->Add((tr = eW->GetTrack(i)));
	}
	int nntr = eW->Nn();
	EdbVTA *vta = 0;
	for(int i=0; i<nntr; i++)
	{
	if ((vta = eW->GetVTn(i)))
	{
	    if (vta->Flag() == 0) //track
	    {
		eArrTr->Add((tr = vta->GetTrack()));
	    }
	    else if (vta->Flag() == 1) //segment
	    {
		eArrSegP->Add((EdbSegP *)(vta->GetTrack()));
	    }
	    else if (vta->Flag() == 3) //vertex
	    {
		eArrV->Add((EdbVertex *)(vta->GetTrack()));
	    }
	}
	}
	eIndVert = eArrV->IndexOf(eW);
    }
    else if (eSegment)
    {
	eIndVert = -1;
	eVerRec->SegmentNeighbor(eSegment, Rmax, Dpat, ImpMax, eSegWmin, eArrSegP, eArrTr, eArrV);
	if (eArrV->GetEntries()) eDrawVertex = 1;
    }
    Draw();
}
//_____________________________________________________________________________
void EdbDisplay::DrawAllObjects()
{
    fTrigPad->cd();
    fTrigPad->GetListOfPrimitives()->Remove(fAllButton);
    if (fStyle/2 == 1)
	fEnvButton->SetFillColor(33);
    else
	fEnvButton->SetFillColor(38);
    fTrigPad->GetListOfPrimitives()->Add(fEnvButton);
    fEnvButton->SetPad(0.05,0.47,0.95,0.56);
    fEnvButton->Draw();
    fTrigPad->Modified(kTRUE);
    fTrigPad->Update();
    fTrigPad->Draw();
    fPad->cd();

    if (eArrVSave || eArrTrSave || eArrSegPSave)
    {
	delete eArrV;
	eArrV = eArrVSave;
	eArrVSave = 0;
	eIndVert = eIndVertSave;
	delete eArrTr;
	eArrTr = eArrTrSave;
	eArrTrSave = 0;
	delete eArrSegP;
	eArrSegP = eArrSegPSave;
	eArrSegPSave = 0;
    }
    if (eArrV && eWorking && (eIndVert >= 0)) 
    {
//	eArrV->RemoveAt(eIndVert);
	eArrV->AddAt(eWorking, eIndVert);
    }
    Draw();    
}
//_____________________________________________________________________________
void EdbDisplay::DialogNeighborParameters()
{
    if (eRadMax == 0)
    {
	printf("No neighborhood in EdbPVGen object!\n");
	fflush(stdout);
	return;
    }

    eWait_Answer = true;

    fMain = new TGMainFrame(gClient->GetRoot(), 260, 340);
    TGMainFrame *fTra = fMain;

//    fTra->Connect("CloseWindow()", "EdbDisplay", this, "CloseDialogModifiedParams()");

   // use hierarchical cleaning
//    fTra->SetCleanup(kDeepCleanup);

    TGHorizontalFrame *fF[4] = {0,0,0,0};
    TGLabel *fLabel[4] = {0,0,0,0};
    Double_t parinit[4] = {eRadMax, eDpat,  eImpMax, eSegWmin};
    Double_t parmax[4] =  { 150000.,   50.,   10000., 1000. };
    Double_t parmin[4] =  {    0.,    0.,      0.,    0. };
    //if (eVerRec) if (eVerRec->ePVR) parmax[1] = ((eVerRec->ePVR)->Npatterns()-1);
    char *parlabel[4] = {"Maximal radius","+/- patterns", "Maximal impact", "Minimal seg W"};

    TGGC myGC = *(gClient->GetResourcePool()->GetFrameGC());
    TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
    if (myfont) myGC.SetFont(myfont->GetFontHandle());

    TGVerticalFrame *fF1 = new TGVerticalFrame(fTra, 240, 250);
    TGLayoutHints *fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
    fTra->AddFrame(fF1, fL1);
    TGLayoutHints *fL2 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2);
    for (int i = 0; i < 4; i++) {
      fF[i] = new TGHorizontalFrame(fF1, 220, 35);
      fF1->AddFrame(fF[i], fL2);
      fNumericEntries[i] = new TGNumberEntry(fF[i], parinit[i], 12, i+1000,
                                 TGNumberFormat::kNESInteger,
				 TGNumberFormat::kNEANonNegative,
				 TGNumberFormat::kNELLimitMinMax,
				 parmin[i], parmax[i]);
      fF[i]->AddFrame(fNumericEntries[i], fL2);
      fLabel[i] = new TGLabel(fF[i], parlabel[i], myGC(), myfont->GetFontStruct());
      fF[i]->AddFrame(fLabel[i], fL2);
    }

    TGHorizontalFrame *fFrame = new TGHorizontalFrame(fTra, 220, 30);
    char cmda[256];
    sprintf(cmda,
    "((EdbDisplay*)((gROOT->GetListOfSpecials())->FindObject(\"%s\")))->AcceptModifiedParams()",fTitle);
    TGTextButton *abut = new TGTextButton(fFrame, "&Accept", cmda);
    fFrame->AddFrame(abut,
		     new TGLayoutHints(kLHintsCenterY | kLHintsLeft,
		     5,5,3,4));

    char cmdc[256];
    sprintf(cmdc,
    "((EdbDisplay*)((gROOT->GetListOfSpecials())->FindObject(\"%s\")))->CancelDialogModifiedParams()",fTitle);
    TGTextButton *cbut = new TGTextButton(fFrame, "&Cancel", cmdc);
    fFrame->AddFrame(cbut,
		     new TGLayoutHints(kLHintsCenterY | kLHintsRight,
    		     5, 5, 4, 4));
//    fFrame->Resize(220, 30);

    fTra->AddFrame(fFrame,
		   new TGLayoutHints(kLHintsCenterY | kLHintsCenterX,
		   5, 5, 4, 4));

    fTra->MapSubwindows();
    fTra->Resize(260, fTra->GetDefaultHeight());
    fTra->SetWindowName("Neighborhood Parameters");
    fTra->MapWindow();
    gClient->WaitFor(fTra);
}
//_____________________________________________________________________________
void EdbDisplay::DialogTrackParameters()
{
    eWait_Answer = true;

    fMain = new TGMainFrame(gClient->GetRoot(), 280, 360);
    TGMainFrame *fTra = fMain;

//    fTra->Connect("CloseWindow()", "EdbDisplay", this, "CloseDialogModifiedParams()");

   // use hierarchical cleaning
//    fTra->SetCleanup(kDeepCleanup);

    TGHorizontalFrame *fF[4] = {0,0,0,0};
    TGLabel *fLabel[4] = {0,0,0,0};
    Double_t parinit[4] = {    eP,    eM, eTImpMax, eTProbMin};
    Double_t parmax[4] =  { 1000.,   10., 1000000.,        1.};
    Double_t parmin[4] =  {    0.,    0.,       0.,        0.};
    char *parlabel[4] = {"Momentum","Mass", "MaxImpact", "MinProb"};

    TGGC myGC = *(gClient->GetResourcePool()->GetFrameGC());
    TGFont *myfont = gClient->GetFont("-adobe-helvetica-bold-r-*-*-12-*-*-*-*-*-iso8859-1");
    if (myfont) myGC.SetFont(myfont->GetFontHandle());

    TGVerticalFrame *fF1 = new TGVerticalFrame(fTra, 260, 250);
    TGLayoutHints *fL1 = new TGLayoutHints(kLHintsTop | kLHintsLeft, 2, 2, 2, 2);
    fTra->AddFrame(fF1, fL1);
    TGLayoutHints *fL2 = new TGLayoutHints(kLHintsCenterY | kLHintsLeft, 2, 2, 2, 2);
    for (int i = 0; i < 4; i++) {
      fF[i] = new TGHorizontalFrame(fF1, 240, 35);
      fF1->AddFrame(fF[i], fL2);
      fNumericEntries[i] = new TGNumberEntry(fF[i], parinit[i], 12, i+2000,
                                 TGNumberFormat::kNESReal,
				 TGNumberFormat::kNEANonNegative,
				 TGNumberFormat::kNELLimitMinMax,
				 parmin[i], parmax[i]);
      fF[i]->AddFrame(fNumericEntries[i], fL2);
      fLabel[i] = new TGLabel(fF[i], parlabel[i], myGC(), myfont->GetFontStruct());
      fF[i]->AddFrame(fLabel[i], fL2);
    }

    TGHorizontalFrame *fFrame = new TGHorizontalFrame(fTra, 240, 30);
    char cmda[256];
    sprintf(cmda,
//    "((EdbDisplay*)(gROOT->GetListOfSpecials()->FindObject((TObject *)0x%08x)))->AcceptModifiedTrackParams()",this);
    "((EdbDisplay*)((gROOT->GetListOfSpecials())->FindObject(\"%s\")))->AcceptModifiedTrackParams()",fTitle);
    TGTextButton *abut = new TGTextButton(fFrame, "&Accept", cmda);
    fFrame->AddFrame(abut,
		     new TGLayoutHints(kLHintsCenterY | kLHintsLeft,
		     5,5,3,4));

    char cmdc[256];
    sprintf(cmdc,
    "((EdbDisplay*)((gROOT->GetListOfSpecials())->FindObject(\"%s\")))->CancelDialogModifiedTrackParams()",fTitle);
    TGTextButton *cbut = new TGTextButton(fFrame, "&Cancel", cmdc);
    fFrame->AddFrame(cbut,
		     new TGLayoutHints(kLHintsCenterY | kLHintsRight,
    		     5, 5, 4, 4));
//    fFrame->Resize(220, 30);

    fTra->AddFrame(fFrame,
		   new TGLayoutHints(kLHintsCenterY | kLHintsCenterX,
		   5, 5, 4, 4));

    fTra->MapSubwindows();
    fTra->Resize(280, fTra->GetDefaultHeight());
    fTra->SetWindowName("Track Parameters");
    fTra->MapWindow();
    gClient->WaitFor(fTra);
}
//_____________________________________________________________________________
void EdbDisplay::AcceptModifiedParams()
{
    eWait_Answer = false;
    if(fNumericEntries[0]) eRadMax = fNumericEntries[0]->GetIntNumber();
    if(fNumericEntries[1]) eDpat   = fNumericEntries[1]->GetIntNumber();
    if(fNumericEntries[2]) eImpMax = fNumericEntries[2]->GetIntNumber();
    if(fNumericEntries[3]) eSegWmin= fNumericEntries[3]->GetIntNumber();

    for (Int_t i = 0; i < 4; i++) SafeDelete(fNumericEntries[i]);

    fMain->SendCloseMessage();
    fMain = 0;
}
//_____________________________________________________________________________
void EdbDisplay::CloseDialogModifiedParams()
{
    eWait_Answer = false;
    fMain = 0;
}
//_____________________________________________________________________________
void EdbDisplay::CancelDialogModifiedParams()
{
    eWait_Answer = false;
    for (Int_t i = 0; i < 3; i++) SafeDelete(fNumericEntries[i]);
    fMain->SendCloseMessage();
    fMain = 0;
}
//_____________________________________________________________________________
void EdbDisplay::AcceptModifiedTrackParams()
{
    eWait_Answer = false;
    if(fNumericEntries[0]) eP = fNumericEntries[0]->GetNumber();
    if(fNumericEntries[1]) eM = fNumericEntries[1]->GetNumber();
    if(fNumericEntries[2]) eTImpMax  = fNumericEntries[2]->GetNumber();
    if(fNumericEntries[3]) eTProbMin = fNumericEntries[3]->GetNumber();
    for (Int_t i = 0; i < 4; i++) SafeDelete(fNumericEntries[i]);
    fMain->SendCloseMessage();
    fMain = 0;
}
//_____________________________________________________________________________
void EdbDisplay::CloseDialogModifiedTrackParams()
{
    eWait_Answer = false;
    fMain = 0;
}
//_____________________________________________________________________________
void EdbDisplay::CancelDialogModifiedTrackParams()
{
    eWait_Answer = false;
    for (Int_t i = 0; i < 4; i++) SafeDelete(fNumericEntries[i]);
    fMain->SendCloseMessage();
    fMain = 0;
}
//_____________________________________________________________________________
void EdbDisplay::DialogModifiedVTX()
{
    eWait_Answer = true;
    fMain = new TGMainFrame(gClient->GetRoot(), 250, 40);
    TGMainFrame *fTra = fMain;
    TGHorizontalFrame *fFrame = new TGHorizontalFrame(fTra, 220, 30);
    char cmda[256];
    sprintf(cmda,
    "((EdbDisplay*)(gROOT->GetListOfSpecials()->FindObject(\"%s\")))->AcceptModifiedVTX()",fTitle);
    TGTextButton *abut = new TGTextButton(fFrame, "&Accept", cmda);
    fFrame->AddFrame(abut,
		     new TGLayoutHints(kLHintsCenterY | kLHintsLeft,
		     5,5,3,4));

    char cmdc[256];
    sprintf(cmdc,
    "((EdbDisplay*)(gROOT->GetListOfSpecials()->FindObject(\"%s\")))->CancelModifiedVTX()",fTitle);
    TGTextButton *cbut = new TGTextButton(fFrame, "&Cancel", cmdc);
    fFrame->AddFrame(cbut,
		     new TGLayoutHints(kLHintsCenterY | kLHintsRight,
    		     5, 5, 4, 4));
    fFrame->Resize(210, 30);

    fTra->AddFrame(fFrame,
		   new TGLayoutHints(kLHintsCenterY | kLHintsCenterX,
		   5, 5, 4, 4));

    fTra->MapSubwindows();
    fTra->Resize(220, fTra->GetDefaultHeight());
    fTra->SetWindowName("Modified Vertex Exist!");
    fTra->MapWindow();
}
//_____________________________________________________________________________
void EdbDisplay::CloseDialogModifiedVTX()
{
    eWait_Answer = false;
    fMain->SendCloseMessage();
    fMain = 0;
}
//=============================================================================
void EdbDisplay::DrawVTXTracks(char *type, EdbVertex *v)
{
  EdbVertex *vv = v;
  if (!vv)
  {
    if (!strcmp(type, "Original") && eVertex)
    {
	vv = eVertex;
    }
    else if (!strcmp(type, "Previous") && ePrevious)
    {
	vv = ePrevious;
    }
    else if (!strcmp(type, "Modified") && eWorking)
    {
	vv = eWorking;
    }
  }
  if (!vv) return;
  bool drawrem = true;
  if (vv == eVertex && eWorking) drawrem = false;
  if (vv == ePrevious) drawrem = false;

  fCanvasVTX->cd();

  if (fVTXTracks)
  {
    fVTXTracks->Clear();
  }
  else
  {
    fVTXTracks = new TPaveText(0.05, 0.05, 0.95, 0.72);
    fVTXTracks->ResetBit(kCanDelete);
  }
  char line[128];
  EdbTrackP *tr = 0;  
  TText *t = 0;
  int ntr = vv->N();
  sprintf(line, "           Tracks parameters for %s vertex", type);
  t = fVTXTracks->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);
  strcpy(line, "    #   ID   Nseg   Mass       P       Chi2/ndf    Prob     Chi2Contrib     Impact");
  t = fVTXTracks->AddText(line);
  t->SetTextColor(kBlue);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);
  for (int i=0; i<ntr; i++)
  {
    tr = vv->GetTrack(i);
    sprintf(line,"%4d  %4d  %4d   %7.4f  %7.2f    %5.2f     %7.4f       %6.3f    %7.2f",
		   i, tr->ID(), tr->N(), tr->M(), tr->P(),
		      tr->Chi2()/tr->N(), tr->Prob(), vv->V()->track_chi2(i), vv->Impact(i));
    t = fVTXTracks->AddText(line);
    t->SetTextColor(kBlack);
    t->SetTextSize(0.03);
    t->SetTextAlign(12);
    t->SetTextFont(102);
  }
  t = fVTXTracks->AddText("");
  t->SetTextColor(kBlack);
  t->SetTextSize(0.03);
  t->SetTextAlign(12);
  t->SetTextFont(102);
  fVTXTracks->Draw();

  float dy = 0.67/(2+ntr+1);
  for (int i=0; i<50; i++)
  {
    if (!drawrem)
    {
	if (fRemBut[i])
	{
	    fRemBut[i]->Delete();
	    fRemBut[i] = 0;
	}
	continue;
    }
    if ( i < ntr )
    {
      if (!fRemBut[i])
      {
	char but[256], tit[32]="Rem";
	sprintf(but,
	"((EdbDisplay*)(gROOT->GetListOfSpecials()->FindObject(\"%s\")))->RemoveTrackFromTable(%d)",fTitle,i);
	fRemBut[i] = new TButton(tit, but, 0.88, 0.635-(i+2)*dy, 0.94, 0.675-(i+2)*dy);
	fRemBut[i]->SetToolTipText("Remove corresponding track");
	fRemBut[i]->ResetBit(kCanDelete);
	fRemBut[i]->SetFillColor(38);
	sprintf(tit,"RemTr%d",i);
	fRemBut[i]->SetName(tit);
	fRemBut[i]->Draw();
      }
      else
      {
        fRemBut[i]->SetPad(0.88,0.635-(i+2)*dy,0.94,0.675-(i+2)*dy);
	fRemBut[i]->Draw();
      }
    }
    else if (fRemBut[i])
    {
	fRemBut[i]->Delete();
	fRemBut[i] = 0;
    }
  }

  fPad->Modified(kTRUE);
  fPad->Update();
}

//=============================================================================
void EdbDisplay::SelectVertexTracks(TObjArray *vtx)
{
  if (!vtx) return;

  if (!eArrTr) eArrTr = new TObjArray();
  else eArrTr->Clear();

  Int_t nv = vtx->GetEntries();

  for (Int_t i = 0; i < nv; i++) {
    EdbVertex *vertex = (EdbVertex*)(vtx->At(i));
    if (!vertex || vertex->Flag() < 0) continue;
    for (Int_t j = 0; j < vertex->N(); j++) eArrTr->Add(vertex->GetTrack(j));
  }
}

//=============================================================================
void EdbDisplay::ClearSegmentEnv()
{
	if (fPad->GetListOfPrimitives()->FindObject(eSegPM))
	{
	    fPad->GetListOfPrimitives()->Remove(eSegPM);
	}
	delete eSegPM;
	eSegPM = 0;
	TList *li = fTrigPad->GetListOfPrimitives();
	if (eArrVSave || eArrSegPSave || eArrTrSave)
	{
	    li->Remove(fAllButton);
	    delete eArrV;
	    eArrV = eArrVSave;
	    eArrVSave = 0;
	    if (eArrTrSave)
	    {
		delete eArrTr;
		eArrTr = eArrTrSave;
		eArrTrSave = 0;
	    }
	    if (eArrSegPSave)
	    {
		delete eArrSegP;
		eArrSegP = eArrSegPSave;
		eArrSegPSave = 0;
	    }
	    Draw();    
	}
	else li->Remove(fEnvButton);
}
