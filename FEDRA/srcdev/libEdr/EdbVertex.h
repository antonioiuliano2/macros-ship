#ifndef ROOT_EdbVertex
#define ROOT_EdbVertex

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbVertex                                                            //
//                                                                      //
// Class to reconstruct vertexes                                        //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include <TObject.h>
#include <TTree.h>
#include <TVector3.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include "EdbPattern.h"
#include "vt++/VtTrack.hh"
#include "vt++/VtVertex.hh"
#include "vt++/VtRelation.hh"
#include "vt++/VtKalman.hh"

class EdbPVRec;

//_________________________________________________________________________
class EdbVTA: public TObject {

 private:

  EdbTrackP *eTrack;          // pointer to track (or segment)
  EdbVertex *eVertex; 	      // pointer to vertex

  Int_t      eZpos;            // 1-track start, 0-track end connect to the vertex
  Int_t	     eFlag;            // 0-neighbor track;
		               // 1-neighbor segment;
			       // 2-direct (attached) track connection
			       // 3-neighbor vertex;
  Float_t    eImp;             // impact parameter
  Float_t    eDist;            // distance from vertex to the nearest track point
  
 public:
  EdbVTA();
  EdbVTA(EdbVTA &vta);
  EdbVTA(EdbTrackP *tr, EdbVertex *v );
  virtual ~EdbVTA();

  Int_t      Zpos()      const  {return eZpos;}
  Int_t      Flag()      const  {return eFlag;}
  Float_t    Imp()       const  {return eImp;}
  Float_t    Dist()      const  {return eDist;}
  EdbTrackP *GetTrack()  const  {return eTrack;}
  EdbVertex *GetVertex() const  {return eVertex;}

  void Set0();
  void SetZpos(int zpos)       {eZpos = zpos;}
  void SetFlag(int flag)       {eFlag = flag;}
  void SetImp(float imp)       {eImp = imp;}
  void SetDist(float dist)     {eDist = dist;}
  void SetTrack(EdbTrackP *tr) {eTrack = tr;}
  void SetVertex(EdbVertex *v) {eVertex = v;}
  void Print();

  void AddVandT();

  ClassDef(EdbVTA,1)  // vertex-track association
};

//_________________________________________________________________________
class EdbVertex: public TObject {

 private:

  TList   eVTn;       // vertex neighborhood tracks and segments
  TList   eVTa;       // attached tracks

  Float_t eX;         // for generated vertexes - the real vertex position 
  Float_t eY;         // for reconstructed ones - average of track connection  
  Float_t eZ;         // points, used as local coordiantes origin (0,0,0) 
		      // to avoid accuracy problem

  Int_t	  eFlag;      // 0 - neutral (tracks starts attached only)
		      // 1 - charge (tracks ends&starts attached)
		      // 2 - back neutral (tracks ends attached only)
		      // 3 - neutral, linked (has common track with other vertex)
		      // 4 - charge, linked
		      // 5 - back neutral, linked
  Int_t   eMCEvt;
  Int_t	  eID;
  Float_t eQuality;   // Probability/(vsigmax**2+vsigmay**2)

  VERTEX::Vertex *eV; // pointer to VtVertex object
  
 public:
  EdbVertex();
  //EdbVertex(EdbVertex &v);
  virtual ~EdbVertex();

  //  int MakeV( bool usemom = true, bool usesegpar=false );
 
  void       Clear();
  void       ClearV();
  void       ClearNeighborhood();
  void       AddVTA(EdbVTA *vta);
  void       RemoveVTA(EdbVTA *vta);
  void       ResetTracks();
  void       Print();
  Int_t      Compare(const TObject *o) const;
  Bool_t     EstimateVertexMath( float& xv, float& yv, float& zv, float& d );
  void       Edb2Vt(const EdbTrackP& tr, VERTEX::Track& t, float X0 = 0., float m = 0.139 );
  void       Edb2Vt(const EdbSegP& s, VERTEX::Track& t,float X0 = 0., float m = 0.139 );
  Float_t    Chi2Track(EdbTrackP *tr, int zpos, float X0 = 0.);

  Float_t    MinDist();
  Float_t    Volume() { if(V()) return V()->vxerr()*V()->vyerr()*V()->vzerr(); else return 0; }
  Float_t    MaxAperture();
  Float_t    MaxImpact() { EdbVTA *vta=GetMaxImpVTA(); return vta? vta->Imp(): 0; }
  EdbVertex *GetConnectedVertex(int nv);
  EdbVertex *GetConnectedVertexForTrack(int it);
  Bool_t     IsEqual(const TObject *o) const;

  Int_t      N()          const  {return eVTa.GetSize();}
  Int_t      Nn()         const  {return eVTn.GetSize();}
  Int_t      Nv();
  Int_t      Flag()       const  {return eFlag;}
  Int_t      MCEvt()      const  {return eMCEvt;}
  Int_t      ID()         const  {return eID;}
  Int_t      Zpos(int i)         {return GetVTa(i)->Zpos();}
  ULong_t    Hash()       const  {return eID;}
  Bool_t     IsSortable() const  {return kTRUE;}
  Float_t    X()          const  {return eX;}
  Float_t    Y()          const  {return eY;}
  Float_t    Z()          const  {return eZ;}
  Float_t    VX()         const  {return eV ? (eV->vx() + eX) : 1000000.;}
  Float_t    VY()         const  {return eV ? (eV->vy() + eY) : 1000000.;}
  Float_t    VZ()         const  {return eV ? (eV->vz() + eZ) : 1000000.;}
  Float_t    Quality()           {return eQuality;}
  TList     *VTa()               {return &eVTa;}
  TList     *VTn()               {return &eVTn;} 
  EdbVTA    *GetVTa(int i)       {return (EdbVTA*)(eVTa.At(i));}
  EdbVTA    *GetVTn(int i)       {return (EdbVTA*)(eVTn.At(i));}
  EdbTrackP *GetTrack(int i)     {return GetVTa(i)->GetTrack();}
  EdbTrackP *GetTrackN(int i)     {return GetVTn(i)->GetTrack();}
  EdbSegP   *GetTrackV(int i, bool usesegpar=false);
  EdbVTA    *GetMaxImpVTA();

  float      CheckImpGeom(const EdbTrackP *tr);
  float      CheckImp(const EdbTrackP *tr);
  EdbVTA    *CheckImp(const EdbTrackP *tr, float ImpMax, int zpos, float dist);
  Float_t    Impact(int i) { return GetVTa(i)?  GetVTa(i)->Imp(): 1000000.; }
  Float_t    DistSeg(EdbSegP *seg, float X0 = 0.);
  Float_t    DistTrack(EdbTrackP *tr, int zpos, float X0 = 0.);
  Float_t    ImpTrack(int i) { return DistSeg( GetTrackV(i) ); }

  VERTEX::Vertex *V()     const  {return eV;}

  void       SetID(int ID = 0)                 {eID = ID;}
  void       SetXYZ(float x, float y, float z) {eX=x; eY=y; eZ=z;} 
  void       SetFlag(int flag = 0 )            {eFlag = flag;}
  void       SetMC(int mEvt=0) 		       {eMCEvt=mEvt;}
  void       SetV(VERTEX::Vertex *v)           {eV=v;}
  void       SetQuality( float q = 0 )         {eQuality = q;}

  Bool_t     TrackInVertex( EdbTrackP *t );
  Int_t      CheckDiscardedTracks();
  Int_t      EstimateVertexFlag();
  EdbTrackP *MeanTrack();

  ClassDef(EdbVertex,3)  // vertex class for OPERA emulsion data
};

//_________________________________________________________________________
class EdbVertexPar: public TObject {
 
 public:
  Float_t    eZbin;         // z- granularity (default is 100 microns)
  Float_t    eAbin;         // safety margin for angular aperture of vertex products
  Float_t    eDZmax;        // maximum z-gap in the track-vertex group
  Float_t    eProbMin;      // minimum acceptable probability for chi2-distance between tracks
  Float_t    eImpMax;       // maximal acceptable impact parameter (preliminary check)
  Float_t    eImpMaxV;      // if the impact is <= eImpMaxV the 2-vertex is accepted disregard to it's probability
  Bool_t     eUseMom;       // use or not track momentum for vertex calculations
  Bool_t     eUseSegPar;    // use only the nearest measured segments for vertex fit (as Neuchatel)
  Int_t      eQualityMode;  // vertex quality estimation method (0:=Prob/(sigVX^2+sigVY^2); 1:= inverse average track-vertex distance)
  Bool_t     eUseKalman;    // use or not Kalman for the vertex fit. Default is true
  Bool_t     eUseLimits;    // if true - look for the vertex only inside limits defined by eVmin:eVmax, default is false
  TVector3   eVmin,eVmax;   // limits for the vertex search
  
  EdbVertexPar();
  virtual ~EdbVertexPar(){}
  ClassDef(EdbVertexPar,1)  // vertex reconstruction parameters
};

//_________________________________________________________________________
class EdbVertexRec: public EdbVertexPar {

 private:

  EdbVertex *eVertex;
  EdbVertex *eWorking;

 public:

  TList      eVTA;          // vertex-track associations
  TObjArray *eEdbTracks;
  TObjArray *eVTX;          // array of vertex
  EdbPVRec  *ePVR;          // patterns volume (optional)

 public:
  EdbVertexRec() { Set0(); }
  EdbVertexRec( EdbVertexPar &vpar ) : EdbVertexPar(vpar) { Set0(); }
  virtual ~EdbVertexRec();
  
  void      Set0();
  void      Reset();
  const TVector3  *Vmin() const {return &eVmin;}
  const TVector3  *Vmax() const {return &eVmax;}
  int        RefitAll();
  void       AcceptPolish();
  void       RejectPolish();
  void       StatVertexN();
  void	     AcceptModifiedVTX(EdbVertex *eVertex, EdbVertex *eWorking);
  void       CancelModifiedVTX(EdbVertex *eVertex, EdbVertex *eWorking);
  void       FillTracksStartEnd(TIndexCell &starts, TIndexCell &ends);
  
  Int_t      MakeV(EdbVertex &edbv, bool isRefit=false);
  Int_t      FindVertex();

  EdbVertex *Make1Vertex(TObjArray &tracks, float zexpected);
  EdbVertex *StripBadTracks( EdbVertex &v,  float impMax, int ntrMin );

  EdbVertex *ProbVertex2(EdbTrackP *tr1, EdbTrackP *tr2, int zpos1, int zpos2 );
  Int_t      ProbVertexN();
  Int_t      ProbVertexN_old();
  Int_t      ProbVertexNpos(int zpos);
  void       CheckVTX();
  EdbVertex *TestVTAGroup(TObjArray &arrvta);
  int        EstimateVertexFlag(int zpos1, int zpos2);
  
  Bool_t     CheckDZ2(float z1, float z2, int zpos1, int zpos2, float z );
  Bool_t     IsInsideLimits(EdbSegP &s);
  Int_t      FindSimilarTracks(EdbTrackP &t, TObjArray &found, int nsegmin=2, float dMax=100., float dTheta=0.01, float dZmax=50000.);
  Int_t      FindSimilarTracksE(EdbSegP &spred, TObjArray &found, bool startend,
				float impact, float dthetaMax, float dxy, float zminT, float zmaxT, float zminV, float zmaxV);
  bool       CompatibleSegments(EdbSegP &spred, EdbSegP &stest,
				float impact, float dthetaMax, float dxy, float zminT, float zmaxT, float zminV, float zmaxV);
  Int_t      FindSimilarSegments(EdbSegP &spred, TObjArray &found, EdbPattern &pat,
				float impact, float dthetaMax, float dxy, float zminT, float zmaxT, float zminV, float zmaxV);

  Int_t	     LinkedVertexes();
  Int_t      LoopVertex(TIndexCell &list1, TIndexCell &list2, int zpos1, 
			int zpos2);
  Int_t      AddSegmentToVertex(EdbSegP *s, float ImpMax = 25., 
				float ProbMin = 0.01, float Mom = 0.3);
  Int_t      VertexPolish(EdbVertex *v, int refill = 0, float RadMax = 1000., 
			  int Dpat = 2, float ImpMax = 25., 
			  float ProbMin = 0.01, float Mom = 0.3);
  Int_t	     VertexTuning(int seltype = 0);
  Int_t	     VertexNeighbor(float RadMax = 1000., int Dpat = 1, 
			    float ImpMax = 1000000.);
  Int_t	     VertexNeighbor(EdbVertex *v, float RadMax = 1000., int Dpat = 1, 
			    float ImpMax = 1000000.);
  Int_t	     SelVertNeighbor(EdbVertex *v, int seltype, float RadMax, 
			     int Dpat, TObjArray *ao);
  Int_t	     SelSegNeighbor(EdbSegP *s, int seltype, float RadMax, int Dpat, TObjArray *ao);
  Int_t	     SegmentNeighbor(EdbSegP *s, float RadMax = 1000., int Dpat = 1, float ImpMax = 1000000., float SegWmin=9, 
                             TObjArray *aseg = 0, TObjArray *atr = 0, TObjArray *arv = 0);
  Float_t    CheckImpact( EdbSegP *s1,   EdbSegP *s2, int zpos1, int zpos2, float pv[3] );
  Float_t    CheckImpactN( EdbSegP *s1,   EdbSegP *s2, float pv[3], bool &parallel, float dzMax=6000. );
  Bool_t     EstimateVertexQuality(EdbVertex &v);
  Bool_t     EstimateVertexPosition(EdbVertex &v);
  Double_t   Tdistance(const VERTEX::Track& t1, const VERTEX::Track& t2);
  Double_t   Tdistance(const EdbSegP& s1, const EdbSegP& s2);
  Double_t   TdistanceChi2(const EdbTrackP& tr1, const EdbTrackP& tr2);
  Double_t   TdistanceChi2(const EdbSegP& s1, const EdbSegP& s2, float m);
  Double_t   MoveTrackToOtherVertex(EdbVertex *v2, int it2max, EdbVertex *v1, 
				    int seltype, EdbVertex **v2n, 
				    EdbVertex **v1n);
  EdbVertex *AddTrackToVertex(EdbVertex *eVertex, EdbTrackP *eTr, int zpos);
  EdbVertex *RemoveTrackFromVertex(EdbVertex *eVertex, int itr);
  EdbVertex *RemoveVTAFromVertex( EdbVertex &vtx, EdbVTA &vta );

  EdbTrackP *GetEdbTrack(const int index);
  EdbVTA    *AddTrack(EdbVertex &edbv, EdbTrackP *track, int zpos);

  void       SetPVRec(EdbPVRec *pvr) {ePVR = pvr;}
  void       AddVTA(EdbVTA *vta)     {eVTA.Add((TObject*)vta);}
  Int_t      Nvtx()           const  {return eVTX ? eVTX->GetEntries() : 0;}
  EdbVertex *GetVertex(Int_t &i)     {return eVTX ? (EdbVertex*)eVTX->At(i):0;}

  void AddVertex(EdbVertex *vtx) {
    if (!eVTX) eVTX = new TObjArray();
    eVTX->Add((TObject*)vtx);
  }

  int CheckTrack(EdbTrackP &track, int zpos);


  ClassDef(EdbVertexRec,3) //reconstruct vertexes in OPERA emulsion data
};

#endif /* ROOT_EdbVertex */
