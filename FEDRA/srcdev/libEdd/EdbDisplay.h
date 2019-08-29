#ifndef ROOT_EdbDisplay
#define ROOT_EdbDisplay

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDisplay                                                           //
//                                                                      //
// Class to display pattern volume in 3D                                //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TGNumberEntry.h"
#include "EdbDisplayBase.h"
#include "EdbPVRec.h"
#include "EdbVertex.h"
#include "EdbBrick.h"

class EdbSegG;
//class EdbBrickP;

//_________________________________________________________________________
class EdbDisplay: public EdbDisplayBase {

 private:

  Int_t      eDrawTracks;  // tracks drawing option

  Int_t      eDrawVertex;  // vertex drawing option

  TArrayI    *eColors;
  TArrayF    *eDZs;

  TGNumberEntry *fNumericEntries[3];

 public:

  Int_t        eFromPlate; // plates range (for colors normalization)
  Int_t        eToPlate;

  //EdbBrickP    *eB;        // brick object with the plates geometry inside
  EdbVertexRec *eVerRec;
  TObjArray *eArrSegP;     // array of segments to be drawn
  TObjArray *eArrTr;       // array of tracks to be drawn
  TObjArray *eArrV;        // array of vertices to be drawn
  TObjArray *eArrSegPSave;     // saved array of segments to be drawn
  TObjArray *eArrTrSave;   // saved array of tracks to be drawn
  TObjArray *eArrVSave;    // saved array of vertices to be drawn
  EdbVertex *eWorking;     // working vertex
  EdbVertex *eVertex;      // current selected vertex
  EdbVertex *ePrevious;    // saved previous vertex modifications
  EdbSegP   *eSegment;     // working segment (for segment neighborhood)
  EdbTrackP *eTrack;       // working intermediate track (track creation)
  EdbTrackP *eTrack1;      // working intermediate track (track splitting)
  EdbTrackP *eTrack2;      // working intermediate track (track splitting)
  TPolyMarker3D *eSegPM;   // green mark for segment selected as working
  Bool_t eWait_Answer;	   // set TRUE when answer received
  Int_t eIndVert;	   // Index of selected vertex in ArrV
  Int_t eIndVertSave;	   // Index of selected vertex in ArrV (seved)
  TList eCreatedTracks;    // list of tracks, created during vertex operations 
  Double_t eRadMax;        // Maximal Radius for neighborhood
  Double_t eDpat;          // +/- patterns for neighborhood
  Double_t eImpMax;        // Maximal impact for neighborhood
  Double_t eP;             // track momentum (creation from segment, propagation)
  Double_t eM;             // track mass (creation from segment, propagation)
  Double_t eTImpMax;       // Maximal impact for interactive add track
  Double_t eTProbMin;      // Minimal probability for interactive add track
  Double_t eSegWmin;       // Minimal segment W for neighbouring selection

  TObjArray *eArrSegG;     // additional array of segments for the presentation purpose only

 public:

  EdbDisplay() : EdbDisplayBase() {Set0();}
  ~EdbDisplay();

  EdbDisplay(const char *title,
	     Float_t x0, Float_t x1,
	     Float_t y0, Float_t y1,
	     Float_t z0, Float_t z1,
	     TCanvas *Canvas = 0) : 
    EdbDisplayBase(title, x0, x1, y0, y1, z0, z1, Canvas) {Set0();}

    EdbDisplay(const char *title, EdbLayer &la, TCanvas *Canvas = 0) : 
      EdbDisplayBase(title, la.Xmin(), la.Xmax(), la.Ymin(), la.Ymax(), 
		     la.Zmin(), la.Zmax(), Canvas)
      {Set0();}

  static EdbDisplay *EdbDisplayExist(const char *title);
  void Delete();
  void Set0();
  void GuessRange( float margZmin=3000,float margZmax=1000,float margR=300 );
  void SetVerRec(EdbVertexRec *evr) { eVerRec = evr; };

  void Refresh();
  void SetArrSegG(TObjArray *arrg) {eArrSegG = arrg;}
  void SetArrSegP(TObjArray *arr);
  void SetArrTr(TObjArray *arr);
  void SetDrawTracks(int opt) {eDrawTracks=opt;}

  EdbVertexRec *VerRec() const {return eVerRec;}
  //void PatternDraw(EdbPattern &pat);
  void TrackDraw(EdbTrackP *tr, Color_t kColor=kWhite);
  EdbSegG *SegLine(const EdbSegP *seg);

  void DrawRef(float start[3], float end[3]);
  void SelectVertexTracks(TObjArray *vtx);
  void SetArrV(TObjArray *arrv);
  void VertexDraw(EdbVertex *v);
  void SetDrawVertex(int opt) {eDrawVertex=opt;}
  void CancelModifiedVTX();
  void DeleteModifiedVTX();
  void AcceptModifiedVTX();
  void DialogModifiedVTX();
  void CloseDialogModifiedVTX();
  void UndoModifiedVTX();
  void DrawVertexEnvironment();
  void DrawAllObjects();
  void DrawVTXTracks(char *type, EdbVertex *v = 0);
  void RemoveTrackFromTable( int ivt = 0 );
  void DialogNeighborParameters();
  void AcceptModifiedParams();
  void CloseDialogModifiedParams();
  void CancelDialogModifiedParams();
  void DialogTrackParameters();
  void AcceptModifiedTrackParams();
  void CloseDialogModifiedTrackParams();
  void CancelDialogModifiedTrackParams();
  void ClearSegmentEnv();
  void DrawSegmentExtrapolationLine(const EdbSegP &s, float zmin, float zmax);
 
  ClassDef(EdbDisplay,1) //FEDRA Event Display
};
//_________________________________________________________________________
class EdbVertexG : public TPolyMarker3D {
 private:

  EdbVertex *eV;
  EdbDisplay *eD;

 public:
  EdbVertexG():TPolyMarker3D(1) {eV=0; eD=0;}
  EdbVertexG(EdbDisplay *D):TPolyMarker3D(1) {eV=0; eD=D;}
  virtual ~EdbVertexG(){}

  void SetVertex(EdbVertex *v) {eV=v;}

  virtual void          DumpVertex();    // *MENU*
  virtual void          InspectVertex(); // *MENU*
  virtual void		SetAsWorking();  // *MENU*
  virtual void          DeleteVertex();  // *MENU*
  virtual void          RemoveKink();    // *MENU*
  virtual void          TestVertex();    // *MENU*
  virtual const char *	GetTitle() const;
  virtual const char *	GetName() const;
  virtual char *	GetObjectInfo(int px, int py) const;

 ClassDef(EdbVertexG,1)  //Vertex
};

//_________________________________________________________________________
class EdbTrackG : public TPolyMarker3D {
 private:

  EdbTrackP *eTr;
  EdbDisplay *eD;
  
 public:
  EdbTrackG() {eTr=0; eD=0;}
  EdbTrackG(EdbDisplay *D) {eTr=0; eD=D;}
  EdbTrackG(Int_t nhits, EdbDisplay *D):TPolyMarker3D(nhits) {eTr=0; eD=D;}
  virtual ~EdbTrackG(){}

  void SetTrack(EdbTrackP *tr) {eTr=tr;}

  virtual void          DumpTrack();       // *MENU*
  virtual void          InspectTrack();    // *MENU*
  virtual void		SetAsWorkingVertex();        // *MENU*
  virtual void		RemoveTrackFromVertex();     // *MENU*
  virtual void		AddTrackToVertex();          // *MENU*
  virtual void		FixNewTrack();     // *MENU*
  virtual void		DeleteTrack();     // *MENU*
  virtual void		UndoNewTrack();    // *MENU*
  virtual void		UndoSplit();       // *MENU*
  virtual void		UndoRemoveKink();  // *MENU*
  virtual void		AddToNewTrack();      // *MENU*
  virtual void		AddToNewTrackAndFit();// *MENU*
  virtual void		InfoTrackVert();      // *MENU*
  virtual void		EstimateMomentum();   // *MENU*
  virtual const char *	GetTitle() const;
  virtual const char *	GetName() const;
  virtual char *	GetObjectInfo(int px, int py) const;

 ClassDef(EdbTrackG,1)  //Track
};

//_________________________________________________________________________
class EdbSegG : public TPolyLine3D {
 private:

  const EdbSegP *eSeg;
  EdbDisplay *eD;

 public:
  EdbSegG() {eSeg=0; eD=0;}
  EdbSegG(EdbSegP &s);
  EdbSegG(EdbDisplay *D) {eSeg=0; eD=D;}
  EdbSegG(Int_t nhits):TPolyLine3D(nhits) {eSeg=0; eD=0;}
  EdbSegG(Int_t nhits, EdbDisplay *D):TPolyLine3D(nhits) {eSeg=0; eD=D;}
  virtual ~EdbSegG(){}

  void SetSeg(const EdbSegP *s) {eSeg=s;}
  float X() { return GetP()[0]; }
  float Y() { return GetP()[1]; }
  float Z() { return GetP()[2]; }

  virtual void          DumpSegment(); 	   // *MENU*
  virtual void          InspectSegment();  // *MENU*
  virtual void		AddAsTrackToVertex();      // *MENU*
  virtual void		AddToNewTrack();      // *MENU*
  virtual void		AddToNewTrackAndFit();// *MENU*
  virtual void		RemoveFromTrack(); // *MENU*
  virtual void		SplitTrack();      // *MENU*
  virtual void		SetAsWorking();    // *MENU*
  virtual void		InfoSegVert();     // *MENU*
  virtual void		InfoSegSeg();      // *MENU*
  virtual const char *	GetTitle() const;
  virtual const char *	GetName() const;
  virtual char *	GetObjectInfo(int px, int py) const;

 ClassDef(EdbSegG,1) //Segment
};


#endif /* ROOT_EdbDisplay */
