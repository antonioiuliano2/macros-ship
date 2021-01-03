#ifndef ROOT_EdbDataSet
#define ROOT_EdbDataSet

//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbDataSet                                                           //
//                                                                      //
// OPERA data set definition&reconstruction                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////
#include "TObjArray.h"
#include "TFile.h"
#include "TCut.h"
#include "EdbMask.h"
#include "EdbAffine.h"
#include "EdbRun.h"
#include "EdbPVRec.h"
#include "EdbLayer.h"
#include "EdbSegmentCut.h"
#include "EdbVertex.h"
#include "TMatrix.h"
#include "TArrayC.h"

//______________________________________________________________________________
class EdbDataPiece : public TNamed {

 public:
  Int_t        eAFID;           // 1-use fiducial marks transformations, 0 - do not
  Int_t        eCLUST;          // 1-use clusters, 0 - do not

 private:
  Int_t        ePlate;          // plate id
  Int_t        ePiece;          // piece id in this plate
  Int_t        eFlag;           // 0-do nothing, 1-do something
  TObjArray    eRunFiles;       //
  //  TString      eFileNameRaw;    // name of the raw data file (run)

  EdbLayer    *eLayers[3];      // base(0),up(1),down(2)  layers
  EdbScanCond *eCond[3];        //
  TIndexCell  *eAreas[3];       // base/up/down  surface areas list
  TObjArray   *eCuts[3];        // array of cuts
  Float_t      eCutCP[6];       //
  Float_t      eCutGR;          // grain cut (chi)
  Int_t        eOUTPUT;         //

  TCut        *eRCuts[3];       //! root-style text cuts

 public: 
  TString      eFileNameCP;     // name of the couples data file
  TString      eFileNamePar;    // name of the parameters file
  TIndexCell  *eCouplesInd;     //

  EdbRun      *eRun;            //!
  TTree       *eCouplesTree;    //!

  EdbMask     *eEraseMask;      // id's (entries) of segments to be erased when read couples tree

 public:
  EdbDataPiece();
  EdbDataPiece(int plate, int piece, char* file, int flag);
  virtual ~EdbDataPiece();

  void Set0();

  int  Plate() const {return ePlate;}
  int  InitCouplesInd();
  int  GetLinkedSegEntr(int side, int aid, int vid, int sid, TArrayI &entr) const;
  void SetVolume0( float x0, float y0, float z0, float tx=0, float ty=0 );
  void SetVolumeA( float dx, float dy ) { GetLayer(0)->SetDXDY( dx, dy); }
  void AddRunFile( const char *name );
  void CloseRun();
  const char *GetRunFile(int i) const;
  const char *MakeName();
  const char *MakeNameCP(const char *dir);
  const char *GetNameCP() const {return eFileNameCP.Data();}
  void  Print();
  void  WriteCuts();
  int   CheckCCD(int maxentr=2000);
  int   RemoveCCDPeak(TMatrix &matr);
  int   UpdateSegmentCut(EdbSegmentCut cut);

  void SetCouplesTree(TTree *tree) {eCouplesTree=tree;}
  int       Nruns() const { return eRunFiles.GetEntriesFast(); }
  int       Flag() const {return eFlag;}
  EdbLayer *GetMakeLayer(int id);
  EdbLayer *GetLayer(int id)
    { if(eLayers[id]) return (EdbLayer *)eLayers[id]; else return 0; }
  EdbScanCond *GetMakeCond(int id);
  EdbScanCond *GetCond(int id)
    { if(eCond[id]) return (EdbScanCond *)eCond[id]; else return 0; }


  void           SetOUTPUT(int out=1) {eOUTPUT=out;}
  void           SetCutGR(float chi) {eCutGR=chi;}
  void           AddCutCP(float var[6]);
  void           AddSegmentCut(int layer, int xi, float var[10]);
  void           AddSegmentCut(int layer, int xi, float min[5], float max[5]);
  int            NCuts(int layer);
  EdbSegmentCut *GetCut(int layer, int i)
    { return (EdbSegmentCut *)(eCuts[layer]->UncheckedAt(i)); }

  void           AddRCut(int layer, TCut &cut);
  TCut           *GetRCut(int layer)  { return eRCuts[layer]; }

  float           GetCutGR() const {return eCutGR;}
  int             GetOUTPUT() const {return eOUTPUT;}

  int  AcceptViewHeader(const EdbViewHeader *head);
  void MakeNamePar(const char *dir);
  int  CorrectAngles();
  int  CorrectAngles(TTree *tree);
  void CorrectShrinkage( int layer, float shr );
  int  UpdateShrPar( int layer );
  int  UpdateAffPar(  int layer, EdbAffine2D &aff);
  int  UpdateAffTPar( int layer, EdbAffine2D &aff);
  int  UpdateZPar( int layer, float z );
  int  TakePiecePar();
  int  ReadPiecePar(const char *file);
  int  MakeLinkListArea(int irun);
  int  MakeLinkListCoord(int irun);
  int  GetAreaData(EdbPVRec *ali, int const area, int const side);
  int  TakeRawSegment(EdbView *view, int id, EdbSegP &segP, int side);
  int  PassCuts(int id, float var[5]);
  int  PassCutCP(float var[6]);

  float GetRawSegmentPix( EdbSegment *seg );
  float CalculateSegmentChi2( EdbSegment *seg, float sx, float sy, float sz );

  int  GetRawData(EdbPVRec *ali);
  int  GetCPData( EdbPattern *pat, EdbPattern *p1=0, EdbPattern *p2=0 );
  int  GetCPData_new( EdbPattern *pat, EdbPattern *p1=0, EdbPattern *p2=0, TIndex2 *trseg=0 );
  int  TakeCPSegment(EdbSegCouple &cp, EdbSegP &segP);

  int   InitCouplesTree( const char *mode="READ" );
  static TTree *InitCouplesTree( const char *file, const char *mode );
  void CloseCPData();

  ClassDef(EdbDataPiece,2)  // Edb raw data unit (scanned plate) associated with run file
};


//______________________________________________________________________________
class EdbDataSet : public TNamed {

 private:

  TString    eInputList;    // list of input data (runs)
  TString    eAnaDir;       // path for analysis data directory
  TString    eParDir;       // path for parameters directory
  TString    eDBFileName;   // root file to keep pieces parameters
  TFile      *eDBFile;      // the file (database) to save all parameters

  TObjArray  ePieces;       // array of runs

 public:
  EdbDataSet();
  EdbDataSet(const char *file);
  virtual ~EdbDataSet();

  void Set0();
  int N() const { return ePieces.GetEntriesFast(); }
  EdbDataPiece *GetPiece(int id) 
    { if(id<ePieces.GetEntriesFast()) return (EdbDataPiece *)ePieces.At(id); else return 0;}

  const char *GetAnaDir() const {return eAnaDir.Data();}
  const char *GetParDir() const {return eParDir.Data();}
  int  ReadDataSetDef(const char *file);
  int  GetRunList(const char *file);
  void PrintRunList();
  void WriteRunList();

  EdbDataPiece *FindPiece(const char *name);

  void Print();

  ClassDef(EdbDataSet,1)  // OPERA emulsion data set
};

//______________________________________________________________________________
class EdbDataProc : public TObject {

 private:

  EdbDataSet *eDataSet;

  EdbPVRec   *ePVR;

  int  eNoUpdate;

 public:
  EdbDataProc();
  EdbDataProc(int npl, TArrayI &ids, TArrayF &zs);
  EdbDataProc(const char *file);
  virtual ~EdbDataProc();

  EdbDataSet *GetDataSet() {return eDataSet;}
  EdbPVRec *PVR() const {return ePVR;}
  EdbPVRec* GetPVR() const {return ePVR;}
  void   SetPVR(EdbPVRec *pvr)  {ePVR=pvr;}

  EdbPVRec *ExtractDataVolume( EdbSegP &v, int plmin, int plmax,
			       float accept[4], 
			       int datatype=0 );
  EdbPVRec *ExtractDataVolume( EdbTrackP &tr, float binx=20, float bint=10,
			       int datatype=0 );
  EdbPVRec *ExtractDataVolumeF( EdbTrackP &tr, float binx=20, float bint=10,
			       int datatype=0 );

  int  InitVolume(int datatype=0, const char *rcut="1");
  int  InitVolume(EdbPVRec *ali, int datatype=0, TIndex2 *trseg=0 );
  int  InitVolumeTracks(EdbPVRec *ali, const char *rcut);
  int  InitVolumeRaw(EdbPVRec *ali);
  int  Process(){ return Link(); }  // to be removed
  int  CheckCCD();
  int  Link();
  int  Link(EdbDataPiece &piece);
  void Align(int alignFlag);

  static int LinkTracksWithFlag( EdbPVRec *ali, float p, float probmin, int nsegmin, int maxgap, int flag, float mass=0.1396 );
  void LinkTracks(int alg=0, float p=-1.);
  void LinkTracksC(int alg=0, float p=-1.);
  void LinkRawTracks(int alg=0);
  void AlignLinkTracks(int alg=0, int alignFlag=0);

  void SetNoUpdate(int nu) { eNoUpdate=nu; }
  int  NoUpdate() const    { return eNoUpdate; }

  int    ShrinkCorr() {return 1;}
  int    CheckShrinkage( EdbPVRec *ali, int couple, float &shr1, float &shr2 );
  void   CorrectAngles();

  void   AjustZ(int doZ);
  void   FineAlignment(int doFine);
  void   FineAlignmentTracks();

  void   FillCouplesTree( TTree *tree, EdbPVRec *al, int fillraw=0 );
  void   CloseCouplesTree( TTree *tree );

  static int MakeTracksTree(EdbPVRec *ali=0, const char *file="linked_tracks.root");
  static int MakeVertexTree(TObjArray &vtxarr, const char *file);
  static int MakeTracksTree(TObjArray &tracks, float xv=0, float yv=0, const char *file="linked_tracks.root");
  static int ReadVertexTree( EdbVertexRec &vertexrec, const char     *fname, const char *rcut, map<int,EdbTrackP*> trackID_map);
  static int ReadVertexTree( EdbVertexRec &vertexrec, const char     *fname, const char *rcut, TObjArray *builttracks = 0);
  EdbVertex* GetVertexFromTree( EdbVertexRec &vertexrec, const char     *fname, const int vertexID );
  static int ReadTracksTree(EdbPVRec &ali,
			    const char *fname="linked_tracks.root",
			    //			    int   nsegMin=3,
			    //			    float probMin=0.01, 
			    const char *rcut="t.eFlag>-1&&nseg>2&&t.eProb>.01" );

  TIndex2 *MakeTracksSegmentsList( EdbPVRec &ali );

  ClassDef(EdbDataProc,1)  // emulsion data processing
};

#endif /* ROOT_EdbDataSet */
