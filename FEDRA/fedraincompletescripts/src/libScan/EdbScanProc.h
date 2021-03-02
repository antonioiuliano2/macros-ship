#ifndef ROOT_EdbScanProc
#define ROOT_EdbScanProc

#include "EdbRunAccess.h"
#include "EdbDataSet.h"
#include "EdbScanClient.h"
#include "EdbScanSet.h"

class TEnv;

class EdbScanProc : public TNamed
{
 public:
 TString eProcDirClient;    // directory path for root data
 TString eProcDirServer;    // directory path for root data
  TString eParDir;           // directory path for off-line processing parameters

private:
  TString eServerCreatedRunName; //EdbRun file name that is created by scanserver side

public: 
  EdbScanProc();
  virtual ~EdbScanProc(){}

  char   *BrickDir(int brick);
  bool    CheckDir(const char *dir, bool create=true);
  bool    CheckDirWritable(const char *dir);
  bool    CheckAFFDir(int brick, bool create=true);
  bool    CheckBrickDir(EdbID id, bool create=true);
  bool    CheckPlateDir(EdbID id, bool create=true);
  bool    CheckProcDir(int id[4], bool create=true);
  bool    CheckProcDir( EdbID id, bool create=true) {int id4[4]; id.Get(id4); return CheckProcDir(id4,create); }
  void    MakeFileName(TString &s, int id[4], const char *suffix, bool inplate=true);
  void    MakeFileName(TString &s, EdbID id, const char *suffix, bool inplate=true) 
    {int id4[4]; id.Get(id4); return MakeFileName(s,id4,suffix,inplate);}
  void    MakeFileNameSrv(TString &s, int id[4], const char *suffix, bool inplate=true);
  void    MakeFileNameSrv(TString &s, EdbID id, const char *suffix, bool inplate=true) 
    {int id4[4]; id.Get(id4); return MakeFileNameSrv(s,id4,suffix,inplate);}
  void    MakeAffName(TString &s, int id1[4], int id2[4], const char *suffix="aff.par");
  void    MakeAffName(TString &s, EdbID id1, EdbID id2, const char *suffix="aff.par") 
    { int id14[4]; id1.Get(id14); int id24[4]; id2.Get(id24); return MakeAffName(s,id14,id24,suffix); }
  bool    GetMap(int brick, TString &map);
  bool    AddParLine(const char *file, const char *line, bool recreate=false );
  bool    MakeInPar(int id[4], const char *option);
  bool    MakeInPar(EdbID id, const char *option)  {int id4[4]; id.Get(id4); return MakeInPar(id4,option);}
  void    MakeInParSet(EdbID id, const char *option);
  int     CopyFile(int id1[4], int id2[4], const char *suffix, bool overwrite);
  int     CopyPar(EdbID id1, EdbID id2, bool overwrite=true) {
    int id14[4], id24[4]; id1.Get(id14); id2.Get(id24); 
    return CopyPar(id14,id24,overwrite);}
  int     CopyPar(int id1[4], int id2[4], bool overwrite=true) {return CopyFile(id1,id2,"par",overwrite);}
  int     CopyPred(int id1[4],int id2[4], bool overwrite=true) {return CopyFile(id1,id2,"pred.root",overwrite);}
  int     CopyAFFPar(int id1c[4],int id2c[4], int id1p[4], int id2p[4], bool overwrite=true);
  int     RemoveFile(EdbID id, const char *suffix);
  void    CopyParSet(EdbID idset1, EdbID idset2);
  bool    ReadPiecePar(EdbID id, EdbPlateP &plate);

  void    MakeScannedIDList( EdbID id0, EdbScanSet &sc, int pl_from, int pl_to, const char *suffix );
  void    CheckFiles( EdbScanSet &sc, const char *suffix );

  int     ReadPatTXT(EdbPattern &pred, EdbID id, const char *suffix, int flag=-1) {int id4[4]; id.Get(id4); return ReadPatTXT(pred,id4,suffix,flag);}
  int     ReadPatTXT(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     ReadPatTXT(const char *file, EdbPattern &pred, int flag=-1);
  int     WritePatTXT(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     WriteSBcndTXT(int id[4], const char *suffix="man.sbt.txt");
  int     ReadPatRoot(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     WritePatRoot(EdbPattern &pred, int id[4], const char *suffix, int flag=-1);
  int     ReadPred(EdbPattern &pred, int id[4], int flag=-1) {return ReadPatRoot(pred,id,"pred.root",flag);}
  int     WritePred(EdbPattern &pred, int id[4], int flag=-1) {return WritePatRoot(pred,id,"pred.root",flag);}
  int     ReadFound(EdbPattern &pred, int id[4], int flag=-1) {return ReadPatRoot(pred,id,"found.root",flag);}
  int     WriteFound(EdbPattern &pred, int id[4], int flag=-1) {return WritePatRoot(pred,id,"found.root",flag);}
  int     ReadFound(EdbPattern &pred, EdbID id, int flag=-1) {int id4[4]; id.Get(id4); return ReadFound(pred,id4,flag);}
  int     ReadPred(EdbPattern &pred, EdbID id, int flag=-1) {int id4[4]; id.Get(id4); return ReadPred(pred,id4,flag);}
  int     WritePred(EdbPattern &pred, EdbID id, int flag=-1) {int id4[4]; id.Get(id4); return WritePred(pred,id4,flag);}
  int     WriteFound(EdbPattern &found, EdbID id, int flag=-1) {int id4[4]; id.Get(id4); return WriteFound(found,id4,flag);}
  int     WritePatTXT(EdbPattern &pred, EdbID id, const char *suffix, int flag=-1) {int id4[4]; id.Get(id4); return WritePatTXT(pred,id4,suffix,flag);}

  bool    WaitFileReady(const char* fname_); //waits file copied/moved, in ready state
  EdbRun *InitRun(int id[4], char* runname_ = NULL, char* runnamesrv_ = NULL, bool createrun_=true);
  bool    FlashRawDir(EdbScanClient &scan, int id[4]);
  int     LoadPlate(EdbScanClient &scan, int id[4], int attempts=1);
  int     ScanAreas(EdbScanClient::ScanType st, EdbScanClient &scan, int id[4], int flag=-1, const char *opt="NOCLCLFRAMESUM");
  int     ScanAreas(EdbScanClient::ScanType st, EdbScanClient &scan, EdbPattern &pred, int id[4], const char *opt="NOCLCLFRAMESUM"); // NEW!!!

  bool    InitPiece(EdbDataPiece &piece, int id[4]);
  bool    InitPiece(EdbDataPiece &piece, EdbID id) {int id4[4]; id.Get(id4); return InitPiece(piece,id4);}
  int     ReadPiece(EdbDataPiece &piece, EdbPattern &pat);
  int     ReadPatCP(EdbPattern &pat, int id[4], TCut c="1");
  int     ReadPatCP(EdbPattern &pat, EdbID id, TCut c="1") {int id4[4]; id.Get(id4); return ReadPatCP(pat,id4,c);}
  int     ReadPatCPnopar(EdbPattern &pat, EdbID id, TCut cut="1", bool do_erase=false, bool read_mt=false);
  int     ReadPatCPnopar(EdbPattern &pat, const char *file, TCut cut="1", EdbMask *erase_mask=0, bool read_mt=false);
  EdbMask* ReadEraseMask(EdbID id);
  void    MakeEraseFile(EdbID id, EdbPattern &pat);
  bool    ApplyAffZ(EdbPattern &pat,int id1[4],int id2[4]);
  bool    GetAffZ(EdbAffine2D &aff, float &z,int id1[4],int id2[4]);
  bool    SetAFFDZ(int id1[4], int id2[4], float dz);
  bool    SetAFF0(int id1[4], int id2[4]);
  bool    MakeAFFSet(EdbScanSet &sc);
  bool    MakeParSet(EdbScanSet &sc);
  bool    PrepareSetStructure(EdbScanSet &sc);

  int     ConvertAreas(EdbScanClient &scan, int id[4], int flag=-1, const char *opt="NOCLCLFRAMESUM");
  int     CorrectAngles(int id[4]);
  int     LinkRun(int id[4], int noUpdate=1);
  int     LinkRunAll(int id[4], int npre=3, int nfull=1, int correct_ang=1);
  int     LinkRunAll(EdbID id, int npre=3, int nfull=1, int correct_ang=1)  {int id4[4]; id.Get(id4); return LinkRunAll(id4,npre,nfull,correct_ang);}
  int     LinkSet(EdbScanSet &sc, int npre=3, int nfull=1, int correct_ang=1);

  void     GetPatternSide( EdbID id, int side, EdbLayer &la, const char *segcut, int afid, EdbPattern &p);
  void     LinkRunTest(EdbID id, EdbPlateP &plate, TEnv &cenv);
  void     LinkRunNew(EdbID id, EdbPlateP &plate, TEnv &cenv);
  void     LinkSetNew(EdbScanSet &sc, TEnv &cenv);
  void     LinkSetNewTest(EdbScanSet &sc, TEnv &cenv);

  //void     WriteSetGeom(EdbScanSet &sc);
  //void     ReadSetGeom(EdbScanSet &sc);

  int     AlignNewNopar(EdbID id1, EdbID id2, TEnv &cenv, EdbAffine2D *aff=0, float dz=0);
  int     AlignNewLayernopar(EdbID id1, EdbID id2, TEnv &cenv, EdbPlateP &plate2, EdbAffine2D *aff=0, float dz=0);
  
  bool    UpdateAFFPar( EdbID id1, EdbID id2, EdbLayer &l, EdbAffine2D *aff0=0);
  bool    UpdatePlatePar( EdbID id, EdbLayer &l);
  int     AlignSetNewNopar(EdbScanSet &sc, TEnv &cenv);
  void    AlignSetNewNopar(EdbID id, TEnv &cenv);

  void    AlignSet(EdbID id, int npre, int nfull, const char *opt="-z"  );
  int     AlignSet( EdbScanSet &sc, int npre=1, int nfull=3, const char *opt="-z");
  int     Align(EdbID id1, EdbID id2, const char *option="") {int id14[4]; id1.Get(id14); int id24[4]; id2.Get(id24); return Align(id14,id24,option);}
  int     Align(int id1[4], int id2[4], const char *option="");
  int     AlignAll(int id1[4], int id2[4], int npre=1, int nfull=3, const char *opt="-z");
  int     AlignAll(EdbID id1, EdbID id2, int npre=1, int nfull=3, const char *opt="-z")
    {int id41[4]; id1.Get(id41); int id42[4]; id2.Get(id42); return AlignAll(id41, id42, npre, nfull, opt);}
  
  int     TrackSetBT( EdbScanSet &sc, TEnv &cenv);
  int     ReadTracksTree(EdbID id, EdbPVRec &ali, TCut cut="1");
  int     ReadTracksTree(const char *name, EdbPVRec &ali, TCut cut="1");

  bool    CorrectPredWithFound(int id1[4], int id2[4], const char *opt="-z", int patmin=6);
  bool    CorrectAffWithPred(int id1[4], int id2[4], const char *opt="-z", int patmin=6, const char *parfile="fullalignment");
  bool    ProjectFound(int id1[4],int id2[4]);
  bool    ProjectFound(EdbID id1,EdbID id2) {int id14[4]; id1.Get(id14); int id24[4]; id2.Get(id24); return ProjectFound(id14, id24); }

  int     FindPredictions(EdbPattern &pred, int id[4], EdbPattern &found, int maxholes=3);
  int     FindPredictions(int id[4], int flag=-1, int maxholes=3);

  bool    InitRunAccessNew(EdbRunAccess &ra, EdbID id, EdbPlateP &plate, bool do_update=false);
  bool    InitRunAccessNew(EdbRunAccess &ra, EdbID idset, int idplate, bool do_update=false);
    
  bool    InitRunAccess(EdbRunAccess &ra, int id[4], bool do_update=false);
  bool    InitRunAccess(EdbRunAccess &ra, EdbID id, bool do_update=false) {int id4[4]; id.Get(id4); return InitRunAccess(ra, id4, do_update); }

  int     FindPredictionsRawSet(EdbID idp, EdbScanSet &ss, int npl);
  int     FindPredictionsRaw(EdbID idp, EdbID idr);
  int     FindPredictionsRaw(EdbPattern &pred, EdbPattern &found, EdbRunAccess &ra, 
			     EdbScanCond &condBT, EdbScanCond &condMT, 
			     float delta_theta=0.1, float puls_min=5., float puls_mt=9., float chi2max=1.6, FILE *out=0 );
  int     FindCompliments( EdbSegP &s, EdbPattern &pat, TObjArray &found, float chi2max, TArrayF &chiarr );
  void    SetDefaultCondBT(EdbScanCond &cond);
  void    SetDefaultCondMT(EdbScanCond &cond);

  void    OptimizeScanPath(EdbPattern &pin, EdbPattern &pout,int brick);
  int     RemoveDublets(EdbPattern &pin, EdbPattern &pout,int brick);

  bool    AddAFFtoScanSet(EdbScanSet &sc, EdbID id1, EdbID id2);
  bool    AddAFFtoScanSet(EdbScanSet &sc, int id1[4], int id2[4]);
  bool    AddAFFtoScanSet(EdbScanSet &sc, int b1, int p1, int s1, int e1,int b2, int p2, int s2, int e2);

  int     AssembleScanSet(  EdbScanSet &ss);
  int     ReadScanSetCP(    EdbScanSet &ss, EdbPVRec &ali, TCut c="1", bool do_erase=true, int minplate=-1000, int maxplate=-1000);
  int     ReadScanSetCP(    EdbID id, EdbPVRec &ali, TCut c="1", bool do_erase=true, bool do_assemble=true, int minplate=-1000, int maxplate=-1000);
  int     ReadFoundSegment( EdbID id,  EdbSegP &s, int flag=-1);
  int     ReadFoundTrack(   EdbScanSet &ss, EdbTrackP &track, int flag=-1);
  int     ReadFoundTracks(   EdbScanSet &ss,  EdbPVRec &ali, int flag=-1);
  int     ReadManFoundTracks(   EdbScanSet &ss,  EdbPVRec &ali, int flag=-1);

  void    CheckRunQualityRaw( EdbID idss ) {}
  void    CheckSetQualityRaw( EdbID idss );
  
  void    AlignOverlaps(EdbID id, EdbPattern &p1,EdbPattern &p2, TEnv &cenv, const char *suff);
  void    CheckViewOverlaps( EdbID id, TEnv &cenv );

  int     WriteScanSet(EdbID id, EdbScanSet &ss);
  EdbScanSet  *ReadScanSet(EdbID id);
  EdbScanSet  *ReadScanSetGlobal(EdbID id, bool x_marks);

  int     WriteSBTrack(EdbTrackP &t, int path, EdbID id);  //to remove?
  int     WriteSBTracks(TObjArray &tracks, EdbID id);
  TObjArray *ReadSBTracks(EdbID id);
  void    MergeSetSBT(EdbID id, EdbScanSet &ss);
  void    MergeSetSBT(EdbID id);

  void    PrepareVolumesPred(int id[4], EdbPattern &points, int before=5, int after=5, 
			     int pmin=1, int pmax=57, EdbScanSet *sc=0);

  int     MakeTracksPred(TObjArray &tracks, EdbID id, EdbLayer &layer);

  int     TestAl(int id1[4], int id2[4]);
  int     TestAl(EdbID id1, EdbID id2)   {int id14[4]; id1.Get(id14); int id24[4]; id2.Get(id24); return TestAl(id14,id24); }
  int     TestAl(EdbPattern &p1, EdbPattern &p2);
  int     TestAl(const char *cpfile1, const char *cpfile2, TCut &cut, float dz, EdbAffine2D *aff=0);

  int     ReadMarksSet(EdbMarksSet &ms, int brick, const char *filename, char spacer='_', char shape='S');
  int     WriteMarksSet(EdbMarksSet &ms, int brick, const char *filename, char spacer='_', char shape='S', int plate=1);

  int     AlignRaw(EdbID id1, EdbID id2, TEnv &cenv, EdbAffine2D *applyAff=0);
  void    AlignRawSet(EdbID id1, EdbID id2, TEnv &cenv);
  void    UpdateSetWithAff(EdbID idset, EdbAffine2D aff);
  void    UpdateSetWithAff( EdbID id, EdbID id1, EdbID id2 );
  void    UpdateSetWithAff( EdbID id,  EdbID idu );
  void    UpdateSetWithPlatePar( EdbID id );
  void    UpdateSetWithPlatePar( EdbScanSet &ss );

  void    MakeLinkSetSummary( EdbID id );
  
  void    MakeAlignSetSummary( EdbID id );
  void    MakeAlignSetSummary( EdbID id1, EdbID id2, const char *fout, const char *opt="UPDATE" );

  int     FindRawTrack( EdbTrackP &pred, EdbTrackP &found, EdbID idset, int plate,  TEnv &cenv);
  int     FindRawTrack( EdbTrackP &pred, EdbTrackP &found, EdbID idset, int plate);

  void    UpdateAlignSummaryTree(EdbID idset1, EdbID idset2, TTree &tree);
  bool    ReadAffToLayer( EdbLayer &la, EdbID id1, EdbID id2 );

  void    ExtractRawVolume(EdbID id, EdbID idnew,EdbSegP pred, int plate, TEnv &cenv);
  void    ExtractRawVolume(EdbScanSet &ss, EdbScanSet &ssnew, EdbSegP &pred, float dR);

  void SetServerRunName(const char* fname_);
  const char* GetServerRunName()const;
  
  void ReadUncorrectedBTforFoundTracks( EdbPVRec &ali );

  void    LogPrint(int brick, int level, const char *rout, const char *msgfmt, ...);
  void    Print();
  
  ClassDef(EdbScanProc,1)  // scanned data processing
};
#endif /* ROOT_EdbScanProc */
