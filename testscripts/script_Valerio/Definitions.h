#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TNamed.h"
#include "EdbSegP.h"
#include "EdbPattern.h"


// FEDRA VARIABLES

   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

// VTX TREE

   // Declaration of leaf types
   Int_t           vID;
   Float_t         vx;
   Float_t         vy;
   Float_t         vz;
   Float_t         vtx_max_aperture;
   Float_t         probability;
   Int_t           n;
   Int_t           TrackID[40];   //[n]
   Float_t         TX[40];   //[n]
   Float_t         TY[40];   //[n]
   Int_t           nseg[40];   //[n]
   Float_t         trackfill[40];   //[n]
   Int_t           incoming[40];   //[n]
   Float_t         impactparameter[40];   //[n]
   Int_t           trk_num_holes[40];   //[n]
   Int_t           trk_max_gap[40];   //[n]
   Int_t           MCEventID[40];   //[n]
   Int_t           MCTrackID[40];   //[n]
   Int_t           MCMotherID[40];   //[n]

   // List of branches
   TBranch        *b_vID;   //!
   TBranch        *b_vx;   //!
   TBranch        *b_vy;   //!
   TBranch        *b_vz;   //!
   TBranch        *b_vtx_max_aperture;   //!
   TBranch        *b_probability;   //!
   TBranch        *b_n;   //!
   TBranch        *b_TrackID;   //!
   TBranch        *b_TX;   //!
   TBranch        *b_TY;   //!
   TBranch        *b_nseg;   //!
   TBranch        *b_trackfill;   //!
   TBranch        *b_incoming;   //!
   TBranch        *b_impactparameter;   //!
   TBranch        *b_trk_num_holes;   //!
   TBranch        *b_trk_max_gap;   //!
   TBranch        *b_MCEventID;   //!
   TBranch        *b_MCTrackID;   //!
   TBranch        *b_MCMotherID;   //!


// TRACKS VARIABLES

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxt = 1;
   static constexpr Int_t kMaxs = 29;
   static constexpr Int_t kMaxsf = 29;

   // Declaration of leaf types
   Int_t           trid;
   Int_t           tnseg;
   Int_t           npl;
   Int_t           n0;
   Float_t         xv;
   Float_t         yv;
   Float_t         w;
   EdbTrackP       *t_;
   Int_t           s_;
   UInt_t          s_fUniqueID[kMaxs];   //[s_]
   UInt_t          s_fBits[kMaxs];   //[s_]
   EdbSegP         *s_EdbTrack2D;
   Int_t           s_ePID[kMaxs];   //[s_]
   Int_t           s_eID[kMaxs];   //[s_]
   Int_t           s_eVid[kMaxs][2];   //[s_]
   Int_t           s_eAid[kMaxs][2];   //[s_]
   Int_t           s_eFlag[kMaxs];   //[s_]
   Int_t           s_eTrack[kMaxs];   //[s_]
   Float_t         s_eX[kMaxs];   //[s_]
   Float_t         s_eY[kMaxs];   //[s_]
   Float_t         s_eZ[kMaxs];   //[s_]
   Float_t         s_eTX[kMaxs];   //[s_]
   Float_t         s_eTY[kMaxs];   //[s_]
   Float_t         s_eSZ[kMaxs];   //[s_]
   Float_t         s_eChi2[kMaxs];   //[s_]
   Float_t         s_eProb[kMaxs];   //[s_]
   Float_t         s_eW[kMaxs];   //[s_]
   Float_t         s_eVolume[kMaxs];   //[s_]
   Float_t         s_eDZ[kMaxs];   //[s_]
   Float_t         s_eDZem[kMaxs];   //[s_]
   Float_t         s_eP[kMaxs];   //[s_]
   Int_t           s_eMCTrack[kMaxs];   //[s_]
   Int_t           s_eMCEvt[kMaxs];   //[s_]
   EdbID           s_eScanID[kMaxs];
   Int_t           sf_;
   UInt_t          sf_fUniqueID[kMaxsf];   //[sf_]
   UInt_t          sf_fBits[kMaxsf];   //[sf_]
   EdbSegP         *sf_EdbTrack2D;
   Int_t           sf_ePID[kMaxsf];   //[sf_]
   Int_t           sf_eID[kMaxsf];   //[sf_]
   Int_t           sf_eVid[kMaxsf][2];   //[sf_]
   Int_t           sf_eAid[kMaxsf][2];   //[sf_]
   Int_t           sf_eFlag[kMaxsf];   //[sf_]
   Int_t           sf_eTrack[kMaxsf];   //[sf_]
   Float_t         sf_eX[kMaxsf];   //[sf_]
   Float_t         sf_eY[kMaxsf];   //[sf_]
   Float_t         sf_eZ[kMaxsf];   //[sf_]
   Float_t         sf_eTX[kMaxsf];   //[sf_]
   Float_t         sf_eTY[kMaxsf];   //[sf_]
   Float_t         sf_eSZ[kMaxsf];   //[sf_]
   Float_t         sf_eChi2[kMaxsf];   //[sf_]
   Float_t         sf_eProb[kMaxsf];   //[sf_]
   Float_t         sf_eW[kMaxsf];   //[sf_]
   Float_t         sf_eVolume[kMaxsf];   //[sf_]
   Float_t         sf_eDZ[kMaxsf];   //[sf_]
   Float_t         sf_eDZem[kMaxsf];   //[sf_]
   Float_t         sf_eP[kMaxsf];   //[sf_]
   Int_t           sf_eMCTrack[kMaxsf];   //[sf_]
   Int_t           sf_eMCEvt[kMaxsf];   //[sf_]
   EdbID           sf_eScanID[kMaxsf];

   // List of branches
   TBranch        *b_trid;   //!
   TBranch        *b_tnseg;   //!
   TBranch        *b_npl;   //!
   TBranch        *b_n0;   //!
   TBranch        *b_xv;   //!
   TBranch        *b_yv;   //!
   TBranch        *b_w;   //!
   TBranch        *b_t_;   //!
   TBranch        *b_s_;   //!
   TBranch        *b_s_fUniqueID;   //!
   TBranch        *b_s_fBits;   //!
   TBranch        *b_s_EdbTrack2D;   //!
   TBranch        *b_s_ePID;   //!
   TBranch        *b_s_eID;   //!
   TBranch        *b_s_eVid;   //!
   TBranch        *b_s_eAid;   //!
   TBranch        *b_s_eFlag;   //!
   TBranch        *b_s_eTrack;   //!
   TBranch        *b_s_eX;   //!
   TBranch        *b_s_eY;   //!
   TBranch        *b_s_eZ;   //!
   TBranch        *b_s_eTX;   //!
   TBranch        *b_s_eTY;   //!
   TBranch        *b_s_eSZ;   //!
   TBranch        *b_s_eChi2;   //!
   TBranch        *b_s_eProb;   //!
   TBranch        *b_s_eW;   //!
   TBranch        *b_s_eVolume;   //!
   TBranch        *b_s_eDZ;   //!
   TBranch        *b_s_eDZem;   //!
   TBranch        *b_s_eP;   //!
   TBranch        *b_s_eMCTrack;   //!
   TBranch        *b_s_eMCEvt;   //!
   TBranch        *b_s_eScanID;   //!
   TBranch        *b_sf_;   //!
   TBranch        *b_sf_fUniqueID;   //!
   TBranch        *b_sf_fBits;   //!
   TBranch        *b_sf_EdbTrack2D;   //!
   TBranch        *b_sf_ePID;   //!
   TBranch        *b_sf_eID;   //!
   TBranch        *b_sf_eVid;   //!
   TBranch        *b_sf_eAid;   //!
   TBranch        *b_sf_eFlag;   //!
   TBranch        *b_sf_eTrack;   //!
   TBranch        *b_sf_eX;   //!
   TBranch        *b_sf_eY;   //!
   TBranch        *b_sf_eZ;   //!
   TBranch        *b_sf_eTX;   //!
   TBranch        *b_sf_eTY;   //!
   TBranch        *b_sf_eSZ;   //!
   TBranch        *b_sf_eChi2;   //!
   TBranch        *b_sf_eProb;   //!
   TBranch        *b_sf_eW;   //!
   TBranch        *b_sf_eVolume;   //!
   TBranch        *b_sf_eDZ;   //!
   TBranch        *b_sf_eDZem;   //!
   TBranch        *b_sf_eP;   //!
   TBranch        *b_sf_eMCTrack;   //!
   TBranch        *b_sf_eMCEvt;   //!
   TBranch        *b_sf_eScanID;   //!


// MC VARIABLES

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxcbmroot_Stack_MCTrack = 2880;
   static constexpr Int_t kMaxcbmroot_Box_BoxPoint = 24529;
   static constexpr Int_t kMaxcbmroot_SciFi_SciFiPoint = 215;
   static constexpr Int_t kMaxcbmroot_PixelModules_PixelModulesPoint = 2326;
   static constexpr Int_t kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint = 1;
   static constexpr Int_t kMaxcbmroot_MuonTagger_MuonTaggerPoint = 642;
   static constexpr Int_t kMaxcbmroot_Event_MCEventHeader = 1;

   // Declaration of leaf types
   Int_t           MCTrack_;
   UInt_t          MCTrack_fUniqueID[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   UInt_t          MCTrack_fBits[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Int_t           MCTrack_fPdgCode[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Int_t           MCTrack_fMotherId[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fPx[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fPy[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fPz[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fM[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fStartX[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fStartY[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fStartZ[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fStartT[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Double32_t      MCTrack_fW[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Int_t           MCTrack_fProcID[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Int_t           MCTrack_fNPoints[kMaxcbmroot_Stack_MCTrack];   //[cbmroot.Stack.MCTrack_]
   Int_t           BoxPoint_;
   UInt_t          BoxPoint_fUniqueID[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   UInt_t          BoxPoint_fBits[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Int_t           BoxPoint_fTrackID[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   UInt_t          BoxPoint_fEventId[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fPx[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fPy[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fPz[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fTime[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fLength[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fELoss[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Int_t           BoxPoint_fDetectorID[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fX[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fY[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Double32_t      BoxPoint_fZ[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Int_t           BoxPoint_fPdgCode[kMaxcbmroot_Box_BoxPoint];   //[cbmroot.Box.BoxPoint_]
   Int_t           SciFiPoint_;
   UInt_t          SciFiPoint_fUniqueID[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   UInt_t          SciFiPoint_fBits[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Int_t           SciFiPoint_fTrackID[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   UInt_t          SciFiPoint_fEventId[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fPx[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fPy[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fPz[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fTime[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fLength[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fELoss[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Int_t           SciFiPoint_fDetectorID[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fX[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fY[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Double32_t      SciFiPoint_fZ[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Int_t           SciFiPoint_fPdgCode[kMaxcbmroot_SciFi_SciFiPoint];   //[cbmroot.SciFi.SciFiPoint_]
   Int_t           PixelModulesPoint_;
   UInt_t          PixelModulesPoint_fUniqueID[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   UInt_t          PixelModulesPoint_fBits[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Int_t           PixelModulesPoint_fTrackID[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   UInt_t          PixelModulesPoint_fEventId[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fPx[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fPy[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fPz[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fTime[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fLength[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fELoss[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Int_t           PixelModulesPoint_fDetectorID[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fX[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fY[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Double32_t      PixelModulesPoint_fZ[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Int_t           PixelModulesPoint_fPdgCode[kMaxcbmroot_PixelModules_PixelModulesPoint];   //[cbmroot.PixelModules.PixelModulesPoint_]
   Int_t           MufluxSpectrometerPoint_;
   UInt_t          MufluxSpectrometerPoint_fUniqueID[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   UInt_t          MufluxSpectrometerPoint_fBits[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Int_t           MufluxSpectrometerPoint_fTrackID[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   UInt_t          MufluxSpectrometerPoint_fEventId[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fPx[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fPy[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fPz[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fTime[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fLength[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fELoss[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Int_t           MufluxSpectrometerPoint_fDetectorID[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fX[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fY[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double32_t      MufluxSpectrometerPoint_fZ[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Int_t           MufluxSpectrometerPoint_fPdgCode[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Double_t        MufluxSpectrometerPoint_fdist2Wire[kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint];   //[cbmroot.MufluxSpectrometer.MufluxSpectrometerPoint_]
   Int_t           MuonTaggerPoint_;
   UInt_t          MuonTaggerPoint_fUniqueID[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   UInt_t          MuonTaggerPoint_fBits[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Int_t           MuonTaggerPoint_fTrackID[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   UInt_t          MuonTaggerPoint_fEventId[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fPx[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fPy[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fPz[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fTime[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fLength[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fELoss[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Int_t           MuonTaggerPoint_fDetectorID[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fX[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fY[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Double32_t      MuonTaggerPoint_fZ[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
   Int_t           MuonTaggerPoint_fPdgCode[kMaxcbmroot_MuonTagger_MuonTaggerPoint];   //[cbmroot.MuonTagger.MuonTaggerPoint_]
/*FairMCEventHeader *MCEventHeader_;
   UInt_t          MCEventHeader_TNamed_fUniqueID;
   UInt_t          MCEventHeader_TNamed_fBits;
   TString         MCEventHeader_TNamed_fName;
   TString         MCEventHeader_TNamed_fTitle;
   UInt_t          MCEventHeader_fRunId;
   UInt_t          MCEventHeader_fEventId;
   Double32_t      MCEventHeader_fX;
   Double32_t      MCEventHeader_fY;
   Double32_t      MCEventHeader_fZ;
   Double32_t      MCEventHeader_fT;
   Double32_t      MCEventHeader_fB;
   Int_t           MCEventHeader_fNPrim;
   Bool_t          MCEventHeader_fIsSet;
   Double32_t      MCEventHeader_fRotX;
   Double32_t      MCEventHeader_fRotY;
   Double32_t      MCEventHeader_fRotZ;
*/
   // List of branches
   TBranch        *b_cbmroot_Stack_MCTrack_;   //!
   TBranch        *b_MCTrack_fUniqueID;   //!
   TBranch        *b_MCTrack_fBits;   //!
   TBranch        *b_MCTrack_fPdgCode;   //!
   TBranch        *b_MCTrack_fMotherId;   //!
   TBranch        *b_MCTrack_fPx;   //!
   TBranch        *b_MCTrack_fPy;   //!
   TBranch        *b_MCTrack_fPz;   //!
   TBranch        *b_MCTrack_fM;   //!
   TBranch        *b_MCTrack_fStartX;   //!
   TBranch        *b_MCTrack_fStartY;   //!
   TBranch        *b_MCTrack_fStartZ;   //!
   TBranch        *b_MCTrack_fStartT;   //!
   TBranch        *b_MCTrack_fW;   //!
   TBranch        *b_MCTrack_fProcID;   //!
   TBranch        *b_MCTrack_fNPoints;   //!
   TBranch        *b_cbmroot_Box_BoxPoint_;   //!
   TBranch        *b_BoxPoint_fUniqueID;   //!
   TBranch        *b_BoxPoint_fBits;   //!
   TBranch        *b_BoxPoint_fTrackID;   //!
   TBranch        *b_BoxPoint_fEventId;   //!
   TBranch        *b_BoxPoint_fPx;   //!
   TBranch        *b_BoxPoint_fPy;   //!
   TBranch        *b_BoxPoint_fPz;   //!
   TBranch        *b_BoxPoint_fTime;   //!
   TBranch        *b_BoxPoint_fLength;   //!
   TBranch        *b_BoxPoint_fELoss;   //!
   TBranch        *b_BoxPoint_fDetectorID;   //!
   TBranch        *b_BoxPoint_fX;   //!
   TBranch        *b_BoxPoint_fY;   //!
   TBranch        *b_BoxPoint_fZ;   //!
   TBranch        *b_BoxPoint_fPdgCode;   //!
   TBranch        *b_cbmroot_SciFi_SciFiPoint_;   //!
   TBranch        *b_SciFiPoint_fUniqueID;   //!
   TBranch        *b_SciFiPoint_fBits;   //!
   TBranch        *b_SciFiPoint_fTrackID;   //!
   TBranch        *b_SciFiPoint_fEventId;   //!
   TBranch        *b_SciFiPoint_fPx;   //!
   TBranch        *b_SciFiPoint_fPy;   //!
   TBranch        *b_SciFiPoint_fPz;   //!
   TBranch        *b_SciFiPoint_fTime;   //!
   TBranch        *b_SciFiPoint_fLength;   //!
   TBranch        *b_SciFiPoint_fELoss;   //!
   TBranch        *b_SciFiPoint_fDetectorID;   //!
   TBranch        *b_SciFiPoint_fX;   //!
   TBranch        *b_SciFiPoint_fY;   //!
   TBranch        *b_SciFiPoint_fZ;   //!
   TBranch        *b_SciFiPoint_fPdgCode;   //!
   TBranch        *b_cbmroot_PixelModules_PixelModulesPoint_;   //!
   TBranch        *b_PixelModulesPoint_fUniqueID;   //!
   TBranch        *b_PixelModulesPoint_fBits;   //!
   TBranch        *b_PixelModulesPoint_fTrackID;   //!
   TBranch        *b_PixelModulesPoint_fEventId;   //!
   TBranch        *b_PixelModulesPoint_fPx;   //!
   TBranch        *b_PixelModulesPoint_fPy;   //!
   TBranch        *b_PixelModulesPoint_fPz;   //!
   TBranch        *b_PixelModulesPoint_fTime;   //!
   TBranch        *b_PixelModulesPoint_fLength;   //!
   TBranch        *b_PixelModulesPoint_fELoss;   //!
   TBranch        *b_PixelModulesPoint_fDetectorID;   //!
   TBranch        *b_PixelModulesPoint_fX;   //!
   TBranch        *b_PixelModulesPoint_fY;   //!
   TBranch        *b_PixelModulesPoint_fZ;   //!
   TBranch        *b_PixelModulesPoint_fPdgCode;   //!
   TBranch        *b_cbmroot_MufluxSpectrometer_MufluxSpectrometerPoint_;   //!
   TBranch        *b_MufluxSpectrometerPoint_fUniqueID;   //!
   TBranch        *b_MufluxSpectrometerPoint_fBits;   //!
   TBranch        *b_MufluxSpectrometerPoint_fTrackID;   //!
   TBranch        *b_MufluxSpectrometerPoint_fEventId;   //!
   TBranch        *b_MufluxSpectrometerPoint_fPx;   //!
   TBranch        *b_MufluxSpectrometerPoint_fPy;   //!
   TBranch        *b_MufluxSpectrometerPoint_fPz;   //!
   TBranch        *b_MufluxSpectrometerPoint_fTime;   //!
   TBranch        *b_MufluxSpectrometerPoint_fLength;   //!
   TBranch        *b_MufluxSpectrometerPoint_fELoss;   //!
   TBranch        *b_MufluxSpectrometerPoint_fDetectorID;   //!
   TBranch        *b_MufluxSpectrometerPoint_fX;   //!
   TBranch        *b_MufluxSpectrometerPoint_fY;   //!
   TBranch        *b_MufluxSpectrometerPoint_fZ;   //!
   TBranch        *b_MufluxSpectrometerPoint_fPdgCode;   //!
   TBranch        *b_MufluxSpectrometerPoint_fdist2Wire;   //!
   TBranch        *b_cbmroot_MuonTagger_MuonTaggerPoint_;   //!
   TBranch        *b_MuonTaggerPoint_fUniqueID;   //!
   TBranch        *b_MuonTaggerPoint_fBits;   //!
   TBranch        *b_MuonTaggerPoint_fTrackID;   //!
   TBranch        *b_MuonTaggerPoint_fEventId;   //!
   TBranch        *b_MuonTaggerPoint_fPx;   //!
   TBranch        *b_MuonTaggerPoint_fPy;   //!
   TBranch        *b_MuonTaggerPoint_fPz;   //!
   TBranch        *b_MuonTaggerPoint_fTime;   //!
   TBranch        *b_MuonTaggerPoint_fLength;   //!
   TBranch        *b_MuonTaggerPoint_fELoss;   //!
   TBranch        *b_MuonTaggerPoint_fDetectorID;   //!
   TBranch        *b_MuonTaggerPoint_fX;   //!
   TBranch        *b_MuonTaggerPoint_fY;   //!
   TBranch        *b_MuonTaggerPoint_fZ;   //!
   TBranch        *b_MuonTaggerPoint_fPdgCode;   //!
/*   TBranch        *b_cbmroot_Event_MCEventHeader_;   //!
TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fUniqueID;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fBits;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fName;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fTitle;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fRunId;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fEventId;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fX;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fY;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fZ;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fT;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fB;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fNPrim;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fIsSet;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotX;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotY;   //!
   TBranch        *b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotZ;   //!
*/



void  vtx_reader_Fedra(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("vID", &vID, &b_vID);
   fChain->SetBranchAddress("vx", &vx, &b_vx);
   fChain->SetBranchAddress("vy", &vy, &b_vy);
   fChain->SetBranchAddress("vz", &vz, &b_vz);
   fChain->SetBranchAddress("vtx_max_aperture", &vtx_max_aperture, &b_vtx_max_aperture);
   fChain->SetBranchAddress("probability", &probability, &b_probability);
   fChain->SetBranchAddress("n", &n, &b_n);
   fChain->SetBranchAddress("TrackID", TrackID, &b_TrackID);
   fChain->SetBranchAddress("TX", TX, &b_TX);
   fChain->SetBranchAddress("TY", TY, &b_TY);
   fChain->SetBranchAddress("nseg", nseg, &b_nseg);
   fChain->SetBranchAddress("trackfill", trackfill, &b_trackfill);
   fChain->SetBranchAddress("incoming", incoming, &b_incoming);
   fChain->SetBranchAddress("impactparameter", impactparameter, &b_impactparameter);
   fChain->SetBranchAddress("trk_num_holes", trk_num_holes, &b_trk_num_holes);
   fChain->SetBranchAddress("trk_max_gap", trk_max_gap, &b_trk_max_gap);
   fChain->SetBranchAddress("MCEventID", MCEventID, &b_MCEventID);
   fChain->SetBranchAddress("MCTrackID", MCTrackID, &b_MCTrackID);
   fChain->SetBranchAddress("MCMotherID", MCMotherID, &b_MCMotherID);
}


void tracks_reader_Fedra(TTree *tree)
{

   // Set object pointer
   t_ = 0;
   s_EdbTrack2D = 0;
   sf_EdbTrack2D = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("trid", &trid, &b_trid);
   fChain->SetBranchAddress("nseg", &tnseg, &b_tnseg);
   fChain->SetBranchAddress("npl", &npl, &b_npl);
   fChain->SetBranchAddress("n0", &n0, &b_n0);
   fChain->SetBranchAddress("xv", &xv, &b_xv);
   fChain->SetBranchAddress("yv", &yv, &b_yv);
   fChain->SetBranchAddress("w", &w, &b_w);
   fChain->SetBranchAddress("t.", &t_, &b_t_);
   fChain->SetBranchAddress("s", &s_, &b_s_);
   fChain->SetBranchAddress("s.fUniqueID", s_fUniqueID, &b_s_fUniqueID);
   fChain->SetBranchAddress("s.fBits", s_fBits, &b_s_fBits);
   fChain->SetBranchAddress("s.EdbTrack2D", s_EdbTrack2D, &b_s_EdbTrack2D);
   fChain->SetBranchAddress("s.ePID", s_ePID, &b_s_ePID);
   fChain->SetBranchAddress("s.eID", s_eID, &b_s_eID);
   fChain->SetBranchAddress("s.eVid[2]", s_eVid, &b_s_eVid);
   fChain->SetBranchAddress("s.eAid[2]", s_eAid, &b_s_eAid);
   fChain->SetBranchAddress("s.eFlag", s_eFlag, &b_s_eFlag);
   fChain->SetBranchAddress("s.eTrack", s_eTrack, &b_s_eTrack);
   fChain->SetBranchAddress("s.eX", s_eX, &b_s_eX);
   fChain->SetBranchAddress("s.eY", s_eY, &b_s_eY);
   fChain->SetBranchAddress("s.eZ", s_eZ, &b_s_eZ);
   fChain->SetBranchAddress("s.eTX", s_eTX, &b_s_eTX);
   fChain->SetBranchAddress("s.eTY", s_eTY, &b_s_eTY);
   fChain->SetBranchAddress("s.eSZ", s_eSZ, &b_s_eSZ);
   fChain->SetBranchAddress("s.eChi2", s_eChi2, &b_s_eChi2);
   fChain->SetBranchAddress("s.eProb", s_eProb, &b_s_eProb);
   fChain->SetBranchAddress("s.eW", s_eW, &b_s_eW);
   fChain->SetBranchAddress("s.eVolume", s_eVolume, &b_s_eVolume);
   fChain->SetBranchAddress("s.eDZ", s_eDZ, &b_s_eDZ);
   fChain->SetBranchAddress("s.eDZem", s_eDZem, &b_s_eDZem);
   fChain->SetBranchAddress("s.eP", s_eP, &b_s_eP);
   fChain->SetBranchAddress("s.eMCTrack", s_eMCTrack, &b_s_eMCTrack);
   fChain->SetBranchAddress("s.eMCEvt", s_eMCEvt, &b_s_eMCEvt);
   fChain->SetBranchAddress("s.eScanID", s_eScanID, &b_s_eScanID);
   fChain->SetBranchAddress("sf", &sf_, &b_sf_);
   fChain->SetBranchAddress("sf.fUniqueID", sf_fUniqueID, &b_sf_fUniqueID);
   fChain->SetBranchAddress("sf.fBits", sf_fBits, &b_sf_fBits);
   fChain->SetBranchAddress("sf.EdbTrack2D", sf_EdbTrack2D, &b_sf_EdbTrack2D);
   fChain->SetBranchAddress("sf.ePID", sf_ePID, &b_sf_ePID);
   fChain->SetBranchAddress("sf.eID", sf_eID, &b_sf_eID);
   fChain->SetBranchAddress("sf.eVid[2]", sf_eVid, &b_sf_eVid);
   fChain->SetBranchAddress("sf.eAid[2]", sf_eAid, &b_sf_eAid);
   fChain->SetBranchAddress("sf.eFlag", sf_eFlag, &b_sf_eFlag);
   fChain->SetBranchAddress("sf.eTrack", sf_eTrack, &b_sf_eTrack);
   fChain->SetBranchAddress("sf.eX", sf_eX, &b_sf_eX);
   fChain->SetBranchAddress("sf.eY", sf_eY, &b_sf_eY);
   fChain->SetBranchAddress("sf.eZ", sf_eZ, &b_sf_eZ);
   fChain->SetBranchAddress("sf.eTX", sf_eTX, &b_sf_eTX);
   fChain->SetBranchAddress("sf.eTY", sf_eTY, &b_sf_eTY);
   fChain->SetBranchAddress("sf.eSZ", sf_eSZ, &b_sf_eSZ);
   fChain->SetBranchAddress("sf.eChi2", sf_eChi2, &b_sf_eChi2);
   fChain->SetBranchAddress("sf.eProb", sf_eProb, &b_sf_eProb);
   fChain->SetBranchAddress("sf.eW", sf_eW, &b_sf_eW);
   fChain->SetBranchAddress("sf.eVolume", sf_eVolume, &b_sf_eVolume);
   fChain->SetBranchAddress("sf.eDZ", sf_eDZ, &b_sf_eDZ);
   fChain->SetBranchAddress("sf.eDZem", sf_eDZem, &b_sf_eDZem);
   fChain->SetBranchAddress("sf.eP", sf_eP, &b_sf_eP);
   fChain->SetBranchAddress("sf.eMCTrack", sf_eMCTrack, &b_sf_eMCTrack);
   fChain->SetBranchAddress("sf.eMCEvt", sf_eMCEvt, &b_sf_eMCEvt);
   fChain->SetBranchAddress("sf.eScanID", sf_eScanID, &b_sf_eScanID);
}


void cbmsim_reader_MC(TTree *tree)
{
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("MCTrack", &MCTrack_, &b_cbmroot_Stack_MCTrack_);
   fChain->SetBranchAddress("MCTrack.fUniqueID", MCTrack_fUniqueID, &b_MCTrack_fUniqueID);
   fChain->SetBranchAddress("MCTrack.fBits", MCTrack_fBits, &b_MCTrack_fBits);
   fChain->SetBranchAddress("MCTrack.fPdgCode", MCTrack_fPdgCode, &b_MCTrack_fPdgCode);
   fChain->SetBranchAddress("MCTrack.fMotherId", MCTrack_fMotherId, &b_MCTrack_fMotherId);
   fChain->SetBranchAddress("MCTrack.fPx", MCTrack_fPx, &b_MCTrack_fPx);
   fChain->SetBranchAddress("MCTrack.fPy", MCTrack_fPy, &b_MCTrack_fPy);
   fChain->SetBranchAddress("MCTrack.fPz", MCTrack_fPz, &b_MCTrack_fPz);
   fChain->SetBranchAddress("MCTrack.fM", MCTrack_fM, &b_MCTrack_fM);
   fChain->SetBranchAddress("MCTrack.fStartX", MCTrack_fStartX, &b_MCTrack_fStartX);
   fChain->SetBranchAddress("MCTrack.fStartY", MCTrack_fStartY, &b_MCTrack_fStartY);
   fChain->SetBranchAddress("MCTrack.fStartZ", MCTrack_fStartZ, &b_MCTrack_fStartZ);
   fChain->SetBranchAddress("MCTrack.fStartT", MCTrack_fStartT, &b_MCTrack_fStartT);
   fChain->SetBranchAddress("MCTrack.fW", MCTrack_fW, &b_MCTrack_fW);
   fChain->SetBranchAddress("MCTrack.fProcID", MCTrack_fProcID, &b_MCTrack_fProcID);
   fChain->SetBranchAddress("MCTrack.fNPoints", MCTrack_fNPoints, &b_MCTrack_fNPoints);
   fChain->SetBranchAddress("BoxPoint", &BoxPoint_, &b_cbmroot_Box_BoxPoint_);
   fChain->SetBranchAddress("BoxPoint.fUniqueID", BoxPoint_fUniqueID, &b_BoxPoint_fUniqueID);
   fChain->SetBranchAddress("BoxPoint.fBits", BoxPoint_fBits, &b_BoxPoint_fBits);
   fChain->SetBranchAddress("BoxPoint.fTrackID", BoxPoint_fTrackID, &b_BoxPoint_fTrackID);
   fChain->SetBranchAddress("BoxPoint.fEventId", BoxPoint_fEventId, &b_BoxPoint_fEventId);
   fChain->SetBranchAddress("BoxPoint.fPx", BoxPoint_fPx, &b_BoxPoint_fPx);
   fChain->SetBranchAddress("BoxPoint.fPy", BoxPoint_fPy, &b_BoxPoint_fPy);
   fChain->SetBranchAddress("BoxPoint.fPz", BoxPoint_fPz, &b_BoxPoint_fPz);
   fChain->SetBranchAddress("BoxPoint.fTime", BoxPoint_fTime, &b_BoxPoint_fTime);
   fChain->SetBranchAddress("BoxPoint.fLength", BoxPoint_fLength, &b_BoxPoint_fLength);
   fChain->SetBranchAddress("BoxPoint.fELoss", BoxPoint_fELoss, &b_BoxPoint_fELoss);
   fChain->SetBranchAddress("BoxPoint.fDetectorID", BoxPoint_fDetectorID, &b_BoxPoint_fDetectorID);
   fChain->SetBranchAddress("BoxPoint.fX", BoxPoint_fX, &b_BoxPoint_fX);
   fChain->SetBranchAddress("BoxPoint.fY", BoxPoint_fY, &b_BoxPoint_fY);
   fChain->SetBranchAddress("BoxPoint.fZ", BoxPoint_fZ, &b_BoxPoint_fZ);
   fChain->SetBranchAddress("BoxPoint.fPdgCode", BoxPoint_fPdgCode, &b_BoxPoint_fPdgCode);
   fChain->SetBranchAddress("SciFiPoint", &SciFiPoint_, &b_cbmroot_SciFi_SciFiPoint_);
   fChain->SetBranchAddress("SciFiPoint.fUniqueID", SciFiPoint_fUniqueID, &b_SciFiPoint_fUniqueID);
   fChain->SetBranchAddress("SciFiPoint.fBits", SciFiPoint_fBits, &b_SciFiPoint_fBits);
   fChain->SetBranchAddress("SciFiPoint.fTrackID", SciFiPoint_fTrackID, &b_SciFiPoint_fTrackID);
   fChain->SetBranchAddress("SciFiPoint.fEventId", SciFiPoint_fEventId, &b_SciFiPoint_fEventId);
   fChain->SetBranchAddress("SciFiPoint.fPx", SciFiPoint_fPx, &b_SciFiPoint_fPx);
   fChain->SetBranchAddress("SciFiPoint.fPy", SciFiPoint_fPy, &b_SciFiPoint_fPy);
   fChain->SetBranchAddress("SciFiPoint.fPz", SciFiPoint_fPz, &b_SciFiPoint_fPz);
   fChain->SetBranchAddress("SciFiPoint.fTime", SciFiPoint_fTime, &b_SciFiPoint_fTime);
   fChain->SetBranchAddress("SciFiPoint.fLength", SciFiPoint_fLength, &b_SciFiPoint_fLength);
   fChain->SetBranchAddress("SciFiPoint.fELoss", SciFiPoint_fELoss, &b_SciFiPoint_fELoss);
   fChain->SetBranchAddress("SciFiPoint.fDetectorID", SciFiPoint_fDetectorID, &b_SciFiPoint_fDetectorID);
   fChain->SetBranchAddress("SciFiPoint.fX", SciFiPoint_fX, &b_SciFiPoint_fX);
   fChain->SetBranchAddress("SciFiPoint.fY", SciFiPoint_fY, &b_SciFiPoint_fY);
   fChain->SetBranchAddress("SciFiPoint.fZ", SciFiPoint_fZ, &b_SciFiPoint_fZ);
   fChain->SetBranchAddress("SciFiPoint.fPdgCode", SciFiPoint_fPdgCode, &b_SciFiPoint_fPdgCode);
   fChain->SetBranchAddress("PixelModulesPoint", &PixelModulesPoint_, &b_cbmroot_PixelModules_PixelModulesPoint_);
   fChain->SetBranchAddress("PixelModulesPoint.fUniqueID", PixelModulesPoint_fUniqueID, &b_PixelModulesPoint_fUniqueID);
   fChain->SetBranchAddress("PixelModulesPoint.fBits", PixelModulesPoint_fBits, &b_PixelModulesPoint_fBits);
   fChain->SetBranchAddress("PixelModulesPoint.fTrackID", PixelModulesPoint_fTrackID, &b_PixelModulesPoint_fTrackID);
   fChain->SetBranchAddress("PixelModulesPoint.fEventId", PixelModulesPoint_fEventId, &b_PixelModulesPoint_fEventId);
   fChain->SetBranchAddress("PixelModulesPoint.fPx", PixelModulesPoint_fPx, &b_PixelModulesPoint_fPx);
   fChain->SetBranchAddress("PixelModulesPoint.fPy", PixelModulesPoint_fPy, &b_PixelModulesPoint_fPy);
   fChain->SetBranchAddress("PixelModulesPoint.fPz", PixelModulesPoint_fPz, &b_PixelModulesPoint_fPz);
   fChain->SetBranchAddress("PixelModulesPoint.fTime", PixelModulesPoint_fTime, &b_PixelModulesPoint_fTime);
   fChain->SetBranchAddress("PixelModulesPoint.fLength", PixelModulesPoint_fLength, &b_PixelModulesPoint_fLength);
   fChain->SetBranchAddress("PixelModulesPoint.fELoss", PixelModulesPoint_fELoss, &b_PixelModulesPoint_fELoss);
   fChain->SetBranchAddress("PixelModulesPoint.fDetectorID", PixelModulesPoint_fDetectorID, &b_PixelModulesPoint_fDetectorID);
   fChain->SetBranchAddress("PixelModulesPoint.fX", PixelModulesPoint_fX, &b_PixelModulesPoint_fX);
   fChain->SetBranchAddress("PixelModulesPoint.fY", PixelModulesPoint_fY, &b_PixelModulesPoint_fY);
   fChain->SetBranchAddress("PixelModulesPoint.fZ", PixelModulesPoint_fZ, &b_PixelModulesPoint_fZ);
   fChain->SetBranchAddress("PixelModulesPoint.fPdgCode", PixelModulesPoint_fPdgCode, &b_PixelModulesPoint_fPdgCode);
   fChain->SetBranchAddress("MufluxSpectrometerPoint", &MufluxSpectrometerPoint_, &b_cbmroot_MufluxSpectrometer_MufluxSpectrometerPoint_);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fUniqueID", &MufluxSpectrometerPoint_fUniqueID, &b_MufluxSpectrometerPoint_fUniqueID);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fBits", &MufluxSpectrometerPoint_fBits, &b_MufluxSpectrometerPoint_fBits);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fTrackID", &MufluxSpectrometerPoint_fTrackID, &b_MufluxSpectrometerPoint_fTrackID);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fEventId", &MufluxSpectrometerPoint_fEventId, &b_MufluxSpectrometerPoint_fEventId);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fPx", &MufluxSpectrometerPoint_fPx, &b_MufluxSpectrometerPoint_fPx);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fPy", &MufluxSpectrometerPoint_fPy, &b_MufluxSpectrometerPoint_fPy);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fPz", &MufluxSpectrometerPoint_fPz, &b_MufluxSpectrometerPoint_fPz);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fTime", &MufluxSpectrometerPoint_fTime, &b_MufluxSpectrometerPoint_fTime);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fLength", &MufluxSpectrometerPoint_fLength, &b_MufluxSpectrometerPoint_fLength);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fELoss", &MufluxSpectrometerPoint_fELoss, &b_MufluxSpectrometerPoint_fELoss);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fDetectorID", &MufluxSpectrometerPoint_fDetectorID, &b_MufluxSpectrometerPoint_fDetectorID);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fX", &MufluxSpectrometerPoint_fX, &b_MufluxSpectrometerPoint_fX);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fY", &MufluxSpectrometerPoint_fY, &b_MufluxSpectrometerPoint_fY);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fZ", &MufluxSpectrometerPoint_fZ, &b_MufluxSpectrometerPoint_fZ);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fPdgCode", &MufluxSpectrometerPoint_fPdgCode, &b_MufluxSpectrometerPoint_fPdgCode);
   fChain->SetBranchAddress("MufluxSpectrometerPoint.fdist2Wire", &MufluxSpectrometerPoint_fdist2Wire, &b_MufluxSpectrometerPoint_fdist2Wire);
   fChain->SetBranchAddress("MuonTaggerPoint", &MuonTaggerPoint_, &b_cbmroot_MuonTagger_MuonTaggerPoint_);
   fChain->SetBranchAddress("MuonTaggerPoint.fUniqueID", MuonTaggerPoint_fUniqueID, &b_MuonTaggerPoint_fUniqueID);
   fChain->SetBranchAddress("MuonTaggerPoint.fBits", MuonTaggerPoint_fBits, &b_MuonTaggerPoint_fBits);
   fChain->SetBranchAddress("MuonTaggerPoint.fTrackID", MuonTaggerPoint_fTrackID, &b_MuonTaggerPoint_fTrackID);
   fChain->SetBranchAddress("MuonTaggerPoint.fEventId", MuonTaggerPoint_fEventId, &b_MuonTaggerPoint_fEventId);
   fChain->SetBranchAddress("MuonTaggerPoint.fPx", MuonTaggerPoint_fPx, &b_MuonTaggerPoint_fPx);
   fChain->SetBranchAddress("MuonTaggerPoint.fPy", MuonTaggerPoint_fPy, &b_MuonTaggerPoint_fPy);
   fChain->SetBranchAddress("MuonTaggerPoint.fPz", MuonTaggerPoint_fPz, &b_MuonTaggerPoint_fPz);
   fChain->SetBranchAddress("MuonTaggerPoint.fTime", MuonTaggerPoint_fTime, &b_MuonTaggerPoint_fTime);
   fChain->SetBranchAddress("MuonTaggerPoint.fLength", MuonTaggerPoint_fLength, &b_MuonTaggerPoint_fLength);
   fChain->SetBranchAddress("MuonTaggerPoint.fELoss", MuonTaggerPoint_fELoss, &b_MuonTaggerPoint_fELoss);
   fChain->SetBranchAddress("MuonTaggerPoint.fDetectorID", MuonTaggerPoint_fDetectorID, &b_MuonTaggerPoint_fDetectorID);
   fChain->SetBranchAddress("MuonTaggerPoint.fX", MuonTaggerPoint_fX, &b_MuonTaggerPoint_fX);
   fChain->SetBranchAddress("MuonTaggerPoint.fY", MuonTaggerPoint_fY, &b_MuonTaggerPoint_fY);
   fChain->SetBranchAddress("MuonTaggerPoint.fZ", MuonTaggerPoint_fZ, &b_MuonTaggerPoint_fZ);
   fChain->SetBranchAddress("MuonTaggerPoint.fPdgCode", MuonTaggerPoint_fPdgCode, &b_MuonTaggerPoint_fPdgCode);
   /*fChain->SetBranchAddress("MCEventHeader.", &MCEventHeader_, &b_cbmroot_Event_MCEventHeader_);
   fChain->SetBranchAddress("MCEventHeader.TNamed.fUniqueID", &MCEventHeader_TNamed_fUniqueID, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fUniqueID);
   fChain->SetBranchAddress("MCEventHeader.TNamed.fBits", &MCEventHeader_TNamed_fBits, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fBits);
   fChain->SetBranchAddress("MCEventHeader.TNamed.fName", &MCEventHeader_TNamed_fName, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fName);
   fChain->SetBranchAddress("MCEventHeader.TNamed.fTitle", &MCEventHeader_TNamed_fTitle, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fTitle);
   fChain->SetBranchAddress("MCEventHeader.fRunId", &MCEventHeader_fRunId, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRunId);
   fChain->SetBranchAddress("MCEventHeader.fEventId", &MCEventHeader_fEventId, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fEventId);
   fChain->SetBranchAddress("MCEventHeader.fX", &MCEventHeader_fX, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fX);
   fChain->SetBranchAddress("MCEventHeader.fY", &MCEventHeader_fY, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fY);
   fChain->SetBranchAddress("MCEventHeader.fZ", &MCEventHeader_fZ, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fZ);
   fChain->SetBranchAddress("MCEventHeader.fT", &MCEventHeader_fT, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fT);
   fChain->SetBranchAddress("MCEventHeader.fB", &MCEventHeader_fB, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fB);
   fChain->SetBranchAddress("MCEventHeader.fNPrim", &MCEventHeader_fNPrim, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fNPrim);
   fChain->SetBranchAddress("MCEventHeader.fIsSet", &MCEventHeader_fIsSet, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fIsSet);
   fChain->SetBranchAddress("MCEventHeader.fRotX", &MCEventHeader_fRotX, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotX);
   fChain->SetBranchAddress("MCEventHeader.fRotY", &MCEventHeader_fRotY, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotY);
   fChain->SetBranchAddress("MCEventHeader.fRotZ", &MCEventHeader_fRotZ, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotZ);
   */
}
