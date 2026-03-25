#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "TClonesArray.h"
#include "TObject.h"
#include "TNamed.h"

// FEDRA VARIABLES

TTree          *tree;   //!pointer to the analyzed TTree or TChain
Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

// Declaration of leaf types
Int_t           vID;
Float_t         vx;
Float_t         vy;
Float_t         vz;
Float_t         vtx_max_aperture;
Float_t         probability;
Int_t           n;
Float_t         TX[54];   //[n]
Float_t         TY[54];   //[n]
Int_t           nseg[54];   //[n]
Float_t         trackfill[54];   //[n]
Int_t           incoming[54];   //[n]
Float_t         impactparameter[54];   //[n]
Int_t           trk_num_holes[54];   //[n]
Int_t           trk_max_gap[54];   //[n]
Int_t           MCEventID[54];   //[n]
Int_t           MCTrackID[54];   //[n]
Int_t           MCMotherID[54];   //[n]

// List of branches
TBranch        *b_vID;   //!
TBranch        *b_vx;   //!
TBranch        *b_vy;   //!
TBranch        *b_vz;   //!
TBranch        *b_vtx_max_aperture;   //!
TBranch        *b_probability;   //!
TBranch        *b_n;   //!
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


// MC VARIABLES

// Fixed size dimensions of array or collections stored in the TTree if any.
   static constexpr Int_t kMaxcbmroot_Stack_MCTrack = 13353;
   static constexpr Int_t kMaxcbmroot_Box_BoxPoint = 44993;
   static constexpr Int_t kMaxcbmroot_PixelModules_PixelModulesPoint = 3024;
   static constexpr Int_t kMaxcbmroot_Spectrometer_SpectrometerPoint = 42;
   static constexpr Int_t kMaxcbmroot_MufluxSpectrometer_MufluxSpectrometerPoint = 678;
   static constexpr Int_t kMaxcbmroot_MuonTagger_MuonTaggerPoint = 1871;
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
   Int_t           SpectrometerPoint_;
   UInt_t          SpectrometerPoint_fUniqueID[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   UInt_t          SpectrometerPoint_fBits[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Int_t           SpectrometerPoint_fTrackID[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   UInt_t          SpectrometerPoint_fEventId[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fPx[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fPy[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fPz[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fTime[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fLength[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fELoss[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Int_t           SpectrometerPoint_fDetectorID[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fX[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fY[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Double32_t      SpectrometerPoint_fZ[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
   Int_t           SpectrometerPoint_fPdgCode[kMaxcbmroot_Spectrometer_SpectrometerPoint];   //[cbmroot.Spectrometer.SpectrometerPoint_]
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
/* FairMCEventHeader *MCEventHeader_;
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
   TBranch        *b_cbmroot_Spectrometer_SpectrometerPoint_;   //!
   TBranch        *b_SpectrometerPoint_fUniqueID;   //!
   TBranch        *b_SpectrometerPoint_fBits;   //!
   TBranch        *b_SpectrometerPoint_fTrackID;   //!
   TBranch        *b_SpectrometerPoint_fEventId;   //!
   TBranch        *b_SpectrometerPoint_fPx;   //!
   TBranch        *b_SpectrometerPoint_fPy;   //!
   TBranch        *b_SpectrometerPoint_fPz;   //!
   TBranch        *b_SpectrometerPoint_fTime;   //!
   TBranch        *b_SpectrometerPoint_fLength;   //!
   TBranch        *b_SpectrometerPoint_fELoss;   //!
   TBranch        *b_SpectrometerPoint_fDetectorID;   //!
   TBranch        *b_SpectrometerPoint_fX;   //!
   TBranch        *b_SpectrometerPoint_fY;   //!
   TBranch        *b_SpectrometerPoint_fZ;   //!
   TBranch        *b_SpectrometerPoint_fPdgCode;   //!
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
/* TBranch        *b_cbmroot_Event_MCEventHeader_;   //!
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

void vtx_reader_Fedra(TTree *tree)
{

  if (!tree) return;
  
  fCurrent = -1;
  tree->SetMakeClass(1);
  // Set object pointer
  //MCEventHeader_ = 0;
  
   tree->SetBranchAddress("vID", &vID, &b_vID);
   tree->SetBranchAddress("vx", &vx, &b_vx);
   tree->SetBranchAddress("vy", &vy, &b_vy);
   tree->SetBranchAddress("vz", &vz, &b_vz);
   tree->SetBranchAddress("vtx_max_aperture", &vtx_max_aperture, &b_vtx_max_aperture);
   tree->SetBranchAddress("probability", &probability, &b_probability);
   tree->SetBranchAddress("n", &n, &b_n);
   tree->SetBranchAddress("TX", TX, &b_TX);
   tree->SetBranchAddress("TY", TY, &b_TY);
   tree->SetBranchAddress("nseg", nseg, &b_nseg);
   tree->SetBranchAddress("trackfill", trackfill, &b_trackfill);
   tree->SetBranchAddress("incoming", incoming, &b_incoming);
   tree->SetBranchAddress("impactparameter", impactparameter, &b_impactparameter);
   tree->SetBranchAddress("trk_num_holes", trk_num_holes, &b_trk_num_holes);
   tree->SetBranchAddress("trk_max_gap", trk_max_gap, &b_trk_max_gap);
   tree->SetBranchAddress("MCEventID", MCEventID, &b_MCEventID);
   tree->SetBranchAddress("MCTrackID", MCTrackID, &b_MCTrackID);
   tree->SetBranchAddress("MCMotherID", MCMotherID, &b_MCMotherID);
}

void vtx_reader_MC(TTree *tree)
{
   
   // Set branch addresses and branch pointers
   if (!tree) return;
   fCurrent = -1;
   tree->SetMakeClass(1);

   tree->SetBranchAddress("MCTrack", &MCTrack_, &b_cbmroot_Stack_MCTrack_);
   tree->SetBranchAddress("MCTrack.fUniqueID", MCTrack_fUniqueID, &b_MCTrack_fUniqueID);
   tree->SetBranchAddress("MCTrack.fBits", MCTrack_fBits, &b_MCTrack_fBits);
   tree->SetBranchAddress("MCTrack.fPdgCode", MCTrack_fPdgCode, &b_MCTrack_fPdgCode);
   tree->SetBranchAddress("MCTrack.fMotherId", MCTrack_fMotherId, &b_MCTrack_fMotherId);
   tree->SetBranchAddress("MCTrack.fPx", MCTrack_fPx, &b_MCTrack_fPx);
   tree->SetBranchAddress("MCTrack.fPy", MCTrack_fPy, &b_MCTrack_fPy);
   tree->SetBranchAddress("MCTrack.fPz", MCTrack_fPz, &b_MCTrack_fPz);
   tree->SetBranchAddress("MCTrack.fM", MCTrack_fM, &b_MCTrack_fM);
   tree->SetBranchAddress("MCTrack.fStartX", MCTrack_fStartX, &b_MCTrack_fStartX);
   tree->SetBranchAddress("MCTrack.fStartY", MCTrack_fStartY, &b_MCTrack_fStartY);
   tree->SetBranchAddress("MCTrack.fStartZ", MCTrack_fStartZ, &b_MCTrack_fStartZ);
   tree->SetBranchAddress("MCTrack.fStartT", MCTrack_fStartT, &b_MCTrack_fStartT);
   tree->SetBranchAddress("MCTrack.fW", MCTrack_fW, &b_MCTrack_fW);
   tree->SetBranchAddress("MCTrack.fProcID", MCTrack_fProcID, &b_MCTrack_fProcID);
   tree->SetBranchAddress("MCTrack.fNPoints", MCTrack_fNPoints, &b_MCTrack_fNPoints);
   tree->SetBranchAddress("BoxPoint", &BoxPoint_, &b_cbmroot_Box_BoxPoint_);
   tree->SetBranchAddress("BoxPoint.fUniqueID", BoxPoint_fUniqueID, &b_BoxPoint_fUniqueID);
   tree->SetBranchAddress("BoxPoint.fBits", BoxPoint_fBits, &b_BoxPoint_fBits);
   tree->SetBranchAddress("BoxPoint.fTrackID", BoxPoint_fTrackID, &b_BoxPoint_fTrackID);
   tree->SetBranchAddress("BoxPoint.fEventId", BoxPoint_fEventId, &b_BoxPoint_fEventId);
   tree->SetBranchAddress("BoxPoint.fPx", BoxPoint_fPx, &b_BoxPoint_fPx);
   tree->SetBranchAddress("BoxPoint.fPy", BoxPoint_fPy, &b_BoxPoint_fPy);
   tree->SetBranchAddress("BoxPoint.fPz", BoxPoint_fPz, &b_BoxPoint_fPz);
   tree->SetBranchAddress("BoxPoint.fTime", BoxPoint_fTime, &b_BoxPoint_fTime);
   tree->SetBranchAddress("BoxPoint.fLength", BoxPoint_fLength, &b_BoxPoint_fLength);
   tree->SetBranchAddress("BoxPoint.fELoss", BoxPoint_fELoss, &b_BoxPoint_fELoss);
   tree->SetBranchAddress("BoxPoint.fDetectorID", BoxPoint_fDetectorID, &b_BoxPoint_fDetectorID);
   tree->SetBranchAddress("BoxPoint.fX", BoxPoint_fX, &b_BoxPoint_fX);
   tree->SetBranchAddress("BoxPoint.fY", BoxPoint_fY, &b_BoxPoint_fY);
   tree->SetBranchAddress("BoxPoint.fZ", BoxPoint_fZ, &b_BoxPoint_fZ);
   tree->SetBranchAddress("BoxPoint.fPdgCode", BoxPoint_fPdgCode, &b_BoxPoint_fPdgCode);
   tree->SetBranchAddress("PixelModulesPoint", &PixelModulesPoint_, &b_cbmroot_PixelModules_PixelModulesPoint_);
   tree->SetBranchAddress("PixelModulesPoint.fUniqueID", PixelModulesPoint_fUniqueID, &b_PixelModulesPoint_fUniqueID);
   tree->SetBranchAddress("PixelModulesPoint.fBits", PixelModulesPoint_fBits, &b_PixelModulesPoint_fBits);
   tree->SetBranchAddress("PixelModulesPoint.fTrackID", PixelModulesPoint_fTrackID, &b_PixelModulesPoint_fTrackID);
   tree->SetBranchAddress("PixelModulesPoint.fEventId", PixelModulesPoint_fEventId, &b_PixelModulesPoint_fEventId);
   tree->SetBranchAddress("PixelModulesPoint.fPx", PixelModulesPoint_fPx, &b_PixelModulesPoint_fPx);
   tree->SetBranchAddress("PixelModulesPoint.fPy", PixelModulesPoint_fPy, &b_PixelModulesPoint_fPy);
   tree->SetBranchAddress("PixelModulesPoint.fPz", PixelModulesPoint_fPz, &b_PixelModulesPoint_fPz);
   tree->SetBranchAddress("PixelModulesPoint.fTime", PixelModulesPoint_fTime, &b_PixelModulesPoint_fTime);
   tree->SetBranchAddress("PixelModulesPoint.fLength", PixelModulesPoint_fLength, &b_PixelModulesPoint_fLength);
   tree->SetBranchAddress("PixelModulesPoint.fELoss", PixelModulesPoint_fELoss, &b_PixelModulesPoint_fELoss);
   tree->SetBranchAddress("PixelModulesPoint.fDetectorID", PixelModulesPoint_fDetectorID, &b_PixelModulesPoint_fDetectorID);
   tree->SetBranchAddress("PixelModulesPoint.fX", PixelModulesPoint_fX, &b_PixelModulesPoint_fX);
   tree->SetBranchAddress("PixelModulesPoint.fY", PixelModulesPoint_fY, &b_PixelModulesPoint_fY);
   tree->SetBranchAddress("PixelModulesPoint.fZ", PixelModulesPoint_fZ, &b_PixelModulesPoint_fZ);
   tree->SetBranchAddress("PixelModulesPoint.fPdgCode", PixelModulesPoint_fPdgCode, &b_PixelModulesPoint_fPdgCode);
   tree->SetBranchAddress("SpectrometerPoint", &SpectrometerPoint_, &b_cbmroot_Spectrometer_SpectrometerPoint_);
   tree->SetBranchAddress("SpectrometerPoint.fUniqueID", SpectrometerPoint_fUniqueID, &b_SpectrometerPoint_fUniqueID);
   tree->SetBranchAddress("SpectrometerPoint.fBits", SpectrometerPoint_fBits, &b_SpectrometerPoint_fBits);
   tree->SetBranchAddress("SpectrometerPoint.fTrackID", SpectrometerPoint_fTrackID, &b_SpectrometerPoint_fTrackID);
   tree->SetBranchAddress("SpectrometerPoint.fEventId", SpectrometerPoint_fEventId, &b_SpectrometerPoint_fEventId);
   tree->SetBranchAddress("SpectrometerPoint.fPx", SpectrometerPoint_fPx, &b_SpectrometerPoint_fPx);
   tree->SetBranchAddress("SpectrometerPoint.fPy", SpectrometerPoint_fPy, &b_SpectrometerPoint_fPy);
   tree->SetBranchAddress("SpectrometerPoint.fPz", SpectrometerPoint_fPz, &b_SpectrometerPoint_fPz);
   tree->SetBranchAddress("SpectrometerPoint.fTime", SpectrometerPoint_fTime, &b_SpectrometerPoint_fTime);
   tree->SetBranchAddress("SpectrometerPoint.fLength", SpectrometerPoint_fLength, &b_SpectrometerPoint_fLength);
   tree->SetBranchAddress("SpectrometerPoint.fELoss", SpectrometerPoint_fELoss, &b_SpectrometerPoint_fELoss);
   tree->SetBranchAddress("SpectrometerPoint.fDetectorID", SpectrometerPoint_fDetectorID, &b_SpectrometerPoint_fDetectorID);
   tree->SetBranchAddress("SpectrometerPoint.fX", SpectrometerPoint_fX, &b_SpectrometerPoint_fX);
   tree->SetBranchAddress("SpectrometerPoint.fY", SpectrometerPoint_fY, &b_SpectrometerPoint_fY);
   tree->SetBranchAddress("SpectrometerPoint.fZ", SpectrometerPoint_fZ, &b_SpectrometerPoint_fZ);
   tree->SetBranchAddress("SpectrometerPoint.fPdgCode", SpectrometerPoint_fPdgCode, &b_SpectrometerPoint_fPdgCode);
   tree->SetBranchAddress("MufluxSpectrometerPoint", &MufluxSpectrometerPoint_, &b_cbmroot_MufluxSpectrometer_MufluxSpectrometerPoint_);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fUniqueID", MufluxSpectrometerPoint_fUniqueID, &b_MufluxSpectrometerPoint_fUniqueID);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fBits", MufluxSpectrometerPoint_fBits, &b_MufluxSpectrometerPoint_fBits);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fTrackID", MufluxSpectrometerPoint_fTrackID, &b_MufluxSpectrometerPoint_fTrackID);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fEventId", MufluxSpectrometerPoint_fEventId, &b_MufluxSpectrometerPoint_fEventId);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fPx", MufluxSpectrometerPoint_fPx, &b_MufluxSpectrometerPoint_fPx);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fPy", MufluxSpectrometerPoint_fPy, &b_MufluxSpectrometerPoint_fPy);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fPz", MufluxSpectrometerPoint_fPz, &b_MufluxSpectrometerPoint_fPz);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fTime", MufluxSpectrometerPoint_fTime, &b_MufluxSpectrometerPoint_fTime);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fLength", MufluxSpectrometerPoint_fLength, &b_MufluxSpectrometerPoint_fLength);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fELoss", MufluxSpectrometerPoint_fELoss, &b_MufluxSpectrometerPoint_fELoss);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fDetectorID", MufluxSpectrometerPoint_fDetectorID, &b_MufluxSpectrometerPoint_fDetectorID);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fX", MufluxSpectrometerPoint_fX, &b_MufluxSpectrometerPoint_fX);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fY", MufluxSpectrometerPoint_fY, &b_MufluxSpectrometerPoint_fY);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fZ", MufluxSpectrometerPoint_fZ, &b_MufluxSpectrometerPoint_fZ);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fPdgCode", MufluxSpectrometerPoint_fPdgCode, &b_MufluxSpectrometerPoint_fPdgCode);
   tree->SetBranchAddress("MufluxSpectrometerPoint.fdist2Wire", MufluxSpectrometerPoint_fdist2Wire, &b_MufluxSpectrometerPoint_fdist2Wire);
   tree->SetBranchAddress("MuonTaggerPoint", &MuonTaggerPoint_, &b_cbmroot_MuonTagger_MuonTaggerPoint_);
   tree->SetBranchAddress("MuonTaggerPoint.fUniqueID", MuonTaggerPoint_fUniqueID, &b_MuonTaggerPoint_fUniqueID);
   tree->SetBranchAddress("MuonTaggerPoint.fBits", MuonTaggerPoint_fBits, &b_MuonTaggerPoint_fBits);
   tree->SetBranchAddress("MuonTaggerPoint.fTrackID", MuonTaggerPoint_fTrackID, &b_MuonTaggerPoint_fTrackID);
   tree->SetBranchAddress("MuonTaggerPoint.fEventId", MuonTaggerPoint_fEventId, &b_MuonTaggerPoint_fEventId);
   tree->SetBranchAddress("MuonTaggerPoint.fPx", MuonTaggerPoint_fPx, &b_MuonTaggerPoint_fPx);
   tree->SetBranchAddress("MuonTaggerPoint.fPy", MuonTaggerPoint_fPy, &b_MuonTaggerPoint_fPy);
   tree->SetBranchAddress("MuonTaggerPoint.fPz", MuonTaggerPoint_fPz, &b_MuonTaggerPoint_fPz);
   tree->SetBranchAddress("MuonTaggerPoint.fTime", MuonTaggerPoint_fTime, &b_MuonTaggerPoint_fTime);
   tree->SetBranchAddress("MuonTaggerPoint.fLength", MuonTaggerPoint_fLength, &b_MuonTaggerPoint_fLength);
   tree->SetBranchAddress("MuonTaggerPoint.fELoss", MuonTaggerPoint_fELoss, &b_MuonTaggerPoint_fELoss);
   tree->SetBranchAddress("MuonTaggerPoint.fDetectorID", MuonTaggerPoint_fDetectorID, &b_MuonTaggerPoint_fDetectorID);
   tree->SetBranchAddress("MuonTaggerPoint.fX", MuonTaggerPoint_fX, &b_MuonTaggerPoint_fX);
   tree->SetBranchAddress("MuonTaggerPoint.fY", MuonTaggerPoint_fY, &b_MuonTaggerPoint_fY);
   tree->SetBranchAddress("MuonTaggerPoint.fZ", MuonTaggerPoint_fZ, &b_MuonTaggerPoint_fZ);
   tree->SetBranchAddress("MuonTaggerPoint.fPdgCode", MuonTaggerPoint_fPdgCode, &b_MuonTaggerPoint_fPdgCode);
   /*tree->SetBranchAddress("MCEventHeader.", &MCEventHeader_, &b_cbmroot_Event_MCEventHeader_);
   tree->SetBranchAddress("MCEventHeader.TNamed.fUniqueID", &MCEventHeader_TNamed_fUniqueID, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fUniqueID);
   tree->SetBranchAddress("MCEventHeader.TNamed.fBits", &MCEventHeader_TNamed_fBits, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fBits);
   tree->SetBranchAddress("MCEventHeader.TNamed.fName", &MCEventHeader_TNamed_fName, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fName);
   tree->SetBranchAddress("MCEventHeader.TNamed.fTitle", &MCEventHeader_TNamed_fTitle, &b_MCEventHeader_cbmroot_Event_MCEventHeader_TNamed_fTitle);
   tree->SetBranchAddress("MCEventHeader.fRunId", &MCEventHeader_fRunId, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRunId);
   tree->SetBranchAddress("MCEventHeader.fEventId", &MCEventHeader_fEventId, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fEventId);
   tree->SetBranchAddress("MCEventHeader.fX", &MCEventHeader_fX, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fX);
   tree->SetBranchAddress("MCEventHeader.fY", &MCEventHeader_fY, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fY);
   tree->SetBranchAddress("MCEventHeader.fZ", &MCEventHeader_fZ, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fZ);
   tree->SetBranchAddress("MCEventHeader.fT", &MCEventHeader_fT, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fT);
   tree->SetBranchAddress("MCEventHeader.fB", &MCEventHeader_fB, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fB);
   tree->SetBranchAddress("MCEventHeader.fNPrim", &MCEventHeader_fNPrim, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fNPrim);
   tree->SetBranchAddress("MCEventHeader.fIsSet", &MCEventHeader_fIsSet, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fIsSet);
   tree->SetBranchAddress("MCEventHeader.fRotX", &MCEventHeader_fRotX, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotX);
   tree->SetBranchAddress("MCEventHeader.fRotY", &MCEventHeader_fRotY, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotY);
   tree->SetBranchAddress("MCEventHeader.fRotZ", &MCEventHeader_fRotZ, &b_MCEventHeader_cbmroot_Event_MCEventHeader_fRotZ);
    */
     }
