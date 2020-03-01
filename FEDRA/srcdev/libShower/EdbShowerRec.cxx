//////////////////////////////////////////////////////////////////////////
//                                                                      //
// EdbShowerRec                                                         //
//                                                                      //
// Base Class for the shower reconstruction in OPERA brick.             //
// This class has been initially developed by F. Juget, G. Lutter       //
// and is ongoingly modified and developed by Frank Meisel, which is    //
// currently the corresponding author: frank.meisel@lhep.unibe.ch       //
//                                                                      //
// The main goal of this library is to make the shower search in the    //
// scanned emulsion data as easy as possible, and as most automatically.//
// There exist also two manuals, to be found under "documentation". In  //
// these some general handling of the package is explained.             //
//                                                                      //
// The EdbShowerRec class makes following things:                       //
// Starting from the emulsion scanned data (represented as an EdbPVRec  //
// object in FEDRA), it does search and reconstruction of showers,      //
// based on different algorithms, which the user can set indivially     //
// (normaly not necessary). After that, we try to identify              //
// characteristics of the shower, as it is: energy and particle ID.     //
//                                                                      //
//                                                                      //
// Nota Bene:   all lines in the code trailing with   "///"             //
//              are todo comments.                                      //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

#include "EdbShowerRec.h"
#include "EdbShowerAlg.h"
#include "EdbPVRQuality.h"

ClassImp(EdbShowerRec)
//----------------------------------------------------------------
EdbShowerRec::EdbShowerRec()
{
    // Default Constructor:
    Init();
    Set0();

}
//----------------------------------------------------------------
EdbShowerRec::EdbShowerRec(TObjArray* InBTArray, int num,int MAXPLATE,  int DATA, int piece2, int piece2par,int UPDOWN,EdbPVRec  *pvr)
{
    // Constructor giving a TObjArray of EdbSegments as Initiator Basetracks.
    // Old style constructor: Fill the array variables (arrays with 10000 indices..)
    // from the TObjArray of EdbSegP.
    // Though this is the "original" reconstruction function, it should not be used anymore
    // since it just does everything mixed together and the different steps are not to
    // clear for the user.

    Init();
    Set0();

    // Interim arrays, code structure is copied from shower_btr.C :
    double x0[10000];
    double y0[10000];
    double z0[10000];
    double tx0[10000];
    double ty0[10000];
    double chi20[10000];
    double P0[10000];
    int   W0[10000];
    int   Flag0[10000];
    int   TRid[10000];
    double Esim[10000];
    int   id[10000];
    int   Pid[10000];
    int Ncand=0;
    float Zoff;

    int nentries = InBTArray->GetEntries();
    EdbSegP *s  = 0;
    for (int i=0; i<nentries; i++) {
        s=(EdbSegP*)InBTArray->At(i);
        x0[Ncand]  =  s->X();
        y0[Ncand]  =  s->Y();
        Zoff       =  s->Z();
        z0[Ncand]  =  Zoff;
        tx0[Ncand] =  s->TX();
        ty0[Ncand] =  s->TY();
        chi20[Ncand] = s->Chi2();
        W0[Ncand] = s->W();
        P0[Ncand] = s->P();
        Flag0[Ncand] = s->Flag();
        TRid[Ncand] = s->MCEvt(); //this has to be a positive number
        Esim[Ncand] = s->P();
        id[Ncand]  = s->ID();
        Pid[Ncand] =    Zoff/1273+1;
        if (gEDBDEBUGLEVEL>2) cout << s->X() << " " << s->Y() << " " << s->Z() << " " << s->TX() << " "  << s->TY() << " " << Pid[Ncand] << " " <<  Esim[Ncand] << "  "  << TRid[Ncand]  << endl;
        cout << chi20[Ncand] << " " <<  W0[Ncand] << "  "  << P0[Ncand] << "  "  << Flag0[Ncand]  << endl;
        Ncand++;
    }

    // Start the actual reconstruction function:
    rec(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,chi20,W0,P0,Flag0,Pid,id,TRid,Esim,piece2,piece2par,UPDOWN,pvr);

}


//----------------------------------------------------------------
EdbShowerRec::EdbShowerRec(EdbPVRec  *pvr)
{
    // Constructor giving a EdbPVRec object directly.
    // We try to do the reconstruction in the following way:
    // a) if EdbPVRec has tracking already done then take eTracks for
    //    Inititator Basetracks
    // b) if EdbPVRec has vertexing already done then take eVertex for
    //    additional cut of the Initiator Basetracks which are in
    //    an IP of less than 250 microns.
    // c) Start BGEstimation, Reconstruction, Energy, ID  directly
    //    from this EdbPVRec object
    //
    // After these steps showers are stored as an TObjarray of EdbTrackP
    // in this class and can be accessed.
    // Showers can also be written in the linked_tracks.root file, or
    // in the Shower.root file. But this file stores the showers just as
    // a tree with leafs and the object structure is lost.

    cout << "----------------------------------------------------------" << endl;
    cout << "EdbShowerRec::EdbShowerRec(EdbPVRec  *pvr) " << endl;
    cout << "  Constructor giving a EdbPVRec object directly. We try to do the reconstruction in the following way: " << endl;
    cout << "  a) if EdbPVRec has tracking already done then take eTracks for Inititator Basetracks."<< endl;
    cout << "  b) if EdbPVRec has vertexing already done then take eVertex for additional cut of the Initiator Basetracks which are in an IP of less than 250 microns." << endl;
    cout << "  c) Start BGEstimation, Reconstruction, Energy, ID  directly from this EdbPVRec object" << endl;
    cout <<"----------  (some of step c are still: )    T O D O -----------------------" << endl;
    cout << endl;
    cout << "After these steps showers are stored as an TObjarray of EdbTrackP "<< endl;
    cout << "in this class and can be accessed. " << endl;
    cout << "Showers can also be written in the linked_tracks.root file, or  " << endl;
    cout << "in the Shower.root file. But this file stores the showers just as " << endl;
    cout << "a tree with leafs and the object structure is lost. " << endl;
    cout << endl;
    cout << "WARNING:: EdbShowerRec::EdbShowerRec(EdbPVRec* pvr)   DO NOTHING YET !!!" << endl;
    cout << "PLEASE DO THESE STEPS MANUALLY UNTIL SHOWER LIBRARY IS FULLY FUNCTIONAL" << endl;
    cout << "By stopping here, we avoid misunderstandings when all in one is done but not yet fully tested." << endl;
    cout << endl;
    cout << "----------------------------------------------------------" << endl;

    Init();
    Set0();
    SetEdbPVRec(pvr);

}

//----------------------------------------------------------------

EdbShowerRec::EdbShowerRec(EdbPVRec  *pvr, Int_t OptionType)
{
    // Constructor giving a EdbPVRec object directly and an option number:
    cout << "----------------------------------------------------------" << endl;
    cout << "EdbShowerRec::EdbShowerRec(EdbPVRec  *pvr, Int_t OptionType) " << endl;
    cout << "   Constructor giving a EdbPVRec object directly and an option number:" << endl;
    cout << "   OptionType = " << OptionType << endl;
    cout << "----------------------------------------------------------" << endl;
    cout << "   0: Set Initial Values, do nothing more." << endl;
    cout << "   1: Check EdbPVRec object (presence, bg dens), do nothing more." << endl;
    cout << "   2: Check other Input objects (linked_tracks.root, feedback file), do nothing more." << endl;
    cout << "   3: Do normal (electron) reconstruction using Cone Tube algorithm (OI implementation)" << endl;
    cout << "   4: Do pair search (gamma) reconstruction using  GS algorithm." << endl;
    cout << "   5: Do both Reco and ID+Energy estimations." << endl;
    cout << "  -1: Do everything possible!" << endl;
    cout << "----------------------------------------------------------" << endl;


    Init();
    Set0();

    if (OptionType==1) {
        EdbPVRQuality* qual = new EdbPVRQuality(pvr);
        EdbPVRec  *pvr_new=qual->GetEdbPVRec();
        SetEdbPVRec(pvr_new);
    }
    if (OptionType==2) {
        // check here if original pvr has tracks and/or vertices
        // check here if linked_tracks.root exists, if feedback file exits
        cout << "if (OptionType==2) .." << endl;
    }
    if (OptionType==3) {
        // start reco for OI algo here...
        cout << "if (OptionType==3) .." << endl;
    }
    if (OptionType==4) {
        // start reco for OI algo here...
        cout << "if (OptionType==4) .." << endl;
    }
    if (OptionType==5) {
        // Invoke ID and Energy estimation here...
        cout << "if (OptionType==5) .." << endl;
    }



}

//----------------------------------------------------------------
EdbShowerRec::~EdbShowerRec()
{
    // Default Destructor.
    delete eInBTArray;
    delete eRecoShowerArray;
}




//----------------------------------------------------------------
void EdbShowerRec::Init()
{
    // Create (one time) needed internal objects

    // Check if eInBTArray array exits, otherwise create.
    if (!eInitiatorBTIsCreated) {
        eInBTArray=new TObjArray();
        eInitiatorBTIsCreated=kTRUE;
        eInBTArrayN=0;
    }

    // Create Array For storing the reconstructed showers.
    eRecoShowerArray=new TObjArray();
    eRecoShowerArrayN=0;

    return;
}



//----------------------------------------------------------------
void EdbShowerRec::Set0()
{
    // Reset all internal variables to default values.

    eUseAliSub=kFALSE;
    // This is new/changed (06.June 2011):
    // If we used eUseAliSub=true, then it results in large memory
    // consumption when having many InitiatorBTs. To avoid this the standard
    // behaviour is reconstruction on normal EdbPVRec object. Only in case
    // of either large MC/Data ratio or only few InitiatorBTs, we do it this
    // way. See also the function: CheckEdbPVRecSetUseAliSub.

    eShowerTreeIsDone=kFALSE;
    eDoEnergyCalculation=kFALSE;
    eEnergyIsDone=kFALSE;
    eShowersN=0;

    eAlgoParameterConeRadius=800;      // Radius of spanning cone from first BT
    eAlgoParameterConeAngle=0.02;       // Opening angle of spanning cone from first BT
    eAlgoParameterConnectionDR=150;      // Connection Criterium: delta R
    eAlgoParameterConnectionDT=0.15;      // Connection Criterium: delta T
    eAlgoParameterNPropagation=3;      // N plates backpropagating

    eFilename_LinkedTracks="linked_tracks.root";

    eQualityPar[0]=0.12;
    eQualityPar[1]=1.0;
    eUseQualityPar=kFALSE;

    mlp1=mlp2=mlp3=mlp4=0;
}

//-------------------------------------------------------------------
void EdbShowerRec::InitPiece(EdbDataPiece &piece, const char *cpfile, const char *parfile)
{
    piece.eFileNameCP  = cpfile;
    piece.eFileNamePar = parfile;
    piece.TakePiecePar();
    //piece.Print();
}
//-------------------------------------------------------------------
void EdbShowerRec::initproc( const char *def, int iopt,   const char *rcut="1" )
{
    dproc = new EdbDataProc(def);
    dproc->InitVolume(iopt, rcut);
    gAli = dproc->PVR();
}
//-------------------------------------------------------------------
void EdbShowerRec::SaveResults()
{
    fileout->cd();
    treesaveb->Write();
    fileout->Close();
    fileout->Delete();
}
//-------------------------------------------------------------------


//-------------------------------------------------------------------
//void EdbShowerRec::rec(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, int *TRid, int piece2, int piece2par,int DOWN, float Rcut=100., float Tcut=0.05)
// void EdbShowerRec::rec(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, int *plate0, int *id0, int *TRid, double *Esim, int piece2, int piece2par, int DOWN)
// {
//     //if (DOWN == 1) {recdown(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,TRid,piece2, piece2par,Rcut,Tcut);}
//   if (DOWN == 1) {recdown(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,plate0,id0,TRid,Esim,piece2, piece2par);}
//   if (DOWN == 0) {recup(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,plate0,id0,TRid,Esim,piece2, piece2par);}
// }


//-------------------------------------------------------------------
void EdbShowerRec::rec(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, double* chi20, int* W0, double* P0, int* Flag0, int *plate0, int *id0, int *TRid, double *Esim, int piece2, int piece2par, int DOWN,EdbPVRec  *pvr)
{
    // General reconstruction intermediate function to differentiate between upstream
    // or downstream reco.
    // This function contains also arrays for Chi2,W,P,Flag of the first BT...

//    if (DOWN == 1) {
//        recdown(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,chi20,W0,P0,Flag0,plate0,id0,TRid,Esim,mcDigitIndexTop0,mcDigitIndexBottom0,piece2, piece2par,pvr);
//    }
//    if (DOWN == 0) {
//        recup(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,chi20,W0,P0,Flag0,plate0,id0,TRid,Esim,mcDigitIndexTop0,mcDigitIndexBottom0,piece2, piece2par,pvr);
//    }
    return;
}



//-------------------------------------------------------------------
void EdbShowerRec::rec(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, int *plate0, int *id0, int *TRid, double *Esim, int *mcDigitIndexTop0, int *mcDigitIndexBottom0,int piece2, int piece2par,int DOWN,EdbPVRec  *pvr)
{
    // Another reconstruction function ...
    double chi20[10000];
    double P0[10000];
    int W0[10000];
    int Flag0[10000];
    chi20[0]=-9999999;
    W0[0]=-9999999;
    P0[0]=-9999999;
    Flag0[0]=-9999999;
    if (DOWN == 1) {
        recdown(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,chi20,W0,P0,Flag0,plate0,id0,TRid,Esim,mcDigitIndexTop0,mcDigitIndexBottom0,piece2, piece2par,pvr);
    }
    if (DOWN == 0) {
        recup(num,MAXPLATE,DATA,Ncand,x0,y0,z0,tx0,ty0,chi20,W0,P0,Flag0,plate0,id0,TRid,Esim,mcDigitIndexTop0,mcDigitIndexBottom0,piece2, piece2par,pvr);
    }
    return;
}

//-------------------------------------------------------------------

void EdbShowerRec::rec(TObjArray *sarr, EdbPVRec  *pvr)
{
    // Another reconstruction function ...
    int Ncand = 0;
    double x0[10000];
    double y0[10000];
    double z0[10000];
    double tx0[10000];
    double ty0[10000];
    double chi20[10000];
    double P0[10000];
    int   W0[10000];
    int   Flag0[10000];
    int   TRid[10000];
    double Esim[10000];
    int   id[10000];
    int   Pid[10000];
    float Zoff;
    int nentries = sarr->GetEntries();
    EdbSegP *s  = 0;
    for (int i=0; i<nentries; i++) {

        s= (EdbSegP*)sarr->At(i);
        x0[Ncand]  =  s->X();
        y0[Ncand]  =  s->Y();
        Zoff       =  s->Z();
        z0[Ncand]  =  Zoff;
        tx0[Ncand] =  s->TX();
        ty0[Ncand] =  s->TY();
        chi20[Ncand] = s->Chi2();
        W0[Ncand] = s->W();
        P0[Ncand] = s->P();
        Flag0[Ncand] = s->Flag();
        TRid[Ncand] = s->MCEvt(); //this has to be a positive number
        Esim[Ncand] = s->P();
        id[Ncand]  = s->ID();
        Pid[Ncand] =    Zoff/1273+1;
        //cout << s->X() << " " << s->Y() << " " << s->Z() << " " << s->TX() << " "  << s->TY() << " " << Pid[Ncand] << " " <<  Esim[Ncand] << "  "  << TRid[Ncand]  << endl;
        //cout << chi20[Ncand] << " " <<  W0[Ncand] << "  "  << P0[Ncand] << "  "  << Flag0[Ncand]  << endl;
        Ncand++;
    }
    rec(0,57,0,Ncand,x0,y0,z0,tx0,ty0,chi20,W0,P0,Flag0,Pid,id,TRid,Esim,1,1,1,pvr);

}

//-------------------------------------------------------------------



//-------------------------------------------------------------------
//-------------------------------------------------------------------
//void EdbShowerRec::recdown(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, int *TRid, int piece2, int piece2par, float Rcut=100., float Tcut=0.05)
void EdbShowerRec::recdown(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, double* chi20, int* W0, double* P0,int* Flag0,  int *plate0, int *id0, int *TRid, double *Esim, int *mcDigitIndexTop0, int *mcDigitIndexBottom0, int piece2, int piece2par,EdbPVRec  *pvr)
{
    // Old style implementation of downstream reconstruction.
    // Starts with arrays that contain Initiator Basetracks,
    // and uses the EdbPVRec object for source data.
    // Reconstruction algorithm is the standard "ConeTube" basetrack
    // connection algorithm.
    // Reconstructed Showers are first saved as Float_t data in a
    // Tree, afterwards they will be converted into an EdbTrackP array.
    //

    // shower tree definition
    treesaveb = new TTree("treebranch","tree of branchtrack");
    treesaveb->Branch("number_eventb",&number_eventb,"number_eventb/I");
    treesaveb->Branch("sizeb",&sizeb,"sizeb/I");
    treesaveb->Branch("sizeb15",&sizeb15,"sizeb15/I");
    treesaveb->Branch("sizeb20",&sizeb20,"sizeb20/I");
    treesaveb->Branch("sizeb30",&sizeb30,"sizeb30/I");
    treesaveb->Branch("E_MC",&E_MC,"E_MC/F");

    treesaveb->Branch("isizeb",&isizeb,"isizeb/I");
    treesaveb->Branch("showerID",&showerID,"showerID/I");
    treesaveb->Branch("idb",idb,"idb[sizeb]/I");
    treesaveb->Branch("plateb",plateb,"plateb[sizeb]/I");
    treesaveb->Branch("xb",xb,"xb[sizeb]/F");
    treesaveb->Branch("yb",yb,"yb[sizeb]/F");
    treesaveb->Branch("zb",zb,"zb[sizeb]/F");
    treesaveb->Branch("txb",txb,"txb[sizeb]/F");
    treesaveb->Branch("tyb",tyb,"tyb[sizeb]/F");
    treesaveb->Branch("nfilmb",nfilmb,"nfilmb[sizeb]/I");
    treesaveb->Branch("ntrace1simub",ntrace1simub,"ntrace1simu[sizeb]/I");
    treesaveb->Branch("ntrace2simub",ntrace2simub,"ntrace2simu[sizeb]/I"); // s->W()
    treesaveb->Branch("ntrace3simub",ntrace3simub,"ntrace3simu[sizeb]/F"); // s->P()
    treesaveb->Branch("ntrace4simub",ntrace4simub,"ntrace4simu[sizeb]/I"); // s->Flag()
    treesaveb->Branch("chi2btkb",chi2btkb,"chi2btkb[sizeb]/F");
    treesaveb->Branch("deltarb",deltarb,"deltarb[sizeb]/F");
    treesaveb->Branch("deltathetab",deltathetab,"deltathetab[sizeb]/F");
    treesaveb->Branch("deltaxb",deltaxb,"deltaxb[sizeb]/F");
    treesaveb->Branch("deltayb",deltayb,"deltayb[sizeb]/F");
    treesaveb->Branch("tagprimary",tagprimary,"tagprimary[sizeb]/F");
    treesaveb->Branch("purityb",&purityb,"purityb/F");
    treesaveb->Branch("mcDigitIndexTop",mcDigitIndexTop,"mcDigitIndexTop[sizeb]/I");
    treesaveb->Branch("mcDigitIndexBottom",mcDigitIndexBottom,"mcDigitIndexBottom[sizeb]/I");

    fileout = new TFile("shower1.root","RECREATE");

    //additionally needed, if not given by constructor.
    // this is only a temporary solution. When giving over the constructor with
    // segments of EdbSegP class, they will contain the needed information then.
    cout << " void EdbShowerRec::recdown()  Additionally needed, if not given by constructor: " << endl;
    cout << " void EdbShowerRec::recdown()  this is only a temporary solution. When giving over the constructor with " << endl;
    cout << " void EdbShowerRec::recdown()  segments of EdbSegP class, they will contain the needed information then." << endl;
    cout << " void EdbShowerRec::recdown()  For the momemt, if Chi2, W,P,Flag of the starting BT is not given..." << endl;
    cout << " void EdbShowerRec::recdown()  It is assumed to be   ...   -9999999  ..." << endl;
//   float chi20[10000];
//   float P0[10000];
//   int W0[10000];
//   int Flag0[10000];
//   chi20[0]=-9999999;
//   W0[0]=-9999999;
//   P0[0]=-9999999;
//   Flag0[0]=-9999999;

    float Z0, SX0,SY0;
    float X0, Y0, Xss,Yss;
    float Delta, rayon,rayon2, diff, delta;
    float Xe[4], Ye[4];
    char fname2[128];

    char nme[64];
    float  a11, a12, a21, a22, bb1, bb2;
    a11=1;
    a12=0;
    a21=0;
    a22=1;
    bb1=0;
    bb2=0;
    float Zoff;
    double ZZ0[60];
    int ind;
    const Int_t eNSM = 10000;
    const Int_t eNTM = 1000;

    TArrayF  xbb(eNTM*eNSM);
    TArrayF  ybb(eNTM*eNSM);
    TArrayF  zbb(eNTM*eNSM);
    TArrayF  txbb(eNTM*eNSM);
    TArrayF  tybb(eNTM*eNSM);
    TArrayF  chi2bb(eNTM*eNSM);
    TArrayI  sigbb(eNTM*eNSM);
    TArrayI  wbb(eNTM*eNSM);
    TArrayI  add(eNTM*eNSM);
    TArrayF  Pbb(eNTM*eNSM);
    TArrayI  Flagbb(eNTM*eNSM);

    for (int i=0; i<eNSM; i++)
    {
        for (int j=0; j<eNTM; j++)
        {
            ind = i*eNTM+j;
            xbb[ind];
            ybb[ind];
            zbb[ind];
            txbb[ind];
            tybb[ind];
            chi2bb[ind];
            wbb[ind];
            sigbb[ind];
            Pbb[ind];
            Flagbb[ind];
        }
    }

    TList  *eS = new TList;
    for (int ii=0; ii<Ncand; ii++)
    {
        EdbShowerRec* shower = new EdbShowerRec();

        shower->SetID(ii);
        shower->SetX0(x0[ii]);
        shower->SetY0(y0[ii]);
        shower->SetZ0(z0[ii]);
        shower->SetTx(tx0[ii]);
        shower->SetTy(ty0[ii]);
        shower->SetTrID(TRid[ii]);
        shower->SetChi20(chi20[ii]);
        shower->SetW0(W0[ii]);
        shower->SetP0(P0[ii]);
        shower->SetFlag0(Flag0[ii]);


        shower->SetDeltarb(200.,0);
        shower->SetDelthetab(0.5,0);
        shower->SetPrimary(0,0);
        shower->SetNFilmb(1,0);
        shower->SetIDb(id0[ii],0);
        shower->SetPlateb(plate0[ii],0);
        shower->SetXb(x0[ii],0);
        shower->SetYb(y0[ii],0);
        shower->SetZb(z0[ii],0);
        shower->SetTXb(tx0[ii],0);
        shower->SetTYb(ty0[ii],0);
        shower->SetChi2b(chi20[ii],0);
        shower->SetWb(W0[ii],0);
        shower->SetPb(P0[ii],0);
        shower->SetFlagb(Flag0[ii],0);
        shower->SetmcDigitIndexTop(mcDigitIndexTop0[ii],0);
        shower->SetmcDigitIndexBottom(mcDigitIndexBottom0[ii],0);


        ind = 0 +  ii*eNTM;

        xbb[ind] =  x0[ii];
        ybb[ind] =  y0[ii];
        zbb[ind] =  z0[ii];
        txbb[ind] =  tx0[ii];
        tybb[ind] =  ty0[ii];
        chi2bb[ind] =  chi20[ii];
        wbb[ind] =  W0[ii];
        Pbb[ind] =  P0[ii];
        Flagbb[ind] =  Flag0[ii];

        add[ii] = 1;
        eS->Add(shower);
    }

    EdbShowerRec* a=0;
// save the Z position of the first sheet
    /*
        if (num<10)
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/0%d_00%d.par","par",num,piece2par);
            }
            else
            {
                sprintf(nme,"%s/0%d_0%d.par","par",num,piece2par);
            }
        }
        else
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/%d_00%d.par","par",num,piece2par);
            }
            else
            {
                sprintf(nme,"%s/%d_0%d.par","par",num,piece2par);
            }
        }
        ReadAffAndZ(nme,&Zoff,&a11,&a12,&a21,&a22,&bb1,&bb2 );
        //  ZZ0[0]= z0[0];
        ZZ0[0]= Zoff;
    */
    int nntr;
    int ok;
    int max=0 ;
    int index = 0;
    double rr;

//    float Rcut = 180.;
//    float  Tcut = 0.18;
//   float Rcut = 150.;
//   float  Tcut = 0.15;
    float Rcut = eAlgoParameterConnectionDR;
    float  Tcut = eAlgoParameterConnectionDT;

    float mini[5] = {-20000.,-20000.,-1,-1,0};
    float maxi[5] = {20000.,20000.,1,1,100};
    /*
    for (int d=num+1; d<MAXPLATE+1; d++) {
        index++;
        printf( "Process plate nr (w.r.t.first BT) : %d\n",d);

        if (d<10)
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/0%d_00%d.par","par",d,piece2par);
                sprintf(fname2,"%s/0%d_00%d.cp.root","data",d,piece2);
            }
            else
            {
                sprintf(nme,"%s/0%d_0%d.par","par",d,piece2par);
                sprintf(fname2,"%s/0%d_0%d.cp.root","data",d,piece2);
            }
        }
        else
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/%d_00%d.par","par",d,piece2par);
                sprintf(fname2,"%s/%d_00%d.cp.root","data",d,piece2);
            }
            else
            {
                sprintf(nme,"%s/%d_0%d.par","par",d,piece2par);
                sprintf(fname2,"%s/%d_0%d.cp.root","data",d,piece2);
            }
        }
    */
    EdbPVRec ali;
    EdbPattern *pat=0;
    EdbPattern *pat2=0;
    EdbDataPiece Piece;
    /*
    InitPiece(Piece,fname2,nme);
    pat2 = new EdbPattern(0.,0., Piece.GetLayer(0)->Z(),100 );
    Z0 =  Piece.GetLayer(0)->Z();
    ReadPiece(Piece, *pat2);

    ReadAffAndZ(nme,&Zoff,&a11,&a12,&a21,&a22,&bb1,&bb2 );


    Z0 =  Zoff;
    ZZ0[index] =  Zoff;
    */
//   EdbPVRec  *pvr;
//   LoadEdbPVRec();
//   pvr= GetEdbPVRec();



    for (int d=0; d<pvr->Npatterns(); d++) {
        index++;
        printf( "Process plate nr (w.r.t.first BT) : %d\r",d+1);
        if (d==pvr->Npatterns()-1) printf("\n");
        pat2 = pvr->GetPattern(d);

        Z0 =  pat2->Z();
        ZZ0[index] = pat2->Z();

        int ii =0;

        TIter next(eS);
        while ((a=(EdbShowerRec*)next()))
        {

            //cout << " DEBUG:: a->GetZ0() =  "  << a->GetZ0() << "   Z0=  " <<   Z0 << endl;

            if (a->GetZ0()>Z0) ii++;
            if (a->GetZ0()>Z0) continue;
            Xss =  a->GetX0() + fabs(Z0-a->GetZ0())*a->GetTx();
            Yss =  a->GetY0() + fabs(Z0-a->GetZ0())*a->GetTy();


            mini[0] = Xss-1250.;
            mini[1] = Yss-1250.;
            mini[2] = -0.5;
            mini[3] = -0.5;
            mini[4] = 0.;
            maxi[0] = Xss+1250.;
            maxi[1] = Yss+1250.;
            maxi[2] = 0.5;
            maxi[3] = 0.5;
            maxi[4] = 100;

            pat = pat2->ExtractSubPattern(mini,maxi);
            ali.AddPattern(pat);
            //ali.Print(); cout << "This was ali.Print(); "<< endl;

            EdbSegP *ss=0;

            int incr=-1;

            for (int ip=0; ip<ali.Npatterns(); ip++) {
                incr++;
            };

            for (int i=0; i<ali.GetPattern(incr)->N(); i++)
            {
                ss = ali.GetPattern(incr)->GetSegment(i);
//         ss->PrintNice();
                ok = 0;
                if (DATA ==1) {
                    X0 = a11*ss->X() + a12*ss->Y() + bb1;
                    Y0 = a21*ss->X() + a22*ss->Y() + bb2;

                    SX0 = a11*ss->TX()+ a12*ss->TY() ;
                    SY0 = a21*ss->TX()+a22*ss->TY();
                }
                else
                {
                    X0 = ss->X() ;
                    Y0 = ss->Y() ;
                    SX0 = ss->TX();
                    SY0 = ss->TY();
                }

                if (fabs(SX0)>0.4) continue;
                if (fabs(SY0)>0.4) continue;

                //if ( a->GetTrID() == ss->MCEvt()|| ss->MCEvt()<0)
                //{
                //Delta = 20.+   fabs(0.02*fabs(Z0-a->GetZ0())); // cone with 20 mrad
                //rayon = 800.; // cylinder with 800 microns radius
                Delta = 20.+   fabs(eAlgoParameterConeAngle*fabs(Z0-a->GetZ0())); // cone with eAlgoParameterConeAngle mrad
                rayon = eAlgoParameterConeRadius; // cylinder with eAlgoParameterConeRadius microns radius

                diff = sqrt((X0-Xss)*(X0-Xss)+(Y0-Yss)*(Y0-Yss));

//           cout << "diff = " << diff << endl;
//           ss->PrintNice();

                if ( diff<rayon&&diff<Delta )
                {
                    if (index==1) max = 1;
                    if (index==2) max = 2;
                    if (index>2) max = 3;

                    for (int el=0; el<max; el++)
                    {
                        //extrapolation in previous plate
                        Xe[el] =  X0 - fabs(Z0-ZZ0[index-el-1])*SX0;
                        Ye[el] =  Y0 - fabs(Z0-ZZ0[index-el-1])*SY0;
                        nntr = add[ii];

                        for (int l2=0; l2<nntr; l2++)
                        {
                            ind = l2 + ii*eNTM;

                            if (fabs(a->GetZb(l2)-(Z0-(el+1)*1273))<500)
                            {
//                   ss->PrintNice();

                                rayon2 = sqrt((Xe[el]-a->GetXb(l2))*(Xe[el]-a->GetXb(l2))+((Ye[el]-a->GetYb(l2))*(Ye[el]-a->GetYb(l2))));
                                delta = sqrt((SX0-a->GetTXb(l2))*(SX0-a->GetTXb(l2))+((SY0-a->GetTYb(l2))*(SY0-a->GetTYb(l2))));
                                if (rayon2<Rcut&&delta<Tcut)
                                {

                                    if (ok==0)
                                    {
                                        ind = add[ii] + ii*eNTM;

                                        xbb[ind] =  X0 ;
                                        ybb[ind] =  Y0 ;
                                        zbb[ind] =  Z0 ;
                                        txbb[ind] =  SX0 ;
                                        tybb[ind] =  SY0 ;
                                        a->SetIDb(ss->ID(),add[ii]);
                                        a->SetPlateb(d,add[ii]);
                                        a->SetXb(X0,add[ii]);
                                        a->SetYb(Y0,add[ii]);
                                        a->SetZb(Z0,add[ii]);
                                        a->SetTXb(SX0,add[ii]);
                                        a->SetTYb(SY0,add[ii]);
                                        a->SetDeltarb(rayon2,add[ii]);
                                        a->SetDelthetab(delta,add[ii]);
                                        //a->SetNFilmb( index+1,add[ii]);
                                        rr =  (Z0-a->GetZ0())/1273;
                                        a->SetNFilmb(int(round(rr)+1),add[ii]);
                                        if (ss->MCEvt()==a->GetTrID()) sigbb[ind] = 1;
                                        if (ss->MCEvt()<0) sigbb[ind] = 0;
                                        chi2bb[ind] = ss->Chi2();
                                        wbb[ind] = int(ss->W());
                                        Pbb[ind] = ss->P();
                                        Flagbb[ind] = ss->Flag();

                                        a->SetmcDigitIndexTop(ss->MCEvt(),add[ii]);
                                        a->SetmcDigitIndexBottom(ss->MCTrack(),add[ii]);
                                        ok =1;

//                       ss->PrintNice();

                                        add[ii]++;
                                    }
                                    else if (ok==1&&rayon2<a->GetDeltarb(add[ii]-1))
                                    {
                                        a->SetDeltarb(rayon2,add[ii]-1);
                                        a->SetDelthetab(delta,add[ii]-1);

                                    }
                                }
                            }
                        }
                    }
                }
                //}
            }
            ii++;
        }
    }


    int  i = 0;
    int loopcount=0;
    int evtid = 0;
    int nMC=0;
    TIter next(eS);
    while ((a=(EdbShowerRec*)next()))
    {
//     cout << "while ((a=(EdbShowerRec*)next()))   loopcount= " << loopcount << endl;

        number_eventb = a->GetTrID();
        sizeb = add[i];
        a->SetSize(add[i]);
        sizeb15 = 0;
        sizeb20 = 0;
        sizeb30 = 0;
        isizeb = 0;
        purityb = 0;
        nMC=0;
        E_MC = Esim[i];

        for (int j=0; j<add[i]; j++)
        {
            if ((a->GetZb(j)>a->GetZ0()+5000)&&isizeb==0) isizeb = 1; //cross at least 4 plates
            if (a->GetNFilmb(j)<16)  sizeb15++;
            if (a->GetNFilmb(j)<21)  sizeb20++;
            if (a->GetNFilmb(j)<31)  sizeb30++;
            if (a->GetNFilmb(j)==15)  isizeb=15;
            if (a->GetNFilmb(j)==20)  isizeb=20;
            if (a->GetNFilmb(j)==30)  isizeb=30;
            if (a->GetNFilmb(j)==50)  isizeb=50;
        }
        a->SetSize15(sizeb15);
        a->SetSize20(sizeb20);
        a->SetSize30(sizeb30);

        for (int j=0; j<add[i]; j++)
        {
            idb[j] = a->GetIDb(j);
            //plateb[j] = a->GetPlateb(j);
            plateb[j] = floor(a->GetZb(j)/1273)+31;// F. Brunet : ePlate setters and getters obsolet with new plate thickness
            xb[j] = a->GetXb(j);
            yb[j] = a->GetYb(j);
            zb[j] = a->GetZb(j);
            txb[j] = a->GetTXb(j);
            tyb[j] = a->GetTYb(j);
            deltarb[j]  =  a->GetDeltarb(j);
            deltathetab[j] = a->GetDeltathetab(j);
            nfilmb[j] = a->GetNFilmb(j);
            tagprimary[j] = 0;

            mcDigitIndexTop[j] = a->GetmcDigitIndexTop(j);
            mcDigitIndexBottom[j] = a->GetmcDigitIndexBottom(j);

            //gEDBDEBUGLEVEL=3;
            if (gEDBDEBUGLEVEL>2) cout <<  xb[j] << "  " << yb[j]  << "  " << zb[j] << "  " << txb[j] << "  " << tyb[j] << "  " << sizeb << endl;

            if (j==0) tagprimary[j] = 1;
            if (j==0) chi2btkb[j] = chi2bb[j+i*eNTM];
            if (j>0)   chi2btkb[j] = chi2bb[j+i*eNTM];
            //      if (j==0) ntrace1simub[j]= 1;
            if (j==0&&a->GetTrID()>=0) ntrace1simub[j]= 1;
            if (j==0&&a->GetTrID()<0) ntrace1simub[j]= 0;

            if (j>0)  ntrace1simub[j] = sigbb[j+i*eNTM];
            if (j==0) ntrace2simub[j]= wbb[j+i*eNTM];
            if (j>0)  ntrace2simub[j] = wbb[j+i*eNTM];
            ntrace3simub[j]= Pbb[j+i*eNTM];
            ntrace4simub[j]= Flagbb[j+i*eNTM];

            //cout << " j=  " <<  j << "  i=  " <<  i  << "  eNTM= " << eNTM <<endl;
            //cout << " j+i*eNTM=  " <<  j+i*eNTM << "  chi2bb[j+i*eNTM];    =  " <<  chi2bb[j+i*eNTM]      << "  eNTM= " << eNTM <<endl;

            if (ntrace1simub[j]>=1) ++nMC;
        }
        if (sizeb>0) purityb=float(nMC)/float(sizeb);
        if (isizeb > 0)//cross at least 4 plates
        {
            showerID = evtid;
            treesaveb->Fill();
            evtid++;
        }
        i++;
    }

    eShowersN=treesaveb->GetEntries();
    treesaveb->Show(0);
    SaveResults();
    NeuralNet();


    //CalculateEnergyValues();
    // deprecated since were going to have the transition to EdbShowerAlg, which does energy...
    // New... still to check:
    //TransferTreebranchShowerTreeIntoShowerObjectArray


    // Open File. Get Tree.
    TFile *_file0 = new TFile("Shower2.root");
    TTree *tr = (TTree*)_file0->Get("treebranch");
    TransferTreebranchShowerTreeIntoShowerObjectArray(tr);
    TObjArray* RecoShowerArray = GetRecoShowerArray() ;
    // Check the number of entries:
    cout << tr->GetEntries() << endl;
    cout << GetRecoShowerArrayN() << endl;
    // Instantate ShowerAlgorithmEnergy Class
    EdbShowerAlgESimple* ShowerAlgEnergyInstance = new EdbShowerAlgESimple();
    // Run Shower Energy Algorith on all shower/tracks
    ShowerAlgEnergyInstance->DoRun(RecoShowerArray);
    // Write the values then into "Shower.root" (old Shower.root will be replaced(and backuped)!");
    ShowerAlgEnergyInstance->WriteNewRootFile("Shower2.root","treebranch");

#ifndef WIN32
    strcpy(cmd,"rm  shower1.root");
    gSystem->Exec(cmd);
    strcpy(cmd,"rm  Shower2.root");
    gSystem->Exec(cmd);
#else
    gSystem->Exec("del shower1.root");
    gSystem->Exec("del Shower2.root");
#endif

    return;
}
//-------------------------------------------------------------------
//-------------------------------------------------------------------
//void EdbShowerRec::recup(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, int *TRid, int piece2, int piece2par,float Rcut=100., float Tcut=0.05)
void EdbShowerRec::recup(int num,int MAXPLATE,  int DATA, int Ncand, double *x0, double *y0, double *z0, double *tx0, double *ty0, double* chi20, int* W0, double* P0,int* Flag0,  int *plate0, int *id0, int *TRid, double *Esim, int *mcDigitIndexTop0, int *mcDigitIndexBottom0, int piece2, int piece2par,EdbPVRec  *pvr)
{

// shower tree definition

    treesaveb = new TTree("treebranch","tree of branchtrack");
    treesaveb->Branch("number_eventb",&number_eventb,"number_eventb/I");
    treesaveb->Branch("sizeb",&sizeb,"sizeb/I");
    treesaveb->Branch("sizeb15",&sizeb15,"sizeb15/I");
    treesaveb->Branch("sizeb20",&sizeb20,"sizeb20/I");
    treesaveb->Branch("sizeb30",&sizeb30,"sizeb30/I");
    treesaveb->Branch("E_MC",&E_MC,"E_MC/F");
    treesaveb->Branch("idb",idb,"idb[sizeb]/I");
    treesaveb->Branch("plateb",plateb,"plateb[sizeb]/I");
    treesaveb->Branch("showerID",&showerID,"showerID/I");
    treesaveb->Branch("isizeb",&isizeb,"isizeb/I");
    treesaveb->Branch("xb",xb,"xb[sizeb]/F");
    treesaveb->Branch("yb",yb,"yb[sizeb]/F");
    treesaveb->Branch("zb",zb,"zb[sizeb]/F");
    treesaveb->Branch("txb",txb,"txb[sizeb]/F");
    treesaveb->Branch("tyb",tyb,"tyb[sizeb]/F");
    treesaveb->Branch("nfilmb",nfilmb,"nfilmb[sizeb]/I");
    treesaveb->Branch("ntrace1simub",ntrace1simub,"ntrace1simu[sizeb]/I");
    treesaveb->Branch("ntrace2simub",ntrace2simub,"ntrace2simu[sizeb]/I");
    treesaveb->Branch("ntrace3simub",ntrace3simub,"ntrace3simu[sizeb]/F");
    treesaveb->Branch("ntrace4simub",ntrace4simub,"ntrace4simu[sizeb]/I");
    treesaveb->Branch("chi2btkb",chi2btkb,"chi2btkb[sizeb]/F");
    treesaveb->Branch("deltarb",deltarb,"deltarb[sizeb]/F");
    treesaveb->Branch("deltathetab",deltathetab,"deltathetab[sizeb]/F");
    treesaveb->Branch("deltaxb",deltaxb,"deltaxb[sizeb]/F");
    treesaveb->Branch("deltayb",deltayb,"deltayb[sizeb]/F");
    treesaveb->Branch("tagprimary",tagprimary,"tagprimary[sizeb]/F");
    treesaveb->Branch("purityb",&purityb,"purityb/F");
    treesaveb->Branch("mcDigitIndexTop",mcDigitIndexTop,"mcDigitIndexTop[sizeb]/I");
    treesaveb->Branch("mcDigitIndexBottom",mcDigitIndexBottom,"mcDigitIndexBottom[sizeb]/I");

    fileout = new TFile("shower1.root","RECREATE");

    //additionally needed, if not given by constructor.
    // this is only a temporary solution. When giving over the constructor with
    // segments of EdbSegP class, they will contain the needed information then.
    cout << " void EdbShowerRec::recdown()  Additionally needed, if not given by constructor: " << endl;
    cout << " void EdbShowerRec::recdown()  this is only a temporary solution. When giving over the constructor with " << endl;
    cout << " void EdbShowerRec::recdown()  segments of EdbSegP class, they will contain the needed information then." << endl;
    cout << " void EdbShowerRec::recdown()  For the momemt, if Chi2, W,P,Flag of the starting BT is not given..." << endl;
    cout << " void EdbShowerRec::recdown()  It is assumed to be   ...   -9999999  ..." << endl;
//   float chi20[10000];
//   float P0[10000];
//   int W0[10000];
//   int Flag0[10000];
//   chi20[0]=-9999999;
//   W0[0]=-9999999;
//   P0[0]=-9999999;
//   Flag0[0]=-9999999;

    float Z0, SX0,SY0;
    float X0, Y0, Xss,Yss;
    float Delta, rayon,rayon2, diff, delta;
    float Xe[4], Ye[4];
    char fname2[128];

    char nme[64];
    float  a11, a12, a21, a22, bb1, bb2;
    a11=1;
    a12=0;
    a21=0;
    a22=1;
    bb1=0;
    bb2=0;
    float Zoff;
    double ZZ0[57];
    int ind;
    const Int_t eNSM = 10000;
    const Int_t eNTM = 1000;

    TArrayF  xbb(eNTM*eNSM);
    TArrayF  ybb(eNTM*eNSM);
    TArrayF  zbb(eNTM*eNSM);
    TArrayF  txbb(eNTM*eNSM);
    TArrayF  tybb(eNTM*eNSM);
    TArrayI  sigbb(eNTM*eNSM);
    TArrayF  chi2bb(eNTM*eNSM);
    TArrayI  wbb(eNTM*eNSM);
    TArrayI add(eNTM*eNSM);
    TArrayF  Pbb(eNTM*eNSM);
    TArrayI  Flagbb(eNTM*eNSM);

    for (int i=0; i<eNSM; i++)
    {
        for (int j=0; j<eNTM; j++)
        {
            ind = i*eNTM+j;
            xbb[ind];
            ybb[ind];
            zbb[ind];
            txbb[ind];
            tybb[ind];
            sigbb[ind];
            wbb[ind];
            chi2bb[ind];
            Pbb[ind];
            Flagbb[ind];
        }
    }

    TList  *eS = new TList;

    for (int ii=0; ii<Ncand; ii++)
    {
        EdbShowerRec* shower = new EdbShowerRec();

        shower->SetID(ii);
        shower->SetX0(x0[ii]);
        shower->SetY0(y0[ii]);
        shower->SetZ0(z0[ii]);
        shower->SetTx(tx0[ii]);
        shower->SetTy(ty0[ii]);
        shower->SetTrID(TRid[ii]);
        shower->SetChi20(chi20[ii]);
        shower->SetW0(W0[ii]);
        shower->SetP0(P0[ii]);
        shower->SetFlag0(Flag0[ii]);

        shower->SetDeltarb(200.,0);
        shower->SetDelthetab(0.5,0);
        shower->SetPrimary(0,0);
        shower->SetNFilmb(1,0);
        shower->SetIDb(id0[ii],0);
        shower->SetPlateb(plate0[ii],0);
        shower->SetXb(x0[ii],0);
        shower->SetYb(y0[ii],0);
        shower->SetZb(z0[ii],0);
        shower->SetTXb(tx0[ii],0);
        shower->SetTYb(ty0[ii],0);
        shower->SetChi2b(chi20[ii],0);
        shower->SetWb(W0[ii],0);
        shower->SetPb(P0[ii],0);
        shower->SetFlagb(Flag0[ii],0);
        shower->SetmcDigitIndexTop(mcDigitIndexTop0[ii],0);
        shower->SetmcDigitIndexBottom(mcDigitIndexBottom0[ii],0);

        ind = 0 +  ii*eNTM;

        xbb[ind] =  x0[ii];
        ybb[ind] =  y0[ii];
        zbb[ind] =  z0[ii];
        txbb[ind] =  tx0[ii];
        tybb[ind] =  ty0[ii];
        chi2bb[ind] =  chi20[ii];
        wbb[ind] =  W0[ii];
        Pbb[ind] =  P0[ii];
        Flagbb[ind] =  Flag0[ii];

        add[ii] = 1;
        eS->Add(shower);


    }

    EdbShowerRec* a=0;
// save the Z position of the first sheet
    if (num<10)
    {
        if (piece2par<10)
        {
            sprintf(nme,"%s/0%d_00%d.par","par",num,piece2par);
        }
        else
        {
            sprintf(nme,"%s/0%d_0%d.par","par",num,piece2par);
        }
    }
    else
    {
        if (piece2par<10)
        {
            sprintf(nme,"%s/%d_00%d.par","par",num,piece2par);
        }
        else
        {
            sprintf(nme,"%s/%d_0%d.par","par",num,piece2par);
        }
    }

    ReadAffAndZ(nme,&Zoff,&a11,&a12,&a21,&a22,&bb1,&bb2 );
    //  ZZ0[0]= z0[0];
    ZZ0[0]= Zoff;

    int nntr;
    int ok;
    int max=0 ;
    int index = 0;
    double rr;

    float Rcut = 150.;
    float  Tcut = 0.15;

    float mini[5] = {-20000.,-20000.,-1,-1,0};
    float maxi[5] = {20000.,20000.,1,1,100};

//    for (int d=num+1; d<MAXPLATE+1; d++){
//  if (DOWN==0) break;
    for (int d=num-1; d>num-MAXPLATE; d--) {
        index++;
        //sprintf( "Process plate number: %d/%d\n",index+1,MAXPLATE-PLATE);

        //if(d<10) {sprintf(nme,"%s/0%d_00%d.par","par",d,piece2par);}
        //else {sprintf(nme,"%s/%d_00%d.par","par",d,piece2par);}

        if (d<10) {
            sprintf(fname2,"%s/0%d_00%d.cp.root","data",d,piece2);
        }
        else {
            sprintf(fname2,"%s/%d_00%d.cp.root","data",d,piece2);
        }

        if (d<10)
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/0%d_00%d.par","par",d,piece2par);
                sprintf(fname2,"%s/0%d_00%d.cp.root","data",d,piece2);
            }
            else
            {
                sprintf(nme,"%s/0%d_0%d.par","par",d,piece2par);
                sprintf(fname2,"%s/0%d_0%d.cp.root","data",d,piece2);
            }
        }
        else
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/%d_00%d.par","par",d,piece2par);
                sprintf(fname2,"%s/%d_00%d.cp.root","data",d,piece2);
            }
            else
            {
                sprintf(nme,"%s/%d_0%d.par","par",d,piece2par);
                sprintf(fname2,"%s/%d_0%d.cp.root","data",d,piece2);
            }
        }

        EdbPVRec ali;
        EdbPattern *pat=0;
        EdbPattern *pat2=0;
        EdbDataPiece Piece;

        InitPiece(Piece,fname2,nme);
        pat2 = new EdbPattern(0.,0., Piece.GetLayer(0)->Z(),100 );
        Z0 =  Piece.GetLayer(0)->Z();
        ReadPiece(Piece, *pat2);

        ReadAffAndZ(nme,&Zoff,&a11,&a12,&a21,&a22,&bb1,&bb2 );
        Z0 =  Zoff;
        ZZ0[index] =  Zoff;

        int ii =0;
        TIter next(eS);
        while ((a=(EdbShowerRec*)next()))
        {
            if (a->GetZ0()>Z0) ii++;
            if (a->GetZ0()>Z0) continue;
            Xss =  a->GetX0() + fabs(Z0-a->GetZ0())*a->GetTx();
            Yss =  a->GetY0() + fabs(Z0-a->GetZ0())*a->GetTy();

            mini[0] = Xss-1250.;
            mini[1] = Yss-1250.;
            mini[2] = -0.5;
            mini[3] = -0.5;
            mini[4] = 0.;
            maxi[0] = Xss+1250.;
            maxi[1] = Yss+1250.;
            maxi[2] = 0.5;
            maxi[3] = 0.5;
            maxi[4] = 100.;

            pat = pat2->ExtractSubPattern(mini,maxi);
            ali.AddPattern(pat);

            EdbSegP *ss=0;

            int incr=-1;

            for (int ip=0; ip<ali.Npatterns(); ip++) {
                incr++;
            };

            for (int i=0; i<ali.GetPattern(incr)->N(); i++)
            {
                ss = ali.GetPattern(incr)->GetSegment(i);
                ok = 0;
                if (DATA ==1) {
                    X0 = a11*ss->X() + a12*ss->Y() + bb1;
                    Y0 = a21*ss->X() + a22*ss->Y() + bb2;

                    SX0 = a11*ss->TX()+ a12*ss->TY() ;
                    SY0 = a21*ss->TX()+a22*ss->TY();
                }
                else
                {
                    X0 = ss->X() ;
                    Y0 = ss->Y() ;
                    SX0 = ss->TX();
                    SY0 = ss->TY();
                }

                if (fabs(SX0)>0.4) continue;
                if (fabs(SY0)>0.4) continue;

                //if ( a->GetTrID() == ss->MCEvt()|| ss->MCEvt()<0)
                //{
                Delta = 20.+   fabs(0.02*fabs(Z0-a->GetZ0())); // cone with 20 mrad
                rayon = 800.; // cylinder with 800 microns radius
                diff = sqrt((X0-Xss)*(X0-Xss)+(Y0-Yss)*(Y0-Yss));

                if ( diff<rayon&&diff<Delta )
                {
                    if (index==1) max = 1;
                    if (index==2) max = 2;
                    if (index>2) max = 3;

                    for (int el=0; el<max; el++)
                    {
                        //extrapolation in previous plate
                        Xe[el] =  X0 - fabs(Z0-ZZ0[index-el-1])*SX0;
                        Ye[el] =  Y0 - fabs(Z0-ZZ0[index-el-1])*SY0;
                        nntr = add[ii];

                        for (int l2=0; l2<nntr; l2++)
                        {
                            ind = l2 + ii*eNTM;
                            if (fabs(zbb[ind]-(Z0-(el+1)*1273))<500)
                            {
                                rayon2 = sqrt((Xe[el]-xbb[ind])*(Xe[el]-xbb[ind])+((Ye[el]-ybb[ind])*(Ye[el]-ybb[ind])));
                                delta = sqrt((SX0-txbb[ind])*(SX0-txbb[ind])+((SY0-tybb[ind])*(SY0-tybb[ind])));
                                if (rayon2<Rcut&&delta<Tcut)
                                {
                                    if (ok==0)
                                    {
                                        ind = add[ii] + ii*eNTM;

                                        xbb[ind] =  X0 ;
                                        ybb[ind] =  Y0 ;
                                        zbb[ind] =  Z0 ;
                                        txbb[ind] =  SX0 ;
                                        tybb[ind] =  SY0 ;

                                        a->SetIDb(ss->ID(),add[ii]);
                                        a->SetPlateb(d,add[ii]);
                                        a->SetXb(X0,add[ii]);
                                        a->SetYb(Y0,add[ii]);
                                        a->SetZb(Z0,add[ii]);
                                        a->SetTXb(SX0,add[ii]);
                                        a->SetTYb(SY0,add[ii]);
                                        a->SetDeltarb(rayon2,add[ii]);
                                        a->SetDelthetab(delta,add[ii]);
                                        //a->SetNFilmb( index+1,add[ii]);
                                        rr =  (Z0-a->GetZ0())/1273;
                                        a->SetNFilmb(int(round(rr)+1),add[ii]);

                                        if (ss->MCEvt()==a->GetTrID()) sigbb[ind] = 1;
                                        if (ss->MCEvt()<0) sigbb[ind] = 0;
                                        chi2bb[ind] = ss->Chi2();
                                        wbb[ind] = int(ss->W());
                                        Pbb[ind] = ss->P();
                                        Flagbb[ind] = ss->Flag();

                                        a->SetmcDigitIndexTop(ss->MCEvt(),add[ii]);
                                        a->SetmcDigitIndexBottom(ss->MCTrack(),add[ii]);
                                        ok =1;
                                        add[ii]++;
                                    }
                                    else if (ok==1&&rayon2<a->GetDeltarb(add[ii]-1))
                                    {
                                        a->SetDeltarb(rayon2,add[ii]-1);
                                        a->SetDelthetab(delta,add[ii]-1);
                                    }
                                }
                            }
                        }
                    }
                }
                //}
            }
            ii++;
        }
    }

    int  i = 0;
    int evtid = 0;
    int nMC=0;
    TIter next(eS);
    while ((a=(EdbShowerRec*)next()))
    {
        number_eventb = a->GetID();
        sizeb = add[i];
        a->SetSize(add[i]);
        sizeb15 = 0;
        sizeb20 = 0;
        sizeb30 = 0;
        isizeb = 0;
        purityb = 0;
        nMC=0;
        E_MC = Esim[i];

        for (int j=0; j<add[i]; j++)
        {
            if ((a->GetZb(j)>a->GetZ0()+5000)&&isizeb==0) isizeb = 1; //cross at least 4 plates
            if (a->GetNFilmb(j)<16)  sizeb15++;
            if (a->GetNFilmb(j)<21)  sizeb20++;
            if (a->GetNFilmb(j)<31)  sizeb30++;
            if (a->GetNFilmb(j)==15)  isizeb=15;
            if (a->GetNFilmb(j)==20)  isizeb=20;
            if (a->GetNFilmb(j)==30)  isizeb=30;
            if (a->GetNFilmb(j)==50)  isizeb=50;
        }
        a->SetSize15(sizeb15);
        a->SetSize20(sizeb20);
        a->SetSize30(sizeb30);

        for (int j=0; j<add[i]; j++)
        {
            idb[j] = a->GetIDb(j);
            plateb[j] = a->GetPlateb(j);
            xb[j] = a->GetXb(j);
            yb[j] = a->GetYb(j);
            zb[j] = a->GetZb(j);
            txb[j] = a->GetTXb(j);
            tyb[j] = a->GetTYb(j);
            deltarb[j]  =  a->GetDeltarb(j);
            deltathetab[j] = a->GetDeltathetab(j);
            nfilmb[j] = a->GetNFilmb(j);
            tagprimary[j] = 0;

            mcDigitIndexTop[j] = a->GetmcDigitIndexTop(j);
            mcDigitIndexBottom[j] = a->GetmcDigitIndexBottom(j);

//      printf("txb[j]  %f\n", a->GetTXb(j));

            if (gEDBDEBUGLEVEL>2) cout <<  xb[j] << "  " << yb[j]  << "  " << zb[j] << "  " << txb[j] << "  " << tyb[j] << "  " << sizeb << endl;

            if (j==0) tagprimary[j] = 1;
            if (j==0) chi2btkb[j] = chi2bb[j+i*eNTM];
            if (j>0)   chi2btkb[j] = chi2bb[j+i*eNTM];
            //      if (j==0) ntrace1simub[j]= 1;
            if (j==0&&a->GetTrID()>=0) ntrace1simub[j]= 1;
            if (j==0&&a->GetTrID()<0) ntrace1simub[j]= 0;

            if (j>0)  ntrace1simub[j] = sigbb[j+i*eNTM];
            if (j==0) ntrace2simub[j]= wbb[j+i*eNTM];
            if (j>0)  ntrace2simub[j] = wbb[j+i*eNTM];
            ntrace3simub[j]= Pbb[j+i*eNTM];
            ntrace4simub[j]= Flagbb[j+i*eNTM];

            //cout << " j=  " <<  j << "  i=  " <<  i  << "  eNTM= " << eNTM <<endl;
            //cout << " j+i*eNTM=  " <<  j+i*eNTM << "  chi2bb[j+i*eNTM];    =  " <<  chi2bb[j+i*eNTM]      << "  eNTM= " << eNTM <<endl;

            if (ntrace1simub[j]>=1) ++nMC;
        }
        if (sizeb>0) purityb=float(nMC)/float(sizeb);
        if (isizeb >0) //cross at least 4 plates
        {
            // number_eventb = evtid;
            showerID = evtid;
            treesaveb->Fill();
            evtid++;
        }
        i++;
    }

    eShowersN=treesaveb->GetEntries();
    //treesaveb->Show(0);

    SaveResults();
    NeuralNet();
    // Attention: NeuralNet() ->Gives Shower2.root as output!

    // Open File. Get Tree.
    TFile *_file0 = new TFile("Shower2.root");
    TTree *tr = (TTree*)_file0->Get("treebranch");
    TransferTreebranchShowerTreeIntoShowerObjectArray(tr);
    TObjArray* RecoShowerArray = GetRecoShowerArray() ;
    // Check the number of entries:
    cout << tr->GetEntries() << endl;
    cout << GetRecoShowerArrayN() << endl;
    // Instantate ShowerAlgorithmEnergy Class
    EdbShowerAlgESimple* ShowerAlgEnergyInstance = new EdbShowerAlgESimple();
    // Run Shower Energy Algorith on all shower/tracks
    ShowerAlgEnergyInstance->DoRun(RecoShowerArray);
    cout << "// Write the values then into Shower.root (old Shower.root will be replaced(and backuped)!);" << endl;

    ShowerAlgEnergyInstance->WriteNewRootFile();


#ifndef WIN32
    strcpy(cmd,"rm  shower1.root");
    gSystem->Exec(cmd);
    strcpy(cmd,"rm  Shower2.root");
    gSystem->Exec(cmd);
#else
    gSystem->Exec("del shower1.root");
    gSystem->Exec("del Shower2.root");
#endif

    return;
}

//-------------------------------------------------------------------
//-------------------------------------------------------------------

//void EdbShowerRec::remove(char *shfname, char *def,float Rcut=100., float Tcut=0.05)
void EdbShowerRec::remove(char *shfname, char *def, int MAXPLATE, int piece2par)
{

    treesaveb = new TTree("treebranch","tree of branchtrack");
    treesaveb->Branch("number_eventb",&number_eventb,"number_eventb/I");
    treesaveb->Branch("sizeb",&sizeb,"sizeb/I");
    treesaveb->Branch("sizeb15",&sizeb15,"sizeb15/I");
    treesaveb->Branch("sizeb20",&sizeb20,"sizeb20/I");
    treesaveb->Branch("sizeb30",&sizeb30,"sizeb30/I");
    treesaveb->Branch("E_MC",&E_MC,"E_MC/F");
    treesaveb->Branch("isizeb",&isizeb,"isizeb/I");
    treesaveb->Branch("xb",xb,"xb[sizeb]/F");
    treesaveb->Branch("yb",yb,"yb[sizeb]/F");
    treesaveb->Branch("zb",zb,"zb[sizeb]/F");
    treesaveb->Branch("txb",txb,"txb[sizeb]/F");
    treesaveb->Branch("tyb",tyb,"tyb[sizeb]/F");
    treesaveb->Branch("nfilmb",nfilmb,"nfilmb[sizeb]/I");
    treesaveb->Branch("ntrace1simub",ntrace1simub,"ntrace1simu[sizeb]/I");
    treesaveb->Branch("ntrace2simub",ntrace2simub,"ntrace2simu[sizeb]/I");
    treesaveb->Branch("chi2btkb",chi2btkb,"chi2btkb[sizeb]/F");
    treesaveb->Branch("deltarb",deltarb,"deltarb[sizeb]/F");
    treesaveb->Branch("deltathetab",deltathetab,"deltathetab[sizeb]/F");
    treesaveb->Branch("deltaxb",deltaxb,"deltaxb[sizeb]/F");
    treesaveb->Branch("deltayb",deltayb,"deltayb[sizeb]/F");
    treesaveb->Branch("tagprimary",tagprimary,"tagprimary[sizeb]/F");
    fileout = new TFile("shower.root","RECREATE");

//    vert(shfname,def,Rcut,Tcut);
    vert(shfname,def,MAXPLATE,piece2par);

    SaveResults();

}
//-------------------------------------------------------------------
//-------------------------------------------------------------------
// void EdbShowerRec::vert(const char *name2,const char *def,float Rcut=100., float Tcut=0.05)
void EdbShowerRec::vert(const char *name2,const char *def, int MAXPLATE, int piece2par)
{
    double ZZ0[57];
    char nme[64];
    int ind;
    float  a11, a12, a21, a22, bb1, bb2;
    float Zoff;
    const Int_t eNSM = 10000;
    const Int_t eNTM = 1000;

    TArrayF  xbb(eNTM*eNSM);
    TArrayF  ybb(eNTM*eNSM);
    TArrayF  zbb(eNTM*eNSM);
    TArrayF  txbb(eNTM*eNSM);
    TArrayF  tybb(eNTM*eNSM);
    TArrayF  deltarbb(eNTM*eNSM);
    TArrayF  deltatbb(eNTM*eNSM);
    TArrayI  sigbb(eNTM*eNSM);
    TArrayI  nfilmbb(eNTM*eNSM);
    TArrayI Flag(eNTM*eNSM);
    TArrayI add(eNTM*eNSM);
    TArrayI EvNum(eNSM);

    for (int i=0; i<eNSM; i++)
    {
        for (int j=0; j<eNTM; j++)
        {
            ind = i*eNTM+j;
            xbb[ind];
            ybb[ind];
            zbb[ind];
            txbb[ind];
            tybb[ind];
            deltarbb[ind];
            deltatbb[ind];
            sigbb[ind];
            nfilmbb[ind];
            Flag[ind];

        }
    }

    file = new TFile(name2 ,"READ");
    TTree *treebranch_e = (TTree *) file->Get("treebranch");
    treebranch_e->SetBranchAddress("number_eventb", &number_eventb);
    treebranch_e->SetBranchAddress("nfilmb", nfilmb);
    treebranch_e->SetBranchAddress("sizeb", &sizeb);
    treebranch_e->SetBranchAddress("tagprimary", tagprimary);
    treebranch_e->SetBranchAddress("txb", txb);
    treebranch_e->SetBranchAddress("tyb", tyb);
    treebranch_e->SetBranchAddress("xb", xb);
    treebranch_e->SetBranchAddress("yb",yb);
    treebranch_e->SetBranchAddress("zb",zb);
    treebranch_e->SetBranchAddress("chi2btkb",chi2btkb);

    int nentries_e = int(treebranch_e->GetEntries());
    TList  *eS = new TList;

    for (Int_t ii = 0; ii<nentries_e; ii++)
    {
        treebranch_e->GetEntry(ii);

        EdbShowerRec* shower = new EdbShowerRec();
        shower->SetID(ii);
        shower->SetX0(xb[0]);
        shower->SetY0(yb[0]);
        shower->SetZ0(zb[0]);
        shower->SetTx(txb[0]);
        shower->SetTy(tyb[0]);
        shower->SetTrID(number_eventb);

        shower->SetDeltarb(200.,0);
        shower->SetDelthetab(0.5,0);
        shower->SetPrimary(0,0);
        shower->SetNFilmb(1,0);
        shower->SetSize(sizeb);

        shower->SetXb(xb[0],0);
        shower->SetYb(yb[0],0);
        shower->SetZb(zb[0],0);
        shower->SetTXb(txb[0],0);
        shower->SetTYb(tyb[0],0);

        for (Int_t ibtke = 1; ibtke < sizeb; ibtke++)
        {
            shower->SetXb(xb[ibtke],ibtke);
            shower->SetYb(yb[ibtke],ibtke);
            shower->SetZb(zb[ibtke],ibtke);
            shower->SetTXb(txb[ibtke],ibtke);
            shower->SetTYb(tyb[ibtke],ibtke);
            shower->SetNFilmb(nfilmb[ibtke],ibtke);
        }
        eS->Add(shower);
    }

    EdbShowerRec* a=0;

    //================================================================
    float Z0, Z1, TX0,TY0  ;
    float X0, Y0, X1, Y1, Xss,Yss;
    float Delta, rayon,rayon2, diff, delta;
    float diff2, Delta2;
    float Xe[4], Ye[4];
    int  nseg;

    for (int jj=0; jj<eNSM; jj++)
    {
        for (int ii=0; ii<eNTM; ii++)
        {
            ind = jj*eNTM +ii;
            Flag[ind] = 0;
        }
    }

    EdbTrackP *tr=0;
//read the linked_tracks.root tree
    initproc(def, 100,"nseg>3");

    //printf("nseg>3&&s.eChi2<(s.eW*0.3-3.)\n");
    gAli->FillCell(30,30,0.009,0.009);

    int ntr = gAli->eTracks->GetEntries();
    printf("ntr = %d\n",ntr);

    EdbSegP *s = new EdbSegP();
    for (int i=0; i<ntr; i++) {
        tr = (EdbTrackP*)(gAli->eTracks->At(i));
        if (i%1000==0)  printf( "Processed tracks %d/%d\n",i,ntr);
        nseg = tr->N();
        if (nseg<4) continue;  //do not take into account tracks with less than 4 segments

        s = tr->GetSegment(0);
        X1 = s->X();
        Y1 = s->Y();
        Z1 = s->Z();
        for (int l=0; l<nseg; l++)
        {
            s = tr->GetSegment(l);
            Z0 = s->Z();
            TX0 =  s->TX();
            TY0 = s->TY();
            X0 =  s->X();
            Y0 =  s->Y();
// loop over the shower
            int  ii=0;
            TIter next(eS);
            while ((a=(EdbShowerRec*)next()))
            {
                for (Int_t ibtke = 1; ibtke < a->GetSize(); ibtke++)
                {
                    if ( a->GetXb(ibtke)==X0&&a->GetYb(ibtke)==Y0&&a->GetZb(ibtke)==Z0 )
                    {
                        rayon = 800.; // cylinder with 400 microns radius

                        // Extrapolation of the first shower basetrack at the Z basetrack position of the first basetrack of the track
                        Xss =  a->GetX0() + (Z1-a->GetZ0())*a->GetTx();
                        Yss =  a->GetY0() + (Z1-a->GetZ0())*a->GetTy();

                        // check if the first  basetrack of the track is inside the cone
                        diff2 = sqrt((X1-Xss)*(X1-Xss)+(Y1-Yss)*(Y1-Yss));
                        Delta2 = 20. +   fabs(0.02*(Z1-a->GetZ0())); // cone with 20 mrad

                        // Extrapolation of the first shower basetrack at the current Z basetrack position
                        Xss =  a->GetX0() + (Z0-a->GetZ0())*a->GetTx();
                        Yss =  a->GetY0() + (Z0-a->GetZ0())*a->GetTy();
                        // check if the current basetrack is inside the cone
                        diff = sqrt((X0-Xss)*(X0-Xss)+(Y0-Yss)*(Y0-Yss));
                        Delta = 20. +   fabs(0.02*(Z0-a->GetZ0())); // cone with 20 mrad

                        if ( (diff<rayon&&diff<Delta)&&(diff2>Delta2||diff2>rayon) )  // (the BT must be in the cone)&&(the first BT of the track must be outside the cone)
                        {
                            printf( "%f %f %f %f %d %f\n", Z1, Delta2,diff2,Z0,ii,chi2btkb[ibtke]);
                            Flag[ibtke+ii*eNTM] = 1;
                        }
                    }
                }
                ii++;
            }
        }
    }
    file->Close();
    gAli->Delete();

// get Z value for each plate
    ZZ0[0] = 0.;
    for (Int_t iii = 1; iii<MAXPLATE; iii++)
    {
        if (iii<10)
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/0%d_00%d.par","par",iii,piece2par);
            }
            else
            {
                sprintf(nme,"%s/0%d_0%d.par","par",iii,piece2par);
            }
        }
        else
        {
            if (piece2par<10)
            {
                sprintf(nme,"%s/%d_00%d.par","par",iii,piece2par);
            }
            else
            {
                sprintf(nme,"%s/%d_0%d.par","par",iii,piece2par);
            }
        }
        //  if(iii<10) {sprintf(nme,"%s/0%d_00%d.par","par",iii,piece2par);}
        //  else {sprintf(nme,"%s/%d_00%d.par","par",iii,piece2par);}

        ReadAffAndZ(nme,&Zoff,&a11,&a12,&a21,&a22,&bb1,&bb2);
        ZZ0[iii] = Zoff;

    }

    float zpos,rr;
    int nntr, d, max =0;
    int ok;
//    float Rcut = 100.;
//    float  Tcut = 0.05;
    float Rcut = 150.;
    float  Tcut = 0.15;

    int ige =0;

    TIter next(eS);
    while ((a=(EdbShowerRec*)next()))
    {
        ind = 0 + ige*eNTM;
        add[ige] = 1;
        EvNum[ige] = a->GetTrID();
        xbb[ind] =  a->GetX0();
        ybb[ind] =  a->GetY0();
        zbb[ind] =  a->GetZ0();
        txbb[ind] = a->GetTx();
        tybb[ind] = a->GetTy();
        deltarbb[ind]= 200. ;
        deltatbb[ind]= 0.5;

        nfilmbb[ind] = 1;
        for ( int ibtke = 1; ibtke < a->GetSize() ; ibtke++)
        {
            ok = 0;
            Xss =  a->GetX0() + (a->GetZb(ibtke)-a->GetZ0())*a->GetTx();
            Yss =  a->GetY0() + (a->GetZb(ibtke)-a->GetZ0())*a->GetTy();

            Delta = 20. + fabs(0.02*(a->GetZb(ibtke)-a->GetZ0())); // cone with 20
            rayon = 800.; // cylinder with 400 microns radius
            diff = sqrt((a->GetXb(ibtke)-Xss)*(a->GetXb(ibtke)-Xss)+(a->GetYb(ibtke)-Yss)*(a->GetYb(ibtke)-Yss));

            if ( diff<rayon&&diff<Delta&&Flag[ibtke+ige*eNTM]==0 )  //BT must be in the cone and does not belong to a track which starts outside the cone volume
            {

                d = a->GetNFilmb(ibtke);
                zpos = a->GetZb(ibtke);
                if (zpos>1000&&zpos<1500.)
                {
                    max = 1;
                } else if (zpos<3000)
                {
                    max = 2;
                }
                else
                {
                    max = 3;
                }
                for (int el=0; el<max; el++)
                {

                    rr =  a->GetZb(ibtke)/1273;
                    d = int(round(rr));

                    Xe[el] =  a->GetXb(ibtke) -(a->GetZb(ibtke)-ZZ0[d-el])*a->GetTXb(ibtke);
                    Ye[el] =  a->GetYb(ibtke) -(a->GetZb(ibtke)-ZZ0[d-el])*a->GetTYb(ibtke);
                    nntr = add[ige];

                    for (int l2=0; l2<nntr; l2++) {

                        if (fabs(a->GetZb(l2)-(a->GetZb(ibtke)-(el+1)*1273))<500)
                        {

                            rayon2 = sqrt((Xe[el]-a->GetXb(l2))*(Xe[el]-a->GetXb(l2))+((Ye[el]-a->GetYb(l2))*(Ye[el]-a->GetYb(l2))));
                            delta = sqrt((a->GetTXb(ibtke)-a->GetTXb(l2))*(a->GetTXb(ibtke)-a->GetTXb(l2))+((a->GetTYb(ibtke)-a->GetTYb(l2))*(a->GetTYb(ibtke)-a->GetTYb(l2))));


                            if (rayon2<Rcut&&delta<Tcut&&Flag[ige+l2*eNSM]==0)  // connection criteria with  a BT which does not belong to a track which starts outside the cone volume
                            {
                                if (ok==0)
                                {
                                    //ind = ige + add[ige]*eNSM;
                                    ind = add[ige] + ige*eNTM;
                                    xbb[ind] =   a->GetXb(ibtke);
                                    ybb[ind] =   a->GetYb(ibtke);
                                    zbb[ind] =   a->GetZb(ibtke);
                                    txbb[ind] =  a->GetTXb(ibtke);
                                    tybb[ind] =  a->GetTYb(ibtke);
                                    deltarbb[ind] =  rayon2 ;
                                    deltatbb[ind] =  delta ;
                                    EvNum[ige]= a->GetTrID();

                                    nfilmbb[ind] =  a->GetNFilmb(ibtke);

                                    sigbb[ind] = int(chi2btkb[ibtke]);
                                    ok =1;
                                    add[ige]++;

                                }
                                else if (ok==1&&rayon2<deltarbb[(add[ige]-1)+ige*eNTM])
                                {
                                    ind = (add[ige]-1) + ige*eNTM;
                                    deltarbb[ind] =  rayon2 ;
                                    deltatbb[ind] =  delta ;
                                }
                            }
                        }
                    } //for (int l2=0;l2<nntr;l2++)
                } // for (int el=0;el<max;el++)
            }//if ( diff<rayon&&diff<Delta )
        }   // for ( ibtke = 1;
        ige++;
    } // while

    for (int i=0; i<nentries_e; i++)
    {
//  number_eventb = i;
        number_eventb = EvNum[i];

        sizeb = add[i];
        sizeb15 = 0;
        sizeb20 = 0;
        sizeb30 = 0;

        isizeb = 0;

        for (int j=0; j<sizeb; j++)
        {
            ind = j + i*eNTM;
            if (zbb[ind]>zb[0]+5000) isizeb = 1;
            if (nfilmbb[ind]<16)  sizeb15++;
            if (nfilmbb[ind]<21)  sizeb20++;
            if (nfilmbb[ind]<31)  sizeb30++;
//  }

//  for (int j=0;j<sizeb; j++)
//  {
            xb[j] = xbb[ind];
            yb[j] = ybb[ind];
            zb[j] = zbb[ind];
            txb[j] = txbb[ind];
            tyb[j] = tybb[ind];
            deltarb[j]  = deltarbb[ind] ;
            deltathetab[j] = deltatbb[ind];

            nfilmb[j] = nfilmbb[ind];
            tagprimary[j] = 0;
            if (j==0) tagprimary[j] = 1;
//      if (j==0) chi2btkb[j] = 1;
//      if (j>0)  chi2btkb[j] = sigbb[ind];
            if (j==0) ntrace1simub[j]= 1;
            if (j>0)  ntrace1simub[j] = sigbb[j+i*eNTM];
        }
        if (isizeb > 0)treesaveb->Fill();
    }

}
//-------------------------------------------------------------------
//-------------------------------------------------------------------

void EdbShowerRec::ReadAffAndZ(char *fname, Float_t *pZoffs,  Float_t *a11, Float_t *a12,  Float_t *a21, Float_t *a22, Float_t *b1, Float_t *b2)
{
    char            buf[256];
    char            key[256];
    Int_t id;
    Float_t fTranslateX,fTranslateY;

    FILE *fp=fopen(fname,"r");
    if (fp==NULL)   {
        printf("ERROR open file: %s n", fname);
        return;
    } else
        // printf( "Read affine transformation from file: %s\n", fname );
        while ( fgets(buf,256,fp)!=NULL ) {

            for ( int i=0; i<(int)strlen(buf); i++ )
                if ( buf[i]=='#' )  {
                    buf[i]='0';                       // cut out comments starting from #
                    break;
                }

            if ( sscanf(buf,"%s",key)!=1 ) continue;
            if      ( !strcmp(key,"INCLUDE")   )
            {
                //        sscanf(buf+strlen(key),"%s",name);
                //ReadAffAndZ(name,pZoffs);
            }
            else if ( !strcmp(key,"AFFXY")  )
            {
                sscanf(buf+strlen(key),"%d %f %f %f %f %f %f",&id,a11,a12,a21,a22,b1,b2);
            }
            else if ( !strcmp(key,"ZLAYER")   )
            {
                sscanf(buf+strlen(key),"%d %f",&id,pZoffs);
            }
            else if ( !strcmp(key,"VOLUME0")   )
            {
                sscanf(buf+strlen(key),"%f %f",&fTranslateX,&fTranslateY);
            }
        }
    fclose(fp);
}

//-------------------------------------------------------------------

int EdbShowerRec::ReadPiece( EdbDataPiece &piece, EdbPattern &pat )
{
    //assuming that piece was initialised correctly

    if (!piece.InitCouplesTree("READ")) return 0;
    piece.GetCPData_new( &pat,0,0,0 );
    pat.SetSegmentsZ();
    pat.Transform(    piece.GetLayer(0)->GetAffineXY()   );
    pat.TransformA(   piece.GetLayer(0)->GetAffineTXTY() );
    pat.TransformShr( piece.GetLayer(0)->Shr()  );
    return 1;
}

//-------------------------------------------------------------------

void EdbShowerRec::NeuralNet() {

    float offset;  // Z position of last (first) plate in the data set (usually = 0.)
    offset = 0.;

//ouput file
    fileout3 = new TFile("Shower2.root","RECREATE");
    //definition of treebranch
    treesaveb3 = new TTree("treebranch","tree of branchtrack");
    treesaveb3->Branch("number_eventb",&number_eventb,"number_eventb/I");
    treesaveb3->Branch("sizeb",&sizeb,"sizeb/I");
    treesaveb3->Branch("sizeb15",&sizeb15,"sizeb15/I");
    treesaveb3->Branch("sizeb20",&sizeb20,"sizeb20/I");
    treesaveb3->Branch("sizeb30",&sizeb30,"sizeb30/I");
    treesaveb3->Branch("output15",&output15,"output15/F");
    treesaveb3->Branch("output20",&output20,"output20/F");
    treesaveb3->Branch("output30",&output30,"output30/F");
    treesaveb3->Branch("output50",&output50,"output50/F");
    treesaveb3->Branch("eProb90",&eProb90,"eProb90/I");
    treesaveb3->Branch("eProb1",&eProb1,"eProb1/I");
    treesaveb3->Branch("E_MC",&E_MC,"E_MC/F");
    treesaveb3->Branch("idb",idb,"idb[sizeb]/I");
    treesaveb3->Branch("plateb",plateb,"plateb[sizeb]/I");
    treesaveb3->Branch("showerID",&showerID,"showerID/I");
    treesaveb3->Branch("isizeb",&isizeb,"isizeb/I");
    treesaveb3->Branch("xb",xb,"xb[sizeb]/F");
    treesaveb3->Branch("yb",yb,"yb[sizeb]/F");
    treesaveb3->Branch("zb",zb,"zb[sizeb]/F");
    treesaveb3->Branch("txb",txb,"txb[sizeb]/F");
    treesaveb3->Branch("tyb",tyb,"tyb[sizeb]/F");
    treesaveb3->Branch("nfilmb",nfilmb,"nfilmb[sizeb]/I");
    treesaveb3->Branch("ntrace1simub",ntrace1simub,"ntrace1simu[sizeb]/I");
    treesaveb3->Branch("ntrace2simub",ntrace2simub,"ntrace2simu[sizeb]/I");
    treesaveb3->Branch("ntrace4simub",ntrace4simub,"ntrace4simu[sizeb]/I");
    treesaveb3->Branch("ntrace3simub",ntrace3simub,"ntrace3simu[sizeb]/F");
    treesaveb3->Branch("chi2btkb",chi2btkb,"chi2btkb[sizeb]/F");
    treesaveb3->Branch("deltarb",deltarb,"deltarb[sizeb]/F");
    treesaveb3->Branch("deltathetab",deltathetab,"deltathetab[sizeb]/F");
    treesaveb3->Branch("deltaxb",deltaxb,"deltaxb[sizeb]/F");
    treesaveb3->Branch("deltayb",deltayb,"deltayb[sizeb]/F");
    treesaveb3->Branch("tagprimary",tagprimary,"tagprimary[sizeb]/F");
    treesaveb3->Branch("purityb",&purityb,"purityb/F");
    treesaveb3->Branch("trackdensb",&trackdensb,"trackdensb/F");
    treesaveb3->Branch("mcDigitIndexTop",mcDigitIndexTop,"mcDigitIndexTop[sizeb]/I");
    treesaveb3->Branch("mcDigitIndexBottom",mcDigitIndexBottom,"mcDigitIndexBottom[sizeb]/I");

    treesaveb3->Branch("Energy",&EnergyCorrectedb,"EnergyCorrectedb/F");
    treesaveb3->Branch("EnergyUnCorrected",&EnergyUnCorrectedb,"EnergyUnCorrectedb/F");
    treesaveb3->Branch("EnergySigma",&EnergySigmaCorrectedb,"EnergySigmaCorrectedb/F");
    treesaveb3->Branch("EnergySigmaUnCorrected",&EnergySigmaUnCorrectedb,"EnergySigmaUnCorrectedb/F");



    char fname_e[128];
    Float_t Rcut30, Tcut30;
    Float_t Rcut20, Tcut20;
    Float_t Rcut15, Tcut15;
    Float_t Rcut50, Tcut50;

    int flag15;
    int flag20;
    int flag30;
    int flag50;

    Double_t params[53];

    Float_t ntrke50;
    Float_t ntrke30;
    Float_t ntrke20;
    Float_t ntrke15;

    TTree *varia = new TTree("MonteCarlo", "Filtered Monte Carlo Events");
    varia->Branch("bin3", &bin3, "bin3/F");
    varia->Branch("bin4", &bin4, "bin4/F");
    varia->Branch("bin5", &bin5, "bin5/F");
    varia->Branch("bin6", &bin6, "bin6/F");
    varia->Branch("bin7", &bin7, "bin7/F");
    varia->Branch("bin8", &bin8, "bin8/F");
    varia->Branch("bin9", &bin9, "bin9/F");
    varia->Branch("bin10", &bin10, "bin10/F");
    varia->Branch("bin11", &bin11, "bin11/F");
    varia->Branch("bin12", &bin12, "bin12/F");
    varia->Branch("bin13", &bin13, "bin13/F");
    varia->Branch("bin14", &bin14, "bin14/F");
    varia->Branch("bin15", &bin15, "bin15/F");
    varia->Branch("bin16", &bin16, "bin16/F");
    varia->Branch("bin17", &bin17, "bin17/F");
    varia->Branch("bin18", &bin18, "bin18/F");
    varia->Branch("bin19", &bin19, "bin19/F");
    varia->Branch("bin20", &bin20, "bin20/F");
    varia->Branch("bin21", &bin21, "bin21/F");
    varia->Branch("bin22", &bin22, "bin22/F");
    varia->Branch("bin23", &bin23, "bin23/F");
    varia->Branch("bin24", &bin24, "bin24/F");
    varia->Branch("bin25", &bin25, "bin25/F");
    varia->Branch("bin26", &bin26, "bin26/F");
    varia->Branch("bin27", &bin27, "bin27/F");
    varia->Branch("bin28", &bin28, "bin28/F");
    varia->Branch("bin29", &bin29, "bin29/F");
    varia->Branch("bin30", &bin30, "bin30/F");
    varia->Branch("bin31", &bin31, "bin31/F");
    varia->Branch("bin32", &bin32, "bin32/F");
    varia->Branch("bin33", &bin33, "bin33/F");
    varia->Branch("bin34", &bin34, "bin34/F");
    varia->Branch("bin35", &bin35, "bin35/F");
    varia->Branch("bin41", &bin41, "bin41/F");
    varia->Branch("bin42", &bin42, "bin42/F");
    varia->Branch("bin43", &bin43, "bin43/F");
    varia->Branch("bin44", &bin44, "bin44/F");
    varia->Branch("bin45", &bin45, "bin45/F");
    varia->Branch("bin46", &bin46, "bin46/F");
    varia->Branch("bin47", &bin47, "bin47/F");
    varia->Branch("bin48", &bin48, "bin48/F");
    varia->Branch("bin49", &bin49, "bin49/F");
    varia->Branch("bin50", &bin50, "bin50/F");
    varia->Branch("bin51", &bin51, "bin51/F");
    varia->Branch("bin52", &bin52, "bin52/F");
    varia->Branch("bin53", &bin53, "bin53/F");
    varia->Branch("bin54", &bin54, "bin54/F");
    varia->Branch("bin55", &bin55, "bin55/F");
    varia->Branch("bin56", &bin56, "bin56/F");
    varia->Branch("bin57", &bin57, "bin57/F");
    varia->Branch("bin58", &bin58, "bin58/F");
    varia->Branch("bin59", &bin59, "bin59/F");
    varia->Branch("bin60", &bin60, "bin60/F");
    varia->Branch("type",  &type,  "type/I");


    // Weightfiles using Expanded Path Name, hopefully it does work for both windows and linux
    TString weightpath_20=TString(gSystem->ExpandPathName("$FEDRA_ROOT"))+TString("/src/libShower/weights_20.txt");
    TString weightpath_15=TString(gSystem->ExpandPathName("$FEDRA_ROOT"))+TString("/src/libShower/weights_15.txt");
    TString weightpath_30=TString(gSystem->ExpandPathName("$FEDRA_ROOT"))+TString("/src/libShower/weights_30.txt");
    TString weightpath_50=TString(gSystem->ExpandPathName("$FEDRA_ROOT"))+TString("/src/libShower/weights_50.txt");

    TMultiLayerPerceptron *mlp1 = new TMultiLayerPerceptron("@bin3,@bin4,@bin5,@bin6,@bin7,@bin8,@bin9,@bin10,@bin11,@bin12,@bin13,@bin14,@bin15,@bin16,@bin17,@bin18,@bin19,@bin20,@bin31,@bin32,@bin33,@bin34,@bin35:63:21:type",varia,"bin31>0","bin31>0");
    mlp1->LoadWeights(weightpath_20);

    TMultiLayerPerceptron *mlp2 = new TMultiLayerPerceptron("@bin3,@bin4,@bin5,@bin6,@bin7,@bin8,@bin9,@bin10,@bin11,@bin12,@bin13,@bin14,@bin15,@bin31,@bin32,@bin33,@bin34,@bin35:63:21:type",varia,"bin31>0","bin31>0");
    mlp2->LoadWeights(weightpath_15);

    TMultiLayerPerceptron *mlp3 = new TMultiLayerPerceptron("@bin3,@bin4,@bin5,@bin6,@bin7,@bin8,@bin9,@bin10,@bin11,@bin12,@bin13,@bin14,@bin15,@bin16,@bin17,@bin18,@bin19,@bin20,@bin21,@bin22,@bin23,@bin24,@bin25,@bin26,@bin27,@bin28,@bin29,@bin30,@bin31,@bin32,@bin33,@bin34,@bin35:63:21:type",varia,"bin31>0","bin31>0");
    mlp3->LoadWeights(weightpath_30);

    TMultiLayerPerceptron *mlp4 = new TMultiLayerPerceptron("@bin3,@bin4,@bin5,@bin6,@bin7,@bin8,@bin9,@bin10,@bin11,@bin12,@bin13,@bin14,@bin15,@bin16,@bin17,@bin18,@bin19,@bin20,@bin21,@bin22,@bin23,@bin24,@bin25,@bin26,@bin27,@bin28,@bin29,@bin30,@bin41,@bin42,@bin43,@bin44,@bin45,@bin46,@bin47,@bin48,@bin49,@bin50,@bin51,@bin52,@bin53,@bin54,@bin55,@bin56,@bin57,@bin58,@bin59,@bin60,@bin31,@bin32,@bin33,@bin34,@bin35:63:21:type",varia,"bin31>0","bin31>0");
    mlp4->LoadWeights(weightpath_50);

    TH1F *histogl_e = new TH1F("longitudinal profile","longitudinal profile",50,0,50);

    TH1F *historcut30 = new TH1F("Delta r 30"," ",50,0,200);
    TH1F *histotcut30 = new TH1F("Delta theta 30"," ",50,0,0.2);
    TH1F *historcut50 = new TH1F("Delta r 50"," ",50,0,200);
    TH1F *histotcut50 = new TH1F("Delta theta 50"," ",50,0,0.2);

    TH1F *historcut20 = new TH1F("Delta r 20 "," ",50,0,200);
    TH1F *histotcut20 = new TH1F("Delta theta 20"," ",50,0,0.2);
    TH1F *historcut15 = new TH1F("Delta r 15"," ",50,0,200);
    TH1F *histotcut15 = new TH1F("Delta theta 15"," ",50,0,0.2);

    float val30,val20,val15,val50;

    sprintf(fname_e,"shower1.root");  //input shower tree file

    TFile *file_e2 = new TFile(fname_e ,"READ");
    TTree *treebranch_e = (TTree *) file_e2->Get("treebranch");
    treebranch_e->SetBranchAddress("nfilmb", nfilmb);
    treebranch_e->SetBranchAddress("sizeb", &sizeb);
    treebranch_e->SetBranchAddress("sizeb15", &sizeb15);
    treebranch_e->SetBranchAddress("sizeb20", &sizeb20);
    treebranch_e->SetBranchAddress("sizeb30", &sizeb30);
    treebranch_e->SetBranchAddress("isizeb", &isizeb);
    treebranch_e->SetBranchAddress("chi2btkb", &chi2btkb);
    treebranch_e->SetBranchAddress("tagprimary", tagprimary);
    treebranch_e->SetBranchAddress("txb", txb);
    treebranch_e->SetBranchAddress("tyb", tyb);
    treebranch_e->SetBranchAddress("xb", xb);
    treebranch_e->SetBranchAddress("yb", yb);
    treebranch_e->SetBranchAddress("zb", zb);
    treebranch_e->SetBranchAddress("showerID", &showerID);
    treebranch_e->SetBranchAddress("plateb", plateb);
    treebranch_e->SetBranchAddress("idb", idb);
    treebranch_e->SetBranchAddress("deltarb", deltarb);
    treebranch_e->SetBranchAddress("deltathetab", deltathetab);
    treebranch_e->SetBranchAddress("number_eventb", &number_eventb);
    treebranch_e->SetBranchAddress("ntrace1simub", &ntrace1simub);
    treebranch_e->SetBranchAddress("ntrace2simub", &ntrace2simub);
    treebranch_e->SetBranchAddress("ntrace3simub", &ntrace3simub);
    treebranch_e->SetBranchAddress("ntrace4simub", &ntrace4simub);
    treebranch_e->SetBranchAddress("purityb", &purityb);
    treebranch_e->SetBranchAddress("E_MC", &E_MC);
    treebranch_e->SetBranchAddress("mcDigitIndexTop",mcDigitIndexTop);
    treebranch_e->SetBranchAddress("mcDigitIndexBottom",mcDigitIndexBottom);


    for (Int_t ik = 0; ik < treebranch_e->GetEntries(); ik++)
    {
        treebranch_e->GetEntry(ik);
        histogl_e->Reset();
        historcut15->Reset();
        histotcut15->Reset();
        historcut20->Reset();
        histotcut20->Reset();
        historcut30->Reset();
        histotcut30->Reset();
        historcut50->Reset();
        histotcut50->Reset();

        output50 = -10;
        output30 = -10;
        output20 = -10;
        output15 = -10;

        if (isizeb>0) {
            ntrke50 = 0;
            Rcut50 = 0.;
            Tcut50 = 0.;
            ntrke30 = 0;
            Rcut30 = 0.;
            Tcut30 = 0.;
            ntrke20 = 0;
            Rcut20 = 0.;
            Tcut20 = 0.;
            ntrke15 = 0;
            Rcut15 = 0.;
            Tcut15 = 0.;
            flag15 = 0;
            flag20 = 0;
            flag30 = 0;
            flag50 = 0;

            if (isizeb==15)flag15=1;
            if (isizeb==20)flag20=1;
            if (isizeb==30)flag30=1;
            if (isizeb==50)flag50=1;


            for (Int_t ibtke = 0; ibtke < sizeb; ibtke++) { // loop on btk

                histogl_e->Fill(nfilmb[ibtke]);
                if (nfilmb[ibtke]<51)ntrke50++;
                if (nfilmb[ibtke]<31)ntrke30++;
                if (nfilmb[ibtke]<21)ntrke20++;
                if (nfilmb[ibtke]<16)ntrke15++;

                if (ibtke>0&&nfilmb[ibtke]<51)
                {
                    Rcut50 += deltarb[ibtke];
                    Tcut50 += deltathetab[ibtke];
                    historcut50->Fill(deltarb[ibtke]);
                    histotcut50->Fill(deltathetab[ibtke]);
                }
                if (ibtke>0&&nfilmb[ibtke]<31)
                {
                    Rcut30 += deltarb[ibtke];
                    Tcut30 += deltathetab[ibtke];
                    historcut30->Fill(deltarb[ibtke]);
                    histotcut30->Fill(deltathetab[ibtke]);
                }
                if (ibtke>0&&nfilmb[ibtke]<21)
                {
                    Rcut20 += deltarb[ibtke];
                    Tcut20 += deltathetab[ibtke];
                    historcut20->Fill(deltarb[ibtke]);
                    histotcut20->Fill(deltathetab[ibtke]);
                }
                if (ibtke>0&&nfilmb[ibtke]<16)
                {
                    Rcut15 += deltarb[ibtke];
                    Tcut15 += deltathetab[ibtke];
                    historcut15->Fill(deltarb[ibtke]);
                    histotcut15->Fill(deltathetab[ibtke]);
                }
            }

            bin1 = histogl_e->GetBinContent(1);
            bin2 = histogl_e->GetBinContent(2);
            bin3 = histogl_e->GetBinContent(3);
            bin4 = histogl_e->GetBinContent(4);
            bin5 = histogl_e->GetBinContent(5);
            bin6 = histogl_e->GetBinContent(6);
            bin7 = histogl_e->GetBinContent(7);
            bin8 = histogl_e->GetBinContent(8);
            bin9 = histogl_e->GetBinContent(9);
            bin10 = histogl_e->GetBinContent(10);
            bin11 = histogl_e->GetBinContent(11);
            bin12 = histogl_e->GetBinContent(12);
            bin13 = histogl_e->GetBinContent(13);
            bin14 = histogl_e->GetBinContent(14);
            bin15 = histogl_e->GetBinContent(15);
            bin16 = histogl_e->GetBinContent(16);
            bin17 = histogl_e->GetBinContent(17);
            bin18 = histogl_e->GetBinContent(18);
            bin19 = histogl_e->GetBinContent(19);
            bin20 = histogl_e->GetBinContent(20);
            bin21 = histogl_e->GetBinContent(21);
            bin22 = histogl_e->GetBinContent(22);
            bin23 = histogl_e->GetBinContent(23);
            bin24 = histogl_e->GetBinContent(24);
            bin25 = histogl_e->GetBinContent(25);
            bin26 = histogl_e->GetBinContent(26);
            bin27 = histogl_e->GetBinContent(27);
            bin28 = histogl_e->GetBinContent(28);
            bin29 = histogl_e->GetBinContent(29);
            bin30 = histogl_e->GetBinContent(30);
            bin41 = histogl_e->GetBinContent(31);
            bin42 = histogl_e->GetBinContent(32);
            bin43 = histogl_e->GetBinContent(33);
            bin44 = histogl_e->GetBinContent(34);
            bin45 = histogl_e->GetBinContent(35);
            bin46 = histogl_e->GetBinContent(36);
            bin47 = histogl_e->GetBinContent(37);
            bin48 = histogl_e->GetBinContent(38);
            bin49 = histogl_e->GetBinContent(39);
            bin50 = histogl_e->GetBinContent(40);
            bin51 = histogl_e->GetBinContent(41);
            bin52 = histogl_e->GetBinContent(42);
            bin53 = histogl_e->GetBinContent(43);
            bin54 = histogl_e->GetBinContent(44);
            bin55 = histogl_e->GetBinContent(45);
            bin56 = histogl_e->GetBinContent(46);
            bin57 = histogl_e->GetBinContent(47);
            bin58 = histogl_e->GetBinContent(48);
            bin59 = histogl_e->GetBinContent(49);
            bin60 = histogl_e->GetBinContent(50);
// input variables for 20 sheets
            bin31 = ntrke20;
            bin32 = Rcut20/ntrke20;
            bin33 = Tcut20/ntrke20;
            bin34 =  historcut20->GetRMS();
            bin35 =  histotcut20->GetRMS();

            params[0] = bin3;
            params[1] = bin4;
            params[2] = bin5;
            params[3] = bin6;
            params[4] = bin7;
            params[5] = bin8;
            params[6] = bin9;
            params[7] = bin10;
            params[8] = bin11;
            params[9] = bin12;
            params[10] = bin13;
            params[11] = bin14;
            params[12] = bin15;
            params[13] = bin16;
            params[14] = bin17;
            params[15] = bin18;
            params[16] = bin19;
            params[17] = bin20;
            params[18] = bin31;
            params[19] = bin32;
            params[20] = bin33;
            params[21] = bin34;
            params[22] = bin35;

            val20=mlp1->Evaluate(0, params);

// input variables for 30 sheets
            bin31 = ntrke30;
            bin32 = Rcut30/ntrke30;
            bin33 = Tcut30/ntrke30;
            bin34 =  historcut30->GetRMS();
            bin35 =  histotcut30->GetRMS();

            params[18] = bin21;
            params[19] = bin22;
            params[20] = bin23;
            params[21] = bin24;
            params[22] = bin25;
            params[23] = bin26;
            params[24] = bin27;
            params[25] = bin28;
            params[26] = bin29;
            params[27] = bin30;
            params[28] = bin31;
            params[29] = bin32;
            params[30] = bin33;
            params[31] = bin34;
            params[32] = bin35;
            val30 = mlp3->Evaluate(0, params);

// input variables for 50 sheets
            bin31 = ntrke50;
            bin32 = Rcut50/ntrke50;
            bin33 = Tcut50/ntrke50;
            bin34 =  historcut50->GetRMS();
            bin35 =  histotcut50->GetRMS();

            params[28] = bin41;
            params[29] = bin42;
            params[30] = bin43;
            params[31] = bin44;
            params[32] = bin45;
            params[33] = bin46;
            params[34] = bin47;
            params[35] = bin48;
            params[36] = bin49;
            params[37] = bin50;
            params[38] = bin51;
            params[39] = bin52;
            params[40] = bin53;
            params[41] = bin54;
            params[42] = bin55;
            params[43] = bin56;
            params[44] = bin57;
            params[45] = bin58;
            params[46] = bin59;
            params[47] = bin60;

            params[48] = bin31;
            params[49] = bin32;
            params[50] = bin33;
            params[51] = bin34;
            params[52] = bin35;
            val50 = mlp4->Evaluate(0, params);

// input variables for 15 sheets

            bin31 = ntrke15;
            bin32 = Rcut15/ntrke15;
            bin33 = Tcut15/ntrke15;
            bin34 =  historcut15->GetRMS();
            bin35 =  histotcut15->GetRMS();
            params[13] = bin31;
            params[14] = bin32;
            params[15] = bin33;
            params[16] = bin34;
            params[17] = bin35;

            val15=mlp2->Evaluate(0, params);

// save final output values

            output50= val50;
            output30= val30;
            output20= val20;
            output15= val15;
            eProb90=0;
            eProb1=0;

            if (output15>0.38&&isizeb==1) eProb90=1;
            if (flag15==1&&output15>0.38) eProb90=1;
            if (flag20==1&&output20>0.37) eProb90=1;
            if (flag30==1&&output30>0.36) eProb90=1;
            if (flag50==1&&output50>0.34)  eProb90=1;

            if (output15>0.78&&isizeb==1) eProb1=1;
            if (flag15==1&&output15>0.78) eProb1=1;
            if (flag20==1&&output20>0.73) eProb1=1;
            if (flag30==1&&output30>0.7) eProb1=1;
            if (flag50==1&&output50>0.73) eProb1=1;
            /*
            if (output15>0.345) eProb90=1;
            if (output15>0.345&&output20>0.355) eProb90=1;
            if (output15>0.345&&output20>0.355&&output30>0.385) eProb90=1;
            if (output15>0.345&&output20>0.355&&output30>0.385&&output50>0.37) eProb90=1;

            if (output15>0.64) eProb1=1;
            if (output15>0.64&&output20>0.65) eProb1=1;
            if (output15>0.64&&output20>0.65&&output30>0.64) eProb1=1;
            if (output15>0.64&&output20>0.65&&output30>0.64&output50>0.645) eProb1=1;
            */
            if (eProb90==1&&flag50==1)
            {
                cout<<  " output50: " << val50  <<  " output20: " << val20 <<  "   output15: " << val15 << endl;
            }
            treesaveb3->Fill();
            //cout << "Do treesaveb3->Fill();"<<endl;
        }
    }

    eShowersN=treesaveb3->GetEntries();

    file_e2->Close();

    fileout3->cd();
    treesaveb3->Write();
    fileout3->Close();
}


//-------------------------------------------------------------------
int EdbShowerRec::AddInitiatorBT(EdbSegP* seg)
{
    // Check if eInBTArray array exits, otherwise create.
    if (!eInitiatorBTIsCreated) {
        eInBTArray=new TObjArray(9999);
        eInitiatorBTIsCreated=kTRUE;
    }
    // Add the segment to the eInBTArray array
    eInBTArray->Add(seg);

    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::AddInitiatorBT(seg) added. Now entries in eInBTArray: ";
    if (gEDBDEBUGLEVEL>2) cout << eInBTArray->GetEntries() << endl;

    eInBTArrayIsFilled=kTRUE;
    return 1;
}
//-------------------------------------------------------------------




//-------------------------------------------------------------------
void EdbShowerRec::PrintInitiatorBTs()
{
    // Print the eInBTArray nice...
    if (!eInBTArrayIsFilled) {
        cout << "void EdbShowerRec::PrintInitiatorBTs()   WARNING: eInBTArrayIsFilled="<<eInBTArrayIsFilled<<" . RETURN. "<< endl;
        return;
    }

    cout << "void EdbShowerRec::PrintInitiatorBTs()   eInBTArray->GetEntries()=" << eInBTArray->GetEntries() << endl;

    EdbSegP* seg;
    for (int i=0; i< eInBTArray->GetEntries(); ++i) {
        seg=(EdbSegP*) eInBTArray->At(i);
        seg->PrintNice();
    }
    return;
}
//-------------------------------------------------------------------







//-------------------------------------------------------------------
void EdbShowerRec::CalculateEnergyValues()
{
    //=C= -----------------will be deprecated soon-------------

    if (gEDBDEBUGLEVEL>2) cout << "------------------------------------------------------"<<endl;
    if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::CalculateEnergyValues()"<<endl;
    if (gEDBDEBUGLEVEL>2) cout << "------------------------------------------------------"<<endl;

//   if (!gROOT->GetClass("TMultiLayerPerceptron")) {
//     cout << "Need to load libMLP seperately"<<endl;
//     gSystem->Load("libMLP");
//   }

    //gEDBDEBUGLEVEL=3;


    ///--------------------------------------------
    ///---CreateANNTree:
    ///--------------------------------------------
    Energy_CreateANNTree();
    ///--------------------------------------------

    ///--------------------------------------------
    ///---CreateEnergyResolutionFitFunction:
    ///--------------------------------------------
    Energy_CreateEnergyResolutionFitFunction();
    ///--------------------------------------------

    ///--------------------------------------------
    ///---Extract the desired profile
    ///--------------------------------------------
    Energy_ExtractShowerParametrisationProfile();
    ///--------------------------------------------

    return;
}




//-------------------------------------------------------------------
void EdbShowerRec::LoadShowerFile(TString ShowerFileName="Shower2.root")
{
    if (gEDBDEBUGLEVEL>2) cout << "------------------------------------------------------"<<endl;
    if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::LoadShowerFile()"<<endl;
    if (gEDBDEBUGLEVEL>2) cout << "------------------------------------------------------"<<endl;

    FileReconstructedShowerTree = new TFile(ShowerFileName,"READ");
    if (FileReconstructedShowerTree==0)  {
        cout << "ERROR: FileReconstructedShowerTree==NULL. ABORTING."<<endl;
        return;
    }

    //  if (gEDBDEBUGLEVEL>2) FileReconstructedShowerTree->ls();
    if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::LoadShowerFile() done."<<endl;
    return;
}


//-------------------------------------------------------------------
void EdbShowerRec::Energy_ExtractShowerParametrisationProfile()
{
    if (gEDBDEBUGLEVEL>2) cout << "------------------------------------------------------"<<endl;
    if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile()"<<endl;
    if (gEDBDEBUGLEVEL>2) cout << "------------------------------------------------------"<<endl;
    //ouput file
    fileout3 = new TFile("Shower.root","RECREATE");
    treesaveb3 = new TTree("treebranch","tree of branchtrack");
    treesaveb3->Branch("number_eventb",&number_eventb,"number_eventb/I");
    treesaveb3->Branch("sizeb",&sizeb,"sizeb/I");
    treesaveb3->Branch("sizeb15",&sizeb15,"sizeb15/I");
    treesaveb3->Branch("sizeb20",&sizeb20,"sizeb20/I");
    treesaveb3->Branch("sizeb30",&sizeb30,"sizeb30/I");
    treesaveb3->Branch("output15",&output15,"output15/F");
    treesaveb3->Branch("output20",&output20,"output20/F");
    treesaveb3->Branch("output30",&output30,"output30/F");
    treesaveb3->Branch("output50",&output50,"output50/F");
    treesaveb3->Branch("eProb90",&eProb90,"eProb90/I");
    treesaveb3->Branch("eProb1",&eProb1,"eProb1/I");
    treesaveb3->Branch("E_MC",&E_MC,"E_MC/F");
    treesaveb3->Branch("idb",idb,"idb[sizeb]/I");
    treesaveb3->Branch("plateb",plateb,"plateb[sizeb]/I");
    treesaveb3->Branch("showerID",&showerID,"showerID/I");
    treesaveb3->Branch("isizeb",&isizeb,"isizeb/I");
    treesaveb3->Branch("xb",xb,"xb[sizeb]/F");
    treesaveb3->Branch("yb",yb,"yb[sizeb]/F");
    treesaveb3->Branch("zb",zb,"zb[sizeb]/F");
    treesaveb3->Branch("txb",txb,"txb[sizeb]/F");
    treesaveb3->Branch("tyb",tyb,"tyb[sizeb]/F");
    treesaveb3->Branch("nfilmb",nfilmb,"nfilmb[sizeb]/I");
    treesaveb3->Branch("ntrace1simub",ntrace1simub,"ntrace1simu[sizeb]/I");
    treesaveb3->Branch("ntrace2simub",ntrace2simub,"ntrace2simu[sizeb]/I");
    treesaveb3->Branch("ntrace3simub",ntrace3simub,"ntrace3simu[sizeb]/F");
    treesaveb3->Branch("ntrace4simub",ntrace4simub,"ntrace4simu[sizeb]/I");
    treesaveb3->Branch("chi2btkb",chi2btkb,"chi2btkb[sizeb]/F");
    treesaveb3->Branch("deltarb",deltarb,"deltarb[sizeb]/F");
    treesaveb3->Branch("deltathetab",deltathetab,"deltathetab[sizeb]/F");
    treesaveb3->Branch("deltaxb",deltaxb,"deltaxb[sizeb]/F");
    treesaveb3->Branch("deltayb",deltayb,"deltayb[sizeb]/F");
    treesaveb3->Branch("purityb",&purityb,"purityb/F");
    treesaveb3->Branch("trackdensb",&trackdensb,"trackdensb/F");
    treesaveb3->Branch("tagprimary",tagprimary,"tagprimary[sizeb]/F");

    treesaveb3->Branch("mcDigitIndexTop",mcDigitIndexTop,"mcDigitIndexTop[sizeb]/I");
    treesaveb3->Branch("mcDigitIndexBottom",mcDigitIndexBottom,"mcDigitIndexBottom[sizeb]/I");

    treesaveb3->Branch("Energy",&EnergyCorrectedb,"EnergyCorrectedb/F");
    treesaveb3->Branch("EnergyUnCorrected",&EnergyUnCorrectedb,"EnergyUnCorrectedb/F");
    treesaveb3->Branch("EnergySigma",&EnergySigmaCorrectedb,"EnergySigmaCorrectedb/F");
    treesaveb3->Branch("EnergySigmaUnCorrected",&EnergySigmaUnCorrectedb,"EnergySigmaUnCorrectedb/F");

    // -------------------------------------------------------------
    //  If no showerTree was created we have to load the ShowerTree
    //  most probably from an already reconstructed file.
    //  Float_t energy_shot_particle;
    if (!eShowerTreeIsDone) LoadShowerFile();
    // Get Treebranch:
    // Either from shower file (then it has to be casted) or from the
    // reconstruction, then its automatically there.
    if (!eShowerTreeIsDone) {
        treesaveb = (TTree*)(FileReconstructedShowerTree->Get("treebranch"));
        treesaveb->SetBranchAddress("nfilmb", nfilmb);
        treesaveb->SetBranchAddress("sizeb", &sizeb);
        treesaveb->SetBranchAddress("sizeb15", &sizeb15);
        treesaveb->SetBranchAddress("sizeb20", &sizeb20);
        treesaveb->SetBranchAddress("sizeb30", &sizeb30);
        treesaveb->SetBranchAddress("output15", &output15);
        treesaveb->SetBranchAddress("output20", &output20);
        treesaveb->SetBranchAddress("output30", &output30);
        treesaveb->SetBranchAddress("output50", &output50);
        treesaveb->SetBranchAddress("eProb1", &eProb1);
        treesaveb->SetBranchAddress("eProb90", &eProb90);
        treesaveb->SetBranchAddress("isizeb", &isizeb);
        treesaveb->SetBranchAddress("chi2btkb", &chi2btkb);
        treesaveb->SetBranchAddress("tagprimary", tagprimary);
        treesaveb->SetBranchAddress("txb", txb);
        treesaveb->SetBranchAddress("tyb", tyb);
        treesaveb->SetBranchAddress("xb", xb);
        treesaveb->SetBranchAddress("yb", yb);
        treesaveb->SetBranchAddress("zb", zb);
        treesaveb->SetBranchAddress("showerID", &showerID);
        treesaveb->SetBranchAddress("plateb", plateb);
        treesaveb->SetBranchAddress("idb", idb);
        treesaveb->SetBranchAddress("deltarb", deltarb);
        treesaveb->SetBranchAddress("deltathetab", deltathetab);
        treesaveb->SetBranchAddress("number_eventb", &number_eventb);
        treesaveb->SetBranchAddress("ntrace1simub", &ntrace1simub);
        treesaveb->SetBranchAddress("ntrace2simub", &ntrace2simub);
        treesaveb->SetBranchAddress("ntrace3simub", &ntrace3simub);
        treesaveb->SetBranchAddress("ntrace4simub", &ntrace4simub);
        treesaveb->SetBranchAddress("purityb", &purityb);
        treesaveb->SetBranchAddress("trackdensb", &trackdensb);
        treesaveb->SetBranchAddress("E_MC", &E_MC);
        treesaveb->SetBranchAddress("Energy", &EnergyCorrectedb);
        treesaveb->SetBranchAddress("EnergyUnCorrected", &EnergyUnCorrectedb);
        treesaveb->SetBranchAddress("EnergySigma", &EnergySigmaCorrectedb);
        treesaveb->SetBranchAddress("EnergySigmaUnCorrected", &EnergySigmaUnCorrectedb);
    }
    else {
        ; /* treesaveb already declared and filled. */
    }
    // -------------------------------------------------------------


    // Just local variables which are neeeded once for the shower parametrisation...
    Float_t Dr[57];
    Float_t X0[57];
    Float_t Y0[57];
    Float_t TX0=0. ,TY0=0.;
    Float_t theta[57];
    Float_t dist;
    Int_t nentries_e;
    Int_t nentries_withisizeb=0;
    Int_t NumberOfPlates;
    Int_t sizeb_for_plates;
    Float_t BT_deltaR_mean;
    Float_t BT_deltaR_rms;
    Float_t BT_deltaT_mean;
    Float_t BT_deltaT_rms;
    Int_t longprofile[57]; // Depens on how much sheets;
    Int_t numberofilms;
    Float_t EnergyMC;
    TString histname;
    TH1D *histo_nbtk;
    TH1D *histo_longprofile;
    TH1D* histo_deltaR;
    TH1D* histo_deltaT;
    // Just local variables which are neeeded once for the shower parametrisation...END

    //gEDBDEBUGLEVEL=3;

    // Create ShowerProfileTree ...

    if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile(): Create ShowerProfileTree ..."<<endl;
    fileout3->cd();
    ShowerParametrisationTree = new TTree("ShowerProfileTree","ShowerProfileTree");

    ShowerParametrisationTree->Branch("sizeb_for_plates",&sizeb_for_plates,"sizeb_for_plates/I");
    ShowerParametrisationTree->Branch("BT_deltaR_mean",&BT_deltaR_mean,"BT_deltaR_mean/F");
    ShowerParametrisationTree->Branch("BT_deltaR_rms",&BT_deltaR_rms,"BT_deltaR_rms/F");
    ShowerParametrisationTree->Branch("BT_deltaT_mean",&BT_deltaT_mean,"BT_deltaT_mean/F");
    ShowerParametrisationTree->Branch("BT_deltaT_rms",&BT_deltaT_rms,"BT_deltaT_rms/F");
    ShowerParametrisationTree->Branch("longprofile",longprofile,"longprofile[58]/I");
    ShowerParametrisationTree->Branch("numberofilms",&numberofilms,"numberofilms/I");
    ShowerParametrisationTree->Branch("EnergyMC",&E_MC,"EnergyMC/F");

    ShowerParametrisationTree->Branch("Energy",&EnergyCorrectedb,"EnergyCorrectedb/F");
    ShowerParametrisationTree->Branch("EnergyUnCorrected",&EnergyUnCorrectedb,"EnergyUnCorrectedb/F");
    ShowerParametrisationTree->Branch("EnergySigma",&EnergySigmaCorrectedb,"EnergySigmaCorrectedb/F");
    ShowerParametrisationTree->Branch("EnergySigmaUnCorrected",&EnergySigmaUnCorrectedb,"EnergySigmaUnCorrectedb/F");

    //=C= -----------------------------------------------------
    histname="histo_nbtk_av";
    histo_nbtk_av   = new TH1D(histname,"Avergage basetrack numbers",21,0.0,210.0);
    histname="histo_longprofile_av";
    histo_longprofile_av  = new TH1D(histname,"Basetracks per emulsion number",56,0.0,56.0);
    histname="histo_deltaR_mean";
    histo_deltaR_mean     = new TH1D(histname,"Mean #deltar of all BTs in one shower",61,0.0,61.0);
    histname="histo_deltaT_mean";
    histo_deltaT_mean       = new TH1D(histname,"Mean #delta#theta of all BTs in one shower",50,0.0,0.05);
    histname="histo_deltaR_rms";
    histo_deltaR_rms        = new TH1D(histname,"RMS #deltar of all BTs in one shower",61,0.0,61.0);
    histname="histo_deltaT_rms";
    histo_deltaT_rms      = new TH1D(histname,"RMS #delta#theta of all BTs in one shower",50,0.0,0.025);
    //=C= These Histos are for each Shower Entrie filled new
    histname="histo_nbtk";
    histo_nbtk              = new TH1D(histname,"Basetracks in the shower",250,0.0,250.0);
    histname="histo_longprofile";
    histo_longprofile     = new TH1D(histname,"Basetracks per emulsion number",56,0.0,56.0);
    histname="histo_deltaR";
    histo_deltaR            = new TH1D(histname,"Mean #deltar of all BTs in Shower",100,0.0,100.0);
    histname="histo_deltaT";
    histo_deltaT          = new TH1D(histname,"Mean #delta#theta of all BTs in Shower",100,0.0,0.05);
    //=C= -----------------------------------------------------
    eHistEnergyNominal      = new TH1F("eHistEnergyNominal","eHistEnergyNominal",100,0,20);
    eHistEnergyCorrectedb   = new TH1F("eHistEnergyCorrectedb","eHistEnergyCorrectedb",100,0,20);
    eHistEnergyUnCorrectedb = new TH1F("eHistEnergyUnCorrectedb","eHistEnergyUnCorrectedb",100,0,20);
    //=C= -----------------------------------------------------
    histo_nbtk_av->Reset();
    histo_nbtk->Reset();
    histo_longprofile_av->Reset();
    histo_deltaR_mean->Reset();
    histo_deltaT_mean->Reset();
    histo_deltaR_rms->Reset();
    histo_deltaT_rms->Reset();
    //=C= -----------------------------------------------------

    //=C= Again Define some Variables used for calculating:
    nentries_e = (Int_t) treesaveb -> GetEntries();
    for (Int_t j = 0; j<57; ++j) {
        Dr[j] = 0.03*j*1273. +20.0;
    }
    Float_t zmax=0;
    Float_t zmin=0;

    //=C= -----------------------------------------------------
    Double_t params[70];
    Double_t val;
    Double_t val_corr;
    Double_t sigma_val_corr;
    Double_t sigma_val;

    Int_t NrOfANNInputNeurons[]= {11,16,21,26,31,36,41,46};
    for (Int_t i=0; i<8; ++i) NrOfANNInputNeurons[i]=NrOfANNInputNeurons[i]+1+4;
    //=C= -----------------------------------------------------

    cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile(): Loop over reconstructed showers (=" << nentries_e << " total)."<< endl;
    //=C= -----------------------------------------------------
    //=C= Loop over reconstructed Showers (Entries):
    for (Int_t ige = 0; ige < nentries_e; ++ige) {

        /// DEBUG -----------------------------------
        //    if (ige>5&&ige<nentries_e-5) continue;
        //    if (ige>900) continue;
        //    if (ige<nentries_e-900) continue;
        /// DEBUG -----------------------------------

        if (ige%50==0) cout<<"===   Doing Entry "<< ige << endl;
        treesaveb->GetEntry(ige);
        if (gEDBDEBUGLEVEL>3) treesaveb->Show(ige);  /// DEBUG OPTION

        // GetEnergyMC:
        EnergyMC=E_MC;
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() === EnergyMC= "<< EnergyMC << endl;

        // Search for local zmax and zmin in the whole shower event to
        // guess from that the number of plates...
        zmax=0;
        zmin=0; // Reinit zmax, min;
        for (Int_t j=0; j<sizeb; ++j) {
            if (zb[j]>zmax) zmax=zb[j];
            if (zb[j]<zmin) zmin=zb[j];
        }
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() === zmax= "<< zmax << endl;
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() === zmin= "<< zmin << endl;

        // Calculate the plates.
        Int_t pl=int(zmax-zmin)/1273+1;
        NumberOfPlates=pl;
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() === pl= "<< pl << endl;

        // Make sure that its at least 1 or more...
        NumberOfPlates=max(NumberOfPlates,1);
        NumberOfPlates=min(NumberOfPlates,57);
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  NumberOfPlates= "<< NumberOfPlates << endl;


        histo_longprofile      ->Reset();
        histo_deltaR           ->Reset();
        histo_deltaT           ->Reset();
        for (int id=0; id<57; id++) {
            theta[id]= 0;
            X0[id] = id*1273.*txb[0] +  xb[0];
            Y0[id] = id*1273.*tyb[0] +  yb[0];
            longprofile[id]=-1;
        }
        if (tagprimary[0]==1) {
            TX0 = txb[0];
            TY0 = tyb[0];
        }
        ++nentries_withisizeb;

        sizeb_for_plates=0;
        numberofilms=NumberOfPlates;
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  numberofilms= " << numberofilms <<endl;

        //=C= -----------------------------------------------------
        //=C= loop over the basetracks in the shower (boucle sur les btk)
        for (Int_t ibtke = 0; ibtke < sizeb; ibtke++) {
            dist = sqrt((xb[ibtke]- X0[nfilmb[ibtke]-1])*(xb[ibtke]- X0[nfilmb[ibtke]-1])+(yb[ibtke]- Y0[nfilmb[ibtke]-1])*(yb[ibtke]- Y0[nfilmb[ibtke]-1]));
            // inside the cone
            if (dist<Dr[nfilmb[ibtke]-1]&&nfilmb[ibtke]<=numberofilms) {

                histo_longprofile         ->Fill(nfilmb[ibtke]);
                histo_longprofile_av      ->Fill(nfilmb[ibtke]);
                theta[nfilmb[ibtke]]+= (TX0-txb[ibtke])*(TX0-txb[ibtke])+(TY0-tyb[ibtke])*(TY0-tyb[ibtke]);

                if (nfilmb[ibtke]<=numberofilms) {
                    ++sizeb_for_plates;
                }
                if (ibtke>0&&nfilmb[ibtke]<=numberofilms) {
                    histo_deltaR            ->Fill(deltarb[ibtke]);
                    histo_deltaT            ->Fill(deltathetab[ibtke]);
                }
            }
        }//=C= END OF loop over the basetracks in the shower
        //=C= -----------------------------------------------------

        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  ige sizeb_for_plates sizeb" << ige << "  " << sizeb_for_plates <<"   " <<  sizeb << endl;
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  histo_deltaR->GetMean() histo_deltaT->GetMean()" << histo_deltaR->GetMean() << "  " << histo_deltaT->GetMean() <<"   " << endl;
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  histo_deltaR->GetRMS() histo_deltaT->GetRMS()" << histo_deltaR->GetRMS() << "  " << histo_deltaT->GetRMS() <<"   " << endl;

        //=C= Fill NumberBT Hist----Just temporary....for all showers:
        histo_nbtk                    ->Fill(sizeb_for_plates);
        histo_nbtk_av                 ->Fill(sizeb_for_plates);
        //=C= Fill dR,dT Mean and RMS Histos for all showers:
        histo_deltaR_mean             ->Fill(histo_deltaR->GetMean());
        histo_deltaT_mean             ->Fill(histo_deltaT->GetMean());
        histo_deltaR_rms              ->Fill(histo_deltaR->GetRMS());
        histo_deltaT_rms              ->Fill(histo_deltaT->GetRMS());

        //=C= -----------------------------------------------------
        // "Variable" 0 describes the (MonteCarlo) energy of the particle
        // energy_shot_particle;  (if there)
        // Variable 1 corresponds to number of Basetracks in the Shower.
        // sizeb;
        // Variable 2,3,4,5 corresponds to dR, mean,rms; dT, mean,rms;
        BT_deltaR_mean=histo_deltaR->GetMean();
        BT_deltaR_rms=histo_deltaR->GetRMS();
        BT_deltaT_mean=histo_deltaT->GetMean();
        BT_deltaT_rms=histo_deltaT->GetRMS();
        // Variables longprofile[...] corresponds to BinEntries Of the longitudialProfile (averaged):
        // Take Care!! BinContent(0) is the UnderflowBin! so it has to go from 0 (plate1) to nbfilm+1 (f.e.plate11->nbfilm+1=12) !
        for (Int_t i=1; i<=numberofilms+1; ++i) {
            longprofile[i-1] = int(histo_longprofile->GetBinContent(i));
        }
        //=C= -----------------------------------------------------


        //=C= -----------------------------------------------------
        //=C= --- Here we differentiate which Platetype we have.
        //=C= --- Can be for each Shower different!
        Int_t type=0;
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  NumberOfPlates    " <<  NumberOfPlates << endl;
        if (NumberOfPlates<16) {
            type=0;
        }
        if (NumberOfPlates>=16 && NumberOfPlates<21) {
            type=1;
        }
        if (NumberOfPlates>=21 && NumberOfPlates<26) {
            type=2;
        }
        if (NumberOfPlates>=26 && NumberOfPlates<31) {
            type=3;
        }
        if (NumberOfPlates>=31 && NumberOfPlates<36) {
            type=4;
        }
        if (NumberOfPlates>=36 && NumberOfPlates<41) {
            type=5;
        }
        if (NumberOfPlates>=41 && NumberOfPlates<46) {
            type=6;
        }
        if (NumberOfPlates>=46) {
            type=7;
        }
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  Type   " << type << endl;
        //=C= -----------------------------------------------------


        //=C= -----------------------------------------------------
        //=C= Fill the ANN Inputvariables
        //=C= For the specific Sample
        //=C= -----------------------------------------------------
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  " <<  ige;
        if (gEDBDEBUGLEVEL>2) cout << "     type=  " << type << endl;
        // Reset All  inANN  Values to Zero:
        for (Int_t j=0; j<70; j++) {
            inANN[j]=0;
        }

        // Variable 0 corresponds to the "type" i.e. to the Energy of the Shower.
        inANN[0]= EnergyMC;
        // Variable 1 corresponds to number of Basetracks in the Shower.
        inANN[1]= sizeb_for_plates;
        // Variable 2,3,4,5 corresponds to dR, mean,rms; dT, mean,rms;
        inANN[2]= BT_deltaR_mean;
        inANN[3]= BT_deltaR_rms;
        inANN[4]= BT_deltaT_mean;
        inANN[5]= BT_deltaT_rms;
        // Variable 6... (end) to the longitudinale profile...
        for (Int_t j=0; j<numberofilms; j++) inANN[6+j]=longprofile[j+1];
        //=C= -----------------------------------------------------


        //---------------Parameters as   param array suited for ANN to evaluate: ---------------
        // Reset All Values to Zero:
        for (Int_t j=0; j<70; j++) {
            params[j]=0;
        }
        //=C= NeuralNetworkLayout: check to which layout this specific shower belongs:
        for (Int_t j=1; j<=NrOfANNInputNeurons[type]; ++j) {
            params[j-1]=inANN[j];
            if (gEDBDEBUGLEVEL>3) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  j  params[j-1]=inANN[j]   " << j<< "  " << params[j-1] << "  " << inANN[j] << endl; ///DEBUG OPTION
        }


        //--- Main Core Part: Evaluate function of Neural Network, specific to the
        //---                 ANN which is usefull for the shower length!
        //------------------------------------------
        TMultiLayerPerceptron* mlp= (TMultiLayerPerceptron*)eEnergyANN->At(type);
        val=mlp->Evaluate(0,params);
        val_corr=Energy_CorrectEnergy(val, type);
        //------------------------------------------
        sigma_val=Energy_CalcSigma(val, type);
        sigma_val_corr=Energy_CalcSigma(val_corr, type);
        //------------------------------------------

        //------------------------------------------
        EnergyCorrectedb=val_corr;
        EnergyUnCorrectedb=val;
        EnergySigmaCorrectedb=sigma_val_corr;
        EnergySigmaUnCorrectedb=sigma_val;
        EnergyMC=E_MC;
        //------------------------------------------

        if (gEDBDEBUGLEVEL>3) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() === " << mlp->GetStructure().Data() <<endl; ///DEBUG OPTION
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  EnergyUnCorrectedb, EnergyCorrectedb,EnergyMC:    " <<  EnergyUnCorrectedb << "    ==  " << EnergyCorrectedb << "      " <<  EnergyMC << endl; ///DEBUG OPTION
        if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_ExtractShowerParametrisationProfile() ===  EnergySigmaCorrectedb, EnergySigmaUnCorrectedb:    " <<  EnergySigmaUnCorrectedb << "    ==  " << EnergySigmaCorrectedb << "      " << endl; ///DEBUG OPTION

        eHistEnergyUnCorrectedb->Fill(val);
        eHistEnergyCorrectedb->Fill(val_corr);
        eHistEnergyNominal->Fill( inANN[0]/1000. );

        ShowerParametrisationTree->Fill();
        treesaveb3->Fill();
        //--------------------------------------------------------------------
    }  //=C= END of Loop over reconstructed Showers (Entries):
    cout << "Loop over reconstructed Showers finished."<< endl;
    
    // Scale histo_longprofile_av to number of entries:
    histo_longprofile_av->Scale(1.0/nentries_withisizeb);

    // Write to File:
    fileout3->cd();
    treesaveb3->Write();
    fileout3->Close();

    eEnergyIsDone=kTRUE;
    return;
}




//-------------------------------------------------------------------
Float_t EdbShowerRec::Energy_CorrectEnergy(Float_t MeasuredEnergy, Int_t NPlatesNr)
{
    //=C= -----------------will be deprecated soon-------------
    //=C= Correct the given ANN  Energy linear according
    //=C= to the formula
    //=C= E_corr = 1/p1 * ( E_meas -p0 )
    //=C= p1,p0 are different for the plate numbers.
    //=C= -----------------------------------------------------
    Float_t p0[8];
    Float_t p1[8];
    Float_t CorrEnergy;
    // New Values
    p0[0]=1.01;
    p1[0]=2.83;  // 11 Plates
    p0[1]=1.07;
    p1[1]=1.48;  // 16 Plates
    p0[2]=0.97;
    p1[2]=1.03;  // 21 Plates
    p0[3]=0.98;
    p1[3]=0.89;  // 26 Plates
    p0[4]=0.92;
    p1[4]=1.05;  // 31 Plates
    p0[5]=0.90;
    p1[5]=1.03;  // 36 Plates
    p0[6]=0.88;
    p1[6]=1.07;  // 41 Plates
    p0[7]=0.87;
    p1[7]=1.10;  // 46 Plates
    if (gEDBDEBUGLEVEL>2) cout << "void EdbShowerRec::Energy_CorrectEnergy()                   === NPlatesNr  p0[NPlatesNr] " << NPlatesNr << "  " << p0[NPlatesNr] << endl;
    CorrEnergy= 1.0/p1[NPlatesNr] * ( MeasuredEnergy-p0[NPlatesNr] );

    return CorrEnergy;
}





//-------------------------------------------------------------------
void EdbShowerRec::Energy_CreateEnergyResolutionFitFunction()
{
    //=C= -----------------will be deprecated soon-------------
    eEnergyResolutionFitFunction_All = new TF1("EnergyResolutionFitFunction_All","TMath::Sqrt([0]*[0]/sqrt(x)/sqrt(x)+[2]*[2]+[1]*[1]/x/x)",0.9,48.0);
    eEnergyResolutionFitFunction_All->SetParameter(0,0.5);  //Initialize with A=50%
    eEnergyResolutionFitFunction_All->SetParameter(1,0.2);  //Initialize with B=20%
    eEnergyResolutionFitFunction_All->SetParameter(2,0.1);  //Initialize with C=10%
    eEnergyResolutionFitFunction_All->SetParName(0,"ResolutionFactor A");
    eEnergyResolutionFitFunction_All->SetParName(1,"ResolutionFactor B");
    eEnergyResolutionFitFunction_All->SetParName(2,"ResolutionFactor C");
    eEnergyResolutionFitFunction_All->SetLineColor(2);
    eEnergyResolutionFitFunction_All->SetParLimits(0,0,2);
    eEnergyResolutionFitFunction_All->SetParLimits(1,0,2);
    eEnergyResolutionFitFunction_All->SetParLimits(2,0,2);

    return;
}


//-------------------------------------------------------------------
Float_t EdbShowerRec::Energy_CalcSigma(Float_t Energy, Int_t type)
{
    //=C= -----------------will be deprecated soon-------------
    // To understand where these values come from, please look in
    // the note of energy measurement by FWM. to be found on the opera doc server.
    //=C= -----------------------------------------------------
    Float_t SigmaParametrizationArray[8][3]= { { 1.27,0.0,0.28}, {0.79,0.0,0.26}, {0.48,0.0,0.23}, {0.50,0.0,0.15} ,{0.49,0.0,0.10} ,{0.50,0.17,0.04}, {0.49,0.24,0.00},{0.45,0.37,0.03} };

    eEnergyResolutionFitFunction_All->SetParameter(0,SigmaParametrizationArray[type][0]);
    eEnergyResolutionFitFunction_All->SetParameter(1,SigmaParametrizationArray[type][1]);
    eEnergyResolutionFitFunction_All->SetParameter(2,SigmaParametrizationArray[type][2]);

    Float_t CalcSigma = eEnergyResolutionFitFunction_All->Eval(Energy);
    return CalcSigma;
}




//-------------------------------------------------------------------
void EdbShowerRec::Energy_CreateANNTree()
{
    //=C= -----------------will be deprecated soon-------------

    if (gEDBDEBUGLEVEL>2) cout  << "------------------------------------------"<<endl;
    if (gEDBDEBUGLEVEL>2) cout  << "void EdbShowerRec::Energy_CreateANNTree()"<<endl;
    if (gEDBDEBUGLEVEL>2) cout  << "------------------------------------------"<<endl;

    //=C= Load MLP Root Class:
//  if (!gROOT->GetClass("TMultiLayerPerceptron")) {
//    gSystem->Load("libMLP");
//  }
    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Energy_CreateANNTree() === Load(libMLP) done." << endl;

    //=C= -----------------------------------------------------
    //=C=
    //=C= -----------------------------------------------------
    //=C= Explanation: ReconstructedPlatesArray[ReconstructedPlatesNr]:
    //=C= for BTs per each EmulsionNr (11,16,21..56)
    //=C= Explanation: 4: dR, dT: mean and rms.
    //=C= Explanation: 1: number of Basetracks
    //=C= Explanation: [0]: type=energy of the Shower
    //=C= For the FJ Profile we have 1+4+(11,16,21=NumberOfPlates)... InputNeurons:
    //=C= -----------------------------------------------------
    //=C= We Create different instances of ANNs for different shower lengths:
    //=C= -----------------------------------------------------
/// DEBUG //  inANN = new Double_t[70];
    ANNInterimTreeInputvariables = new TTree("ANNInterimTreeInputvariables","ANNInterimTreeInputvariables");
    ANNInterimTreeInputvariables  ->  Branch("inANN", inANN, "inANN[70]/D");
    //=C= -----------------------------------------------------

    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Energy_CreateANNTree() === CreateANNTree()  is done."<<endl;
    //gEDBDEBUGLEVEL=3;

    Int_t ReconstructedPlatesArray[]  = {11,16,21,26,31,36,41,46};
    Int_t NrOfANNInputNeurons[]       = {11,16,21,26,31,36,41,46};
    eEnergyANN = new TObjArray(9);

    //=C= -----------------------------------------------------
    //=C= This ANN Layout has the following form:
    //=C= inputNeurons            :     hiddenLayer1          :   hiddenLayer2          :   type
    //=C= inANN[1],...,inANN[x]   :     NrOfANNInputNeurons+1 :   NrOfANNInputNeurons   :   inANN[0]
    //=C= -----------------------------------------------------

    for (Int_t i=0; i<8; ++i) {
        NrOfANNInputNeurons[i]=NrOfANNInputNeurons[i]+1+4;
        TString layout="";
        //=C= -----------------------------------------------------
        for (Int_t j=1; j<NrOfANNInputNeurons[i]; ++j) layout += "@inANN["+TString(Form("%d",j))+"],";
        layout += "@inANN["+TString(Form("%d",NrOfANNInputNeurons[i]))+"]:"+TString(Form("%d",NrOfANNInputNeurons[i]+1))+":"+TString(Form("%d",NrOfANNInputNeurons[i]));
        layout+=":inANN[0]/1000"; //if (Geant4) ([E]=MeV) //else layout+=":inANN[0]"; (Geant3) ([E]=GeV)
        //=C= -----------------------------------------------------
//     eEnergyANN[i]
        TMultiLayerPerceptron* eEANN =  new TMultiLayerPerceptron(layout,ANNInterimTreeInputvariables,"","(Entry$)%2");

        //=C= -----------------------------------------------------
        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Energy_CreateANNTree() === ANN_MLP->GetStructure() : " << eEANN->GetStructure().Data() << endl ;
        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Energy_CreateANNTree() === LoadANNWeightFile():"<<endl;

        TString weightfilestring=TString(gSystem->ExpandPathName("$FEDRA_ROOT"))+TString("/src/libShower/")+TString(Form("EnergyMeasurement_WeightFile_%d_Plates.txt",ReconstructedPlatesArray[i]));

        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Energy_CreateANNTree() === weightfilestring=  " << weightfilestring.Data() << endl;


        eEANN->LoadWeights(weightfilestring);
        eEnergyANN->Add(eEANN);
    }
    //=C= -----------------------------------------------------

    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Energy_CreateANNTree() === is done."<<endl;
    return;
}

//-------------------------------------------------------------------
int EdbShowerRec::SaveLNK()
{
    const char *fname="linked_tracks.root";
    const char *rcut="";

    FILE *f = fopen(fname, "r");
    if (!f)
    {
        printf("Cannot find linked_tracks.root file !\n");
        return -1;
    }
    fclose(f);

    EdbPVRec *ali = new EdbPVRec();
    EdbDataProc *dproc = new EdbDataProc();
    dproc->ReadTracksTree(*ali, "linked_tracks.root",rcut);

    sprintf(cmd,"cp -v linked_tracks.root linked_tracks.backup.root");
    gSystem->Exec(cmd);

    int ntracks=ali->Ntracks();
    cout << "nb tracks: " << ntracks << endl;

    ali->Link();
    ali->FillTracksCell();
    ali->MakeTracks(2,0);

    // --------------------------------------------------------------

    EdbSegP *segS;
    EdbTrackP *trackS;

    char fname_e[128];
    float x,y,tx,ty,z;

    Int_t sizeb;
    Int_t isizeb;
    Float_t xb[10000];
    Float_t yb[10000];
    Float_t zb[10000];
    Float_t txb[10000];
    Float_t tyb[10000];
    Int_t nfilmb[10000];
    //   Int_t ngrainb[10000];
    //   Int_t ntrace1simub[10000];
    Int_t ntrace2simub[10000];
    Float_t chi2btkb[10000];
    int eProb90;
    int eProb1;

    int number_eventb,w2,plate[10000];

    sprintf(fname_e,"Shower.root");  //input shower tree file

    f = fopen(fname_e, "r");
    if (!f)
    {
        printf("Cannot find Shower.root file !\n");
        return -2;
    }
    fclose(f);


    TFile *file_e2 = new TFile(fname_e ,"READ");
    TTree *treebranch_e = (TTree*)file_e2->Get("treebranch");
    treebranch_e->SetBranchAddress("nfilmb", nfilmb);   // max(nfilmb)=npl
    treebranch_e->SetBranchAddress("sizeb", &sizeb);   // nseg
    treebranch_e->SetBranchAddress("isizeb", &isizeb);
    treebranch_e->SetBranchAddress("chi2btkb", &chi2btkb);
    treebranch_e->SetBranchAddress("txb", txb);
    treebranch_e->SetBranchAddress("tyb", tyb);
    treebranch_e->SetBranchAddress("xb", xb);
    treebranch_e->SetBranchAddress("yb", yb);
    treebranch_e->SetBranchAddress("zb", zb);
    treebranch_e->SetBranchAddress("number_eventb", &number_eventb);
    treebranch_e->SetBranchAddress("ntrace2simub", &ntrace2simub); // grain
    treebranch_e->SetBranchAddress("eProb1", &eProb1);
    treebranch_e->SetBranchAddress("eProb90", &eProb90);
    int nentries =  treebranch_e->GetEntries();

    for (int i=0; i<nentries; i++)
    {
        trackS = new EdbTrackP();
        treebranch_e->GetEntry(i);

        if (eProb1==0||eProb90==0) continue; // remove non electron!!!

        for (int j=0; j<sizeb; j++)
        {
            x=xb[j];
            y=yb[j];
            tx=txb[j];
            ty=tyb[j];
            z=zb[j];
            // plate[j] = abs(57+(int)(z/1273.));  // compiler instability with gcc4.3
            plate[j] = TMath::Abs(57+(int)(z/1273.));  // compiler instability with gcc4.3
            // cout << plate[j] << endl;
            w2=ntrace2simub[j];
            segS = new EdbSegP();
            segS->Set(j, x, y, tx, ty, w2, 0);
            segS->SetZ(z);
            segS->SetDZ(300.);
            //segS->SetPID(plate[j]); // DISPLAY PROBLEM IF RIGHT PID IS USED !!!!!!!!!!!!!!
            segS->SetPID(j);
            trackS->AddSegment(segS);
            trackS->AddSegmentF(segS);
//    cout << i<< "     " << x << " " << y << " " << z << endl;
        }
        trackS->SetSegmentsTrack(i);
        trackS->SetPID(plate[0]);
        trackS->SetID(i);
        trackS->SetNpl(nfilmb[sizeb-1]);
        trackS->SetProb(1.);
        trackS->SetFlag(-100);
        trackS->FitTrack();
        ali->AddTrack(trackS);
    }

    file_e2->Close();

    // --------------------------------------------------------------

    dproc->MakeTracksTree(ali,"linked_tracks.root");

    cout << endl << "----------------------------------------------" << endl;
    cout << "Nb of tracks before adding showers :  " << ntracks << endl;
    cout << "Nb of tracks after adding showers  :  " << ali->Ntracks() << endl;
    cout << "Nb of stored showers               :  " << ali->Ntracks()-ntracks << endl;
    cout << "----------------------------------------------" << endl << endl;

    return 1;
}





//-------------------------------------------------------------------
TObjArray* EdbShowerRec::ShowerToEdbSegPArray()
{
    // Very simple function to transfer all treebranch entries to an
    // TObjArray of EdbSegPs. Mainly used for eda.C drawing framework

    cout << "-------------------EdbShowerRec::ShowerToEdbSegPArray()---------------------------" << endl << endl;
    TObjArray* segArray=new TObjArray(999999);

    float x,y,tx,ty,z;
    Int_t sizeb;
    Float_t xb[10000];
    Float_t yb[10000];
    Float_t zb[10000];
    Float_t txb[10000];
    Float_t tyb[10000];
    Int_t nfilmb[10000];
    Int_t ntrace2simub[10000];
    int eProb90;
    int eProb1;
    int w2,plate[10000];

    if (!eShowerTreeIsDone) {
        TTree * treesaveb = (TTree*)(FileReconstructedShowerTree->Get("treebranch"));
        treesaveb->SetBranchAddress("sizeb", &sizeb);
        treesaveb->SetBranchAddress("txb", txb);
        treesaveb->SetBranchAddress("tyb", tyb);
        treesaveb->SetBranchAddress("xb", xb);
        treesaveb->SetBranchAddress("yb", yb);
        treesaveb->SetBranchAddress("zb", zb);
        treesaveb->SetBranchAddress("eProb90", &eProb90);
        treesaveb->SetBranchAddress("eProb1", &eProb1);
        treesaveb->SetBranchAddress("nfilmb", nfilmb);

    }
    int nentries =  treesaveb->GetEntries();
    if (gEDBDEBUGLEVEL>2) treesaveb->Print();

    cout << "nentries= " << nentries << endl;

    for (int i=0; i<nentries; i++)
    {
        treesaveb->GetEntry(i);
        treesaveb->Show(i);

        for (int j=0; j<sizeb; j++)
        {
            x=xb[j];
            y=yb[j];
            tx=txb[j];
            ty=tyb[j];
            z=zb[j];
            // This is to have the absolute plate reference frame where Plate 57 is at Z==0
            // plate[j] = abs(57+(int)(z/1273.)); // compiler instability with gcc4.3
            plate[j] = TMath::Abs(57+(int)(z/1273.)); // fix.
            cout << "DOING j= " << j << "    Z==  " << z << "    Platee[j]=   " << plate[j] << endl;
            w2=ntrace2simub[j];
            EdbSegP* segS = new EdbSegP();
            segS->Set(j, x, y, tx, ty, w2, 0);
            segS->SetZ(z);
            segS->SetDZ(300.);
            //segS->SetPID(plate[j]); // DISPLAY PROBLEM IF RIGHT PID IS USED !!!!!!!!!!!!!!
            segS->SetPID(j);
            segArray->Add(segS);
//    cout << i<< "     " << x << " " << y << " " << z << endl;
        }
    }
    cout << "segArray->GetEntries();= " << segArray->GetEntries() << endl;

    return segArray;
}





//----------------------------------------------------------------



Int_t EdbShowerRec::GetTreeBranchEntryNr(EdbSegP* seg)
{
    // Returns the number of the treebranch entry for which the segment is starting point of a shower!
    // if the segment is not starting point at all, return -1;

    cout << "EdbShowerRec::GetTreeBranchEntryNr()"<< endl;

    if (!eShowerTreeIsDone) {
        LoadShowerFile("Shower.root");
        treesaveb = (TTree*)(FileReconstructedShowerTree->Get("treebranch"));
        treesaveb->SetBranchAddress("txb", txb);
        treesaveb->SetBranchAddress("tyb", tyb);
        treesaveb->SetBranchAddress("xb", xb);
        treesaveb->SetBranchAddress("yb", yb);
        treesaveb->SetBranchAddress("zb", zb);
    }
    else {
        cout << "EdbShowerRec::GetTreeBranchEntryNr    treesaveb already declared and filled"<<endl;
        /* treesaveb already declared and filled. */
    }

    cout << "EdbShowerRec::GetTreeBranchEntryNr    Looping over treesaveb now..."<<endl;

    for (int i=0; i<treesaveb->GetEntries(); i++) {
        treesaveb->GetEntry(i);
        if (TMath::Abs(seg->X()-xb[0])>1.0) continue;
        if (TMath::Abs(seg->Y()-yb[0])>1.0) continue;
        if (TMath::Abs(seg->Z()-zb[0])>1.0) continue;
        if (TMath::Abs(seg->TX()-txb[0])>0.01) continue;
        if (TMath::Abs(seg->TY()-tyb[0])>0.01) continue;
        treesaveb->Show(i);
        return i;
    }
    return -1;
}



//----------------------------------------------------------------



Float_t* EdbShowerRec::GetEnergyValues(Int_t TreeEntry)
{
    Float_t* dummy = new Float_t[5];
    dummy[0]=-999.;
    dummy[1]=-999.;
    dummy[2]=-999.;
    dummy[3]=-999.;
    dummy[4]=-999.;

    cout << "eShowerTreeIsDone = " << eShowerTreeIsDone << endl;

    if (!eShowerTreeIsDone) {
        cout << "EdbShowerRec::GetEnergyValues    !eShowerTreeIsDone" << endl;
        LoadShowerFile("Shower.root");
    }
    // Get Treebranch:
    // Either from shower file (then it has to be casted) or from the
    // reconstruction, then its automatically there.
    if (!eShowerTreeIsDone) {
        LoadShowerFile("Shower.root");
        treesaveb = (TTree*)(FileReconstructedShowerTree->Get("treebranch"));
        treesaveb->SetBranchAddress("sizeb", &sizeb);
        treesaveb->SetBranchAddress("E_MC", &E_MC);
        treesaveb->SetBranchAddress("Energy", &EnergyCorrectedb);
        treesaveb->SetBranchAddress("EnergyUnCorrected", &EnergyUnCorrectedb);
        treesaveb->SetBranchAddress("EnergySigma", &EnergySigmaCorrectedb);
        treesaveb->SetBranchAddress("EnergySigmaUnCorrected", &EnergySigmaUnCorrectedb);
    }
    else {
        cout << "EdbShowerRec::GetEnergyValues    treesaveb already declared and filled"<<endl;
        /* treesaveb already declared and filled. */
    }
    cout << "EdbShowerRec::GetEnergyValues    LoadShowerFile  SetBranchAddress  done" << endl;

    treesaveb->GetEntry(TreeEntry);
    if (gEDBDEBUGLEVEL>2) treesaveb->Show(TreeEntry);

    dummy[0]=E_MC;
    dummy[1]=EnergyCorrectedb;
    dummy[2]=EnergyUnCorrectedb;
    dummy[3]=EnergySigmaCorrectedb;
    dummy[4]=EnergySigmaUnCorrectedb;

    cout << "E_MC  EnergyCorrectedb  EnergyUnCorrectedb  EnergySigmaCorrectedb  EnergySigmaUnCorrectedb:" << endl;
    cout << "  " <<  dummy[0] << "  " << dummy[1] << "  " << dummy[2] << "  " << dummy[3] <<  "  " << dummy[4] << "  " << endl;

    return dummy;
}



//----------------------------------------------------------------



Float_t* EdbShowerRec::GetShowerValues(Int_t TreeEntry)
{
    Float_t* dummy = new Float_t[10];
    dummy[0]=-999.;
    dummy[1]=-999.;
    dummy[2]=-999.;
    dummy[3]=-999.;
    dummy[4]=-999.;
    dummy[5]=-999.;
    dummy[6]=-999.;
    dummy[7]=-999.;

    if (!eShowerTreeIsDone) {
        cout << "EdbShowerRec::GetEnergyValues    !eShowerTreeIsDone" << endl;
        LoadShowerFile("Shower.root");
    }
    // Get Treebranch:
    // Either from shower file (then it has to be casted) or from the
    // reconstruction, then its automatically there.
    if (!eShowerTreeIsDone) {
        LoadShowerFile("Shower.root");
        treesaveb = (TTree*)(FileReconstructedShowerTree->Get("treebranch"));
        treesaveb->SetBranchAddress("sizeb", &sizeb);
        treesaveb->SetBranchAddress("txb", txb);
        treesaveb->SetBranchAddress("tyb", tyb);
        treesaveb->SetBranchAddress("xb", xb);
        treesaveb->SetBranchAddress("yb", yb);
        treesaveb->SetBranchAddress("zb", zb);
        treesaveb->SetBranchAddress("eProb90", &eProb90);
        treesaveb->SetBranchAddress("eProb1", &eProb1);
        treesaveb->SetBranchAddress("nfilmb", nfilmb);

    }
    else {
        ; /* treesaveb already declared and filled. */
    }
    cout << "EdbShowerRec::GetEnergyValues    LoadShowerFile  SetBranchAddress  done" << endl;

    treesaveb->GetEntry(TreeEntry);
    if (gEDBDEBUGLEVEL>2) treesaveb->Show(TreeEntry);

    dummy[0]=sizeb;
    dummy[1]=xb[0];
    dummy[2]=yb[0];
    dummy[3]=zb[0];
    dummy[4]=txb[0];
    dummy[5]=tyb[0];
    dummy[6]=nfilmb[sizeb-1]-nfilmb[0]+1; // number of passed plates [Firstplate..Lastplate]
    dummy[7]=eProb90;
    dummy[8]=eProb1;
    /*
    cout << "sizeb xb[0] yb[0] zb[0] txb[0] tyb[0] nplates eProb90 eProb1:" << endl;
    cout << "  " <<  dummy[0] << "  " << dummy[1] << "  " << dummy[2] << "  " << dummy[3] <<  "  " << dummy[4] << "  ";
    cout << "  " <<  dummy[5] << "  " << dummy[6] << "  " << dummy[7] << "  " << dummy[8] <<  "  " << dummy[9] << "  "<<endl;;
    */
    return dummy;
}




//_____________________________________________________________________________________________
void EdbShowerRec::PrintRecoShowerArray()
{
    // Print the array for the Reconstructed Showers.
    Log(2,"EdbShowerRec::PrintRecoShowerArray","EdbShowerRec::PrintRecoShowerArray");
    EdbTrackP* show=0;
    for (int i=0; i<GetRecoShowerArrayN(); ++i) {
        show=(EdbTrackP*)eRecoShowerArray->At(i);
        show->PrintNice();
    }
    cout << "Printed all " <<  GetRecoShowerArrayN()  << " showers of the eRecoShowerArray" << endl;
    Log(2,"EdbShowerRec::PrintRecoShowerArray","EdbShowerRec::PrintRecoShowerArray...done.");
    return;
}


//_____________________________________________________________________________________________

Double_t EdbShowerRec::DeltaThetaComponentwise(EdbSegP* s1,EdbSegP* s2)
{
    // Calculate Delta Theta of two segments.
    // Be aware that this DeltaTheta function returns the difference between the
    // component values of dTheta!!!
    // Acutally this function should be the normal way to calculate dTheta correctly...
    Double_t tx1,tx2,ty1,ty2;
    tx1=s1->TX();
    tx2=s2->TX();
    ty1=s1->TY();
    ty2=s2->TY();
    Double_t dt= TMath::Sqrt( (tx1-tx2)*(tx1-tx2) + (ty1-ty2)*(ty1-ty2) );
    return dt;
}

//_____________________________________________________________________________________________

Double_t EdbShowerRec::DeltaR_WithPropagation(EdbSegP* s,EdbSegP* stest)
{
    // Calculate Delta R of two segments.
    // If they have different Z positions, then first segment in the argument of the bracket
    // is propagated to the position of the second segment.
    if (s->Z()==stest->Z()) return TMath::Sqrt((s->X()-stest->X())*(s->X()-stest->X())+(s->Y()-stest->Y())*(s->Y()-stest->Y()));
    Double_t zorig;
    Double_t dR;
    zorig=s->Z();
    s->PropagateTo(stest->Z());
    dR=TMath::Sqrt( (s->X()-stest->X())*(s->X()-stest->X())+(s->Y()-stest->Y())*(s->Y()-stest->Y()) );
    s->PropagateTo(zorig);
    return dR;
}

//______________________________________________________________________________

Double_t EdbShowerRec::DeltaR_WithoutPropagation(EdbSegP* s,EdbSegP* stest)
{
    // Calculate Delta R of two segments using no z propagation.
    return TMath::Sqrt((s->X()-stest->X())*(s->X()-stest->X())+(s->Y()-stest->Y())*(s->Y()-stest->Y()));
}

//_____________________________________________________________________________________________

Bool_t EdbShowerRec::IsInConeTube(EdbSegP* TestingSegment, EdbSegP* StartingSegment, Double_t CylinderRadius, Double_t ConeAngle)
{
    // General Function which returns kTRUE if the Testing BaeTrack is in a cone defined
    // by the StartingBaseTrack. In case of same starting Z position, a max distance cut ( in Delta R)
    // of 20microns is assumed (i.e. the cone is not totally closed at the start).
    // In case of TestingSegment==StartingSegment this function should also correctly return kTRUE.

    if (gEDBDEBUGLEVEL>3) cout << "Bool_t EdbShowerAlg::IsInConeTube("<<TestingSegment<< ","<< StartingSegment <<","<< CylinderRadius << "," << ConeAngle << endl;
    if (gEDBDEBUGLEVEL>3) cout << "Bool_t EdbShowerAlg::IsInConeTube() Test Segment " << TestingSegment << " vs. Starting Segment " << StartingSegment << endl;

    // if (TestingSegment->Z()-StartingSegment->Z()<2.0) return kFALSE; //wrong...
    if (StartingSegment->Z()>TestingSegment->Z() ) return kFALSE;

    TVector3 x1(StartingSegment->X(),StartingSegment->Y(),StartingSegment->Z());
    TVector3 x2(TestingSegment->X(),TestingSegment->Y(),TestingSegment->Z());
    TVector3 direction_x1(StartingSegment->TX()*1273,StartingSegment->TY()*1273,1273);

    // u1 is the difference vector of the position!
    TVector3 u1=x2-x1;

    Double_t direction_x1_norm= direction_x1.Mag();
    Double_t cosangle=  (direction_x1*u1)/(u1.Mag()*direction_x1_norm);
    Double_t angle = TMath::ACos(cosangle);
    // This is the old version of angle calculation. It does not give the same results as in ROOT
    // when use TVector3.Angle(&TVector3). // For this IsInConeTube() we use therefore the ROOT calculation.
    angle=u1.Angle(direction_x1);

    // For the case where the two basetracks have same z position
    // the angle is about 90 degree so it makes no sense to calculate it...
    // therefore we set it artificially to zero:
    // if (TMath::Abs(StartingSegment->Z()-TestingSegment->Z())<2.0) {
    if (StartingSegment->Z()==TestingSegment->Z() ) {
        if (gEDBDEBUGLEVEL>3) cout << "Bool_t EdbShowerAlg::IsInConeTube()   Same Z position of TestingSegment and StartingSegment, Set angle artificially to zero" << endl;
        angle=0.0;

        // Check for position distance for 20microns if
        // Testing Segment is in same Z as StartingSegment
        if (gEDBDEBUGLEVEL>3) cout << "Bool_t EdbShowerAlg::IsInConeTube()   Check for position distance for 20microns if Testing Segment is in same Z as StartingSegment" << endl;
        if (gEDBDEBUGLEVEL>3) cout << "Bool_t EdbShowerAlg::IsInConeTube()   DeltaR_WithoutPropagation(StartingSegment,TestingSegment) = "<< DeltaR_WithoutPropagation(StartingSegment,TestingSegment) << endl;
        if (DeltaR_WithoutPropagation(StartingSegment,TestingSegment)<20) return kTRUE;
        if (DeltaR_WithoutPropagation(StartingSegment,TestingSegment)>=20) return kFALSE;
    }


    /// Outside if angle greater than ConeAngle (to be fulfilled for Cone and Tube in both cases)
    if (gEDBDEBUGLEVEL>3) cout << "Bool_t EdbShowerAlg::IsInConeTube()   Check if AngleVector now within the ConeAngleVector (<"<< ConeAngle<<"): " <<   angle << endl;
    if (angle>ConeAngle) {
        return kFALSE;
    }


    /// if angle smaller than ConeAngle, then you can differ between Tuberadius and CylinderRadius
    Double_t TubeDistance = 1.0/direction_x1_norm  *  ( (x2-x1).Cross(direction_x1) ).Mag();

    if (gEDBDEBUGLEVEL>3) cout << "Bool_t EdbShowerAlg::IsInConeTube()   Check if TestingSegment is now within the Tube (<"<< CylinderRadius<<"): " <<   TubeDistance << endl;

    if (TubeDistance>CylinderRadius) {
        return kFALSE;
    }

    return kTRUE;
}


//_____________________________________________________________________________________________


void EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks()
{
    // Fill the Initiator Basetrack Array by the linked_tracks file,
    // that is given in the  eFilename_LinkedTracks  String.
    // We check, if the file exists, and if there are tracks inside.
    // If so, we fill Initiator Basetrack Array with (a copy of ) these
    // tracks.

    Log(2,"EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks","EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks.");

    // -----------------
    // Check if file  eFilename_LinkedTracks  exits and if it has tracks inside....
    // -----------------  tODO
    cout << "EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks()   TODO::CHECK IF " << eFilename_LinkedTracks   <<" exists and if it has tracks inside!!!!" << endl;
    /// -----------------

    EdbSegP* segment=0;
    EdbSegP*    s2=0;
    EdbTrackP*  t=0;

    TFile * fil = new TFile(eFilename_LinkedTracks);
    fil->ls();
    TTree* tr= (TTree*)fil->Get("tracks");
    int nentr = 0;
    nentr =int(tr->GetEntries());

    //   check if tracks has entries: if so then ok, otherwise return directly:
    if (nentr>0) {
        ;
    }
    else {
        cout << "EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks()   tracks file has  NO  entries...Return." << endl;
        return;
    }

    TClonesArray *seg= new TClonesArray("EdbSegP",60);

    //   check if   eInBTArray  was created, if not do it now.:
    if (!eInBTArray) eInBTArray=new TObjArray(nentr);

    int nseg,n0,npl;
    cout << nentr << "  entries Total"<<endl;
    tr->SetBranchAddress("t."  , &t  );
    tr->SetBranchAddress("s"  , &seg  );
    tr->SetBranchAddress("nseg"  , &nseg  );
    tr->SetBranchAddress("n0"  , &n0  );
    tr->SetBranchAddress("npl"  , &npl  );

    for (int i=0; i<nentr; i++ ) {
        tr->GetEntry(i);

        // Take first BT of the tracks:
        s2=(EdbSegP*) seg->At(0);
        // Clone it (otherwise it will be overtaken by the next one...)
        segment=(EdbSegP*)s2->Clone();
        // Add it to array
        eInBTArray->Add(segment);
    }
    delete tr;
    delete fil;

    // Update Entry numbers.
    SetInBTArrayN(eInBTArray->GetEntries());

    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks   Printing InBTs" << endl;
    EdbSegP* testseg=0;
    for (int i=0; i<eInBTArrayN; i++ ) {
        if (gEDBDEBUGLEVEL>2) {
            testseg=(EdbSegP*)eInBTArray->At(i);
            testseg->PrintNice();
        }
    }

    eInBTArrayIsFilled=kTRUE;

    Log(2,"EdbShowerRec::Fill_eInBTArray_ByBaseTracksOf_eAli","EdbShowerRec::Fill_eInBTArray_ByBaseTracksOf_eAli.   After TCut: InBT N=%d", eInBTArrayN );
    Log(2,"EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks","EdbShowerRec::Fill_eInBTArray_ByLinkTracks_eFilename_LinkedTracks...done.");
    return;
}



//_____________________________________________________________________________________________



void EdbShowerRec::LoadEdbPVRec()
{
    // Function to load the EdbPatterVolumeReconstruction object into memory.
    // Shower reco needs to have full basetrack data to improve performance.
    // Because full basetrack data (corresponds to cp.root files) is memory extensive,
    // standard operations rely on linked tracks only.
    // Warning: some hardcoded strings are inside herebject for analysis use!

    Log(2,"EdbShowerRec::LoadEdbPVRec","Loads the full cp Basetracks into the  eAli  object.");
    // Loads the full cp Basetracks into the  eAli  object:
    if (!eAli) eAli=new EdbPVRec();

    // **** Data set initialization
    EdbDataProc* dset = new EdbDataProc("lnk.def");
    // **** Volume initialization
    dset->InitVolume(0);
    // **** Get Pattern Volume Reconstruction object
    eAli = dset->PVR();

    if (gEDBDEBUGLEVEL>2) eAli->Print();

    // Set private variables for the EdbPVrec Object now:
    eAliNpat=eAli->Npatterns();
    eAliNpatM1=eAliNpat-1;
    eNumberPlate_eAliPID=eAliNpatM1;

    eAliLoaded=kTRUE;
    Log(2,"EdbShowerRec::LoadEdbPVRec","Loads the full cp Basetracks into the  eAli  object.Done.");
}



//_____________________________________________________________________________________________


void EdbShowerRec::Transform_eAli( EdbSegP* InitiatorBT, Float_t ExtractSize=1000)
{
    //    Transform an EdbPVRec Object into a smaller EdbPVRec object.
    //    The length (number of plates/patterns) of eAli_Sub is not changed.
    //    Only XY-size (+- ExtractSize) and MC cuts (same MC or BG) around
    //    the InitiatorBT are applied.
    // 	  This method increases speed of processing, useful for large simulation MC
    //    datasets.
    //	  ( Code Detail Warning:
    //	    In case there are tooo many InitiatorBTs (order of 10k), the memory
    //      consumption be very high, the reason is that we cannot delete the
    //      small EdbPVRec objects build for the previous InitiatorBT, since the
    //      addresses of the segments will be deleted as well, segfault will appear.
    //      In this case the  EdbPattern::ExtractSubPattern function will not only
    //      clone the InitiatorBTs, but also the EdbPatterns, and the existing
    //      EdbSegPs in it.
    //      This method is more intended to use for the standalone program
    //      __ShowReco__  because in this program the EdbSegP objects are not used
    //      in its following code way, and then one can easyly delete the
    //      eAli_Sub object.
    //
    //      In this case, use the SetUseAliSub(0), then reconstruction will operate
    //      on the original EdbPVRec object.
    //    )
    //
    //    Test have shown that, the  Transform_eAli  with ExtractSubPattern  method
    //    gains highly in speed the more MC BTs are in the volume. For a
    //    sample consisting only of data, the only advantage is the loop over a
    //    smaller reconstruction area. If many MC BTs are in the volume, the advantage
    //    by cutting out also the other MC for one eAli_Sub object becomes really
    //    visible in time (but unfortunately also in memory).
    //
    //    So maybe a good tradeoff is like this:
    //      If ratio of MC/Data basetracks in the volume is significant then
    //      SetUseAliSub(1).
    //      If number of InitiatorBTs exceeds 5000 (number choosen by me after few
    //      tests) then we automatically switch to SetUseAliSub(0).
    //      For volumes containing only Data basetracks (i.e. MCEvt==-999) we
    //      do right from the start the operation on SetUseAliSub(0).
    /// IMPLEMENT THIS !!!
    //
    //    -----------------------------------------------------------------------------

    Log(3,"EdbShowerRec::Transform_eAli","Transform eAli to eAli_Sub.");

    // eUseAliSub=kTRUE; // Now possible to choose between use of full ali, or extracted ali.
    eAli_Sub_halfsize=ExtractSize;

    // IF TO RECREATE THE gALI_SUB in RECONSTRUCTION  or if to use gAli global (slowlier but maybe memory leak safe).
    if (eUseAliSub==0)
    {
        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Transform_eAli   UseAliSub==0  No new eAli_Sub created. Use eAli instead. "<<endl;
        eAli_Sub=eAli;
        if (gEDBDEBUGLEVEL>3) eAli_Sub->Print();
        return;
    }

    Int_t npat;
    npat = eAli->Npatterns();  //number of plates

    // has to be deleted in some part of the script outside this function...
    // Dont forget , otherwise memory heap overflow!
    /// DEBUG       if (eAli_Sub) { delete eAli_Sub;eAli_Sub=0;}
    /// original, but keeps not adresses of segment in eAli; if it is deleted, then segfault will appear...
    /// But if not removed, memory can increase very quickly.
    /// In this case, use the SetUseAliSub(0), then reconstruction will operate on the
    /// on the original EdbPVRec Object.
    if (NULL!=eAli_Sub) {
        ;
    } /// do nothing now... let it live... delete eAli_Sub;eAli_Sub=0;}
    eAli_Sub = new EdbPVRec();

    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Transform_eAli   eAli_Sub  created."<<endl;

    // Create SubPattern objects
    EdbSegP* ExtrapolateInitiatorBT=0;
    ExtrapolateInitiatorBT = (EdbSegP*)InitiatorBT->Clone();

    Int_t InitiatorBTMCEvt=InitiatorBT->MCEvt();

    // Create Variables For ExtractSubpattern boundaries
    Float_t mini[5];
    Float_t maxi[5];
    //Float_t ExtractSize=1000;// now in fucntion header
    mini[0]=ExtrapolateInitiatorBT->X()-ExtractSize;
    mini[1]=ExtrapolateInitiatorBT->Y()-ExtractSize;
    maxi[0]=ExtrapolateInitiatorBT->X()+ExtractSize;
    maxi[1]=ExtrapolateInitiatorBT->Y()+ExtractSize;
    mini[2]=-0.7;
    mini[3]=-0.7;
    mini[4]=0.0;
    maxi[2]=0.7;
    maxi[3]=0.7;
    maxi[4]=100.0;

    EdbPattern* singlePattern;
    Float_t ExtrapolateInitiatorBT_zpos_orig=ExtrapolateInitiatorBT->Z();

    //InitiatorBT->PrintNice();

    // Add the subpatterns in a loop for the plates:
    // in reverse ordering.due to donwstream behaviour (!):
    // (Only downstream is supported now...)
    for (Int_t ii=0; ii<npat; ++ii) {

        Float_t zpos=eAli->GetPattern(ii)->Z();
        //cout << " Doing pattern " << ii << " at Z postion = "  << zpos << endl;

        ExtrapolateInitiatorBT->PropagateTo(zpos);

        mini[0]=ExtrapolateInitiatorBT->X()-ExtractSize;
        mini[1]=ExtrapolateInitiatorBT->Y()-ExtractSize;
        maxi[0]=ExtrapolateInitiatorBT->X()+ExtractSize;
        maxi[1]=ExtrapolateInitiatorBT->Y()+ExtractSize;

        singlePattern=(EdbPattern*)eAli->GetPattern(ii)->ExtractSubPattern(mini,maxi,InitiatorBTMCEvt);
        // This sets PID() analogue to (upstream), nut not PID of the BTs !
        singlePattern-> SetID(eAli->GetPattern(ii)->ID());
        // This sets PID() analogue to (upstream), nut not PID of the BTs !
        singlePattern-> SetPID(eAli->GetPattern(ii)->PID());

        eAli_Sub->AddPattern(singlePattern);

        // Propagate back...! (in order not to change the original BT)
        ExtrapolateInitiatorBT->PropagateTo(ExtrapolateInitiatorBT_zpos_orig);

        //cout << "singlePattern: N segments filled:  " << singlePattern->N() << endl;
    }

    delete ExtrapolateInitiatorBT;
    eAli_SubNpat=eAli_Sub->Npatterns();  //number of plates

    if (gEDBDEBUGLEVEL>2) eAli_Sub->Print();
    return;
}

//_____________________________________________________________________________________________

Float_t EdbShowerRec::GetEdbPVRecNSeg(EdbPVRec* ali, Int_t mode)
{
    // Return Number of Segments (=BTs) in the EdbPVRec volume.
    // The mode number is used as following:
    //    mode = 0 : return ALL BTs.
    //    mode = 1 : return BTs with s->MCEvt() > 0    (MC)
    //    mode = 2 : return BTs with s->MCEvt() < 0    (DATA)
    //    mode = 3 : return fraction nTotalMC/(nTotal);

    if (NULL==ali) return 0;
    Int_t nTotal=0;
    Int_t nTotalMC=0;
    Int_t nTotalDATA=0;
    Int_t npat=ali->Npatterns();
    for (Int_t ii=0; ii<npat; ++ii) {
        Int_t nseg=ali->GetPattern(ii)->N();
        nTotal+=nseg;
        for (Int_t jj=0; jj<nseg; ++jj) {
            EdbSegP* seg=ali->GetPattern(ii)->GetSegment(jj);
            if (seg->MCEvt()>0) ++nTotalMC;
            if (seg->MCEvt()<0) ++nTotalDATA;
        }
    }
    Float_t ratio=(Float_t)(nTotalMC)/(Float_t)(nTotal);
//   cout << "nTotalDATA basetracks in volume: " << nTotalDATA << endl;
//   cout << "nTotalMC basetracks in volume: " << nTotalMC << endl;
//   cout << "nTotal basetracks in volume: " << nTotal << endl;
//   cout << "ratio basetracks in volume: " << ratio << endl;

    if (mode == 1) return (Float_t)nTotalMC;
    if (mode == 2) return (Float_t)nTotalDATA;
    if (mode == 3) return ratio;
    return (Float_t)nTotal;

}


//_____________________________________________________________________________________________


void EdbShowerRec::CheckEdbPVRecSetUseAliSub()
{
    // Small function that test the ratio of MC/Data basetracks in the volume.
    // If ratio is high (>20%), and number of InitiatorBTN is mediate (<1000),
    // then we will use the SetUseAliSub(1) option, which means:
    //   -faster reconstruction, but higher memory consumption.
    // If not we switch to the standard reco method SetUseAliSub(0):
    //   -lower reconstruction, but stable memory consumption.

    cout << "void EdbShowerRec::CheckEdbPVRecSetUseAliSub()" << endl;
    cout << "void EdbShowerRec::CheckEdbPVRecSetUseAliSub() GetUseAliSub() = " << GetUseAliSub() << endl;

    if (gEDBDEBUGLEVEL>2) {
        cout << "void EdbShowerRec::CheckEdbPVRecSetUseAliSub() " <<  eAli << endl;
        cout << "void EdbShowerRec::CheckEdbPVRecSetUseAliSub() " <<  GetEdbPVRecNSeg(eAli,3) << endl;
        cout << "void EdbShowerRec::CheckEdbPVRecSetUseAliSub() " <<  eInBTArrayN << endl;
    }

    if (NULL==eAli) return;
    if (GetEdbPVRecNSeg(eAli,3)<0.20) return;
    if (eInBTArrayN>1000) return;

    // If we arrived here, it might be time to switch to the fast reconstruction!
    SetUseAliSub(1);

    cout << "EdbShowerRec::CheckEdbPVRecSetUseAliSub()   We arrived here, so it might be time to switch to the fast reconstruction!" << endl;
    cout << "EdbShowerRec::CheckEdbPVRecSetUseAliSub()   GetUseAliSub() = " << GetUseAliSub() << endl;

    return;
}


//_____________________________________________________________________________________________


void EdbShowerRec::Execute()
{
    cout << "EdbShowerRec::Execute()" << endl;
    cout << "EdbShowerRec::Execute   DOING MAIN SHOWER RECONSTRUCTION HERE" << endl;
    cout << "EdbShowerRec::Execute   USE Standard CONETUBE ALGO with Standard Parameters." << endl;

    // Check if the eRecoShowerArray exists, if not create it;
    if (!eRecoShowerArray) eRecoShowerArray=new TObjArray();
    //cout << "eRecoShowerArray->GetEntries() " << eRecoShowerArray->GetEntries()  << endl;

    // Check the use of eAliSub or not:
    CheckEdbPVRecSetUseAliSub();


    EdbSegP* InBT;
    EdbSegP* Segment;
    EdbTrackP* RecoShower;

    Bool_t    StillToLoop=kTRUE;
    Int_t     ActualPID;
    Int_t     newActualPID;
    Int_t     STEP=-1;
    Int_t     NLoopedPattern=0;

    Int_t     eLastPlate_eAliPID=0;
    Int_t     eFirstPlate_eAliPID=0;
    Float_t   eLastPlate_eAliPIDZ=0;
    Float_t   eFirstPlate_eAliPIDZ=0;

    // Already defined as private variable:
    eAliNpatM1=eAliNpat-1;
    eNumberPlate_eAliPID=eAliNpatM1;
    if (gEDBDEBUGLEVEL>2) {
        cout << "  eAliNpat  = " << eAliNpat << endl;
    }


    // "First Plate" means always the one most   upstream, i.e. the lowest  Z position
    // "Last  Plate" means always the one most downstream, i.e. the highest Z position
    if  (eAli->GetPattern(0)->Z() > eAli->GetPattern(eAliNpatM1)->Z()) {
        eLastPlate_eAliPID=eAli->GetPattern(0)->PID();
        eFirstPlate_eAliPID=eAli->GetPattern(eAliNpatM1)->PID();
        eLastPlate_eAliPIDZ=eAli->GetPattern(0)->Z();
        eFirstPlate_eAliPIDZ=eAli->GetPattern(eAliNpatM1)->Z();
    }
    else {
        eLastPlate_eAliPID=eAli->GetPattern(eAliNpatM1)->PID();
        eFirstPlate_eAliPID=eAli->GetPattern(0)->PID();
        eLastPlate_eAliPIDZ=eAli->GetPattern(eAliNpatM1)->Z();
        eFirstPlate_eAliPIDZ=eAli->GetPattern(0)->Z();
    }

    if (gEDBDEBUGLEVEL>1) {
        cout << "EdbShowerRec::Execute()  eAli->GetPattern(0)->Z() > eAli->GetPattern(eAliNpatM1)->Z() " << eAli->GetPattern(0)->Z()<< " " << eAli->GetPattern(eAliNpatM1)->Z() << endl;
        cout << "EdbShowerRec::Execute()  eLastPlate_eAliPID  = " << eLastPlate_eAliPID << endl;
        cout << "EdbShowerRec::Execute()  eFirstPlate_eAliPID  = " << eFirstPlate_eAliPID << endl;
        cout << "EdbShowerRec::Execute()  eLastPlate_eAliPIDZ  = " << eLastPlate_eAliPIDZ << endl;
        cout << "EdbShowerRec::Execute()  eFirstPlate_eAliPIDZ  = " << eFirstPlate_eAliPIDZ << endl;
        cout << endl;
    }
// DEBUG!!!!

    if (eFirstPlate_eAliPID<eLastPlate_eAliPID) STEP=-1;
    if (eFirstPlate_eAliPID>eLastPlate_eAliPID) STEP=1;
    if (gEDBDEBUGLEVEL>1) cout << "EdbShowerRec::Execute   STEP for patternloop direction =  " << STEP << endl;

    //--- Loop over InBTs:
    if (gEDBDEBUGLEVEL>1) {
        cout << "EdbShowerRec::Execute   Loop over InBTs:" << endl;
        cout << "EdbShowerRec::Execute   ("<< eInBTArrayN << "):" << endl;
    }



    // Since eInBTArray is usually filled in ascending ordering by zpositon
    // We use the descending loop to begin with BT with lowest z first.
    for (Int_t i=eInBTArrayN-1; i>=0; --i) {
//  for (Int_t i=eInBTArrayN-1; i==eInBTArrayN-1; --i) {

        // CounterOutPut
        if (gEDBDEBUGLEVEL==2) if ((i%1)==0) cout << eInBTArrayN <<" InBT in total, still to do:"<<Form("%4d",i)<< "\r\r\r\r"<<flush;
        if (gEDBDEBUGLEVEL==3) cout << eInBTArrayN <<" InBT in total, still to do:"<< i << endl;

        //-----------------------------------
        // 1) Make local_gAli with cut parameters:
        //-----------------------------------

        // Create new Shower Object for storage;
        // For the moment this is treated as an "EdbTrackP" object...
        RecoShower = new EdbTrackP(*(EdbSegP*)eInBTArray->At(i));

        // Get InitiatorBT from eInBTArray
        InBT=(EdbSegP*)eInBTArray->At(i);
        if (gEDBDEBUGLEVEL>2) {
            cout << " " << endl;
            cout << "EdbShowerRec::Execute() Start reconstruction for Initiator BaseTrack:" <<  endl;
            InBT->PrintNice();
        }

        // Add InBT to RecoShower:
        // This has to be done, since by definition the first BT in the RecoShower is the InBT.
        // Otherwise, also the definition of shower axis and transversal profiles is wrong!
        RecoShower -> AddSegment(InBT);
        // Attention: If you have a EdbTrackP object and you add a segment, then
        // you have to manually set Counters for N(), NPl(). It is not done auto.
        RecoShower -> SetCounters();
        //cout << "Segment  (InBT) " << Segment << " was added to RecoShower." <<  endl;

        // Transform (make size smaller, extract only events having same MC) the  eAli  object:
        Transform_eAli(InBT);
        if (gEDBDEBUGLEVEL>3) eAli_Sub->Print();
        Int_t npat_int=0;
        Int_t npat_total=0;
        Int_t npatN=0;
        // Get greatest Z Value:
        Double_t greatestZValueofeAliSub=TMath::Max(eAli_Sub->GetPattern(0)->Z(),eAli_Sub->GetPattern(eAli_Sub->Npatterns()-1)->Z());
        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Execute--- --- greatest Z value of eAliSub= " << greatestZValueofeAliSub << endl;


        //-----------------------------------
        // 2) Loop over (whole) eAli, check BT for Cuts
        // eAli_Sub
        // 				We have to sort the Initiator PID and pattern PID,
        // 				since the sample of Bari has shown that these two PIDs maybe
        // 				not necessaryly be the same!!!
        //-----------------------------------


        EdbPattern* patternActualPID = eAli_Sub->GetPatternByZ(InBT->Z());
        ActualPID     = patternActualPID->PID();
        newActualPID   = ActualPID ;
        //InBT->PrintNice();
        //cout << endl;

        while (StillToLoop) {
            if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Execute--- --- Doing patterloop for PID: " << ActualPID << ", for patterns Z position: " << eAli_Sub->GetPatternByPID(ActualPID)->Z() << endl;

            patternActualPID = eAli_Sub->GetPatternByPID(ActualPID);

            // This reco Alg we take only BTs that come downstream the InBT:
            // (should now never hapen anymore).
            if (patternActualPID->Z() < InBT->Z() ) {
                cout << "WARNING : patternActualPID->Z() < InBT->Z() " << endl;
                continue;
            }

            for (Int_t btloop_cnt=0; btloop_cnt<patternActualPID->GetN(); ++btloop_cnt) {

                Segment = (EdbSegP*)patternActualPID->GetSegment(btloop_cnt);

                if (gEDBDEBUGLEVEL>4) {
                    cout << "Checking Segment " << btloop_cnt << endl;
                    Segment->PrintNice();
                }


                // Now apply cut conditions: Cut for Quality: --------------------
                if ( (eUseQualityPar==kTRUE) && (Segment->Chi2()>Segment->W()*eQualityPar[0]-eQualityPar[1]) ) {
                    continue;
                }


                // Now apply cut conditions: OI OfficialImplementation Alg  --------------------
                if ( !IsInConeTube(Segment, InBT, eAlgoParameterConeRadius , eAlgoParameterConeAngle) ) continue;
                //cout << "Segment " << Segment << " passed IsInConeTube"<<endl;
                if ( !FindPrecedingBTs(Segment, InBT, eAli_Sub, RecoShower)) continue;
                if (gEDBDEBUGLEVEL>4) cout << "Segment " << Segment << " passed FindPrecedingBTs"<<endl;
                // end of    cut conditions: OI OfficialImplementation Alg  --------------------

                // If we arrive here, Basetrack  Segment  has passed criteria
                // and is then added to the RecoShower:
                // Check if its not the InBT which is already added:
                if (Segment->X()==InBT->X()&&Segment->Y()==InBT->Y()) {
                    ;    // is InBT, do nothing;
                }
                else {
                    RecoShower -> AddSegment(Segment);
                    // Attention: If you have a EdbTrackP object and you add a segment, then
                    // you have to manually set Counters for N(), NPl(). It is not done auto.
                    RecoShower -> SetCounters();
                }
                if (gEDBDEBUGLEVEL>4) cout << "Segment  " << Segment << " was added at  &Segment : " << &Segment  <<  endl;

            } // of btloop_cnt


            if (gEDBDEBUGLEVEL>2) {
                cout << "EdbShowerRec::Execute--- --- ActualPID= " << ActualPID << "  done. Reconstructed shower has up to now: " << RecoShower->N()  << " Segments." << endl;
                cout << "EdbShowerRec::Execute--- --- NLoopedPattern = " << NLoopedPattern << "  until now. " <<  endl;
            }


            // Calc BT density around shower:
            EdbPattern* pat_interim = patternActualPID;
            CalcTrackDensity(pat_interim,eAli_Sub_halfsize,npat_int,npat_total,npatN);

            //------------
            /// CALCULATE NEWACTUALPID NUMBER BY GETPATTERNSUCCEDING:
            EdbPattern* pat_ActualSucceeding;
            pat_ActualSucceeding = eAli_Sub->GetPatternSucceding( patternActualPID ) ;


            if (pat_ActualSucceeding!=NULL) {

                // Set New Values:
                newActualPID = pat_ActualSucceeding->PID();
                ++NLoopedPattern;

                //cout <<"patternActualPID          : at,PID,Z " << patternActualPID << " " << patternActualPID->PID() << "  " << patternActualPID->Z() << endl;
                //cout <<"pat_ActualSucceeding: at,PID,Z " << pat_ActualSucceeding << " " << pat_ActualSucceeding->PID() << "  " << pat_ActualSucceeding->Z() << endl;
            }


            // DAU user solution for stopping the loop:
            if (pat_ActualSucceeding==NULL) {
                StillToLoop=kFALSE;
                if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Execute--- ---Stopp Loop since: pat_ActualSucceeding==NULL"<<endl;
            }

            if (gEDBDEBUGLEVEL>3) {
                cout << "EdbShowerRec::Execute--- --- StillToLoop= " << StillToLoop << endl;
                cout << "EdbShowerRec::Execute--- ---ActualPID= " << ActualPID <<" Finished. Take newActualPID = " << newActualPID << " now. (this is end of  while (StillToLoop) ) " << endl;
            }

            ActualPID=newActualPID;
            //cout << endl << endl;
        } // of // while (StillToLoop)


        if (gEDBDEBUGLEVEL>2) cout << "Finshed  __while (StillToLoop)__  . Now Adding reconstructed shower to the array (only in case it has two or more Basetracks) of reconstructed showers ..." << endl;
        if (gEDBDEBUGLEVEL>2) cout << "Shower Has  ..." << RecoShower->N() << "  basetracks in it." << endl;

        // Add Shower Object to Shower Reco Array.
        // Not, if its empty or containing only one BT:
        if (RecoShower->N()>1) eRecoShowerArray->Add(RecoShower);
        // cout << "DEBUG_TEST    eRecoShowerArray->Add(RecoShower);" << endl;
        // cout << "is eRecoShowerArray existing anyway??? " << eRecoShowerArray << endl;

        // Set back loop values:
        StillToLoop=kTRUE;
        NLoopedPattern=0;

    } // of  //   for (Int_t i=eInBTArrayN-1; i>=0; --i) {


    // Set new value for  eRecoShowerArrayN  (may now be < eInBTArrayN).
    SetRecoShowerArrayN(eRecoShowerArray->GetEntries());

    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::eRecoShowerArray() Entries: " << eRecoShowerArray->GetEntries() << endl;
    if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::Execute()...done." << endl;

    if (gEDBDEBUGLEVEL>2) PrintRecoShowerArray();


    ///========================   WRITE TO FILE..... ======================================================
    ///========================   TEMPORARY SOLUTION TO ACCES FUNCTION NEURALNET AND ENERGY .....==========
    ///========================  Transfer ..... ======================================================
    TransferShowerObjectArrayIntoEntryOfTreebranchShowerTree(eRecoShowerArray);
    NeuralNet();
    CalculateEnergyValues();



    cout << "======================== " << eQualityPar[0]  << endl;
    cout << "======================== " << eQualityPar[1]  << endl;
    return;
}

//______________________________________________________________________________




//______________________________________________________________________________

Bool_t EdbShowerRec::FindPrecedingBTs(EdbSegP* s, EdbSegP* InBT, EdbPVRec *gAli, EdbTrackP* shower)
{
    //cout << "EdbShowerRec::FindPrecedingBTs   Find BTs to be connected with the Test BT:" << endl;

    EdbSegP* s_TestBT;
    Int_t nentries=shower->N();
    Double_t dZ;

    // In case the shower is empty we do not have to search for a preceeding BT:
    if (nentries==0) return kTRUE;

    // We do not test BaseTracks which are before the initiator BaseTrack!
    if (s->Z()<InBT->Z()) {
        cout << "--- --- EdbShowerRec::FindPrecedingBTs(): s->Z()<InBT->Z()..: We do not test BaseTracks which are before the initiator BaseTrack!  return kFALSE;" << endl;
        return kFALSE;
    }

    // For the very first Z position we do not test
    // if testBT has Preceeders, only if it it has a BT around (case for e+e- coming from gammma):
    // Take 50microns and 80mrad in (dR/dT) around.
    // This does not affect the normal results, but helps for
    // events which may have a second BT close to InBT (like in e+e-)
    if (TMath::Abs(s->Z()-InBT->Z())<2.0) {
        //cout << "--- --- EdbShowerRec::FindPrecedingBTs(): s->Z()~InBT->Z()..: DeltaThetaComponentwise(s, InBT)= " << DeltaThetaComponentwise(s, InBT) <<  " <<  DeltaR_WithoutPropagation(s, InBT) = " << DeltaR_WithoutPropagation(s, InBT) << endl;
        if (DeltaThetaComponentwise(s, InBT) < 0.08 && DeltaR_WithoutPropagation(s, InBT) < 50.0 ) return kTRUE;
    }

    if (gEDBDEBUGLEVEL>3) cout << "--- --- EdbShowerRec::FindPrecedingBTs(): Testing " << nentries << " entries of the shower." << endl;

    for (Int_t i=nentries-1; i>=0; --i) {
        s_TestBT = (EdbSegP*)( shower->GetSegment(i) );

        if (gEDBDEBUGLEVEL>3) cout << "--- --- EdbShowerRec::FindPrecedingBTs(): Do  s_TestBT->ID()  s->ID()  s_TestBT->MCEvt() s_TestBT->Z()  s->Z() "<< s_TestBT->ID() << " " << s->ID() << " " << s_TestBT->MCEvt() <<"  " << s_TestBT->Z() << "  " <<  s->Z()  <<endl;
        if (gEDBDEBUGLEVEL>3) s_TestBT->PrintNice();
        if (gEDBDEBUGLEVEL>3) s->PrintNice();


        dZ=TMath::Abs(s_TestBT->Z()-s->Z());
        if (dZ<30) continue;                  // Exclude the case of same Zpositions...
        if (dZ>(eAlgoParameterNPropagation*1273.0)+650.0) continue;     // Exclude the case of more than eAlgoParameterNPropagation plates before...

        if (gEDBDEBUGLEVEL>3) cout << "--- --- EdbShowerRec::FindPrecedingBTs(): Checking dT,dR and dZ for i:  " << i << "  " << DeltaThetaComponentwise(s, s_TestBT)  << "  " << DeltaR_WithPropagation(s, s_TestBT) << "  "<<dZ << endl;

        if (DeltaThetaComponentwise(s, s_TestBT) >  eAlgoParameterConnectionDT ) continue;
        if (DeltaR_WithPropagation(s, s_TestBT) > eAlgoParameterConnectionDR ) continue;

        if (gEDBDEBUGLEVEL>3) cout << "--- --- EdbShowerRec::FindPrecedingBTs(): Checking dT,dR and dZ for i:  " << i << "  " << DeltaThetaComponentwise(s, s_TestBT)  << "  " << DeltaR_WithPropagation(s, s_TestBT) << "  "<<dZ << "   ok!"<<endl;
        return kTRUE;
    }
    //---------------------------------------------
    return kFALSE;
}




//------------------------------------------------------------------------------------------------------

void EdbShowerRec::ClearInBTArray() {
    if (eInBTArrayN>0) eInBTArray->Clear();
    eInBTArrayN=0;
    return;
}
void EdbShowerRec::ClearRecoShowerArray() {
    if (eRecoShowerArrayN>0) eRecoShowerArray->Clear();
    eRecoShowerArrayN=0;
    return;
}

//------------------------------------------------------------------------------------------------------


//-------------------------------------------------------------------------------------------
void EdbShowerRec::TransferShowerObjectArrayIntoEntryOfTreebranchShowerTree(TObjArray* showrecoarray)
{
    Log(2, "EdbShowerRec", "--- void* TransferShowerObjectArrayIntoEntryOfTreebranchShowerTree() ---Interim SOLUTIION !!!! ");
    TFile* interim_fileout = new TFile("shower1.root","RECREATE");
    // EdbShowerRec::T-S-O-A-I-E-O-T-S-T

    // Interim SOLUTIION !!!!  shower tree definition

    //---------------------------------------------------------------------------------
    //-     VARIABLES: shower_  "treebranch"  reconstruction
    //---------------------------------------------------------------------------------
    Int_t shower_number_eventb, shower_sizeb, shower_isizeb,shower_showerID;
    Int_t shower_sizeb15, shower_sizeb20, shower_sizeb30;
    Float_t shower_energy_shot_particle;
    Float_t shower_xb[9999];
    Float_t shower_yb[9999];
    Float_t shower_zb[9999];
    Float_t shower_txb[9999];
    Float_t shower_tyb[9999];
    Float_t shower_deltarb[9999];
    Float_t shower_deltathetab[9999];
    Float_t shower_deltaxb[9999];
    Float_t shower_deltayb[9999];
    Int_t   shower_nfilmb[9999];
    //Int_t   shower_ngrainb[9999];
    Float_t shower_chi2btkb[9999];
    Int_t shower_ntrace1simub[9999]; // MCEvt
    Int_t shower_ntrace2simub[9999]; // s->W()
    Float_t shower_ntrace3simub[9999]; // s->P()
    Int_t shower_ntrace4simub[9999]; // s->Flag()
    Float_t shower_tagprimary[9999];
    Int_t   shower_idb[9999];
    Int_t   shower_plateb[9999];
    Float_t shower_deltasigmathetab[58];
    Int_t   shower_numberofilms;
    Float_t shower_purityb,shower_E_MC,shower_trackdensb;
    //---------------------------------------------------------------------------------

    TTree* interim_treesaveb = new TTree("treebranch","tree of branchtrack");
    interim_treesaveb->Branch("number_eventb",&shower_number_eventb,"number_eventb/I");
    interim_treesaveb->Branch("sizeb",&shower_sizeb,"sizeb/I");
    interim_treesaveb->Branch("sizeb15",&shower_sizeb15,"sizeb15/I");
    interim_treesaveb->Branch("sizeb20",&shower_sizeb20,"sizeb20/I");
    interim_treesaveb->Branch("sizeb30",&shower_sizeb30,"sizeb30/I");
    interim_treesaveb->Branch("E_MC",&shower_E_MC,"E_MC/F");
    interim_treesaveb->Branch("isizeb",&shower_isizeb,"isizeb/I");
    interim_treesaveb->Branch("showerID",&shower_showerID,"showerID/I");
    interim_treesaveb->Branch("idb",shower_idb,"idb[sizeb]/I");
    interim_treesaveb->Branch("plateb",shower_plateb,"plateb[sizeb]/I");
    interim_treesaveb->Branch("xb",shower_xb,"xb[sizeb]/F");
    interim_treesaveb->Branch("yb",shower_yb,"yb[sizeb]/F");
    interim_treesaveb->Branch("zb",shower_zb,"zb[sizeb]/F");
    interim_treesaveb->Branch("txb",shower_txb,"txb[sizeb]/F");
    interim_treesaveb->Branch("tyb",shower_tyb,"tyb[sizeb]/F");
    interim_treesaveb->Branch("nfilmb",shower_nfilmb,"nfilmb[sizeb]/I");
    interim_treesaveb->Branch("ntrace1simub",shower_ntrace1simub,"ntrace1simu[sizeb]/I");
    interim_treesaveb->Branch("ntrace2simub",shower_ntrace2simub,"ntrace2simu[sizeb]/I"); // s->W()
    interim_treesaveb->Branch("ntrace3simub",shower_ntrace3simub,"ntrace3simu[sizeb]/F"); // s->P()
    interim_treesaveb->Branch("ntrace4simub",shower_ntrace4simub,"ntrace4simu[sizeb]/I"); // s->Flag()
    interim_treesaveb->Branch("chi2btkb",shower_chi2btkb,"chi2btkb[sizeb]/F");
    interim_treesaveb->Branch("deltarb",shower_deltarb,"deltarb[sizeb]/F");
    interim_treesaveb->Branch("deltathetab",shower_deltathetab,"deltathetab[sizeb]/F");
    interim_treesaveb->Branch("deltaxb",shower_deltaxb,"deltaxb[sizeb]/F");
    interim_treesaveb->Branch("deltayb",shower_deltayb,"deltayb[sizeb]/F");
    interim_treesaveb->Branch("tagprimary",shower_tagprimary,"tagprimary[sizeb]/F");
    interim_treesaveb->Branch("purityb",&shower_purityb,"purityb/F");
    interim_treesaveb->Branch("trackdensb",&shower_trackdensb,"trackdensb/F");

    Int_t showrecoarrayN=showrecoarray->GetEntries();
    cout << "showrecoarray->GetEntries(); " <<   showrecoarray->GetEntries() << endl;


    // LOOP OVER ENTRIES FROM showrecoarray:
    for (Int_t ooo=0; ooo< showrecoarrayN; ++ooo) {
        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::T-S-O-A-I-E-O-T-S-T: Doing entry " << ooo << endl;

        // Get Track from  showrecoarray
        EdbTrackP* segarray = (EdbTrackP*)showrecoarray->At(ooo);
        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::T-S-O-A-I-E-O-T-S-T: segarray = " << segarray  << endl;
        if (gEDBDEBUGLEVEL>2) cout << "EdbShowerRec::T-S-O-A-I-E-O-T-S-T: segarray->N(),N0(): " << segarray->N() << " " << segarray->N0() << endl;
        if (gEDBDEBUGLEVEL>2) segarray->PrintNice();

        EdbSegP* seg;
        EdbSegP* Inseg;
        Int_t diff_pid;
        Float_t min_shower_deltathetab=99999;
        Float_t min_shower_deltar=99999;
        Float_t test_shower_deltathetab=99999;
        Float_t test_shower_deltar=99999;
        Float_t test_shower_deltax,test_shower_deltay;

        Float_t shower_sizebNHELP=0;
        Float_t shower_sizebMCNHELP=0;

        Float_t extrapol_x,extrapol_y, extrapo_diffz;

        // Initialize arrays...
        shower_sizeb15=0;
        shower_sizeb20=0;
        shower_sizeb30=0;
        shower_sizeb=0;
        shower_energy_shot_particle=0.0;
        shower_numberofilms=0;
        for (int ii=0; ii<9999; ii++)  {
            shower_xb[ii]=0;
            shower_yb[ii]=0;
            shower_zb[ii]=0;
            shower_txb[ii]=0;
            shower_tyb[ii]=0;
            shower_nfilmb[ii]=0;
            shower_tagprimary[ii]=0;
            shower_ntrace1simub[ii]=0;
            shower_ntrace2simub[ii]=0;
            shower_ntrace3simub[ii]=0;
            shower_ntrace4simub[ii]=0;
            shower_deltaxb[ii]=0;
            shower_deltayb[ii]=0;
            shower_chi2btkb[ii]=0;
            shower_idb[ii]=0;
            shower_plateb[ii]=0;
        }
        for (int i=1; i<59; ++i) {
            shower_deltasigmathetab[i]=0;
        }

        // Part To calculate the TransfereedVariables....
        shower_sizeb=segarray->N();
        Inseg=(EdbSegP*)segarray->GetSegment(0);
        shower_energy_shot_particle=Inseg->P();
        if (gEDBDEBUGLEVEL>2)  cout << "EdbShowerRec::T-S-O-A-I-E-O-T-S-T   shower_sizeb= "<< shower_sizeb << endl;
        if (gEDBDEBUGLEVEL>2)  cout << "EdbShowerRec::T-S-O-A-I-E-O-T-S-T   Inseg->PrintNice()= " << endl;
        if (gEDBDEBUGLEVEL>2)  Inseg->PrintNice();



        if (gEDBDEBUGLEVEL>3) cout << "--- --- ---------------------"<<endl;
        //-------------------------------------
        for (int ii=0; ii<shower_sizeb; ii++)  {
            if (gEDBDEBUGLEVEL>3) cout << "EdbShowerRec::T-S-O-A-I-E-O-T-S-T   Doing Segment "<< ii << endl;

            if (ii>=9999) {
                cout << "EdbShowerRec::T-S-O-A-I-E-O-T-S-T: WARNING: shower_sizeb ( " << shower_sizeb<< ") greater than SHOWERARRAY.   Set sizeb to 9999 and  Stop filling!."<<endl;
                shower_sizeb=9999;
                break;
            }
//       seg=(EdbSegP*)segarray->At(ii); // works when TObjArray
            seg=(EdbSegP*)segarray->GetSegment(ii);  // works when EdbTrackP
            if(!seg) cout<<"ERRORE SEGMENTO NON TROVATO"<<endl;
            //-------------------------------------
            shower_xb[ii]=seg->X();
            shower_yb[ii]=seg->Y();
            shower_txb[ii]=seg->TX();
            shower_tyb[ii]=seg->TY();
            shower_zb[ii]=seg->Z();
            shower_chi2btkb[ii]=seg->Chi2();
            shower_deltathetab[ii]=0.5;
            shower_deltarb[ii]=200;
            shower_tagprimary[ii]=0;
            if (ii==0) shower_tagprimary[ii]=1;
            shower_isizeb=1; // always 1, not needed anymore
            if (seg->MCEvt()>0) {
                shower_ntrace1simub[ii]=seg->MCEvt();
                shower_ntrace2simub[ii]=seg->W();
                shower_ntrace3simub[ii]=seg->P();
                shower_ntrace4simub[ii]=seg->Flag();
            }
            else {
                // keep the seg->BT settings for BG:
                // that come out of "normal" scanned data from fedra:
                // shower_ntrace1simub=-999
                // shower_ntrace2simub=seg->W();
                // shower_ntrace3simub=-999
                // shower_ntrace4simub=0
                shower_ntrace1simub[ii]=-999;
                shower_ntrace2simub[ii]=seg->W();
                shower_ntrace3simub[ii]=-999;
                shower_ntrace4simub[ii]=0;
            }
            shower_idb[ii]=seg->ID();
            shower_plateb[ii]=seg->PID();
            shower_E_MC=seg->P();
            //-------------------------------------
            // PUT HERE:  deltarb,deltarb, nflimb, sizeb15......
            diff_pid=TMath::Abs( Inseg->PID()-seg->PID() )+1;
            // (does this work for up/downsream listing??)
            // (yes, since InBT->PID is also changed.)
            // but works only if the gAli Object has no missing plates
            // otherwise f.e. PID(1) and PID(2) are not necessaryly abay by dZ=1273
            // (could be Z(1)=1273 and Z(2)=3900...)

            // Calc pur:
            shower_sizebNHELP++;
            // Increase the shower_sizebMCNHELP variable ONLY and ONLY if there is MC segment
            if (seg->MCEvt()>0) shower_sizebMCNHELP++;


            // InBT:
            if (ii==0) {
                shower_deltathetab[0]=0.5;
                shower_deltarb[0]=200;
                shower_nfilmb[0]=1;
            }
            // All other BTs:
            if (ii>0) {
                // its correct like this, since this is the way it is done in
                // the official FJ-Algorithm:
                shower_nfilmb[ii]=diff_pid;
                if (gEDBDEBUGLEVEL>3) cout << "--- ---Inseg->PID() seg->PID() ii diif_pid shower_nfilmb[ii]  " << Inseg->PID()<< "   "  <<  seg->PID() << "   " << ii<< "  " << diff_pid<<"  "<< shower_nfilmb[ii]<<"  " << endl;

                if (diff_pid >= 15 ) shower_sizeb15++;
                if (diff_pid >= 20 ) shower_sizeb20++;
                if (diff_pid >= 30 ) shower_sizeb30++;

                // PUT HERE:  calculation routine for shower_deltasigmathetab
                // see referenc in thesis of Luillo Esposito, page 109.
                shower_deltasigmathetab[diff_pid]=shower_deltasigmathetab[diff_pid]+(TMath::Power(shower_txb[ii]-shower_txb[0],2)+TMath::Power(shower_tyb[ii]-shower_tyb[0],2));

                // PUT HERE:  calculation routine for shower_deltathetab, shower_deltarb
                // ExSetTreebranchNametrapolate the BT [ii] to the position [jj] and then calc the
                // position and slope differences;
                // For the backward extrapolation of the   shower_deltathetab and shower_deltarb
                // calulation for BaseTrack(ii), Basetrack(jj)->Z() hast to be smaller.
                min_shower_deltathetab=99999;   // Reset
                min_shower_deltar=99999;        // Reset

                for (int jj=0; jj<shower_sizeb; jj++)  {
                    if (ii==jj) continue;

                    // since we do not know if BTs are ordered by their Z positions:
                    // and cannot cut directly on the number in the shower entry:
                    if (shower_zb[ii]<shower_zb[jj]) continue;



                    extrapo_diffz=shower_zb[ii]-shower_zb[jj];
                    if (TMath::Abs(extrapo_diffz)>4*1273+1.0) continue;
                    if (TMath::Abs(extrapo_diffz)<1.0) continue; // remove same positions.

                    extrapol_x=shower_xb[ii]-shower_txb[ii]*extrapo_diffz; // minus, because its ii after jj.
                    extrapol_y=shower_yb[ii]-shower_tyb[ii]*extrapo_diffz; // minus, because its ii after jj.

                    // Delta radius we need to extrapolate.
                    test_shower_deltax=extrapol_x;//shower_txb[ii]*(shower_zb[ii]-shower_zb[jj])+shower_xb[ii];
                    test_shower_deltay=extrapol_y;//shower_tyb[ii]*(shower_zb[ii]-shower_zb[jj])+shower_yb[ii];
                    test_shower_deltax=test_shower_deltax-shower_xb[jj];
                    test_shower_deltay=test_shower_deltay-shower_yb[jj];
                    test_shower_deltar=TMath::Sqrt(test_shower_deltax*test_shower_deltax+test_shower_deltay*test_shower_deltay);

                    // Delta theta we do not need to extrapolate. (old version...)
                    //test_shower_deltathetab=TMath::Sqrt(shower_txb[ii]*shower_txb[ii]+shower_tyb[ii]*shower_tyb[ii]);
                    //test_shower_deltathetab=test_shower_deltathetab-TMath::Sqrt(shower_txb[jj]*shower_txb[jj]+shower_tyb[jj]*shower_tyb[jj]);
                    //test_shower_deltathetab=TMath::Abs(test_shower_deltathetab);
                    //----
                    // As before in ShowRec this way of calculation is not equivalent as calculating
                    // DeltaTheta domponentwise:
                    // Code from libShower:
                    // delta = sqrt((SX0-a->GetTXb(l2))*(SX0-a->GetTXb(l2))+((SY0-a->GetTYb(l2))*(SY0-a->GetTYb(l2))));
                    test_shower_deltathetab=TMath::Sqrt(TMath::Power(shower_txb[ii]-shower_txb[jj],2)+TMath::Power(shower_tyb[ii]-shower_tyb[jj],2));

                    // Check if both dr,dt match parameter criteria and then just take these values.....
                    // Maybe a change is necessary because it is not exactly the same as in the off. algorithm.
                    //if (test_shower_deltar<400 && test_shower_deltathetab<0.8 ) {
                    if (test_shower_deltar<150 && test_shower_deltathetab<0.15 ) {
                        if (test_shower_deltar<min_shower_deltar) {
                            min_shower_deltathetab=test_shower_deltathetab;
                            min_shower_deltar=test_shower_deltar;
                            shower_deltathetab[ii]=min_shower_deltathetab;
                            shower_deltarb[ii]=min_shower_deltar;
                        }
                    }
                }
            }
            //-------------------------------------

            shower_purityb=shower_sizebMCNHELP/shower_sizebNHELP;
            shower_trackdensb=trackdensb;
        }
        //-------------------------------------
        for (int i=1; i<59; ++i) {
            shower_deltasigmathetab[i]=TMath::Sqrt(shower_deltasigmathetab[i]);
        }
        shower_numberofilms=shower_nfilmb[shower_sizeb-1]; // there we assume (this is correct always?) when
        // the last shower BT is in the last film...(otherwise we would again have to loop on sizeb array);


        // Fill Tree:
        interim_treesaveb->Fill();
        if (gEDBDEBUGLEVEL>2) interim_treesaveb->Show(interim_treesaveb->GetEntries()-1);

        if (gEDBDEBUGLEVEL>2) cout << "Now we have interim_treesaveb  Entries: "  << interim_treesaveb->GetEntries() << endl;
        if (gEDBDEBUGLEVEL>3) interim_treesaveb->Show(interim_treesaveb->GetEntries()-1);


    } // LOOP OVER ENTRIES FROM showrecoarray: END.
    //-----------------------------------------------------------------------------------------------



    ///------------  ONLY FOR EDBSHOWERREC.... NOT FOR MY OWN STUFF!!---------------------------
    //TFile* interim_fileout = new TFile("shower1.root","RECREATE");
    interim_fileout->cd();
    interim_treesaveb->Write();
    interim_fileout->Close();
    //delete interim_fileout;
    //delete interim_treesaveb;
    ///------------  ONLY FOR EDBSHOWERREC.... NOT FOR MY OWN STUFF!!---------------------------


    return;
}

//-------------------------------------------------------------------------------------------




void EdbShowerRec::TransferTreebranchShowerTreeIntoShowerObjectArray(TTree* treebranch) {
    Log(2, "EdbShowerRec", "--- void* TransferTreebranchShowerTreeIntoShowerObjectArray() ---Interim SOLUTIION !!!! ");

    // Decaring ObjArray which are needed. Use out of private member variables of this class.
    TObjArray* showarray;
    showarray = new TObjArray();
    TTree* eShowerTree = treebranch;

    // EdbShowerP* show=0; // libShowRec
    EdbTrackP* show=0; // libShower

    EdbSegP*    showseg=0;

    //-     VARIABLES: shower_  "treebranch"  reconstruction
    Int_t shower_number_eventb, shower_sizeb, shower_isizeb,shower_showerID;
    Int_t shower_sizeb15, shower_sizeb20, shower_sizeb30;
    Float_t shower_energy_shot_particle;
    Float_t shower_xb[10000];
    Float_t shower_yb[10000];
    Float_t shower_zb[10000];
    Float_t shower_txb[10000];
    Float_t shower_tyb[10000];
    Float_t shower_deltarb[10000];
    Float_t shower_deltathetab[10000];
    Float_t shower_deltaxb[10000];
    Float_t shower_deltayb[10000];
    Int_t   shower_nfilmb[10000];
    Float_t shower_chi2btkb[10000];
    Int_t shower_ntrace1simub[10000]; // MCEvt
    Int_t shower_ntrace2simub[10000]; // s->W()
    Float_t shower_ntrace3simub[10000]; // s->P()
    Int_t shower_ntrace4simub[10000]; // s->Flag()
    Float_t shower_tagprimary[10000];
    Int_t   shower_idb[10000];
    Int_t   shower_plateb[10000];
    Float_t shower_deltasigmathetab[58];
    Int_t   shower_numberofilms;
    Float_t shower_purb; // purity of shower
    Int_t shower_mcDigitIndexTop[1000]; //
    Int_t shower_mcDigitIndexBottom[1000]; //


    // Set Addresses of treebranch tree:
    eShowerTree->SetBranchAddress("number_eventb",&shower_number_eventb);
    eShowerTree->SetBranchAddress("sizeb",&shower_sizeb);
    eShowerTree->SetBranchAddress("sizeb15",&shower_sizeb15);
    eShowerTree->SetBranchAddress("sizeb20",&shower_sizeb20);
    eShowerTree->SetBranchAddress("sizeb30",&shower_sizeb30);
    eShowerTree->SetBranchAddress("isizeb",&shower_isizeb);
    eShowerTree->SetBranchAddress("xb",shower_xb);
    eShowerTree->SetBranchAddress("yb",shower_yb);
    eShowerTree->SetBranchAddress("zb",shower_zb);
    eShowerTree->SetBranchAddress("txb",shower_txb);
    eShowerTree->SetBranchAddress("tyb",shower_tyb);
    eShowerTree->SetBranchAddress("nfilmb",shower_nfilmb);
    eShowerTree->SetBranchAddress("ntrace1simub",shower_ntrace1simub);  // s.eMCEvt
    eShowerTree->SetBranchAddress("ntrace2simub",shower_ntrace2simub); // s.eW
    eShowerTree->SetBranchAddress("ntrace3simub",shower_ntrace3simub); // s.eP
    eShowerTree->SetBranchAddress("ntrace4simub",shower_ntrace4simub); // s.eFlag
    eShowerTree->SetBranchAddress("chi2btkb",shower_chi2btkb);
    eShowerTree->SetBranchAddress("deltarb",shower_deltarb);
    eShowerTree->SetBranchAddress("deltathetab",shower_deltathetab);
    eShowerTree->SetBranchAddress("deltaxb",shower_deltaxb);
    eShowerTree->SetBranchAddress("deltayb",shower_deltayb);
    eShowerTree->SetBranchAddress("tagprimary",shower_tagprimary);
    eShowerTree->SetBranchAddress("energy_shot_particle",&shower_energy_shot_particle);
    eShowerTree->SetBranchAddress("E_MC",&shower_energy_shot_particle);
    eShowerTree->SetBranchAddress("showerID",&shower_showerID);
    eShowerTree->SetBranchAddress("idb",shower_idb);
    eShowerTree->SetBranchAddress("plateb",shower_plateb);
    eShowerTree->SetBranchAddress("deltasigmathetab",shower_deltasigmathetab);
    eShowerTree->SetBranchAddress("lenghtfilmb",&shower_numberofilms);
    eShowerTree->SetBranchAddress("purityb",&shower_purb); // shower purity
    eShowerTree->SetBranchAddress("mcDigitIndexTop",shower_mcDigitIndexTop);
    eShowerTree->SetBranchAddress("mcDigitIndexBottom",shower_mcDigitIndexBottom);


    Int_t nent=eShowerTree->GetEntries();
    cout << "nent = " << nent << endl;

    Int_t nenttotalcount=0;
    cout << "nenttotalcount = " << nenttotalcount << endl;

    Bool_t QUICK_READ=kFALSE;
// 	QUICK_READ=kTRUE; cout << "QUICK_READ=kTRUE; "<< endl;
    if (QUICK_READ) nent=1; // read only first shower... (mainly for debug purposes).

    // Loop over Tree Entries (==different showers):
    for (int i=0; i<nent; ++i) {
        eShowerTree->GetEntry(i);
        //cout << "i = " << i << endl;
        //cout << "shower_tagprimary[0] = " << shower_tagprimary[0] << endl;
        if (gEDBDEBUGLEVEL>3) eShowerTree->Show(i);


        // Since we do not take showers with below 1segment
        //if (shower_sizeb<=1) continue;

//      show = new EdbShowerP(shower_sizeb); // libShowRec
        show = new EdbTrackP(shower_sizeb); // libShower

        // For Frederics algorithm, shower_number_eventb gives NOT -999 as for a BG track,
        // but "0". So we need to disticnct also for "0" cases:
        show ->SetMC(shower_number_eventb,shower_number_eventb);
        if (shower_number_eventb==0) show ->SetMC(-999,999);

        // To Have unique profile histograms in the memory:
        // show = new EdbShowerP(0,eRecoShowerArrayN,shower_sizeb);

        //cout << "------------------- shower_tagprimary[0]  shower_tagprimary[shower_sizeb-1]  " << shower_tagprimary[0] << "  " <<  shower_tagprimary[shower_sizeb-1] << endl;

        int nplmax=0;

        // Loop over shower Entries (==different basetracks):
        for (int j=0; j<shower_sizeb; ++j) {

            //cout << "------------------ --- j= " << j << endl;

            // Keep Notice of the order!!!
            // shower_tagprimary  has to be added first!!
            // temporary solution: if that is not the case as below, no segment will be added....
            // ...
            if (shower_tagprimary[0]!=1) continue;

            // int id, float x, float y, float tx, float ty, float w = 0, int flag = 0
            showseg = new EdbSegP(j,shower_xb[j],shower_yb[j],shower_txb[j],shower_tyb[j],shower_ntrace2simub[j],0);

            showseg->SetZ(shower_zb[j]);
            showseg->SetChi2(shower_chi2btkb[j]);
            showseg->SetP(shower_ntrace3simub[j]);
            showseg->SetFlag(shower_ntrace4simub[j]);
            showseg->SetMC(shower_ntrace1simub[j],shower_ntrace1simub[j]);
            // For Frederics algorithm, ntrace1simub gives NOT -999 as for a BG track,
            // but "0". So we need to disticnct also for "0" cases:
            if (shower_ntrace1simub[j]==0) showseg->SetMC(-999,-999);

            showseg->SetPID(shower_plateb[j]);
            // in real data plateb is 57 for emulsion 57...(sim may be different, depending on orfeo...)

            /// DEBUGTEST
            /// DEBUG TEST FOR FREDERIC !!! (FWM, 27.08.2010)
            showseg->SetPID(shower_nfilmb[j]);
            /// showseg->SetPID(shower_nfilmb[j]);
            // in real data plateb is 57 for emulsion 57...(sim may be different, depending on orfeo...)

            showseg->SetID(shower_idb[j]);
            showseg->SetProb(1);

            //Set Npl right. In libShowRec we dont have to to this since we call the EdbShowerP::Update() function
            //int npl = TMath::Abs(shower_plateb[j]-shower_plateb[0]);
            // The first plate of shower_plateb[0] when used in ::rec() mode is somehow wrong,
            // therefore we rather take nfilmb for plate references...
            int npl = TMath::Abs(shower_nfilmb[j]-shower_nfilmb[0]);
            if (npl>nplmax) {
                show->SetNpl(npl);
            }

            //showseg->PrintNice();

            show->AddSegment(showseg);
            // Attention: If you have a EdbTrackP object and you add a segment, then
            // you have to manually set Counters for N(), NPl(). It is not done auto.
            show -> SetCounters();
        } // Loop over shower Entries (==different basetracks):


        // This ID should be unique!
        // (Means every shower has one different!)
        show->SetID(nenttotalcount);
        nenttotalcount++;

        // Add Shower Object to Shower Reco Array.
        // Not, if its empty:
        // Not, if its containing only one BT:
        if (show->N()>=1) {
            showarray->Add(show);
            ++eRecoShowerArrayN;
            if (gEDBDEBUGLEVEL>3) cout << "EdbShowRec::Treebranch_To_RecoShowerArray   showarray->Add(show) ...   done." << endl;
        }

        if (gEDBDEBUGLEVEL>3) show->PrintNice();

    }// Loop over Tree Entries (==different showers):
    //---------------------------------------------------

    if (gEDBDEBUGLEVEL>2) eShowerTree->Show(0);
    if (gEDBDEBUGLEVEL>2) eShowerTree->Show(nent-1);

    cout << "EdbShowRec::Treebranch_To_RecoShowerArray   Loop over Tree Entries (==different showers)  done." << endl;
    cout << "EdbShowRec::Treebranch_To_RecoShowerArray   showarray->N()" << showarray->GetEntries() <<endl;

    SetRecoShowerArray(showarray);
    SetRecoShowerArrayN(showarray->GetEntries());

    Log(2,"EdbShowRec::Treebranch_To_RecoShowerArray","EdbShowRec::Treebranch_To_RecoShowerArray...done.");
    return;
}






//------------------------------------------------------------------------------------------------------



void EdbShowerRec::ResetAlgoParameters()
{
    eAlgoParameterConeRadius=800;      // Radius of spanning cone from first BT
    eAlgoParameterConeAngle=0.015;       // Opening angle of spanning cone from first BT
    eAlgoParameterConnectionDR=150;      // Connection Criterium: delta R
    eAlgoParameterConnectionDT=0.15;      // Connection Criterium: delta T
    eAlgoParameterNPropagation=3;      // N plates backpropagating
    return;
}

//------------------------------------------------------------------------------------------------------

void EdbShowerRec::SetAlgoParameter(Double_t paravalue, Int_t paranr)
{
    if (paranr==0) eAlgoParameterConeRadius=paravalue;
    else if (paranr==1) eAlgoParameterConeAngle=paravalue;
    else if (paranr==2) eAlgoParameterConnectionDR=paravalue;
    else if (paranr==3) eAlgoParameterConnectionDT=paravalue;
    else if (paranr==4) eAlgoParameterNPropagation=paravalue;
    else {
        return;
    }
}

//------------------------------------------------------------------------------------------------------

Double_t EdbShowerRec::GetAlgoParameter(Int_t paranr)
{
    if (paranr==0) return eAlgoParameterConeRadius;
    else if (paranr==1) return eAlgoParameterConeAngle;
    else if (paranr==2) return eAlgoParameterConnectionDR;
    else if (paranr==3) return eAlgoParameterConnectionDT;
    else if (paranr==4) return eAlgoParameterNPropagation;
    else {
        return 0;
    }
}

//------------------------------------------------------------------------------------------------------

void EdbShowerRec::PrintParameters() {
    cout << "void            EdbShowerRec::PrintParameters() " << endl;
    cout << "eAlgoParameterConeRadius = " << eAlgoParameterConeRadius << endl;
    cout << "eAlgoParameterConeAngle = " << eAlgoParameterConeAngle << endl;
    cout << "eAlgoParameterConnectionDR = " << eAlgoParameterConnectionDR << endl;
    cout << "eAlgoParameterConnectionDT = " << eAlgoParameterConnectionDT << endl;
    cout << "eAlgoParameterNPropagation = " << eAlgoParameterNPropagation << endl;
    return;
}

//------------------------------------------------------------------------------------------------------

void EdbShowerRec::CalcTrackDensity(EdbPattern* pat_interim,Int_t pat_interim_halfsize,Int_t& npat_int,Int_t& npat_total,Int_t& npatN)
{
    if (gEDBDEBUGLEVEL>3) cout << "-------------void EdbShowerRec::CalcTrackDensity(&pat_interim,pat_interim_halfsize,&npat_int,&npat_total)"<<endl;
    npat_int=pat_interim->GetN();
    if (npat_int<=0) return;
    npat_total+=npat_int;
    ++npatN;
    if (npatN>0) trackdensb=(Float_t)npat_total/(Float_t)npatN/pat_interim_halfsize/pat_interim_halfsize/4.0*1000.0*1000.0; // BT/mm2
    if (gEDBDEBUGLEVEL>3) {
        cout << "npat_int = " << npat_int << endl;
        cout << "npat_total = " << npat_total << endl;
        cout << "npatN = " << npatN << endl;
        cout << "shower_trackdensb = " << trackdensb << endl;
    }
    return;
}


//------------------------------------------------------------------------------------------------------

Double_t EdbShowerRec::InvariantMass(EdbSegP* s1, EdbSegP* s2, Double_t Momentum1,Double_t Momentum2, Double_t ID1, Double_t ID2)
{
    cout << "EdbShowerRec::InvariantMass()" << endl;
    cout << "s1->P()  " << s1->P() << endl;
    cout << "s2->P()  " << s2->P() << endl;

    // If no momentum is specified, then it will be taken from segments momentum:
    if (Momentum1<0) Momentum1=s1->P();
    if (Momentum2<0) Momentum2=s2->P();

    //mass reconstruction for Pi0 using 2 tracks
    //you have to know the momentum of the 2 tracks and their slopes (Tx,Ty)
    //first implementation: frederic juget
    TLorentzVector v1;
    TLorentzVector v2;
    Double_t energy1=0;
    Double_t energy2=0;
    v1.SetPxPyPzE(s1->TX(),s1->TY(),1,Momentum1);
    v2.SetPxPyPzE(s2->TX(),s2->TY(),1,Momentum2);
    TLorentzVector vsum=v1+v1;
    Double_t invM=vsum.Mag2();
    cout << v1.Mag() << endl;
    cout << v2.Mag() << endl;
    cout << vsum.Mag() << endl;
    return invM;
}


//------------------------------------------------------------------------------------------------------


void EdbShowerRec::Help()
{
    // Help Function
    cout << "--------------------------------------------------------------------------" << endl;
    cout << "--                                                                      --" << endl;
    cout << "-- EdbShowerRec                                                         --" << endl;
    cout << "--                                                                      --" << endl;
    cout << "-- Base Class for the shower reconstruction in OPERA brick.             --" << endl;
    cout << "-- This class has been initially developed by F. Juget, G. Lutter       --" << endl;
    cout << "-- and is ongoingly modified and developed by Frank Meisel, which is    --" << endl;
    cout << "-- currently the corresponding author: frank.meisel@gmx.de              --" << endl;
    cout << "--                                                                      --" << endl;
    cout << "-- The main goal of this library is to make the shower search in the    --" << endl;
    cout << "-- scanned emulsion data as easy as possible, and as most automatically.--" << endl;
    cout << "-- There exist also two manuals, to be found under 'documentation'. In  --" << endl;
    cout << "-- these some general handling of the package is explained.             --" << endl;
    cout << "--                                                                      --" << endl;
    cout << "-- The EdbShowerRec class makes following things:                       --" << endl;
    cout << "-- Starting from the emulsion scanned data (represented as an EdbPVRec  --" << endl;
    cout << "-- object in FEDRA), it does search and reconstruction of showers,      --" << endl;
    cout << "-- based on different algorithms, which the user can set indivially     --" << endl;
    cout << "-- (normaly not necessary). After that, we try to identify              --" << endl;
    cout << "-- characteristics of the shower, as it is: energy and particle ID.     --" << endl;
    cout << "--                                                                      --" << endl;
    cout << "--------------------------------------------------------------------------" << endl;
    return;
}


//------------------------------------------------------------------------------------------------------


void EdbShowerRec::Print() {
    // Print Function
    cout << "-- EdbShowerRec   Print                                                     --" << endl;
}


//------------------------------------------------------------------------------------------------------

void EdbShowerRec::SetInBTArray( EdbPVRec* Ali, Int_t mode ) {

    // SetInBTArray( EdbPVRec* Ali , Int_t mode )
    // Set the Initiator BaseTrack Array From the EdbPVRec
    // object. Using a mode of:
    // -1) if it has tracks inside, we take first segment of the
    // -1) tracks              (TO BE IMPLEMENTED...mode=-1)
    //  0) otherwise (mode is 0) ALL Basetracks in the volume.
    //  0) ATTENTION...TAKES  VERY  LONG
    // 	1) If  mode  is 1, then we take only the
    // 	1) first Z plate            ... Mainly DEBUG purpose ....
    // 	2) If  mode  is 2, then we take only 5 random
    // 	2) basetracksfirst Z plate  ... Mainly DEBUG purpose ....

    if (gEDBDEBUGLEVEL>1) {
        cout << "-----------------------------------------------------------------------------" << endl;
        cout << "-- EdbShowerRec   SetInBTArray( EdbPVRec* Ali , Int_t mode )               --" << endl;
        cout << "-- EdbShowerRec   Set the Initiator BaseTrack Array From the EdbPVRec      --" << endl;
        cout << "-- EdbShowerRec   object.                                                  --" << endl;
        cout << "-- EdbShowerRec   A) if it has tracks inside, we take first segment of the --" << endl;
        cout << "-- EdbShowerRec   A) tracks              (TO BE IMPLEMENTED...mode=-1)     --" << endl;
        cout << "-- EdbShowerRec   B) otherwise (mode is 0) ALL Basetracks in the volume.   --" << endl;
        cout << "-- EdbShowerRec   B) ATTENTION...TAKES  VERY  LONG                         --" << endl;
        cout << "-- EdbShowerRec   C) If  mode  is 1, then we take only the                 --" << endl;
        cout << "-- EdbShowerRec   C) first Z plate            ... Mainly DEBUG purpose ....--" << endl;
        cout << "-- EdbShowerRec   D) If  mode  is 2, then we take only (1,5) random            --" << endl;
        cout << "-- EdbShowerRec   D) basetracksfirst Z plate  ... Mainly DEBUG purpose ....--" << endl;
    }

    TObjArray* InBTArray = new TObjArray();

    Int_t nBT=1;
    if (mode>10)  {
        // help construction, i hope noone ever needs to look at this ugly part of code :)
        nBT=mode;
        mode=2;
    }

    if (1==mode) {
        InBTArray->Clear();
        EdbPattern* pat = Ali->GetPatternZLowestHighest();
        int nseg= pat->N();
        for (int j=0; j<nseg; j++) {
            EdbSegP* seg = pat->GetSegment(j);
            InBTArray->Add(seg);
        }
    }
    else if (2==mode) {
        for (int k=0; k<nBT; k++) {
            EdbPattern* pat = Ali->GetPattern(gRandom->Uniform(0,Ali->Npatterns()));
            EdbSegP* seg = pat->GetSegment(gRandom->Uniform(0,pat->N()));
            //seg->PrintNice();
            InBTArray->Add(seg);
        }
        cout << "-- EdbShowerRec   D) Set " << nBT << " Random BTs for shower reco. --" << endl;
    }
    else {
        int npat=Ali->Npatterns();
        for (int i=0; i<npat; i++) {
            EdbPattern* pat = Ali->GetPattern(i);
            int nseg= pat->N();
            for (int j=0; j<nseg; j++) {
                EdbSegP* seg = pat->GetSegment(j);
                InBTArray->Add(seg);
            }
        }
    }


    if (gEDBDEBUGLEVEL>1) cout << "-- EdbShowerRec   InBTArray->GetEntries(): " << setw(20) << InBTArray->GetEntries()  << "       --" << endl;
    eInBTArray = InBTArray;
    eInBTArrayN=eInBTArray->GetEntries();
    if (eInBTArrayN>0) {
        eInBTArrayIsFilled=kTRUE;
    }

    if (gEDBDEBUGLEVEL>1) cout << "-- EdbShowerRec   SetInBTArray( EdbPVRec* Ali )                       done.--" << endl;
    if (gEDBDEBUGLEVEL>1) cout << "-----------------------------------------------------------------------------" << endl;
    return;
}

