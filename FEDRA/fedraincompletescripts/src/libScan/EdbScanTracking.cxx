  //-- Author :  Valeri Tioukov   19/02/2011
  //////////////////////////////////////////////////////////////////////////
  //                                                                      //
  // EdbScanTracking                                                      //
  //                                                                      //
  // To handle tracking in the scanset                                    //
  //                                                                      //
  //////////////////////////////////////////////////////////////////////////
#include "EdbLog.h"
#include "EdbDataSet.h"
#include "EdbScanTracking.h"
#include "EdbAlignmentV.h"
#include "EdbPlateAlignment.h"
#include "EdbTrackFitter.h"
#include "EdbPhys.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TROOT.h"
  
using namespace TMath;
  
ClassImp(EdbTrackAssembler)
ClassImp(EdbScanTracking)
  
//--------------------------------------------------------------------------------------
EdbTrackAssembler::~EdbTrackAssembler()
{
  SafeDelete(eHistNcnd);
  SafeDelete(eHistProbAll);
  SafeDelete(eHistProbBest);
  SafeDelete(eHistThetaBest);
  SafeDelete(eHistThetaAll);
}
  
  //--------------------------------------------------------------------------------------
  EdbTrackAssembler::EdbTrackAssembler()
{
  eMapMarg = 50.; // [microns]
  eZ = 0;
  eCellN=10;    //mean n/cell
  eDTmax=0.07;
  eDRmax=45.;
  eDZGapMax = 5000;
  eProbMin = 0.001;
  eCollisionsRate=0;
  
  eHistNcnd     = new TH1F("Ncnd","number of candidates after preliminary selection", 20,0.5,20.5);
  eHistProbBest = new TH1F("ProbBest","prob for best selected candidate", 250,0,1);
  eHistProbAll  = new TH1F("ProbAll","prob for all candidates", 250,0,1);
  eHistThetaBest = new TH1F("ThetaBest","angle theta for best selected candidate", 180,0,TMath::PiOver2());
  eHistThetaAll  = new TH1F("ThetaAll","angle theta for all candidates", 180,0,TMath::PiOver2());

    // for basetracks:
  eCond.SetDefault();
  eCond.SetSigma0( 4, 4, 0.005, 0.005 );
  eCond.SetPulsRamp0(14., 21.);
  eCond.SetPulsRamp04(14., 21.);
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::DoubletsFilterOut(EdbPattern &p)
{
  EdbAlignmentV adup;
  adup.eDVsame[0]=adup.eDVsame[1]=10.;
  adup.eDVsame[2]=adup.eDVsame[3]=0.08;
  adup.FillGuessCell(p,p,1.);
  adup.FillCombinations();
  adup.DoubletsFilterOut(0);   // assign flag -10 to the duplicated segments
}
  
  //--------------------------------------------------------------------------------------
void EdbTrackAssembler::CheckPatternAlignment(EdbPattern &p, EdbPlateP &plate,  int nsegmin)
{
  Log(0,"EdbTrackAssembler::CheckPatternAlignment","");
  
  ExtrapolateTracksToZ( p.Z(), nsegmin);
  int ntr = eTrZ.GetEntriesFast();
  EdbPattern ptr( 0, 0, p.Z(), ntr ); 
  for(int i=0; i<ntr; i++) ptr.AddSegment( *((EdbSegP*)(eTrZ.UncheckedAt(i))) );
    
  EdbPlateAlignment al;
  al.SetSigma( 25, 0.015 );
  al.eOffsetMax = 500.;
  al.eDZ        = 0;
  al.eDPHI      = 0.00;
    //al.eDoCoarse=1;
  al.Align(ptr,p,0);
  
  EdbAffine2D *aff = al.eCorrL[0].GetAffineXY();
  aff->Invert();
  aff->Print();
  p.Transform(aff);
  plate.GetAffineXY()->Transform(aff);
  
  EdbAffine2D *afftxty = al.eCorrL[0].GetAffineTXTY();
  afftxty->Invert();
  afftxty->Print();
  p.TransformA(afftxty);
  plate.GetAffineTXTY()->Transform(afftxty);
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::AddPattern(EdbPattern &p)
{
    //int ntrBefore=eTracks.GetEntriesFast();
    //if(ntrBefore>0) ExtrapolateTracksToZ(p.Z());
    
    //DoubletsFilterOut();
  int nseg = p.N();
  Log(3,"EdbTrackAssembler::AddPattern","try to add %d segments",p.N());
  int attached=0;
  for(int j=0; j<nseg; j++) {
    EdbSegP *s = p.GetSegment(j);
    if(s->Flag()==-10)   continue;
    s->SetErrors();
    eCond.FillErrorsCov(s->TX(),s->TY(),s->COV());
    if( !AddSegment( *(p.GetSegment(j)) ) )
      AddSegmentAsTrack( *(p.GetSegment(j)) );
    else attached++;
  }
  int ntrAfter = eTracks.GetEntriesFast();
    //int totSegTr=0;
    //for(int i=0; i<ntrAfter; i++) totSegTr += ((EdbTrackP*)(eTracks.At(i)))->N();
  Log(2,"EdbTrackAssembler::AddPattern","with z=%10.2f   %d/%d attached/tried; total collisions: %d;  tracks: %d",
      p.Z(), attached, nseg, eCollisionsRate, ntrAfter );
}
  
//--------------------------------------------------------------------------------------
EdbTrackP *EdbTrackAssembler::AddSegment(EdbSegP &s)
{
  TObjArray trsel;
  float v[2] = { s.X(), s.Y() };
  int nsel =  eTrZMap.SelectObjectsC( v, eDRmax+50 , trsel );
  Log(3,"EdbTrackAssembler::AddSegment", "nsel = %d",nsel);
  if(!nsel) { 
    return 0;  }
    float prob, probmax = eProbMin;
    EdbSegP *ssbest = 0;
    int ncnd = 0;
    for(int i=0; i<nsel; i++) {
      EdbSegP *ss = (EdbSegP*)(trsel.At(i));
      prob = ProbSeg( *ss, s );
      Log(3,"EdbTrackAssembler::AddSegment", "prob(probmin) = %f (%f) ",prob, eProbMin);
      if(prob<eProbMin) continue;
      ncnd++;
      if(eHistProbAll) eHistProbAll->Fill(prob);
      if(eHistThetaAll) eHistThetaAll->Fill(ATan(ss->Theta()));
      if( prob > probmax ) { ssbest = ss; probmax=prob; }
    }
    if(!ssbest)  return 0;
    s.SetProb(probmax);
    if(eHistNcnd)     eHistNcnd->Fill(ncnd);
    if(eHistProbBest) eHistProbBest->Fill(probmax);
    if(eHistThetaBest) eHistThetaBest->Fill(ATan(ssbest->Theta()));
    EdbTrackP *t = (EdbTrackP*)(ssbest);
    EdbSegP *sz = t->GetSegmentWithClosestZ( t->Z(), 45. );
    if(!sz) t->AddSegment( eSegments.AddSegment(s) );
    else {
      if( !SameSegment(s,*sz) ) {
        if( s.Prob() > sz->Prob() )  t->SubstituteSegment( sz ,  eSegments.AddSegment(s) );
        eCollisionsRate++;
      }
    }
    return t;
}

  //--------------------------------------------------------------------------------------
  bool EdbTrackAssembler::SameSegment( EdbSegP &s1, EdbSegP &s2 )
{
  if( Abs( s1.X() - s2.X() )   <0.000001 &&
      Abs( s1.Y() - s2.Y() )   <0.000001 &&
      Abs( s1.TX()- s2.TX() )  <0.000001 &&
      Abs( s1.TY()- s2.TY() )  <0.000001 &&
      Abs( s1.W() - s2.W() )   <0.000001  )    return true;
  return false;
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::RecalculateSegmentsProb( EdbTrackP &tr )
{
    // assumed that track is fitted: reset the segmets probabilities
  int n=tr.N();
  for(int i=0; i<n; i++) 
    tr.GetSegment(i)->SetProb( ProbSeg( *(tr.GetSegmentF(i)),  *(tr.GetSegment(i)) ) ); 
}
  
  //--------------------------------------------------------------------------------------
  float EdbTrackAssembler::ProbSeg( EdbSegP &s1, EdbSegP &s2 )
{
  // return the probability that the second segment can belong to track defined by s1
  float dtx = s1.TX() - s2.TX();
  if( Abs( dtx ) > eDTmax )    return 0;
  float dty = s1.TY() - s2.TY();
  if( Abs( dty ) > eDTmax )    return 0;
  double dt2 = dtx*dtx +  dty*dty;
  if(dt2>eDTmax*eDTmax)        return 0;
  
  float dz = s2.Z()-s1.Z();
  float dx = s2.X() - (s1.X() + dz*s1.TX());
  if( Abs( dx ) > eDRmax )     return 0;
  float dy = s2.Y() - (s1.Y() + dz*s1.TY());
  if( Abs( dy ) > eDRmax )     return 0;
  double dr2 = dx*dx +  dy*dy;
  if(dr2>eDRmax*eDRmax)        return 0;
  
  float prob=0;
  if(eDoUseMCS==1){
    ///workaround to get previous(not propagated) segment
    EdbTrackP* t = dynamic_cast<EdbTrackP*> (&s1);
    EdbSegP* seg = t?(const_cast<EdbSegP*>(t->TrackEnd())):0;
    prob = seg?(eFitter.ProbSegMCS(seg, &s2)):0;
  }
  else if(eDoUseMCS==2)
  {
    EdbTrackP* t = dynamic_cast<EdbTrackP*> (&s1);
    EdbSegP* seg = t?(const_cast<EdbSegP*>(t->TrackEnd())):0;
    if(seg) {
      float chi = eFitter.Chi2Seg( seg, &s2 );
      prob = (float)TMath::Prob( chi*chi, 4);
    }
    if(prob<0.001&&gEDBDEBUGLEVEL>2) {
      Log(3,"EdbTrackAssembler::ProbSeg", "dx,dy,dz,dtx,dty: %f %f %f %f %f  prob= %f",dx,dy,dz,dtx,dty, prob);
      seg->Print();
      s2.Print();
      eCond.Print();
    }

   }
   else
   {
     EdbSegP s;
     float chi = eFitter.Chi2SegM( s1, s2, s, eCond, eCond );
     prob = (float)TMath::Prob( chi*chi, 4);
   }
   

  prob *= eCond.ProbSeg( s2.Theta(), s2.W() );            // the proability component depending on the grains number
  prob *= (float)TMath::Prob( s2.Chi2()*s2.Chi2(), 4 );   // the proability component depending on the segment strength

  return prob;
}
  
  //--------------------------------------------------------------------------------------
  EdbTrackP *EdbTrackAssembler::AddSegmentAsTrack(EdbSegP &s)
{
  if(s.W()<16    )  return 0;
  if(s.Chi2()>2.5)  return 0;
  EdbTrackP *t = new EdbTrackP( eSegments.AddSegment(s), 0.139);    // EdbTrackAssembler is owner of segments 
  eTracks.Add(t);
    //EdbSegP ss;
    //t->MakePredictionTo(eZ,ss);
    //eTrZ.AddSegment(ss);
  return t;
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::ExtrapolateTracksToZ(float z, int nsegmin)
{
  eTrZ.Clear();
  eTrZMap.CleanCells();
    
  eZ=z;
  int n=eTracks.GetEntriesFast();
  for(int i=0; i<n; i++) {
    EdbTrackP *t = (EdbTrackP*)(eTracks.At(i));
      
    if( t->N() < nsegmin )      continue;
    if( !AcceptDZGap(*t, z) )   continue;
    t->MakePredictionTo(eZ,*t);                        // extrapolation of tracks itself
      //((EdbSegP *)t)->PrintNice();
    eTrZ.Add(t);
  }
    
    //FillTrZMap();
    //if(gEDBDEBUGLEVEL>2) eTrZMap.PrintStat();
}
  
  //--------------------------------------------------------------------------------------
  bool EdbTrackAssembler::AcceptDZGap(EdbTrackP &t, float z)
{
  float z1 = t.GetSegmentFirst()->Z();
  float z2 = t.GetSegmentLast()->Z();
  if(Min( Abs(z1-z), Abs(z2-z)) > eDZGapMax ) return false;
  return true;
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::FillTrZMap()
{
  int n=eTrZ.GetEntriesFast();
  Log(2,"EdbTrackAssembler::FillTrZMap", "with %d tracks",n);
  for(int i=0; i<n; i++) {
    EdbSegP *s = (EdbSegP*)(eTrZ.At(i));
    eTrZMap.AddObject( s->X(), s->Y(), s );
  }
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::InitTrZMap( const char *str )
{
  int   nx=0, ny=0, ncell=0;
  float xmi,xma, ymi, yma;
  sscanf(str,"%d %f %f %d %f %f %d",&nx,&xmi,&xma,&ny,&ymi,&yma,&ncell);
  InitTrZMap( nx,xmi,xma,ny,ymi,yma,ncell );
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::InitTrZMap( int nx, float xmi, float xma, 
                                      int ny, float ymi, float yma,  int ncell)
{
  float  mi[2] = {  xmi, ymi };
  float  ma[2] = {  xma, yma };
  int     n[2] = { nx, ny };
  eTrZMap.InitCell(ncell, n, mi, ma);
}
  
  //--------------------------------------------------------------------------------------
  void EdbTrackAssembler::InitTrZMap()
{
  /*  float  mi[2] = {  eTrZ.Xmin()-eMapMarg, eTrZ.Ymin()-eMapMarg };
  float  ma[2] = {  eTrZ.Xmax()+eMapMarg, eTrZ.Ymax()+eMapMarg };
  float  dens = eTrZ.N()/( (ma[0]-mi[0])*(ma[1]-mi[1]));
  float  step=10000;
  if(dens>0.00000001) step = Sqrt( eCellN/dens );
  int    n[2] = { int((ma[0]-mi[0])/step)+1, int((ma[1]-mi[1])/step)+1 };
  float stepX = (ma[0]-mi[0])/n[0];
  float stepY = (ma[1]-mi[1])/n[1];
  n[0] = int((ma[0]-mi[0]+1.)/stepX);
  n[1] = int((ma[1]-mi[1]+1.)/stepY);
  eTrZMap.InitCell(3*eCellN, n, mi, ma);*/
}

//--------------------------------------------------------------------------------------
void EdbTrackAssembler::SetSegmentsErrors()
{
  int ntr = eTracks.GetEntriesFast();
    for( int i=0; i<ntr; i++ )     {
      EdbTrackP *t = (EdbTrackP*)(eTracks.At(i));
      if(t->Flag()==-10) continue;
      int nseg=t->N();
      if(nseg>0)  { 
        for(int j=0; j<nseg; j++) {
          EdbSegP *s = t->GetSegment(j);
          s->SetErrors();
          eCond.FillErrorsCov(s->TX(),s->TY(),s->COV());
        }
      }
    }
}

//--------------------------------------------------------------------------------------
void EdbTrackAssembler::FitTracks()
{
  EdbTrackFitter fit;

  int ntr = eTracks.GetEntriesFast();
  for( int i=0; i<ntr; i++ )     {
    EdbTrackP *t = (EdbTrackP*)(eTracks.At(i));
    if(t->Flag()==-10) continue;
    int nseg=t->N();
    t->FitTrackKFS(0,10000);
        //fit.FitTrackLine(*t);
    if(nseg>1) RecalculateSegmentsProb(*t);
  }
}

//______________________________________________________________________________
void EdbTrackAssembler::CombTracks( TObjArray &selected )
{
    // eliminate crossing&overlapping tracks with multiple segments usage
  
  int nsegMin=2;
  int nGapMax=50;
  
  int ntr = eTracks.GetEntriesFast();
  Log(3,"EdbTrackAssembler::CombTracks","Comb %d tracks");
  
    // *** sort tracks by quality
  
  TIndexCell cn;  //"nseg:prob:entry"
  Long_t v[3];
  
  int nsegtot=0;
  EdbTrackP *tr=0;
  for(int i=0; i<ntr; i++) {
    tr = (EdbTrackP*)(eTracks.At(i));
    if( tr->Flag() == -10 ) continue;
    if( tr->N() < nsegMin ) continue;
    tr->SetID(i);
    tr->SetCounters();
    nsegtot += tr->SetSegmentsTrack(-1);
    v[0]= -(tr->N());
    v[1]= (Long_t)((1.-tr->Prob())*100);
    v[2]= i;
    cn.Add(3,v);
  }
  cn.Sort();
  
  Log(3,"EdbTrackAssembler::CombTracks","%d tracks with %d segments for processing...",ntr,nsegtot);
  
    // *** set track ID for segments attached to
  
  TIndexCell *cp=0, *c=0;
  int nn=cn.GetEntriesFast();
  for(int i=nn-1; i>=0; i--) {
    cp = cn.At(i);                              // tracks with fixed npl
    int np = cp->GetEntriesFast();
    for(int ip=np-1; ip>=0; ip--) {
      c = cp->At(ip);                           // tracks with fixed Npl & Prob
      int nt = c->GetEntriesFast();
      for(int it=0; it<nt; it++) {
        tr = (EdbTrackP*)(eTracks.At( c->At(it)->Value() ) );
        tr->SetSegmentsTrack();
      }
    }
  }
  
  
  cp=0;   c=0;
  nn=cn.GetEntriesFast();
  for(int i=0; i<nn; i++) {
    cp = cn.At(i);                              // tracks with fixed npl

    int np = cp->GetEntriesFast();
    for(int ip=0; ip<np; ip++) {
      c = cp->At(ip);                           // tracks with fixed Npl & Prob

      int nt = c->GetEntriesFast();
      for(int it=0; it<nt; it++) {
  
        tr = (EdbTrackP*)(eTracks.At( c->At(it)->Value() ) );

        if(tr->RemoveAliasSegments()>0){
          if(tr->N()<nsegMin)             tr->SetFlag(-10);
          if(tr->CheckMaxGap()>nGapMax)   tr->SetFlag(-10);
        }

        if( tr->Flag() != -10 ) selected.Add(tr);
      }
    }
  }

}

//=======================================================================================
  EdbScanTracking::EdbScanTracking()
{
  eNsegMin=2;
  eNgapMax=50;
}

//--------------------------------------------------------------------------------------
void EdbScanTracking::TrackSetBT(EdbID idset, TEnv &env)
{
  
  // read scanset object
  EdbScanSet *ss = eSproc->ReadScanSet(idset);
  if(!ss) { Log(1,"EdbScanTracking::TrackSetBT",
    "Error! set for %s do not found",idset.AsString()); return; }
  
    int npl = ss->eIDS.GetSize();
    if(npl<2) { Log(1,"EdbScanTracking::TrackSetBT", "Warning! npl<2 : %d stop tracking!",npl); return; }
  
  // create and init tracking object
    EdbTrackAssembler etra;
  
    etra.eCond.SetSigma0(        env.GetValue("fedra.track.Sigma0"         , "3 3 0.005 0.005") );
    etra.eCond.SetPulsRamp0(     env.GetValue("fedra.track.PulsRamp0"      , "15 20") );
    etra.eCond.SetPulsRamp04(    env.GetValue("fedra.track.PulsRamp04"     , "15 20") );
    etra.eCond.SetDegrad(        env.GetValue("fedra.track.Degrad"         , 4) );
    etra.SetRadLength(           env.GetValue("fedra.track.RadX0"          , 5810.) );
    etra.eDoUseMCS              = env.GetValue("fedra.track.do_use_mcs"    , 0 );
      
    etra.eDTmax                 = env.GetValue("fedra.track.DTmax"          ,     0.07 );
    etra.eDRmax                 = env.GetValue("fedra.track.DRmax"          ,    45.   );
    etra.eDZGapMax              = env.GetValue("fedra.track.DZGapMax"       ,  5000.   );
    etra.eProbMin               = env.GetValue("fedra.track.probmin"        ,  0.001   );
    bool        do_erase        = env.GetValue("fedra.track.erase"          ,  false   );
    const char  *cut            = env.GetValue("fedra.readCPcut"            ,     "1"  );
    bool        do_misalign     = env.GetValue("fedra.track.do_misalign"    ,      0   );
    int         npass           = env.GetValue("fedra.track.npass"          ,      1   );
    float       misalign_offset = env.GetValue("fedra.track.misalign_offset",    500.  );
    bool        do_local_corr   = env.GetValue("fedra.track.do_local_corr"  ,      1   );
    int        NcpMin_local_corr= env.GetValue("fedra.track.NcpMin_local_corr"  ,  0   );
    bool        eDoRealign      = env.GetValue("fedra.track.do_realign"     ,      0   );
    bool        do_comb         = env.GetValue("fedra.track.do_comb"        ,      0   );
    eNsegMin                    = env.GetValue("fedra.track.NsegMin"        ,      2   );
    float       momentum        = env.GetValue("fedra.track.momentum"       ,      2.  );  
    etra.SetMomentum (momentum);
    
    etra.InitTrZMap(  env.GetValue("fedra.track.TrZmap", "2400 0 120000   2000 0 100000   30" ) );
  
    //EdbPattern p;
    
    EdbAffine2D misalign[60];
    if(do_misalign) {
        //           1 2 3 4 5 6  7  8  9
      int dx[9] = {0,0,1,1,1,0,-1,-1,-1};
      int dy[9] = {0,1,1,0,-1,-1,-1,0,1};
      for(int i=0; i<60; i++) {
        misalign[i].ShiftX( dx[i%9] * misalign_offset );
        misalign[i].ShiftY( dy[i%9] * misalign_offset );
        if(gEDBDEBUGLEVEL>1) printf("%d |  %d  %d\n",i, dx[i%9], dy[i%9]);
      }
    }
  
  // read segments and use them for tracking
    for(int ipass=0; ipass<npass; ipass++) {
      if(gEDBDEBUGLEVEL>1) printf("\n\n*************** ipass=%d ************\n",ipass);
      etra.eCollisionsRate=0;
      for(int i=0; i<npl; i++) {
        EdbID *id = ss->GetID(i);
      
        EdbPlateP *plate = ss->GetPlate(id->ePlate);
      
        EdbPattern p;
        eSproc->ReadPatCPnopar(p,*id, cut, do_erase);
        p.SetZ(plate->Z());
        p.SetSegmentsZ();
        p.SetID(i);
        p.SetPID(i);
        p.SetSegmentsPID();
      //plate->Print();
        p.Transform(    plate->GetAffineXY()   );
        p.TransformShr( plate->Shr() );
        p.TransformA(   plate->GetAffineTXTY() );
        p.SetSegmentsPlate(id->ePlate);
        Log(1,"EdbScanTracking::TrackSetBT",
    "Starting local correction for plate %i",i); 
        if(do_local_corr && plate->Z()<0.) {
          //combine maps from the various plates
          EdbCorrectionMap corrmap = EdbCorrectionMap();
          corrmap.Init(2,0,120000,2,0,100000); 
          corrmap.ApplyCorrections(plate->Map());
          
           //invert map we want everything to plate 25, and this maps goes from plate whichplate+1 to plate whichplate
          for(int i=0; i<corrmap.Ncell(); i++) {
            EdbLayer *loc  = corrmap.GetLayer(i);
            loc->Invert();
          }

          //combining the map with all the map from the other plates up to the most downstream one
          for (int iplate = 1; iplate< npl - (i+1); iplate++){
  
           EdbCorrectionMap *othermap = new EdbCorrectionMap();

           othermap->Init(2,0,120000,2,0,100000); 
  
           othermap->ApplyCorrections(ss->GetPlate(iplate+i)->Map());
  
           for(int i=0; i<othermap->Ncell(); i++) {
            EdbLayer *loc  = othermap->GetLayer(i);
            loc->Invert();
           }

           corrmap.ApplyCorrections(*othermap); //apply inverse corrections of layer
          }
          int nseg = p.N();
          for(int j=0; j<nseg; j++) {
            EdbSegP *s = p.GetSegment(j);
            //int ncp_maplayer = -1;
            //if (plate->Map().GetLayer(s->X(),s->Y()))
            // ncp_maplayer = plate->Map().GetLayer(s->X(),s->Y())->Ncp();
            //if (ncp_maplayer > NcpMin_local_corr)
            plate->CorrectSegLocal(*s);
            //else 
            //  Log(3, "EdbScanTracking::TrackSetBT","Only %i couples in layer, not applying correction.", ncp_maplayer);          
          }
        }
      
      
        if(do_misalign) {
          p.Transform(&misalign[i]);
          Log(2,"EdbScanTracking::TrackSetBT","apply misalignment of %f",misalign_offset);
        //misalign[i].Print();
        }
  
        if(i>0) etra.ExtrapolateTracksToZ(p.Z());
        if(eDoRealign) 
        {
          if( i==1 ) etra.CheckPatternAlignment(p,*plate,1);
          if( i>1  ) etra.CheckPatternAlignment(p,*plate,2);
        }
        etra.FillTrZMap();
        etra.AddPattern(p);
      }
    }
    
    if(eDoRealign) eSproc->WriteScanSet(idset,*ss);
    
    int ntr = etra.Tracks().GetEntriesFast();
    
    for(int i=0; i<ntr; i++) {
      EdbTrackP *t = (EdbTrackP *)(etra.Tracks().At(i));
      if(t->N()<eNsegMin)  t->SetFlag(-10);
    }
    
    etra.SetSegmentsErrors();
    etra.FitTracks();

    TObjArray selectedTracks(ntr);
    if(do_comb) {
      etra.CombTracks(selectedTracks);
    } else {
      int cnt=0;
      for(int i=0; i<ntr; i++) {
        EdbTrackP *t = (EdbTrackP *)(etra.Tracks().At(i));
        if(t->Flag()!=-10)  {
          t->SetID(cnt++);
          t->SetCounters();
          t->SetSegmentsTrack();
          t->SetP(momentum);
          selectedTracks.Add(t);
        }
      }
    }
    EdbDataProc::MakeTracksTree( selectedTracks, 0., 0., Form("b%s.trk.root", idset.AsString()) );
    TFile f( Form("b%s.trk.root", idset.AsString()) ,"UPDATE");
    env.Write();
    f.Close();
    
    SaveHist(idset,etra);
}

//--------------------------------------------------------------------------------------
void EdbScanTracking::SaveHist(EdbID idset, EdbTrackAssembler &etra)
{
  TFile f( Form("b%s.trk.root", idset.AsString()) ,"UPDATE");
  if(etra.eHistNcnd) etra.eHistNcnd->Write();
  if(etra.eHistProbAll) etra.eHistProbAll->Write();
  if(etra.eHistProbBest) etra.eHistProbBest->Write();
  TH1F *probrest = 0, *probPurity=0;
  if(etra.eHistProbAll&&etra.eHistProbBest) {
    probrest = (TH1F*)(etra.eHistProbAll->Clone("ProbRest"));
    probrest->SetTitle("prob for the other candidates");
    probrest->Add(etra.eHistProbBest,-1);
    probrest->Write();
    probPurity = (TH1F*)(probrest->Clone("ProbPurity"));
    probPurity->SetTitle("Nother/Nall vs prob");
    probPurity->Divide(etra.eHistProbAll);
    probPurity->Write();
  }
  if(etra.eHistThetaAll) etra.eHistThetaAll->Write();
  if(etra.eHistThetaBest) etra.eHistThetaBest->Write();
  TH1F *thetarest = 0, *thetaPurity=0;
  if(etra.eHistThetaAll&&etra.eHistThetaBest) {
    thetarest = (TH1F*)(etra.eHistThetaAll->Clone("ThetaRest"));
    thetarest->SetTitle("theta for other candidates");
    thetarest->Add(etra.eHistThetaBest,-1);
    thetarest->Write();
    thetaPurity = (TH1F*)(thetarest->Clone("ThetaPurity"));
    thetaPurity->SetTitle("Nother/Nall vs theta");
    thetaPurity->Divide(etra.eHistThetaAll);
    thetaPurity->Write();
  }
  
  gStyle->SetPalette(1);
  gStyle->SetOptStat(1);
  bool batch = gROOT->IsBatch();
  gROOT->SetBatch();
  
  TCanvas *c = new TCanvas("purity","tracking purity",900,800);
  c->Divide(2,3);

  c->cd(1)->SetLogy(); 
  etra.eHistNcnd->SetAxisRange(0,10);
  etra.eHistNcnd->Draw();

  c->cd(2); probrest->SetLineColor(kBlue); probrest->Draw();

  c->cd(3)->SetLogy();
  etra.eHistProbAll->Draw();
  etra.eHistProbBest->SetLineColor(kRed); etra.eHistProbBest->Draw("same");
  probrest->SetLineColor(kBlue); probrest->Draw("same");

  c->cd(5);
  probPurity->Draw();

  c->cd(4)->SetLogy();
  etra.eHistThetaAll->Draw();
  etra.eHistThetaBest->SetLineColor(kRed); etra.eHistThetaBest->Draw("same");
  thetarest->SetLineColor(kBlue); thetarest->Draw("same");

  c->cd(6);
  thetaPurity->Draw();

  c->Write();
  SafeDelete(c);
  gROOT->SetBatch(batch);
  f.Close();
}

//--------------------------------------------------------------------------------------
void EdbScanTracking::TrackAli(EdbPVRec &ali, TEnv &env)
{
    EdbTrackAssembler etra;
  
    etra.eCond.SetSigma0(        env.GetValue("fedra.track.Sigma0"         , "3 3 0.005 0.005") );
    etra.eCond.SetPulsRamp0(     env.GetValue("fedra.track.PulsRamp0"      , "15 20") );
    etra.eCond.SetPulsRamp04(    env.GetValue("fedra.track.PulsRamp04"     , "15 20") );
    etra.eCond.SetDegrad(        env.GetValue("fedra.track.Degrad"          , 4) );
    etra.eCond.SetRadX0(         env.GetValue("fedra.track.RadX0"          , 5810.) );
      
    etra.eDTmax                 = env.GetValue("fedra.track.DTmax"          ,     0.07 );
    etra.eDRmax                 = env.GetValue("fedra.track.DRmax"          ,    45.   );
    etra.eDZGapMax              = env.GetValue("fedra.track.DZGapMax"       ,  5000.   );
    etra.eProbMin               = env.GetValue("fedra.track.probmin"        ,  0.001   );
    
    bool        do_misalign     = env.GetValue("fedra.track.do_misalign"    ,      0   );
    int         npass           = env.GetValue("fedra.track.npass"          ,      1   );
    float       misalign_offset = env.GetValue("fedra.track.misalign_offset",    500.  );
    //bool        do_local_corr   = env.GetValue("fedra.track.do_local_corr"  ,      1   );
    bool        eDoRealign      = env.GetValue("fedra.track.do_realign"     ,      0   );
    bool        do_comb         = env.GetValue("fedra.track.do_comb"        ,      0   );
    eNsegMin                    = env.GetValue("fedra.track.NsegMin"        ,      2   );
    float       momentum        = env.GetValue("fedra.track.momentum",     2. );
    etra.InitTrZMap(  env.GetValue("fedra.track.TrZmap", "2400 0 120000   2000 0 100000   30" ) );

    EdbAffine2D misalign[60];
    if(do_misalign) {
        //           1 2 3 4 5 6  7  8  9
      int dx[9] = {0,0,1,1,1,0,-1,-1,-1};
      int dy[9] = {0,1,1,0,-1,-1,-1,0,1};
      for(int i=0; i<60; i++) {
        misalign[i].ShiftX( dx[i%9] * misalign_offset );
        misalign[i].ShiftY( dy[i%9] * misalign_offset );
        if(gEDBDEBUGLEVEL>1) printf("%d |  %d  %d\n",i, dx[i%9], dy[i%9]);
      }
    }
  
    int npl = ali.Npatterns();
  
    // read segments and use them for tracking
    for(int ipass=0; ipass<npass; ipass++) {
      if(gEDBDEBUGLEVEL>1) printf("\n\n*************** ipass=%d ************\n",ipass);
      etra.eCollisionsRate=0;
      for(int i=0; i<npl; i++) {
        
        EdbPattern &p = *ali.GetPattern(i);
        
        //p.SetZ(plate->Z());
        p.SetSegmentsZ();
        p.SetID(i);
        p.SetPID(i);
        p.SetSegmentsPID();
//        p.SetSegmentsPlate(id->ePlate);
        if(gEDBDEBUGLEVEL>1) printf("pattern with z: %f\n", p.Z());
      
        if(do_misalign) {
          p.Transform(&misalign[i]);
          Log(2,"EdbScanTracking::TrackSetBT","apply misalignment of %f",misalign_offset);
        }
  
        if(i>0) etra.ExtrapolateTracksToZ(p.Z());
        //if( eDoRealign && i==1 ) etra.CheckPatternAlignment(p,1);
        //if( eDoRealign && i>1  ) etra.CheckPatternAlignment(p,2);
        etra.FillTrZMap();
        etra.AddPattern(p);
      }
    }

    int ntr = etra.Tracks().GetEntriesFast();
    
    for(int i=0; i<ntr; i++) {
      EdbTrackP *t = (EdbTrackP *)(etra.Tracks().At(i));
      if(t->N()<eNsegMin)  t->SetFlag(-10);
    }
    
    etra.SetSegmentsErrors();
    etra.FitTracks();

    TObjArray selectedTracks(ntr);
    if(do_comb) {
      etra.CombTracks(selectedTracks);
    } else {
      int cnt=0;
      for(int i=0; i<ntr; i++) {
        EdbTrackP *t = (EdbTrackP *)(etra.Tracks().At(i));
        if(t->Flag()!=-10)  {
          t->SetID(cnt++);
          t->SetCounters();
          t->SetSegmentsTrack();
          t->SetP(momentum);
          selectedTracks.Add(t);
        }
      }
    }
    
    EdbID idset;
    EdbDataProc::MakeTracksTree( selectedTracks, 0., 0., Form("b%s.trk.root", idset.AsString()) );
    TFile f( Form("b%s.trk.root", idset.AsString()) ,"UPDATE");
    env.Write();
    f.Close();
    
    SaveHist(idset,etra);
}

