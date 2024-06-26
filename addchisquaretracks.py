import numpy as np
import pandas as pd
import ROOT as r
from ROOT.TMath import Sqrt, Cos, Sin, ATan2

#degradation functions
def SigmaX(ax):
  '''computing angular degradation of error in x position'''
  eSigmaX0 = 5
  eDegrad = 5    
  return eSigmaX0*(1. + np.abs(ax)*eDegrad)

def SigmaY(ay):
  '''computing angular degradation of error in y position'''
  eSigmaY0 = 5
  eDegrad = 5  
  return eSigmaY0*(1. + np.abs(ay)*eDegrad)

def SigmaTX(ax):
  '''computing angular degradation of error in y position'''
  eSigmaTX0 = 0.003
  eDegrad = 5  
  return eSigmaTX0*(1. + np.abs(ax)*eDegrad)

def SigmaTY(ay):
  '''computing angular degradation of error in y position'''
  eSigmaTY0 = 0.003
  eDegrad = 5    
  return eSigmaTY0*(1. + np.abs(ay)*eDegrad)

def getcovmatrix(tx,ty):
  '''trying to get covariance matrix'''  
  theta = Sqrt( tx*tx + ty*ty )
  sx    = SigmaX(theta)      
  sy    = SigmaY(0)
  stx   = SigmaTX(theta)
  sty   = SigmaTY(0)

  cov = r.TMatrixD(5,5)
  cov[0,0] = sx*sx
  cov[1,1] = sy*sy
  cov[2,2] = stx*stx
  cov[3,3] = sty*sty
  cov[4,4] = 1             
     
  Phi = ATan2(ty,tx)
  t = r.TMatrixD(5,5)
  tt = r.TMatrixD(5,5)
 
  t[0,0] =  Cos(Phi)
  t[0,1] = -Sin(Phi)
  t[1,0] =  Sin(Phi)
  t[1,1] =  Cos(Phi)
  t[2,2] =  Cos(Phi)
  t[2,3] = -Sin(Phi)
  t[3,2] =  Sin(Phi)
  t[3,3] =  Cos(Phi)
  t[4,4] =  1.
  tt[0,0] =  t(0,0)
  tt[1,0] =  t(0,1)
  tt[0,1] =  t(1,0)
  tt[1,1] =  t(1,1)
  tt[2,2] =  t(2,2)
  tt[3,2] =  t(2,3)
  tt[2,3] =  t(3,2)
  tt[3,3] =  t(3,3)
  tt[4,4] =  t(4,4)

  #covfinal = t*(cov*tt)
  covtemp = r.TMatrixD(5,5)
  covfinal = r.TMatrixD(5,5)

  covtemp.Mult(cov,tt)
  covfinal.Mult(t,covtemp)

  return covfinal

def CHI2(itrack, df,trackdf):
  '''compute CHI2 for track itrack in dataframe df'''
  chi2 = 0

  dftracked = df.query("FEDRATrackID=={}".format(itrack))
  dftracked = dftracked.sort_values("PID",ascending=False)
  dftracked = dftracked.reset_index()
  #computing chisquare with respect to each base-track
  for tx,ty in zip(dftracked["TX"],dftracked["TY"]):
   segmentcov = getcovmatrix(tx,ty)
   #retrieving sigma for TX and TY
   stx = segmentcov[2,2]
   sty = segmentcov[3,3]

   dtx = tx-trackdf.loc[itrack,"TX"]
   dty = ty-trackdf.loc[itrack,"TY"]
   #NOTICE that after computation of the chi square, a square root is applied.
   #We are computing the average of the chi for each base-track
   chi2 += np.sqrt(dtx*dtx/stx + dty*dty/sty )
  #averaging over all base-tracks
  
  return chi2/len(dftracked)

df = pd.read_csv("b000001_withvertices.csv")
df = df.query("FEDRATrackID>=0") #only segments associated with volume tracks

quarterfilenames = ["firstquarter","secondquarter","thirdquarter","fourthquarter"]
trackdflist = []
for filename in quarterfilenames:
 trackdflist.append(pd.read_csv("trackdf_"+filename+".csv"))

def CHI2_quarter(quarter, itrack, df, trackdflist):
    '''
    retrieve chi2 of a track with id itrack
    taking into account different quarters.
    quarter: quarter you are referring to (1,2,3,4)
    itrack: Track reconstruction index
    df: big dataframe with all segments associated to a volume track
    trackdflist: list of dataframes with track reconstruction
    '''
    trackdf = trackdflist[quarter-1]
    return CHI2(itrack, df, trackdf)
#example

examplechi = CHI2_quarter(1,61,df,trackdflist)
