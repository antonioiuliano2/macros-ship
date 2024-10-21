import ROOT as r
import time,os,sys,argparse
from array import array

import gzip
import pickle
import numpy as np
import os

def get_data(file):
    with gzip.open(file,'rb') as f:
        data = pickle.load(f)
    return data

masssq = {}
PDG = r.TDatabasePDG.Instance()

def getMasssq(pid):

    apid = abs(int(pid))
    if not apid in masssq:
        masssq[apid] = PDG.GetParticle(apid).Mass()**2
    return masssq[apid]


def rotate(ctheta,stheta,cphi,sphi,px,py,pz):
  #rotate around y-axis
     px1=ctheta*px+stheta*pz
     pzr=-stheta*px+ctheta*pz
     #rotate around z-axis
     pxr=cphi*px1-sphi*py
     pyr=sphi*px1+cphi*py
     return pxr,pyr,pzr



def makeMuonDIS():

    
    parser= argparse.ArgumentParser(description='Script to generate DIS events')


    parser.add_argument(                                                 
        '-f',
        '--inputFile',
        help='''Input file to use ''')
    #parser.add_argument(
    #    '-l',
    #    '--muonsLocation',
    #    help='''Location of the simulated muons''')
    parser.add_argument(                                                 
        '-nPerJobs', 
        '--nPerJobs',
        type=int,
        help='''The number of the job ''')
    parser.add_argument(
        '-nJobs', 
        '--nJob',
        type=int,
        help=''' Number of Jobs to generate''')
    parser.add_argument(
        '-nDISPerMuon',
        '--nDIS',
        type=int,
        default=10000,
        help=''' Number of DIS per muon to generate''')



    args = parser.parse_args()
    # read file with muons hitting concrete wall/SBT and different parts of the set up
    #muonIn  = r.TFile.Open(args.inputFile, 'read')
    #muonChoice=args.muonsLocation
    #if muonChoice=="SBT":
    #    sTree=muonIn.muonsSBT
    #else:
    #    sTree=muonIn.muonsTr1
        
    sTree = get_data(args.inputFile)
    #sTree   == muonIn.muonsSBT
    #sTree=muonIn.muonsSBT
    nJob    = args.nJob
    nPerJob = args.nPerJobs
    nMult   = args.nDIS
    print(nPerJob, nJob, nMult)
    #define the start and the end event per job, the total number of events
    #nTOT = sTree.GetEntries()
    nTOT = sTree.shape[1]
    nStart = nPerJob*nJob
    nEnd   = min(nTOT,nStart+nPerJob)
    print(nTOT)
    print("The name of the file is:", args.inputFile, "The job number is", nJob, "The number per Job", nPerJob, "The number DIS", nMult)
    print("The start event is=", nStart, "The End events is=", nEnd)
    #prepare the outputFile
    # DIS event
    # incoming muon,      id:px:py:pz:x:y:z:w
    # outgoing particles, id:px:py:pz
    fout  = r.TFile('muonDis_'+str(nJob)+'.root','recreate')
    dTree = r.TTree('DIS','muon DIS')
    iMuon       = r.TClonesArray("TVectorD") 
    iMuonBranch = dTree.Branch("InMuon",iMuon,32000,-1)
    dPart       = r.TClonesArray("TVectorD") 
    dPartBranch = dTree.Branch("Particles",dPart,32000,-1)
    fcross = open('sigmadata_'+str(nJob)+'.txt', "w")



    myPythia = r.TPythia6()
     # msel 2 includes diffractive parts
    myPythia.SetMSEL(2)      
    # To get below 10 GeV, you have to change PARP(2)
    myPythia.SetPARP(2,2)     
    for kf in [211,321,130,310,3112,3122,3222,3312,3322,3334]:
       kc = myPythia.Pycomp(kf) 
       myPythia.SetMDCY(kc,1,0)



     

    R = int(time.time()%900000000)
    myPythia.SetMRPY(1,R)
    mutype = {-13:'gamma/mu+',13:'gamma/mu-'}

    # stop pythia printout during loop
    myPythia.SetMSTU(11, 11)
    print("start production ",nStart,nEnd)
    nMade = 0
    intime = time.perf_counter()    
    #for k in range(nStart, nEnd):
    imuon = nStart
    for px,py,pz,x,y,z,pid,w in zip(*sTree[:,nStart:nEnd]):
        #conversion to FairShip Coordinates
        z -= 50
        z *= 100
        x*=100
        y*=100
        pid=int(pid)
        #rc = sTree.GetEvent(k)
        #print("The number of event is",k,", time = ", time.perf_counter() - intime, " s, job = ", nJob)
        # make n events / muon
        #px, py, pz = sTree.px, sTree.py, sTree.pz
        #x, y, z = sTree.x, sTree.y, sTree.z
        #pid, w = sTree.id, sTree.w
        p = r.TMath.Sqrt(px ** 2 + py ** 2 + pz ** 2)
        #print("muon",imuon,"pID",pid,"weight",w," the momentum is=", p,". The components (px,py,pz) are", px,",",py,",",pz)
        E = r.TMath.Sqrt(getMasssq(pid) + p ** 2)
        Ecm = r.TMath.Sqrt(2*E*PDG.GetParticle(2212).Mass()) #proton mass is lower of neutron, so I use it for the cut
        if (Ecm<2.): #cut due to SetPARP(2,2)
            continue
        # px=p*sin(theta)cos(phi),py=p*sin(theta)sin(phi),pz=p*cos(theta)
        theta = r.TMath.ACos( pz / p)
        phi = r.TMath.ATan2(py, px)
        ctheta, stheta = r.TMath.Cos(theta), r.TMath.Sin(theta)
        cphi, sphi = r.TMath.Cos(phi), r.TMath.Sin(phi)
        isProton=1 #flag 
        xsec=0
        mu = array('d', [pid, px, py, pz, E, x, y, z, w,isProton,xsec])        
        muPart = r.TVectorD(11, mu)
        myPythia.Initialize('FIXT', mutype[pid],'p+', p)
        myPythia.Pylist(1)
        for a in range(nMult):
            if a==nMult/2:
                myPythia.Initialize('FIXT', mutype[pid],'n0', p)
                isProton=0
            dPart.Clear()
            iMuon.Clear()
            muPart[9]=isProton
            iMuon[0] =muPart
            myPythia.GenerateEvent()
            myPythia.Pyedit(1)
        	#xsec=myPythia.GetPARI(1)
        	#fcross.write(f"{xsec}\n")
            for itrk in range(1,myPythia.GetN()+1):
                xsec=myPythia.GetPARI(1)
                #fcross.write(f"{xsec}\n") 
                muPart[10]=xsec
                did = myPythia.GetK(itrk,2)
                dpx,dpy,dpz = rotate(ctheta,stheta,cphi,sphi,myPythia.GetP(itrk,1),myPythia.GetP(itrk,2),myPythia.GetP(itrk,3))
        	    #print(did,dpx,dpy,dpz)
                psq =   dpx**2+dpy**2+dpz**2
                E = r.TMath.Sqrt(getMasssq(did)+psq)
                m = array('d',[did,dpx,dpy,dpz,E])
                part = r.TVectorD(5,m)
                nPart = dPart.GetEntries()
                #if itrk == 1: 
                #    print("sigma_DIS = ",muPart[10])
                #print("Track ", itrk, ": ID = ",did,". E, px, py, pz = ",E,",",dpx,",",dpy,",",dpz)
                if dPart.GetSize()==nPart: dPart.Expand(nPart+10)
                dPart[nPart] = part
                if itrk == 1: fcross.write(f"{xsec}\n")
            nMade+=1
            if nMade%10000==0: print('made so far ',nMade)
            dTree.Fill()
        imuon+=1 
    fout.cd()  
    dTree.Write()
    myPythia.SetMSTU(11, 6)
    #myPythia.Pystat(2)
    
    print("created nJob ",nJob,':',nStart,' - ',nEnd," events")
if __name__ == '__makeMuonDIS__':
    makeMuonDIS()
    fcross.close()
makeMuonDIS()

