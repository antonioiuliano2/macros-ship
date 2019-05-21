const TRotation& opera_to_beam();
const TRotation& beam_to_opera();
const double TXfascio = -6.7E-3;
const double TYfascio = 58E-3;
///

TRotation* _opera_to_beam = NULL;
TRotation* _beam_to_opera = NULL;

const TRotation& opera_to_beam() {
    
    if (_opera_to_beam) return *_opera_to_beam;
    
    _opera_to_beam = new TRotation();
    _opera_to_beam -> RotateX(TYfascio);
    double phi = TMath::ATan(TMath::Tan(TXfascio)*TMath::Cos(TYfascio));
    _opera_to_beam -> RotateY(-phi);
    
    return *_opera_to_beam;
    
}


Int_t GetDirectionCosine(RParticle *part, Float_t cosdir[3]){
    
    Float_t p = ptrue(part);
    
    cosdir[0]=part->Px()/p;
    cosdir[1]=part->Py()/p;
    cosdir[2]=part->Pz()/p;
    
    return 0;
    
}

Int_t GetDirectionCosineCorrected(RParticle *part, Float_t cosdir[3]){
    
    GetDirectionCosine(part, cosdir);
    
    TVector3 cosdirV;
    cosdirV[0]=cosdir[0];
    cosdirV[1]=cosdir[1];
    cosdirV[2]=cosdir[2];
    
    TVector3 cosdircorrV;
    
    cosdircorrV=opera_to_beam()*cosdirV;
    
    cosdir[0]=cosdircorrV[0];
    cosdir[1]=cosdircorrV[1];
    cosdir[2]=cosdircorrV[2];
    
    return 0;
}

Float_t GetpT2ry(RParticle *parent, RParticle *daughter, Float_t psmeardau){
    
    Float_t parenttx=parent->Px()/parent->Pz();
    Float_t parentty=parent->Py()/parent->Pz();
    Float_t dautx=daughter->Px()/daughter->Pz();
    Float_t dauty=daughter->Py()/daughter->Pz();
    
    Float_t costh=(1+parenttx*dautx+parentty*dauty)/sqrt((1+parenttx*parenttx+parentty*parentty)*(1+dautx*dautx+dauty*dauty));
    
    Float_t pT2= psmeardau*TMath::Sin(TMath::ACos(costh));
    
    return pT2;
}

Float_t GetpTtrue(RParticle *parent, RParticle *daughter){
    
    Float_t psmeardau=ptrue(daughter);
    
    Float_t parenttx=parent->Px()/parent->Pz();
    Float_t parentty=parent->Py()/parent->Pz();
    Float_t dautx=daughter->Px()/daughter->Pz();
    Float_t dauty=daughter->Py()/daughter->Pz();
    
    Float_t costh=(1+parenttx*dautx+parentty*dauty)/sqrt((1+parenttx*parenttx+parentty*parentty)*(1+dautx*dautx+dauty*dauty));
    
    Float_t pT2= psmeardau*TMath::Sin(TMath::ACos(costh));
    
    // cout << LOUT(pT2) << endl;
    
    return pT2;
}

Float_t GetPTmiss(Float_t Px_had_1ry, Float_t Py_had_1ry, Float_t Px_had_2ry, Float_t Py_had_2ry){
    
    return sqrt((Px_had_1ry+Px_had_2ry)*(Px_had_1ry+Px_had_2ry)+(Py_had_1ry+Py_had_2ry)*(Py_had_1ry+Py_had_2ry));
}

Float_t GetPhi(RParticle *parent, Float_t Px_had_1ry, Float_t Py_had_1ry){
    
    Float_t cosdirparent[3] = {0};
    GetDirectionCosineCorrected(parent, cosdirparent);
    
    Float_t Norma1ry = sqrt(Px_had_1ry*Px_had_1ry+Py_had_1ry*Py_had_1ry);
    Float_t Norma2ry = sqrt(cosdirparent[0]*cosdirparent[0]+cosdirparent[1]*cosdirparent[1]);
    if (Norma1ry==0) {
        return -99;
    }
    else{
        Float_t phi = TMath::ACos((Px_had_1ry*cosdirparent[0]+Py_had_1ry*cosdirparent[1])/(Norma1ry*Norma2ry));
        phi=phi*180/TMath::ACos(-1);
        
        return phi;
    }
}


Float_t Kinkangle(RParticle *parent, RParticle *daughter[3], int ndau){
    
    Float_t kink=0., mtx=0.,mty=0.,ttx=0.,tty=0.;
    
    ttx = parent->Px()/parent->Pz();
    tty = parent->Py()/parent->Pz();
    int j=0;
    for (int i=0; i<ndau; i++) {
        if (daughter[i]) {
            mtx = daughter[i]->Px()/daughter[i]->Pz();
            mty = daughter[i]->Py()/daughter[i]->Pz();
            j++;
        }
        kink += sqrt((ttx-mtx)*(ttx-mtx)+(tty-mty)*(tty-mty));
    }
    return kink/j;
}


Float_t MinimaMassaInvariante(RParticle *tau, RParticle *daughter[3], Float_t psmeardau[3]){
    
    Float_t Mmin=0;
    
    Float_t mpi=0.1396; //GeV/c^2
    
    Float_t Etot=0;
    
    Float_t Pxtot=0,Pytot=0,Pztot=0;
    
    Float_t pTtot=0, pLtot=0;
    Float_t Ptot=0;
    
    Float_t costh=0;
    
    Float_t E=0;
    
    Float_t tautx=tau->Px()/tau->Pz();
    Float_t tauty=tau->Py()/tau->Pz();
    
    Float_t Ptau=ptrue(tau);
    
    for (int i=0; i<3; i++) {
        Float_t p=psmeardau[i];
        E=sqrt(mpi*mpi+p*p);
        Etot+=E;
        
        Float_t psmearc[3]={0};
        
        SmearMomentumComponentsNotCorrected(daughter[i], p, psmearc);
        
        Pxtot+=psmearc[0];
        Pytot+=psmearc[1];
        Pztot+=psmearc[2];
    }
    
    Ptot=sqrt(Pxtot*Pxtot+Pytot*Pytot+Pztot*Pztot);
    
    costh=(Pxtot*tau->Px()+Pytot*tau->Py()+Pztot*tau->Pz())/(Ptot*Ptau);
    
    pTtot=Ptot*(TMath::Sin(TMath::ACos(costh)));
    pLtot=Ptot*costh;
    
    cout << LOUT(pTtot) << endl;
    
    //cout << LOUT(costh) << LOUT(p) << LOUT(E) << LOUT(pT2ry) << LOUT(pL2ry) << endl;
    
    
    Float_t pTnu=-pTtot;
    
    Mmin=sqrt((sqrt(Etot*Etot-pLtot*pLtot)+abs(pTnu))*(sqrt(Etot*Etot-pLtot*pLtot)+abs(pTnu)) - (pTtot+pTnu)*(pTtot+pTnu));
    
    
    Float_t z0=(pLtot*abs(pTnu))/(sqrt(Etot*Etot-pLtot*pLtot)); //momento longitudinale del neutrino minimizzato
    
    Float_t Mmin2=sqrt((Etot+sqrt(pTnu*pTnu+z0*z0))*(Etot+sqrt(pTnu*pTnu+z0*z0))-(pLtot+z0)*(pLtot+z0)-(pTtot+pTnu)*(pTtot+pTnu));
    
    
    return Mmin;
    
}



Float_t MassaInvariante(RParticle *daughter[3], Float_t psmear[3]){
    
    Float_t mpi=0.1396; //GeV/c^2
    
    Float_t M=0;
    
    Float_t pcomp[3] = {0};
    Float_t p=0;
    Float_t E;
    Float_t Etot=0;
    Float_t ptot=0, pxtot=0, pytot=0, pztot=0;
    
    for(int i=0; i<3; i++){
        SmearMomentumComponentsNotCorrected(daughter[i],psmear[i],pcomp);
        p=sqrt(pcomp[0]*pcomp[0]+pcomp[1]*pcomp[1]+pcomp[2]*pcomp[2]);
        
        E=sqrt(mpi*mpi+p*p);
        Etot+=E;
        pxtot+=pcomp[0];
        pytot+=pcomp[1];
        pztot+=pcomp[2];
    }
    
    
    ptot=sqrt(pxtot*pxtot+pytot*pytot+pztot*pztot);
    
    if (Etot!=ptot) {
        M=sqrt(Etot*Etot-ptot*ptot);
    }
    
    return M;
}

