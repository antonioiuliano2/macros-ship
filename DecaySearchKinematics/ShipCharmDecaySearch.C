#include "ShipCharmDecaySearch.h"
#include "TLorentzVector.h"
//***************Track kinematics which require only angles and positions information********

float ShipCharmDecaySearch::AverageKinkAngle(float parenttx, float parentty,const float* daughterstx, const float* daughtersty, int ndaughters){
	float kink=0.;

	//loop on daughters
	for (int i = 0; i < ndaughters; i++){
 
		float daughtertx = daughterstx[i];
		float daughterty = daughtersty[i];

  		kink += TMath::Sqrt(pow(parenttx - daughtertx,2) + pow(parentty - daughterty,2));
 }
 
 	return kink/ndaughters;
}

float ShipCharmDecaySearch::KinkAngle(float parenttx, float parentty, float daughtertx, float daughterty){
    //only one daughter, average kinkangle equal to kink angle
    return AverageKinkAngle(parenttx,parentty,&daughtertx,&daughterty,1);

}

std::vector<float> ShipCharmDecaySearch::FedraTrackKink(EdbTrackP* mytrack){
 std::vector<float> kinksearchresult; //first is index, second maximum kink, third rmax
 int nseg = mytrack->N();
 float kinkangles[nseg-1];
 //loop on subsequent segment pairs of the tracks
 for (int i = 0; i < nseg-1; i++){
  EdbSegP *firstseg = (EdbSegP*) mytrack->GetSegment(i);
  EdbSegP *secondseg = (EdbSegP*) mytrack->GetSegment(i+1);

  kinkangles[i]=(KinkAngle(firstseg->TX(), firstseg->TY(), secondseg->TX(), secondseg->TY()));
 }
 //getting maximum and rms
 float deltathetamax = TMath::MaxElement(nseg-1, kinkangles);
 int locmax = TMath::LocMax(nseg-1,kinkangles);
 float deltathetarms = TMath::RMS(nseg-1, kinkangles);
 float rmax = deltathetamax/deltathetarms;

 float kinkrmin = 3.;
 if (rmax > kinkrmin){ //kink found, return kink location
     kinksearchresult.push_back(locmax);
     kinksearchresult.push_back(deltathetamax);
     kinksearchresult.push_back(rmax);
     return kinksearchresult;
 }
 else{ //kink not found, loc is -1
     kinksearchresult.push_back(-1);
     kinksearchresult.push_back(deltathetamax);
     kinksearchresult.push_back(rmax);
     return kinksearchresult;
 }
}

float ShipCharmDecaySearch::IPtoVertex(TVector3 vertexpos, TVector3 trackstartpos, float tracktx, float trackty){
 //'''Impact parameter of track with respect to primary vertex'''
//from EdbEDAUtil CalcIP
 float ax = tracktx;
 float ay = trackty;

 float x = vertexpos(0);
 float y = vertexpos(1);
 float z = vertexpos(2);
 
 float bx = trackstartpos(0) - ax * trackstartpos(2) ;
 float by = trackstartpos(1) - ay * trackstartpos(2);

 
 float a;
 float r;
 float xx,yy,zz;
	
 a = (ax*(x-bx)+ay*(y-by)+1.*(z-0.))/(ax*ax+ay*ay+1.);
 xx = bx +ax*a;
 yy = by +ay*a;
 zz = 0. +1.*a;
 r = sqrt((xx-x)*(xx-x)+(yy-y)*(yy-y)+(zz-z)*(zz-z));

 return r;
}

float ShipCharmDecaySearch::TransverseIPtoVertex(TVector3 vertexpos, TVector3 trackstartpos, float tracktx, float trackty){
 //'''Impact parameter of track with respect to primary vertex'''
 
 float dz = vertexpos(2) - trackstartpos(2);
 float ipx = tracktx * dz + trackstartpos(0) -vertexpos(0);
 float ipy = trackty * dz + trackstartpos(1) - vertexpos(1);

 float ip = TMath::Sqrt(pow(ipx,2)+pow(ipy,2));

 return ip;
}

//************particle kinematics which require momentum information**********************

float ShipCharmDecaySearch::GetpT(float parenttx, float parentty, float daughtertx, float daughterty, float daughtermomentum){
//transverse momentum of daughter with respect to mother direction
 
 float costh = (1+parenttx*daughtertx + parentty*daughterty)/TMath::Sqrt((1+parenttx*parenttx+parentty*parentty)*(1+daughtertx*daughtertx+daughterty*daughterty));

 float pT = daughtermomentum * TMath::Sin(TMath::ACos(costh));
 
 return pT;

}

float ShipCharmDecaySearch::InvariantMass(std::vector<TVector3> momenta, std::vector<float> masses){
    int nparticles = momenta.size();
    TLorentzVector ptot = TLorentzVector(0.,0.,0.,0.);

    for (int i = 0; i<nparticles;i++){
        TLorentzVector p_i = TLorentzVector(momenta[i],masses[i]);
        ptot += p_i;
    }
    float invmass = ptot.Mag();
    return invmass;
}

float ShipCharmDecaySearch::InvariantMass(TVector3 momentum1, TVector3 momentum2, float mass1, float mass2){
 //now simplified by using TLorentzVectors
 TLorentzVector p1 = TLorentzVector(momentum1, mass1);
 TLorentzVector	p2 = TLorentzVector(momentum2, mass2);

 TLorentzVector ptot = p1 + p2;

 float invmass = ptot.Mag();

 return invmass;
}

TRotation* ShipCharmDecaySearch::rotation_to_beam() {
    //if already done, do not repeat it
    if (fbeamrotation) return fbeamrotation;
    //if attribute do not exists, create it and apply the rotations
    fbeamrotation = new TRotation();
    fbeamrotation -> RotateX(fbeamTY);
    double phi = TMath::ATan(TMath::Tan(fbeamTX)*TMath::Cos(fbeamTY));
    fbeamrotation -> RotateY(-phi);
    
    return fbeamrotation;
    
}
//TRotation* ShipCharmDecaySearch::fbeamrotation = nullptr; //it has to exist

ClassImp(ShipCharmDecaySearch)
