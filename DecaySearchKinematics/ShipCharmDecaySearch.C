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

float ShipCharmDecaySearch::IPtoVertex(TVector3 vertexpos, TVector3 trackstartpos, float tracktx, float trackty){
 //'''Impact parameter of track with respect to primary vertex'''
 
 float dz = vertexpos(2) - trackstartpos(2);
 float ipx = tracktx * dz + trackstartpos(0);
 float ipy = trackty * dz + trackstartpos(1);

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
