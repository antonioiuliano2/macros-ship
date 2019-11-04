//always useful to collect some formulas here
//vector = TLorentzVector(0,0,400.,400.)

/////////////////OPERATIONS WITH TRIVECTORS/////////////////////////////////

double IPtoTarget(TVector3 target, TVector3 vertex, TVector3 trackmom){

double p = trackmom.Mag();

TVector3 normmom;
TVector3 verticesdist;

normmom = trackmom * (1./p);
verticesdist = target - vertex;

double delta = verticesdist.Dot(normmom);

TVector3 ipvector;
ipvector = verticesdist - (delta * normmom);

return ipvector.Mag();
}

/////////////////OPERATIONS WITH LORENTZ VECTORS////////////////////////////

double KineticEnergy(TLorentzVector *vector){
 double energy = vector->E();
 double momentum = vector->P();
 double mass = TMath::Sqrt(pow(energy,2)-pow(momentum,2));

 double kineticenergy = energy - mass;

 return kineticenergy;
}


void Boost(TLorentzVector &vector, double vx, double vy, double vz){
 
  vector.Boost(vx,vy,vz); //from the rod frame to the original frame

 }
