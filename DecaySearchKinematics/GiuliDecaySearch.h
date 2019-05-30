#include "TObject.h"
#include "TVector3.h"
#include "TRotation.h"
#include <vector>

/*
Class created by Antonio to do kinematic operations for decay search. 
Structure is a simple ROOT class, which can be compiled with ACLIC:

.L GiuliDecaySearch.C+

then loaded with

gSystem->Load("GiuliDecaySearch_C.so")

and used


*/
class GiuliDecaySearch:public TObject {

	public:

	GiuliDecaySearch( float beamTX = 0., float beamTY = 0.) { fbeamTX = beamTX; fbeamTY = beamTY; fbeamrotation = NULL;} 
    ~GiuliDecaySearch(){  if (fbeamrotation) delete fbeamrotation;}
    //static functions are not connected to a particular object, can be launched simply by using GiuliDecaySearch::Myfunc(argument)
	static float AverageKinkAngle(float parenttx, float parentty, const float* daughterstx, const float* daughtersty, int ndaughters);

	static float GetpT(float parenttx, float parentty, float daughterttx, float daughtertty, float daughtermomentum);

	static float IPtoVertex(TVector3 vertexpos, TVector3 trackstartpos, float tracktx, float trackty);
     
    static float InvariantMass(std::vector<TVector3> momenta, std::vector<float> masses);

    static float InvariantMass(TVector3 momentum1, TVector3 momentum2, float mass1, float mass2);

    //functions to take into account beam angle need the angles, so they are not static
    TRotation* rotation_to_beam(); 

	ClassDef(GiuliDecaySearch,1) // List of functions for Decay Search
	private:
    float fbeamTX;
    float fbeamTY;
    TRotation* fbeamrotation; 

	};