#include "EfficiencyCut.h"

//reading the file, prepare TTreeReader
EfficiencyCut::EfficiencyCut(const char* inputfilename, const char* geofilename):
 simfile(),
 reader(),
 tracks(reader,"MCTrack"),
 rpcpoints(reader,"ShipRpcPoint"),
 targetpoints(reader,"TargetPoint")
{
 simfile = TFile::Open(inputfilename);
 reader.SetTree("cbmsim",simfile);    
 //importing geofile:
 TGeoManager * tgeom = new TGeoManager("Geometry", "Geane geometry");
 tgeom->Import(geofilename);
 //fucntion for mcs
 fOPERA_mcsres = new TF2("fOPERA_mcsres","([0]+[1]*x)/TMath::Sqrt(y)+ ([2]+[3]*x)+([4]+[5]*x)*TMath::Sqrt(y)");
 fOPERA_mcsres->SetParameters(0.397,0.019,0.176,0.042,-0.014,-0.003);
 //initializing histograms
 hnuP = new TH1D("hnuP","Neutrino momentum;P[GeV/C]",400,0,400);
 hq2_x = new TH2D("hq2_x",";log10(x);log10(Q2)",15,-3,0,10,-1.5,3.5);
}

//find brick from position in the geometry
Bool_t FindBrick(Double_t x, Double_t y, Double_t z, Int_t &NWall,  Int_t &NRow, Int_t &NColumn, Int_t &NPlate)
{
 gGeoManager->FindNode(x,y,z);
 if (gGeoManager->GetLevel() == 0) return kFALSE;//we have in the cave, no mother volume present
 const char *name = gGeoManager->FindNode(x,y,z)->GetMotherVolume()->GetName(); //go there
 if(strcmp(name, "Brick") == 0){
  NPlate = gGeoManager->GetMother(0)->GetNumber()+1; //volumes numbers start from 0, instead we are used of thinking as 1 is the first plate
  NColumn = gGeoManager->GetMother(2)->GetNumber();
  NRow = gGeoManager->GetMother(3)->GetNumber();
  NWall = gGeoManager->GetMother(4)->GetNumber();

  return kTRUE;
 }
 else return kFALSE;
}

bool EfficiencyCut::GeometricalEfficiency(double offsetxy, int Nminplates){

    double vx = tracks[0].GetStartX();
    double vy = tracks[1].GetStartY();
    double vz = tracks[2].GetStartZ();

    const int Ntotplates = 59;
    double startbrickx = 0.5;
    double startbricky = 0.5;

    int InteractionWall;
    int whichrow, whichcolumn; //for findbrick function
    int whichplate = 60;

    bool foundbrick = FindBrick(vx, vy, vz, InteractionWall, whichrow, whichcolumn, whichplate);
    if (!foundbrick){
        cout<<"No brick found!, please check"<<endl;
        return false;
    }
    TGeoVolume *volEmulsion = gGeoManager->GetVolume("Emulsion");
    TGeoBBox *Emulsion = (TGeoBBox*) volEmulsion->GetShape();
    double dxtarget = Emulsion->GetDX();
    double dytarget = Emulsion->GetDY();
    //double dxtarget = 20.05;
    //double dytarget = 20.05;
    cout<<" Wall "<<InteractionWall<<" row "<< whichrow<<" column "<<whichcolumn<<" plate "<<whichplate<<endl;
    cout<<" vx "<<vx<<" vy "<<vy<<endl;
    //longitudinal edge
    if (whichplate > (Ntotplates - Nminplates +1)) return false;
    //upper transverse edges
    if (TMath::Abs(vx) > ((2*dxtarget) + startbrickx - offsetxy)) return false;
    if (TMath::Abs(vy) > ((2*dytarget) + startbricky - offsetxy)) return false;
    //lower transverse edges
    if (TMath::Abs(vx) < (startbrickx + offsetxy)) return false;
    if (TMath::Abs(vy) < (startbricky + offsetxy)) return false;
        
    else return true;

}

bool EfficiencyCut::VisibleVertexLocation(){
     //********************************SECOND CONDITION: VERTEX LOCATION: at least one track with mom > 1 GeV and tantheta<1********/
    const double maxtantheta = 1.;
    const double minmomentum = 1.;
    TDatabasePDG *pdg = TDatabasePDG::Instance();
    //getting particle charge and momentum
    int trackID = 0;
    for (auto &track:tracks){ //loop over MCTracks, find neutrino daughters
     double momentum = track.GetP();
     double charge = 0.;
     double pdgcode = track.GetPdgCode();

     double tx = track.GetPx()/track.GetPz();
     double ty = track.GetPy()/track.GetPz();
     double tantheta = TMath::Sqrt(tx * tx + ty * ty);
     if (pdg->GetParticle(pdgcode)) charge = pdg->GetParticle(pdgcode)->Charge(); 
     //is it from the primary neutrino interaction?
     if(track.GetMotherId()==0 && TMath::Abs(charge)>0){ //charged neutrino daughter  
          if(tantheta<maxtantheta && momentum > minmomentum){ //clearly visible
            primaryvisible.push_back(trackID); //I do not add the tau lepton or charmed hadron since it decays too soon
            npl_nudaughterID[trackID] = 0;
            //cout<<ientry<<" "<<pdgcode<<" "<<momentum<<" "<<tantheta<<endl;
            } //end tantheta condition
      } //end charged daugter condition
      trackID++;
     }//end over track loop
    if (primaryvisible.size()>0) return true;
    else return false;
}

int EfficiencyCut::MCSmeasurement(double maxres, int nplmin){
    //this function must be called after VisibleVertexLocation
    //looping over targetpoints
    if (targetpoints.GetSize() == 0) return -1;
    int vertexbrickID = targetpoints[0].GetDetectorID()/1e+4; //the first hit stored is the closest temporarily to nu interaction
    for (const TargetPoint& targetpoint: targetpoints){
        int detectorID = targetpoint.GetDetectorID();
        int trackID = targetpoint.GetTrackID();
        int brickID = detectorID/1e+4;         
        //selection to accept hit for counting
        if (detectorID%2 == 1.) continue; //only accept one side of emulsion film      
        if (brickID != vertexbrickID) continue; //hit does not came from the same vertex as neutrino
        if (npl_nudaughterID.find(trackID) == npl_nudaughterID.end()) continue;//hit is not from a neutrino daughter

        npl_nudaughterID[trackID]++; //all conditions are accepted, increasing counter
    }
    //checking targetpoints
    for ( const auto &myPair : npl_nudaughterID ) {
         int npl = myPair.second;
         int trackID = myPair.first;

         double momentum = tracks[trackID].GetP();

         if(npl<nplmin) continue;

         double res = fOPERA_mcsres->Eval(momentum, npl); //sigma(1/P)/(1/P)
         cout<<"test "<<res<<" momentum "<<momentum<<" npl "<<npl<<endl;   
         if (res < maxres) measured_momentum_nudaughterID[trackID] = true;
         else measured_momentum_nudaughterID[trackID] = false;
   }
   return 0;
}