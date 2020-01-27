//Simple interface to ShowerRec library. It allows to search for shower looking for segments stored in a EdbPVRec (Created on Thursday 23 January 2020 by A.Iuliano)
#ifndef SimpleShowerRecInterface_h
#define SimpleShowerRecInterface_h
class SimpleShowerRecInterface:public TObject {

  private:
    EdbPVRec* eEdbPVRec; //container of segments and patterns of all involved plates

  public:
   SimpleShowerRecInterface(){eEdbPVRec=0;};
   ~SimpleShowerRecInterface(){};
   
   //void SetAlgoParameter(Double_t paravalue, Int_t paranr) {eEdbPVRec->SetAlgoParameter(Double_t paravalue, Int_t paranr);}; 
   /*setting parameters of EdbShowerRec
   if (paranr==0) eAlgoParameterConeRadius=paravalue;
   else if (paranr==1) eAlgoParameterConeAngle=paravalue;
   else if (paranr==2) eAlgoParameterConnectionDR=paravalue;
   else if (paranr==3) eAlgoParameterConnectionDT=paravalue;
   else if (paranr==4) eAlgoParameterNPropagation=paravalue;
   */
   void BuildPVRec(const int ibrick, const int nplates, int* PID,  TString *selections, char* outputfilename); //build a new PVRec according to given selection, write it to file
   void loadcouples(const int ibrick, const int nplates, int* PID, TString *selections); //load couples according to given selection, called by BuildPVRec, not to be called externally (private?)
   void LoadPVRec(TFile *inputfile);// load an EdbPVRec saved into a file

   void RecoFromTrack(int ntracks, int *tracklist, int* iseglist, char * filename = "linked_tracks.root"); //start reconstruction from a given track
   TObjArray* DrawShower(int ishower, char * showerfilename = "Shower.root");//draw only one entry from the tree
   void DrawAllShowers(char * showerfilename = "Shower.root"); //draw all reconstructed showers


   ClassDef(SimpleShowerRecInterface,1)

};

#endif
