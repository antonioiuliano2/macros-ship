//class to have all selections as functions, so we can generalize the nutau, numu, nue, charm, etc. cases

class EfficiencyCut: public TObject{
    public:
        EfficiencyCut(const char* inputfilename, const char* geofilename);
        
        void GetEntry(int ientry){
            reader.SetEntry(ientry);
            primaryvisible.clear();
            npl_nudaughterID.clear();
            measured_momentum_nudaughterID.clear();}
        int GetEntries(){return reader.GetEntries();}
        double GetEventWeight(){return tracks[0].GetWeight();} 
        //efficiency functions
        int GeometricalEfficiency(double offsetxy = 0.1, int Nminplates = 4); //is the vertex well contained in the brick?
        bool VisibleVertexLocation(double maxtantheta = 1., double minmomentum = 1.);//checking charge and momentum of neutrino daughters
        int MCSmeasurement(double maxres = 0.3, int nplmin = 3);
        bool SpectrometerAcceptance(double posres = 100.*1e-4, double sagittares = 0.02122);
        //fill and drawing
        void FillHistograms(TH1D *hnuP, TH2D *hq2_x); //fill histograms for this event


        TFile* simfile;
        TTreeReader reader;
        TTreeReaderArray<ShipMCTrack> tracks;
        TTreeReaderArray<TargetPoint> targetpoints;
        TTreeReaderArray<ShipRpcPoint> rpcpoints;
 
        ROOT::RVec<int> primaryvisible;
        map<int,int> npl_nudaughterID;
        map<int,bool> measured_momentum_nudaughterID;        
        TF2 *fOPERA_mcsres;

};