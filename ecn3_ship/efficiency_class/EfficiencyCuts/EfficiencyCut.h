//class to have all selections as functions, so we can generalize the nutau, numu, nue, charm, etc. cases

class EfficiencyCut: public TObject{
    public:
        EfficiencyCut(const char* inputfilename, const char* geofilename);
        
        void GetEntry(int ientry){
            reader.SetEntry(ientry);
            fmuonID = 1; //1 by default, will be overwritten if other muons need to be followed
            //clearing all lists
            primaryvisible.clear();
            daughters.clear();
            visibledaughters.clear();
            npl_nudaughterID.clear();
            measured_momentum_nudaughterID.clear();
        }
        int GetEntries(){return reader.GetEntries();}
        double GetEventWeight(){return tracks[0].GetWeight();} 
        //efficiency functions
        int GeometricalEfficiency(double offsetxy = 0.1, int Nminplates = 4); //is the vertex well contained in the brick?
        bool VisibleVertexLocation(double maxtantheta = 1., double minmomentum = 1.);//checking charge and momentum of neutrino daughters
        bool DecaySearch(bool tausim, double maxdl = 0.4, double minkinkangle = 0.02, double minip = 10e-4, double mindaumomentum = 0.1, double maxdautantheta = 1.);
        int MCSmeasurement(double maxres = 0.3, int nplmin = 3);
        bool SpectrometerAcceptance(double posres = 100.*1e-4, double sagittares = 0.02122);
        bool DecaySpectrometerAcceptance(double posres, double sagittares);
        int DecayChannel();
        //fill and drawing
        void FillHistograms(TH1D *hnuP, TH2D *hq2_x, TH1D* hlP); //fill histograms for this event


        TFile* simfile;
        TTreeReader reader;
        TTreeReaderArray<ShipMCTrack> tracks;
        TTreeReaderArray<TargetPoint> targetpoints;
        TTreeReaderArray<ShipRpcPoint> rpcpoints;
        TTreeReaderArray<strawtubesPoint> strawtubespoints;

        int fmuonID;

        map<int,ROOT::RVecD> energies_q2x;
        ROOT::RVecI primaryvisible; //from neutrino vertex
        ROOT::RVecI daughters; //from decay (of charmed hadron or tau lepton)
        ROOT::RVecI visibledaughters; //from decay (of charmed hadron or tau lepton) visible daughters
        map<int,int> npl_nudaughterID;
        map<int,bool> measured_momentum_nudaughterID;        
        TF2 *fOPERA_mcsres;

        TH1D *hdeltaTX;
        TH1D *hdeltaTY;

};