using namespace ROOT;
void read_eduardmuonfile(){
    //first, accessing the file from EOS
    TFile *inputfile = TFile::Open("root:://eosuser.cern.ch//eos/user/e/edursov/ship_data/snd_ship_ecn3/muons_3_planes.root");
    TTree *cbmsim = (TTree*) inputfile->Get("cbmsim");
    cout<<"Reading a tree with "<<cbmsim->GetEntries()<<endl;
    TCanvas *cxy = new TCanvas();
    cbmsim->Draw("sco_1Point.fY:sco_1Point.fX>>hxy(140,-600,800,140,-600,800)","MCTrack[0].fW","COLZ");
    TCanvas *cxy_zoomed = new TCanvas();
    cbmsim->Draw("sco_1Point.fY:sco_1Point.fX>>hxy_zoomed(400,-200,200,400,-200,200)","MCTrack[0].fW","COLZ");

    TCanvas *cP = new TCanvas();
    cbmsim->Draw("TMath::Sqrt(sco_1Point.fPx*sco_1Point.fPx+sco_1Point.fPy*sco_1Point.fPy+sco_1Point.fPz*sco_1Point.fPz)","MCTrack[0].fW");

    TCanvas *cpdgcode = new TCanvas();
    cbmsim->Draw("sco_1Point.fPdgCode","MCTrack[0].fW");
}