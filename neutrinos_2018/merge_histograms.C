void merge_histograms(){

TFile *f10 = TFile::Open("root://eospublic.cern.ch//eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_10.0_withCharm_nu.root"); //10 GeV cut
TFile *f1 = TFile::Open("root://eospublic.cern.ch//eos/experiment/ship/data/Mbias/background-prod-2018/pythia8_Geant4_1.0_withCharm_nu.root"); //1 GeV cut

Int_t neutrinos[6] = {12,14,16,-12,-14,-16}; //same numbers as the names of the histograms to get
TH1D *h10[6];
TH1D *h1[6];

TFile *outfile = new TFile("/home/utente/Lavoro/Analisi/macro_SHiP/temp_outputs/neutrinos_merged.root","RECREATE");

//TFile *outfile = new TFile("/afs/cern.ch/work/a/aiuliano/public/macros-ship/temp_outputs/neutrinos_merged.root","RECREATE");
TH1D *hsum[6];

TH2D *h2D_10[6];
TH2D *h2D_1 [6];
TH2D *h2D_sum[6];

TCanvas *c[6];
TCanvas *c_2D[6];

for (int ineutrino = 0; ineutrino < 6; ineutrino++){ //loop on neutrino flavours

 if (neutrinos[ineutrino]>0){
  h10[ineutrino] = (TH1D*) f10->Get(Form("10%i", neutrinos[ineutrino])); //getting the histogram from thomas file
  h1[ineutrino] = (TH1D*) f1->Get(Form("10%i", neutrinos[ineutrino]));
  hsum[ineutrino] = new TH1D(Form("10%i",neutrinos[ineutrino]),Form("Merge of histograms for neutrino %i",neutrinos[ineutrino]),400,0,400);
 }
 else {
  h10[ineutrino] = (TH1D*) f10->Get(Form("20%i",-1* neutrinos[ineutrino]));
  h1[ineutrino] = (TH1D*) f1->Get(Form("20%i",-1* neutrinos[ineutrino]));
  hsum[ineutrino] =  new TH1D(Form("20%i",-1 * neutrinos[ineutrino]),Form("Merge of histograms for neutrino %i",neutrinos[ineutrino]),400,0,400);
 }  

 if (neutrinos[ineutrino]>0){
  h2D_10[ineutrino] = (TH2D*) f10->Get(Form("12%i", neutrinos[ineutrino]));
  h2D_1[ineutrino] = (TH2D*) f1->Get(Form("12%i", neutrinos[ineutrino]));
  h2D_sum[ineutrino] = new TH2D(Form("12%i",neutrinos[ineutrino]), Form("log10pt vs log10p for neutrino %i", neutrinos[ineutrino]),100,-0.3,1.7,100,-2,1);
 }
 else{
  h2D_10[ineutrino] = (TH2D*) f10->Get(Form("22%i",-1* neutrinos[ineutrino]));
  h2D_1[ineutrino] = (TH2D*) f1->Get(Form("22%i",-1* neutrinos[ineutrino]));
  h2D_sum[ineutrino] = new TH2D(Form("22%i",-1 *neutrinos[ineutrino]), Form("log10pt vs log10p for neutrino %i", neutrinos[ineutrino]),100,-0.3,1.7,100,-2,1);
 } 
 
 //hsum->Add(h14_10, h14_1, 0.5, 0.5);
 for (int i = 1; i <= 400; i++){ //merging the histograms
   if (i <= 10) hsum[ineutrino]->SetBinContent(i, h1[ineutrino]->GetBinContent(i)); //below 10 GeV, only one histogram has values
   else hsum[ineutrino]->SetBinContent(i, h10[ineutrino]->GetBinContent(i));
 }

 for (int i = 1; i <= h2D_10[ineutrino]->GetNbinsX() * h2D_10[ineutrino]->GetNbinsY(); i++){
   Int_t binx = 0, biny = 0, binz = 0; 
   h2D_sum[ineutrino]->GetBinXYZ(i, binx, biny, binz); //to decompose the bin in the three axis components;

   if (binx <= 65) h2D_sum[ineutrino]->SetBinContent(i, h2D_1[ineutrino]->GetBinContent(i)); //below 10 GeV, only one histogram has values
   else h2D_sum[ineutrino]->SetBinContent(i, h2D_10[ineutrino]->GetBinContent(i));
 }
 //cout<<"Number of muon neutrinos in 5*10^13 pot: "<<h14_1->Integral(1,9) + hsum->Integral(10,400)<<endl;
 cout<<Form("Number of neutrinos of PDG %i in 5*10^13 pot: ", neutrinos[ineutrino])<<hsum[ineutrino]->Integral()<<endl;
 gStyle->SetOptStat(1111111);
 c[ineutrino] = new TCanvas();
 hsum[ineutrino]->Draw();
 c[ineutrino]->SetLogy();
 c[ineutrino]->SetLogx();
 c[ineutrino]->Print(Form("/afs/cern.ch/work/a/aiuliano/public/macros-ship/temp_outputs/neutrinos_plots/%i.png",neutrinos[ineutrino]));

 c_2D[ineutrino] = new TCanvas();
 h2D_sum[ineutrino]->Draw();
 c_2D[ineutrino]->Print(Form("/afs/cern.ch/work/a/aiuliano/public/macros-ship/temp_outputs/neutrino_plots/2D_%i.png",neutrinos[ineutrino]));
 }
outfile->Write();
outfile->Close();
}
