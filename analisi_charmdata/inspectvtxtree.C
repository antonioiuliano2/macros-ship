using namespace ROOT;
void createchain(){

  //create the tchain

  TChain chain("vtx");
  chain.Add("firstquarter/vertextree_test.root");
  chain.Add("secondquarter/vertextree_test.root");
  chain.Add("thirdquarter/vertextree_test.root");
  chain.Add("fourthquarter/vertextree_test.root");

  TH1I *hquarter = new TH1I("hquarter","Quarter of reconstruction",5,0,5);
  TH2F *hxy = new TH2F("hxy","Position of reconstructed vertices;x[#mum];y[#mum]",120,0,120000,100,0,100000);

  float vx, vy;

  chain.SetBranchAddress("vx",&vx);
  chain.SetBranchAddress("vy",&vy);

  int nevents = chain.GetEntries();

  for (int i = 0; i< nevents; i++){
    chain.GetEvent(i);
    hxy->Fill(vx,vy);
    hquarter->Fill(chain.GetTreeNumber());
  }
  TCanvas *c = new TCanvas(); 
  hxy->Draw("COLZ");
  TCanvas *c2 = new TCanvas();
  hquarter->Draw();
}

void inspectvtxtree(){
  //importing tree as a dataframe
  RDataFrame dvertices("vtx","mergedvertextree.root");
  //function to associate quarter (needed for input tracks)
  auto recognizequarter = [](float x, float y){
   if (x < 62500 && y <50000) return 1;
   else if (x > 62500 && y < 50000) return 2;
   else if (x < 62500 && y > 50000) return 3;
   else if (x > 62500 && y > 50000) return 4;
   //it should never return 0
   else return 0; 

  };

  auto d1 = dvertices.Define("quarter",recognizequarter, {"vx","vy"});
  //plotting information
  //position of reconstructed vertices
  auto hxy = d1.Histo2D({"hxy", "Reconstructed vertices position;x[#mum];y[#mum]",125,0,125000,100,0,100000},"vx","vy");
  auto hz = d1.Histo1D({"hz","Reconstructed vertex position;z[#mum]",110,-100000,10000},"vz");
  //molteplicity
  auto hn = d1.Histo1D({"hn","Number of tracks in vertex;ntracks",40,0,40},"n");
  //quarter according to vertex position (due to reconstruction separated in 4 files) 
  auto hquarter = d1.Histo1D("quarter");
  //drawing the canvases
  TCanvas *c1 = new TCanvas();
  hquarter->DrawClone();
  TCanvas *c2 = new TCanvas();
  c2->Divide(1,2);
  c2->cd(1);
  hxy->DrawClone("COLZ");
  c2->cd(2);
  hz->DrawClone();
  TCanvas *c3 = new TCanvas();
  hn->DrawClone();

  //d1.Snapshot("vtx","mergedvertextree_labelled.root"); //not working, RDataFrames do not like TClonesArray, same as FairShip
}
