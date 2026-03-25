//una minimacro che uso per studiare le distribuzioni dei vertici nel bersaglio di SHiP (modificata il 17 Giugno 2017)

void dafile();
void dafunzione();

void studio_vertici_charm(){
  dafile();
  // dafunzione();
}

void dafile(){
  //TFile *inputfile = TFile::Open("~/Lavoro/Analisi/Archivio_cronologico/Marzo_2018/cascade_vs_primary.root");
  TFile *inputfile = TFile::Open("/eos/user/a/aiuliano/public/sims_FairShip/sim_charm/bigtarget_v2/ship.conical.Pythia8CharmOnly-TGeant4.root");
 // TCanvas *c1 = (TCanvas*) inputfile->Get("c1");
  //TH1D* hprimary =  (TH1D*) c1->GetPrimitive("hprimary");
  //TH1D* hsecondary = (TH1D*) c1->GetPrimitive("hsecondary");
  TTree *tree = (TTree*) inputfile->Get("cbmsim");
  tree->Draw("MCTrack[0].fStartZ+230.8>>hprimary(1000)","MCEventHeader.fRunId == 1");
  tree->Draw("MCTrack[0].fStartZ+230.8>>hsecondary(1000)","MCEventHeader.fRunId > 1");
  tree->Draw("MCTrack[0].fStartZ+230.8>>hall(1000)");   

  tree->Draw("TMath::ATan(TMath::Sqrt(pow(MCTrack[0].fPx/MCTrack[0].fPz,2)+pow(MCTrack[0].fPy/MCTrack[0].fPz,2)))>>hthetaprimary","MCEventHeader.fRunId == 1");
  tree->Draw("TMath::ATan(TMath::Sqrt(pow(MCTrack[0].fPx/MCTrack[0].fPz,2)+pow(MCTrack[0].fPy/MCTrack[0].fPz,2)))>>hthetacascade(300,0,0.03)","MCEventHeader.fRunId > 1");
 
  TH1D* hthetaprimary = (TH1D*)gDirectory->Get("hthetaprimary");
  TH1D* hthetacascade = (TH1D*)gDirectory->Get("hthetacascade");

  TCanvas *ctheta = new TCanvas();
  ctheta->Divide(1,2);
  ctheta->cd(1);
  hthetaprimary->Draw();
  ctheta->cd(2);
  hthetacascade->Draw();
  TCanvas *cvertex = new TCanvas();
  TH1D* hall = (TH1D*) gDirectory->Get("hall");
  TH1D* hprimary = (TH1D*) gDirectory->Get("hprimary");
  TH1D* hsecondary = (TH1D*) gDirectory->Get("hsecondary");

  hsecondary->SetLineColor(kRed);
  hsecondary->SetLineWidth(2);
  hsecondary->SetLineStyle(2);
  //c1->Clear();
  hprimary->Draw();
  hprimary->GetXaxis()->SetTitle("cm");
  hprimary->Fit("expo");
  hsecondary->Draw("SAMES");
  
  TLine *l0 = new TLine(2.8,0,2.8,580);
  TLine *l1 = new TLine(5.6,0,5.6,580);
  TLine *l2 = new TLine(11.2,0,11.2,580);
  TLine *l3 = new TLine(16.8,0,16.8,580);
  TLine *l4 = new TLine(22.4,0,22.4,580);
  TLine *l5 = new TLine(28.0,0,28.0,580);

  l0->Draw();
  l1->Draw();
  l2->Draw();
  l3->Draw();
  l4->Draw();
  l5->Draw();
  
  TLegend *leg = new TLegend(0.3,0.7,0.48,0.9);
  leg->AddEntry("hprimary", "Primaries");
  leg->AddEntry("hsecondary", "Secondaries");
  leg->Draw();
/*
   TH1D *hprimary0 = (TH1D*) hprimary->Clone("hprimary0"); //we want to highlight a specific area of the histogram
  hprimary0->GetXaxis()->SetRange(0,hprimary->FindBin(2.8)-1);
  hprimary0->SetFillColor(kYellow);
  hprimary0->Draw("SAME");
  TH1D *hsecondary0 = (TH1D*) hsecondary->Clone("hsecondary0");
  hsecondary0->GetXaxis()->SetRange(0,hsecondary->FindBin(2.8));
  hsecondary0->SetFillColor(kYellow);
  hsecondary0->Draw("SAME");

  TImage *img = TImage::Create();
  cvertex->Update();
  img->FromPad(cvertex);  
  img->WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/ECC0.png");  
  hprimary0->GetXaxis()->SetRange(hprimary->FindBin(2.8), hprimary->FindBin(5.6));
  hsecondary0->GetXaxis()->SetRange(hsecondary->FindBin(2.8), hsecondary->FindBin(5.6));
  cvertex->Draw();
  cvertex->Update();
  img->FromPad(cvertex);
  img->WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/ECC1.png");
  hprimary0->GetXaxis()->SetRange(hprimary->FindBin(5.6)+1, hprimary->FindBin(5.6*2));
  hsecondary0->GetXaxis()->SetRange(hsecondary->FindBin(5.6)+1, hsecondary->FindBin(5.6*2));
  cvertex->Draw();
  cvertex->Update();
  img->FromPad(cvertex);
  img->WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/ECC2.png");
  hprimary0->GetXaxis()->SetRange(hprimary->FindBin(5.6*2)+1, hprimary->FindBin(5.6*3));
  hsecondary0->GetXaxis()->SetRange(hsecondary->FindBin(5.6*2)+1, hsecondary->FindBin(5.6*3));
  cvertex->Draw();
  cvertex->Update();
  img->FromPad(cvertex);
  img->WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/ECC3.png");
  hprimary0->GetXaxis()->SetRange(hprimary->FindBin(5.6*3)+1, hprimary->FindBin(5.6*4));
  hsecondary0->GetXaxis()->SetRange(hsecondary->FindBin(5.6*3)+1, hsecondary->FindBin(5.6*4));
  cvertex->Draw();
  cvertex->Update();
  img->FromPad(cvertex);
  img->WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/ECC4.png");
  hprimary0->GetXaxis()->SetRange(hprimary->FindBin(5.6*4)+1, hprimary->FindBin(5.6*5));
  hsecondary0->GetXaxis()->SetRange(hsecondary->FindBin(5.6*4)+1, hsecondary->FindBin(5.6*5));
  cvertex->Draw();
  cvertex->Update();
  img->FromPad(cvertex);
  img->WriteImage("/home/utente/Lavoro/Analisi/temporanei/vertici_produzione_charm/ECC5.png");
  */
  TCanvas *c = new TCanvas();
  hall->Draw();
  hall->Fit("expo");
  hall->GetXaxis()->SetTitle("cm");


  cout<<"Primari charm da 0 a Run1/2: "<<hprimary->Integral(hprimary->FindBin(0), hprimary->FindBin(2.8)-1)/hprimary->GetEntries()<<endl; //I was 'asked' to put  the first half
  cout<<"Secondari charm da 0 a Run1/2: "<<hsecondary->Integral(hsecondary->FindBin(0), hsecondary->FindBin(2.8)-1)/hsecondary->GetEntries()<<endl; //I was 'asked' to put the first half

  cout<<"Primari charm da 0 a Run2/2: "<<hprimary->Integral(hprimary->FindBin(2.8), hprimary->FindBin(5.6)-1)/hprimary->GetEntries()<<endl; //Second half
  cout<<"Secondari charm da 0 a Run2/2: "<<hsecondary->Integral(hsecondary->FindBin(2.8), hsecondary->FindBin(5.6)-1)/hsecondary->GetEntries()<<endl; //Second half

  cout<<"Primari charm da 0 a Run1: "<<hprimary->Integral(hprimary->FindBin(0),hprimary->FindBin(5.6)-1)/hprimary->GetEntries()<<endl;
  cout<<"Secondari da 0 a Run1: "<<hsecondary->Integral(hsecondary->FindBin(0),hsecondary->FindBin(5.6)-1)/hsecondary->GetEntries()<<endl;
  cout<<"Primari charm da Run1 a Run2: "<<hprimary->Integral(hprimary->FindBin(5.6*1),hprimary->FindBin(5.6*2)-1)/hprimary->GetEntries()<<endl;
  cout<<"Secondari da Run1 a Run2: "<<hsecondary->Integral(hsecondary->FindBin(5.6*1),hsecondary->FindBin(5.6*2)-1)/hsecondary->GetEntries()<<endl;
  cout<<"Primari charm da Run2 a Run3: "<<hprimary->Integral(hprimary->FindBin(5.6*2),hprimary->FindBin(5.6*3)-1)/hprimary->GetEntries()<<endl;
  cout<<"Secondari da Run2 a Run3: "<<hsecondary->Integral(hsecondary->FindBin(5.6*2),hsecondary->FindBin(5.6*3)-1)/hsecondary->GetEntries()<<endl;
  cout<<"Primari charm da Run3 a Run4: "<<hprimary->Integral(hprimary->FindBin(5.6*3),hprimary->FindBin(5.6*4)-1)/hprimary->GetEntries()<<endl;
  cout<<"Secondari da Run3 a Run4: "<<hsecondary->Integral(hsecondary->FindBin(5.6*3),hsecondary->FindBin(5.6*4)-1)/hsecondary->GetEntries()<<endl;
  cout<<"Primari charm da Run4 a Run5: "<<hprimary->Integral(hprimary->FindBin(5.6*4),hprimary->FindBin(5.6*5)-1)/hprimary->GetEntries()<<endl;
  cout<<"Secondari da Run4 a Run5: "<<hsecondary->Integral(hsecondary->FindBin(5.6*4),hsecondary->FindBin(5.6*5)-1)/hsecondary->GetEntries()<<endl;
  cout<<"Primari charm da 0 a Run5:"<<hprimary->Integral(hprimary->FindBin(0), hprimary->FindBin(5.6*5)-1)/hprimary->GetEntries()<<endl;
  cout<<"Primari charm da 0 a Run12:"<<hprimary->Integral(hprimary->FindBin(0), hprimary->FindBin(5.6*12)-1)/hprimary->GetEntries()<<endl;
  cout<<"Secondari da 0 a Run5:"<<hsecondary->Integral(hsecondary->FindBin(0), hsecondary->FindBin(5.6*5)-1)/hsecondary->GetEntries()<<endl;

 cout<<"Primari charm da Run5 a Run6: "<<hprimary->Integral(hprimary->FindBin(5.6*5),hprimary->FindBin(5.6*6)-1)/hprimary->GetEntries()<<endl;
 cout<<"Secondari da Run5 a Run6: "<<hsecondary->Integral(hsecondary->FindBin(5.6*5),hsecondary->FindBin(5.6*6)-1)/hsecondary->GetEntries()<<endl;
 cout<<"Primari charm da Run6 a Run7: "<<hprimary->Integral(hprimary->FindBin(5.6*6),hprimary->FindBin(5.6*7)-1)/hprimary->GetEntries()<<endl;
 cout<<"Secondari da Run6 a Run7: "<<hsecondary->Integral(hsecondary->FindBin(5.6*6),hsecondary->FindBin(5.6*7)-1)/hsecondary->GetEntries()<<endl;
 cout<<"Primari charm da Run7 a Run8: "<<hprimary->Integral(hprimary->FindBin(5.6*7),hprimary->FindBin(5.6*8)-1)/hprimary->GetEntries()<<endl;
 cout<<"Secondari da Run7 a Run8: "<<hsecondary->Integral(hsecondary->FindBin(5.6*7),hsecondary->FindBin(5.6*8)-1)/hsecondary->GetEntries()<<endl;
 cout<<"Primari charm da Run8 a Run9: "<<hprimary->Integral(hprimary->FindBin(5.6*8),hprimary->FindBin(5.6*9)-1)/hprimary->GetEntries()<<endl;
 cout<<"Secondari da Run8 a Run9: "<<hsecondary->Integral(hsecondary->FindBin(5.6*8),hsecondary->FindBin(5.6*9)-1)/hsecondary->GetEntries()<<endl;
 cout<<"Primari charm da Run9 a Run10: "<<hprimary->Integral(hprimary->FindBin(5.6*9),hprimary->FindBin(5.6*10)-1)/hprimary->GetEntries()<<endl;
 cout<<"Secondari da Run9 a Run10: "<<hsecondary->Integral(hsecondary->FindBin(5.6*9),hsecondary->FindBin(5.6*10)-1)/hsecondary->GetEntries()<<endl;
 cout<<"Primari charm da Run10 a Run11: "<<hprimary->Integral(hprimary->FindBin(5.6*10),hprimary->FindBin(5.6*11)-1)/hprimary->GetEntries()<<endl;
 cout<<"Secondari da Run10 a Run11: "<<hsecondary->Integral(hsecondary->FindBin(5.6*10),hsecondary->FindBin(5.6*11)-1)/hsecondary->GetEntries()<<endl;
 cout<<"Primari charm da Run11 a Run12: "<<hprimary->Integral(hprimary->FindBin(5.6*11),hprimary->FindBin(5.6*12)-1)/hprimary->GetEntries()<<endl;
 cout<<"Secondari da Run11 a Run12: "<<hsecondary->Integral(hsecondary->FindBin(5.6*11),hsecondary->FindBin(5.6*12)-1)/hsecondary->GetEntries()<<endl;

  cout<<"Tutti i charm da 0 a Run1/2: "<<hall->Integral(hall->FindBin(0), hall->FindBin(2.8)-1)/hall->GetEntries()<<endl; //I was 'asked' to put the first half
  cout<<"Tutti i charm da 0 a Run2/2: "<<hall->Integral(hall->FindBin(2.8), hall->FindBin(5.6)-1)/hall->GetEntries()<<endl; //second half
  cout<<"Tutti i charm da 0 a Run1: "<<hall->Integral(hall->FindBin(0), hall->FindBin(5.6)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run1 a Run2: "<<hall->Integral(hall->FindBin(5.6*1),hall->FindBin(5.6*2)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run2 a Run3: "<<hall->Integral(hall->FindBin(5.6*2),hall->FindBin(5.6*3)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run3 a Run4: "<<hall->Integral(hall->FindBin(5.6*3),hall->FindBin(5.6*4)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run4 a Run5: "<<hall->Integral(hall->FindBin(5.6*4),hall->FindBin(5.6*5)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da 0 a Run5:"<<hall->Integral(hall->FindBin(0), hall->FindBin(5.6*5)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run5 a Run6: "<<hall->Integral(hall->FindBin(5.6*5),hall->FindBin(5.6*6)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run6 a Run7: "<<hall->Integral(hall->FindBin(5.6*6),hall->FindBin(5.6*7)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run7 a Run8: "<<hall->Integral(hall->FindBin(5.6*7),hall->FindBin(5.6*8)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run8 a Run9: "<<hall->Integral(hall->FindBin(5.6*8),hall->FindBin(5.6*9)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run9 a Run10: "<<hall->Integral(hall->FindBin(5.6*9),hall->FindBin(5.6*10)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run10 a Run11: "<<hall->Integral(hall->FindBin(5.6*10),hall->FindBin(5.6*11)-1)/hall->GetEntries()<<endl;
  cout<<"Tutti i charm da Run11 a Run12: "<<hall->Integral(hall->FindBin(5.6*11),hall->FindBin(5.6*12)-1)/hall->GetEntries()<<endl;
/*  cout<<"Primari da 0 a Run1: "<<hprimary->Integral(hprimary->FindBin(0),hprimary->FindBin(13.99))/hprimary->GetEntries()<<endl;
  cout<<"Secondari da 0 a Run1: "<<hsecondary->Integral(hsecondary->FindBin(0),hsecondary->FindBin(13.99))/hsecondary->GetEntries()<<endl;
  cout<<"Primari da Run1 a Run2: "<<hprimary->Integral(hprimary->FindBin(14.01),hprimary->FindBin(28.99))/hsecondary->GetEntries()<<endl;
  cout<<"Secondari da Run1 a Run2: "<<hsecondary->Integral(hsecondary->FindBin(14.01),hsecondary->FindBin(28.99))/hsecondary->GetEntries()<<endl;
  cout<<"Primari da 0 a Run2:"<<hprimary->Integral(hprimary->FindBin(0), hprimary->FindBin(28.99))/hprimary->GetEntries()<<endl;
  cout<<"Secondari da 0 a Run2:"<<hsecondary->Integral(hsecondary->FindBin(0), hsecondary->FindBin(28.99))/hsecondary->GetEntries()<<endl;*/
 // hprimary->Rebin(10);
 // hsecondary->Rebin(10);
  //hall->Rebin(10);
}
void dafunzione(){
  TF1 *f1= new TF1("f1", "(x + [0] < 58.0) * exp(-x/15.25) + (x + [0] > 58.0) * exp(-x/9.95)", 0, 1000); 

  TH1D * hvz = new TH1D("hvz", "Vertex distribution",300,0,150);
  TH1D * hvz1 = new TH1D("hvz1", "Vertex distribution",300,0,150);
  
  Double_t zminimum = 116.2;
  Double_t checkmin = -14.6;
  Double_t checkmax = 0;
  Double_t zv1;
  Double_t zv2;
  
  ifstream file;
  file.open("secondari200000.txt");
  
  Int_t check[200000];
  
  for (int i = 0; i < 200000; i++){
  file >> check[i];
  }
  
  int ncharm = 0;
  for (int i = 0; i < 200000; i++){
    ncharm++;
    f1->SetParameter(0,0);
    if (i%2 == 0){ //genero 100000 charm primari e 100000 secondari
    //if (check[ncharm-1] == 1){
    //zv1 = speranza->Exp(16) - zminimum; //only the first interaction length
      zv1 = f1->GetRandom(); //only the first interaction length
      zv2 = 0.;
      //cout<<"primario"<<endl;
      //}
      hvz->Fill(zv1+zv2);
    }
    else{
      //zv1 = speranza->Exp(16) - zminimum;
      zv1 = f1->GetRandom(); //only the first interaction length
      f1->SetParameter(0,zv1);
      //zv2 = speranza->Exp(16);
      zv2 = f1->GetRandom();
      //cout<<"secondario"<<endl;
      hvz1->Fill(zv1+zv2);
}
    
    //if ((zv1 + zv2 > checkmin) && (zv1 + zv2 < checkmax)) cout<<zv1<<" "<<zv2<<" "<<zv1+zv2<<endl;
    
    if ((i % 10000) == 0) cout<<i<<endl;
  }
  
  //hprimary.Draw();
  TCanvas *c = new TCanvas();
  hvz->Draw();
  hvz->GetXaxis()->SetTitle("cm");
  //hvz->Scale(1./hvz->Integral());
  hvz1->SetLineColor(kRed);
  hvz1->SetLineWidth(2);
  hvz1->SetLineStyle(2);
  //hvz1->Scale(1./hvz1->Integral());
  cout<<endl;
  cout<<"Integrale totale primari: "<<hvz->Integral()<<endl;
  cout<<"Integrale totale cascata: "<<hvz1->Integral()<<endl;
  //Run1: 1 blocco da 8 cm e 2 blocchi da 2.5 cm. Run2: 5 blocchi da 2.5 cm di TZM (a cui nella realtà va sottratto 0.1 cm di tantalio da ambo i lati)
  cout<<"Primari da 0 a Run1: "<<hvz->Integral(hvz->FindBin(0),hvz->FindBin(13.0))<<endl;
  cout<<"Secondari da 0 a Run1: "<<hvz1->Integral(hvz1->FindBin(0),hvz1->FindBin(13.0))<<endl;
  cout<<"Primari da Run1 a Run2: "<<hvz->Integral(hvz->FindBin(13.0),hvz->FindBin(25.5))<<endl;
  cout<<"Secondari da Run1 a Run2: "<<hvz1->Integral(hvz1->FindBin(13.0),hvz1->FindBin(25.5))<<endl;
  cout<<"Frazione primari nei primi 2 run: "<<(double) (hvz->Integral(hvz->FindBin(0), hvz->FindBin(25.5)))/(hvz->Integral())<<endl;
  cout<<"Frazione secondari nei primi 2 run: "<<(double) (hvz1->Integral(hvz1->FindBin(0),hvz1->FindBin(25.5)))/(hvz1->Integral())<<endl;
  hvz1->Draw("SAMES");
  Double_t ymax = hvz->GetMaximum();
  TLine *l1 = new TLine(13.0,0,13.0,ymax);
  TLine *l2 = new TLine(25.5,0,25.5,ymax);
  TLegend *leg = new TLegend(0.1,0.7,0.48,0.9);
  l1->Draw();
  l2->Draw();
  leg->AddEntry("hvz", "Primaries");
  leg->AddEntry("hvz1", "Secondaries");
  leg->Draw();
  c->Print("./cascade/vertex_distributions.root");
  c->Print("./cascade/vertex_distributions.png");
  }
//11.5 + 12.4
