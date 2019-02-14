//macro realizzata per stimare la velocità con cui dovrà viaggiare il piano per avere una distribuzione spaziale uniforme. (creata il 15 Marzo 2017)

void velocity_estimator(){
  const Double_t vymin = 2.0;
  const Double_t vymax = 2.0;
  const Double_t xstepmin = 1.05;
  const Double_t xstepmax = 1.05;
  
  bool ciclox = false;
  bool cicloy = false;

  const Double_t xlength = 12.5; //in cm
  const Double_t ylength = 10.;
  const Double_t tspill = 4.; //in secondi 

  TRandom3 *rangen = new TRandom3();

  //parametri distribuzione spaziale originaria
  const Double_t mean = 0.;
  const Double_t sigmax = 1./2.35;
  const Double_t sigmay = 1./2.35;
  
  //const Int_t nevents = 1000000 * 10; //numero eventi da generare per spill (1 mil)
  const Int_t nevents = 1e+5; 
  //  const Int_t nevents = 100;
  Double_t x; //posizione iniziale (distribuzione gaussiana)
  Double_t y;
  Double_t t; //istante di tempo

  Double_t xstep;
  Double_t vy;
  //Double_t vy = 0.2.;
  
  //Double_t x1[nevents * 1000]; //posizione finale
  //Double_t y1[nevents * 1000];
  Double_t x1;
  Double_t y1;
  TH1D *hx = new TH1D("hx", "Istogramma posizioni iniziali", 100, -3, 3);
  TH1D *ht = new TH1D("ht", "Istogramma istanti di tempo", 100, 0, 5);
  TH1D *hx1 = new TH1D("hx1", "Istogramma posizioni finali", 100, -5, 5);
  
  // TH2D *hxy1 = new TH2D("hxy1", "Istogramma posizioni finali nel piano xy", 100, -5, 5, 10, 0, 1); //solo il primo cm (per test veloci)
  TH2D *hxy1 = new TH2D("hxy1", "Istogramma posizioni finali nel piano xy", 100, -6.25, 6.25, 100, -5, 5); //tutto il bersaglio 
  //TH2D *hxy1 = new TH2D("hxy1", "Istogramma posizioni finali nel piano xy", 80, -4, 4, 80, -4, 4); //tutto il bersaglio bordi esclusi
  TH1D *hdist = new TH1D("hdist", "Istogramma distanze", 100, 0, 100);
  
  Double_t offsetx,offsety;
  Int_t tmax;
  Int_t contatore;
  Int_t nspill;
  Int_t nspilltot, neventsperspill; //numero di spill effettuate
  
  Int_t numero;
  Int_t numeromedio = 0;
  Int_t nbins = 0; //numero medio di particelle nei bin diversi da zero, di numero nbins.
  Double_t sommaquadra = 0;

  TGraph2D *numbergraph = new TGraph2D();
  TGraph2D *chigraph = new TGraph2D();
  TGraph *graph = new TGraph();
  //TGraph2D *graph = new TGraph2D();
  TGraph *graph1 = new TGraph();
  Int_t npoints = 0; //tiene conto dei punti nel grafico
  Int_t npoints1 = 0;

  //const Double_t xstep = 0.3;
  //Double_t tstep;
  xstep = xstepmin;
  while (xstep < xstepmax + 0.1){
    vy = vymin;
    //tstep = xstep/vy;
    while (vy < vymax + 0.1){
      nspill = 0;
      nspilltot = 2 * TMath::Ceil((ylength / (2 * xstep))) * TMath::Ceil((xlength / (vy*tspill)));
      neventsperspill = (Int_t) nevents/(nspilltot);
      contatore = 0;
      hx->Reset();
      ht->Reset();
      hx1->Reset();
      hxy1->Reset();
      hdist->Reset();
      nbins = 0;
      numeromedio = 0;
      sommaquadra = 0.;
      offsetx = -xlength/2. + 0.5;
      
      while (offsetx < 6.){ //includere il bordo
	cout<<"Il percorso parte da x = "<<offsetx<<endl;
	offsety = -xlength/2.;
	tmax = 0;
	while (vy * tmax < xlength){ //voglio coprire tutto il bersaglio, anche se l'ultima spill è in parte 'sprecata',quando vy * 5 non è un divisore di 10
	  cout<<"OFFSET:"<<offsetx<<endl;
	  nspill++;
	  for (int i = 0; i < neventsperspill; i++){
	    x = rangen->Gaus(mean,sigmax);
	    y = rangen->Gaus(mean,sigmay); //rimane ferma in ogni spill
	    t = rangen->Uniform(0,tspill);
	    //t = (Double_t) rangen->Integer((tspill/tstep)) * tstep;
	    //cout<<t<<endl;
	    y1 = y + vy * t + offsety; //in uno spill non deve coprire tutta la lunghezza	
	    x1 = x + offsetx + (nspill % 2)  * xstep;	  
	    
	    hx1->Fill(x1);
	    hxy1->Fill(x1,y1);
	    contatore++;
	  }
      tmax += tspill;
      offsety += vy * tspill;
      
	}
	
	tmax = 0;
	
	while (vy * tmax < xlength){ 
	  cout<<"OFFSET:"<<offsetx<<endl;
	  nspill++;
	  for (int i = 0; i < neventsperspill; i++){
	    x = rangen->Gaus(mean,sigmax);
	    y = rangen->Gaus(mean,sigmay); //rimane ferma in ogni spill
	    t = rangen->Uniform(0,tspill);
	   // t = (Double_t) rangen->Integer((tspill/tstep)) * tstep;
	    
	    y1 = y - vy * t + offsety; //in uno spill non deve coprire tutta la lunghezza
	    x1 = x + offsetx + (nspill % 2)  * xstep;	   
	    // cout<<t<<" "<<vy*t<<" "<<x1<<" "<<y1<<endl;
	    hx1->Fill(x1);
	    hxy1->Fill(x1,y1);
	    
	    contatore++;
	  }
	  tmax += tspill;
	  offsety -= vy * tspill;
	  
	}
	offsetx += 2 * xstep;
      }
      
	if ((ciclox == false) && (cicloy == false)){
       	TCanvas *c0 = new TCanvas();
	hx->Draw();
	
	TCanvas *c = new TCanvas();
	ht->Draw();
	
	TCanvas *c1 = new TCanvas();
	hx1->Draw();
	TCanvas *c2 = new TCanvas();
	hxy1->GetXaxis()->SetTitle("cm");
	hxy1->GetYaxis()->SetTitle("cm");
	hxy1->Draw("COLZ");
	}
      for (int nx = 0; nx < hxy1->GetNbinsX(); nx++){
	for (int ny = 0; ny < hxy1->GetNbinsY(); ny++){
	  if(hxy1->GetBinContent(nx+1,ny+1) > 0){//bin non nulli
	    numero = hxy1->GetBinContent(nx+1,ny+1);
	    numeromedio += numero;
	    nbins++;
	  }
	}
  }
      numeromedio = Double_t (numeromedio/nbins);
      //calcolo il 'chi quadro', dato dalla somma (N-<N>)**2/sigma**2, con sigma = sqrt(N)
      for (int nx = 0; nx < hxy1->GetNbinsX(); nx++){
	for (int ny = 0; ny < hxy1->GetNbinsY(); ny++){
	  if(hxy1->GetBinContent(nx+1,ny+1)> 0){//bin non nulli
	    numero = hxy1->GetBinContent(nx+1,ny+1);
	    sommaquadra += pow((numero - numeromedio),2)/numero;
	  }
	}
      }
      
      cout<<"Numero bins non nulli:"<<nbins<<" vy = "<<vy<<" xstep = "<<xstep<<endl;
      cout<<"Numero medio: "<<numeromedio<<" chi quadrato "<<sommaquadra<<" chi quadrato diviso gradi di libertà "<<sommaquadra/(nbins-1)<<endl;
      cout<<"Numero spill: "<<nspill<<" "<<nspilltot<<endl;
      cout<<"Numero protoni integrati per spill "<<neventsperspill<<" mentre in totale "<<neventsperspill * nspill<< endl;
      npoints++;
      numbergraph->SetPoint(npoints, vy, xstep, nspill);
      chigraph->SetPoint(npoints, vy, xstep, sommaquadra/(nbins-1));
      if (ciclox == true) vy = vy + 0.1;
      else break;
    }//fine primo ciclo while
   
    if (cicloy == true) {
      xstep = xstep + 0.1;
    }
    else break;
    graph->SetPoint(npoints,vy, sommaquadra/(nbins-1));   
  }//fine secondo ciclo while
  graph->SetMarkerStyle(7);
  TCanvas *cgraph = new TCanvas();
  //graph->Draw("AP");
  graph->GetXaxis()->SetTitle("v[cm/s]");
  graph->GetYaxis()->SetTitle("#chi ^{2}/(ndf)");
  TCanvas *cnumbergraph = new TCanvas();
  numbergraph->GetXaxis()->SetTitle("v[cm/s]");
  numbergraph->GetXaxis()->SetRangeUser(vymin,vymax);
  numbergraph->GetYaxis()->SetRangeUser(xstepmin,xstepmax);
  //numbergraph->Draw("surf1");
  //TCanvas *cchigraph = new TCanvas();
  //chigraph->Draw("surf1");
}
