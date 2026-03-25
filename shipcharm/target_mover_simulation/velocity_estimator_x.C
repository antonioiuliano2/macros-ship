//macro realizzata per stimare la velocità con cui dovrà viaggiare il piano per avere una distribuzione spaziale uniforme. (creata il 15 Marzo 2017)

void velocity_estimator_x(){
  const Double_t vxmin = 2.6;
  const Double_t vxmax = 2.6;
  const Double_t ystepmin = 1.0;
  const Double_t ystepmax = 1.0;
  
  bool ciclox = false;
  bool cicloy = false;

  const Double_t xlength = 12.5; //in cm
  const Double_t ylength = 10.;
  const Double_t tspill = 4.8; //in secondi 

  TRandom3 *rangen = new TRandom3();

  //parametri distribuzione spaziale originaria
  const Double_t mean = 0.;
  const Double_t sigmax = 1/2.35;
  const Double_t sigmay = 1./2.35;
  
  //const Int_t nevents = 1000000 * 10; //numero eventi da generare per spill (1 mil)
  const Int_t nevents = 1e+5; 
  //  const Int_t nevents = 100;
  Double_t x; //posizione iniziale (distribuzione gaussiana)
  Double_t y;
  Double_t t; //istante di tempo

  Double_t ystep;
  Double_t vx;
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
  Double_t tmax;
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
  ystep = ystepmin;
  while (ystep < ystepmax + 0.1){
    vx = vxmin;
    //tstep = xstep/vx;
    while (vx < vxmax + 0.1){
      nspill = 0;
      nspilltot = 2 * TMath::Ceil((ylength / (2 * ystep))) * TMath::Ceil((xlength / (vx*tspill)));
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
      offsety = -ylength/2.;
      
      while (offsety < 5.){ //includere il bordo
	cout<<"Il percorso parte da y = "<<offsety<<endl;
	offsetx = -xlength/2.;
	tmax = 0;
	while (vx * tmax < xlength){ //voglio coprire tutto il bersaglio, anche se l'ultima spill è in parte 'sprecata',quando vx * 5 non è un divisore di 10
	  cout<<"OFFSET:"<<offsetx<<endl;
	  nspill++;
	  for (int i = 0; i < neventsperspill; i++){
	    x = rangen->Gaus(mean,sigmax);
	    y = rangen->Gaus(mean,sigmay); //rimane ferma in ogni spill
	    t = rangen->Uniform(0,tspill);
	    //t = (Double_t) rangen->Integer((tspill/tstep)) * tstep;
	    //cout<<t<<endl;
	    x1 = x + vx * t + offsetx; //in uno spill non deve coprire tutta la lunghezza	
	    y1 = y + offsety + (nspill % 2)  * ystep;	  
	    
	    hx1->Fill(x1);
	    hxy1->Fill(x1,y1);
	    contatore++;
	  }
      tmax += tspill;
      cout<<tmax<<" "<<vx*tmax<<endl;
      offsetx += vx * tspill;
      
	}
	
	tmax = 0;
	
	while (vx * tmax < xlength){ 
	  cout<<"OFFSET:"<<offsetx<<endl;
	  nspill++;
	  for (int i = 0; i < neventsperspill; i++){
	    x = rangen->Gaus(mean,sigmax);
	    y = rangen->Gaus(mean,sigmay); //rimane ferma in ogni spill
	    t = rangen->Uniform(0,tspill);
	   // t = (Double_t) rangen->Integer((tspill/tstep)) * tstep;
	    
	    x1 = x - vx * t + offsetx; //in uno spill non deve coprire tutta la lunghezza
	    y1 = y + offsety + (nspill % 2)  * ystep;	   
	    // cout<<t<<" "<<vx*t<<" "<<x1<<" "<<y1<<endl;
	    hx1->Fill(x1);
	    hxy1->Fill(x1,y1);
	    
	    contatore++;
	  }
	  tmax += tspill;
	  offsetx -= vx * tspill;
	  
	}
	offsety += 2 * ystep;
        cout<<"PROVA:"<<offsety<<endl;
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
      
      cout<<"Numero bins non nulli:"<<nbins<<" vx = "<<vx<<" ystep = "<<ystep<<endl;
      cout<<"Numero medio: "<<numeromedio<<" chi quadrato "<<sommaquadra<<" chi quadrato diviso gradi di libertà "<<sommaquadra/(nbins-1)<<endl;
      cout<<"Numero spill: "<<nspill<<" "<<nspilltot<<endl;
      cout<<"Numero protoni integrati per spill "<<neventsperspill<<" mentre in totale "<<neventsperspill * nspill<< endl;
      npoints++;
      numbergraph->SetPoint(npoints, vx, ystep, nspill);
      chigraph->SetPoint(npoints, vx, ystep, sommaquadra/(nbins-1));
      if (ciclox == true) vx = vx + 0.1;
      else break;
    }//fine primo ciclo while
   
    if (cicloy == true) {
      ystep = ystep + 0.1;
    }
    else break;
    graph->SetPoint(npoints,vx, sommaquadra/(nbins-1));   
  }//fine secondo ciclo while
  graph->SetMarkerStyle(7);
  TCanvas *cgraph = new TCanvas();
  //graph->Draw("AP");
  graph->GetXaxis()->SetTitle("v[cm/s]");
  graph->GetYaxis()->SetTitle("#chi ^{2}/(ndf)");
  TCanvas *cnumbergraph = new TCanvas();
  numbergraph->GetXaxis()->SetTitle("v[cm/s]");
  numbergraph->GetXaxis()->SetRangeUser(vxmin,vxmax);
  numbergraph->GetYaxis()->SetRangeUser(ystepmin,ystepmax);
  //numbergraph->Draw("surf1");
  //TCanvas *cchigraph = new TCanvas();
  //chigraph->Draw("surf1");
}
