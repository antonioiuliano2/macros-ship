// // // ------------------------->
// // // // ----------------------->
// // // ------------------------->
// // --------------------------->
// ----------------------------->

void vz_plot(){

   TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   dir.ReplaceAll("vz_plot.C","");
   dir.ReplaceAll("/./","/");
   ifstream in;
   in.open(Form("%serrors.dat",dir.Data()));

   Float_t x,y,z,u,v;
   Int_t nlines = 0;
   Float_t err_MC[50]={};
   Float_t err_MC_pr[50]={};
   Float_t err_MC_hd[50]={};
   Float_t err_DT[50]={};

   while (1) {
     in >> x >> y >> z >> u >> v;
     err_DT[nlines]=y;
     err_MC[nlines]=z;
     err_MC_pr[nlines]=u;
     err_MC_hd[nlines]=v;
     if (!in.good()) break;
     if (nlines < 5) printf("x=%8f, y=%8f, z=%8f, u=%8f, v=%8f\n",x,y,z,u,v);
     nlines++;
   }
   printf(" found %d points\n",nlines);

   in.close();

  TFile *f = new TFile("full_vz.root","READ"); //simulazione
  TFile *fb = new TFile("with_bkg_full_vz.root","READ"); //dati

  TH1F* hvz1 = (TH1F*)f->Get("vz_ch1");
  TH1F* hvz2 = (TH1F*)f->Get("vz_ch2");
  TH1F* hvz3 = (TH1F*)f->Get("vz_ch3");
  TH1F* hvz4 = (TH1F*)f->Get("vz_ch4");
  TH1F* hvz5 = (TH1F*)f->Get("vz_ch5");
  TH1F* hvz6 = (TH1F*)f->Get("vz_ch6");

  TH1F* hvz1MCpr = (TH1F*)f->Get("vz_ch1_pr");
  TH1F* hvz2MCpr = (TH1F*)f->Get("vz_ch2_pr");
  TH1F* hvz3MCpr = (TH1F*)f->Get("vz_ch3_pr");
  TH1F* hvz4MCpr = (TH1F*)f->Get("vz_ch4_pr");
  TH1F* hvz5MCpr = (TH1F*)f->Get("vz_ch5_pr");
  TH1F* hvz6MCpr = (TH1F*)f->Get("vz_ch6_pr");

  TH1F* hvz1MChd = (TH1F*)f->Get("vz_ch1_hd");
  TH1F* hvz2MChd = (TH1F*)f->Get("vz_ch2_hd");
  TH1F* hvz3MChd = (TH1F*)f->Get("vz_ch3_hd");
  TH1F* hvz4MChd = (TH1F*)f->Get("vz_ch4_hd");
  TH1F* hvz5MChd = (TH1F*)f->Get("vz_ch5_hd");
  TH1F* hvz6MChd = (TH1F*)f->Get("vz_ch6_hd");

  
  TH1F* hvz1wb = (TH1F*)fb->Get("vz_ch1");
  TH1F* hvz2wb = (TH1F*)fb->Get("vz_ch2");
  TH1F* hvz3wb = (TH1F*)fb->Get("vz_ch3");
  TH1F* hvz4wb = (TH1F*)fb->Get("vz_ch4");
  TH1F* hvz5wb = (TH1F*)fb->Get("vz_ch5");
  TH1F* hvz6wb = (TH1F*)fb->Get("vz_ch6");

  ofstream log_file("vz_global.txt");

  //hvz2MCpr->Scale(4); // only 1 quarter
  //hvz2MChd->Scale(4); // only 1 quarter
  
  hvz1->Scale(0.25);
  hvz2->Scale(0.33);
  hvz3->Scale(3.86);
  hvz4->Scale(4.91);
  hvz5->Scale(3.75);
  hvz6->Scale(6.0);

  hvz1->Rebin(3);
  hvz2->Rebin(3);
  hvz3->Rebin(3);
  hvz4->Rebin(3);
  hvz5->Rebin(3);
  hvz6->Rebin(3); 
  // rapportare il numero di protoni alla configurazione CH1
  hvz1MCpr->Scale(1);    // C_mc = MCPOT_CH1 / MCPOT_CHX    135000
  hvz2MCpr->Scale(1);                              //135000
  hvz3MCpr->Scale(3.86);                           //40000
  hvz4MCpr->Scale(4.91);
  hvz5MCpr->Scale(3.75);           // err MC = C_mc * rad(N)
  hvz6MCpr->Scale(6.0);             
  
  hvz1MCpr->Rebin(3); //5 bin
  hvz2MCpr->Rebin(3); //5 bin
  hvz3MCpr->Rebin(3); //10 bin
  hvz4MCpr->Rebin(3); //10 bin
  hvz5MCpr->Rebin(3); //10 bin
  hvz6MCpr->Rebin(3); //10 bin
  
  hvz1MChd->Scale(1);
  hvz2MChd->Scale(1);
  hvz3MChd->Scale(3.86);
  hvz4MChd->Scale(4.91);
  hvz5MChd->Scale(3.75);
  hvz6MChd->Scale(6.0);
  
  
  hvz1MChd->Rebin(3);
  hvz2MChd->Rebin(3);
  hvz3MChd->Rebin(3);
  hvz4MChd->Rebin(3);
  hvz5MChd->Rebin(3);
  hvz6MChd->Rebin(3);

  hvz1wb->Scale(0.25); //abbiamo solo 4 run di CH1 e 3 run di CH2
  hvz2wb->Scale(0.33);
  hvz3wb->Scale(3.86); //abbiamo 1 run per CH3, CH4, CH5, CH6
  hvz4wb->Scale(4.91);
  hvz5wb->Scale(3.75);
  hvz6wb->Scale(6.0);

  // AGGIUNGERE SCALE FACTOR
  

  hvz1wb->Rebin(3);
  hvz2wb->Rebin(3);
  hvz3wb->Rebin(3);
  hvz4wb->Rebin(3);
  hvz5wb->Rebin(3);
  hvz6wb->Rebin(3);

  TH1F *hvz1MCfull = new TH1F("hMC1","hMC1",5,-35,-5);
  TH1F *hvz2MCfull = new TH1F("hMC2","hMC2",5,-35,-5);
  TH1F *hvz3MCfull = new TH1F("hMC3","hMC3",10,-65,-5);
  TH1F *hvz4MCfull = new TH1F("hMC4","hMC4",10,-65,-5);
  TH1F *hvz5MCfull = new TH1F("hMC5","hMC5",10,-65,-5);
  TH1F *hvz6MCfull = new TH1F("hMC6","hMC6",10,-65,-5);

  hvz1MCfull->Add(hvz1MCpr);
  hvz1MCfull->Add(hvz1MChd);
  hvz2MCfull->Add(hvz2MCpr);
  hvz2MCfull->Add(hvz2MChd);
  hvz3MCfull->Add(hvz3MCpr);
  hvz3MCfull->Add(hvz3MChd);
  hvz4MCfull->Add(hvz4MCpr);
  hvz4MCfull->Add(hvz4MChd);
  hvz5MCfull->Add(hvz5MCpr);
  hvz5MCfull->Add(hvz5MChd);
  hvz6MCfull->Add(hvz6MCpr);
  hvz6MCfull->Add(hvz6MChd);


  //graph calcolato a mano, tenendo conto delle posizioni dei bin (il rebin)
  TGraphErrors *grMCpr = new TGraphErrors();

  grMCpr->SetPoint(0,3,hvz1MCpr->GetBinContent(1));
  grMCpr->SetPoint(1,9,hvz1MCpr->GetBinContent(2));
  grMCpr->SetPoint(2,15,hvz1MCpr->GetBinContent(3));
  grMCpr->SetPoint(3,21,hvz1MCpr->GetBinContent(4));
  grMCpr->SetPoint(4,27,hvz1MCpr->GetBinContent(5));

  grMCpr->SetPoint(5,38,hvz2MCpr->GetBinContent(1));
  grMCpr->SetPoint(6,44,hvz2MCpr->GetBinContent(2));
  grMCpr->SetPoint(7,50,hvz2MCpr->GetBinContent(3));
  grMCpr->SetPoint(8,56,hvz2MCpr->GetBinContent(4));
  grMCpr->SetPoint(9,62,hvz2MCpr->GetBinContent(5));

  grMCpr->SetPoint(10,83,hvz3MCpr->GetBinContent(1));
  grMCpr->SetPoint(11,89,hvz3MCpr->GetBinContent(2));
  grMCpr->SetPoint(12,95,hvz3MCpr->GetBinContent(3));
  grMCpr->SetPoint(13,101,hvz3MCpr->GetBinContent(4));
  grMCpr->SetPoint(14,107,hvz3MCpr->GetBinContent(5));
  grMCpr->SetPoint(15,113,hvz3MCpr->GetBinContent(6));
  grMCpr->SetPoint(16,119,hvz3MCpr->GetBinContent(7));
  grMCpr->SetPoint(17,125,hvz3MCpr->GetBinContent(8));
  grMCpr->SetPoint(18,131,hvz3MCpr->GetBinContent(9));
  grMCpr->SetPoint(19,137,hvz3MCpr->GetBinContent(10));

  grMCpr->SetPoint(20,158,hvz4MCpr->GetBinContent(1));
  grMCpr->SetPoint(21,164,hvz4MCpr->GetBinContent(2));
  grMCpr->SetPoint(22,170,hvz4MCpr->GetBinContent(3));
  grMCpr->SetPoint(23,176,hvz4MCpr->GetBinContent(4));
  grMCpr->SetPoint(24,182,hvz4MCpr->GetBinContent(5));
  grMCpr->SetPoint(25,188,hvz4MCpr->GetBinContent(6));
  grMCpr->SetPoint(26,194,hvz4MCpr->GetBinContent(7));
  grMCpr->SetPoint(27,200,hvz4MCpr->GetBinContent(8));
  grMCpr->SetPoint(28,206,hvz4MCpr->GetBinContent(9));
  grMCpr->SetPoint(29,212,hvz4MCpr->GetBinContent(10));

  grMCpr->SetPoint(30,233,hvz5MCpr->GetBinContent(1));
  grMCpr->SetPoint(31,239,hvz5MCpr->GetBinContent(2));
  grMCpr->SetPoint(32,245,hvz5MCpr->GetBinContent(3));
  grMCpr->SetPoint(33,251,hvz5MCpr->GetBinContent(4));
  grMCpr->SetPoint(34,257,hvz5MCpr->GetBinContent(5));
  grMCpr->SetPoint(35,263,hvz5MCpr->GetBinContent(6));
  grMCpr->SetPoint(36,269,hvz5MCpr->GetBinContent(7));
  grMCpr->SetPoint(37,275,hvz5MCpr->GetBinContent(8));
  grMCpr->SetPoint(38,281,hvz5MCpr->GetBinContent(9));
  grMCpr->SetPoint(39,287,hvz5MCpr->GetBinContent(10));
  
  grMCpr->SetPoint(40,308,hvz6MCpr->GetBinContent(1));
  grMCpr->SetPoint(41,314,hvz6MCpr->GetBinContent(2));
  grMCpr->SetPoint(42,320,hvz6MCpr->GetBinContent(3));
  grMCpr->SetPoint(43,326,hvz6MCpr->GetBinContent(4));
  grMCpr->SetPoint(44,332,hvz6MCpr->GetBinContent(5));
  grMCpr->SetPoint(45,338,hvz6MCpr->GetBinContent(6));
  grMCpr->SetPoint(46,344,hvz6MCpr->GetBinContent(7));
  grMCpr->SetPoint(47,350,hvz6MCpr->GetBinContent(8));
  grMCpr->SetPoint(48,356,hvz6MCpr->GetBinContent(9));
  grMCpr->SetPoint(49,362,hvz6MCpr->GetBinContent(10));

  // ERRORS (letti da un altro file)
  for(int i=0;i<50;i++){
    grMCpr->SetPointError(i,0,err_MC_pr[i]);
  }

  TGraphErrors *grMChd = new TGraphErrors();

  grMChd->SetPoint(0,3,hvz1MChd->GetBinContent(1));
  grMChd->SetPoint(1,9,hvz1MChd->GetBinContent(2));
  grMChd->SetPoint(2,15,hvz1MChd->GetBinContent(3));
  grMChd->SetPoint(3,21,hvz1MChd->GetBinContent(4));
  grMChd->SetPoint(4,27,hvz1MChd->GetBinContent(5));

  grMChd->SetPoint(5,38,hvz2MChd->GetBinContent(1));
  grMChd->SetPoint(6,44,hvz2MChd->GetBinContent(2));
  grMChd->SetPoint(7,50,hvz2MChd->GetBinContent(3));
  grMChd->SetPoint(8,56,hvz2MChd->GetBinContent(4));
  grMChd->SetPoint(9,62,hvz2MChd->GetBinContent(5));

  grMChd->SetPoint(10,83,hvz3MChd->GetBinContent(1));
  grMChd->SetPoint(11,89,hvz3MChd->GetBinContent(2));
  grMChd->SetPoint(12,95,hvz3MChd->GetBinContent(3));
  grMChd->SetPoint(13,101,hvz3MChd->GetBinContent(4));
  grMChd->SetPoint(14,107,hvz3MChd->GetBinContent(5));
  grMChd->SetPoint(15,113,hvz3MChd->GetBinContent(6));
  grMChd->SetPoint(16,119,hvz3MChd->GetBinContent(7));
  grMChd->SetPoint(17,125,hvz3MChd->GetBinContent(8));
  grMChd->SetPoint(18,131,hvz3MChd->GetBinContent(9));
  grMChd->SetPoint(19,137,hvz3MChd->GetBinContent(10));

  grMChd->SetPoint(20,158,hvz4MChd->GetBinContent(1));
  grMChd->SetPoint(21,164,hvz4MChd->GetBinContent(2));
  grMChd->SetPoint(22,170,hvz4MChd->GetBinContent(3));
  grMChd->SetPoint(23,176,hvz4MChd->GetBinContent(4));
  grMChd->SetPoint(24,182,hvz4MChd->GetBinContent(5));
  grMChd->SetPoint(25,188,hvz4MChd->GetBinContent(6));
  grMChd->SetPoint(26,194,hvz4MChd->GetBinContent(7));
  grMChd->SetPoint(27,200,hvz4MChd->GetBinContent(8));
  grMChd->SetPoint(28,206,hvz4MChd->GetBinContent(9));
  grMChd->SetPoint(29,212,hvz4MChd->GetBinContent(10));

  grMChd->SetPoint(30,233,hvz5MChd->GetBinContent(1));
  grMChd->SetPoint(31,239,hvz5MChd->GetBinContent(2));
  grMChd->SetPoint(32,245,hvz5MChd->GetBinContent(3));
  grMChd->SetPoint(33,251,hvz5MChd->GetBinContent(4));
  grMChd->SetPoint(34,257,hvz5MChd->GetBinContent(5));
  grMChd->SetPoint(35,263,hvz5MChd->GetBinContent(6));
  grMChd->SetPoint(36,269,hvz5MChd->GetBinContent(7));
  grMChd->SetPoint(37,275,hvz5MChd->GetBinContent(8));
  grMChd->SetPoint(38,281,hvz5MChd->GetBinContent(9));
  grMChd->SetPoint(39,287,hvz5MChd->GetBinContent(10));
  
  grMChd->SetPoint(40,308,hvz6MChd->GetBinContent(1));
  grMChd->SetPoint(41,314,hvz6MChd->GetBinContent(2));
  grMChd->SetPoint(42,320,hvz6MChd->GetBinContent(3));
  grMChd->SetPoint(43,326,hvz6MChd->GetBinContent(4));
  grMChd->SetPoint(44,332,hvz6MChd->GetBinContent(5));
  grMChd->SetPoint(45,338,hvz6MChd->GetBinContent(6));
  grMChd->SetPoint(46,344,hvz6MChd->GetBinContent(7));
  grMChd->SetPoint(47,350,hvz6MChd->GetBinContent(8));
  grMChd->SetPoint(48,356,hvz6MChd->GetBinContent(9));
  grMChd->SetPoint(49,362,hvz6MChd->GetBinContent(10));

   // ERRORS

  for(int i=0;i<50;i++){
    grMChd->SetPointError(i,0,err_MC_hd[i]);
  }

  TGraphErrors *grMC = new TGraphErrors();

  grMC->SetPoint(0,3,hvz1MCfull->GetBinContent(1));
  grMC->SetPoint(1,9,hvz1MCfull->GetBinContent(2));
  grMC->SetPoint(2,15,hvz1MCfull->GetBinContent(3));
  grMC->SetPoint(3,21,hvz1MCfull->GetBinContent(4));
  grMC->SetPoint(4,27,hvz1MCfull->GetBinContent(5));

  grMC->SetPoint(5,38,hvz2MCfull->GetBinContent(1));
  grMC->SetPoint(6,44,hvz2MCfull->GetBinContent(2));
  grMC->SetPoint(7,50,hvz2MCfull->GetBinContent(3));
  grMC->SetPoint(8,56,hvz2MCfull->GetBinContent(4));
  grMC->SetPoint(9,62,hvz2MCfull->GetBinContent(5));

  grMC->SetPoint(10,83,hvz3MCfull->GetBinContent(1));
  grMC->SetPoint(11,89,hvz3MCfull->GetBinContent(2));
  grMC->SetPoint(12,95,hvz3MCfull->GetBinContent(3));
  grMC->SetPoint(13,101,hvz3MCfull->GetBinContent(4));
  grMC->SetPoint(14,107,hvz3MCfull->GetBinContent(5));
  grMC->SetPoint(15,113,hvz3MCfull->GetBinContent(6));
  grMC->SetPoint(16,119,hvz3MCfull->GetBinContent(7));
  grMC->SetPoint(17,125,hvz3MCfull->GetBinContent(8));
  grMC->SetPoint(18,131,hvz3MCfull->GetBinContent(9));
  grMC->SetPoint(19,137,hvz3MCfull->GetBinContent(10));

  grMC->SetPoint(20,158,hvz4MCfull->GetBinContent(1));
  grMC->SetPoint(21,164,hvz4MCfull->GetBinContent(2));
  grMC->SetPoint(22,170,hvz4MCfull->GetBinContent(3));
  grMC->SetPoint(23,176,hvz4MCfull->GetBinContent(4));
  grMC->SetPoint(24,182,hvz4MCfull->GetBinContent(5));
  grMC->SetPoint(25,188,hvz4MCfull->GetBinContent(6));
  grMC->SetPoint(26,194,hvz4MCfull->GetBinContent(7));
  grMC->SetPoint(27,200,hvz4MCfull->GetBinContent(8));
  grMC->SetPoint(28,206,hvz4MCfull->GetBinContent(9));
  grMC->SetPoint(29,212,hvz4MCfull->GetBinContent(10));

  grMC->SetPoint(30,233,hvz5MCfull->GetBinContent(1));
  grMC->SetPoint(31,239,hvz5MCfull->GetBinContent(2));
  grMC->SetPoint(32,245,hvz5MCfull->GetBinContent(3));
  grMC->SetPoint(33,251,hvz5MCfull->GetBinContent(4));
  grMC->SetPoint(34,257,hvz5MCfull->GetBinContent(5));
  grMC->SetPoint(35,263,hvz5MCfull->GetBinContent(6));
  grMC->SetPoint(36,269,hvz5MCfull->GetBinContent(7));
  grMC->SetPoint(37,275,hvz5MCfull->GetBinContent(8));
  grMC->SetPoint(38,281,hvz5MCfull->GetBinContent(9));
  grMC->SetPoint(39,287,hvz5MCfull->GetBinContent(10));
  
  grMC->SetPoint(40,308,hvz6MCfull->GetBinContent(1));
  grMC->SetPoint(41,314,hvz6MCfull->GetBinContent(2));
  grMC->SetPoint(42,320,hvz6MCfull->GetBinContent(3));
  grMC->SetPoint(43,326,hvz6MCfull->GetBinContent(4));
  grMC->SetPoint(44,332,hvz6MCfull->GetBinContent(5));
  grMC->SetPoint(45,338,hvz6MCfull->GetBinContent(6));
  grMC->SetPoint(46,344,hvz6MCfull->GetBinContent(7));
  grMC->SetPoint(47,350,hvz6MCfull->GetBinContent(8));
  grMC->SetPoint(48,356,hvz6MCfull->GetBinContent(9));
  grMC->SetPoint(49,362,hvz6MCfull->GetBinContent(10));

  // ERRORS

  for(int i=0;i<50;i++){
    grMC->SetPointError(i,0,err_MC[i]);
  }
  /*
  grMC->SetPointError(0,0,51);
  grMC->SetPointError(1,0,52);
  grMC->SetPointError(2,0,51);
  grMC->SetPointError(3,0,53);
  grMC->SetPointError(4,0,52);

  grMC->SetPointError(5,0,193);
  grMC->SetPointError(6,0,192);
  grMC->SetPointError(7,0,184);
  grMC->SetPointError(8,0,191);
  grMC->SetPointError(9,0,191);

  grMC->SetPointError(10,0,22);
  grMC->SetPointError(11,0,23);
  grMC->SetPointError(12,0,22);
  grMC->SetPointError(13,0,23);
  grMC->SetPointError(14,0,23);
  grMC->SetPointError(15,0,22);
  grMC->SetPointError(16,0,23);
  grMC->SetPointError(17,0,23);
  grMC->SetPointError(18,0,22);
  grMC->SetPointError(19,0,21);

  grMC->SetPointError(20,0,19);
  grMC->SetPointError(21,0,19);
  grMC->SetPointError(22,0,19);
  grMC->SetPointError(23,0,20);
  grMC->SetPointError(24,0,19);
  grMC->SetPointError(25,0,19);
  grMC->SetPointError(26,0,19);
  grMC->SetPointError(27,0,19);
  grMC->SetPointError(28,0,19);
  grMC->SetPointError(29,0,18);

  grMC->SetPointError(30,0,19);
  grMC->SetPointError(31,0,20);
  grMC->SetPointError(32,0,19);
  grMC->SetPointError(33,0,19);
  grMC->SetPointError(34,0,19);
  grMC->SetPointError(35,0,19);
  grMC->SetPointError(36,0,18);
  grMC->SetPointError(37,0,19);
  grMC->SetPointError(38,0,19);
  grMC->SetPointError(39,0,17);

  
  grMC->SetPointError(40,0,16);
  grMC->SetPointError(41,0,15);
  grMC->SetPointError(42,0,15);
  grMC->SetPointError(43,0,15);
  grMC->SetPointError(44,0,14);
  grMC->SetPointError(45,0,14);
  grMC->SetPointError(46,0,14);
  grMC->SetPointError(47,0,13);
  grMC->SetPointError(48,0,15);
  grMC->SetPointError(49,0,14);
  */


  TGraphErrors *grwb = new TGraphErrors();
  
  grwb->SetPoint(0,3,hvz1wb->GetBinContent(1));
  grwb->SetPoint(1,9,hvz1wb->GetBinContent(2));
  grwb->SetPoint(2,15,hvz1wb->GetBinContent(3));
  grwb->SetPoint(3,21,hvz1wb->GetBinContent(4));
  grwb->SetPoint(4,27,hvz1wb->GetBinContent(5));

  grwb->SetPoint(5,38,hvz2wb->GetBinContent(1));
  grwb->SetPoint(6,44,hvz2wb->GetBinContent(2));
  grwb->SetPoint(7,50,hvz2wb->GetBinContent(3));
  grwb->SetPoint(8,56,hvz2wb->GetBinContent(4));
  grwb->SetPoint(9,62,hvz2wb->GetBinContent(5));

  grwb->SetPoint(10,83,hvz3wb->GetBinContent(1));
  grwb->SetPoint(11,89,hvz3wb->GetBinContent(2));
  grwb->SetPoint(12,95,hvz3wb->GetBinContent(3));
  grwb->SetPoint(13,101,hvz3wb->GetBinContent(4));
  grwb->SetPoint(14,107,hvz3wb->GetBinContent(5));
  grwb->SetPoint(15,113,hvz3wb->GetBinContent(6));
  grwb->SetPoint(16,119,hvz3wb->GetBinContent(7));
  grwb->SetPoint(17,125,hvz3wb->GetBinContent(8));
  grwb->SetPoint(18,131,hvz3wb->GetBinContent(9));
  grwb->SetPoint(19,137,hvz3wb->GetBinContent(10));

  grwb->SetPoint(20,158,hvz4wb->GetBinContent(1));
  grwb->SetPoint(21,164,hvz4wb->GetBinContent(2));
  grwb->SetPoint(22,170,hvz4wb->GetBinContent(3));
  grwb->SetPoint(23,176,hvz4wb->GetBinContent(4));
  grwb->SetPoint(24,182,hvz4wb->GetBinContent(5));
  grwb->SetPoint(25,188,hvz4wb->GetBinContent(6));
  grwb->SetPoint(26,194,hvz4wb->GetBinContent(7));
  grwb->SetPoint(27,200,hvz4wb->GetBinContent(8));
  grwb->SetPoint(28,206,hvz4wb->GetBinContent(9));
  grwb->SetPoint(29,212,hvz4wb->GetBinContent(10));

  grwb->SetPoint(30,233,hvz5wb->GetBinContent(1));
  grwb->SetPoint(31,239,hvz5wb->GetBinContent(2));
  grwb->SetPoint(32,245,hvz5wb->GetBinContent(3));
  grwb->SetPoint(33,251,hvz5wb->GetBinContent(4));
  grwb->SetPoint(34,257,hvz5wb->GetBinContent(5));
  grwb->SetPoint(35,263,hvz5wb->GetBinContent(6));
  grwb->SetPoint(36,269,hvz5wb->GetBinContent(7));
  grwb->SetPoint(37,275,hvz5wb->GetBinContent(8));
  grwb->SetPoint(38,281,hvz5wb->GetBinContent(9));
  grwb->SetPoint(39,287,hvz5wb->GetBinContent(10));
  
  grwb->SetPoint(40,308,hvz6wb->GetBinContent(1));
  grwb->SetPoint(41,314,hvz6wb->GetBinContent(2));
  grwb->SetPoint(42,320,hvz6wb->GetBinContent(3));
  grwb->SetPoint(43,326,hvz6wb->GetBinContent(4));
  grwb->SetPoint(44,332,hvz6wb->GetBinContent(5));
  grwb->SetPoint(45,338,hvz6wb->GetBinContent(6));
  grwb->SetPoint(46,344,hvz6wb->GetBinContent(7));
  grwb->SetPoint(47,350,hvz6wb->GetBinContent(8));
  grwb->SetPoint(48,356,hvz6wb->GetBinContent(9));
  grwb->SetPoint(49,362,hvz6wb->GetBinContent(10));

  // ERRORS

  for(int i=0;i<50;i++){
    grwb->SetPointError(i,0,err_DT[i]);
  }
  /*
  grwb->SetPointError(0,0,96);
  grwb->SetPointError(1,0,96);
  grwb->SetPointError(2,0,99);
  grwb->SetPointError(3,0,95);
  grwb->SetPointError(4,0,96);

  grwb->SetPointError(5,0,165);
  grwb->SetPointError(6,0,161);
  grwb->SetPointError(7,0,156);
  grwb->SetPointError(8,0,168);
  grwb->SetPointError(9,0,156);

  grwb->SetPointError(10,0,178);
  grwb->SetPointError(11,0,167);
  grwb->SetPointError(12,0,169);
  grwb->SetPointError(13,0,166);
  grwb->SetPointError(14,0,163);
  grwb->SetPointError(15,0,163);
  grwb->SetPointError(16,0,167);
  grwb->SetPointError(17,0,164);
  grwb->SetPointError(18,0,171);
  grwb->SetPointError(19,0,165);

  grwb->SetPointError(20,0,239);
  grwb->SetPointError(21,0,244);
  grwb->SetPointError(22,0,206);
  grwb->SetPointError(23,0,236);
  grwb->SetPointError(24,0,229);
  grwb->SetPointError(25,0,225);
  grwb->SetPointError(26,0,216);
  grwb->SetPointError(27,0,251);
  grwb->SetPointError(28,0,215);
  grwb->SetPointError(29,0,212);

  grwb->SetPointError(30,0,161);
  grwb->SetPointError(31,0,155);
  grwb->SetPointError(32,0,173);
  grwb->SetPointError(33,0,164);
  grwb->SetPointError(34,0,162);
  grwb->SetPointError(35,0,173);
  grwb->SetPointError(36,0,184);
  grwb->SetPointError(37,0,169);
  grwb->SetPointError(38,0,163);
  grwb->SetPointError(39,0,168);

  
  grwb->SetPointError(40,0,173);
  grwb->SetPointError(41,0,174);
  grwb->SetPointError(42,0,174);
  grwb->SetPointError(43,0,172);
  grwb->SetPointError(44,0,166);
  grwb->SetPointError(45,0,163);
  grwb->SetPointError(46,0,173);
  grwb->SetPointError(47,0,163);
  grwb->SetPointError(48,0,166);
  grwb->SetPointError(49,0,162);
  */
  //// log_file
  
  log_file << 0 << " " << 0 << " " << 3 << " " << hvz1MCfull->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz1MCfull->GetBinContent(1)) << endl;
  log_file << 0 << " " << 1 << " " << 9 << " " << hvz1MCfull->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz1MCfull->GetBinContent(2)) << endl;
  log_file << 0 << " " << 2 << " " << 15 << " " << hvz1MCfull->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz1MCfull->GetBinContent(3)) << endl;
  log_file << 0 << " " << 3 << " " << 21 << " " << hvz1MCfull->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz1MCfull->GetBinContent(4)) << endl;
  log_file << 0 << " " << 4 << " " << 27 << " " << hvz1MCfull->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz1MCfull->GetBinContent(5)) << endl;

  log_file << 0 << " " << 5 << " " << 38 << " " << hvz2MCfull->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz2MCfull->GetBinContent(1)) << endl;
  log_file << 0 << " " << 6 << " " << 44 << " " << hvz2MCfull->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz2MCfull->GetBinContent(2)) << endl;
  log_file << 0 << " " << 7 << " " << 50 << " " << hvz2MCfull->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz2MCfull->GetBinContent(3)) << endl;
  log_file << 0 << " " << 8 << " " << 56 << " " << hvz2MCfull->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz2MCfull->GetBinContent(4)) << endl;
  log_file << 0 << " " << 9 << " " << 62 << " " << hvz2MCfull->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz2MCfull->GetBinContent(5)) << endl;

  log_file << 0 << " " << 10 << " " << 83 << " " << hvz3MCfull->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(1)) << endl;
  log_file << 0 << " " << 11 << " " << 89 << " " << hvz3MCfull->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(2)) << endl;
  log_file << 0 << " " << 12 << " " << 95 << " " << hvz3MCfull->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(3)) << endl;
  log_file << 0 << " " << 13 << " " << 101 << " " << hvz3MCfull->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(4)) << endl;
  log_file << 0 << " " << 14 << " " << 107 << " " << hvz3MCfull->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(5)) << endl;
  log_file << 0 << " " << 15 << " " << 113 << " " << hvz3MCfull->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(6)) << endl;
  log_file << 0 << " " << 16 << " " << 119 << " " << hvz3MCfull->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(7)) << endl;
  log_file << 0 << " " << 17 << " " << 125 << " " << hvz3MCfull->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(8)) << endl;
  log_file << 0 << " " << 18 << " " << 131 << " " << hvz3MCfull->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(9)) << endl;
  log_file << 0 << " " << 19 << " " << 137 << " " << hvz3MCfull->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz3MCfull->GetBinContent(10)) << endl;

  log_file << 0 << " " << 20 << " " << 158 << " " << hvz4MCfull->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(1)) << endl;
  log_file << 0 << " " << 21 << " " << 164 << " " << hvz4MCfull->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(2)) << endl;
  log_file << 0 << " " << 22 << " " << 170 << " " << hvz4MCfull->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(3)) << endl;
  log_file << 0 << " " << 23 << " " << 176 << " " << hvz4MCfull->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(4)) << endl;
  log_file << 0 << " " << 24 << " " << 182 << " " << hvz4MCfull->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(5)) << endl;
  log_file << 0 << " " << 25 << " " << 188 << " " << hvz4MCfull->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(6)) << endl;
  log_file << 0 << " " << 26 << " " << 194 << " " << hvz4MCfull->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(7)) << endl;
  log_file << 0 << " " << 27 << " " << 200 << " " << hvz4MCfull->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(8)) << endl;
  log_file << 0 << " " << 28 << " " << 206 << " " << hvz4MCfull->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(9)) << endl;
  log_file << 0 << " " << 29 << " " << 212 << " " << hvz4MCfull->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz4MCfull->GetBinContent(10)) << endl;

  log_file << 0 << " " << 30 << " " << 233 << " " << hvz5MCfull->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(1)) << endl;
  log_file << 0 << " " << 31 << " " << 239 << " " << hvz5MCfull->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(2)) << endl;
  log_file << 0 << " " << 32 << " " << 245 << " " << hvz5MCfull->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(3)) << endl;
  log_file << 0 << " " << 33 << " " << 251 << " " << hvz5MCfull->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(4)) << endl;
  log_file << 0 << " " << 34 << " " << 257 << " " << hvz5MCfull->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(5)) << endl;
  log_file << 0 << " " << 35 << " " << 263 << " " << hvz5MCfull->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(6)) << endl;
  log_file << 0 << " " << 36 << " " << 269 << " " << hvz5MCfull->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(7)) << endl;
  log_file << 0 << " " << 37 << " " << 275 << " " << hvz5MCfull->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(8)) << endl;
  log_file << 0 << " " << 38 << " " << 281 << " " << hvz5MCfull->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(9)) << endl;
  log_file << 0 << " " << 39 << " " << 287 << " " << hvz5MCfull->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz5MCfull->GetBinContent(10)) << endl;
  
  log_file << 0 << " " << 40 << " " << 308 << " " << hvz6MCfull->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(1)) << endl;
  log_file << 0 << " " << 41 << " " << 314 << " " << hvz6MCfull->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(2)) << endl;
  log_file << 0 << " " << 42 << " " << 320 << " " << hvz6MCfull->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(3)) << endl;
  log_file << 0 << " " << 43 << " " << 326 << " " << hvz6MCfull->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(4)) << endl;
  log_file << 0 << " " << 44 << " " << 332 << " " << hvz6MCfull->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(5)) << endl;
  log_file << 0 << " " << 45 << " " << 338 << " " << hvz6MCfull->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(6)) << endl;
  log_file << 0 << " " << 46 << " " << 344 << " " << hvz6MCfull->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(7)) << endl;
  log_file << 0 << " " << 47 << " " << 350 << " " << hvz6MCfull->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(8)) << endl;
  log_file << 0 << " " << 48 << " " << 356 << " " << hvz6MCfull->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(9)) << endl;
  log_file << 0 << " " << 49 << " " << 362 << " " << hvz6MCfull->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz6MCfull->GetBinContent(10)) << endl;

  //////////////

  log_file << 1 << " " << 0 << " " << 3 << " " << hvz1MChd->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz1MChd->GetBinContent(1)) << endl;
  log_file << 1 << " " << 1 << " " << 9 << " " << hvz1MChd->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz1MChd->GetBinContent(2)) << endl;
  log_file << 1 << " " << 2 << " " << 15 << " " << hvz1MChd->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz1MChd->GetBinContent(3)) << endl;
  log_file << 1 << " " << 3 << " " << 21 << " " << hvz1MChd->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz1MChd->GetBinContent(4)) << endl;
  log_file << 1 << " " << 4 << " " << 27 << " " << hvz1MChd->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz1MChd->GetBinContent(5)) << endl;

  log_file << 1 << " " << 5 << " " << 38 << " " << hvz2MChd->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz2MChd->GetBinContent(1)) << endl;
  log_file << 1 << " " << 6 << " " << 44 << " " << hvz2MChd->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz2MChd->GetBinContent(2)) << endl;
  log_file << 1 << " " << 7 << " " << 50 << " " << hvz2MChd->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz2MChd->GetBinContent(3)) << endl;
  log_file << 1 << " " << 8 << " " << 56 << " " << hvz2MChd->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz2MChd->GetBinContent(4)) << endl;
  log_file << 1 << " " << 9 << " " << 62 << " " << hvz2MChd->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz2MChd->GetBinContent(5)) << endl;

  log_file << 1 << " " << 10 << " " << 83 << " " << hvz3MChd->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(1)) << endl;
  log_file << 1 << " " << 11 << " " << 89 << " " << hvz3MChd->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(2)) << endl;
  log_file << 1 << " " << 12 << " " << 95 << " " << hvz3MChd->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(3)) << endl;
  log_file << 1 << " " << 13 << " " << 101 << " " << hvz3MChd->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(4)) << endl;
  log_file << 1 << " " << 14 << " " << 107 << " " << hvz3MChd->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(5)) << endl;
  log_file << 1 << " " << 15 << " " << 113 << " " << hvz3MChd->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(6)) << endl;
  log_file << 1 << " " << 16 << " " << 119 << " " << hvz3MChd->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(7)) << endl;
  log_file << 1 << " " << 17 << " " << 125 << " " << hvz3MChd->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(8)) << endl;
  log_file << 1 << " " << 18 << " " << 131 << " " << hvz3MChd->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(9)) << endl;
  log_file << 1 << " " << 19 << " " << 137 << " " << hvz3MChd->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz3MChd->GetBinContent(10)) << endl;

  log_file << 1 << " " << 20 << " " << 158 << " " << hvz4MChd->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(1)) << endl;
  log_file << 1 << " " << 21 << " " << 164 << " " << hvz4MChd->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(2)) << endl;
  log_file << 1 << " " << 22 << " " << 170 << " " << hvz4MChd->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(3)) << endl;
  log_file << 1 << " " << 23 << " " << 176 << " " << hvz4MChd->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(4)) << endl;
  log_file << 1 << " " << 24 << " " << 182 << " " << hvz4MChd->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(5)) << endl;
  log_file << 1 << " " << 25 << " " << 188 << " " << hvz4MChd->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(6)) << endl;
  log_file << 1 << " " << 26 << " " << 194 << " " << hvz4MChd->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(7)) << endl;
  log_file << 1 << " " << 27 << " " << 200 << " " << hvz4MChd->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(8)) << endl;
  log_file << 1 << " " << 28 << " " << 206 << " " << hvz4MChd->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(9)) << endl;
  log_file << 1 << " " << 29 << " " << 212 << " " << hvz4MChd->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz4MChd->GetBinContent(10)) << endl;

  log_file << 1 << " " << 30 << " " << 233 << " " << hvz5MChd->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(1)) << endl;
  log_file << 1 << " " << 31 << " " << 239 << " " << hvz5MChd->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(2)) << endl;
  log_file << 1 << " " << 32 << " " << 245 << " " << hvz5MChd->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(3)) << endl;
  log_file << 1 << " " << 33 << " " << 251 << " " << hvz5MChd->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(4)) << endl;
  log_file << 1 << " " << 34 << " " << 257 << " " << hvz5MChd->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(5)) << endl;
  log_file << 1 << " " << 35 << " " << 263 << " " << hvz5MChd->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(6)) << endl;
  log_file << 1 << " " << 36 << " " << 269 << " " << hvz5MChd->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(7)) << endl;
  log_file << 1 << " " << 37 << " " << 275 << " " << hvz5MChd->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(8)) << endl;
  log_file << 1 << " " << 38 << " " << 281 << " " << hvz5MChd->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(9)) << endl;
  log_file << 1 << " " << 39 << " " << 287 << " " << hvz5MChd->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz5MChd->GetBinContent(10)) << endl;
  
  log_file << 1 << " " << 40 << " " << 308 << " " << hvz6MChd->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(1)) << endl;
  log_file << 1 << " " << 41 << " " << 314 << " " << hvz6MChd->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(2)) << endl;
  log_file << 1 << " " << 42 << " " << 320 << " " << hvz6MChd->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(3)) << endl;
  log_file << 1 << " " << 43 << " " << 326 << " " << hvz6MChd->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(4)) << endl;
  log_file << 1 << " " << 44 << " " << 332 << " " << hvz6MChd->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(5)) << endl;
  log_file << 1 << " " << 45 << " " << 338 << " " << hvz6MChd->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(6)) << endl;
  log_file << 1 << " " << 46 << " " << 344 << " " << hvz6MChd->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(7)) << endl;
  log_file << 1 << " " << 47 << " " << 350 << " " << hvz6MChd->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(8)) << endl;
  log_file << 1 << " " << 48 << " " << 356 << " " << hvz6MChd->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(9)) << endl;
  log_file << 1 << " " << 49 << " " << 362 << " " << hvz6MChd->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz6MChd->GetBinContent(10)) << endl;


  //////////////

  log_file << 2 << " " << 0 << " " << 3 << " " << hvz1MCpr->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz1MCpr->GetBinContent(1)) << endl;
  log_file << 2 << " " << 1 << " " << 9 << " " << hvz1MCpr->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz1MCpr->GetBinContent(2)) << endl;
  log_file << 2 << " " << 2 << " " << 15 << " " << hvz1MCpr->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz1MCpr->GetBinContent(3)) << endl;
  log_file << 2 << " " << 3 << " " << 21 << " " << hvz1MCpr->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz1MCpr->GetBinContent(4)) << endl;
  log_file << 2 << " " << 4 << " " << 27 << " " << hvz1MCpr->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz1MCpr->GetBinContent(5)) << endl;

  log_file << 2 << " " << 5 << " " << 38 << " " << hvz2MCpr->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz2MCpr->GetBinContent(1)) << endl;
  log_file << 2 << " " << 6 << " " << 44 << " " << hvz2MCpr->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz2MCpr->GetBinContent(2)) << endl;
  log_file << 2 << " " << 7 << " " << 50 << " " << hvz2MCpr->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz2MCpr->GetBinContent(3)) << endl;
  log_file << 2 << " " << 8 << " " << 56 << " " << hvz2MCpr->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz2MCpr->GetBinContent(4)) << endl;
  log_file << 2 << " " << 9 << " " << 62 << " " << hvz2MCpr->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz2MCpr->GetBinContent(5)) << endl;

  log_file << 2 << " " << 10 << " " << 83 << " " << hvz3MCpr->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(1)) << endl;
  log_file << 2 << " " << 11 << " " << 89 << " " << hvz3MCpr->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(2)) << endl;
  log_file << 2 << " " << 12 << " " << 95 << " " << hvz3MCpr->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(3)) << endl;
  log_file << 2 << " " << 13 << " " << 101 << " " << hvz3MCpr->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(4)) << endl;
  log_file << 2 << " " << 14 << " " << 107 << " " << hvz3MCpr->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(5)) << endl;
  log_file << 2 << " " << 15 << " " << 113 << " " << hvz3MCpr->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(6)) << endl;
  log_file << 2 << " " << 16 << " " << 119 << " " << hvz3MCpr->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(7)) << endl;
  log_file << 2 << " " << 17 << " " << 125 << " " << hvz3MCpr->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(8)) << endl;
  log_file << 2 << " " << 18 << " " << 131 << " " << hvz3MCpr->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(9)) << endl;
  log_file << 2 << " " << 19 << " " << 137 << " " << hvz3MCpr->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz3MCpr->GetBinContent(10)) << endl;

  log_file << 2 << " " << 20 << " " << 158 << " " << hvz4MCpr->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(1)) << endl;
  log_file << 2 << " " << 21 << " " << 164 << " " << hvz4MCpr->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(2)) << endl;
  log_file << 2 << " " << 22 << " " << 170 << " " << hvz4MCpr->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(3)) << endl;
  log_file << 2 << " " << 23 << " " << 176 << " " << hvz4MCpr->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(4)) << endl;
  log_file << 2 << " " << 24 << " " << 182 << " " << hvz4MCpr->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(5)) << endl;
  log_file << 2 << " " << 25 << " " << 188 << " " << hvz4MCpr->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(6)) << endl;
  log_file << 2 << " " << 26 << " " << 194 << " " << hvz4MCpr->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(7)) << endl;
  log_file << 2 << " " << 27 << " " << 200 << " " << hvz4MCpr->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(8)) << endl;
  log_file << 2 << " " << 28 << " " << 206 << " " << hvz4MCpr->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(9)) << endl;
  log_file << 2 << " " << 29 << " " << 212 << " " << hvz4MCpr->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz4MCpr->GetBinContent(10)) << endl;

  log_file << 2 << " " << 30 << " " << 233 << " " << hvz5MCpr->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(1)) << endl;
  log_file << 2 << " " << 31 << " " << 239 << " " << hvz5MCpr->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(2)) << endl;
  log_file << 2 << " " << 32 << " " << 245 << " " << hvz5MCpr->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(3)) << endl;
  log_file << 2 << " " << 33 << " " << 251 << " " << hvz5MCpr->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(4)) << endl;
  log_file << 2 << " " << 34 << " " << 257 << " " << hvz5MCpr->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(5)) << endl;
  log_file << 2 << " " << 35 << " " << 263 << " " << hvz5MCpr->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(6)) << endl;
  log_file << 2 << " " << 36 << " " << 269 << " " << hvz5MCpr->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(7)) << endl;
  log_file << 2 << " " << 37 << " " << 275 << " " << hvz5MCpr->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(8)) << endl;
  log_file << 2 << " " << 38 << " " << 281 << " " << hvz5MCpr->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(9)) << endl;
  log_file << 2 << " " << 39 << " " << 287 << " " << hvz5MCpr->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz5MCpr->GetBinContent(10)) << endl;
  
  log_file << 2 << " " << 40 << " " << 308 << " " << hvz6MCpr->GetBinContent(1) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(1)) << endl;
  log_file << 2 << " " << 41 << " " << 314 << " " << hvz6MCpr->GetBinContent(2) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(2)) << endl;
  log_file << 2 << " " << 42 << " " << 320 << " " << hvz6MCpr->GetBinContent(3) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(3)) << endl;
  log_file << 2 << " " << 43 << " " << 326 << " " << hvz6MCpr->GetBinContent(4) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(4)) << endl;
  log_file << 2 << " " << 44 << " " << 332 << " " << hvz6MCpr->GetBinContent(5) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(5)) << endl;
  log_file << 2 << " " << 45 << " " << 338 << " " << hvz6MCpr->GetBinContent(6) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(6)) << endl;
  log_file << 2 << " " << 46 << " " << 344 << " " << hvz6MCpr->GetBinContent(7) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(7)) << endl;
  log_file << 2 << " " << 47 << " " << 350 << " " << hvz6MCpr->GetBinContent(8) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(8)) << endl;
  log_file << 2 << " " << 48 << " " << 356 << " " << hvz6MCpr->GetBinContent(9) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(9)) << endl;
  log_file << 2 << " " << 49 << " " << 362 << " " << hvz6MCpr->GetBinContent(10) << " " << 0 << " " << TMath::Sqrt(hvz6MCpr->GetBinContent(10)) << endl;

   //////////////

  log_file << 3 << " " << 0 << " " << 3 << " " << hvz1wb->GetBinContent(1) << " " << 0 << " " << 29 << endl;
  log_file << 3 << " " << 1 << " " << 9 << " " << hvz1wb->GetBinContent(2) << " " << 0 << " " << 29 << endl;
  log_file << 3 << " " << 2 << " " << 15 << " " << hvz1wb->GetBinContent(3) << " " << 0 << " " << 32 << endl;
  log_file << 3 << " " << 3 << " " << 21 << " " << hvz1wb->GetBinContent(4) << " " << 0 << " " << 29 << endl;
  log_file << 3 << " " << 4 << " " << 27 << " " << hvz1wb->GetBinContent(5) << " " << 0 << " " << 30 << endl;

  log_file << 3 << " " << 5 << " " << 38 << " " << hvz2wb->GetBinContent(1) << " " << 0 << " " << 40 << endl;
  log_file << 3 << " " << 6 << " " << 44 << " " << hvz2wb->GetBinContent(2) << " " << 0 << " " << 38 << endl;
  log_file << 3 << " " << 7 << " " << 50 << " " << hvz2wb->GetBinContent(3) << " " << 0 << " " << 36 << endl;
  log_file << 3 << " " << 8 << " " << 56 << " " << hvz2wb->GetBinContent(4) << " " << 0 << " " << 40 << endl;
  log_file << 3 << " " << 9 << " " << 62 << " " << hvz2wb->GetBinContent(5) << " " << 0 << " " << 36 << endl;

  log_file << 3 << " " << 10 << " " << 83 << " " << hvz3wb->GetBinContent(1) << " " << 0 << " " << 117 << endl;
  log_file << 3 << " " << 11 << " " << 89 << " " << hvz3wb->GetBinContent(2) << " " << 0 << " " << 104 << endl;
  log_file << 3 << " " << 12 << " " << 95 << " " << hvz3wb->GetBinContent(3) << " " << 0 << " " << 106 << endl;
  log_file << 3 << " " << 13 << " " << 101 << " " << hvz3wb->GetBinContent(4) << " " << 0 << " " << 103 << endl;
  log_file << 3 << " " << 14 << " " << 107 << " " << hvz3wb->GetBinContent(5) << " " << 0 << " " << 99 <<  endl;
  log_file << 3 << " " << 15 << " " << 113 << " " << hvz3wb->GetBinContent(6) << " " << 0 << " " << 98 << endl;
  log_file << 3 << " " << 16 << " " << 119 << " " << hvz3wb->GetBinContent(7) << " " << 0 << " " << 104 <<  endl;
  log_file << 3 << " " << 17 << " " << 125 << " " << hvz3wb->GetBinContent(8) << " " << 0 << " " << 100 <<  endl;
  log_file << 3 << " " << 18 << " " << 131 << " " << hvz3wb->GetBinContent(9) << " " << 0 << " " << 109 <<  endl;
  log_file << 3 << " " << 19 << " " << 137 << " " << hvz3wb->GetBinContent(10) << " " << 0 << " " << 101 << endl;

  log_file << 3 << " " << 20 << " " << 158 << " " << hvz4wb->GetBinContent(1) << " " << 0 << " " << 170 << endl;
  log_file << 3 << " " << 21 << " " << 164 << " " << hvz4wb->GetBinContent(2) << " " << 0 << " " << 176 << endl;
  log_file << 3 << " " << 22 << " " << 170 << " " << hvz4wb->GetBinContent(3) << " " << 0 << " " << 128 << endl;
  log_file << 3 << " " << 23 << " " << 176 << " " << hvz4wb->GetBinContent(4) << " " << 0 << " " << 166 << endl;
  log_file << 3 << " " << 24 << " " << 182 << " " << hvz4wb->GetBinContent(5) << " " << 0 << " " << 157 << endl;
  log_file << 3 << " " << 25 << " " << 188 << " " << hvz4wb->GetBinContent(6) << " " << 0 << " " << 152 << endl;
  log_file << 3 << " " << 26 << " " << 194 << " " << hvz4wb->GetBinContent(7) << " " << 0 << " " << 141 << endl;
  log_file << 3 << " " << 27 << " " << 200 << " " << hvz4wb->GetBinContent(8) << " " << 0 << " " << 184 << endl;
  log_file << 3 << " " << 28 << " " << 206 << " " << hvz4wb->GetBinContent(9) << " " << 0 << " " << 140 << endl;
  log_file << 3 << " " << 29 << " " << 212 << " " << hvz4wb->GetBinContent(10) << " " << 0 << " " << 135 << endl;

  log_file << 3 << " " << 30 << " " << 233 << " " << hvz5wb->GetBinContent(1) << " " << 0 << " " << 112 << endl;
  log_file << 3 << " " << 31 << " " << 239 << " " << hvz5wb->GetBinContent(2) << " " << 0 << " " << 103 << endl;
  log_file << 3 << " " << 32 << " " << 245 << " " << hvz5wb->GetBinContent(3) << " " << 0 << " " << 128 << endl;
  log_file << 3 << " " << 33 << " " << 251 << " " << hvz5wb->GetBinContent(4) << " " << 0 << " " << 116 << endl;
  log_file << 3 << " " << 34 << " " << 257 << " " << hvz5wb->GetBinContent(5) << " " << 0 << " " << 114 << endl;
  log_file << 3 << " " << 35 << " " << 263 << " " << hvz5wb->GetBinContent(6) << " " << 0 << " " << 128 << endl;
  log_file << 3 << " " << 36 << " " << 269 << " " << hvz5wb->GetBinContent(7) << " " << 0 << " " << 142 << endl;
  log_file << 3 << " " << 37 << " " << 275 << " " << hvz5wb->GetBinContent(8) << " " << 0 << " " << 123 << endl;
  log_file << 3 << " " << 38 << " " << 281 << " " << hvz5wb->GetBinContent(9) << " " << 0 << " " << 115 << endl;
  log_file << 3 << " " << 39 << " " << 287 << " " << hvz5wb->GetBinContent(10) << " " << 0 << " " << 122 << endl;
  
  log_file << 3 << " " << 40 << " " << 308 << " " << hvz6wb->GetBinContent(1) << " " << 0 << " " << 122 << endl;
  log_file << 3 << " " << 41 << " " << 314 << " " << hvz6wb->GetBinContent(2) << " " << 0 << " " << 123 << endl;
  log_file << 3 << " " << 42 << " " << 320 << " " << hvz6wb->GetBinContent(3) << " " << 0 << " " << 124 << endl;
  log_file << 3 << " " << 43 << " " << 326 << " " << hvz6wb->GetBinContent(4) << " " << 0 << " " << 120 << endl;
  log_file << 3 << " " << 44 << " " << 332 << " " << hvz6wb->GetBinContent(5) << " " << 0 << " " << 110 << endl;
  log_file << 3 << " " << 45 << " " << 338 << " " << hvz6wb->GetBinContent(6) << " " << 0 << " " << 105 << endl;
  log_file << 3 << " " << 46 << " " << 344 << " " << hvz6wb->GetBinContent(7) << " " << 0 << " " << 122 << endl;
  log_file << 3 << " " << 47 << " " << 350 << " " << hvz6wb->GetBinContent(8) << " " << 0 << " " << 105 << endl;
  log_file << 3 << " " << 48 << " " << 356 << " " << hvz6wb->GetBinContent(9) << " " << 0 << " " << 110 << endl;
  log_file << 3 << " " << 49 << " " << 362 << " " << hvz6wb->GetBinContent(10) << " " << 0 << " " << 103 << endl;


  

  gStyle->SetOptFit(1111);
  
  

  grMCpr->SetLineWidth(2);
  grMCpr->SetLineColor(4);
  //grMCpr->SetMarkerStyle(21);
  grMCpr->SetMarkerColor(4);
  //grMCpr->Draw("APL && same");
  grMCpr->Fit("expo");
  grMCpr->SetTitle("MC protons");

  TF1 *hd_fit = new TF1("hd_fit","[0]*TMath::Log(x)",0,365);
  hd_fit->SetParameter(0, 180);
  //hd_fit->SetParameter(1, 0.79);
  //hd_fit->SetParameter(2, -170);
  
  grMChd->SetLineWidth(2);
  grMChd->SetLineColor(3);
  //grMChd->SetMarkerStyle(22);
  grMChd->SetMarkerColor(3);
  //grMChd->Fit("expo");
  grMChd->Fit("hd_fit");
  grMChd->SetTitle("MC hadrons");

  TF1 *fit = new TF1("fit","expo + [2]*TMath::Log(x)",0,365);
  fit->SetParameter(0, 7.8);
  fit->SetParameter(1, -0.0058);
  fit->SetParLimits(2, 135,145);
  //fit->SetParLimits(3, -1.62,-1.22);
  //fit->FixParameter(4, -170);

  grMC->SetLineWidth(2);
  grMC->SetLineColor(1);
  grMC->SetMarkerStyle(20);
  grMC->SetMarkerColor(1);
  //grMChd->Fit("expo");
  grMC->Fit(fit);
  grMC->SetTitle("Full MC");

  fit->SetParameter(0, 7.8);
  fit->SetParameter(1, -0.0058);
  fit->SetParLimits(2, 135,145);
  
  grwb->SetLineWidth(2);
  //grwb->SetLineColor(4);
  grwb->SetMarkerStyle(20);
  grwb->SetMarkerColor(3);
  //grwb->Draw("AP");
  //grwb->Fit("expo");
  grwb->Fit(fit);
  grwb->SetTitle("Full data");

  
  
 
  TMultiGraph *mgr = new TMultiGraph();
  mgr->Add(grwb);
  mgr->Add(grMC);
  mgr->GetYaxis()->SetRangeUser(0,3000);
  mgr->GetXaxis()->SetLimits(0,365);
  mgr->GetXaxis()->SetTitle("vz [mm]");
  mgr->GetYaxis()->SetTitle("FULL DATA");
  mgr->Draw("AP");

  /*
  TCanvas *c1 = new TCanvas();
  grMCpr->GetYaxis()->SetRangeUser(0,3500);
  grMCpr->GetXaxis()->SetLimits(0,365);
  grMCpr->GetXaxis()->SetTitle("vz [mm]");
  grMCpr->GetYaxis()->SetTitle("Proton vertices");
  grMCpr->Draw("APB");
  gr->Draw("P && same");
  //grwb->Draw("P && sames");
  */
  TCanvas *c2 = new TCanvas();
  grMC->GetYaxis()->SetRangeUser(0,3500);
  grMC->GetXaxis()->SetLimits(0,365);
  grMC->GetXaxis()->SetTitle("vz [mm]");
  grMC->GetYaxis()->SetTitle("Selected vertices");
  grMC->Draw("AP");
  grMCpr->Draw("P && same");
  grMChd->Draw("P && same");
  // grwb->Draw("P && same");

log_file.close();
  
}
