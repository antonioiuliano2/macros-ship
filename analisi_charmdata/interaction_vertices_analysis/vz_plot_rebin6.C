// // // ------------------------->
// // // // ----------------------->
// // // ------------------------->
// // --------------------------->
// ----------------------------->

void vz_plot_rebin6(){

   //TString dir = gSystem->UnixPathName(gInterpreter->GetCurrentMacroName());
   //dir.ReplaceAll("vz_plot_rebin6.C","");
   //dir.ReplaceAll("/./","/");
   TString dir = TString("/home/utente/Lavoro/BDT_vertices_Valerio/afterBDT_plots/");
   ifstream in;
   in.open(Form("%serrors.dat",dir.Data()));

   Int_t iline;
   Float_t x,y,z,u,v;
   Int_t nlines = 0;
   Float_t err_MC[30]={};
   Float_t err_MC_pr[30]={};
   Float_t err_MC_hd[30]={};
   Float_t err_DT[30]={};

   while (1) {
     in >> x >> y >> z >> u >> v;
     err_DT[nlines]=y;
     err_MC[nlines]=z;
     err_MC_pr[nlines]=u;
     err_MC_hd[nlines]=v;
     if (!in.good()) break;
     printf("x=%8f, y=%8f, z=%8f, u=%8f, v=%8f\n",x,y,z,u,v);
     nlines++;
   }
   printf(" found %d points\n",nlines);

   in.close();

  TFile *f = new TFile((dir+TString("full_vz.root")).Data(),"READ"); //full_vz_1quarter.root
  //TFile *fb = new TFile("with_bkg_full_vz_1quarter.root","READ");
  TFile *fb = new TFile((dir+TString("with_bkg_full_vz.root")).Data(),"READ"); //dati

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
  /*
  hvz1MCpr->Scale(4); 
  hvz1MChd->Scale(4); 
  hvz2MCpr->Scale(4); 
  hvz2MChd->Scale(4);
  hvz3MCpr->Scale(4); 
  hvz3MChd->Scale(4); 
  hvz4MCpr->Scale(4); 
  hvz4MChd->Scale(4);
  hvz5MCpr->Scale(4); 
  hvz5MChd->Scale(4); 
  hvz6MCpr->Scale(4); 
  hvz6MChd->Scale(4);

  hvz1wb->Scale(4);
  hvz2wb->Scale(4);
  hvz3wb->Scale(4);
  hvz4wb->Scale(4);
  hvz5wb->Scale(4);
  hvz6wb->Scale(4);
  */
  
  // ARTEFACTS
  //hvz2MCpr->Scale(2); // per rebin6
  //hvz2MChd->Scale(2); // per rebin6
  //hvz1MCpr->Scale(2); // per rebin6
  //hvz1MChd->Scale(2); // per rebin6
  //hvz1wb->Scale(2); // per rebin6
  //hvz2wb->Scale(2); // per rebin6

  /*float eff_pr_ch1 = 0.896;
  float eff_hd_ch1 = 0.691;
  float eff_pr_ch2 = 0.795;
  float eff_hd_ch2 = 0.528;
  float eff_pr_ch3 = 0.751;
  float eff_hd_ch3 = 0.394;
  float eff_pr_ch4 = 0.728;
  float eff_hd_ch4 = 0.365;
  float eff_pr_ch5 = 0.744;
  float eff_hd_ch5 = 0.348;
  float eff_pr_ch6 = 0.722;
  float eff_hd_ch6 = 0.295;*/
    
  hvz1->Scale(0.25);
  hvz2->Scale(0.33);
  hvz3->Scale(3.86);
  hvz4->Scale(4.91);
  hvz5->Scale(3.75);
  hvz6->Scale(6.0);

  hvz1->Rebin(3);
  hvz2->Rebin(3);
  hvz3->Rebin(6);
  hvz4->Rebin(6);
  hvz5->Rebin(6);
  hvz6->Rebin(6);

  /*
  hvz1MCpr->Scale(1./eff_pr_ch1);
  hvz2MCpr->Scale(1./eff_pr_ch2);
  hvz3MCpr->Scale(1./eff_pr_ch3);
  hvz4MCpr->Scale(1./eff_pr_ch4);
  hvz5MCpr->Scale(1./eff_pr_ch5);
  hvz6MCpr->Scale(1./eff_pr_ch6);*/
  
  hvz1MCpr->Scale(1);    // C_mc = MCPOT_CH1 / MCPOT_CHX    135000
  hvz2MCpr->Scale(1);                              //135000
  hvz3MCpr->Scale(3.86);                           //40000
  hvz4MCpr->Scale(4.91);
  hvz5MCpr->Scale(3.375);           // err MC = C_mc * rad(N)
  hvz6MCpr->Scale(6.0);             
  
  hvz1MCpr->Rebin(3);
  hvz2MCpr->Rebin(3);
  hvz3MCpr->Rebin(6);
  hvz4MCpr->Rebin(6);
  hvz5MCpr->Rebin(6);
  hvz6MCpr->Rebin(6);

  /*
  hvz1MChd->Scale(1./eff_hd_ch1);
  hvz2MChd->Scale(1./eff_hd_ch2);
  hvz3MChd->Scale(1./eff_hd_ch3);
  hvz4MChd->Scale(1./eff_hd_ch4);
  hvz5MChd->Scale(1./eff_hd_ch5);
  hvz6MChd->Scale(1./eff_hd_ch6);*/
  
  hvz1MChd->Scale(1);
  hvz2MChd->Scale(1);
  hvz3MChd->Scale(3.86);
  hvz4MChd->Scale(4.91);
  hvz5MChd->Scale(3.375);
  hvz6MChd->Scale(6.0);
  
  
  hvz1MChd->Rebin(3);
  hvz2MChd->Rebin(3);
  hvz3MChd->Rebin(6);
  hvz4MChd->Rebin(6);
  hvz5MChd->Rebin(6);
  hvz6MChd->Rebin(6);

  hvz1wb->Scale(0.25);
  hvz2wb->Scale(0.33);
  hvz3wb->Scale(3.86);
  hvz4wb->Scale(4.91);
  hvz5wb->Scale(3.375);
  hvz6wb->Scale(6.0);

  // AGGIUNGERE SCALE FACTOR
  

  hvz1wb->Rebin(3);
  hvz2wb->Rebin(3);
  hvz3wb->Rebin(6);
  hvz4wb->Rebin(6);
  hvz5wb->Rebin(6);
  hvz6wb->Rebin(6);

  TH1F *hvz1MCfull = new TH1F("hMC1","hMC1",5,-35,-5);
  TH1F *hvz2MCfull = new TH1F("hMC2","hMC2",5,-35,-5);
  TH1F *hvz3MCfull = new TH1F("hMC3","hMC3",5,-65,-5);
  TH1F *hvz4MCfull = new TH1F("hMC4","hMC4",5,-65,-5);
  TH1F *hvz5MCfull = new TH1F("hMC5","hMC5",5,-65,-5);
  TH1F *hvz6MCfull = new TH1F("hMC6","hMC6",5,-65,-5);

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

  int rebin_factor = 2;

 

  TGraphErrors *grMCpr = new TGraphErrors();

  grMCpr->SetPoint(0,3,rebin_factor*hvz1MCpr->GetBinContent(1));
  grMCpr->SetPoint(1,9,rebin_factor*hvz1MCpr->GetBinContent(2));
  grMCpr->SetPoint(2,15,rebin_factor*hvz1MCpr->GetBinContent(3));
  grMCpr->SetPoint(3,21,rebin_factor*hvz1MCpr->GetBinContent(4));
  grMCpr->SetPoint(4,27,rebin_factor*hvz1MCpr->GetBinContent(5));

  grMCpr->SetPoint(5,38,rebin_factor*hvz2MCpr->GetBinContent(1));
  grMCpr->SetPoint(6,44,rebin_factor*hvz2MCpr->GetBinContent(2));
  grMCpr->SetPoint(7,50,rebin_factor*hvz2MCpr->GetBinContent(3));
  grMCpr->SetPoint(8,56,rebin_factor*hvz2MCpr->GetBinContent(4));
  grMCpr->SetPoint(9,62,rebin_factor*hvz2MCpr->GetBinContent(5));


  grMCpr->SetPoint(10,89,hvz3MCpr->GetBinContent(1));
  grMCpr->SetPoint(11,101,hvz3MCpr->GetBinContent(2));
  grMCpr->SetPoint(12,113,hvz3MCpr->GetBinContent(3));
  grMCpr->SetPoint(13,125,hvz3MCpr->GetBinContent(4));
  grMCpr->SetPoint(14,137,hvz3MCpr->GetBinContent(5));


  grMCpr->SetPoint(15,164,hvz4MCpr->GetBinContent(1));
  grMCpr->SetPoint(16,176,hvz4MCpr->GetBinContent(2));
  grMCpr->SetPoint(17,188,hvz4MCpr->GetBinContent(3));
  grMCpr->SetPoint(18,200,hvz4MCpr->GetBinContent(4));
  grMCpr->SetPoint(19,212,hvz4MCpr->GetBinContent(5));
  
  grMCpr->SetPoint(20,239,hvz5MCpr->GetBinContent(1));
  grMCpr->SetPoint(21,251,hvz5MCpr->GetBinContent(2));
  grMCpr->SetPoint(22,263,hvz5MCpr->GetBinContent(3));
  grMCpr->SetPoint(23,275,hvz5MCpr->GetBinContent(4));
  grMCpr->SetPoint(24,287,hvz5MCpr->GetBinContent(5));
  
  grMCpr->SetPoint(25,314,hvz6MCpr->GetBinContent(1));
  grMCpr->SetPoint(26,326,hvz6MCpr->GetBinContent(2));
  grMCpr->SetPoint(27,338,hvz6MCpr->GetBinContent(3));
  grMCpr->SetPoint(28,350,hvz6MCpr->GetBinContent(4));
  grMCpr->SetPoint(29,362,hvz6MCpr->GetBinContent(5));


  
  
  TGraphErrors *grMChd = new TGraphErrors();

  grMChd->SetPoint(0,3,rebin_factor*hvz1MChd->GetBinContent(1));
  grMChd->SetPoint(1,9,rebin_factor*hvz1MChd->GetBinContent(2));
  grMChd->SetPoint(2,15,rebin_factor*hvz1MChd->GetBinContent(3));
  grMChd->SetPoint(3,21,rebin_factor*hvz1MChd->GetBinContent(4));
  grMChd->SetPoint(4,27,rebin_factor*hvz1MChd->GetBinContent(5));

  grMChd->SetPoint(5,38,rebin_factor*hvz2MChd->GetBinContent(1));
  grMChd->SetPoint(6,44,rebin_factor*hvz2MChd->GetBinContent(2));
  grMChd->SetPoint(7,50,rebin_factor*hvz2MChd->GetBinContent(3));
  grMChd->SetPoint(8,56,rebin_factor*hvz2MChd->GetBinContent(4));
  grMChd->SetPoint(9,62,rebin_factor*hvz2MChd->GetBinContent(5));

  grMChd->SetPoint(10,89,hvz3MChd->GetBinContent(1));
  grMChd->SetPoint(11,101,hvz3MChd->GetBinContent(2));
  grMChd->SetPoint(12,113,hvz3MChd->GetBinContent(3));
  grMChd->SetPoint(13,125,hvz3MChd->GetBinContent(4));
  grMChd->SetPoint(14,137,hvz3MChd->GetBinContent(5));


  grMChd->SetPoint(15,164,hvz4MChd->GetBinContent(1));
  grMChd->SetPoint(16,176,hvz4MChd->GetBinContent(2));
  grMChd->SetPoint(17,188,hvz4MChd->GetBinContent(3));
  grMChd->SetPoint(18,200,hvz4MChd->GetBinContent(4));
  grMChd->SetPoint(19,212,hvz4MChd->GetBinContent(5));
  
  grMChd->SetPoint(20,239,hvz5MChd->GetBinContent(1));
  grMChd->SetPoint(21,251,hvz5MChd->GetBinContent(2));
  grMChd->SetPoint(22,263,hvz5MChd->GetBinContent(3));
  grMChd->SetPoint(23,275,hvz5MChd->GetBinContent(4));
  grMChd->SetPoint(24,287,hvz5MChd->GetBinContent(5));
  
  grMChd->SetPoint(25,314,hvz6MChd->GetBinContent(1));
  grMChd->SetPoint(26,326,hvz6MChd->GetBinContent(2));
  grMChd->SetPoint(27,338,hvz6MChd->GetBinContent(3));
  grMChd->SetPoint(28,350,hvz6MChd->GetBinContent(4));
  grMChd->SetPoint(29,362,hvz6MChd->GetBinContent(5));


  TGraphErrors *grMC = new TGraphErrors();

  grMC->SetPoint(0,3,rebin_factor*hvz1MCfull->GetBinContent(1));
  grMC->SetPoint(1,9,rebin_factor*hvz1MCfull->GetBinContent(2));
  grMC->SetPoint(2,15,rebin_factor*hvz1MCfull->GetBinContent(3));
  grMC->SetPoint(3,21,rebin_factor*hvz1MCfull->GetBinContent(4));
  grMC->SetPoint(4,27,rebin_factor*hvz1MCfull->GetBinContent(5));

  grMC->SetPoint(5,38,rebin_factor*hvz2MCfull->GetBinContent(1));
  grMC->SetPoint(6,44,rebin_factor*hvz2MCfull->GetBinContent(2));
  grMC->SetPoint(7,50,rebin_factor*hvz2MCfull->GetBinContent(3));
  grMC->SetPoint(8,56,rebin_factor*hvz2MCfull->GetBinContent(4));
  grMC->SetPoint(9,62,rebin_factor*hvz2MCfull->GetBinContent(5));

  grMC->SetPoint(10,89,hvz3MCfull->GetBinContent(1));
  grMC->SetPoint(11,101,hvz3MCfull->GetBinContent(2));
  grMC->SetPoint(12,113,hvz3MCfull->GetBinContent(3));
  grMC->SetPoint(13,125,hvz3MCfull->GetBinContent(4));
  grMC->SetPoint(14,137,hvz3MCfull->GetBinContent(5));


  grMC->SetPoint(15,164,hvz4MCfull->GetBinContent(1));
  grMC->SetPoint(16,176,hvz4MCfull->GetBinContent(2));
  grMC->SetPoint(17,188,hvz4MCfull->GetBinContent(3));
  grMC->SetPoint(18,200,hvz4MCfull->GetBinContent(4));
  grMC->SetPoint(19,212,hvz4MCfull->GetBinContent(5));
  
  grMC->SetPoint(20,239,hvz5MCfull->GetBinContent(1));
  grMC->SetPoint(21,251,hvz5MCfull->GetBinContent(2));
  grMC->SetPoint(22,263,hvz5MCfull->GetBinContent(3));
  grMC->SetPoint(23,275,hvz5MCfull->GetBinContent(4));
  grMC->SetPoint(24,287,hvz5MCfull->GetBinContent(5));
  
  grMC->SetPoint(25,314,hvz6MCfull->GetBinContent(1));
  grMC->SetPoint(26,326,hvz6MCfull->GetBinContent(2));
  grMC->SetPoint(27,338,hvz6MCfull->GetBinContent(3));
  grMC->SetPoint(28,350,hvz6MCfull->GetBinContent(4));
  grMC->SetPoint(29,362,hvz6MCfull->GetBinContent(5));


  TGraphErrors *grwb = new TGraphErrors();

  grwb->SetPoint(0,3,rebin_factor*hvz1wb->GetBinContent(1));
  grwb->SetPoint(1,9,rebin_factor*hvz1wb->GetBinContent(2));
  grwb->SetPoint(2,15,rebin_factor*hvz1wb->GetBinContent(3));
  grwb->SetPoint(3,21,rebin_factor*hvz1wb->GetBinContent(4));
  grwb->SetPoint(4,27,rebin_factor*hvz1wb->GetBinContent(5));

  grwb->SetPoint(5,38,rebin_factor*hvz2wb->GetBinContent(1));
  grwb->SetPoint(6,44,rebin_factor*hvz2wb->GetBinContent(2));
  grwb->SetPoint(7,50,rebin_factor*hvz2wb->GetBinContent(3));
  grwb->SetPoint(8,56,rebin_factor*hvz2wb->GetBinContent(4));
  grwb->SetPoint(9,62,rebin_factor*hvz2wb->GetBinContent(5));

  grwb->SetPoint(10,89,hvz3wb->GetBinContent(1));
  grwb->SetPoint(11,101,hvz3wb->GetBinContent(2));
  grwb->SetPoint(12,113,hvz3wb->GetBinContent(3));
  grwb->SetPoint(13,125,hvz3wb->GetBinContent(4));
  grwb->SetPoint(14,137,hvz3wb->GetBinContent(5));


  grwb->SetPoint(15,164,hvz4wb->GetBinContent(1));
  grwb->SetPoint(16,176,hvz4wb->GetBinContent(2));
  grwb->SetPoint(17,188,hvz4wb->GetBinContent(3));
  grwb->SetPoint(18,200,hvz4wb->GetBinContent(4));
  grwb->SetPoint(19,212,hvz4wb->GetBinContent(5));
  
  grwb->SetPoint(20,239,hvz5wb->GetBinContent(1));
  grwb->SetPoint(21,251,hvz5wb->GetBinContent(2));
  grwb->SetPoint(22,263,hvz5wb->GetBinContent(3));
  grwb->SetPoint(23,275,hvz5wb->GetBinContent(4));
  grwb->SetPoint(24,287,hvz5wb->GetBinContent(5));
  
  grwb->SetPoint(25,314,hvz6wb->GetBinContent(1));
  grwb->SetPoint(26,326,hvz6wb->GetBinContent(2));
  grwb->SetPoint(27,338,hvz6wb->GetBinContent(3));
  grwb->SetPoint(28,350,hvz6wb->GetBinContent(4));
  grwb->SetPoint(29,362,hvz6wb->GetBinContent(5));


  float xgraph[30] = {3,9,15,21,27,38,44,50,56,62,89,101,113,125,137,164,176,188,200,212,239,251,263,275,287,314,326,338,350,362};

  float yMC[30]= {};
  float yDT[30]= {};
  
  yMC[0] = rebin_factor*hvz1MCfull->GetBinContent(1);
  yMC[1] = rebin_factor*hvz1MCfull->GetBinContent(2);
  yMC[2] = rebin_factor*hvz1MCfull->GetBinContent(3);
  yMC[3] = rebin_factor*hvz1MCfull->GetBinContent(4);
  yMC[4] = rebin_factor*hvz1MCfull->GetBinContent(5);
  yMC[5] = rebin_factor*hvz2MCfull->GetBinContent(1);
  yMC[6] = rebin_factor*hvz2MCfull->GetBinContent(2);
  yMC[7] = rebin_factor*hvz2MCfull->GetBinContent(3);
  yMC[8] = rebin_factor*hvz2MCfull->GetBinContent(4);
  yMC[9] = rebin_factor*hvz2MCfull->GetBinContent(5);
  yMC[10] = hvz3MCfull->GetBinContent(1);
  yMC[11] = hvz3MCfull->GetBinContent(2);
  yMC[12] = hvz3MCfull->GetBinContent(3);
  yMC[13] = hvz3MCfull->GetBinContent(4);
  yMC[14] = hvz3MCfull->GetBinContent(5);
  yMC[15] = hvz4MCfull->GetBinContent(1);
  yMC[16] = hvz4MCfull->GetBinContent(2);
  yMC[17] = hvz4MCfull->GetBinContent(3);
  yMC[18] = hvz4MCfull->GetBinContent(4);
  yMC[19] = hvz4MCfull->GetBinContent(5);
  yMC[20] = hvz5MCfull->GetBinContent(1);
  yMC[21] = hvz5MCfull->GetBinContent(2);
  yMC[22] = hvz5MCfull->GetBinContent(3);
  yMC[23] = hvz5MCfull->GetBinContent(4);
  yMC[24] = hvz5MCfull->GetBinContent(5);
  yMC[25] = hvz6MCfull->GetBinContent(1);
  yMC[26] = hvz6MCfull->GetBinContent(2);
  yMC[27] = hvz6MCfull->GetBinContent(3);
  yMC[28] = hvz6MCfull->GetBinContent(4);
  yMC[29] = hvz6MCfull->GetBinContent(5);

  yDT[0] = rebin_factor*hvz1wb->GetBinContent(1);
  yDT[1] = rebin_factor*hvz1wb->GetBinContent(2);
  yDT[2] = rebin_factor*hvz1wb->GetBinContent(3);
  yDT[3] = rebin_factor*hvz1wb->GetBinContent(4);
  yDT[4] = rebin_factor*hvz1wb->GetBinContent(5);
  yDT[5] = rebin_factor*hvz2wb->GetBinContent(1);
  yDT[6] = rebin_factor*hvz2wb->GetBinContent(2);
  yDT[7] = rebin_factor*hvz2wb->GetBinContent(3);
  yDT[8] = rebin_factor*hvz2wb->GetBinContent(4);
  yDT[9] = rebin_factor*hvz2wb->GetBinContent(5);
  yDT[10] = hvz3wb->GetBinContent(1);
  yDT[11] = hvz3wb->GetBinContent(2);
  yDT[12] = hvz3wb->GetBinContent(3);
  yDT[13] = hvz3wb->GetBinContent(4);
  yDT[14] = hvz3wb->GetBinContent(5);
  yDT[15] = hvz4wb->GetBinContent(1);
  yDT[16] = hvz4wb->GetBinContent(2);
  yDT[17] = hvz4wb->GetBinContent(3);
  yDT[18] = hvz4wb->GetBinContent(4);
  yDT[19] = hvz4wb->GetBinContent(5);
  yDT[20] = hvz5wb->GetBinContent(1);
  yDT[21] = hvz5wb->GetBinContent(2);
  yDT[22] = hvz5wb->GetBinContent(3);
  yDT[23] = hvz5wb->GetBinContent(4);
  yDT[24] = hvz5wb->GetBinContent(5);
  yDT[25] = hvz6wb->GetBinContent(1);
  yDT[26] = hvz6wb->GetBinContent(2);
  yDT[27] = hvz6wb->GetBinContent(3);
  yDT[28] = hvz6wb->GetBinContent(4);
  yDT[29] = hvz6wb->GetBinContent(5);

  
  


  ofstream out(Form("%serrors_filips.dat",dir.Data()));

  for(int i=0;i<30;i++){

    out << i << " " << xgraph[i] << " " << yMC[i] << " " << err_MC[i] << " " << yDT[i] << " " << 4*err_DT[i] << endl;

  }

  out.close();
  

  // ERRORS
  for(int i=0;i<30;i++){
    if(i<10)grMCpr->SetPointError(i,0,2*err_MC_pr[i]);
    else grMCpr->SetPointError(i,0,err_MC_pr[i]);
  }

   // ERRORS

  for(int i=0;i<30;i++){
    if(i<10)grMChd->SetPointError(i,0,2*err_MC_hd[i]);
    else grMChd->SetPointError(i,0,err_MC_hd[i]);
  }

  // ERRORS

  for(int i=0;i<30;i++){
    if(i<10)grMC->SetPointError(i,0,2*err_MC[i]);
    else grMC->SetPointError(i,0,err_MC[i]);
  }
 

  // ERRORS

  for(int i=0;i<30;i++){
    if(i<10)grwb->SetPointError(i,0,2*err_DT[i]);
    else grwb->SetPointError(i,0,err_DT[i]);
  }

  
  

  gStyle->SetOptFit(1111);
  
  grMCpr->SetLineWidth(2);
  grMCpr->SetLineColor(4);
  //grMCpr->SetMarkerStyle(21);
  grMCpr->SetMarkerColor(4);
  //grMCpr->Draw("APL && same");
  grMCpr->Fit("expo");
  grMCpr->SetTitle("MC protons");

 /*TF1 *hd_fit = new TF1("hd_fit","[0]*TMath::Log([1]*x)",0,365);
  hd_fit->SetParameter(0, 180);
  hd_fit->SetParameter(1, 0.79);
 // hd_fit->SetParameter(2, -170);*/

  TF1 *hd_fit = new TF1("hd_fit","pol3",0,365);
  hd_fit->SetParameter(0, 155);
  hd_fit->SetParameter(1, 23);
  hd_fit->SetParameter(2, -0.109);
  hd_fit->SetParameter(3, -0.00016);
  
  grMChd->SetLineWidth(2);
  grMChd->SetLineColor(3);
  //grMChd->SetMarkerStyle(22);
  grMChd->SetMarkerColor(3);
  //grMChd->Fit("expo");
  grMChd->Fit("hd_fit");
  grMChd->SetTitle("MC hadrons");

/*  TF1 *fit = new TF1("fit","expo(0) + [2]*TMath::Log([3]*x)",0,365);
  fit->SetParameter(0, grMCpr->GetFunction("expo")->GetParameter(0));
  fit->SetParameter(1, grMCpr->GetFunction("expo")->GetParameter(1));;
  fit->FixParameter(2, hd_fit->GetParameter(0));
  fit->FixParameter(3, hd_fit->GetParameter(1));
*/
  //fit->SetParLimits(3, -1.62,-1.22);
  //fit->FixParameter(4, -170);

  TF1 *fit = new TF1("fit","expo(0) + pol3(2)",0,365);
  fit->SetParameter(0, grMCpr->GetFunction("expo")->GetParameter(0));
  fit->SetParameter(1, grMCpr->GetFunction("expo")->GetParameter(1));;
  fit->FixParameter(2, hd_fit->GetParameter(0));
  //fit->SetParameter(2,155);
  fit->FixParameter(3,hd_fit->GetParameter(1));
  fit->FixParameter(4,hd_fit->GetParameter(2));
  fit->FixParameter(5,hd_fit->GetParameter(3));

  grMC->SetLineWidth(2);
  grMC->SetLineColor(1);
  grMC->SetMarkerStyle(20);
  grMC->SetMarkerColor(1);
  //grMChd->Fit("expo");
  grMC->Fit(fit);
  grMC->SetTitle("Full MC");

  //fit->SetParameter(0, 7.8);
  //fit->SetParameter(1, -0.0058);
  //fit->SetParLimits(2, 135,145);
  //fit->SetParameter(2,155);
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
  mgr->GetYaxis()->SetRangeUser(0,6000);
  mgr->GetXaxis()->SetLimits(0,365);
  mgr->GetXaxis()->SetTitle("vz [mm]");
  mgr->GetYaxis()->SetTitle("FULL DATA");
  mgr->Draw("AP");

  
  /*TCanvas *c1 = new TCanvas();
  grMCpr->GetYaxis()->SetRangeUser(0,3500);
  grMCpr->GetXaxis()->SetLimits(0,365);
  grMCpr->GetXaxis()->SetTitle("vz [mm]");
  grMCpr->GetYaxis()->SetTitle("Proton vertices");
  grMCpr->Draw("APB");
  gr->Draw("P && same");
  //grwb->Draw("P && sames");
  */
  TCanvas *c2 = new TCanvas();
  grMC->GetYaxis()->SetRangeUser(0,6000);
  grMC->GetXaxis()->SetLimits(0,365);
  grMC->GetXaxis()->SetTitle("vz [mm]");
  grMC->GetYaxis()->SetTitle("Selected vertices");
  grMC->Draw("AP");
  grMCpr->Draw("P && same");
  grMChd->Draw("P && same");
  // grwb->Draw("P && same");

log_file.close();
  
}
