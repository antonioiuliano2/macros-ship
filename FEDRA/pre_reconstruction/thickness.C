struct HXY{
  int   xb,yb;
  float xmin,xmax,ymin,ymax;
};

HXY hxy;

//float xmin=4000, ymin=4000;
float xb=202,yb=182;

HXY hxy1={xb,yb,
  -155000,-155000+770*xb,
  4000,   4000+yb*851
};
        
HXY hxy2={xb,yb,
  4000,   4000+770*xb,
  4000,   4000+yb*851
};
        
HXY hxy3={xb,yb,
  4000,4000+770*xb,
  4000,4000+yb*565
};
        
HXY hxy4={xb,yb,
  -155000,-155000+770*xb,
  4000,4000+yb*565
};

void thickness(int q)
{
  if(q==1) hxy=hxy1;
  if(q==2) hxy=hxy2;
  if(q==3) hxy=hxy3;
  if(q==4) hxy=hxy4;
  //draw_all();
  //draw_diff();
  draw_diff6();
  printf("\n %d \n", hxy.xb);
}

draw_diff6()
{
  gROOT->SetBatch(1);
  new TCanvas();
  Views->Draw( Form("eZ1:eYview:eXview>>hz1(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ1>1","prof colz");
  Views->Draw( Form("eZ2:eYview:eXview>>hz2(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ2>1","prof colz");
  Views->Draw( Form("eZ3:eYview:eXview>>hz3(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ3>1","prof colz");
  Views->Draw( Form("eZ4:eYview:eXview>>hz4(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ4>1","prof colz");
  gROOT->SetBatch(0);
  
  TH2D *hdz_up = hz1->ProjectionXY("hdz_up");  // this "Projection" is necessary to make Add working correctly
  TH2D *jz2    = hz2->ProjectionXY("jz2");
  hdz_up->Add(jz2,-1);
  
  TH2D *hdz_down = hz3->ProjectionXY("hdz_down");  // this "Projection" is necessary to make Add working correctly
  TH2D *jz4    = hz4->ProjectionXY("jz4");
  hdz_down->Add(jz4,-1);
  
  TH2D *hdz_base = hz2->ProjectionXY("hdz_base");  // this "Projection" is necessary to make Add working correctly
  TH2D *jz3 = hz3->ProjectionXY("jz3");
  hdz_base->Add(jz3,-1);
  
  //printf("Mean thickness: %6.2f %6.2f %6.2f \n",hdz_up->GetMean(3),hdz_base->GetMean(3),hdz_down->GetMean(3));
  
  TCanvas *c = new TCanvas("diff", gDirectory->GetFile()->GetName(), 1200,700);
  c->Divide(3,2);
  //const char *opt ="lego2 z";
  const char *opt ="colz";
  gStyle->SetOptStat("ne");
  
  //c->cd(1)->SetGrid();  hz1->Draw(opt);
  //c->cd(2)->SetGrid();  hz2->Draw(opt);
  //c->cd(3)->SetGrid();  hz3->Draw(opt);
  //c->cd(4)->SetGrid();  hz4->Draw(opt);

  c->cd(4)->SetGrid();  hdz_down->SetMinimum(10);  hdz_down->SetMaximum(100); hdz_down->Draw(opt);
  c->cd(5)->SetGrid();  hdz_up->SetMinimum(10);    hdz_up->SetMaximum(100);   hdz_up->Draw(opt);
  c->cd(6)->SetGrid();  hdz_base->SetMinimum(150); hdz_base->SetMaximum(180); hdz_base->Draw(opt);
  
  TH1F *su = Spectrum(*hdz_up,"up");
  TH1F *sd = Spectrum(*hdz_down,"down");
  TH1F *sb = Spectrum(*hdz_base,"base");
  c->cd(3)->SetGrid();
  sb->SetLineColor(kBlack);
  sb->Draw("h");
  su->SetLineColor(kRed);
  su->Draw("h same");
  sd->SetLineColor(kBlue);
  sd->Draw("h same");
  
  /*
  c->cd(5)->SetGrid();
  Views->Draw("eNcl:eZframe>>h(500,100,600)");
  Views->Draw("eZ3","","same");
  Views->Draw("eZ2","","same");
  
  c->cd(6)->SetGrid();
  Views->SetLineColor(kBlue); Views->Draw("eZ3>>hzlayer(500,100,600)","","");
  Views->SetLineColor(kBlack); Views->Draw("eNcl:eZframe","","same");
  Views->SetLineColor(kRed); Views->Draw("eZ1","","same");
  Views->SetLineColor(kRed); Views->Draw("eZ2","","same");
  Views->SetLineColor(kBlue); Views->Draw("eZ3","","same");
  Views->SetLineColor(kBlue); Views->Draw("eZ4","","same");
  Views->SetLineColor(kBlack);
  */
  c->cd(1)->SetGrid();
  Views->Draw( Form("eNsegments:eYview:eXview>>hxytop(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eNframesTop==0","prof colz");
  c->cd(2)->SetGrid();
  Views->Draw( Form("eNsegments:eYview:eXview>>hxybot(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eNframesTop!=0","prof colz");
 
}

draw_diff()
{
  gROOT->SetBatch(1);
  new TCanvas();
  Views->Draw( Form("eZ1:eYview:eXview>>hz1(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ1>1","prof colz");
  Views->Draw( Form("eZ2:eYview:eXview>>hz2(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ2>1","prof colz");
  Views->Draw( Form("eZ3:eYview:eXview>>hz3(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ3>1","prof colz");
  Views->Draw( Form("eZ4:eYview:eXview>>hz4(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eZ4>1","prof colz");
  gROOT->SetBatch(0);
  
  TH2D *hdz_up = hz1->ProjectionXY("hdz_up");  // this "Projection" is necessary to make Add working correctly
  TH2D *jz2    = hz2->ProjectionXY("jz2");
  hdz_up->Add(jz2,-1);
  
  TH2D *hdz_down = hz3->ProjectionXY("hdz_down");  // this "Projection" is necessary to make Add working correctly
  TH2D *jz4    = hz4->ProjectionXY("jz4");
  hdz_down->Add(jz4,-1);
  
  TH2D *hdz_base = hz2->ProjectionXY("hdz_base");  // this "Projection" is necessary to make Add working correctly
  TH2D *jz3 = hz3->ProjectionXY("jz3");
  hdz_base->Add(jz3,-1);
  
  //printf("Mean thickness: %6.2f %6.2f %6.2f \n",hdz_up->GetMean(3),hdz_base->GetMean(3),hdz_down->GetMean(3));
  
  TCanvas *c = new TCanvas("diff", gDirectory->GetFile()->GetName(), 1200,700);
  c->Divide(4,2);
  //const char *opt ="lego2 z";
  const char *opt ="colz";
  gStyle->SetOptStat("ne");
  
  //c->cd(1)->SetGrid();  hz1->Draw(opt);
  //c->cd(2)->SetGrid();  hz2->Draw(opt);
  //c->cd(3)->SetGrid();  hz3->Draw(opt);
  //c->cd(4)->SetGrid();  hz4->Draw(opt);

  c->cd(1)->SetGrid();  hdz_up->SetMinimum(30);    hdz_up->SetMaximum(100);   hdz_up->Draw(opt);
  c->cd(2)->SetGrid();  hdz_down->SetMinimum(30);  hdz_down->SetMaximum(100); hdz_down->Draw(opt);
  c->cd(3)->SetGrid();  hdz_base->SetMinimum(150); hdz_base->SetMaximum(190); hdz_base->Draw(opt);
  
  TH1F *su = Spectrum(*hdz_up,"up");
  TH1F *sd = Spectrum(*hdz_down,"down");
  TH1F *sb = Spectrum(*hdz_base,"base");
  c->cd(4)->SetGrid();
  sb->SetLineColor(kBlack);
  sb->Draw("h");
  su->SetLineColor(kRed);
  su->Draw("h same");
  sd->SetLineColor(kBlue);
  sd->Draw("h same");
  
  c->cd(5)->SetGrid();
  Views->Draw("eNcl:eZframe>>h(500,100,600)");
  Views->Draw("eZ3","","same");
  Views->Draw("eZ2","","same");
  
  c->cd(6)->SetGrid();
  Views->SetLineColor(kBlue); Views->Draw("eZ3>>hzlayer(500,100,600)","","");
  Views->SetLineColor(kBlack); Views->Draw("eNcl:eZframe","","same");
  Views->SetLineColor(kRed); Views->Draw("eZ1","","same");
  Views->SetLineColor(kRed); Views->Draw("eZ2","","same");
  Views->SetLineColor(kBlue); Views->Draw("eZ3","","same");
  Views->SetLineColor(kBlue); Views->Draw("eZ4","","same");
  Views->SetLineColor(kBlack);
   
  c->cd(7)->SetGrid();
  Views->Draw( Form("eNsegments:eYview:eXview>>hxytop(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eNframesTop==0","prof colz");
  c->cd(8)->SetGrid();
  Views->Draw( Form("eNsegments:eYview:eXview>>hxybot(%d,%f,%f,%d,%f,%f)",hxy.xb,hxy.xmin,hxy.xmax,hxy.yb,hxy.ymin,hxy.ymax),"eNframesTop!=0","prof colz");
 
}

draw_all()
{
  TCanvas *c = new TCanvas("entr","entr",900,700);
  c->Divide(2,2);
  
  c->cd(1)->SetGrid();
  Views->Draw("eZ1:Entry$");
  Views->Draw("eZ2:Entry$","", "same");
  Views->Draw("eZ3:Entry$","", "same");
  Views->Draw("eZ4:Entry$","", "same");
  
  if( gDirectory->Get("PinViews") )
  {
    c->cd(2)->SetGrid();
    PinViews->SetMarkerStyle(6);
    PinViews->Draw("eZ1:Entry$");
    PinViews->Draw("eZ2:Entry$","", "same");
    PinViews->Draw("eZ3:Entry$","", "same");
    PinViews->Draw("eZ4:Entry$","", "same");
  }
  
  c->cd(3)->SetGrid();
  Views->Draw("eZ1:eXview");
  Views->Draw("eZ2:eXview","", "same");
  Views->Draw("eZ3:eXview","", "same");
  Views->Draw("eZ4:eXview","", "same");
 
  c->cd(4)->SetGrid();
  Views->Draw("eZ1:eYview");
  Views->Draw("eZ2:eYview","", "same");
  Views->Draw("eZ3:eYview","", "same");
  Views->Draw("eZ4:eYview","", "same");
 }

TH1F *Spectrum(TH2D &h2, const char *name)
{
  TH1F *h = new TH1F(name,"Spectrum",300,0,300);
  for(int i=1; i<h2.GetNbinsX()+1; i++ )
    for(int j=1; j<h2.GetNbinsY()+1; j++ )
      h->Fill( h2.GetBinContent(i,j) );
  return h;
}
