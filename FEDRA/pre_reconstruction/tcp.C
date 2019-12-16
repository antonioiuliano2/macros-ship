//TCut c_p0("c_p0","abs(s.eTY-0.11)<0.02&&abs(s.eTX-0.02)<0.02&&eN1==1&&eN2==1&&s1.eChi2+s2.eChi2+1>s.eChi2");
TCut c_pall("c_pall","eN1==1&&eN2==1");
TCut c_p0("c_p0","abs(s.eTY-0.1)<0.02&&abs(s.eTX-0.07)<0.02&&eN1==1&&eN2==1");

void tcp()
{
  // draw_peak(c_p0);
  //draw_peak_xy(c_p0);
  draw_distortion(c_pall);
}

void draw_distortion(TCut &cut)
{
  TCanvas *c = new TCanvas( Form("dist_%s",cut.GetName()), cut.GetName(),1200,1000 );
  //defining histograms with axis labels
  TH2F *s0_xy = new TH2F("s0_xy","XY positions of segments;x[#mum];y[#mum]",120,0,120000,100,0,100000);
  
  TProfile2D *s1tx_xy = new TProfile2D("s1tx_xy","Residual map (s1.eTX - s.eTX);x[#mum];y[#mum]",120,0,120000,100,0,100000);
  TProfile2D *s1ty_xy = new TProfile2D("s1ty_xy","Residual map (s1.eTY - s.eTY);x[#mum];y[#mum]",120,0,120000,100,0,100000);
  TProfile2D *s2tx_xy = new TProfile2D("s2tx_xy","Residual map (s2.eTX - s.eTX);x[#mum];y[#mum]",120,0,120000,100,0,100000);
  TProfile2D *s2ty_xy = new TProfile2D("s2ty_xy","Residual map (s2.eTY - s.eTY);x[#mum];y[#mum]",120,0,120000,100,0,100000);
  
  TProfile2D *chi2p_xy = new TProfile2D("chi2p_xy","Chi square map;x[#mum];y[#mum]",120,0,120000,100,0,100000);
  c->Divide(3,2);

  c->cd(1);
  couples->Draw("s.eY:s.eX>>s0_xy", cut ,"colz");
  s0_xy->Draw("COLZ");
  
  c->cd(2);
  couples->Draw("s1.eTX-s.eTX:s.eY:s.eX>>s1tx_xy", cut ,"prof colz");
  s1tx_xy->Draw("COLZ");
  c->cd(3);
  couples->Draw("s1.eTY-s.eTY:s.eY:s.eX>>s1ty_xy", cut ,"prof colz");
  s1ty_xy->Draw("COLZ");
  c->cd(4);
  couples->Draw("eCHI2P:s.eY:s.eX>>chi2p_xy", cut ,"prof colz");
  chi2p_xy->Draw("COLZ");
  c->cd(5);
  couples->Draw("s2.eTX-s.eTX:s.eY:s.eX>>s2tx_xy", cut ,"prof colz");
  s2tx_xy->Draw("COLZ");
  c->cd(6);
  couples->Draw("s2.eTY-s.eTY:s.eY:s.eX>>s2ty_xy", cut ,"prof colz0");
  s2ty_xy->Draw("COLZ");
}

void draw_peak(TCut &cut)
{
  TCanvas *c = new TCanvas( cut.GetName(), cut.GetName(),1200,1000 );
  
  c->Divide(3,3);

  c->cd(1);
  couples->Draw("s.eTY:s.eTX", cut ,"");
  c->cd(2);
  couples->Draw("s1.eTY:s1.eTX", cut ,"");
  c->cd(3);
  couples->Draw("s2.eTY:s2.eTX", cut ,"");
  
  c->cd(4);
  couples->Draw("s.eY:s.eX", cut ,"");
  
  c->cd(5);
  couples->Draw("s1.eChi2+s2.eChi2:s.eChi2", cut ,"");
  
  c->cd(6);
  couples->Draw("s1.eTX", cut ,"");
  c->cd(9);
  couples->Draw("s2.eTX", cut ,"");
  
  c->cd(7);
  couples->Draw("s.eTX", cut ,"");
  
  c->cd(8);
  couples->Draw("s.eTY", cut ,"");
  
  
}

void draw_peak_xy(TCut &cut)
{
  TCanvas *c1 = new TCanvas( cut.GetName(), cut.GetName(),1200,1000 );
  
  c1->Divide(4,3);

  c1->cd(1);
  couples->Draw("s.eTX:s.eX", cut ,"");
  c1->cd(5);
  couples->Draw("s1.eTX:s.eX", cut ,"");
  c1->cd(9);
  couples->Draw("s2.eTX:s.eX", cut ,"");
 
  c1->cd(2);
  couples->Draw("s.eTX:s.eY", cut ,"");
  c1->cd(6);
  couples->Draw("s1.eTX:s.eY", cut ,"");
  c1->cd(10);
  couples->Draw("s2.eTX:s.eY", cut ,"");
  
  c1->cd(3);
  couples->Draw("s.eTY:s.eX", cut ,"");
  c1->cd(7);
  couples->Draw("s1.eTY:s.eX", cut ,"");
  c1->cd(11);
  couples->Draw("s2.eTY:s.eX", cut ,""); 

  c1->cd(4);
  couples->Draw("s.eTY:s.eY", cut ,"");
  c1->cd(8);
  couples->Draw("s1.eTY:s.eY", cut ,"");
  c1->cd(12);
  couples->Draw("s2.eTY:s.eY", cut ,"");
   
}
