//usage: rootprint -f png rootfile.root:hname -S statsoutoftheway.C -D COLZ.
//A simple and effectively way to plot histograms from root file

{
  gStyle->SetStatX(0.3);
  gStyle->SetStatY(0.9);  
  gStyle->SetOptStat(11111);
  
}
