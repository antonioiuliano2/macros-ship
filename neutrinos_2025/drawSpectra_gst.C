using namespace ROOT;
//draw neutrino spectra from gst file (A. Iuliano, 18 March 2026)
void drawSpectra_gst(){

 const double normsim = 4e+19; //reference of simulation weights (one year of data taking)
 const double normship = 5e+13; //reference for five years of data taking
 //open input file with RDataFrame
 RDataFrame df("gst","/home/antonio/Simulations/GENIE_sims/2026_03_16_1year_allflavours/nu_1year_fluxhanae34.gst.root");

 auto df1 = df.Filter("cc && dis"); //interaction type filter (CC DIS only)
 //define histograms
 auto hE_nu_e = df1.Filter("abs(neu)==12").Histo1D({"hE_nu_e","Electron neutrino energy spectrum;E[GeV]",400,0.,400.},"Ev");
 auto hE_nu_mu = df1.Filter("abs(neu)==14").Histo1D({"hE_nu_mu","Muon neutrino energy spectrum;E[GeV]",400,0.,400.},"Ev");
 auto hE_nu_tau = df1.Filter("abs(neu)==16").Histo1D({"hE_nu_tau","Tau neutrino energy spectrum;E[GeV]",400,0.,400.},"Ev");

 hE_nu_e->Scale(normship/normsim);
 hE_nu_mu->Scale(normship/normsim);
 hE_nu_tau->Scale(normship/normsim);

    
 TCanvas *c1 = new TCanvas();
 hE_nu_mu->SetLineColor(kRed);
 hE_nu_tau->SetLineColor(kGreen);

 hE_nu_mu->DrawClone("histo");
 hE_nu_e->DrawClone("histo&&SAMES");
 hE_nu_tau->DrawClone("histo&&SAMES");

 c1->BuildLegend();

}