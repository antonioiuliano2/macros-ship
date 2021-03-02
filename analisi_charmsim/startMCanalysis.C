//perform Valerio's code for MC training
void startMCanalysis(){
  gROOT->ProcessLine(".L /afs/cern.ch/work/a/aiuliano/public/macros-ship/analisi_charmsim/vtx_MC_analysis.C+");
  gROOT->ProcessLine("Loop()");
}
