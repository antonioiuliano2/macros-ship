import ROOT as r

def openP8cascade(nfiles_p8):
 xrootpublic = "root://eospublic.cern.ch/"
 P8CascadePath = xrootpublic+"/eos/experiment/ship/simulation/bkg/Cascade2026/charm/cascade/"

 P8CascadeChain = r.TChain("pythia6")

 for ifile in range(nfiles_p8):
    P8CascadeChain.AddFile(P8CascadePath+"cascade_msel4_W_FTFT_{}.root".format(ifile))
 print("Total Number of entries in Pythia8 {} input files is {}".format(nfiles_p8,P8CascadeChain.GetEntries()))
 p8_df = r.RDataFrame(P8CascadeChain)
 return p8_df

def openP6cascade(nfiles_p6):
 xrootuser = "root://eoshome-a.cern.ch/"
 P6CascadePath = xrootuser+"/eos/user/a/aiuliano/public/sims_FairShip/cascadeproduction_Wtarget_2026/"

 P6CascadeChain = r.TChain("pythia6")
 for ifile in range(nfiles_p6):
    P6CascadeChain.AddFile(P6CascadePath+"Run_{}_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root".format(ifile))
 print("Total Number of entries in Pythia6 {} input files is {}".format(nfiles_p6,P6CascadeChain.GetEntries()))
 p6_df = r.RDataFrame(P6CascadeChain)
 return p6_df

def process_charm_dataframe(df, outputfilename):
 nbins = 400
 Emin = 0
 Emax = 400

 hE_all = df.Histo1D(("hE_all","Energy of all charmed hadrons",nbins,0,Emax),"E")
 
 df_dpm = df.Filter("TMath::Abs(id)==411")
 hE_dpm = df_dpm.Histo1D(("hE_dpm","Energy of charged mesons",nbins,0,Emax),"E")

 df_d0 = df.Filter("TMath::Abs(id)==421")
 hE_d0 = df_d0.Histo1D(("hE_d0","Energy of D0 mesons",nbins,0,Emax),"E")

 df_ds = df.Filter("TMath::Abs(id)==431")
 hE_ds = df_ds.Histo1D(("hE_ds","Energy of Ds mesons",nbins,0,Emax),"E")

 df_lambdac = df.Filter("TMath::Abs(id)==4122")
 hE_lambdac = df_lambdac.Histo1D(("hE_lambdac","Energy of Lambdac baryons",nbins,0,Emax),"E")

 outputfile = r.TFile(outputfilename,"RECREATE")

 hE_all.Write()
 hE_d0.Write()
 hE_ds.Write()
 hE_dpm.Write()
 hE_lambdac.Write()

 outputfile.Close()
#Pythia8
#nfiles_p8 = 215
#p8_df = openP8cascade(nfiles_p8)
#process_charm_dataframe(p8_df,"histos_P8_{}_files.root".format(nfiles_p8))
#Pythia6
#nfiles_p6 = 1
#p6_df = openP6cascade(nfiles_p6)
#process_charm_dataframe(p6_df,"histos_P6_{}_files.root".format(nfiles_p6))

#read histofiles, compare plots

p6histos = r.TFile.Open("histos_P6_1_files.root")
p8histos = r.TFile.Open("histos_P8_215_files.root")

def plot_comparison(histoname):

 c = r.TCanvas("c_"+histoname)
 hE_p6 = p6histos.Get(histoname)
 hE_p8 = p8histos.Get(histoname)

 if histoname != "hE_all":
  #fraction over all simulated charmed hadrons
  hE_p6.Scale(1./p6histos.Get("hE_all").Integral())
  hE_p8.Scale(1./p8histos.Get("hE_all").Integral())

 hE_p6.SetName(hE_p8.GetName()+"_pythia6")
 hE_p8.SetName(hE_p8.GetName()+"_pythia8")
 hE_p6.SetTitle("Pythia6 Charm Cascade Production")
 hE_p8.SetTitle("Pythia8 Charm Cascade Production")

 hE_p8.SetLineColor(r.kRed)
 hE_p6.Draw("histo")
 hE_p8.Draw("histo&&SAMES")
 c.BuildLegend()
 hE_p6.SetTitle("Charm Cascade Comparison "+histoname+";E[GeV];fraction over all charmed hadrons [%]")
 c.Draw()

r.gStyle.SetOptStat("nemri")
plot_comparison("hE_all")
plot_comparison("hE_d0")
plot_comparison("hE_dpm")
plot_comparison("hE_ds")
plot_comparison("hE_lambdac")