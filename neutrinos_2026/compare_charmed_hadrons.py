import ROOT as r

def openP8cascade(nfiles_p8):
 xrootpublic = "root://eospublic.cern.ch/"
 P8CascadePath = xrootpublic+"/eos/experiment/ship/simulation/bkg/Cascade2026/charm/cascade/"

 P8CascadeChain = r.TChain("pythia6")

 for ifile in range(nfiles_p8):
    P8CascadeChain.AddFile(P8CascadePath+"cascade_msel4_W_FTFT_{}.root".format(ifile))
 p8_df = r.RDataFrame(P8CascadeChain)
 return p8_df

def openP6cascade(nfiles_p6):
 xrootuser = "root://eoshome-a.cern.ch/"
 P6CascadePath = xrootuser+"/eos/user/a/aiuliano/public/sims_FairShip/cascadeproduction_Wtarget_2026/"

 P6CascadeChain = r.TChain("pythia6")
 for ifile in range(nfiles_p6):
    P6CascadeChain.AddFile(P6CascadePath+"Run_{}_Cascade1000k-parp16-MSTP82-1-MSEL4-ntuple.root".format(ifile))
 p6_df = r.RDataFrame(P6CascadeChain)
 return p6_df

def process_charm_dataframe(df, outputfilename):
 nbins = 400
 Emin = 0
 Emax = 400

 df_dpm = df.Filter("TMath::Abs(id)==411")
 hE_dpm = df_dpm.Histo1D(("hE_dpm","Energy of charged mesons",nbins,0,Emax),"E")

 df_d0 = df.Filter("TMath::Abs(id)==421")
 hE_d0 = df_d0.Histo1D(("hE_d0","Energy of D0 mesons",nbins,0,Emax),"E")

 df_ds = df.Filter("TMath::Abs(id)==431")
 hE_ds = df_ds.Histo1D(("hE_ds","Energy of Ds mesons",nbins,0,Emax),"E")

 df_lambdac = df.Filter("TMath::Abs(id)==4122")
 hE_lambdac = df_lambdac.Histo1D(("hE_lambdac","Energy of Lambdac baryons",nbins,0,Emax),"E")

 outputfile = r.TFile(outputfilename,"RECREATE")

 hE_d0.Write()
 hE_ds.Write()
 hE_dpm.Write()
 hE_lambdac.Write()

 outputfile.Close()
#Pythia8
nfiles_p8 = 400
p8_df = openP8cascade(nfiles_p8)
process_charm_dataframe(p8_df,"histos_P8_{}_files.root".format(nfiles_p8))
Pythia6
nfiles_p6 = 1
p6_df = openP6cascade(nfiles_p6)
process_charm_dataframe(p6_df,"histos_P6_{}_files.root".format(nfiles_p6))