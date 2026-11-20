import ROOT as r
import numpy as np
df = r.RDataFrame("flux","/eos/experiment/ship/user/aiuliano/inputs_gevgen_fnal/flux_gsimple_Hanae_makeCascade_merged_2026_06/Shuffled_Mbias2026noCharm_makeCascade_merged_nuflux_gsimple.root")
normpot = 6e+20
dfmeta = r.RDataFrame("meta","/eos/experiment/ship/user/aiuliano/inputs_gevgen_fnal/flux_gsimple_Hanae_makeCascade_merged_2026_06/Shuffled_Mbias2026noCharm_makeCascade_merged_nuflux_gsimple.root")

protons_np = dfmeta.AsNumpy(columns=["protons"])
pot = protons_np["protons"][0]
df_nu_tau = df.Filter("abs(pdg)==16")
hE = df_nu_tau.Histo1D("E","wgt")

c = r.TCanvas()
hE.Scale(normpot/pot)
hE.DrawClone("histo")
print("Number of tau neutrinos (and anti-neutrinos) in "," POT is: ",hE.Integral())
c.Draw()