#try to load Thomas tree as a RDataFrame
import ROOT as r

r.gInterpreter.Declare(	
"""	
 #include <map>
 int pdg2value(int id){
  std::map<int, int> mymap;
    mymap [421] = 1;
    mymap [-421] = -1;
    mymap [411] = 2;
    mymap [-411] = -2;
    mymap [431] = 3;
    mymap [-431] = -3;
    mymap [4122] = 4;
    mymap [-4122] = -4;
    mymap [4132] = 5;
    mymap [-4132] = -5;
    mymap [4232] = 6;
    mymap [-4232] = -6;
    mymap [4332] = 7;
    mymap [-4332] = -7;



  std::map<int,int>::iterator it;
  it = mymap.find(id);
  if (it == mymap.end()) return 0;
  else return mymap[id];
 }
""")
df = r.RDataFrame("pythia6", "/home/utente/Scrivania/SHIPBuild/only_charm_Lead.root") #build the dataframe (file with only charm hadrons)

df_new = df.Define("P","sqrt(pow(px,2) + pow(py,2) + pow(pz,2))")
df_new2 = df.Define("name","pdg2value(id)")

df_primary = df_new.Filter('k<2', 'Looking for primary charm production') #count how many protons are primary
df_cascade = df_new.Filter('k>=2', 'Looking for charm cascade production')

report = df_primary.Report() #builds the report

report.Print() #checks the report

#drawing the histograms for primary vs cascade
cpdg = r.TCanvas()
hpdg = df_new2.Histo1D(("hpdg","Produced charm hadron",15, -7,8), "name")
hpdg.Draw()

pdgaxis = hpdg.GetXaxis()
#Setting labels
pdg = r.TDatabasePDG.Instance()
pdgaxis.SetBinLabel(1,pdg.GetParticle(-4332).GetName())
pdgaxis.SetBinLabel(2,pdg.GetParticle(-4232).GetName())
pdgaxis.SetBinLabel(3,pdg.GetParticle(-4132).GetName())
pdgaxis.SetBinLabel(4,pdg.GetParticle(-4122).GetName())
pdgaxis.SetBinLabel(5,pdg.GetParticle(-431).GetName())
pdgaxis.SetBinLabel(6,pdg.GetParticle(-411).GetName())
pdgaxis.SetBinLabel(7,pdg.GetParticle(-421).GetName())
pdgaxis.SetBinLabel(9,pdg.GetParticle(421).GetName())
pdgaxis.SetBinLabel(10,pdg.GetParticle(411).GetName())
pdgaxis.SetBinLabel(11,pdg.GetParticle(431).GetName())
pdgaxis.SetBinLabel(12,pdg.GetParticle(4122).GetName())
pdgaxis.SetBinLabel(13,pdg.GetParticle(4132).GetName())
pdgaxis.SetBinLabel(14,pdg.GetParticle(4232).GetName())
pdgaxis.SetBinLabel(15,pdg.GetParticle(4332).GetName())

c = r.TCanvas()
l = r.TLegend(0.2, 0.2, .8, .8)
hprimary_momentum = df_primary.Histo1D(("hp_primary","Momentum of charm from primary production",400,0,400),"P")
hcascade_momentum = df_cascade.Histo1D(("hp_cascade","Momentum of charm from cascade production",400,0,400),"P")

hcascade_momentum.Draw()
hcascade_momentum.SetTitle("")
hcascade_momentum.GetXaxis().SetTitle("GeV/c")
hprimary_momentum.Draw("SAMES")
hcascade_momentum.SetLineColor(r.kRed)

l.AddEntry('hp_primary','Momentum of charm hadrons from primary production')
l.AddEntry('hp_cascade','Momentum of charm hadrons from cascade production')
l.Draw()
c.Print('charm_momenta.png','png')
