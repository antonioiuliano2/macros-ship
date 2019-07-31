'''tutorial su come disegnare un grafico da un file di testo con colonne
   lancialo con python -i miofile.dat
'''

import ROOT
import sys #per fornire il file direttamente da linea di comando

mygraph = ROOT.TGraph(sys.argv[1]) #sostituire TGraph con TGraphErrors se ci sono altre due colonne con gli errori

def setlabels():
    mygraph.SetTitle("My data")
    mygraph.GetXaxis().SetTitle("x")
    mygraph.GetYaxis().SetTitle("y")

setlabels()
mygraph.Draw("AP*")