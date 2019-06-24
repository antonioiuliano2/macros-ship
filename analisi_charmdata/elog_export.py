#Just discover that pandas can import excel, let's try it (26 February 2019)
import os #to import my $HOME PATH
import pandas as pd
import numpy as np
import ROOT
from root_numpy import fill_hist
from pandas.plotting import register_matplotlib_converters
import matplotlib.dates as mdates

homepath = os.environ['HOME']

register_matplotlib_converters() #for datetimes

df = pd.read_excel(homepath+"/Dropbox/Archivio_cronologico/Giugno_2019/elog_export_shipcharm_21_giugno.xlsx")

# ROOT histograms and Canvases are built outside functions, otherwise they are killed when function ends
croot = ROOT.TCanvas()

hbaseslavich = ROOT.TH1D("hbaseslavich", "Base plastic from Slavich emulsions", 13, 147.5, 212.5)
hbasenagoya = ROOT.TH1D("hbasenagoya", "Base plastic from Nagoya emulsions", 13, 147.5, 212.5)

croot1 = ROOT.TCanvas()

hemuslavich = ROOT.TH1D("hemuslavich", "Slavich emulsions", 14, 45, 115)
hemunagoya = ROOT.TH1D("hemunagoya", "Nagoya emulsions", 14, 45, 115)

def measuredthicknesses():
 ''' plot thicknesses of plastic base and emulsion layer '''
 croot.cd()
 dfslavich = df.query('EmuType == "Slavich"')
 dfnagoya = df.query('EmuType=="Nagoya"')
   
 fill_hist(hbaseslavich, dfslavich['Base'])
 fill_hist(hbasenagoya, dfnagoya['Base'])

 hbaseslavich.Draw()
 hbaseslavich.GetXaxis().SetTitle("thickness[#mum]")
 hbasenagoya.SetLineColor(ROOT.kRed)
 hbasenagoya.Draw("SAMES")

 croot.BuildLegend() # I trust you will build it correctly

 croot1.cd()

 emuthickslavich = np.concatenate((dfslavich['Bottom'],dfslavich['Top']),axis = 0) #merge arrays
 emuthicknagoya = np.concatenate((dfnagoya['Bottom'],dfnagoya['Top']),axis = 0) #merge arrays 

 fill_hist(hemuslavich, emuthickslavich)
 fill_hist(hemunagoya, emuthicknagoya)

 hemuslavich.Draw()
 hemuslavich.GetXaxis().SetTitle("thickness[#mum]")
 hemunagoya.SetLineColor(ROOT.kRed)
 hemunagoya.Draw("SAMES")

 croot1.BuildLegend() # I trust you will build it correctly
 '''Leaving this commented part as an example of how to do hist plots with matplotlib'''
 #fig1 = plt.figure()
 #plt.hist(emuthickslavich, bins=14 , range=(45,115) ,label = 'Slavich')
 #plt.hist(emuthicknagoya, bins=14 , range=(45,115) ,label = 'Nagoya')
 #plt.xlabel('$\mu$m') 
 #plt.title("Emulsion layer thickness")
 #plt.legend()
 #plt.show()

def scanningtime():
 from matplotlib import pyplot as plt #ROOT Canvases do not like matplotlib
 '''check number of plate scanned per mic in each microscope'''
 #getting lists
 df.index = pd.to_datetime(df['Date'])
 dfnomic5 = df.query('Microscope != "mic5"')
 dfweekly = dfnomic5.resample('W').sum()
 #date = dfweekly['Date']
 totalplates = dfweekly['Nplates']
 #print (dfweekly)

 dfmic3 = df.query('Microscope == "mic3"')
 dfweekly = dfmic3.resample('W').sum()
 mic3plates = dfweekly['Nplates']
 #print (dfweekly)

 dfmic2 = df.query('Microscope == "mic2"')
 dfweekly = dfmic2.resample('W').sum()
 mic2plates = dfweekly['Nplates']
 #print (dfweekly)
 #doing the plots, date is automatically displayed as the index of the dataframe
 fig, ax = plt.subplots()
 plt.suptitle('Monthly scan report ')
 ax.plot(mic3plates,'.r', markersize=12, label='mic3')
 ax.plot(mic2plates,'+b', markersize=12, label='mic2')
 ax.plot(totalplates,'*y',markersize=12, label='all')
 ax.legend(loc='upper left',fontsize='large') #borderpad makes the legend larger

 #plt.ylim(0,8)
 fig.show()

def shifters():
 from matplotlib import pyplot as plt #ROOT Canvases do not like matplotlib
 '''reporting the shifter hard work'''
 nAntonio = df['Author'].str.contains('Iuliano').value_counts()[True]
 nAntonia = df['Author'].str.contains('Crescenzo').value_counts()[True]
 nValerio = df['Author'].str.contains('Gentile').value_counts()[True]
 nGiuliana = df['Author'].str.contains('Galati').value_counts()[True]
 nAdele = df['Author'].str.contains('Lauria').value_counts()[True]
 nCristina = df['Author'].str.contains('Montesi').value_counts()[True]
 nBrenda = df['Author'].str.contains('Capone').value_counts()[True]
 nArtem = df['Author'].str.contains('Golovatiuk').value_counts()[True]

 ntotal = df.count().tolist()[0]

 ratios = [] #percentages to fill the pie plot
 ratios.append(float(nAntonio)/ntotal*100)
 ratios.append(float(nAntonia)/ntotal*100)
 ratios.append(float(nValerio)/ntotal*100)
 ratios.append(float(nGiuliana)/ntotal*100)
 ratios.append(float(nAdele)/ntotal*100)
 ratios.append(float(nCristina)/ntotal*100)
 ratios.append(float(nBrenda)/ntotal*100)
 ratios.append(float(nArtem)/ntotal*100)

 print ntotal
 print ratios

 #explode = [0., 0,] #we want to 'extract' the slice of already scanned emulsions 
 colors = ['g','m', 'y','r','b','c','r','m'] 
 labels = ['Antonio', 'Antonia','Valerio','Giuliana','Adele','Cristina','Brenda','Artem']
 fig1,ax1 = plt.subplots()
 ax1.pie(ratios, colors = colors, autopct='%1.1f%%', labels=labels, shadow=True, startangle=90)

 #plt.pie(ratios, explode=explode, colors = colors, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
 ax1.axis('equal') #ensures that the pie is drawn as a circle
 fig1.show()
