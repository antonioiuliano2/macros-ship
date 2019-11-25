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

#df = pd.read_excel(homepath+"/Dropbox/Archivio_cronologico/Giugno_2019/elog_export_shipcharm_21_giugno.xlsx")
df = pd.read_csv("GSIscans_22_11_19.csv")
# ROOT histograms and Canvases are built outside functions, otherwise they are killed when function ends
croot = ROOT.TCanvas()

hbaseslavich = ROOT.TH1D("hbaseslavich", "Base plastic from Slavich emulsions;thickness[#mum]", 14, 147.5, 217.5)
hbasenagoya = ROOT.TH1D("hbasenagoya", "Base plastic from Nagoya emulsions;thickness[#mum]", 14, 147.5, 217.5)

croot1 = ROOT.TCanvas()

hemuslavich = ROOT.TH1D("hemuslavich", "Slavich emulsions;thickness[#mum]", 20, 45, 145)
hemunagoya = ROOT.TH1D("hemunagoya", "Nagoya emulsions;thickness[#mum]", 20, 45, 145)

def measuredthicknesses():
 ''' plot thicknesses of plastic base and emulsion layer '''
 croot.cd()
 dfslavich = df.query('EmuType == "Slavich"')
 dfnagoya = df.query('EmuType=="Nagoya"')
   
 fill_hist(hbaseslavich, dfslavich['Base'])
 fill_hist(hbasenagoya, dfnagoya['Base'])


 hbasenagoya.SetLineColor(ROOT.kRed)
 hbasenagoya.Draw()
 hbaseslavich.Draw("SAMES")
 croot.BuildLegend() # I trust you will build it correctly

 croot1.cd()

 emuthickslavich = np.concatenate((dfslavich['Bottom'],dfslavich['Top']),axis = 0) #merge arrays
 emuthicknagoya = np.concatenate((dfnagoya['Bottom'],dfnagoya['Top']),axis = 0) #merge arrays 

 fill_hist(hemuslavich, emuthickslavich)
 fill_hist(hemunagoya, emuthicknagoya)

 hemunagoya.SetLineColor(ROOT.kRed)
 hemunagoya.Draw()
 hemuslavich.Draw("SAMES")
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
 '''check number of plate scanned per mic in each microscope, need to edit time format with Excel beforehand'''
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
 Authors = ['Emil','Iuliano','Crescenzo','Gentile','Galati','Lauria','Montesi','Capone','Golovatiuk'] #as inserted by shifters, not whole name needs to be inserted
 labels = ['Emil','Antonio', 'Antonia','Valerio','Giuliana','Adele','Cristina','Brenda','Artem'] #names to be displayed in the plot
 colors = ['b','g','m', 'y','r','b','c','r','m'] 
 ratios = [] #percentages to fill the pie plot
 ntotal = len(df)

 for author in Authors:
     nscans = df['Author'].str.contains(author).value_counts()[True] 
     ratios.append(float(nscans)/ntotal*100)

 #sanity check
 if (sum(ratios)<100):
     print ("Warning: sum of scanning percentage lower than 100%: maybe missing some shifters?")
 elif (sum(ratios)>100):
     print ("Error: sum of scanning percentage higher than 100%: this is weird, stopping execution")
     return 2

 if (len(labels) != len(Authors)):
     print ("Error: size of labels of pie chart does not match size of shifter list")
     return 1
 if (len(colors) != len(Authors)):
     print ("Error: size of colors of pie chart does not match size of shifter list")
     return 1

 fig1,ax1 = plt.subplots()
 ax1.pie(ratios, colors = colors, autopct='%1.1f%%', labels=labels, shadow=True, startangle=90)

 #plt.pie(ratios, explode=explode, colors = colors, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
 ax1.axis('equal') #ensures that the pie is drawn as a circle
 fig1.show()

def processing():
    '''checking how many plates have been processed and linked'''
    dfGSI1 = df[df['Scanning Path'].str.contains("GSI1")]
    dfGSI2 = df[df['Scanning Path'].str.contains("GSI2")]
    dfGSI3 = df[df['Scanning Path'].str.contains("GSI3")]
    #removing RE-SCANS and taking only one value
    print(dfGSI2.sort_values("PlateID"))
    print("Scanned a total of {} films from GSI1, {} films from GSI2 and {} films from GSI3".format(len(dfGSI1),len(dfGSI2),len(dfGSI3)))
    dfGSI1 = dfGSI1[dfGSI1["Scanning"]!="TO RE-SCAN"]
    dfGSI2 = dfGSI2[dfGSI2["Scanning"]!="TO RE-SCAN"]
    dfGSI3 = dfGSI3[dfGSI3["Scanning"]!="TO RE-SCAN"]
  
    dfGSI1 = dfGSI1.drop_duplicates(subset=["PlateID"]) 
    dfGSI2 = dfGSI2.drop_duplicates(subset=["PlateID"]) 
    dfGSI3 = dfGSI3.drop_duplicates(subset=["PlateID"]) 
    print(dfGSI2.sort_values("PlateID"))
    print("Excluding duplicates: {} films from GSI1, {} films from GSI2 and {} films from GSI3".format(len(dfGSI1),len(dfGSI2),len(dfGSI3)))

    #plotting linking and alignment vs plateID

    import matplotlib.pyplot as plt
    plt.figure()
    plt.plot(dfGSI1["PlateID"],dfGSI1["Processing"],"b*",label = "GSI1")
    plt.plot(dfGSI2["PlateID"],dfGSI2["Processing"],"r*",label = "GSI2")
    plt.xlabel("PlateID")
    plt.ylabel("processing")
    plt.legend()

    plt.figure()
    plt.plot(dfGSI1["PlateID"],dfGSI1["Linking"],"b*",label = "GSI1")
    plt.plot(dfGSI2["PlateID"],dfGSI2["Linking"],"r*",label = "GSI2")
    plt.xlabel("PlateID")
    plt.ylabel("Linking")
    plt.legend()
    plt.show()