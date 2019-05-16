#Just discover that pandas can import excel, let's try it (26 February 2019)
import os #to import my $HOME PATH
import pandas as pd
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters
import matplotlib.dates as mdates

homepath = os.environ['HOME']

register_matplotlib_converters() #for datetimes

df = pd.read_excel(homepath+"/Dropbox/Archivio_cronologico/Maggio_2019/16May_SHiP_all.xlsx")

#getting lists
df.index = pd.to_datetime(df['Date'])
dfweekly = df.resample('W').sum()
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

#reporting the shifter hard work
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