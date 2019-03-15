#Just discover that pandas can import excel, let's try it (26 February 2019)
import os #to import my $HOME PATH
import pandas as pd
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters
import matplotlib.dates as mdates

homepath = os.environ['HOME']

register_matplotlib_converters() #for datetimes
months = ['January','February','weekly']
for month in months:

    df = pd.read_excel(homepath+"/Dropbox/Charmdata/monthly_scan_report.xls",sheet_name=month)

#getting lists
    date = df[u'Date']
    mic2 = df[u'Mic2']
    mic3 = df[u'Mic3']
    total = df[u'Total']

    fig, ax = plt.subplots()
    plt.suptitle('Monthly scan report ' + month)
    ax.plot(date,mic3,'.r', markersize=12)
    ax.plot(date,mic2,'+b', markersize=12)
    ax.plot(date,total,'*y',markersize=12)
    ax.legend(loc='upper left',fontsize='large') #borderpad makes the legend larger

    ax.set_xticks(date)
    plt.xlabel('day')
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%m-%d-%y"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%m-%d-%y"))
    _=plt.xticks(rotation=90)
    #plt.ylim(0,8)
    fig.show()