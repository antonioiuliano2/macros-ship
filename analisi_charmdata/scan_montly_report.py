#Just discover that pandas can import excel, let's try it (26 February 2019)
import pandas as pd
from matplotlib import pyplot as plt
from pandas.plotting import register_matplotlib_converters
import matplotlib.dates as mdates

register_matplotlib_converters() #for datetimes
#importing Excel worksheets
months = ['January','February']
for month in months:

    df = pd.read_excel("/home/antonio/Dropbox/Charmdata/monthly_scan_report.xls",sheet_name=month)

#getting lists
    date = df[u'Date']
    mic2 = df[u'Mic2']
    mic3 = df[u'Mic3']
    total = df[u'Total']
#doing the plots
    fig, ax = plt.subplots()
    plt.suptitle('Monthly scaning report ' + month)
    ax.plot(date,mic3,'.r')
    ax.plot(date,mic2,'+b')
    ax.plot(date,total,'*y')
    ax.legend(loc='upper right')
#setting labels for date
    ax.set_xticks(date)
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%d-%m"))
    ax.xaxis.set_minor_formatter(mdates.DateFormatter("%d-%m"))
    _=plt.xticks(rotation=90)
    plt.ylim(0,8)
    fig.show()