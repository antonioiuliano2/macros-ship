import numpy as np
import pandas as pd
import sys
'''to replace topology as required by Valerio'''

dfall = pd.read_csv(sys.argv[1])
dfprimaryvertices =dfall[dfall["topology"]==1]
indexes = dfprimaryvertices.groupby("MCEventID")['ntracks'].idxmax()#returning indexes with maximum value
dfprimaryvertices = dfprimaryvertices.loc[indexes] #splicing to keep only the vertices with more tracks

for vtx in dfprimaryvertices["ivtx"]:
    dfall.loc[(dfall.ivtx == vtx),'topology']=1

dfall = dfall.sort_values("ivtx")
dfall = dfall.to_csv("edited_vtx_list.csv",index=False)

