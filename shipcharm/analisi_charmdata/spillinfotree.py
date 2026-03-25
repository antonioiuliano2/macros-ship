import numpy as np
import pandas as pd

df = pd.read_csv(r"/eos/experiment/ship/data/charmxsec/bookkeeping/spillinfo.csv")

print df["RunCode"]
print df["CYCLEID10"]
