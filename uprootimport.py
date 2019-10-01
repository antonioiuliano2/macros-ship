#Testing branch uprooting, to load them into keras, created on 26 September 2019


import uproot

mytree = uproot.open("vertices_firstquarter.root")["vtx"]

probability = mytree.array("probability")
molteplicity = mytree.array("n")

TX = mytree.array("TX")
TY = mytree.array("TY")

