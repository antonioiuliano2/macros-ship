'''test if points are correctly stored'''
import matplotlib.pyplot as plt
import numpy as np
filipsfile = open("errors_filips.dat","r")

npoints = 26
x = np.zeros(npoints)
yMCpr =  np.zeros(npoints)#MonteCarlo protons
yMChd =  np.zeros(npoints)#MonteCarlo hadrons
yMC =  np.zeros(npoints)#MonteCarlo

errMC = np.zeros(npoints)
errMCpr = np.zeros(npoints)
errMChd = np.zeros(npoints)

yDT =  np.zeros(npoints) #data
errDT =  np.zeros(npoints)
#start reading file
for iline, line in enumerate(filipsfile):
 if(iline == 0): #headers
      continue
 points = line.split()
 ipoint = int(points[0])
 
 x[ipoint] = float(points[1])

 yMCpr[ipoint] = float(points[2])
 errMCpr[ipoint] = float(points[3])

 yMChd[ipoint] = float(points[4])
 errMChd[ipoint] = float(points[5])

 yMC[ipoint] = float(points[6])
 errMC[ipoint] = float(points[7])

 yDT[ipoint] = float(points[8])
 errDT[ipoint] = float(points[9])

print("end loop in file, plotting")
fig = plt.figure()
plt.errorbar(x,yMCpr,yerr=errMCpr,fmt="ro")
plt.errorbar(x,yMChd,yerr=errMChd,fmt="bo")
plt.errorbar(x,yMC,yerr=errMC,fmt="ko")
plt.errorbar(x,yDT,yerr=errDT,fmt="go")

plt.show()