'''test if points are correctly stored'''
import matplotlib.pyplot as plt
import numpy as np
filipsfile = open("errors_filips.dat","r")

npoints = 26
x = np.zeros(npoints)
yMC =  np.zeros(npoints)#MonteCarlo
errMC = np.zeros(npoints)
yDT =  np.zeros(npoints) #data
errDT =  np.zeros(npoints)
#start reading file
for line in filipsfile:
 points = line.split()
 ipoint = int(points[0])
 
 x[ipoint] = float(points[1])
 yMC[ipoint] = float(points[2])
 errMC[ipoint] = float(points[3])
 yDT[ipoint] = float(points[4])
 errDT[ipoint] = float(points[5])

print("end loop in file, plotting")
fig = plt.figure()
plt.errorbar(x,yMC,yerr=errMC,fmt="bo")
plt.errorbar(x,yDT,yerr=errDT,fmt="go")

plt.show()