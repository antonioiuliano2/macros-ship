#Script for computing the parameters needed to operate the target mover during July measurement 
#when not specified, units are in mm and s
import math
import matplotlib.pyplot as plt

Ax_rpss = 60
Vx_rpm = 360. #to be converted in mm/s

RpmtoSi = 4. / 60.
Ax = Ax_rpss * 4.
Vx = Vx_rpm * RpmtoSi

emux = 125.
emuy = 100.
beamwidth = 10. #FWHM of beam

ystep = beamwidth

targetenddx = 10 #end thickness of brick container
targetenddy = 10

nspills =int(math.ceil(emuy / ystep))
print(nspills)
#1) place the mover to a fixed distance from beam position
targettobeamx = 15
targettobeamy = 15

#2) set x0 and y0 to place the start of the emulsion at the beam center

tspill = 5
#start computing movement distance
#prespill of 
dxstart = Vx * 0.9 + Ax * 0.1 * 0.1 * 0.5

print "Acceleration length:", dxstart - Vx *0.9

#during a spill
xmovement = Vx * tspill
#after a spill
dxend = Vx * 0.9 + Ax * 0.1 * 0.1 * 0.5#as dx start

x0 = dxstart - (targetenddx + targettobeamx)
y0 = targettobeamy + targetenddy

print 'Posizione iniziale consigliata: ', x0, y0

totaldeltax = dxstart + xmovement + dxend

print 'Durante una spill percorre in x:', xmovement, "nel secondo prima e dopo una spill: ", dxend, 'mentre in totale percorre', totaldeltax

#starting movement pattern
xline = []
yline = []
x = 0.
y = 0.
for i in range(nspills):
 xline.append(x)
 yline.append(y)
 if (i % 2 == 0):
  x += totaldeltax
 else:
  x -= totaldeltax
 xline.append(x)
 yline.append(y)
 y += ystep
 
print(xline)
print(yline)

plt.plot(xline,yline,'r')
plt.show()
 

