import math
from matplotlib import pyplot as plt
#inserting measurements
mover = [] #in the TM reference
survey_original = [] #in the survey reference
survey_transformed = [] #to use TM axis system

nmovements = 4

mover.append([0,0,0]) #units in mm
mover.append([150,0,0])
mover.append([150,100,0])
mover.append([0,100,0])
mover.append([0,0,0])

survey_original.append([9.9201,-0.0915,0.0752]) #units in m
survey_original.append([9.9196,0.0583,0.0754])
survey_original.append([9.9193,0.0590,-0.0245])
survey_original.append([9.9198,-0.0910, -0.0247])
survey_original.append([9.9200,-0.0915, 0.0752])

xmover=[]
ymover=[]
xsurvey_transformed = []
ysurvey_transformed = []
#xtm = ysurvey, ytm = -zsurvey, also using first measurement as origin
for step in survey_original:
	step[0] = (step[0] - 9.9201) *1000
	step[1] = (step[1] + 0.0915)*1000
	step[2] = (step[2] - 0.0752)*1000
	survey_transformed.append([step[1],-step[2],step[0]])
	xsurvey_transformed.append(step[1])
	ysurvey_transformed.append(-step[2])

for step in mover:
	xmover.append(step[0])
	ymover.append(step[1])

plt.figure()

plt.errorbar(xsurvey_transformed,ysurvey_transformed,0.5,0.5,'g.',label='measured')
plt.plot(xsurvey_transformed,ysurvey_transformed,'r',label='measured')
plt.plot(xmover,ymover,'b',label='expected')
plt.xlabel('x[mm]')
plt.ylabel('y[mm]')

def anglexy(index, reference = survey_transformed): 
	if (index>nmovements):
		print ('Error: no more than {} movements have been measured'.format(nmovements))
		1/0
	'''compute angle with respect to the x of the movement number index in the xy plane. Use mover for mover reference'''
	start = reference[index]
	end = reference[index+1]

	return math.atan2(end[1]-start[1],end[0]-start[0])

def xyz(index, reference = survey_transformed):
	if (index>nmovements):
		 print ('Error: no more than {} movements have been measured'.format(nmovements))
		 1/0	
	start = reference[index]
	end = reference[index+1]

	return (end[0]-start[0],end[1]-start[1],end[2]-start[2])


for i in range(nmovements):
 x,y,z = xyz(i)
 xTM,yTM,zTM = xyz(i,mover)
 deltax = x - xTM
 deltay = y -yTM

 print('{number:.3f}th movement: ({xstep:.3f},{ystep:.3f},{zstep:.3f}). Variation: ({deltaxstep:.3f}, {deltaystep:.3f},{zstep:.3f})'.format(number=i+1, xstep=x, ystep=y, zstep=z, deltaxstep=deltax, deltaystep=deltay ))
 print('xy angle during the {number:.4f} movement: {angle:.4f}. Angle variation:{deltatheta:.4f}'.format(number= i+1, angle=anglexy(i),deltatheta = anglexy(i) -anglexy(i,mover)))

 #showing the points where the TM passed
plt.legend(loc='upper center')
plt.show()