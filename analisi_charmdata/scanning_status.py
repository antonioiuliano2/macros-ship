from __future__ import division
from matplotlib import pyplot as plt
'''Report status of scanning of SHIP-charm emulsions. How many emulsions are left to scan?'''
#build lists where booleans will be contained. An array of booleans for the emulsions and a string for location (Naples or Zurich)
bricks = []
bademulsions = 0 #emulsion hard to scan due to problems
#CHARM1 RUNS
bricks.append([[False] * 29, 'Naples'])
bricks.append([[False] * 29, 'Naples'])
bricks.append([[False] * 29, 'Zurich'])
bricks.append([[False] * 29, 'Naples'])
bricks.append([[False] * 29, 'Zurich'])
bricks.append([[False] * 29, 'Naples'])
#CHARM2 RUNS
bricks.append([[False] * 29, 'Naples'])
bricks.append([[False] * 29, 'Zurich'])
bricks.append([[False] * 29, 'Naples'])
bricks.append([[False] * 29, 'Naples'])
bricks.append([[False] * 29, 'Naples'])
bricks.append([[False] * 29, 'Naples'])
#CHARM3 RUNS
bricks.append([[False] * 57, 'Zurich'])
bricks.append([[False] * 57, 'Naples'])
bademulsions = bademulsions + 2 # layer damaged by wrong oil
bricks.append([[False] * 57, 'Naples'])
bademulsions = bademulsions + 57 #Slavich development
#CHARM4 RUNS
bricks.append([[False] * 57, 'Zurich'])
bricks.append([[False] * 57, 'Naples'])
bricks.append([[False] * 57, 'Naples'])
bademulsions = bademulsions + 57 #Slavich development
#CHARM5 RUNS
bricks.append([[False] * 57, 'Naples'])
bricks.append([[False] * 57, 'Naples'])
bricks.append([[False] * 57, 'Naples'])
#CHARM6 RUNS
bricks.append([[False] * 57, 'Naples'])
bricks.append([[False] * 57, 'Naples'])
bricks.append([[False] * 57, 'Naples'])

#Creating the dictionary with our naming convention
run = {'ch1r1':bricks[0],'ch1r2':bricks[1],'ch1r3':bricks[2],'ch1r4':bricks[3],'ch1r5':bricks[4],'ch1r6':bricks[5],\
 'ch2r1':bricks[6],'ch2r2':bricks[7],'ch2r3':bricks[8],'ch2r4':bricks[9],'ch2r5':bricks[10],'ch2r6':bricks[11],\
 'ch3r1':bricks[12],'ch3r2':bricks[13],'ch3r3':bricks[14],'ch4r1':bricks[15],'ch4r2':bricks[16],'ch4r3':bricks[17],\
 'ch5r1':bricks[18],'ch5r2':bricks[19],'ch5r3':bricks[20],'ch6r1':bricks[21],'ch6r2':bricks[22],'ch6r3':bricks[23]}

 

#how many emulsions we need to scan?
ntotalemulsions = 0 #total emulsion to be scanned
ntotalemulsionsNaples = 0
ntotalemulsionsZurich = 0
for brick in bricks:
    ntotalemulsions += len(brick[0])
    if (brick[1] == 'Naples'):        
     ntotalemulsionsNaples += len(brick[0])
    else:
     ntotalemulsionsZurich += len(brick[0])
print 'Total number of emulsion exposed in July 2018: ',ntotalemulsions

# ******************************************Filling the scanned bricks************************************#
def allscanned(name):
    '''all the brick was scanned. Put all elements to True'''
    for i in range(len(run[name][0])):
         (run[name][0])[i] = True
def scanned(name, emulsion):
    '''Scanned emulsion of brick name'''    
    if (emulsion > len(run[name][0])):
        print ('Input emulsion over range: Number of emulsions in brick {} is {}'.format(name, len(run[name][0])))
        return None
    (run[name][0])[emulsion-1] = True

allscanned('ch1r1')
allscanned('ch1r2')
allscanned('ch1r4')
allscanned('ch1r6')

allscanned('ch2r1')
allscanned('ch2r2')
allscanned('ch2r3')
allscanned('ch2r4')
allscanned('ch2r5')
allscanned('ch2r6')

allscanned('ch3r2')
bademulsions = bademulsions + 2 #CH3-R2 P42 and P41 ruined by fixer
allscanned('ch3r3')
bademulsions = bademulsions + 57 #CH3-R3 developed Slavich

allscanned('ch4r2')
allscanned('ch4r3') #just because i need to subtract the bad ones, I know we have not scanned them yet
bademulsions = bademulsions + 57 #CH4-R3 developed Slavich

allscanned('ch5r2')
allscanned('ch5r3')

allscanned('ch6r1')
allscanned('ch6r2')
allscanned('ch6r3')

for index in range(30):
    scanned('ch5r1',index+1)

# ********************************************REPORT******************************************************#    
#counting how many emulsions are left to scan
nscannedemulsions = 0
nscannedemulsionsNaples = 0
nscannedemulsionsZurich = 0


for brick in bricks:
    for emulsion in brick[0]:
        if (emulsion): 
            nscannedemulsions = nscannedemulsions + 1
            if (brick[1] == 'Naples'):
             nscannedemulsionsNaples = nscannedemulsionsNaples + 1
            if (brick[1] == 'Zurich'):
             nscannedemulsionsZurich = nscannedemulsionsZurich + 1

print ('Total number of scanned emulsion so far: {}. Still to be scanned: {}'.format(nscannedemulsions, ntotalemulsions - nscannedemulsions))
print ('We have accomplished: {:.1%} of the scanning workload'.format(nscannedemulsions/ntotalemulsions)) #.2% is the format with percentage with 2 decimals afterwards (like .2f)

def generalreport():
   figure = plt.figure()
   generalratio = nscannedemulsions/ntotalemulsions
   #inserting elements of the graph
   ratios= []
   ratios.append(generalratio*100)
   ratios.append((1-generalratio)*100)

   explode = [0.1,0]
   colors = ['g','y']
   labels = ['Scanned', 'To be done']
   #ready to plot
   plt.pie(ratios, explode=explode, colors=colors, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
   plt.axis('equal')
   plt.show()
def localreport():
   figure = plt.figure()
   napoliratio = float(nscannedemulsionsNaples)/ntotalemulsionsNaples
   scanrate = 24 #assuming 5 emulsions scanned per day and 4 on Friday
   remainingweeks = (ntotalemulsionsNaples - nscannedemulsionsNaples)/scanrate
   print ('Total emulsion films in Naples: {}'.format(ntotalemulsionsNaples))
   print ('Total number of scanned emulsion so far in Naples: {}. Still to be scanned: {}'.format(nscannedemulsionsNaples, ntotalemulsionsNaples - nscannedemulsionsNaples))
   print ('We have accomplished: {:.1%} of the Naples scanning workload'.format(napoliratio)) #.2% is the format with percentage with 2 decimals afterwards (like .2f) 
   print ('Weeks required for completion: {:.2f}'.format(remainingweeks))
   #let's bake a cake
   #damaged = 2 #number of damaged emulsions
   ratios = [] #percentages to fill the pie plot
   ratios.append(napoliratio*100)
   #ratios.append(damaged/ntotalemulsionsNaples * 100)
   ratios.append((1- napoliratio)*100)

   explode = [0.1, 0] #we want to 'extract' the slice of already scanned emulsions 
   colors = ['g', 'y'] 
   labels = ['Scanned', 'To be done']

   plt.pie(ratios, explode=explode, colors = colors, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
   plt.axis('equal') #ensures that the pie is drawn as a circle
   #adding a red layer for bad emulsions
   figure2 = plt.figure()
   completeratios = [] 
   badratio = float(bademulsions)/ntotalemulsionsNaples*100
   goodratio = float(nscannedemulsionsNaples - bademulsions)/ntotalemulsionsNaples*100
   remainingratio = float(100 - badratio - goodratio)
   completeratios.append(goodratio)
   completeratios.append(badratio)
   completeratios.append(remainingratio)
   print "Test", goodratio, badratio, remainingratio
   explode = [0.1, 0, 0] #we want to 'extract' the slice of already scanned emulsions 
   colors = ['g', 'y', 'r']
   labels = ['Scanned normally', 'Less sensitive films', 'To be done'] 
   plt.pie(completeratios, explode=explode, colors = colors, labels=labels, autopct='%1.1f%%', shadow=True, startangle=90)
   plt.axis('equal') #ensures that the pie is drawn as a circle
   plt.show()

   
