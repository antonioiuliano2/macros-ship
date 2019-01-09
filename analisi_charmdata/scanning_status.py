from __future__ import division
'''Report status of scanning of SHIP-charm emulsions. How many emulsions are left to scan?'''
#build lists where booleans will be contained. An array of booleans for the emulsions and a string for location (Naples or Zurich)
bricks = []
#CHARM1 RUNS
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Zurich'])
bricks.append([[False] * 29, 'Napoli'])
#CHARM2 RUNS
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Zurich'])
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Napoli'])
bricks.append([[False] * 29, 'Napoli'])
#CHARM3 RUNS
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])
#CHARM4 RUNS
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])
#CHARM5 RUNS
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])
#CHARM6 RUNS
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])
bricks.append([[False] * 57, 'Napoli'])

#how many emulsions we need to scan?
ntotalemulsions = 0 #total emulsion to be scanned
ntotalemulsionsNapoli = 0
ntotalemulsionsZurich = 0
for brick in bricks:
    ntotalemulsions += len(brick[0])
    if (brick[1] == 'Napoli'):        
     ntotalemulsionsNapoli += len(brick[0])
    else:
     ntotalemulsionsZurich += len(brick[0])
print 'Total number of emulsion exposed in July 2018: ',ntotalemulsions

#Creating the dictionary with our naming convention
run = {'ch1r1':bricks[0],'ch1r2':bricks[1],'ch1r3':bricks[2],'ch1r4':bricks[3],'ch1r5':bricks[4],'ch1r6':bricks[5],\
 'ch2r1':bricks[6],'ch2r2':bricks[7],'ch2r3':bricks[8],'ch2r4':bricks[9],'ch2r5':bricks[10],'ch2r6':bricks[11],\
 'ch3r1':bricks[12],'ch3r2':bricks[13],'ch3r3':bricks[14],'ch4r1':bricks[15],'ch4r2':bricks[16],'ch4r3':bricks[17],\
 'ch5r1':bricks[18],'ch5r2':bricks[19],'ch5r3':bricks[20],'ch6r1':bricks[21],'ch6r2':bricks[22],'ch6r3':bricks[23]}

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

allscanned('ch1r4')
allscanned('ch1r6')
allscanned('ch2r2')
allscanned('ch2r3')
allscanned('ch2r4')
allscanned('ch2r5')

for index in range(22):
    scanned('ch2r1',index+1)
for index in range(10):
    scanned('ch2r6',index+1)
for index in range(15):
    scanned('ch3r2',index+1)

# ********************************************REPORT******************************************************#    
#counting how many emulsions are left to scan
nscannedemulsions = 0
nscannedemulsionsNapoli = 0
nscannedemulsionsZurich = 0


for brick in bricks:
    for emulsion in brick[0]:
        if (emulsion): 
            nscannedemulsions = nscannedemulsions + 1
            if (brick[1] == 'Napoli'):
             nscannedemulsionsNapoli = nscannedemulsionsNapoli + 1
            if (brick[1] == 'Zurich'):
             nscannedemulsionsZurich = nscannedemulsionsZurich + 1

print ('Total number of scanned emulsion so far: {}. Still to be scanned: {}'.format(nscannedemulsions, ntotalemulsions - nscannedemulsions))
print ('We have accomplished: {:.2%} of the scanning workload'.format(nscannedemulsions/ntotalemulsions)) #.2% is the format with percentage with 2 decimals afterwards (like .2f)

def localreport():
   print ('Total number of scanned emulsion so far in Naples: {}. Still to be scanned: {}'.format(nscannedemulsionsNapoli, ntotalemulsionsNapoli - nscannedemulsionsNapoli))
   print ('We have accomplished: {:.2%} of the Naples scanning workload'.format(nscannedemulsionsNapoli/ntotalemulsionsNapoli)) #.2% is the format with percentage with 2 decimals afterwards (like .2f) 