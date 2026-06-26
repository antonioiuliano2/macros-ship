import ROOT as r
import sys
import os

mbiasfile = str(sys.argv[1])
outFileName = str(sys.argv[2])
pot_number = int(sys.argv[3])
cascadefile = str(sys.argv[4])
pot_number_cascade = int(sys.argv[5])

norm_pot = pot_number #normalization pot, same as mbias by default

#Initialization, create gsimple ntuple and metadata ntuple
gsimple_entry = r.genie.flux.GSimpleNtpEntry()
meta_entry = r.genie.flux.GSimpleNtpMeta()
aux_entry = r.genie.flux.GSimpleNtpAux()

#To generate a random seed which is stored in the metadata. (has no meaning so far)
ran = r.TRandom3()
#outputfile setting branches
fOut = r.TFile.Open(outFileName, "RECREATE")
tOut = r.TTree("flux", "a simple flux n-tuple")
metaOut = r.TTree("meta", "metadata for flux n-tuple")

metakey = r.TString(outFileName).Hash()
metakey &= 0x7FFFFFFF

tOut.Branch("gsimple", gsimple_entry)
metaOut.Branch("meta", meta_entry)
tOut.Branch("aux", aux_entry)

# Scoring plane info
# Find a Minimum and Maximum of the x, y, z axis cm

min_z = 999
max_z = -999
# Max x, min x
min_x = 999
max_x = -999
#Max y, min y
min_y = 999
max_y = -999

pdglist = set() #extract all pdg codes of neutrinos in the scoring plane, to be used for metadata
min_weight = 1e10
max_weight = -1e10
max_energy = 0

charmExtern = [4332, 4232, 4132, 4232, 4122, 431, 411, 421, 15] #exclude charmed hadrons and tau leptons
neutrinos = [-12, 12, -14, 14, -16, 16]

def read_mbias_flux(filename, Emin = 1.):
    '''read mbias flux file
        filename: path to the mbias flux file
        Emin: minimum energy of neutrinos to be included in the output file
    '''
    global min_x, max_x, min_y, max_y, min_z, max_z
    global min_weight, max_weight, max_energy

    with open(filename) as f:
     for line in f:
      f = r.TFile.Open(os.environ["EOSSHIP"] + line.rstrip())  
      print("opened file ", line)
      sTree = f.Get("cbmsim")
      for n in range(sTree.GetEntries()):
        sTree.GetEntry(n)
        for v in sTree.PlaneHAPoint:
            nu = v.GetTrackID()
            t = sTree.MCTrack[nu]
            pdgcode = t.GetPdgCode()
            if abs(pdgcode) not in neutrinos:
                continue
            moID = abs(sTree.MCTrack[t.GetMotherId()].GetPdgCode())
            if moID in charmExtern:
                continue  # take heavy flavours from separate production
  
            px = v.GetPx()
            py = v.GetPy()
            pz = v.GetPz()

            E = r.TMath.Sqrt(px*px + py*py + pz*pz)
            if E < Emin:
                continue #skip neutrinos below Emin
            
            x = v.GetX()
            y = v.GetY()
            z = v.GetZ()
            
            gsimple_entry.Reset()
            aux_entry.Reset()
            gsimple_entry.metakey = metakey

            min_x = min(min_x, x)
            max_x = max(max_x, x)
            min_y = min(min_y, y)
            max_y = max(max_y, y)
            min_z = min(min_z, z)
            max_z = max(max_z, z)

            gsimple_entry.pdg = pdgcode
            gsimple_entry.wgt = norm_pot/pot_number
            gsimple_entry.vtxx = x * 1 / 100 # convert from cm to m
            gsimple_entry.vtxy = y * 1 / 100
            gsimple_entry.vtxz = z * 1 / 100
            gsimple_entry.dist = 0. # Distance from hadron decay point to neutrino "vertex", to use for oscillations,
                                      # for example. Don't use.
            gsimple_entry.px = px # in GeV/c
            gsimple_entry.py = py
            gsimple_entry.pz = pz
            gsimple_entry.E = E
            # Accumulate metadata
            pdglist.add(gsimple_entry.pdg) #add pdgcode to the set of pdgcodes
            
            min_weight = min(min_weight, gsimple_entry.wgt)
            max_weight = max(max_weight, gsimple_entry.wgt)
            max_energy = max(max_energy, gsimple_entry.E)
            # All done!
            tOut.Fill()

def read_cascadedecay_file(filename, Emin = 1.):
    '''read cascade decay file
        filename: path to the cascade decay file
        Emin: minimum energy of neutrinos to be included in the output file
    '''
    global min_x, max_x, min_y, max_y, min_z, max_z
    global min_weight, max_weight, max_energy

    f = r.TFile.Open(os.environ["EOSSHIP"] + filename)
    print("opened file ", filename)
    DecayTree = f.Get("Decay")

    for Entry in DecayTree:

        pdgcode = Entry.id
        if abs(pdgcode) not in neutrinos:
            continue

        if Entry.E < Emin:
            continue #skip neutrinos below Emin
        
        gsimple_entry.Reset()
        aux_entry.Reset()
        gsimple_entry.metakey = metakey

        gsimple_entry.pdg = pdgcode
        gsimple_entry.wgt = Entry.weight * norm_pot / pot_number_cascade #normalize to the same POT as the main file

        gsimple_entry.px = Entry.px # in GeV/c
        gsimple_entry.py = Entry.py
        gsimple_entry.pz = Entry.pz
        gsimple_entry.E = Entry.E

        tx = Entry.px / Entry.pz
        ty = Entry.py / Entry.pz

        zstart = 0; #starting point from zero;
        z = min_z;  #at the scoring plane (cm)
        x = tx * (z - zstart) #propagate to the scoring plane
        y = ty * (z - zstart)
    
        min_x = min(min_x, x)
        max_x = max(max_x, x)
        min_y = min(min_y, y)
        max_y = max(max_y, y)

        gsimple_entry.vtxx = x * 1 / 100 # convert from cm to m
        gsimple_entry.vtxy = y * 1 / 100
        gsimple_entry.vtxz = z * 1 / 100
        gsimple_entry.dist = 0. # Distance from hadron decay point to neutrino "vertex", to use for oscillations,
                                      # for example. Don't use.

        # Accumulate metadata
        pdglist.add(gsimple_entry.pdg) #add pdgcode to the set of pdgcodes
            
        min_weight = min(min_weight, gsimple_entry.wgt)
        max_weight = max(max_weight, gsimple_entry.wgt)
        max_energy = max(max_energy, gsimple_entry.E)
        # All done!
        tOut.Fill()

#reading the two input files and filling tree
read_mbias_flux(mbiasfile, Emin = 1.)
read_cascadedecay_file(cascadefile, Emin = 1.)    

# plane corner and plane direction of the scoring plane

plane_corner = [min_x, min_y, min_z]
plane_dir1 = [max_x - min_x, 0, max_z - min_z]
plane_dir2 = [0, max_y - min_y, max_z - min_z]

# sort out metadata

# copy pdg list
meta_entry.pdglist.clear()
for pdgcode in pdglist:
    meta_entry.pdglist.push_back(pdgcode)

meta_entry.maxEnergy = max_energy
meta_entry.maxWgt = max_weight
meta_entry.minWgt = min_weight

meta_entry.protons = pot_number # Number of protons on target

for i in range(3):
     meta_entry.windowBase[i] = plane_corner[i] * 1 / 100 # convert from cm to m
     meta_entry.windowDir1[i] = plane_dir1[i] * 1 / 100
     meta_entry.windowDir2[i] = plane_dir2[i] * 1 / 100

print("Start converting mbias production")
meta_entry.infiles.push_back(mbiasfile)
print("Start converting cascade production")
meta_entry.infiles.push_back(cascadefile)

meta_entry.seed = ran.GetSeed()
meta_entry.metakey = metakey

metaOut.Fill()

fOut.cd()
metaOut.Write()
tOut.Write()
fOut.Close()

print("Done converting")