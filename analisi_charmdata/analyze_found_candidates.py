import ROOT as r
import fedrarootlogon
import os
from argparse import ArgumentParser 
from IPython import get_ipython
ipython = get_ipython()

'''Script to read output of decay search from Valerio, find charm candidates
usage: ipython -i analyze_found_candidates.py ds_file.root vertexfile.root outputhistofile.root
'''

parser = ArgumentParser()
parser.add_argument("-f", "--fedra", dest="vertexfilename", help="file with fedra vertices",
                    required=True)
parser.add_argument("-t", "--tracks", dest="tracksfilename", help="file with fedra tracks",
                    required=True)
parser.add_argument("-d", "--decaysearch", dest="decaysearchfilename", help="decay search output file",
                    required=True)
parser.add_argument("-o", "--outputfile", dest="histofile", help="output file with histogram",default=None)

options = parser.parse_args()
#decay search from Valerio
dsfile = r.TFile.Open(options.decaysearchfilename)
dstree = dsfile.Get("ds")

#reconstructed vertextree
vtxfile = r.TFile.Open(options.vertexfilename)
vtxtree = vtxfile.Get("vtx")

#which vertices I want to study
candidateselection =  "dsvtx.nfound>=2"

dstree.Draw(">>candidates", candidateselection)
dstree.Draw("")
candidates = r.gDirectory.GetList().FindObject("candidates")
ncandidates = dstree.GetEntries(candidateselection)

#setting aliases, python does not likes points
dstree.SetAlias("vid","vtx.fe_id")
dstree.SetAlias("vtx2id","dsvtx.vtx2_vid")
dstree.SetAlias("vtx2_tracks","dsvtx.vtx2_tid")
dstree.SetAlias("oneprongs","dsvtx.prg1_id")


if options.histofile is not None:
 outputhistofile = r.TFile(options.histofile,"RECREATE")
 hmolt = r.TH1I("hmolt","molteplicity of secondary;ntracks",10,0,10)
 hxy = r.TH2D("hxy","xy position of primary vertices;x[#mum];y[#mum]",120,0,120000,100,0,100000)
 hz = r.TH1D("hz","z position of primary vertices;z[#mum]",30,-30000,0,)

candidatesentries = []
ngoodcandidates = 0
for icandidate in range(ncandidates):
 entr = candidates.GetEntry(icandidate)
 #check how many secondaries
 dstree.GetEntry(entr)
 primaryvid = dstree.vid
 vertices = dstree.vtx2id
 oneprongs = dstree.oneprongs
 if vertices.size() > 0.:
  for secondaryvertexid in vertices:
   vtxtree.GetEntry(secondaryvertexid)
   if options.histofile is not None:
    hmolt.Fill(vtxtree.n)
    hxy.Fill(vtxtree.vx, vtxtree.vy)
    hz.Fill(vtxtree.vz)
  ngoodcandidates = ngoodcandidates + 1
  candidatesentries.append(entr)
  print ("Candidates events {}: the primary vertex id is {}, nsecondaries: {}, oneprongtracks: {}".format(entr,primaryvid,vertices.size(),oneprongs.size()))
icandidate = 0
print ("Found candidates {}".format(ngoodcandidates))
if options.histofile is not None:
 cmolt = r.TCanvas()
 hmolt.Draw()
 icandidate = 0
 cxy = r.TCanvas()
 cxy.Divide(1,2)
 cxy.cd(1)
 hxy.Draw()
 cxy.cd(2)
 hz.Draw()

 hmolt.Write()
 hxy.Write()
 hz.Write()

macrospath=os.getenv("MACROSSHIP")
def event(ievent,firstevent = False):
  if not firstevent:
   close() #close previous event
  dstree.GetEntry(ievent)
  primaryvid = dstree.vid
  print ("now inspecting vertex {}".format(primaryvid))
  vtxtree.GetEntry(primaryvid)

  prim_vx = vtxtree.vx

  vertices = dstree.vtx2id
  oneprongs = dstree.oneprongs
  if vertices.size() > 0.:
   vertexids="-nv {}"
   tracksids=" "
   for ivtx in vertices: #adding a {} for each secondary vertex to draw
    vertexids= vertexids +" {}"
   if oneprongs.size() > 0.: 
    tracksids=tracksids+"-nt {}"
   # for trackid in oneprongs:
    # tracksids=tracksids+" {}"
  
    vertices.insert(vertices.end(), oneprongs.begin(), oneprongs.end() );   
    ipython.magic(("run -i "+macrospath+\
    "/FEDRA/VertexTrackDisplay.py -f {} -t {} "+\
     vertexids+tracksids+" -new").format(options.vertexfilename,options.tracksfilename,primaryvid,*vertices))
   else:
  #  ipython.magic(("run -i "+macrospath+\
  #  "/FEDRA/VertexTrackDisplay.py -f {} -t {} "+\
  #   vertexids+tracksids+" -new").format(options.vertexfilename,options.tracksfilename,primaryvid,*vertices)) 
    ipython.magic(("run -i "+macrospath+\
     "/FEDRA/VertexTrackDisplay.py -f {} -t {} "+\
                  vertexids+tracksids+" -new").format(options.vertexfilename,options.tracksfilename,primaryvid,*vertices)) 
  else:
   print ("No secondary vertex")

event(candidatesentries[0],True)

def nextcandidate():
  global icandidate
  icandidate = icandidate + 1
  event(candidatesentries[icandidate])
