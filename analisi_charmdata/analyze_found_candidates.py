import ROOT as r
import fedrarootlogon
import sys
from IPython import get_ipython
ipython = get_ipython()

'''Script to read output of decay search from Valerio, find charm candidates
usage: ipython -i analyze_found_candidates.py ds_file.root vertexfile.root outputhistofile.root
'''

#decay search from Valerio
dsfile = r.TFile.Open(sys.argv[1])
dstree = dsfile.Get("ds")

#reconstructed vertextree
vtxfile = r.TFile.Open(sys.argv[2])
vtxtree = vtxfile.Get("vtx")

#which vertices I want to study
candidateselection =  "dsvtx.nfound>=2"

dstree.Draw(">>candidates", candidateselection)
dstree.Draw("")
candidates = r.gDirectory.GetList().FindObject("candidates")
dstree.SetAlias("vid","vtx.fe_id")
dstree.SetAlias("vtx2id","dsvtx.vtx2_vid")
dstree.SetAlias("vtx2_tracks","dsvtx.vtx2_tid")
ncandidates = dstree.GetEntries(candidateselection)


outputhistofile = r.TFile(sys.argv[3],"RECREATE")
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
 
 if vertices.size() > 0.:
  for secondaryvertexid in vertices:
   vtxtree.GetEntry(secondaryvertexid)
   hmolt.Fill(vtxtree.n)
   hxy.Fill(vtxtree.vx, vtxtree.vy)
   hz.Fill(vtxtree.vz)
  ngoodcandidates = ngoodcandidates + 1
  candidatesentries.append(entr)
  print ("Candidates events {}: the primary vertex id is {}, nsecondaries: {}".format(entr,primaryvid,vertices.size()))

print ("Found candidates {}".format(ngoodcandidates))
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
def event(ievent,firstevent = False):
  if not firstevent:
   close() #close previous event
  dstree.GetEntry(ievent)
  primaryvid = dstree.vid
  print ("now inspecting vertex {}".format(primaryvid))
  vtxtree.GetEntry(primaryvid)

  prim_vx = vtxtree.vx

  vertices = dstree.vtx2id
  if vertices.size() > 0.:
   if vertices.size()==1:
    ipython.magic("run -i ~/Scrivania/macros-ship/FEDRA/VertexTrackDisplay.py -f vertextree_secondquarter.root -t linked_tracks.root -nv {} {} -new".format(primaryvid,vertices[0]))
   elif vertices.size()==2:
    ipython.magic("run -i ~/Scrivania/macros-ship/FEDRA/VertexTrackDisplay.py -f vertextree_secondquarter.root -t linked_tracks.root -nv {} {} {} -new".format(primaryvid,*vertices))
   elif vertices.size()==3:
    ipython.magic("run -i ~/Scrivania/macros-ship/FEDRA/VertexTrackDisplay.py -f vertextree_secondquarter.root -t linked_tracks.root -nv {} {} {} {} -new".format(primaryvid,*vertices))    
   elif vertice.size()==4:
    ipython.magic("run -i ~/Scrivania/macros-ship/FEDRA/VertexTrackDisplay.py -f vertextree_secondquarter.root -t linked_tracks.root -nv {} {} {} {} {} -new".format(primaryvid,*vertices))
  else:
   print ("No secondary vertex")

event(candidatesentries[icandidate],True)

def nextcandidate():
  global icandidate
  icandidate = icandidate + 1
  event(candidatesentries[icandidate])
