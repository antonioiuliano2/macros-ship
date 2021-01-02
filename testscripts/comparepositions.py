import ROOT as r
import fedrarootlogon

ali = r.EdbPVRec()
vrec = r.EdbVertexRec()
dproc = r.EdbDataProc()
scancond = r.EdbScanCond()
ali.SetScanCond(scancond)
#namespace
class NameSpace:
    def __init__(self, **kwargs):
        self.__dict__.update(kwargs)
VERTEX_PAR = NameSpace(DZmax = 3000.,ProbMinV=0.01,ImpMax=15.,UseMom=False,UseSegPar=True,QualityMode=0)

#setting VertexRec parameters
vrec.eDZmax=VERTEX_PAR.DZmax
vrec.eProbMin=VERTEX_PAR.ProbMinV
vrec.eImpMax=VERTEX_PAR.ImpMax
vrec.eUseMom=VERTEX_PAR.UseMom
vrec.eUseSegPar=VERTEX_PAR.UseSegPar
vrec.eQualityMode=VERTEX_PAR.QualityMode

vrec.SetPVRec(ali)
dproc.ReadVertexTree(vrec,"vertextree_test.root","1")

vertices = vrec.eVTX

inputfile = r.TFile.Open("vertextree_test.root")
inputtree = inputfile.Get("vtx")

inputtree.BuildIndex("vID")

hdR = r.TH1D("hdR","Distance dR of vertices, found vs expected;dR[#mu m]",100,0,1)
for vertex in vertices:
    vID = vertex.ID()
    inputtree.GetEntryWithIndex(vID)
    dX = vertex.VX() - inputtree.vx
    dY = vertex.VY() - inputtree.vy
    dZ = vertex.VZ() - inputtree.vz
    dR = r.TVector3(dX,dY,dZ).Mag()
    hdR.Fill(dR)
    if (dR > 10):
        print ("More than 10 microns! ", vID)

hdR.Draw()
