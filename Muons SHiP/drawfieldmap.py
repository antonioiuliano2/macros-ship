import ROOT as r
import numpy

inputfile = r.TFile.Open('$FAIRSHIP/files/nuTauDetField.root')
tree = inputfile.Get('Data')
#defining graphs
outputfile = r.TFile('field_maps.root','RECREATE')
Bxfromxy = r.TGraph2D()
Bxfromxy.SetTitle("Bx map (xy);x[cm];y[cm];Bx[T]")
Bxfromxz = r.TGraph2D()
Bxfromxz.SetTitle("Bx map (xz);z[cm];x[cm];Bx[T]")
Bxfromyz = r.TGraph2D()
Bxfromyz.SetTitle("Bx map (yz);z[cm];y[cm];Bx[T]")

Byfromxy = r.TGraph2D()
Byfromxy.SetTitle("By map (xy);x[cm];y[cm];By[T]")
Byfromxz = r.TGraph2D()
Byfromxz.SetTitle("By map (xz);z[cm];x[cm];By[T]")
Byfromyz = r.TGraph2D()
Byfromyz.SetTitle("By map (yz);z[cm];y[cm];By[T]")

Bzfromxy = r.TGraph2D()
Bzfromxy.SetTitle("Bz map (xy);x[cm];y[cm];Bz[T]")
Bzfromxz = r.TGraph2D()
Bzfromxz.SetTitle("Bz map (xz);z[cm];x[cm];Bz[T]")
Bzfromyz = r.TGraph2D()
Bzfromyz.SetTitle("Bz map (yz);z[cm];y[cm];Bz[T]")


for i in range (tree.GetEntries()):
 tree.GetEntry(i)
 n = Byfromxy.GetN()
 #getting coordinates
 x = tree.x #beam direction
 y = tree.y #width
 z = tree.z #height
 #getting field values
 Bx = tree.Bx
 By = tree.By
 Bz = tree.Bz

 #saving info in graphs
 Bxfromxy.SetPoint(n, x, y, Bx)
 Bxfromxz.SetPoint(n, z, x, Bx)
 Bxfromyz.SetPoint(n, z, y, Bx)
 Byfromxy.SetPoint(n, x, y, By)
 Byfromxz.SetPoint(n, z, x, By)
 Byfromyz.SetPoint(n, z, y, By)
 Bzfromxy.SetPoint(n, x, y, Bz)
 Bzfromxz.SetPoint(n, z, x, Bz)
 Bzfromyz.SetPoint(n, z, y, Bz)

#drawing all the plots and saving them to file
cx = r.TCanvas()
cx.Divide(1,2)
cx.cd(1)
Bxfromxy.GetXaxis().SetTitle('x[cm]')
Bxfromxy.GetYaxis().SetTitle('y[cm]')
Bxfromxy.Draw('COLZ')
Bxfromxy.Write()
cx.Update()
cx.cd(2)
Bxfromxz.GetXaxis().SetTitle('z[cm]')
Bxfromxz.GetYaxis().SetTitle('x[cm]')
Bxfromxz.Draw('COLZ')
Bxfromxy.Write()
cx.Update()

cy = r.TCanvas()
cy.Divide(1,2)
cy.cd(1)
Byfromxy.GetXaxis().SetTitle('x[cm]')
Byfromxy.GetYaxis().SetTitle('y[cm]')
Byfromxy.Draw('COLZ')
Byfromxy.Write()
cy.Update()
cy.cd(2)
Byfromxz.GetXaxis().SetTitle('z[cm]')
Byfromxz.GetYaxis().SetTitle('x[cm]')
Byfromxz.Draw('COLZ')
Byfromxz.Write()
cy.Update()

cz = r.TCanvas()
cz.Divide(1,2)
cz.cd(1)
Bzfromxy.GetXaxis().SetTitle('x[cm]')
Bzfromxy.GetYaxis().SetTitle('y[cm]')
Bzfromxy.Draw('COLZ')
Bzfromxy.Write()
cz.Update()
cz.cd(2)
Bzfromxz.GetXaxis().SetTitle('z[cm]')
Bzfromxz.GetYaxis().SetTitle('x[cm]')
Bzfromxz.Draw('COLZ')
Bzfromxz.Write()
cz.Update()

cyz = r.TCanvas()
cyz.Divide(1,2)
cyz.cd(1)
Bxfromyz.GetXaxis().SetTitle('z[cm]')
Bxfromyz.GetYaxis().SetTitle('y[cm]')
Bxfromyz.Draw('COLZ')
Bxfromyz.Write()
cyz.Update()
cyz.cd(2)
Byfromyz.GetXaxis().SetTitle('z[cm]')
Byfromyz.GetYaxis().SetTitle('y[cm]')
Byfromyz.Draw('COLZ')
Byfromyz.Write()
cyz.Update()

cyz2 = r.TCanvas()
Bzfromyz.GetXaxis().SetTitle('z[cm]')
Bzfromyz.GetYaxis().SetTitle('y[cm]')
Bzfromyz.Draw('COLZ')
Bzfromyz.Write()
cyz2.Update()



outputfile.Close()
