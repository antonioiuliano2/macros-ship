import ROOT as r

inputfile = r.TFile.Open('$FAIRSHIP/files/nuTauDetField.root')
tree = inputfile.Get('Data')
#defining graphs

Bxfromxy = r.TGraph2D()
Bxfromxy.SetTitle("Bx map (xy);x[cm];y[cm];Bx[T]")
Bxfromxz = r.TGraph2D()
Bxfromxz.SetTitle("Bx map (xz);x[cm];z[cm];Bx[T]")
Bxfromyz = r.TGraph2D()
Bxfromyz.SetTitle("Bx map (yz);y[cm];z[cm];Bx[T]")

Byfromxy = r.TGraph2D()
Byfromxy.SetTitle("By map (xy);x[cm];y[cm];By[T]")
Byfromxz = r.TGraph2D()
Byfromxz.SetTitle("By map (xz);x[cm];z[cm];By[T]")
Byfromyz = r.TGraph2D()
Byfromyz.SetTitle("By map (yz);y[cm];z[cm];By[T]")

Bzfromxy = r.TGraph2D()
Bzfromxy.SetTitle("Bz map (xy);x[cm];y[cm];Bz[T]")
Bzfromxz = r.TGraph2D()
Bzfromxz.SetTitle("Bz map (xz);x[cm];z[cm];Bz[T]")
Bzfromyz = r.TGraph2D()
Bzfromyz.SetTitle("Bz map (yz);y[cm];z[cm];Bz[T]")


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
 Bxfromxz.SetPoint(n, x, z, Bx)
 Bxfromyz.SetPoint(n, y, z, Bx)
 Byfromxy.SetPoint(n, x, y, By)
 Byfromxz.SetPoint(n, x, z, By)
 Byfromyz.SetPoint(n, y, z, By)
 Bzfromxy.SetPoint(n, x, y, Bz)
 Bzfromxz.SetPoint(n, x, z, Bz)
 Bzfromyz.SetPoint(n, y, z, Bz)

#drawing plots
cx = r.TCanvas()
cx.Divide(1,2)
cx.cd(1)
Bxfromxy.GetXaxis().SetTitle('x[cm]')
Bxfromxy.GetYaxis().SetTitle('y[cm]')
Bxfromxy.Draw('COLZ')
cx.Update()
cx.cd(2)
Bxfromxz.GetXaxis().SetTitle('x[cm]')
Bxfromxz.GetYaxis().SetTitle('z[cm]')
Bxfromxz.Draw('COLZ')
cx.Update()

cy = r.TCanvas()
cy.Divide(1,2)
cy.cd(1)
Byfromxy.GetXaxis().SetTitle('x[cm]')
Byfromxy.GetYaxis().SetTitle('y[cm]')
Byfromxy.Draw('COLZ')
cy.Update()
cy.cd(2)
Byfromxz.GetXaxis().SetTitle('x[cm]')
Byfromxz.GetYaxis().SetTitle('z[cm]')
Byfromxz.Draw('COLZ')
cy.Update()

cz = r.TCanvas()
cz.Divide(1,2)
cz.cd(1)
Bzfromxy.GetXaxis().SetTitle('x[cm]')
Bzfromxy.GetYaxis().SetTitle('y[cm]')
Bzfromxy.Draw('COLZ')
cz.Update()
cz.cd(2)
Bzfromxz.GetXaxis().SetTitle('x[cm]')
Bzfromxz.GetYaxis().SetTitle('z[cm]')
Bzfromxz.Draw('COLZ')
cz.Update()

cyz = r.TCanvas()
cyz.Divide(1,2)
cyz.cd(1)
Bxfromyz.GetXaxis().SetTitle('y[cm]')
Bxfromyz.GetYaxis().SetTitle('z[cm]')
Bxfromyz.Draw('COLZ')
cyz.Update()
cyz.cd(2)
Byfromyz.GetXaxis().SetTitle('y[cm]')
Byfromyz.GetYaxis().SetTitle('z[cm]')
Byfromyz.Draw('COLZ')
cyz.Update()

cyz2 = r.TCanvas()
Bzfromyz.GetXaxis().SetTitle('y[cm]')
Bzfromyz.GetYaxis().SetTitle('z[cm]')
Bzfromyz.Draw('COLZ')
cyz2.Update()
