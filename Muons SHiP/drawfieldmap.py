import ROOT as r
import numpy

#inputfile = r.TFile.Open('$FAIRSHIP/files/nuTauDetField.root')
inputfile = r.TFile.Open('$HOME/Lavoro/SHIPNuTauBField/nuTauDetField_FairShipRef.root')
tree = inputfile.Get('Data')
#defining graphs
outputfile = r.TFile('field_maps.root','RECREATE')
    
range_tree = inputfile.Get('Range')   
range_tree.GetEntry(0)

#getting range information for the 3 coordinates
xMin = range_tree.xMin
xMax = range_tree.xMax
dx = range_tree.dx
nbinsx = int((xMax - xMin)/(dx))
yMin = range_tree.yMin
yMax = range_tree.yMax
dy = range_tree.dy
nbinsy = int((yMax - yMin)/(dy))
zMin = range_tree.zMin
zMax = range_tree.zMax
dz = range_tree.dz
nbinsz = int((zMax - zMin)/(dz))
#I need to save the maps for the three B coordinates
hBx = r.TProfile3D("hBx","Bx in space", nbinsx+1,xMin, xMax, nbinsy+1,yMin, yMax, nbinsz+1, zMin, zMax)
hBy = r.TProfile3D("hBy","By in space", nbinsx+1,xMin, xMax, nbinsy+1,yMin, yMax, nbinsz+1, zMin, zMax)
hBz = r.TProfile3D("hBz","Bz in space", nbinsx+1,xMin, xMax, nbinsy+1,yMin, yMax, nbinsz+1, zMin, zMax)
print('Numero bin nei 3 assi:',nbinsx, nbinsy, nbinsz)

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

Bfromxz = r.TGraph2D()
Bfromxz.SetTitle("B map(xz);z[cm];x[cm];B[T]")
Bfromyz = r.TGraph2D()
Bfromyz.SetTitle("B map(yz);z[cm];y[cm];B[T]")

for i in range (tree.GetEntries()):
 tree.GetEntry(i)
 n = Byfromxy.GetN()
 #getting coordinates
 x = tree.x #width
 y = tree.y #height
 z = tree.z #beam direction
 #getting field values
 Bx = tree.Bx
 By = tree.By
 Bz = tree.Bz

 B = pow(Bx*Bx + By*By + Bz*Bz,0.5)

 hBx.Fill(x,y,z,Bx)
 hBy.Fill(x,y,z,By)
 hBz.Fill(x,y,z,Bz)
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

 Bfromyz.SetPoint(n,z,y,B)
 Bfromxz.SetPoint(n,z,x,B)

#drawing all the plots and saving them to file
cB = r.TCanvas()
cB.Divide(1,2)
cB.cd(1)
Bfromxz.GetXaxis().SetTitle('z[cm]')
Bfromxz.GetYaxis().SetTitle('x[cm]')
Bfromxz.Draw('COLZ')
Bfromxz.Write()
cB.Update()
cB.cd(2)
Bfromyz.GetXaxis().SetTitle('z[cm]')
Bfromyz.GetYaxis().SetTitle('y[cm]')
Bfromyz.Draw('COLZ')
Bfromyz.Write()
cB.Update()

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


print hBx.GetBinContent(2,1,1)
import numpy as np

outFile = open('mymap.txt', 'w')
outFile.write('#        x(m)          y(m)            z(m)           Bx(T)           By(T)           Bz(T)\n')
    #Expected BIN ORDERING: z -> y -> x
for ix,x in enumerate(np.linspace(xMin,xMax + dx, nbinsx)):
        for iy,y in enumerate(np.linspace(yMin, yMax + dy, nbinsy)):
            for iz,z in enumerate(np.linspace(zMin, zMax + dz, nbinsz)):
                #Getting field value                
                Bx = hBx.GetBinContent(ix+1,iy+1,iz+1) #bin ids start from 1
                By = hBy.GetBinContent(ix+1,iy+1,iz+1)
                Bz = hBz.GetBinContent(ix+1,iy+1,iz+1)
                #Printing the output
                outFile.write('{0:<15.7e} {1:>15.7e} {2:>15.7e} {3:>15.7e} {4:>15.7e} {5:>15.7e}\n'.format(x,y,z,Bx,By,Bz))
                    

outputfile.Close()
