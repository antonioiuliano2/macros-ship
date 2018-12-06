import ROOT as r
f = r.TFile('./studionoise/daanalizzare/particles_number_array.root')
t_film = f.Get('tree')

f_detectors = r.TFile('./studionoise/temporanei/particles_charm_detectors_centralbeam.root')
#f_detectors = r.TFile('./studionoise/temporanei/particles_number_detectors.root')
t_detectors = f_detectors.Get('tree')

hdistanza = r.TH1D("hdistanza","Distanza fra i punti nell'ultimo film di emulsione", 100, 0, 2);
hdistanza1 = r.TH1D("hdistanza1","Distanza fra i punti nel primo piano dello spettrometro", 100, 0, 2);
hdistanza2 = r.TH1D("hdistanza2","Distanza fra i punti nel secondo RPC", 80, 0, 8);

hdistanzaall = r.TH1D("hdistanzaall","Distanza fra i punti nell'ultimo film di emulsione", 100, 0, 2);
hdistanza1all = r.TH1D("hdistanza1all","Distanza fra i punti nel primo piano dello spettrometro", 100, 0, 2);
hdistanza2all = r.TH1D("hdistanza2all","Distanza fra i punti nel secondo RPC", 80, 0, 8);

c = r.TCanvas()
t_film.Draw("xPoint:yPoint>>h(50,-5,5,50,-5,5)","nfilm == 45","COLZ")
xaxis = r.h.GetXaxis()
xaxis.SetTitle("cm")
yaxis = r.h.GetYaxis()
yaxis.SetTitle("cm")
c.Print("./studionoise/temporanei/immaginidistanze/posizioni_ultimofilm.png")
c.Print("./studionoise/temporanei/immaginidistanze/posizioni_ultimofilm.root")

c_zoom = r.TCanvas()
t_film.Draw("xPoint:yPoint>>h1(50,-5,-1,50,2,5)","nfilm == 45","COLZ")
xaxis = r.h1.GetXaxis()
xaxis.SetTitle("cm")
yaxis = r.h1.GetYaxis()
yaxis.SetTitle("cm")
c_zoom.Print("./studionoise/temporanei/immaginidistanze/posizioni_ultimofilm_zoom.png")
c_zoom.Print("./studionoise/temporanei/immaginidistanze/posizioni_ultimofilm_zoom.root")

c0 = r.TCanvas()
t_detectors.Draw("xPoint:yPoint>>h2(80,-8,8,80,-8,8)","nplane == 1", "COLZ")
xaxis = r.h2.GetXaxis()
xaxis.SetTitle("cm")
yaxis = r.h2.GetYaxis()
yaxis.SetTitle("cm")
c0.Print("./studionoise/temporanei/immaginidistanze/posizioni_primopiano_spettrometro.png")
c0.Print("./studionoise/temporanei/immaginidistanze/posizioni_primopiano_spettrometro.root")

c0_zoom = r.TCanvas()
t_detectors.Draw("xPoint:yPoint>>h3(50,-5,-1,50,2,5)","nplane == 1","COLZ")
xaxis = r.h3.GetXaxis()
xaxis.SetTitle("cm")
yaxis = r.h3.GetYaxis()
yaxis.SetTitle("cm")
c0_zoom.Print("./studionoise/temporanei/immaginidistanze/posizioni_spettrometro_zoom.png")
c0_zoom.Print("./studionoise/temporanei/immaginidistanze/posizioni_spettrometro_zoom.root")

c0_1 = r.TCanvas()
t_detectors.Draw("xPoint:yPoint>>h4(100,-50,50,100,-50,50)","nplane == 6", "COLZ")
xaxis = r.h4.GetXaxis()
xaxis.SetTitle("cm")
yaxis = r.h4.GetYaxis()
yaxis.SetTitle("cm")
c0_1.Print("./studionoise/temporanei/immaginidistanze/posizioni_secondorpc.png")
c0_1.Print("./studionoise/temporanei/immaginidistanze/posizioni_secondorpc.root")

#calcolo le distanze fra i punti nell'ultimo film
#for i in range(t_film.GetEntries()):
for i in range(1):
  if (i % 10 == 0):
      print i
  t_film.GetEntry(i)
  xPoint1 = t_film.xPoint
  yPoint1 = t_film.yPoint
  zPoint1 = t_film.zPoint
  npoints = t_film.npoints
  nfilm1 = t_film.nfilm
  for count1 in range(npoints):
   mindistanza = 1000. #parto da un valore minimo molto alto
   if (nfilm1[count1] == 47):
    for count2 in range(npoints):
          if (nfilm1[count2] == 47):
            #print i,j
            distanza = pow(pow(xPoint1[count1] - xPoint1[count2],2) + pow(yPoint1[count1] - yPoint1[count2],2),0.5)
            if ((distanza < mindistanza) and (count1 != count2)):
              mindistanza = distanza
            if (count1 != count2):
              hdistanzaall.Fill(distanza)
    hdistanza.Fill(mindistanza)     

#rivelatori a valle
for i in range (t_detectors.GetEntries()):
#for i in range (100):  
  if (i % 10 == 0):
      print i
 
  t_detectors.GetEntry(i)
  npoints = t_detectors.npoints 
  xPoint1 = t_detectors.xPoint
  yPoint1 = t_detectors.yPoint
  zPoint1 = t_detectors.zPoint
  nplane1 = t_detectors.nplane
  for count1 in range(npoints):
   mindistanza = 1000. #la minima distanza per ogni punto da tutti gli altri
   if (nplane1[count1] == 1):
#    for j in range(t_detectors.GetEntries()):
    for count2 in range(npoints):
          #xPoint2 = t_detectors.xPoint
          #yPoint2 = t_detectors.yPoint
          #zPoint2 = t_detectors.zPoint
          #nplane2 = t_detectors.nplane
          if (nplane1[count2] == 1):
            #print i,j
            distanza = pow(pow(xPoint1[count1] - xPoint1[count2],2) + pow(yPoint1[count1] - yPoint1[count2],2),0.5)
            if ((distanza < mindistanza) and (count1 != count2)):
              mindistanza = distanza
            if (count1 != count2):
              hdistanza1all.Fill(distanza)
   hdistanza1.Fill(mindistanza) #dopo aver calcolato la minima distanza con tutti gli altri punti, riempo l'istogramma
   if (nplane1[count1] == 6):
    for count2 in range(npoints):
          #xPoint2 = t_detectors.xPoint
          #yPoint2 = t_detectors.yPoint
          #zPoint2 = t_detectors.zPoint
          #nplane2 = t_detectors.nplane
          if (nplane1[count2] == 6):
            #print i,j
            distanza = pow(pow(xPoint1[count1] - xPoint1[count2],2) + pow(yPoint1[count1] - yPoint1[count2],2),0.5)
            if ((distanza < mindistanza) and (count1 != count2)):
              mindistanza = distanza
            if (count1 != count2):
              hdistanza2all.Fill(distanza)
   hdistanza2.Fill(mindistanza) #non considero lo stesso punto
    
c1 = r.TCanvas()
hdistanza.Draw()
xaxis = hdistanza.GetXaxis()
xaxis.SetTitle("cm")
c1.Print("./studionoise/temporanei/immaginidistanze/minima_distanza.png")
c1.Print("./studionoise/temporanei/immaginidistanze/minima_distanza.root")

c1_all = r.TCanvas()
hdistanzaall.Draw()
xaxis = hdistanzaall.GetXaxis()
xaxis.SetTitle("cm")
c1_all.Print("./studionoise/temporanei/immaginidistanze/tutte_distanze.png")
c1_all.Print("./studionoise/temporanei/immaginidistanze/tutte_distanze.root")

c2 = r.TCanvas()
hdistanza1.Draw()
xaxis = hdistanza1.GetXaxis()
xaxis.SetTitle("cm")
c2.Print("./studionoise/temporanei/immaginidistanze/minima_distanza_spettrometro.png")
c2.Print("./studionoise/temporanei/immaginidistanze/minima_distanza_spettrometro.root")

c2_all = r.TCanvas()
hdistanza1all.Draw()
xaxis = hdistanza1all.GetXaxis()
xaxis.SetTitle("cm")
c2_all.Print("./studionoise/temporanei/immaginidistanze/tutte_distanze_spettrometro.png")
c2_all.Print("./studionoise/temporanei/immaginidistanze/tutte_distanze_spettrometro.root")

c3 = r.TCanvas()
hdistanza2.Draw()
xaxis = hdistanza2.GetXaxis()
xaxis.SetTitle("cm")
c3.Print("./studionoise/temporanei/immaginidistanze/minima_distanza_rpc.png")
c3.Print("./studionoise/temporanei/immaginidistanze/minima_distanza_rpc.root")

c3_all = r.TCanvas()
hdistanza2all.Draw()
xaxis = hdistanza2all.GetXaxis()
xaxis.SetTitle("cm")
c3_all.Print("./studionoise/temporanei/immaginidistanze/tutte_distanze_rpc.png")
c3_all.Print("./studionoise/temporanei/immaginidistanze/tutte_distanze_rpc.root")


