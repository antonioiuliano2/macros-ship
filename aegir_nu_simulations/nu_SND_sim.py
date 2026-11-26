#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 CERN for the benefit of the SHiP Collaboration
# SPDX-License-Identifier: LGPL-3.0-or-later
"""A first look at an aegir file (booklet ch. 10): hits per detector, and
where along the beam the straw hits are.

    cp $CODE/analysis/first_sim.py .
    python first_sim.py muon_gun.root

Copy it to work/ and extend it in the exercises.
"""

import collections
import os
import sys

import ROOT
import SND_SubDetInfo as SND

def finish(status=0):
    """Flush and exit without PyROOT teardown (see _alive)."""
    sys.stdout.flush()
    sys.stderr.flush()
    os._exit(status)

#sys.path.append(os.path.join(os.environ.get("CODE", ""), "analysis"))
#import shipfiles as sf  # noqa: E402

path = sys.argv[1] if len(sys.argv) > 1 else "muon_gun.root"

# An aegir file holds one RNTuple, "events", with one entry per event. A view
# reads one field; calling it with an entry number gives that entry's value.
reader = ROOT.RNTupleReader.Open("events", path)
hits = reader.GetView["std::vector<SHiP::SimHit>"]("sim_hits")

per_detector = collections.Counter()
h_z = ROOT.TH1D("SND_z", "SND hits;z [m];hits", 6000, 26, 32)
h_xy = ROOT.TH2D("SND_xy", "SND hits;x [m];y [m]; hits", 100, -0.5, 0.5,100,-0.5,0.5)
h_zy = ROOT.TH2D("SND_zy", "SND hits;z [m];y [m]; hits", 100, 26, 32,100,-0.5,0.5)
h_ID = ROOT.TH1D("SND_detID","SND Detector ID (here being main subdetector)",10,0,10)

h_HGC_ids = ROOT.TH1D("HGC_detIDs","HGC Plane ID Hits",125,0,125)
h_MGC_ids = ROOT.TH1D("MGC detIDs","MGC Plane ID Hits",10,0,10)
h_MTC_ids = ROOT.TH1D("MTC detIDs","MTC Plane ID Hits",43,0,43)

h_subids = {"HGC":h_HGC_ids,"MGC":h_MGC_ids,"MTC":h_MTC_ids}

ROOT.gStyle.SetOptStat(111111)

SND_ZSTART = 26425.000 #units are in mm
SND_ZMTCSTART = 28000.000
SND_ZEND = 31475.000
for i in range(reader.GetNEntries()):
    nhits_muMTC = 0
    for h in hits(i):  # the hits of event i
        if (h.position[2] > SND_ZSTART and h.position[2] < SND_ZEND):
         #if(h.pdg_code == 13 and h.track_id == 1):
         #print(i,h.position[0],h.position[1],h.position[2])
         h_z.Fill(h.position[2] / 1000.0)
         h_xy.Fill(h.position[0] / 1000.0, h.position[1] / 1000.0)
         h_zy.Fill(h.position[2] / 1000.0, h.position[1] / 1000.0)
         h_ID.Fill(h.detector_id)
         #now, get subdetector_id, according to detector_id
         subdetectorID, detname = SND.get_subdetID(h.detector_id,h.position[2])
         per_detector[detname] += 1
         h_subids[detname].Fill(subdetectorID)

print(f"{path}: {reader.GetNEntries()} events")
for name, n in per_detector.most_common():
    print(f"  {name!s:12} {n:7d} hits")

canvas = ROOT.TCanvas("c", "", 800, 600)
h_z.Draw("hist")
#canvas.SaveAs("snd_z.pdf")

c2D = ROOT.TCanvas("c2D", "", 800, 800)
h_xy.Draw("COLZ")
#c2D.SaveAs("SND_XY.pdf")

c2D_zy = ROOT.TCanvas("c2D_zy", "", 800, 800)
h_zy.Draw("COLZ")

cID = ROOT.TCanvas("cID","",800,800)
h_ID.Draw()
#c2D_zy.SaveAs("SND_ZY.pdf")
#sf.finish()  # leave without PyROOT's teardown, which can crash

cDetIds = ROOT.TCanvas("cDetIDs")
cDetIds.Divide(2,2)
cDetIds.cd(1)
h_HGC_ids.Draw()
cDetIds.cd(2)
h_MGC_ids.Draw()
cDetIds.cd(3)
h_MTC_ids.Draw()