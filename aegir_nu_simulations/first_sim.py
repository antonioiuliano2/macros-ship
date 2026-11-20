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
h_z = ROOT.TH1D("SND_z", "SND hits;z [m];hits", 100, 26, 32)
h_xy = ROOT.TH2D("SND_xy", "SND hits;x [m];y [m]; hits", 100, -0.5, 0.5,100,-0.5,0.5)
h_zy = ROOT.TH2D("SND_zy", "SND hits;z [m];y [m]; hits", 100, 26, 32,100,-0.5,0.5)

SND_ZSTART = 26425.000
SND_ZMTCSTART = 28000.000
SND_ZEND = 31475.000
for i in range(reader.GetNEntries()):
    nhits_muMTC = 0
    for h in hits(i):  # the hits of event i
        if (h.position[2] > SND_ZMTCSTART and h.position[2] < SND_ZEND):
         if(h.pdg_code == 13 and h.track_id == 1):
          #print(i,h.position[0],h.position[1],h.position[2])
          h_z.Fill(h.position[2] / 1000.0)
          h_xy.Fill(h.position[0] / 1000.0, h.position[1] / 1000.0)
          h_zy.Fill(h.position[2] / 1000.0, h.position[1] / 1000.0)
          nhits_muMTC = nhits_muMTC + 1
    print("WE GOT ",nhits_muMTC, " MTC Hits")

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
#c2D_zy.SaveAs("SND_ZY.pdf")
#sf.finish()  # leave without PyROOT's teardown, which can crash
