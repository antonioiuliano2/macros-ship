#!/bin/bash

filepath="$HOME/cernbox/Synched/Charmdata/BDT_vertices_Valerio/afterBDT_plots"

hadd $filepath/full_vz.root $filepath/full_vz_CH*.root

hadd $filepath/with_bkg_full_vz.root $filepath/with_bkg_full_vz_CH*.root