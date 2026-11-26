'''information about position of SND subdetectors, to be used in the analysis of the hits in the SND'''

def get_Si_X_number(hit_z):
    """among the 120 stations, which one does the hit belong to"""
    Si_Zmin = 26488.7
    Si_pitch = 11.3 #actually needed the distance between two consecutive Si stations
    #there are 3 timing station, 18 mm each
    timing_dz = 18
    if hit_z > 27407:
        ztimingoffset = timing_dz * 2
    elif hit_z > 26937:
        ztimingoffset = timing_dz
    else:
        ztimingoffset = 0
    return int((hit_z - Si_Zmin - ztimingoffset) / Si_pitch)

def get_SiPad_number(hit_z):
    """among the 10 SiPads, which one does the hit belong to"""
    SiPad_Zmin = 27923.5
    SiPad_pitch = 20.0 #20.0 mm
    return int((hit_z - SiPad_Zmin) / SiPad_pitch)

def get_MTC_Station_number(hit_z):
    """among the 42 MTC stations, which one does the hit belong to"""
    MTC_Zmin = 28180.5 
    MTC_pitch = 80.0 #80.0 mm
    return int((hit_z - MTC_Zmin) / MTC_pitch)

def get_subdetID(hit_detectorID, hit_z):
    if (hit_detectorID == 0):
        return get_MTC_Station_number(hit_z), "MTC"
    elif (hit_detectorID == 1):
        return get_SiPad_number(hit_z), "MGC"
    elif (hit_detectorID == 4):
        return get_Si_X_number(hit_z), "HGC"
    else:
        print("No Valid Detector Number {}! Please check".format(hit_detectorID))
        return -1, "Unknown"