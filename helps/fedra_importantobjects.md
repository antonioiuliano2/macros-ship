# FREQUENTLY USED OBJECTS IN FEDRA AND CONCERNING METHODS
Detailed infromation can be found on website http://operaweb.lngs.infn.it/operawiki/index.php/FEDRA or directly in fedra code
# Tracking and Vertexing
## Segment: EdbSegP
Location in fedra: /fedra/src/libEbase/EdbSegP.h 
Class where base tracks are stored. Segment ID can be accessed with ID() 
* To access kinematics: X(),Y(),Z(),TX(),TY() (SX(),SY(),SZ(),STX() and STY() for sigmas) 
* Spherical angle coordinates: Phi() and Theta() (actually Theta() gives TanTheta(), due to historical bug) 
* Emulsion information: Plate(), DZ(), DZem(), W() (the latter weigth w() returns the number of clusters) 
* Fit information (from linking): Chi2(), Prob() 
* MC information (if available): MCEvt(), MCTrack(). Set with   SetMC( int mEvt, int mTrack )
## Track: EdbTrackP
Location in fedra: /fedra/src/libEdr/EdbPattern.h
Class for volume track, inherits from EdbSegP. ID() returns ID of Track.
In general, segments are the data after linking, fitted segmented are the reconstructed positions of the track in the plates after track fit
* To access segments, GetSegmentFirst(), GetSegmentLast() and GetSegment(int i). N() returns number of segments associated to this track 
* Similar methods are used to access fitted segment. Usually are written the same, with an F after Segment (example. GetSegmentFFirst())
* Npl() is the number of plates expected. Different from N(), which the number of plates where the segment has been found.

Important note: linked_tracks.root files do NOT contain EdbTrackP objects, but EdbSegP objects, and track information are stored in separated branches. Probably for memory saving.

## Vertex: EdbVertex
Location in fedra: /fedra/src/libEdr/EdbVertex.h
Class for vertex reconstruction. 
* Vertex Position: X(), Y(), Z() . Also present VX(), VY() and VZ(), which seem to return the same. The latter are the ones suggested in the Fedra Wiki to obtain vertex position;
* To obtain information of vertex fit: chi2(), prob(), ndf(); 
* To access associated tracks: N() gives the number of tracks. GetTrack(int i) gives the i-th track; 
* Impact(i) returns distance of i-th track from the vertex; 
* MaxAperture() returns the maximum angular distance between two tracks associated to this vertex (as in tan(theta)); 
* GetVTa(i) returns VTA object for i-th track; 
* Flag() returns a flag describing vertex: 0 only begin tracks, 1 both begin and end present, 2 only end track present, -10 all tracks associated to other, higher-ranked, vertices

EdbVTA is the object returned by GetVTA(i). Declared also in EdbVertex.h. It describes Vertex-track association:
* Zpos() returns 0 for end track, 1 for begin track; 
* Imp() returns impact parameter; 
* Dist() returns distance from vertex to nearest track point; 
* GetTrack() and GetVertex() return track and vertex objects from the association; 

## EdbVertexComb
Location in fedra: /fedra/src/libEdr/EdbVertexComb.h
Prepare track combinations to study 'anomalous OPERA event in Bari'. Still check for info