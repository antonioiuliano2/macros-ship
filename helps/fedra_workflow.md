# FEDRA analysis workflow for SHiP-charm (last modification 7 december)

I will list here all the procedure starting from LASSO processing:

Connec to scanner@nusrv9.na.infn.it (password is 'neutrino.01')
Path for analysis: /ship/CHARM2018

## New folder
Create a folder like this /ship/CHARM2018/CH1R1/b000001
Copy thickness.C (optional)

## Add a plate
* mkdir p001

* Create a symbolic link to the data folder (actual path depends from scanned brick):
`ln -s /mnt/data/CHARM2018/.../P01/tracks.raw.root 1.0.0.0.raw.root`

* Check thickness with 'thickness.C'

## Create a ScanSet object
'EdbScanSet' is an object which contains all the identificatives of the plates where you want to work, along with affine transformations currently known between the plates. We need to create an EdbScanSet object before every operation of libScan applications (emlink, emalign)  

We can create an EdbScanSet with the command:  

`makescanset -set=1.0.0.0 -dzbase=175 -dz=-1300 -from-plate=29 -to_plate=1`  

It could be useful to save the default configuration in a scanset.txt file, customized with the plates to be analized.  

Warning: Using -reset to reset parameters after an ill-done alignment resets also shrinkage parameters! This has no consequences after linking, but we lose information about how was the shrinkage (we can always find it in the report though). I usually prefer to manually insert the identity transformation in the AFF folder. 

## Linking
Linking is made to obtain basetracks from microtracks. Basetracks are chosen according to a chi-square minimization, using coordinates, angles and cluster number as input (reference FEDRA 2006).  
Linking is done with the command:

`emlink -set=1.0.0.0 -new`

After the linking, check the report b000001.0.0.0.link.pdf and verify that the shrinkage plots show some peaks.  

## Alignment

It is one of the most delicate and important operations, because it allows to connect the different plates. Alignment is done with the command:  

`emalign -set=1.0.0.0 -new`

In order to avoid too long waits, it is convenient to optimize the alignment parameters on a small area, then we can move to the whole scanned area (usually 3 iterations). Aligning parameters are to be inserted in 'align.rootrc'.    

As with the linking, check the report b000001.0.0.0.align.pdf to verify the presence of peaks in xy residuals and angle. If alignment is bad, reset the affine transformations by adding a row in AFF/*aff.par with the identity transformation (1.,0.,0.,1.,0.,0.) to reset the alignment iterations.

