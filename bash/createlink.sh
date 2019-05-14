#!bin/bash
#a small script to make the links to the data folder for a lazy Antonio
#does not work from 1 to 9 due to the need of adding a 0 in the P$i folder

for i in $(seq -w 10 14);
do
 mkdir p0$i
 ln -s /mnt/data/FOOT/2019_GSI/mic5/GSI2/P0$i/tracks.raw.root /home/scanner/foot/2019_GSI/GSI2/b000002/p0$i/2.$i.0.0.raw.root
done
 
