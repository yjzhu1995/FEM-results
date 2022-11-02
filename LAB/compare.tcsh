#!/bin/tcsh
#

# rm .gmtdefaults4 .gmt_bb_info .gmtcommands4

set scale = 2 
gmt gmtset MAP_FRAME_TYPE FANCY

#Bounding box
set w_lon = 231
set e_lon = 240
set s_lat = 48 
set n_lat = 52

#map region box
set map_box = -R$w_lon/$e_lon/$s_lat/$n_lat

#center of projection
set pro_lon = `echo $w_lon $e_lon | awk '{printf"%.2f",($1+$2)/2.0}'`
set pro_lat = `echo $n_lat $s_lat | awk '{printf"%.2f",($1+$2)/2.0}'`

#output file
set filename = modelcompare.ps

##PROJECTIONS (USE ONE COMMENT OUT REST)
########################################

#Cassini
#set map_pro = -Jc$pro_lon/$pro_lat/$scale

#Mercator Projection
set map_pro = -JM7i

#Lambert Projection
#set map_pro = -Ja$pro_lon/$pro_lat/1:50000000

##SETUP 
########################################

# shift origins
set xorigin = -Xf1.5c
set yorigin = -Yf4i

#create colormap for grid
gmt makecpt -Cglobe -T-10000/10000/100 -Z > etopo5_topo.cpt

#simply creates a frame and starts the file
gmt psbasemap  $map_pro $map_box $xorigin $yorigin -B0 -P -K  > $filename

##MAP MAKING -- PLAY AWAY!
########################################

##start with topography
#gmt grdimage /mnt/readynas/yjiang/gpslab/dems/etopo1/etopo1_bed.grd -Cetopo5_topo.cpt $map_pro $map_box -X -Y -P -O -K >> $filename
gmt grdimage ../etopo1_bed.grd -Cetopo5_topo.cpt $map_pro $map_box -X -Y -P -O -K >> $filename
##add lakes, rivers etc ... (use one of them)
# pscoast $map_pro $map_box -Df -Ia -N1 -C0/0/255 -O -K >> $filename
gmt pscoast $map_pro $map_box -X -Y -Df -A10 -Ia -N1 -O -K >> $filename

##add a basegrid - keep this uncommented at the end, to make sure PS-trailer is added.
gmt psbasemap  $map_pro $map_box -X -Y -B3/3 -Lfx1i/1i/0/150 -P -O -K >> $filename

##add tremor/earthquakes for 2011
#awk '{print $7,$6}' neic_2011_earthquakes.txt | psxy $map_pro $map_box -Sc0.2c -G0 -O -K >> $filename
# awk '{print $3,$2}' ETS_tremor_catalog.txt_nv35 | psxy $map_pro $map_box -Sc0.2c -G0 -O -K >> $filename

##how about faultlines on top of this? The '-M' is important as it then understands the file as 'multi-segment' since it contains many faults
#psxy alaska_faults.gmt $map_pro $map_box -M -W10/255/0/0 -O -K >> $filename

gmt psxy $map_pro $map_box  ../bird_pb.gmt -W2p -K -O -G0/0/0 -Sf5/3p  >> $filename
#awk '{print $1,$2}' fault.1 | gmt psxy $map_pro $map_box -W2p,green -K -O >> $filename
#awk '{print $1,$2}' fault.2 | gmt psxy $map_pro $map_box -W2p,green -K -O >> $filename
#awk '{print $2,$3}' QCF_coseismic_2 | gmt psxy $map_pro $map_box -W2p,yellow -K -O >> $filename
#awk '{print $2,$3}' QCF_postseismic_2 | gmt psxy $map_pro $map_box -W2p,purple -K -O >> $filename
##or some velocities
#  awk '{print $1,$2,$17,$16,$19,$18,$20,$3}' vel.txt  | psvelo $map_pro $map_box -X -Y -L -Gred -A+a40+e+p2.5 -W.4,red -Se0.2/0.65/4 -O -K >> $filename
awk '{print $1,$2,$3,$4,$5,$6,$7,$8}' na.gps  | gmt psvelo $map_pro $map_box -X -Y -L -Gred -A+a40+e+p2.5 -W.4,red -Se0.2/0.65/4 -O -K >> $filename
awk '{print $1,$2,$3,$4}' model.gps  | gmt psvelo $map_pro $map_box -X -Y -L -Gblue -A+a40+e+p2.5 -W.4,blue -Se0.2/0.65/4 -O -K >> $filename
#awk '{print $1,$2,$3*1000,$4*1000}' fault.codis  | gmt psvelo $map_pro $map_box -X -Y -L -Ggreen -A+a40+e+p2.5 -W.4,green -Se0.2/0.65/4 -O -K >> $filename
#awk '{print $1,$2,$3*1000,-$4*1000}' trace.gps  | gmt psvelo $map_pro $map_box -X -Y -L -Ggreen -A+a40+e+p2.5 -W.4,green -Se0.2/0.65/4 -O -K >> $filename

# awk '{print $1,$2,$17,$16,$19,$18,$20,$3}' vel.lisa2  | psvelo $map_pro $map_box -X -Y -L -Gpurple -A+a40+e+p2.5 -W.4,purple -Se0.1/0.65/4 -O -K >> $filename
# awk '{print $1,$2,$5-$3,$6-$4,$7,$8,$9,$10}' diff_vel.dat  | psvelo $map_pro $map_box -X -Y -L -G0/255/255 -A+a40+e+p1.5 -W.4 -Sr0.2/0.65/3 -O -K >> $filename

## add velocity scale label
echo '-132  49  10  0  2  2   0   10+/-2 mm/a' | gmt psvelo $map_pro $map_box -X -Y -L -Gred -A+a40+e+p2.5 -W.4,red -Sr0.2/0.65/10 -O -K >> $filename
echo '-132  48.5 20  0  4  4   0   20+/-4 mm/a' | gmt psvelo $map_pro $map_box -X -Y -L -Gpurple -A+a40+e+p2.5 -W.4,purple -Sr0.1/0.65/10 -O -K >> $filename
# echo '-128.75 51.35 12,Helvetica,black 0 BR 5mm' | pstext $map_pro $map_box -F+f+a+j -O -K >> $filename

## PS TRAILER DO NOT COMMENT OUT!
########################################

#simply adds a frame around the map and finishes the file
# psbasemap  $map_pro $map_box -Y -B0 -P -O  >> $filename

#convert to pdf
 ps2pdf $filename

rm gmt.*
