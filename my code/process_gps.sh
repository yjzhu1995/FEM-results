#! /bin/bash
R=99/108/32/38
J=M2i
pi=3.1415927
#gmt set MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE FANCY FORMAT_GEO_MAP ddd:mm:ssF
#gmt gmtset MAP_TICK_LENGTH -0.1
#gmt pscoast -R$R -J$J -Bagf -BWSen -Dh -A -Glightbrown -Wthinnest -P -Slightblue -K -Y5i> $ps
#make gps file follow this format: Lon Lat Ve Vn dVe dVn
location=Maduo
awk '{if ($1 >96 && $1 < 100 && $2 > 33 && $2 < 36)  print $1,$2,$3,$5,$4,$6}' wang20.gps > ${location}.gps
alon=97.8139
alat=34.7463
blon=98.4431
blat=34.5606
C=${blon}/${blat}
E=${alon}/${alat}

gmt project -C$C -E$E -G10 -Q > trace_${location}
awk '{print $1,$2,$3,$4,$5,$6}' ${location}.gps | gmt project -C$C -E$E -Q > vprofile_${location}1_temp
theta=$(echo "$alon $alat $blon $blat" | awk '{print atan2(($3-$1)*cos($2*'$pi'/360+$4*'$pi'/360),$4-$2)}')
echo "$theta"
awk '{if ($8 <200 && $8>-200 && $7<100 && $7>=-100) print $s}' vprofile_${location}1_temp >vprofile_${location}1

awk '{print $8,$4*cos('$theta')+$3*sin('$theta'),$6*cos('$theta')+$5*sin('$theta')}' vprofile_${location}1 > vprofile_${location}1_rotated

######
awk '{print $1,$2,$3}' profile.insar | gmt project -C$C -E$E -Q > vprofile_${location}1_insar_temp



# elon=98.4074
# elat=39.9767

# C2=${elon}/${elat}
# E2=${dlon}/${dlat}
# awk '{print $1,$2,$3,$4,$5,$6}' ${location}.gps | gmt project -C$C2 -E$E2 -Q > vprofile_${location}_temp3
# awk '{if ($8 <200 && $8>-200 && $7<200 && $7>=00) print $s}' vprofile_${location}_temp3 >vprofile_${location}4
# awk '{print $8,$4*cos('$alpha')+$3*sin('$alpha'),$6*cos('$alpha')+$5*sin('$alpha')}' vprofile_${location}4 > vprofile_${location}4_rotated

#rm vprofile_${location}_temp*

