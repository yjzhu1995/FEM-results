#!/usr/bin/env -S bash -e
# GMT modern mode bash template
# Date:    2020-07-24T12:48:37
# User:    yzhu
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
R1=0/15/0/0.6
R2=-5/10/-0.1/0.9
R3=-1/100/0/20
pi=3.1415926
blue=black
red=black,4_3:0
ss=2D_vis_eq_W50D_MaxwellM1e19_dep0_total
T1=100
T2=500
adis=/e/strike-slip/half2Dvis_MaxwellM1e19_sampledbyW20D/adis/
avel=/e/strike-slip/half2Dvis_MaxwellM1e19_sampledbyW20D/avel/
gmt begin dis_tp_100_500 png
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 2x3 -Fs9c/4.8c -M-0.1c/-0.2c -A+o0.2c/0.2c

		gmt subplot set 0,0 -A'(A) t = 0.1T'
			gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1+l"Displacement (U/U@-o@-)" -BWsen	
			awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t30_T300 | gmt plot -JX? -R$R1 -W1,	
#		awk '{print $1*1110/15,atan2($1*1110,15)/3.1415/10}' ./adis/nadis_${ss}_t30_T300 | gmt plot -JX? -R$R1 -W1,gray		
			awk '{print $1*1110/15,$2/10}' ${adis}nadis_${ss}_t30_T1${T1}_T2${T2}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,red
			awk '{print $1*1110/15,$2/50}' ${adis}nadis_${ss}_t30_T1${T2}_T2${T1}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,green

		gmt subplot set 0,1 -A'(b) t = 0.5T'
			# gmt basemap -JX? -R$R2 -Bxa5 -Bya0.2+l"Fault-normal Velocity (V/V@-o@-)" -BWSen
			gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1+l" " -Bwsen
			awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t150_T300 | gmt plot -JX? -R$R1 -W1,	
#			awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t150_T200 | gmt plot -JX? -R$R1 -W1,gray
			# awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t150_T200 | gmt plot -JX? -R$R1 -W2,4:3_0		
			# awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t150_T400 | gmt plot -JX? -R$R1 -W2,8:6_0
			# awk '{print $1*1110/15,atan2($1*1110,15)/3.1415/2}' ./adis/nadis_${ss}_t30_T300 | gmt plot -JX? -R$R1 -W1,gray
#			awk '{print $1*1110/15,$2/30}' ${adis}nadis_${ss}_t150_T1${T1}_T2${T2}_ | gmt plot -JX? -R$R1 -W1,red
			awk '{print $1*1110/15,$2/50}' ${adis}nadis_${ss}_t150_T1${T2}_T2${T1}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,green
			awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t150_T500 | gmt plot -JX? -R$R1 -W1,gray

		gmt subplot set 0,2 -A'(c) t = T'
			# gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1+l"Fault-parallel Displacement (U/U@-o@-)" -BWSen
			gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1+l" " -BwsEn
			awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t300_T300 | gmt plot -JX? -R$R1 -W1,
			awk '{print $1*1110/15,$2/50}' ${adis}nadis_${ss}_t300_T1${T2}_T2${T1}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,green
			awk '{print $1*1110/15,$2/20}' ${adis}nadis_${ss}_t300_T500 | gmt plot -JX? -R$R1 -W1,gray


		gmt subplot set 1,0 -A'(d) t = 0.1T'
			gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1+l"Velocity (V/V@-o@-)" -BWSen
			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t30_T300 | gmt plot  -W1,
			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t30_T1${T1}_T2${T2}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,red
			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t30_T1${T2}_T2${T1}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,green
		gmt subplot set 1,1 -A'(e) t = 0.5T'
			# gmt basemap -JX? -R$R2 -Bxa5 -Bya0.2+l"Fault-normal Velocity (V/V@-o@-)" -BWSen
			gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1 -BwSen
			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t150_T300 | gmt plot  -W1,
#			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t150_T1${T1}_T2${T2}_ | gmt plot -JX? -R$R1 -W1,red
			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t150_T1${T2}_T2${T1}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,green
		gmt subplot set 1,2 -A'(f) t = T'
			# gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1+l"Fault-parallel Displacement (U/U@-o@-)" -BWSen
			gmt basemap -JX? -R$R1 -Bxa5+l"Distance (x/D)" -Bya0.1 -BwSEn
			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t300_T300 | gmt plot  -W1,
			awk '{print $1*1110/15,$2*10}' ${avel}navel_${ss}_t300_T1${T2}_T2${T1}_irr2eq_tp | gmt plot -JX? -R$R1 -W1,green
	gmt subplot end 
gmt end 

