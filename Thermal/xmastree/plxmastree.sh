#!/usr/bin/env -S bash -e
# GMT modern mode bash template
# Date:    2020-07-24T12:48:37
# User:    yzhu
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
R=-1000/2000/0/75
J=X10c/-6c
gmt begin figure8-4 png
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 1x1 -Fs12c/8c -M0.2c/0.2c 
	#gmt subplot begin 2x2 -Fs8c/6c -M0.6c/0.2c -A
		gmt subplot set 0 
			# gmt basemap -JX? -R$R2 -Bxa5+l"Distance (x/D)" -Bya0.2+l"Fault-normal Displacement (U/U@-o@-)" -BWSen
			gmt basemap -J$J -R$R -Bxa1000f200+l"Differential Stress (MPa)" -Bya10f5+l"Depth (km)" -BWseN
			gmt plot middle -J$J -R$R  -W0.5,. 

			# awk '{print -$2,$1}' b.1 | gmt plot -J$J -R$R  -W1,red
			# awk '{print $2,$1}' b.2 | gmt plot -J$J -R$R  -W1,red
			# awk '{print -$2,$1}' d.1 | gmt plot -J$J -R$R  -W1,green
			# awk '{print $2,$1}' d.1 | gmt plot -J$J -R$R  -W1,green
			awk '{print -$2,$1}' bd1 | gmt plot -J$J -R$R  -W1,
			awk '{print $2,$1}' bd2 | gmt plot -J$J -R$R  -W1,
			# awk '{print -$2*2/1000000,$1}' dislocation.dat | gmt plot -J$J -R$R  -W1,blue			
			# awk '{print $2*2/1000000,$1}' dislocation.dat | gmt plot -J$J -R$R  -W1,blue		
			# awk '{print $2*2/1000000,$1}' dislocation_new.dat | gmt plot -J$J -R$R  -W1,gray	
			# awk '{print $2/1000000,$1}' dorn.dat | gmt plot -J$J -R$R  -W1,yellow
	gmt subplot end 
gmt end 

