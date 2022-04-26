#!/usr/bin/env -S bash -e

export GMT_SESSION_NAME=123	# Set a unique session name
R=-134/-120/47/55
J=M3i
arrow=0.3c+ea+h0+p1.5p,black+n30
gmt begin GPS_loc png
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12,
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 1x1 -Fs10c/10c -M0.3c/0.5c 
		gmt subplot set 0,0
			gmt coast -JM? -R$R  -Dh -A100 -G244/243/239 -Wthinnest,244/243/239 -S167/194/223   -Bxa5 -Bya5 -BWSen
            awk '{print $3,$2}' HaidaGwaii.loc | gmt plot -Ss0.2c -Gblue 
			awk '{print $3,$2,$6,$5,$8,$7,0}' disp_VI.txt | gmt velo -Se0.1c/1/0 -A$arrow -Gblack -W0.5,black
			echo "-133 48 10 0 0 0 0" | gmt velo -Se0.1c/1/0 -A$arrow -Gblack -W0.5,black
			echo "-132.5 48.5 10 mm" | gmt text 
            #awk '{ print $3,$2,$1}'  ssite.txt | gmt text -D-0.5/-0.5+v1p,red
	gmt subplot end
gmt end 