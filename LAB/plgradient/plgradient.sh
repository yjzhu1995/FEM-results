#! /bin/bash 
# GMT modern mode bash template
# Date:    2020-07-24T12:48:37
# User:    yzhu
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name

pro=/home1/yzhu/xin_code/2term/Haida_Gwaii_3layer_E20km_H420M1e19K8e17rig1.6_M1e20/disp/aveloc.dat
# pro2=/home1/yzhu/xin_code/2term/Haida_Gwaii_3layer_E20km_H420M1e19K8e17rig1.28_M1e20/disp/aveloc.dat
# pro3=/home1/yzhu/xin_code/2term/Haida_Gwaii_3layer_E20km_H200M1e19K4e17rig0.96_M1e20/disp/aveloc.dat
pro2=/home1/yzhu/xin_code/2term/Haida_Gwaii_3layer_E20km_H200M1e19K4e17rig0.8_M1e20/disp/aveloc.dat
pro3=/home1/yzhu/xin_code/2term/Haida_Gwaii_3layer_E20km_H100M1e19K2e17rig0.25_M1e20/disp/aveloc.dat
color=(green lightblue red orange purple )
gmt begin  gradient_profile png A+m0.5c
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 1x1 -Fs5i/5i -M0.3c/0.3c -A+o0.2c/0.2c
    R=0/500/0/15
        gmt subplot set 0

            awk '{print $2,$3,$4}'  $pro | gmt plot  -R$R -JX? -Sc0.2c -Gblack -Ey -W1,black+s  -Bxa100+l"Diostance along profile (km)" -Bya5+l"Average Dis (16~20)(mm)" -BWSen
            awk '{print $2,$5}'  $pro | gmt plot -Sc0.2c -Gred -W1,red+s  
            awk '{print $2,$5}'  $pro3 | gmt plot  -Sc0.2c -Gdarkblue -W1,darkblue+s 
            awk '{print $2,$5}'  $pro2 | gmt plot  -Sc0.2c -Gorange -W1,orange+s 

		gmt subplot end 
gmt end 
