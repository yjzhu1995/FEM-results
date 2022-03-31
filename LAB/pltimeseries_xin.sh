#! /bin/bash 
# GMT modern mode bash template
# Date:    2020-07-24T12:48:37
# User:    yzhu
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
ori=/home1/yzhu/xin_code/2term/benchmark4_2d_TH500km_L2000km/timeseries/
# ori_maxwell=/home1/yzhu/thrust2d/LAB/p_thrust15_W50km_E30km_1e18/
lab1=/home1/yzhu/xin_code/2term/thrust15_fw50km_E30km_TH500km_L2000km_lab1e16/timeseries/
lab2=/home1/yzhu/xin_code/2term/thrust15_fw50km_E30km_TH500km_L2000km_lab1e15/timeseries/
lab3=/home1/yzhu/xin_code/2term/thrust15_fw50km_E30km_TH500km_L2000km_lab1e14/timeseries/
fem_lab2=/home1/yzhu/thrust2d/LAB/p_thrust15_W50km_E30km_M1e19K1e18_lab10km_1e15/
fem_lab1=/home1/yzhu/thrust2d/LAB/p_thrust15_W50km_E30km_M1e19K1e18_lab10km_1e16/
# fem2d=/home1/yzhu/thrust2d/benchmark/Benchmark_pollitz_L1000km_nongrav/
# fem2d_grav=/home1/yzhu/thrust2d/benchmark/Benchmark_pollitz_L1000km_grav/

factor=111.1949266
ntime=3
times=(0.5 5 50 )
R1=0/5/-3/1
R2=0/5/-15/1
R3=0/5/-30/1
R4=0/5/-5/1
pi=3.1415926


gmt begin LAB_2d_timeseries_Me19K1e18_xin png
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 2x2 -Fs5i/2i -M0.3c/0.3c -A+o0.2c/4c
		gmt subplot set 0,0 -A'(a) Near-field (x=50km)'
                        gmt basemap -JX? -R$R1 -Bxa1 -Bya1+l"East displacement (cm)" -BWSen
                        awk '{print $1,$4*100}' ${ori}/near-field.ts | gmt plot -W1,gray+s           
                        awk '{print $1,$4*100}' ${lab1}/near-field.ts| gmt plot -W1,red+s      
                        awk '{print $1,$4*100}' ${fem_lab1}/near-field.ts| gmt plot -W1,red,4:3_0 
                        awk '{print $1,$4*100}' ${lab2}/near-field.ts | gmt plot -W1,green+s  
                        awk '{print $1,$4*100}' ${fem_lab2}/near-field.ts| gmt plot -W1,green,4:3_0 
                        awk '{print $1,$4*100}' ${lab3}/near-field.ts | gmt plot -W1,blue+s                          

		gmt subplot set 0,1 -A'(b) Far-field (x=500km)'
                        gmt basemap -JX? -R$R2 -Bxa1 -Bya5+l"East displacement (cm)" -BwSEn
                        awk '{print $1,$4*100}' ${ori}/far-field.ts | gmt plot -W1,gray+s           
                        awk '{print $1,$4*100}' ${lab1}/far-field.ts| gmt plot -W1,red+s      
                        awk '{print $1,$4*100}' ${fem_lab1}/far-field.ts| gmt plot -W1,red,4:3_0 
                        awk '{print $1,$4*100}' ${lab2}/far-field.ts | gmt plot -W1,green+s  
                        awk '{print $1,$4*100}' ${fem_lab2}/far-field.ts| gmt plot -W1,green,4:3_0 
                        awk '{print $1,$4*100}' ${lab3}/far-field.ts | gmt plot -W1,blue+s     
                        echo "2  -2.5  No LAB" | gmt text -F+f10p,gray+jML
                        # echo "3  -12.5   No LAB (1e18)" | gmt text -F+f10p,black+jML 
                        echo "3  -2.5   LAB (1e16)"| gmt text -F+f10p,red+jML 
                        echo "2  -5   LAB (1e15)"| gmt text -F+f10p,green+jML 
                        echo "3  -5   LAB (1e14)"| gmt text -F+f10p,blue+jML 
		gmt subplot set 1,0 -A'(c) Near-field (x=50km)'
                        gmt basemap -JX? -R$R3 -Bxa1 -Bya5+l"Vertical displacement (cm)" -BWSen
                        awk '{print $1,$2*100}' ${ori}/near-field.ts | gmt plot -W1,gray+s           
                        awk '{print $1,$2*100}' ${lab1}/near-field.ts| gmt plot -W1,red+s      
                        awk '{print $1,$2*100}' ${fem_lab1}/near-field.ts| gmt plot -W1,red,4:3_0 
                        awk '{print $1,$2*100}' ${lab2}/near-field.ts | gmt plot -W1,green+s  
                        awk '{print $1,$2*100}' ${fem_lab2}/near-field.ts| gmt plot -W1,green,4:3_0 
                        awk '{print $1,$2*100}' ${lab3}/near-field.ts | gmt plot -W1,blue+s     
		gmt subplot set 1,1 -A'(d) Far-field (x=500km)'
                        gmt basemap -JX? -R$R4 -Bxa1 -Bya1+l"Vertical displacement (cm)" -BwSEn
                        awk '{print $1,$2*100}' ${ori}/far-field.ts | gmt plot -W1,gray+s           
                        awk '{print $1,$2*100}' ${lab1}/far-field.ts| gmt plot -W1,red+s      
                        awk '{print $1,$2*100}' ${fem_lab1}/far-field.ts| gmt plot -W1,red,4:3_0 
                        awk '{print $1,$2*100}' ${lab2}/far-field.ts | gmt plot -W1,green+s  
                        awk '{print $1,$2*100}' ${fem_lab2}/far-field.ts| gmt plot -W1,green,4:3_0 
                        awk '{print $1,$2*100}' ${lab3}/far-field.ts | gmt plot -W1,blue+s   
		gmt subplot end 
gmt end 
