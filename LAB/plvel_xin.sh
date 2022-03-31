#! /bin/bash 
# GMT modern mode bash template
# Date:    2020-07-24T12:48:37
# User:    yzhu
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
ori=/home1/yzhu/xin_code/2term/benchmark4_2d_TH500km_L2000km/disp/
lab1=/home1/yzhu/xin_code/2term/thrust15_fw50km_E30km_TH500km_L2000km_lab1e16/disp/
lab2=/home1/yzhu/xin_code/2term/thrust15_fw50km_E30km_TH500km_L2000km_lab1e15/disp/
lab3=/home1/yzhu/xin_code/2term/thrust15_fw50km_E30km_TH500km_L2000km_lab1e14/disp/
fem_lab1=/home1/yzhu/thrust2d/LAB/p_thrust15_W50km_E30km_M1e19K1e18_lab10km_1e16/
fem_lab2=/home1/yzhu/thrust2d/LAB/p_thrust15_W50km_E30km_M1e19K1e18_lab10km_1e15/
# fem2d=/home1/yzhu/thrust2d/benchmark/Benchmark_pollitz_L1000km_nongrav/
# fem2d_grav=/home1/yzhu/thrust2d/benchmark/Benchmark_pollitz_L1000km_grav/

factor=111.1949266
ntime=3
times=(0.5 5 50 )
R1=-150/250/-80/40
R2=-150/250/-50/50
R3=-150/250/-60/30
R4=-150/250/-30/30
pi=3.1415926
blue=black
red=black,4_3:0
ss=3D_vis_eq_W50D_L5D

gmt begin LAB_2d_effects_xin_M1e19K1e18 png
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 2x2 -Fs5i/2i -M0.3c/0.3c -A+o0.2c/4c
		gmt subplot set 0,0 -A'(a) Vertical'
                gmt basemap -JX? -R$R1 -Bxa50 -Bya20+l"Velocity (mm/yr)" -BWsen
                for(( i=1;i<="$ntime";i++ ))
                do
                        a=${times[i-1]};

                        awk '{print $1*'$factor',$3*1000}' ${ori}/xin_po${a}.vel | gmt plot -W1,gray+s            
                        awk '{print $1*'$factor',$3*1000}' ${lab1}/xin_po${a}.vel | gmt plot -W0.5,red  
                        awk '{print $1*'$factor',$3*1000}' ${lab2}/xin_po${a}.vel | gmt plot -W0.5,green  
                        awk '{print $1*'$factor',$3*1000}' ${lab3}/xin_po${a}.vel | gmt plot -W0.5,blue      
                        awk '{print $1*'$factor',$3*1000}' ${fem_lab1}/po${a}yr.vel | gmt plot -W1,red,4:3_0
                        awk '{print $1*'$factor',$3*1000}' ${fem_lab2}/po${a}yr.vel | gmt plot -W1,green,4:3_0
                done
                ###plot legend 
                # echo "50  5  500yr" | gmt text  
                echo "50  -3    50yrs" | gmt text  
                echo "50  -8     5yrs" | gmt text  
                echo "60  -30     0.5yr" | gmt text  
gmt plot -W1,gray << END 
120 -15
125 -15
130 -15
END

gmt plot -W1,red << END 
120 -30
125 -30
130 -30
END

gmt plot -W1,green << END 
120 -45
125 -45
130 -45
END
gmt plot -W1,blue << END 
120 -60
130 -60
END

echo "135  -15   Reference" | gmt text -F+f10p,gray+jML 
echo "135  -30   LAB (1e16)"| gmt text -F+f10p,red+jML 
echo "135  -45   LAB (1e15)"| gmt text -F+f10p,green+jML 
echo "135  -60   LAB (1e14)"| gmt text -F+f10p,blue+jML 
		gmt subplot set 1,0 -A'(b) Horizontal'
			# gmt basemap -JX? -R$R2 -Bxa5 -Bya0.2+l"Fault-normal Velocity (V/V@-o@-)" -BWSen
			gmt basemap -JX? -R$R2 -Bxa50+l"Distance (km)" -Bya10+l"Velocity (mm/yr)" -BWSen+t"Velocity (mm/yr)"

                for(( i=1;i<="$ntime";i++ ))
                do
                a=${times[i-1]};
                        # awk '{print $1*'$factor',$2*100}' po${a}.dis | gmt plot -W1,   
                        awk '{print $1*'$factor',$2*1000}' ${ori}xin_po${a}.vel | gmt plot -W1,gray+s  
                        awk '{print $1*'$factor',$2*1000}' ${lab1}xin_po${a}.vel | gmt plot -W0.5,red
                        awk '{print $1*'$factor',$2*1000}' ${lab2}xin_po${a}.vel | gmt plot -W0.5,green
                        awk '{print $1*'$factor',$2*1000}' ${lab3}xin_po${a}.vel | gmt plot -W0.5,blue
                        awk '{print $1*'$factor',$2*1000}' ${fem_lab1}/po${a}yr.vel | gmt plot -W1,red,4:3_0
                        awk '{print $1*'$factor',$2*1000}' ${fem_lab2}/po${a}yr.vel | gmt plot -W1,green,4:3_0
                done

		gmt subplot set 0,1 -A'(c) Vertical'
                gmt basemap -JX? -R$R3 -Bxa50 -Bya20+l"Displacement (cm)" -BwsEn
    
                for(( i=1;i<="$ntime";i++ ))
                do
                        a=${times[i-1]};
                        awk '{print $1*'$factor',$3*100}' ${ori}xin_po${a}.dis | gmt plot -W1,gray+s
                        awk '{print $1*'$factor',$3*100}' ${lab1}xin_po${a}.dis | gmt plot -W0.5,red
                        awk '{print $1*'$factor',$3*100}' ${lab2}xin_po${a}.dis | gmt plot -W0.5,green
                        awk '{print $1*'$factor',$3*100}' ${lab3}xin_po${a}.dis | gmt plot -W0.5,blue
                done

		gmt subplot set 1,1 -A'(d) Horizontal'
			# gmt basemap -JX? -R$R2 -Bxa5 -Bya0.2+l"Fault-normal Velocity (V/V@-o@-)" -BWSen
			gmt basemap -JX? -R$R4 -Bxa50+l"Distance (km)" -Bya10+l"Displacement (cm)" -BwSEn+t"Displacement (cm)"
                for(( i=1;i<="$ntime";i++ ))
                do
                        a=${times[i-1]}; 
                        awk '{print $1*'$factor',$2*100}' ${ori}xin_po${a}.dis | gmt plot -W1,gray+s
                        awk '{print $1*'$factor',$2*100}' ${lab1}xin_po${a}.dis | gmt plot -W0.5,red
                        awk '{print $1*'$factor',$2*100}' ${lab2}xin_po${a}.dis | gmt plot -W0.5,green
                        awk '{print $1*'$factor',$2*100}' ${lab3}xin_po${a}.dis | gmt plot -W0.5,blue
                done
		gmt subplot end 
gmt end 
