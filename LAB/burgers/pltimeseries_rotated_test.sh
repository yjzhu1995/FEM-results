#! /bin/bash 
# GMT modern mode bash template
# Date:    2020-07-24T12:48:37
# User:    yzhu
# Purpose: Purpose of this script
export GMT_SESSION_NAME=$$	# Set a unique session name
lab=Haida_Gwaii_3layer_E7km_H300M1e20K1e16lowrid_M1e20
lab2=Haida_Gwaii_3layer_E7km_H300M1e20K1e16lowrid20_M1e20
lab3=Haida_Gwaii_3layer_E7km_H300M1e20K1e16_M1e20
ori=/home1/yzhu/xin_code/2term/benchmark4_2d_TH500km_L2000km/timeseries/
# ori_maxwell=/home1/yzhu/thrust2d/LAB/p_thrust15_W50km_E30km_1e18/
xin=/home1/yzhu/xin_code/2term/${lab}/disp/
xin2=/home1/yzhu/xin_code/2term/${lab2}/disp/
xin3=/home1/yzhu/xin_code/2term/${lab3}/disp/
GPS=/home1/yzhu/yzhu/FEM-results/VIsland/VI_timeseriesplot/

# fem2d=/home1/yzhu/thrust2d/benchmark/Benchmark_pollitz_L1000km_nongrav/
# fem2d_grav=/home1/yzhu/thrust2d/benchmark/Benchmark_pollitz_L1000km_grav/
#sta=( QUAD  TFNO  UCLU ) #NANO ALBH  PGC5  PTRF    BAMF  BCOV   BCPH  ELIZ  GLDR  HOLB  NTKA  PTAL   ###not available stations: NEAH BCES P415 SC02  SEAT
sta=(HOLB BCPH BCOV ELIZ BCOV NTKA GLDR QUAD TFNO UCLU BAMF PTAL NANO PTRF CLRS JORD SC04 PGC5 BCNS ALBH ESQM WSLR BCVC BCSF BCMR BCLC CHWK WILL)
# ALBH BAMF BCLC BCMR BCNS BCOV BCPH   BCSF BCVC CHWK CLRS ELIZ ESQM GLDR   HOLB JORD NANO NTKA PGC5 PTAL PTRF    QUAD SC04 TFNO UCLU WILL WOST WSLR
comp=(3 4 2)
gpscomp=(2 3 4)
factor=111.1949266
ntime=3
times=(0.5 5 50 )
R2=2012/2020/-1/3
R1=2012/2020/-2/1
R3=2012/2020/0/1
R4=2012/2020/-5/1
pi=3.1415926
#R=(2010/2020/-15/15 2010/2020/-15/15 2010/2020/-10/10)
R=2010/2020/-10/20
B=(ya5 ya5 ya2)
year=2012.82
gmt begin Timeseries_rotated_3layer_Burgers_difflowrid png A+m0.5c
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 4x7 -Fs5i/2i -M0.3c/0.3c -A+o0.2c/0.2c
        for((i=0;i<4;i++))
        do
            for((j=0;j<7;j++))
            do
                let n=${i}*7+${j}
                name=''${sta[n]}''
                gmt subplot set ${i},${j} -A$name
                        gmt basemap -JX? -R$R -Bxa2 -Bya5 -BWSen
                        awk 'NR>1{print $1,$'${gpscomp[0]}'}' ${GPS}/${sta[n]}.rotate | gmt plot -Sc0.1c -W1,blue -Gblue
                        awk '{print $1+'$year',$2*1000}' ${xin}/${sta[n]}_total.rotate | gmt plot -W2,red+s 
                        awk '{print $1+'$year',$2*1000}' ${xin2}/${sta[n]}_total.rotate | gmt plot -W2,green+s 
                        awk '{print $1+'$year',$2*1000}' ${xin3}/${sta[n]}_total.rotate | gmt plot -W2,black+s                            
            done
        done 
		gmt subplot end 
gmt end 
