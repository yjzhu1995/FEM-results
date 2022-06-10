#!/usr/bin/env -S bash -e

export GMT_SESSION_NAME=123	# Set a unique session name
model=Haida_Gwaii_3layer_E10km_H10kmM1e16_M1e20
R=-129/-122/48/51
J=M3i
gps=/home1/yzhu/AGU2019/GPS/disp_VI.txt
xin1=/home1/yzhu/xin_code/2term/${model}/disp/
xin2=/home1/yzhu/xin_code/2term/Haida_Gwaii_3layer_E10km_H200kmM2e17_M1e20/disp/
xin3=/home1/yzhu/xin_code/2term/Haida_Gwaii_3layer_E10km_H400kmM4e17_M1e20/disp/
arrow1=0.2c+ea+h0.25+p1.5p,black+n30 #0.2c represent arrow length
arrow2=0.2c+ea+h0.25+p1.5p,red+n30
arrow3=0.2c+ea+h0.25+p1.5p,green+n30
arrow4=0.2c+ea+h0.25+p1.5p,yellow+n30
R2=(2010/2020/-25/5 2010/2020/-5/25)
B=(ya5 ya5 ya2)
year=2012.82
GPSts=/home1/yzhu/yzhu/FEM-results/VIsland/VI_timeseriesplot/
sta=(BCPH NTKA ESQM)
comp=(3 4 2)
gpscomp=(2 3 4)
factor=111.1949266
ntime=3
times=(0.5 5 50 )
n=${#sta[*]}

gmt begin Quicklook_tradeoff png A+m0.5c
	gmt gmtset MAP_GRID_CROSS_SIZE_PRIMARY 0.05i MAP_FRAME_TYPE plain MAP_TICK_LENGTH_PRIMARY -0.15 MAP_FRAME_PEN 1p  MAP_TICK_PEN 1p
	gmt gmtset FONT_LABEL 12 FONT_ANNOT_PRIMARY 11 MAP_ANNOT_OFFSET_PRIMARY 0.2 
	gmt gmtset FONT_TITLE 12 FONT_TAG 12,
	gmt gmtset MAP_ANNOT_ORTHO we MAP_LABEL_OFFSET 0.3
	gmt subplot begin 3x3 -Fs5i/2i -M0.3c/0.3c -A+o0.2c/3c
		gmt subplot set 0,0 -A'Mapview of HG'
			# gmt basemap -JM5i -R$R  -Bxa5 -Bya2 -BWSen
			gmt coast -JM7i -R$R  -Dh -A100 -G244/243/239 -Wthinnest,244/243/239 -S167/194/223   -Bxa5 -Bya2 -BWSen 
			# awk '{print $3,$2,$6,$5,$8,$7,0,$1}' $gps | gmt velo   -JM7i -R$R -Se0.05c/0.95/8 -A$arrow -Gblue -W1,blue
			awk '{print $1,$2,$3*1000,$4*1000,0,0,0}' ${xin1}/xin_po8.dis | gmt velo -JM7i -R$R -Se0.05c/1/0 -A$arrow1 -Gblack -W1.5,black
			awk '{print $1,$2,$3*1000,$4*1000,0,0,0}' ${xin2}/xin_po8.dis | gmt velo -JM7i -R$R -Se0.05c/1/0 -A$arrow2 -Gred -W1.5,red
			awk '{print $1,$2,$3*1000,$4*1000,0,0,0}' ${xin3}/xin_po8.dis | gmt velo -JM7i -R$R -Se0.05c/1/0 -A$arrow3 -Ggreen -W1.5,green
		
			echo "-133 48 10 0 0 0 0" | gmt velo -Se0.1c/1/0 -A$arrow -Gblue -W0.5,blue -JM5i -R$R
			echo "-132.5 48.5 10 mm" | gmt text  -JM5i -R$R
		for((i=1;i<3;i++))
        do
            for((j=0;j<n;j++))
            do
                name=''${sta[j]}''
                gmt subplot set ${i},${j} -A$name
                        gmt basemap -JX? -R${R2[i-1]} -Bxa5 -B${B[i-1]} -BWSen
                        awk 'NR>1{print $1,$'${gpscomp[i-1]}'}' ${GPSts}/${sta[j]}.txt | gmt plot -Sc0.1c -W1,blue -Gblue
                        awk '{print $1+'$year',$'${comp[i-1]}'*1000}' ${xin1}/${sta[j]}_total.ts | gmt plot -W2,black+s 
                        awk '{print $1+'$year',$'${comp[i-1]}'*1000}' ${xin2}/${sta[j]}_total.ts | gmt plot -W2,red+s 
                        awk '{print $1+'$year',$'${comp[i-1]}'*1000}' ${xin3}/${sta[j]}_total.ts | gmt plot -W2,green+s 
            done
        done 
            #awk '{ print $3,$2,$1}'  ssite.txt | gmt text -D-0.5/-0.5+v1p,red
	gmt subplot end
gmt end 