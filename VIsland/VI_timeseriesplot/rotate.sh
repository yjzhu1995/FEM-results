#! /bin/bash/

sta=(ALBH BAMF BCLC BCMR BCNS BCOV BCPH BCSF BCVC CHWK CLRS ELIZ ESQM GLDR HOLB JORD NANO NTKA PGC5 PTAL PTRF QUAD SC04 TFNO UCLU WILL WOST WSLR)
n=${#sta[*]}

alon=-122
alat=48
blon=-128
blat=51
C=${blon}/${blat}
E=${alon}/${alat}

project -C$C -E$E -G10 -Q > projection_trace
theta=$(echo "$alon $alat $blon $blat" | awk '{print atan2(($3-$1)*cos($2*'$pi'/360+$4*'$pi'/360),$4-$2)}')
echo "$theta"

for((i=0;i<n;i++))
do
	tsfile=${sta[i]}.txt
	awk '{print $1,-$2*cos('$theta')+$3*sin('$theta'),$3*cos('$theta')+$2*sin('$theta')}' $tsfile > ${sta[i]}.rotate
done

