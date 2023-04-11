#!/usr/bin/bash
grep " TE" log.out | awk '{print $8}' > tmp.L
grep " TE" log.out | awk '{print $4}' > tmp.time
grep " TE" log.out | awk '{print $7}' > tmp.tmax
paste tmp.time tmp.tmax > time_tmax.dat
paste tmp.time tmp.L > time_L.dat

grep "Results have been moved" log.out |awk -F 'KMC_t_' '{print $2}' > tmp.time
grep -a -A 2 "I00000009" log.mulskips.out |grep Iter |awk '{print $NF}' > tmp.zint
grep "Solid phase atoms in MulSKIPS" log.out |awk -F ': ' '{print $2}' > tmp.solidKMC
grep "Solid phase atoms:" log.out |awk -F ': ' '{print $2}' > tmp.solidFEM

paste tmp.time tmp.zint > time_zint.dat
paste tmp.time tmp.solidKMC > time_solidKMC.dat
paste tmp.time tmp.solidFEM > time_solidFEM.dat

grep -a "b(1,1)" log.mulskips.out |awk '{print $2}' > tmp.tempkmc
paste tmp.time tmp.tempkmc > time_tempminkmc.dat

grep -a "Iter" log.mulskips.out |awk '{print $8}' > time_zint_all.dat
grep -a occupied log.mulskips.out |awk '{print $NF}' > time_solidKMC_all.dat
grep -a "^ CountCrystal " log.mulskips.out |awk '{print $2}' > time_solidKMC_all_Si.dat
grep -a "^ CountCrystal " log.mulskips.out |awk '{print $3}' > time_solidKMC_all_Ge.dat
grep -a "^ xLiquid " log.mulskips.out |awk '{print $2}' > time_liquidKMC_Si.dat
grep -a "^ xLiquid " log.mulskips.out |awk '{print $3}' > time_liquidKMC_Ge.dat
grep -a "Atom fraction from last check:" log.mulskips.out |awk '{print $6}' > time_x_all_Si.dat
grep -a "Atom fraction from last check:" log.mulskips.out |awk '{print $7}' > time_x_all_Ge.dat
python3 smooth.py 51 # smooth xSi and xGe

rm tmp*

#gplot -o time_tmax.png -x "Time (ns)" -y "Max T (K)" w lp ::: time_tmax.dat
#gplot -o time_zint.png -x "Time (ns)" -y "S-L interface position (alat)" w lp ::: time_zint.dat

#gplot -x "Time (ns)" -y "Max T (K)" w lp ::: time_tmax.dat
#gplot -x "Time (ns)" -y "S-L interface position (alat)" w lp ::: time_zint.dat

tar cvfz todownload.tar.gz KMC_t_*/I*000.xyz log* Mul*py time_* zint_*dat XGeS_init* #KMC_t_*/InterpXGeS* 

