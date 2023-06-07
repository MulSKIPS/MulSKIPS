#!/usr/bin/bash
grep -a "Iter," log.out |awk '{print $8}' > time_zint.dat
grep -a "Iter, Time, Site" log.out |awk '{print $5}' > time.dat
grep -a "^ CountCrystal " log.out |awk '{print $2}' > time_countSi.dat
grep -a "^ CountCrystal " log.out |awk '{print $3}' > time_countGe.dat
grep -a "^ CountCrystal " log.out |awk '{print $4}' > time_countB.dat
grep -a "^ CountCov " log.out |awk '{print $2}' > time_countCl.dat
grep -a "^ CountCov " log.out |awk '{print $3}' > time_countH.dat
grep -a "^ xGesolid" log.out |awk '{print $2}' > time_x_Si.dat
grep -a "^ xGesolid" log.out |awk '{print $3}' > time_x_Ge.dat
grep -a "^ xGesolid" log.out |awk '{print $4}' > time_x_B.dat

paste time.dat time_zint.dat > time_z.dat

paste time_zint.dat time_x_Si.dat > zint_x_Si.dat
paste time_zint.dat time_x_Ge.dat > zint_x_Ge.dat
paste time_zint.dat time_x_B.dat > zint_x_B.dat
gplot -y 'Atomic fraction' -x 'Grown thickness (nm)' -Y 0:1 u '($1*0.04525):2' w l lw 2 ::: zint_x_*.dat
gplot -x 'Time (ns)' -y 'Grown thickness (nm)' u '1:($2*0.04525-2.172)' w l ::: time_z.dat
