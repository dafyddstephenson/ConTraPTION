#! /bin/bash

#Creates several consecutive job submission scripts to spin model up from 0. Each run (50 years) produces a restart file which is used by the next run. 
#The jobs can be submitted with a dependency, or the SPINUP_template.sh file can be adapted to submit the next job.
#Make sure cn_ocerst_in, nn_it000 and nn_date0 in namelist.TRAJ match the final output of the spinup

for number in {2..19};
do
t0=$((1+(273750*(number-1))))
tend=$((273750*number))
yr=$(((number-1)*50))
cp namelist.SPINUP_template namelist.SPINUP${number}
cp namelist_ice.SPINUP_template namelist_ice.SPINUP${number}
cp SPINUP_template.sh SPINUP${number}.sh
sed -i -e "s#nn_it000    = ???????#nn_it000    = ${t0}#g" namelist.SPINUP${number};
sed -i -e "s#nn_itend    = ???????#nn_itend    = ${tend}#g" namelist.SPINUP${number};
sed -i -e "s#nn_date0    = YYYY0101#nn_date0    = 0${yr}0101#g" namelist.SPINUP${number};
restartnum=`printf '%08i' $((273750*(number-1)))`
sed -i -e "s#SPINUP_IN_restart#SPINUP_${restartnum}_restart#g" namelist.SPINUP${number};
sed -i -e "s#SPINUP_IN_restart_ice#SPINUP_${restartnum}_restart_ice#g" namelist_ice.SPINUP${number};
sed -i -e "s#spinnum=NUM#spinnum=${number}#g" SPINUP${number}.sh;
done