#!/bin/bash

# Generate namelists for a 400-year trajectory run in 50 year increments.

# Rename restart file
if [ -e SPINUP_05201250_restart.nc ];then mv SPINUP_05201250_restart.nc TRAJ_05201250_restart.nc;fi
if [ -e SPINUP_05201250_restart_ice.nc ];then mv SPINUP_05201250_restart_ice.nc TRAJ_05201250_restart_ice.nc;fi


for Y in {000..350..050};do
    
    fname=TRAJ_${Y}_`printf "%03d" $((10#${Y}+50))`.namelist
    fname_ice=TRAJ_${Y}_`printf "%03d" $((10#${Y}+50))`.namelist_ice

    RSTRTSTP=`printf "%08d" $((10#05201250    + $((10#${Y} * 15*365)) ))`
    FRST_STP=`printf "%08d" $((10#${RSTRTSTP} + 1           ))`
    LAST_STP=`printf "%08d" $((10#${RSTRTSTP} + (50*15*365) ))`
    N_OFFSET=`printf "%08d" $((10#${Y}*15*365))`
    STRTDATE=`printf "%04d" $((10#951+10#${Y}))`0101

    echo "################################################################################"
    echo ${Y}
    echo "Restart file    : TRAJ_${RSTRTSTP}_restart.nc"
    echo "Ice restart file: TRAJ_${RSTRTSTP}_restart_ice.nc"
    echo "namelist file   : ${fname}"
    echo "Start date      : ${STRTDATE}"
    echo "First time step : ${FRST_STP}"
    echo "Final time step : ${LAST_STP}"
    echo "First traj step : ${N_OFFSET}"
    
    cp namelist.TRAJ_TEMPLATE ${fname}
    cp namelist_ice.TRAJ_TEMPLATE ${fname_ice}
    sed -i -e "s#STRTDATE#${STRTDATE}#g" ${fname}
    sed -i -e "s#FRST_STP#${FRST_STP}#g" ${fname}
    sed -i -e "s#LAST_STP#${LAST_STP}#g" ${fname}
    sed -i -e "s#RSTRTSTP#${RSTRTSTP}#g" ${fname}
    sed -i -e "s#RSTRTSTP#${RSTRTSTP}#g" ${fname_ice}
    sed -i -e "s#N_OFFSET#${N_OFFSET}#g" ${fname}
done
