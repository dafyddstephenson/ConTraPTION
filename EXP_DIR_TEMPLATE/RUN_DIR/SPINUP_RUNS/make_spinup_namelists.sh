#!/bin/bash

for Y in {000..900..050};do
    
    fname=SPINUP_${Y}_`printf "%03d" $((10#${Y}+50))`.namelist
    fname_ice=SPINUP_${Y}_`printf "%03d" $((10#${Y}+50))`.namelist_ice
    
    if [ ${Y} == 000 ];then RSTRT_TF=.false.;else RSTRT_TF=.true.;fi
    
    RSTRTSTP=`printf "%08d" $((10#${Y} * 15*365))`
    FRST_STP=`printf "%08d" $((10#${RSTRTSTP} + 1           ))`
    LAST_STP=`printf "%08d" $((10#${RSTRTSTP} + (50*15*365) ))`
    STRTDATE=`printf "%04d" $((10#0001+10#${Y}))`0101

    echo "################################################################################"
    echo ${Y}
    if [ ${Y} -gt 0 ];then
	echo "Restart file    : SPINUP_${RSTRTSTP}_restart.nc"
	echo "Ice restart file: SPINUP_${RSTRTSTP}_restart_ice.nc"
    fi 
    echo "From restart    : ${RSTRT_TF}"
    echo "namelist file   : ${fname}"
    echo "Start date      : ${STRTDATE}"
    echo "First time step : ${FRST_STP}"
    echo "Final time step : ${LAST_STP}"
    
    cp namelist.SPINUP_TEMPLATE ${fname}
    cp namelist_ice.SPINUP_TEMPLATE ${fname_ice}
    sed -i -e "s#STRTDATE#${STRTDATE}#g" ${fname}
    sed -i -e "s#FRST_STP#${FRST_STP}#g" ${fname}
    sed -i -e "s#LAST_STP#${LAST_STP}#g" ${fname}
    sed -i -e "s#RSTRT_TF#${RSTRT_TF}#g" ${fname}
    sed -i -e "s#RSTRTSTP#${RSTRTSTP}#g" ${fname}
    sed -i -e "s#RSTRTSTP#${RSTRTSTP}#g" ${fname_ice}
done
