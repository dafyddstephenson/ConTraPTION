#!/bin/bash

ORCA2INPUT_dir=../ORCA2_INPUT
BLD_dir=../../BLD/bin

ln -s ${ORCA2INPUT_dir}/ahmcoef .
ln -s ${ORCA2INPUT_dir}/bathy_level.nc .
ln -s ${ORCA2INPUT_dir}/bathy_meter.nc .
ln -s ${ORCA2INPUT_dir}/bathy_updated.nc .
ln -s ${ORCA2INPUT_dir}/chlorophyll.nc .
ln -s ${ORCA2INPUT_dir}/coordinates.nc .
ln -s ${ORCA2INPUT_dir}/data_1m_potential_temperature_nomask.nc .
ln -s ${ORCA2INPUT_dir}/data_1m_salinity_nomask.nc .
ln -s ${ORCA2INPUT_dir}/geothermal_heating.nc .
ln -s ${ORCA2INPUT_dir}/K1rowdrg.nc .
ln -s ${ORCA2INPUT_dir}/M2rowdrg.nc .
ln -s ${ORCA2INPUT_dir}/mask_itf.nc .
ln -s ${ORCA2INPUT_dir}/ncar_precip.15JUNE2009_orca2.nc .
ln -s ${ORCA2INPUT_dir}/ncar_rad.15JUNE2009_orca2.nc .
ln -s ${ORCA2INPUT_dir}/q_10.15JUNE2009_orca2.nc . 
ln -s ${ORCA2INPUT_dir}/runoff_core_monthly.nc .
ln -s ${ORCA2INPUT_dir}/t_10.15JUNE2009_orca2.nc .
ln -s ${ORCA2INPUT_dir}/u_10.15JUNE2009_orca2.nc .
ln -s ${ORCA2INPUT_dir}/v_10.15JUNE2009_orca2.nc .
ln -s ${ORCA2INPUT_dir}/sss_data.nc .
ln -s ${ORCA2INPUT_dir}/sst_data.nc .
ln -s ${ORCA2INPUT_dir}/subbasins.nc .

ln -s ${BLD_dir}/nemo.exe .
ln -s ${BLD_dir}/nemo_tam.exe .

