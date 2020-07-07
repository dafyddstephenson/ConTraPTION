#!/bin/bash

ORCA2INPUT_dir=/wherever/you/keep/your/ORCA2INPUT
bin_dir=/your_NEMO_3.4_source_code_directory/NEMOGCM/CONFIG/TAM_ORCA2/BLD/bin

################################################################################
# NOTE: can find ORCA2INPUT files at
#    https://zenodo.org/record/1471702#.Xu9ZemozZcA
# Under ORCA2_LIM_v3.4.tar

echo "ORCA2INPUT is located at ${ORCA2INPUT_dir}"
echo "BLD/bin is located at ${bin_dir}"
echo "continue? link NEMO files in working directory? (y/n)"
read linkswitch

if [ ${linkswitch} != y ]
then exit 1
else

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
    ln -s ${ORCA2INPUT_dir}/sss_data.nc .
    ln -s ${ORCA2INPUT_dir}/sst_data.nc .

    ################################################################################
    # CORE Normal-year forcing files

    ln -s ${ORCA2INPUT_dir}/ncar_precip.15JUNE2009_orca2.nc .
    ln -s ${ORCA2INPUT_dir}/ncar_rad.15JUNE2009_orca2.nc .
    ln -s ${ORCA2INPUT_dir}/q_10.15JUNE2009_orca2.nc . 
    ln -s ${ORCA2INPUT_dir}/runoff_core_monthly.nc .
    ln -s ${ORCA2INPUT_dir}/t_10.15JUNE2009_orca2.nc .
    ln -s ${ORCA2INPUT_dir}/u_10.15JUNE2009_orca2.nc .
    ln -s ${ORCA2INPUT_dir}/v_10.15JUNE2009_orca2.nc
    ################################################################################
    #Executables
    ln -s ${bin_dir}/nemo.exe .
    ln -s ${bin_dir}/nemo_tam.exe .
    ln -s ${bin_dir}/../../../../TOOLS/REBUILD_NEMO/rebuild_nemo .
fi

