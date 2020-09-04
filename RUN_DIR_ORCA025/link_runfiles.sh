#!/bin/bash


ORCA025INPUT_dir=/wherever/you/keep/your/ORCA025INPUT
bin_dir=/your_NEMO_3.4_source_code_directory/NEMOGCM/CONFIG/ConTraPTION/BLD/bin

echo "ORCA025INPUT is located at ${ORCA025INPUT_dir}"
echo "BLD/bin is located at ${bin_dir}"
echo "continue? link NEMO files in working directory? (y/n)"
read linkswitch

if [ ${linkswitch} != y ]
then exit 1
else

    ln -s ${ORCA025INPUT_dir}/bathy_level.nc .
    ln -s ${ORCA025INPUT_dir}/bathy_meter.nc .
    ln -s ${ORCA025INPUT_dir}/bfr_coef.nc .
    ln -s ${ORCA025INPUT_dir}/coordinates.nc .
    ln -s ${ORCA025INPUT_dir}/data_1m_potential_temperature_nomask.nc .
    ln -s ${ORCA025INPUT_dir}/data_1m_salinity_nomask.nc .
    ln -s ${ORCA025INPUT_dir}/domain_def.xml . #
    ln -s ${ORCA025INPUT_dir}/field_def.xml . #
    ln -s ${ORCA025INPUT_dir}/geothermal_heating.nc .
    ln -s ${ORCA025INPUT_dir}/ice_in . #
    ln -s ${ORCA025INPUT_dir}/iodef.xml . #
    ln -s ${ORCA025INPUT_dir}/K1rowdrg.nc .
    ln -s ${ORCA025INPUT_dir}/M2rowdrg.nc .
    ln -s ${ORCA025INPUT_dir}/mask_itf.nc .
    ln -s ${ORCA025INPUT_dir}/runoff_1m_nomask.nc . #
    ln -s ${ORCA025INPUT_dir}/sss_1m.nc . #
    
    ################################################################################
    # Weighting files to interpolate forcing
    ln -s ${ORCA025INPUT_dir}/weights_grid03_bilinear_orca025.nc . #
    
    ################################################################################
    # CORE Normal-year forcing files
    ln -s ${ORCA025INPUT_dir}/LWDN_MOD_15JUNE2009.nc . #
    ln -s ${ORCA025INPUT_dir}/Q_10_MOD_15JUNE2009.nc . #
    ln -s ${ORCA025INPUT_dir}/SNOW_15JUNE2009.nc . #
    ln -s ${ORCA025INPUT_dir}/SWDN_MOD_15JUNE2009.nc . #
    ln -s ${ORCA025INPUT_dir}/TPRECIP_15JUNE2009.nc . #
    ln -s ${ORCA025INPUT_dir}/T_10_MOD_15JUNE2009.nc . #
    ln -s ${ORCA025INPUT_dir}/U_10_MOD_15JUNE2009.nc . #
    ln -s ${ORCA025INPUT_dir}/V_10_MOD_15JUNE2009.nc . #
    
    ################################################################################
    # Executables 
    ln -s ${bin_dir}/nemo.exe .
    ln -s ${bin_dir}/nemo_tam.exe .
    ln -s ${bin_dir}/../../../../TOOLS/REBUILD_NEMO/rebuild_nemo .
    ln -s ${bin_dir}/../../../../TOOLS/REBUILD_NEMO/rebuild_nemo.exe .
fi

