#! /bin/bash

# Note, options '--report-bindings -x MALLOC_MMAP_MAX_=-1 -x MALLOC_TRIM_THRESHOLD_=33554432 --mca mpi_paffinity_alone' used for 'mpirun' command obtained from A. Coward

#SBATCH -J PT_TL
#SBATCH -o PT_TL.out
#SBATCH -N 4
#SBATCH -n 64
#SBATCH -t 23:59:59
#SBATCH -p low

module purge
module load shared
module load gcc
module load intel/compiler/64/14.0/2013_sp1.2.144
module load openmpi/intel/64/1.6.5
module load netcdf/intel/64/4.3.0
module load slurm/16.05.8

runname=PT_TL
TAM=tl

ln -sf namelist.${runname} namelist

if [ ! -d ${runname}_output ]; then mkdir ${runname}_output; fi

time `which mpirun` --report-bindings -x MALLOC_MMAP_MAX_=-1 -x MALLOC_TRIM_THRESHOLD_=33554432 --mca mpi_paffinity_alone 1 -np 64 ./nemo_tam.exe
rm 000000_${TAM}_trajectory_????.nc

############################ stitch outputs together in space ######################################
############################    and move to output folder     ######################################
for filename in PTTAM_output_????????_0000.nc; do \ 
    ../../../../TOOLS/REBUILD_NEMO/rebuild_nemo ${filename/_0000.nc} 64; \
	if [ -e ${filename/_0000.nc}.nc ]; #checks stitch has been successful before blindly deleting everything
    then rm ${filename/_0000.nc/}_????.nc;
    fi;
done;
mv PTTAM_output_????????.nc ${runname}_output

echo "stitch done! ncrcat-ing..."


############################## stitch outputs together in time with ncrcat ########################
ncrcat -d t,0,,1 ${runname}_output/PTTAM_output_*.nc ${runname}_output/${runname}_output.nc 
#Andrew Coward's mobilis ncrcat setup to combine outputs
if [ -e ${runname}_output/${runname}_output.nc ];
then rm ${runname}_output/PTTAM_output_*.nc
fi










