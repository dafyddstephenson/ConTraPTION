import netCDF4 as nc
import numpy as np

################################################################################
'''
Python script to generate an initial passive tracer distribution to be read
by NEMOTAM.

The user should provide the location of the "subbasins.nc" file, which provides 
masks for the ocean basins and is found inside the ORCA2_LIM_v3.4.tar file, 
which should have been unpacked to the directory /[PATH]/ORCA2_INPUT.

"subbasins.nc" contains the variables
atlmsk       (Atlantic ocean mask)
atlmsk_nomed (Atlantic ocean mask without Mediterranean Sea)
pacmsk       (Pacific  ocean mask)
indmsk       (Indian   ocean mask)
navlat       (Latitude)
navlon       (Longitude)

As an example, the script creates an initial distribution injecting tracer 
into the surface layer of the Atlantic between the latitudes of 40 and 60N
and longitudes of 20 and 60W.

NOTE: lat and lon here are approximate, and more advanced use should include
T-point locations from the mesh_mask.nc produced by the model with nn_msh = 1
in the namelist.
'''
################################################################################
ORCA2INPUT_dir='../ORCA2_INPUT'
BSNS=nc.Dataset(ORCA2INPUT_dir+'subbasins.nc')

atlmsk=BSNS.variables['atlmsk' ][:]
lon   =BSNS.variables['navlon'][:]
lat   =BSNS.variables['navlat'][:]

Tinit=np.zeros((31,149,182))
Tinit[0,:,:][(lat>40) & (lat<60) & (lon<-20) & (lon>-60)]=1
Tinit*=atlmsk

with nc.Dataset('PT_init.nc','w') as PTIC:
    PTIC.createDimension('z',31)
    PTIC.createDimension('y',149)
    PTIC.createDimension('x',182)
    
    PT0=PTIC.createVariable('Tinit',np.float64,('z','y','x'))
    PT0[:]=Tinit
