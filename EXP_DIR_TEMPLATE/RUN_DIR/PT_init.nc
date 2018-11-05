This should be a netCDF file with the following characteristics:

netcdf PT_init {
dimensions:
	x = 182 ;
	y = 149 ;
	z = 31 ;
variables:
	double Tinit(z, y, x) ;
}

Tinit represents tracer concentration (i.e. 1 if a grid cell is saturated with tracer and 0 if it contains no tracer)