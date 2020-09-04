# ConTraPTION: Configuration for Tracking Passive Tracers Injected into the Ocean with NEMOTAM
## Overview
This is a NEMO *configuration* which can be used to track passive tracers in NEMO with NEMOTAM (see *Stephenson et al.*, 2020). 

NEMOTAM is the tangent-linear and adjoint model (TAM) counterpart to the nonlinear NEMO primitive equation solving ocean circulation model (c.f. *Vidard et al.*, 2011). NEMO is run first to generate a "trajectory" (containing complete information about the ocean state), then (in the passive mode of ConTraPTION) NEMOTAM can be run offline, using this information to track the probable origins (backward/adjoint mode) or destinations (forward/tangent-linear mode) of a passive tracer injected into the NEMO simulation. The tracer can be optionally removed at the surface to produce a record of ventilation.

## Installation
### Installing NEMO and ConTraPTION
This repository requires existing NEMO v3.4 (and NEMOTAM) installs. These can be installed using:

`svn co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-3.4`

This provides the source code for the model and reference configurations, on which the configuration used here should be based.

For a first-time NEMO install, the user will need to set up their machine's architecture. This will either already exist somewhere inside `NEMOGCM/ARCH` (e.g. `ARCH/<DIR>/arch-<ARCHITECTURE>.fcm`) or will have to be created. Instructions for creating/modifying this file are found [here](http://forge.ipsl.jussieu.fr/nemo/wiki/Users/ModelInstall#Setupyourarchitectureconfigurationfile).

At this point, the passive configuration can be created. To do this, run the following command in the `CONFIG` directory:

`./makenemo -d "OPATAM_SRC LIM_SRC_2 OPA_SRC" -n ConTraPTION -m <ARCHITECTURE> add_key "key_mpp_mpi  key_mpp_rep key_nosignedzero key_tam key_diainstant" del_key "key_zdfddm key_iomput"`

This creates a new ORCA2 (two degree resolution) configuration with ocean (`OPA_SRC`) and sea-ice (`LIM_SRC_2`) dynamics, and TAM (`OPATAM_SRC`, `key_tam`) compatibility. 
(NOTE: customise line as required, e.g. regarding parallel processing - `key_mpp_mpi` and `key_mpp_rep`. Replace <ARCHITECTURE> with the relevant substring from the arch. filename)

If you are interested in running the model at eddy-permitting resolution (ORCA025), compile again using:

`./makenemo -n ConTraPTION add_key "key_orca_r025=75 key_dynldf_c2d" del_key "key_orca_r2 key_diaeiv key_traldf_eiv key_dynldf_c3d`
 
This configuration can now be modified to include passive tracer-related subroutines. To do so, move the contents of `MY_SRC` from this repository into the `NEMOGCM/CONFIG/ConTraPTION/MY_SRC` directory and recompile using

`./makenemo -n ConTraPTION`

If running in parallel, the `rebuild_nemo` tool is also required, as the global domain is split into tiles and distributed among the CPUs. Each CPU produces a distinct output file, and these need to be stitched together.

Go to 
`NEMOGCM/TOOLS`

and run the make file:
`maketools -n REBUILD_NEMO -m <ARCHITECTURE>`

### Obtaining and linking forcing and other model input files

To run the model, additional files are required (for surface forcing etc.)

CORE normal-year (repeat annual) forcing (*Large & Yeager*, 2004) is provided by the NEMO Consortium for the ORCA2 configuration [here](https://doi.org/10.5281/zenodo.1471702).

(I have equivalent files for ORCA025 and will archive them online in the coming months. Until then, please contact me directly.)

Unpack this archive into a directory named `ORCA2_INPUT` (or `ORCA025_INPUT`, if you have the files for ORCA025). This will be at the same level as your experiment directory and can be anywhere you choose.

### Creating an experiment directory
At the same level as `ORCA2_INPUT` (or `ORCA025_INPUT`), create an experiment directory (e.g. `RUN_DIR`). Copy the contents of the corresponding template experiment directory from this repository into this directory. The bash script `link_runfiles.sh` creates softlinks to all files necessary to run NEMO and NEMOTAM. Edit it to provide the location of your `ORCA2_INPUT` (or `ORCA025_INPUT`) and `BLD/bin`.

### Directory structure

You should now have a NEMO configuration with the following structure

- `/home/username/NEMO/dev_v3_4_STABLE_2012/NEMOGCM/CONFIG/`
  - `ConTraPTION/`
    - `BLD/`
    - `MY_SRC/`
    - `WORK/`
    - `EXP00/`
    
and an experiment area with the structure    

- `/<PATH TO EXPERIMENT AREA>/`
  - `ORCA2_INPUT`
  - `RUN_DIR`

## Running passive tracer experiments


### Spinning up
To begin without a restart file, the model should first be spun up from rest to provide a start point for all experiments. A template namelist for ORCA2 (`namelist.SPINUP_TEMPLATE`) is included in `RUN_DIR_ORCA2`. However, as it is necessary to spin up for many hundreds of years, it is assumed ORCA025 users will already have access to a restart file (again, I will archive one online in the coming months and update this document. If one is required before this, please contact me).

The ORCA2 spinup run should produce the files `SPINUP_????????_restart.nc` and `SPINUP_????????_restart_ice.nc`, used as start points for the trajectory.

### Running a trajectory
After spinning up, the nonlinear model (NEMO, rather than NEMOTAM) should be run for the desired length of the experiment. Namelists can be found in 
`RUN_DIR/namelist.TRAJ_TEMPLATE` and `RUN_DIR/namelist_ice.TRAJ_TEMPLATE`.

**At the beginning of the namelist file:**

Be sure to set `ln_rstart = .true.`, and set `nn_it000`, `nn_date0` and `cn_ocerst_in` in accordance with the `SPINUP_*.nc` files generated above (`nn_it000` will be the number in the restart filename plus 1; `nn_date0` will be given by the variable `ndastp` in the restart file).

**At the end of the namelist file**:

Be sure TAM trajectory output is produced with `ln_trjhand = .true.`. The trajectory output directory / filename prefix are set by `cn_dirtrj`. All other output files (outside this directory) are not required by `NEMOTAM`. 

Trajectory outputs can be very large. If running a passive simulation in a limited area such that a global trajectory is not required, set `ln_pt_regional` to `.true.` and set the northeast and southwest corners of a quadrilateral domain of interest using `rn_{NE,SW}pt{lat,lon}`. Trajectory output files will only be produced in this region. For the TAM run, they will only be read in this region. (Outside of this region, the initial time step is read by NEMOTAM repeatedly for the entire run)

**Other parameters**:
Other options are detailled in [the NEMO 3.4 manual](http://forge.ipsl.jussieu.fr/little_nemo/export/44/vendor/nemo/current/DOC/NEMO_book.pdf).

### Running the TAM:
**Initial passive tracer distribution**:

The initial passive tracer distribution for a forward run should be saved (as concentrations between 0 and 1) to the variable `pt0_tl` in a netCDF file (e.g. `PT_init.nc`), which the template namelist points to (under `cn_tam_input`).

In the tangent-linear model, concentrations in `PT_init.nc` are propagated forward to produce a future **concentration** distribution. In the adjoint model, the concentrations in `PT_init.nc` are propagated backward to produce a past **volume** distribution (the distinction is inherent to the TAM framework, see *Stephenson et al., 2020*). 

**At the beginning of the namelist file**:

Set `nn_it000 = 1` (this is safe even if beginning later in the trajectory). 
For the tangent-linear model, `nn_itend` should then be set to the desired number of time steps in the run (up to the length of the trajectory). In the adjoint model, `nn_itend` determines the point from which model will run backwards.

**Surface ventilation**:

- **namsbc\_ssr**
  - The namelist parameter `nn_sstr` determines whether the passive tracer is removed at the surface
  - The namelist parameter `rn_dqdt` sets the time scale of the surface restoring scheme. The default, -40Wm^-2K^-1, corresponds to 60 days over a 50 m mixed layer.

**At the end of the namelist file:**

- **namtrj** (trajectory reading options)
  - `cn_dirtrj` should match that set in the trajectory namelist, specifying the location of the trajectory output.
  - `rn_rdttrj` is the timestep length used in the trajectory. This allows NEMO and NEMOTAM to be run using different timesteps, if necessary.
  - `nn_ittrjoffset` determines how many trajectory time-steps to "skip", for example to start a tangent-linear run later.
- **namtst\_tam** (general NEMOTAM options)
  - `ln_swi_opatam` determines TAM mode. `4` for passive forward and `5` for passive backward.
  - `cn_tam_input` is the initial tracer distribution file (in passive mode)
  - `ln_tl_eiv` activates eddy-induced velocities in NEMOTAM in the ORCA2 configuration (following *Gent & McWilliams*, 1990)
- **namtl_trj** (NEMOTAM output options)
  - `ln_trjwri_tan` whether output is desired or not from NEMOTAM (should always be true)
  - `nn_ittrjfrq_tan` the frequency (in timesteps) of NEMOTAM output. A snapshot of the passive tracer fields will be provided this often.
  - `cn_tantrj` the directory / filename structure of NEMOTAM output. Note that this will be ignored in passive mode; all outputs follow the structure `PTTAM_output_????????` in the experiment directory.
  - `ln_tam_out_{t,s,u,v,ssh,hdiv,rot,rhd,rhop}` these are the active variables in NEMOTAM. If NEMOTAM is **not** run in passive mode, these switches determine whether these variables will be present in the NEMOTAM output.
  - `ln_tam_in_{t,s,u,v,ssh,hdiv,rot,rhd,rhop}` these are again the active variables in NEMOTAM. If NEMOTAM is **not** run in passive mode, it will look for these variables in the NEMOTAM input file.
  - `ln_trj_out_{t,s,u,v,hdiv,rot,rhd,rhop}` these switches determine whether the **background** values of these variables (i.e. those in the trajectory, such as the underlying temperature distribution) will be returned alongside the passive tracer output in the NEMOTAM output.
- **nampttam** (passive-tracer-only options)
  - `ln_pt_regional` To save storage space when writing a trajectory. If set to true, trajectory files are only read within the domain specified by the following parameters: 
  - `rn_{NE,SW}pt{lat,lon}` the "corners" (northeast, southwest) of the quadrilateral region read by NEMOTAM. (Outside of this region, the trajectory output at timestep 1 is read by NEMOTAM repeatedly for the entire run)

- **namtra\_adv\_tam** (Options for tracer advection scheme)
  - `ln_traadv_cen2` 2nd order centred scheme
  - `ln_traadv_tvd` total variance diminishing scheme (*)
  - `rn_traadv_weight_h` balance between upwind and centred scheme (lateral advection)
  - `rn_traadv_weight_v` balance between upwind and centred scheme (vertical advection)

(*) If `ln_traadv_tvd = .true.` then the advection scheme is the nonlinear model default. In this case, the adjoint model is no longer a true adjoint of the tangent-linear model.

For `rn_traadv_weight_{h,v}`, a purely centred scheme is determined by value `0` and a purely upwind scheme by value `1`. If negative, then an automatic parameterisation decided by the _weighted mean_ scheme of **Fiadeiro & Veronis** (1977) is implemented.

### TAM output
If run in parallel, the model delegates parts of the grid to different CPUs. The output is thus returned in "tiles" as `PTTAM_output_<TIMESTEP>_????.nc`, which can be pieced together using `TOOLS/rebuild_nemo PTTAM_<TIMESTEP>_output <NO. OF TILES>`. Where the `TOOLS` directory is found at the same level as `CONFIG` (`../../TOOLS` from the install location of this repository). Each file contains the following variables:

- `nav_lon` and `nav_lat` (2D longitude and latitude arrays for the grid)
- `pt_conc_tl` (tangent linear) or `pt_vol_ad` (adjoint): a 4D (x,y,z,t) array describing the spatiotemporal distribution of passive tracer concentrations (tangent-linear) or volumes (adjoint).


**NOTE**: the tracer _concentration_ distribution is provided by the tangent-linear for mathematical consistency. The volume can be calculated by producing a `mesh_mask` file (`nn_msh = 1` in the namelist file) and multiplying the concentrations by `e1t`,`e2t` and `e3t` (the grid dimensions). The adjoint outputs volume automatically.

- `pt_vent_tl` (tangent-linear) or `pt_vent_ad` (adjoint) a 3D (x,y,t) array describing the **cumulative** removal of tracer at the surface during the run.
- any requested variables of the ocean background state at the point in the trajectory corresponding to the output (chosen using `ln_trj_out_` switches in the TAM namelist)

## Modifications to model defaults

Modifications are marked in the model source code using three exclamation marks, the date they were incorporated, and a letter (e.g. `!!!20200211A`) the end of the modification is marked by the same tag with a slash before the date (e.g. `!!!/20200211A`)

A list of modifications and details follows (All changes added by S. Mueller or myself. Asterisks indicate modifications which are essential to run NEMOTAM):

- 20191004B : single adjoint output variables (e.g. tn+tb)
- 20191004C : Add "sshn" to nonlinear trajectory
- 20191004D: Expanding trajectory and TAM filenames to allow at least 8-digit time-steps
- 20191004E*: correct transition of b->n for kdir==1 [see here](http://forge.ipsl.jussieu.fr/nemo/attachment/ticket/1443/trj_tam.F90.diff)
- 20191004F: corrected 'stpr2 - zstp -1' to 'stpr2 - zstp +1' in trj_tam.F90
- 20191004G: adjust output of interpolation coefficients 
- 20191004H: switch to allow adjoint output writing 
- 20191004I: essential modifications to dynzdf_imp_tam.F90 [see here](http://forge.ipsl.jussieu.fr/nemo/attachment/ticket/1362/dynzdf_imp_tam.F90.diff)
- 20191004J*: addition of adjoint time-stepping loop 
- 20191004K*: proper initialisation of TAM variables 
- 201910004L: ability to read cost function/perturbation from netCDF file
- 20191004N*: building of {t,u,v,f}msk_i variables using dom_uniq 
- 20191004O*: force flux SBC (http://forge.ipsl.jussieu.fr/nemo/attachment/ticket/1738/sbcmod_tam.F90.diff)[20160524]
- 20191004P: Passive tracer module + subroutines
- 20191004Q: Modifications to output ventilation record of passive tracer
- 20191004R: trajectory offsetting
- 20191004S: Introduction of weighted-mean avection scheme 
- 20191004T: Introduction of trajectory-upstream scheme and manual weightings
- 20191004U: Introduction of TVD scheme to tangent-linear and reversed fields to adjoint
- 20191004W: Option to include EIV
- 20191013A: Switch and parameters to only write and read trajectory in a lat/lon defined region. If a CPU has no points in this region, it doesn't write trajectory tiles, except on the first time-step of the run. When reading, this CPU simply reads the first time step repeatedly. Any passive tracer concentration outside of this region is set to 0 to prevent instabilities. Dramatically reduces storage space required to run the trajectory for localised passive tracer studies.
- 20200617A: Add switches to output trajectory variables in TAM output files
- 20200622A: Allow NEMO and NEMOTAM to run with differing timesteps
- 20200623A: Add switches to determine which (active) TAM variables are read
- 20200630A: (for ORCA025) bugfix to avoid negative runoff [by J.M. Molines](http://forge.ipsl.jussieu.fr/nemo/browser/trunk/NEMOGCM/NEMO/OPA_SRC/TRA/trazdf.F90?rev=5385)



## References
Fiadeiro, M.E. and Veronis, G., 1977. On weighted-mean schemes for the finite-difference approximation to the advection-diffusion equation. Tellus, 29(6), pp.512-522.

Large, W.G. and Yeager, S.G., 2004. Diurnal to decadal global forcing for ocean and sea-ice models: The data sets and flux climatologies.

Stephenson, D., Müller, S.A. and Sévellec, F., 2020. Tracking water masses using passive-tracer transport in NEMO v3. 4 with NEMOTAM: application to North Atlantic Deep Water and North Atlantic Subtropical Mode Water. Geoscientific Model Development, 13(4), pp.2031-2050.

Vidard, A., Vigilant, F., Benshila, R. and Deltel, C., 2011. NEMO Tangent & Adjoint Models (NemoTam) Reference Manual & User's Guide (Doctoral dissertation, INRIA).


