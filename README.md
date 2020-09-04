# ConTraPTION: Configuration for Tracking Passive Tracers Injected into the Ocean with NEMOTAM
## Overview
**2020-06-23: README UPDATE DUE**
This is a NEMO *configuration* which can be used to track passive tracers in NEMO ORCA2 with NEMOTAM.  

NEMOTAM is a tangent-linear and adjoint model (TAM) counterpart to the nonlinear NEMO primitive equation solving ocean circulation model (c.f. *Vidard et al.*, 2011). As such, it details the evolution of small perturbations to a known background state (the "trajectory") obtained by running the nonlinear model (tangent-linear mode). In adjoint mode, "cost functions" of the ocean state are run backwards along the trajectory to determine their sensitivity to earlier perturbations.

In `ConTraPTION`, passive tracer is injected into the model as a "perturbation" (/"cost function") , and tracked forward (/backward) using the tangent-linear (/adjoint) model. At the surface, the tracer can be removed using a restoring mechanism. A record of tracer removal is returned along with a record of tracer distribution as part of NEMOTAM's output.

## Installation
### Installing NEMO and ConTraPTION
This repository requires existing NEMO v3.4 (and NEMOTAM) installs. These can be installed using:

`svn co https://forge.ipsl.jussieu.fr/nemo/svn/NEMO/releases/release-3.4`

This provides the source code for the model and reference configurations, on which the PT_TAM configuration should be based.

For a first-time NEMO install, the user will need to set up their machine's architecture. This will either already exist somewhere inside `NEMOGCM/ARCH` (e.g. `ARCH/<DIR>/arch-<ARCHITECTURE>.fcm`) or will have to be created. Instructions for creating/modifying this file are found [here](http://forge.ipsl.jussieu.fr/nemo/wiki/Users/ModelInstall#Setupyourarchitectureconfigurationfile).

At this point, the PT_TAM configuration can be created. To do this, run the following command in the `CONFIG` directory:

`./makenemo -d "OPATAM_SRC LIM_SRC_2 OPA_SRC" -n ConTraPTION -m <ARCHITECTURE> add_key "key_mpp_mpi  key_mpp_rep key_nosignedzero key_tam key_diainstant" del_key "key_zdfddm key_iomput"`

This creates a new configuration with ocean (`OPA_SRC`) and sea-ice (`LIM_SRC_2`) dynamics, and TAM (`OPATAM_SRC`, `key_tam`) compatibility. 
(NOTE: customise line as required, e.g. regarding parallel processing - `key_mpp_mpi` and `key_mpp_rep`. Replace <ARCHITECTURE> with the relevant substring from the arch. filename)
 
This configuration can now be modified to include passive tracer-related subroutines. To do so, move the contents of `MY_SRC` from this repository into the `NEMOGCM/CONFIG/ConTraPTION/MY_SRC` directory and recompile using

`./makenemo -n ConTraPTION`

### Obtaining and linking forcing and other model input files
In order to run NEMO-ORCA2, additional files are required, which can be found [here](https://doi.org/10.5281/zenodo.1471702).

Unpack this archive into a directory named `ORCA2_INPUT`. This will be at the same level as your experiment directory and can be anywhere you choose.

### Creating an experiment directory
At the same level as `ORCA2_INPUT`, create an experiment directory (e.g. `RUN_DIR`). Copy the contents of our template experiment directory into this directory. The bash script `link_runfiles.sh` creates softlinks to all files necessary to run NEMO and NEMOTAM. Edit it to provide the location of your `ORCA2_INPUT` and `BLD/bin`.

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
To set up, the model should first be spun up from rest (ideally for around a thousand years) to provide a start point for all experiments.
Inside `RUN_DIR/SPINUP_RUNS` is a bash script `make_spinup_namelists.sh`, which produces namelist files for consecutive runs of 50 years (273750 model time steps) each. They can be tuned by editing this script along with `namelist.SPINUP_template`.

The end result should be the files `SPINUP_????????_restart.nc` and `SPINUP_????????_restart_ice.nc`, used as start points for the trajectory.

### Running a trajectory
After spinning up, the non-linear model should be run for the desired length of the experiment. A script is provided (`RUN_DIR/TRAJECTORY_RUNS/make_trajectory_namelists.sh`) to produce ocean and ice namelists for consecutive 50 year trajectories adding up to 400 years. The namelist parameters themselves can be adjusted in `RUN_DIR/TRAJECTORY_RUNS/namelist.TRAJ_TEMPLATE` and `RUN_DIR/TRAJECTORY_RUNS/namelist_ice.TRAJ_TEMPLATE`.

**At the beginning of the namelist file:**

Be sure to set `ln_rstart = .true.`, and set `nn_it000`, `nn_date0` and `cn_ocerst_in` in accordance with the `SPINUP_*.nc` files generated above.

**At the end of the namelist file**:

Be sure TAM trajectory output is produced with `ln_trjhand = .true.`. The trajectory output location and filename prefix is set by `cn_dirtrj`. All other output files are not required by `NEMOTAM`. 

**Other parameters**:
Other options are detailled in [the NEMO 3.4 manual](http://forge.ipsl.jussieu.fr/little_nemo/export/44/vendor/nemo/current/DOC/NEMO_book.pdf).

### Running the TAM:
**Initial tracer distribution**:

The provided file `PT_init.py` is a Python script for generating initial passive tracer distributions and saving them to the netCDF file `PT_init.nc`, which the template namelist points to.

In the tangent-linear model, concentrations in `PT_init.nc` are propagated forward. In the adjoint model, the concentrations in `PT_init.nc` are automatically multiplied by the local domain volume to produce the volumetric cost function when the model is run backward. 

**At the beginning of the namelist file**:

Ensure `nn_it000 = 1` (even if beginning later in the trajectory). 
For the tangent-linear model, `nn_itend` should be set to the desired number of time steps in the run (up to the length of the trajectory). In the adjoint model, `nn_itend` determines the "start point" relative to the end of the trajectory, from which model will run backwards.

**Surface ventilation**:

- **namsbc\_ssr**
  - The namelist parameter `nn_sstr` determines whether the passive tracer is removed at the surface
  - The namelist parameter `rn_dqdt` sets the time scale of the surface restoring scheme. The default, -40Wm^-2K^-1, corresponds to 60 days over a 50 m mixed layer.

**At the end of the namelist file:**

- **namtrj**
  - `cn_dirtrj` should match that set in the trajectory namelist, specifying the location of the trajectory output.
  - `nn_ittrjoffset` determines the start (end) point of tangent-linear (adjoint) runs, if it is desired to begin the run at a later point in the trajectory.
- **namtst\_tam**
  - `ln_swi_opatam` determines TAM mode. `200` for tangent-linear and `201` for adjoint.
- **namtl_trj** These settings are not seen by PT\_TAM
- **nampttam**
  - `cn_pttam_init` is the initial tracer distribution file
  - `nn_pttam_out_freq` is the desired frequency of output (15 time steps = 1d)
- **namtra\_adv\_tam** Options for passive tracer advection scheme
  - `ln_tl_eiv` include eddy-induced velocity (T) or not (F)
  - `ln_traadv_cen2` 2nd order centred scheme
  - `ln_traadv_tvd` total variance diminishing scheme (**NOTE: NONLINEAR**)
  - `rn_traadv_weight_h` balance between upwind and centred scheme (lateral advection)
  - `rn_traadv_weight_v` balance between upwind and centred scheme (centred advection)

If `ln_traadv_tvd = .true.` then the advection scheme is the nonlinear model default. In this case, the adjoint model is no longer a true adjoint of the tangent-linear model.

For `rn_traadv_weight_{h,v}`, a purely centred scheme is determined by value `0` and a purely upwind scheme by value `1`. If negative, then an automatic parameterisation decided by the _weighted mean_ scheme of Fiadeiro and Veronis (1977) is implemented.

### TAM output
If run in parallel, the model delegates parts of the grid to different CPUs. The output is thus returned in "tiles" as `PTTAM_output_<TIMESTEP>_????.nc`, which can be pieced together using `TOOLS/rebuild_nemo PTTAM_<TIMESTEP>_output <NO. OF TILES>`. Where the `TOOLS` directory is found at the same level as `CONFIG` (`../../TOOLS` from the install location of this repository). Each file contains the following variables:

- `nav_lon` and `nav_lat` (2D longitude and latitude arrays for the grid)
- `pt_conc_tl` (tangent linear) or `pt_vol_ad` (adjoint): a 4D (x,y,z,t) array describing the spatiotemporal distribution of passive tracer concentrations (tangent-linear) or volumes (adjoint).


**NOTE**: the tracer _concentration_ distribution is provided by the tangent-linear for mathematical consistency. The volume can be calculated by producing a `mesh_mask` file (`nn_msh = 1` in the namelist file) and multiplying the concentrations by `e1t`,`e2t` and `e3t` (the grid dimensions). The adjoint outputs volume automatically.

- `pt_vent_tl` (tangent-linear) or `pt_vent_ad` (adjoint) a 3D (x,y,t) array describing the _cumulative_ removal of tracer at the surface during the run (needs to be differenced to produce a time series).
- `tn` and `sn`: the temperature and salinity of the ocean background state at the point in the trajectory corresponding to the output. May be used, for example, to produce TS plots of tracer transformation.

## Modifications to model defaults

There are three types of modifications present: bug fixes to model defaults (e.g. patching undeclared variables in `dynzdf_imp_tam.F90`), general updates for convenience (e.g. extending maximum run length) and modifications specific to passive tracer transport (e.g. bespoke passive tracer advection scheme)

Passive tracer mode (`ln_swi_opatam` >200) is kept separate from the model's standard tangent-linear (`ln_swi_opatam = 2`) and adjoint (`ln_swi_opatam = 3`) modes, using a different set of routines, defined in `MY_SRC/pt_tam.f90`.

### NEMOTAM bug fixes:
- `domwri.F90`: argument `kindic` added to subroutine `dom_uniq`. `kindic` is used by the `oce_tam_init` routine to initialise TAM variables (see following point)
- `oce_tam.F90` building of `{t,u,v,f}msk_i` variables using `dom_uniq` (see [here](http://forge.ipsl.jussieu.fr/nemo/ticket/1499) ).
- `dynzdf_imp_tam.F90` (see [here](http://forge.ipsl.jussieu.fr/nemo/attachment/ticket/1362/dynzdf_imp_tam.F90.diff) ).
- `trj_tam.F90`: 
  - Addition of routines `ad_trj_ini`, `ad_trj_wri` and variables for writing the adjoint trajectory 
  - Corrected trajectory time step reading in adjoint (see [here](http://forge.ipsl.jussieu.fr/nemo/attachment/ticket/1443/trj_tam.F90.diff) ).
- `sbcmod_tam.F90`: removal of bug that forces default (GYRE configuration) surface boundary conditions, enabling flux boundary conditions (see [here](http://forge.ipsl.jussieu.fr/nemo/attachment/ticket/1738/sbcmod_tam.F90.diff) ).
  
### General updates to NEMOTAM
- `tamtrj.F90`
  - `nn_ittrjoffset` parameter added, allowing TL runs to start later than the beginning of the trajectory
  - `cl_dirtrj` variable, which sets output filenames and location is modified to support up to 100 million time steps (18ky) with a consistent filename structure (`PT_TAM_output_????????.nc`), or higher with modified filename structure.
 - `trj_tam.F90`
  - `cl_dirtrj` variable modified as in `tamtrj.F90` and similarly `cl_tantrj` and `cl_adjtrj`
  - added SSH (`sshn`) to trajectory outputs
  - Addition of `lreset` argument to `trj_rea` which "forgets" about previous internal state. This is used if moving to a point in the trajectory not anticipated by the model (e.g. a jump)
- `traadv_tam.F90`
  - `ln_tl_eiv` parameter added, allowing eddy-induced velocity for tracer advection

### Modifications for passive tracer model:

- `sbcmod_tam.F90` handles surface boundary conditions. For `PT_TAM`
  - A new variable `sbc_tmp_rm`is defined **in the adjoint**, which holds the cumulative tracer removed at the surface (summed from the instantaneous variable `sst_m_ad` defined in the routine `sbc_ssr_adj`). In the tangent-linear, the surface removal scheme is handled by routine `pt_tan` in `pt_tam.F90`.
- `nemogcm_tam.F90`
  - Added additional cases for `ln_swi_opatam` variable: values between 200 and 249 signal passive tracer mode, and the routines of `pt_tam.F90` take over.
- `traadv_tam.F90`
  - Commented options for all nonlinear advection schemes except TVD
  - If any scheme other than TVD is selected, the 2nd order linear scheme is forced 
  - For `ln_traadv_cen2`, variables `rn_traadv_weight_h` and `rn_traadv_weight_v` are defined. Their values are written to the `ocean.output` file during the run. Routines `tra_adv_cen2_tan` or `tra_adv_cen2_tan` are called with their values.
- `traadv_cen2_tam.F90`
  - horizontal and vertical diffusion coefficients (from `ldftra_oce` and  `zdf_oce` are now seen by subroutines)
  - added variables `{p,z}r_ups_{h,v}` to determine weighting between centred and upstream schemes
  - Defined `zsgm{u,v,w}` and `ztht{u,v,w}`, the sigma and theta parameters of the Fiadeiro and Veronis (1977, p6) weighted mean advection scheme 
  - Added calculation of theta in u and v
  - Added calculation of sigma from theta. Sigma becomes the weighting parameter between centred and upwind.
  - Added similar calculations for vertical advection
  - Added similar calculations in the adjoint routine
  - Implemented advection scheme choices for adjoint test routine
- **`pt_tam.F90` (new file)**, contains subroutines:
  - `pt_init`
    - reads `cn_pttam_init` from namelist to determine input file
    - reads `nn_pttam_out_freq` from namelist to determine frequency of TAM output
    - allocates variables `tmsk_region`, `tmsk_nasmw`, `tmp_rm` and `sal_rm`
  - `pt_finalise`
    - deallocates variables  `tmsk_region`, `tmsk_nasmw`, `tmp_rm` and `sal_rm` 
  - `pt_run`
    - sees `ln_swi_opatam` and calls either `pt_tan` (200) or `pt_adj` (201) 
  - `pt_tan`
    - called in place of `stp_tan` in `nemogcm_tam.F90`. `pt_tan` also calls `stp_tan` but performs additional steps to reset dynamic feedbacks between time steps (making tracers passive) and to calculate `tmp_rm` surface removal volumes.
  - `pt_adj`
    - called in place of `stp_adj` in `nemogcm_tam.F90`. `pt_tan` also calls `stp_adj` but performs additional steps to reset dynamic feedbacks between time steps (making tracers passive). `tmp_rm` is obtained from `sbc_tmp_rm` here, calculated in our modified`sbcmod_tam.F90`.
  - `pt_tam_wri`
    - This subroutine sets the filename structure and variables of the output files of a passive tracer run. The variable `wri_swi` is used in the call to pass on whether to write tangent-linear or adjoint variables.
## References
Fiadeiro, M.E. and Veronis, G., 1977. On weighted-mean schemes for the finite-difference approximation to the advection-diffusion equation. Tellus, 29(6), pp.512-522.

Vidard, A., Vigilant, F., Benshila, R. and Deltel, C., 2011. NEMO Tangent & Adjoint Models (NemoTam) Reference Manual & User's Guide (Doctoral dissertation, INRIA).
