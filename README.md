# PT\_TAM\_ORCA2
## Overview
This is a NEMO *configuration* which can be used to track passive tracers in NEMO ORCA2 with NEMOTAM.  This repository requires existing NEMO v3.4 and NEMOTAM installs, and should sit in the `CONFIG` folder, e.g.

`/home/username/NEMO/dev_v3_4_STABLE_2012/NEMOGCM/CONFIG`

## Basics
NEMOTAM is a tangent-linear and adjoint model (TAM) counterpart to the nonlinear NEMO primitive equation solving ocean circulation model (c.f. *Vidard et al.*, 2011). As such, it details the evolution of small perturbations to a known background state (the "trajectory") obtained by running the nonlinear model (tangent-linear mode). In adjoint mode, "cost functions" of the ocean state are run backwards along the trajectory to determine their sensitivity to earlier perturbations.

In `PT_TAM_ORCA2`, passive tracer is injected into the model as a "perturbation" (/"cost function") , and tracked forward (/backward) using the tangent-linear (/adjoint) model. At the surface, the tracer is removed using a restoring mechanism. A record of tracer removal is returned along with a record of tracer distribution as part of NEMOTAM's output.

## Directory structure
The default structure of a NEMO configuration is 


- `/home/username/NEMO/dev_v3_4_STABLE_2012/NEMOGCM/CONFIG/`

 - `CONFIG_NAME/`
    - `BLD/`
    - `MY_SRC/`
    - `WORK/`
    - `EXP00/`


This repository contains additional directories, relevant to running passive tracer experiments, which serve as templates. The contents of these can either be copied or linked to from wherever you intend to run your experiments, as long as symbolic links to files in `../../NEMO` remain intact. The directories are:

- `PT_TAM_ORCA2/`
  - `EXP_DIR_TEMPLATE/`
    - `RUN_DIR/`
    - `ORCA2_INPUT/`

`RUN_DIR` can be renamed as desired, and is where the executables and namelists for a run are contained.

`ORCA2_INPUT` contains files necessary for running NEMO in this configuration. The file `README_3.4` describes the meaning of each file. COREII normal-year forcing files are not included in this repository as they are too large. Details on the missing files and how to acquire them are found in `MISSING_FILES`.

## Running passive tracer experiments


### Spinning up
To set up, the model should first be spun up from rest (ideally for around a thousand years) to provide a start point for all experiments.
Inside `RUN_DIR` is a bash script `SPINUP_scriptmaker.sh`, which produces namelist files (and job submission scripts for the *mobilis* HPC system) for consecutive runs of 50 years (273750 model time steps) each. They can be tuned by editing `namelist.SPINUP_template` (and `SPINUP_template.sh` for *mobilis*).

The end result should be the files `SPINUP_????????_restart.nc` and `SPINUP_????????_restart_ice.nc`, used as start points for the trajectory.

### Running a trajectory
After spinning up, the non-linear model should be run for the desired length of the experiment. Parameters can be set in the file `namelist.TRAJ`.

**At the beginning of the namelist file:**

Be sure to set `ln_rstart = .true.`, and set `nn_it000`, `nn_date0` and `cn_ocerst_in` in accordance with the `SPINUP_*.nc` files generated above.

**At the end of the namelist file**:

Be sure TAM trajectory output is produced with `ln_trjhand = .true.`. The trajectory output location and filename prefix is set by `cn_dirtrj`. All other output files are not required by `NEMOTAM`. 

**Other parameters**:
Other options are detailled in [the NEMO 3.4 manual](http://forge.ipsl.jussieu.fr/little_nemo/export/44/vendor/nemo/current/DOC/NEMO_book.pdf).

### Running the TAM:
**Initial tracer distribution**:

The file `PT_init.nc` is a plain text document describing the structure of the input NetCDF file for the model.

In the tangent-linear model, concentrations in `PT_init.nc` are propagated forward. In the adjoint model, the concentrations in `PT_init.nc` are automatically multiplied by the local domain volume to produce the volumetric cost function when the model is run backward. 

**At the beginning of the namelist file**:

Ensure `nn_it000 = 1` (even if beginning later in the trajectory). 
For the tangent-linear model, `nn_itend` should be set to the desired number of time steps in the run (up to the length of the trajectory). In the adjoint model, `nn_itend` determines the "start point" relative to the end of the trajectory, from which model will run backwards.

**At the end of the namelist file:**

- **namtrj**
 - `cn_dirtrj` should match that set in the trajectory namelist, specifying the location of the trajectory output.
 - `nn_ittrjoffset`determines the start (end) point of tangent-linear (adjoint) runs, if it is desired to begin the run at a later point in the trajectory.
- **namtst\_tam**
 - `ln_swi_opatam` determines TAM mode. `200` for tangent-linear and `201` for adjoint.
- **namtl_trj** These settings are not seen by PT\_TAM
- **nampttam**
 - `cn_pttam_init`Initial tracer distribution file
 - `nn_pttam_out_freq`Frequency of output (15 time steps = 1d)
- **namtra\_adv\_tam** Options for passive tracer advection scheme
 - `ln_traadv_cen2` 2nd order centred
 - `ln_traadv_tvd` total variance diminishing (**NONLINEAR**)
 - `rn_traadv_weight_h` balance between upwind and centred scheme (lateral advection)
 - `rn_traadv_weight_v` balance between upwind and centred scheme (centred advection)

If `ln_traadv_tvd = .true.` then the advection scheme is the nonlinear model default. In this case, the adjoint model is no longer a true adjoint of the tangent-linear model.

For `rn_traadv_weight_?`, a purely centred scheme is determined by value `0` and a purely upwind scheme by value `1`. If negative, then an automatic parameterisation decided by the _weighted mean_ scheme of Fiadeiro and Veronis (1977) is implemented.

### TAM output
_not yet complete_

## Modifications to model defaults
### _not yet complete_

## References
Fiadeiro, M.E. and Veronis, G., 1977. On weighted-mean schemes for the finite-difference approximation to the advection-diffusion equation. Tellus, 29(6), pp.512-522.

Vidard, A., Vigilant, F., Benshila, R. and Deltel, C., 2011. NEMO Tangent & Adjoint Models (NemoTam) Reference Manual & User's Guide (Doctoral dissertation, INRIA).