MODULE trdmld_oce
   !!======================================================================
   !!                   ***  MODULE trdmld_oce  ***
   !! Ocean trends :   set tracer and momentum trend variables
   !!======================================================================
   !! History :  1.0  ! 2004-08  (C. Talandier)  New trends organization
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters

   IMPLICIT NONE
   PRIVATE

   PUBLIC   trdmld_oce_alloc    ! Called in trdmld.F90

   LOGICAL, PUBLIC, PARAMETER ::   lk_trdmld = .FALSE.   !: ML trend flag
   !!* mixed layer trends indices
   INTEGER, PARAMETER, PUBLIC ::   jpltrd = 11      !: number of mixed-layer trends arrays
   INTEGER, PUBLIC            ::   jpktrd           !: max level for mixed-layer trends diag.
   !
   INTEGER, PUBLIC, PARAMETER ::   jpmld_xad =  1   !:  zonal      
   INTEGER, PUBLIC, PARAMETER ::   jpmld_yad =  2   !:  meridonal   > advection
   INTEGER, PUBLIC, PARAMETER ::   jpmld_zad =  3   !:  vertical   
   INTEGER, PUBLIC, PARAMETER ::   jpmld_ldf =  4   !:  lateral diffusion (geopot. or iso-neutral)
   INTEGER, PUBLIC, PARAMETER ::   jpmld_for =  5   !:  forcing 
   INTEGER, PUBLIC, PARAMETER ::   jpmld_zdf =  6   !:  vertical diffusion (TKE)
   INTEGER, PUBLIC, PARAMETER ::   jpmld_bbc =  7   !:  geothermal flux
   INTEGER, PUBLIC, PARAMETER ::   jpmld_bbl =  8   !:  bottom boundary layer (advective/diffusive)
   INTEGER, PUBLIC, PARAMETER ::   jpmld_dmp =  9   !:  internal restoring trend
   INTEGER, PUBLIC, PARAMETER ::   jpmld_npc = 10   !:  non penetrative convective adjustment
!! INTEGER, PUBLIC, PARAMETER ::   jpmld_xxx = xx   !:  add here any additional trend (add change jpltrd)
   INTEGER, PUBLIC, PARAMETER ::   jpmld_atf = 11   !:  asselin trend (**MUST BE THE LAST ONE**)

   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: trdmld_oce.F90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

  INTEGER FUNCTION trdmld_oce_alloc()
     !!----------------------------------------------------------------------
     !!                 ***  FUNCTION trdmld_oce_alloc   ***
     !!----------------------------------------------------------------------
     USE lib_mpp
     INTEGER :: ierr(5)
     !!----------------------------------------------------------------------

     ! Initialise jpktrd here as can no longer do it in MODULE body since
     ! jpk is now a variable.
     jpktrd = jpk   !: max level for mixed-layer trends diag.

     ierr(:) = 0

      !
      trdmld_oce_alloc = MAXVAL( ierr )
      IF( lk_mpp                )   CALL mpp_sum ( trdmld_oce_alloc )
      IF( trdmld_oce_alloc /= 0 )   CALL ctl_warn('trdmld_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION trdmld_oce_alloc

   !!======================================================================
END MODULE trdmld_oce
