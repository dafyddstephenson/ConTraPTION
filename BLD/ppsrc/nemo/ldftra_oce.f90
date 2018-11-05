MODULE ldftra_oce
   !!=====================================================================
   !!                      ***  MODULE  ldftra_oce  ***
   !! Ocean physics :  lateral tracer mixing coefficient defined in memory 
   !!=====================================================================
   !! History :  9.0  !  2002-11  (G. Madec)  Original code
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp         ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC ldftra_oce_alloc ! called by nemo_init->nemo_alloc, nemogcm.F90

   !!----------------------------------------------------------------------
   !! Lateral eddy diffusivity coefficients (tracers)
   !!----------------------------------------------------------------------
   !                                                !!* Namelist namtra_ldf : lateral mixing *
   LOGICAL , PUBLIC ::   ln_traldf_lap   = .TRUE.    !: laplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_bilap = .FALSE.   !: bilaplacian operator
   LOGICAL , PUBLIC ::   ln_traldf_level = .FALSE.   !: iso-level direction
   LOGICAL , PUBLIC ::   ln_traldf_hor   = .FALSE.   !: horizontal (geopotential) direction
   LOGICAL , PUBLIC ::   ln_traldf_iso   = .TRUE.    !: iso-neutral direction
   LOGICAL , PUBLIC ::   ln_traldf_grif  = .FALSE.   !: griffies skew flux
   LOGICAL , PUBLIC ::   ln_traldf_gdia  = .FALSE.   !: griffies skew flux streamfunction diagnostics
   REAL(wp), PUBLIC ::   rn_aht_0        = 2000._wp  !: lateral eddy diffusivity (m2/s)
   REAL(wp), PUBLIC ::   rn_ahtb_0       =    0._wp  !: lateral background eddy diffusivity (m2/s)
   REAL(wp), PUBLIC ::   rn_aeiv_0       = 2000._wp  !: eddy induced velocity coefficient (m2/s)
   REAL(wp), PUBLIC ::   rn_slpmax       = 0.01_wp   !: slope limit

   REAL(wp), PUBLIC ::   aht0, ahtb0, aeiv0         !!: OLD namelist names
   LOGICAL , PUBLIC ::   ln_triad_iso    = .FALSE.   !: calculate triads twice
   LOGICAL , PUBLIC ::   ln_botmix_grif  = .FALSE.   !: mixing on bottom
   LOGICAL , PUBLIC ::   l_grad_zps      = .FALSE.   !: special treatment for Horz Tgradients w partial steps 

   REAL(wp), PUBLIC ::   rldf                        !: multiplicative factor of diffusive coefficient
                                                     !: Needed to define the ratio between passive and active tracer diffusion coef. 

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   ahtt, ahtu, ahtv, ahtw   !: ** 2D coefficients ** at T-,U-,V-,W-points

   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'                              eddy induced velocity
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER               ::   lk_traldf_eiv   = .TRUE.   !: eddy induced velocity flag
   
   !                                                                              !!! eddy coefficients at U-, V-, W-points  [m2/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   aeiu , aeiv , aeiw   !: ** 2D coefficients **
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   u_eiv, v_eiv, w_eiv   !: eddy induced velocity [m/s]


   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_oce.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION ldftra_oce_alloc()
     !!----------------------------------------------------------------------
      !!                 ***  FUNCTION ldftra_oce_alloc  ***
     !!----------------------------------------------------------------------
     INTEGER, DIMENSION(3) :: ierr
     !!----------------------------------------------------------------------
     ierr(:) = 0

      ALLOCATE( ahtt(jpi,jpj    ) , ahtu(jpi,jpj    ) , ahtv(jpi,jpj    ) , ahtw(jpi,jpj    ) , STAT=ierr(1) )
      !
      ALLOCATE( aeiu(jpi,jpj    ) , aeiv(jpi,jpj    ) , aeiw(jpi,jpj    ) , STAT=ierr(2) )
      ALLOCATE( u_eiv(jpi,jpj,jpk), v_eiv(jpi,jpj,jpk), w_eiv(jpi,jpj,jpk), STAT=ierr(3))
      ldftra_oce_alloc = MAXVAL( ierr )
      IF( ldftra_oce_alloc /= 0 )   CALL ctl_warn('ldftra_oce_alloc: failed to allocate arrays')
      !
   END FUNCTION ldftra_oce_alloc

   !!=====================================================================
END MODULE ldftra_oce
