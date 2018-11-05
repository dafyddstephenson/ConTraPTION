MODULE sbc_ice
   !!======================================================================
   !!                 ***  MODULE  sbc_ice  ***
   !! Surface module - LIM-3: parameters & variables defined in memory
   !!======================================================================
   !! History :  3.0  ! 2006-08  (G. Madec)  Surface module
   !!            3.2  ! 2009-06  (S. Masson) merge with ice_oce
   !!            3.3.1! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !!            3.4  ! 2011-11  (C. Harris) CICE added as an option
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_lim2' or 'key_lim3' :             LIM-2 or LIM-3 sea-ice model
   !!----------------------------------------------------------------------
   USE par_oce          ! ocean parameters
   USE par_ice_2        ! LIM-2 parameters
   USE lib_mpp          ! MPP library
   USE in_out_manager   ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_ice_alloc ! called in iceini(_2).F90

   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim2    = .TRUE.   !: LIM-2 ice model
   LOGICAL         , PUBLIC, PARAMETER ::   lk_lim3    = .FALSE.  !: no LIM-3
   LOGICAL         , PUBLIC, PARAMETER ::   lk_cice    = .FALSE.  !: no CICE 
   CHARACTER(len=1), PUBLIC, PARAMETER ::   cp_ice_msh = 'C'      !: EVP: 'C'-grid ice-velocity

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qns_ice   !: non solar heat flux over ice                  [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qsr_ice   !: solar heat flux over ice                      [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   qla_ice   !: latent flux over ice                          [W/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dqla_ice  !: latent sensibility over ice                 [W/m2/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   dqns_ice  !: non solar heat flux over ice (LW+SEN+LA)    [W/m2/K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   tn_ice    !: ice surface temperature                          [K]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   alb_ice   !: albedo of ice

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   utau_ice  !: atmos-ice u-stress. VP: I-pt ; EVP: U,V-pts   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   vtau_ice  !: atmos-ice v-stress. VP: I-pt ; EVP: U,V-pts   [N/m2]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr1_i0    !: 1st Qsr fraction penetrating inside ice cover    [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   fr2_i0    !: 2nd Qsr fraction penetrating inside ice cover    [-]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   emp_ice   !: sublimation-snow budget over ice             [kg/m2]



   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: sbc_ice.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_ice_alloc()
      !!----------------------------------------------------------------------
      !!                     ***  FUNCTION sbc_ice_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( qns_ice (jpi,jpj,jpl) , qsr_ice (jpi,jpj,jpl) ,     &
         &      qla_ice (jpi,jpj,jpl) , dqla_ice(jpi,jpj,jpl) ,     &
         &      dqns_ice(jpi,jpj,jpl) , tn_ice  (jpi,jpj,jpl) ,     &
         &      alb_ice (jpi,jpj,jpl) ,                             &
         &      utau_ice(jpi,jpj)     , vtau_ice(jpi,jpj)     ,     &
         &      fr1_i0  (jpi,jpj)     , fr2_i0  (jpi,jpj)     ,     &
         &      emp_ice(jpi,jpj)                              , STAT= sbc_ice_alloc )
         !
      IF( lk_mpp            )   CALL mpp_sum ( sbc_ice_alloc )
      IF( sbc_ice_alloc > 0 )   CALL ctl_warn('sbc_ice_alloc: allocation of arrays failed')
   END FUNCTION sbc_ice_alloc


   !!======================================================================
END MODULE sbc_ice
