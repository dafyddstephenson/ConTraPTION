MODULE sbc_oce_tam
   !!----------------------------------------------------------------------
   !!    This software is governed by the CeCILL licence (Version 2)
   !!----------------------------------------------------------------------
   !!======================================================================
   !!                       ***  MODULE  sbc_oce_tam  ***
   !! Surface module :   variables defined in core memory
   !!                    Tangent and Adjoint Module
   !!======================================================================
   !! History of the direct module:
   !!            3.0   !  2006-06  (G. Madec)  Original code
   !!             -    !  2008-08  (G. Madec)  namsbc moved from sbcmod
   !! History of the TAM module:
   !!            3.0   !  2008-11  (A. Vidard)  Original code
   !!            3.0   !  2009-03  (A. Weaver)  Allocate/initialization routine
   !!            3.2   !  2010-04  (A. Vidard)  3.2 update
   !!            3.4   !  2012-07  (P.-A. Bouttier) 3.4 update
   !!             3.4  ! 2012-09   (A. Vidard) Deallocating and initialising options
   !!----------------------------------------------------------------------
   USE par_oce
   USE sbc_oce
   USE in_out_manager
   USE lib_mpp

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC &
      & sbc_oce_alloc_tam,   & !: Allocate teh TAM fields
      & sbc_oce_dealloc_tam, & !: Deallocate the TAM fields
      & sbc_oce_tam_init       !: Initialize the TAM fields

   !!----------------------------------------------------------------------
   !!              Ocean Surface Boundary Condition fields
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: &
      & qsr_hc_tl,  &
      & sbc_tsc_b_tl,  &
      & sbc_tsc_tl,  &
      & qsr_hc_b_tl
   REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE :: &
      & utau_tl, & !: Tangent linear of sea surface i-stress (ocean referential)     [N/m2]
      & vtau_tl, & !: Tangent linear of sea surface j-stress (ocean referential)     [N/m2]
      & utau_b_tl, & !: Tangent linear of sea surface i-stress (ocean referential)     [N/m2]
      & vtau_b_tl, & !: Tangent linear of sea surface j-stress (ocean referential)     [N/m2]
      & taum_tl, & !: Tangent linear of module of sea surface stress (at T-point)    [N/m2]
      & wndm_tl, & !: Tangent linear of wind speed module at T-point (=|U10m-Uoce|)  [m/s
      & qns_tl,  & !: Tangent linear of sea heat flux: non solar                     [W/m2]
      & qns_b_tl,  & !: Tangent linear of sea heat flux: non solar                     [W/m2]
      & rnf_tl,  & !: Tangent linear of sea heat flux: non solar                     [W/m2]
      & rnf_b_tl,  & !: Tangent linear of sea heat flux: non solar                     [W/m2]
      & qsr_tl,  & !: Tangent linear of sea heat flux:     solar                     [W/m2]
      & qns_tot_tl,  & !: Tangent linear of total sea heat flux: non solar                     [W/m2]
      & qsr_tot_tl,  & !: Tangent linear of total sea heat flux:     solar                     [W/m2]
      & emp_tl,  & !: Tangent linear of freshwater budget: volume flux               [Kg/m2/s]
      & emp_tot_tl,  & !: Tangent linear of freshwater budget: volume flux               [Kg/m2/s]
      & emps_tl, & !: Tangent linear of freshwater budget: concentration/dillution   [Kg/m2/s]
      & emp_b_tl,  & !: Tangent linear of freshwater budget: volume flux               [Kg/m2/s]
      & emps_b_tl, & !: Tangent linear of freshwater budget: concentration/dillution   [Kg/m2/s]
      & fr_i_tl    !: Tangent linear of ice fraction  (between 0 to 1)               -
   LOGICAL, SAVE, PRIVATE :: ll_alloctl = .FALSE.

   REAL(wp), PUBLIC, DIMENSION(:,:,:), ALLOCATABLE :: &
      & qsr_hc_ad,  &
      & sbc_tsc_b_ad,  &
      & sbc_tsc_ad,  &
      & qsr_hc_b_ad

   REAL(wp), PUBLIC, DIMENSION(:,:), ALLOCATABLE :: &
      & utau_ad, & !: Adjoint of sea surface i-stress (ocean referential)     [N/m2]
      & vtau_ad, & !: Adjoint of sea surface j-stress (ocean referential)     [N/m2]
      & utau_b_ad, & !: Adjoint of sea surface i-stress (ocean referential)     [N/m2]
      & vtau_b_ad, & !: Adjoint of sea surface j-stress (ocean referential)     [N/m2]
      & taum_ad, & !: Adjoint of module of sea surface stress (at T-point)    [N/m2]
      & wndm_ad, & !: Adjoint of wind speed module at T-point (=|U10m-Uoce|)  [m/s]
      & qns_ad,  & !: Adjoint of sea heat flux: non solar                     [W/m2
      & qns_b_ad,  & !: Adjoint of sea heat flux: non solar                     [W/m2
      & rnf_ad,  & !: Adjoint of sea heat flux: non solar                     [W/m2
      & rnf_b_ad,  & !: Adjoint of sea heat flux: non solar                     [W/m2
      & qsr_ad,  & !: Adjoint of sea heat flux:     solar                     [W/m2]
      & qns_tot_ad,  & !: Adjoint of total sea heat flux: non solar                     [W/m2
      & qsr_tot_ad,  & !: Adjoint of total sea heat flux:     solar                     [W/m2]
      & emp_ad,  & !: Adjoint of freshwater budget: volume flux               [Kg/m2/s]
      & emp_tot_ad,  & !: Adjoint of freshwater budget: volume flux               [Kg/m2/s]
      & emps_ad, & !: Adjoint of freshwater budget: concentration/dillution   [Kg/m2/s]
      & emp_b_ad,  & !: Adjoint of freshwater budget: volume flux               [Kg/m2/s]
      & emps_b_ad, & !: Adjoint of freshwater budget: concentration/dillution   [Kg/m2/s]
      & fr_i_ad    !: Adjoint of ice fraction  (between 0 to 1)               -

   !!----------------------------------------------------------------------
   !!                     Sea Surface Mean fields
   !!----------------------------------------------------------------------
   REAL(wp), PUBLIC,  DIMENSION(:,:), ALLOCATABLE :: &
      & ssu_m_tl, & !: mean (nn_fsbc time-step) surface sea i-current (U-point) [m/s]
      & ssv_m_tl, & !: mean (nn_fsbc time-step) surface sea j-current (V-point) [m/s]
      & sst_m_tl, & !: mean (nn_fsbc time-step) surface sea temperature     [Celsius]
      & sss_m_tl, & !: mean (nn_fsbc time-step) surface sea salinity            [psu]
      & ssh_m_tl    !: mean (nn_fsbc time-step) surface sea height                [m]

   REAL(wp), PUBLIC,  DIMENSION(:,:), ALLOCATABLE :: &
      & ssu_m_ad, & !: mean (nn_fsbc time-step) surface sea i-current (U-point) [m/s]
      & ssv_m_ad, & !: mean (nn_fsbc time-step) surface sea j-current (V-point) [m/s]
      & sst_m_ad, & !: mean (nn_fsbc time-step) surface sea temperature     [Celsius]
      & sss_m_ad, & !: mean (nn_fsbc time-step) surface sea salinity            [psu]
      & ssh_m_ad    !: mean (nn_fsbc time-step) surface sea height                [m]

   LOGICAL, SAVE, PRIVATE :: ll_allocad = .FALSE.
   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop starting and inding indices 
   !!      to allow unrolling of do-loop using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: vectopt_loop_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   INTEGER FUNCTION sbc_oce_alloc_tam( kmode )
      !! * Arguments
      INTEGER, OPTIONAL :: kmode
      INTEGER :: imode
      INTEGER :: ierr(4)
      IF ( PRESENT(kmode) ) THEN
         imode = kmode
      ELSE
         imode = 0
      END IF
      !------------------------------------------------------------------------
      ierr(:) = 0
      !
      IF ( ( imode == 0 ) .OR. ( imode == 1 ) .AND. ( .NOT. ll_alloctl ) ) THEN
      ALLOCATE( utau_tl(jpi,jpj) , utau_b_tl(jpi,jpj) , taum_tl(jpi,jpj) ,     &
         &      vtau_tl(jpi,jpj) , vtau_b_tl(jpi,jpj) , wndm_tl(jpi,jpj) , fr_i_tl(jpi,jpj), &
         &      STAT=ierr(1) )
         !
      ALLOCATE( qns_tot_tl(jpi,jpj) , qns_tl   (jpi,jpj) , qns_b_tl(jpi,jpj),     &
         &      qsr_tot_tl(jpi,jpj) , qsr_tl   (jpi,jpj) ,                        &
         &      emp_tl    (jpi,jpj) , emp_b_tl (jpi,jpj) ,                        &
         &      rnf_tl    (jpi,jpj) , rnf_b_tl (jpi,jpj) ,                        &
         &      emps_tl   (jpi,jpj) , emps_b_tl(jpi,jpj) , emp_tot_tl(jpi,jpj) , STAT=ierr(2) )
         !
      ALLOCATE(                                                                 &
         &      ssu_m_tl  (jpi,jpj) , sst_m_tl(jpi,jpj) ,                       &
         &      ssv_m_tl  (jpi,jpj) , sss_m_tl  (jpi,jpj), ssh_m_tl(jpi,jpj) , STAT=ierr(3) )
         !
      ALLOCATE( qsr_hc_tl   (jpi,jpj,jpk)  , qsr_hc_b_tl(jpi,jpj,jpk) ,         &
         &      sbc_tsc_b_tl(jpi,jpj,jpts) , sbc_tsc_tl(jpi,jpj,jpts), STAT=ierr(4) )

      sbc_oce_alloc_tam = MAXVAL( ierr )
      IF( lk_mpp            )   CALL mpp_sum ( sbc_oce_alloc_tam )
      IF( sbc_oce_alloc_tam > 0 )   CALL ctl_warn('sbc_oce_alloci_tam: allocation of tangent linear arrays failed')

         ll_alloctl = .TRUE.
         !
      END IF

      IF ( ( imode == 0 ) .OR. ( imode == 2 ) .AND. ( .NOT. ll_allocad ) ) THEN
         ierr(:) = 0
         !
         ALLOCATE( utau_ad(jpi,jpj) , utau_b_ad(jpi,jpj) , taum_ad(jpi,jpj) ,      &
            &      vtau_ad(jpi,jpj) , vtau_b_ad(jpi,jpj) , wndm_ad(jpi,jpj) , fr_i_ad(jpi,jpj),    &
            & STAT=ierr(1) )
         !
         ALLOCATE( qns_tot_ad(jpi,jpj) , qns_ad   (jpi,jpj) , qns_b_ad(jpi,jpj),   &
            &      qsr_tot_ad(jpi,jpj) , qsr_ad   (jpi,jpj) ,                      &
            &      emp_ad    (jpi,jpj) , emp_b_ad (jpi,jpj) ,                      &
            &      rnf_ad    (jpi,jpj) , rnf_b_ad (jpi,jpj) ,                      &
            &      emps_ad   (jpi,jpj) , emps_b_ad(jpi,jpj) , emp_tot_ad(jpi,jpj) , STAT=ierr(2) )
         !
         ALLOCATE(                                                                 &
            &      ssu_m_ad  (jpi,jpj) , sst_m_ad(jpi,jpj) ,                       &
            &      ssv_m_ad  (jpi,jpj) , sss_m_ad  (jpi,jpj), ssh_m_ad(jpi,jpj) , STAT=ierr(3) )
         !
         ALLOCATE( qsr_hc_ad(jpi,jpj,jpk), qsr_hc_b_ad(jpi,jpj,jpk),               &
            &      sbc_tsc_b_ad(jpi,jpj,jpts), sbc_tsc_ad(jpi,jpj,jpts), STAT=ierr(4) )

         sbc_oce_alloc_tam = MAXVAL( ierr )
         IF( lk_mpp            )   CALL mpp_sum ( sbc_oce_alloc_tam )
         IF( sbc_oce_alloc_tam > 0 )   CALL ctl_warn('sbc_oce_alloc_tam: allocation of adjoint arrays failed')
      !
         ll_allocad = .TRUE.
      END IF
   END FUNCTION sbc_oce_alloc_tam
   !
   INTEGER FUNCTION sbc_oce_dealloc_tam( kmode )
      !! * Arguments
      INTEGER, OPTIONAL :: kmode
      INTEGER :: ierr(4)
      INTEGER :: imode
      IF ( PRESENT(kmode) ) THEN
         imode = kmode
      ELSE
         imode = 0
      END IF
      !------------------------------------------------------------------------
      ierr(:) = 0
      !
      IF ( ( imode == 0 ) .OR. ( imode == 1 ) .AND. ( ll_alloctl ) ) THEN
         DEALLOCATE( utau_tl , utau_b_tl , taum_tl ,     &
            &        vtau_tl , vtau_b_tl , wndm_tl , STAT=ierr(1) )
         !
         DEALLOCATE( qns_tot_tl , qns_tl    , qns_b_tl,              &
            &        qsr_tot_tl , qsr_tl    ,                        &
            &        emp_tl     , emp_b_tl  ,                        &
            &        rnf_tl     , rnf_b_tl  ,                        &
            &        emps_tl    , emps_b_tl , emp_tot_tl , STAT=ierr(2) )
         !
         DEALLOCATE(                                                 &
            &        ssu_m_tl   , sst_m_tl ,                         &
            &        ssv_m_tl   , sss_m_tl  , ssh_m_tl , STAT=ierr(3) )
         !
         DEALLOCATE( qsr_hc_tl, qsr_hc_b_tl, sbc_tsc_b_tl , sbc_tsc_tl, STAT=ierr(4) )

         sbc_oce_dealloc_tam = MAXVAL( ierr )
         IF( lk_mpp                  )   CALL mpp_sum ( sbc_oce_dealloc_tam )
         IF( sbc_oce_dealloc_tam > 0 )   CALL ctl_warn('sbc_oce_dealloci_tam: allocation of tangent linear arrays failed')

         ll_alloctl = .FALSE.
         !
      END IF

      IF ( ( imode == 0 ) .OR. ( imode == 2 ) .and. ( ll_allocad ) ) THEN
         ierr(:) = 0
         !
         DEALLOCATE( utau_ad , utau_b_ad , taum_ad ,                 &
            &        vtau_ad , vtau_b_ad , wndm_ad , STAT=ierr(1) )
         !
         DEALLOCATE( qns_tot_ad , qns_ad    , qns_b_ad,              &
            &        qsr_tot_ad , qsr_ad    ,                        &
            &        emp_ad     , emp_b_ad  ,                        &
            &        rnf_ad     , rnf_b_ad  ,                        &
            &        emps_ad    , emps_b_ad , emp_tot_ad , STAT=ierr(2) )
         !
         DEALLOCATE(                                                 &
            &        ssu_m_ad   , sst_m_ad ,                         &
            &        ssv_m_ad   , sss_m_ad  , ssh_m_ad , STAT=ierr(3) )
         !
         DEALLOCATE( qsr_hc_ad, qsr_hc_b_ad, sbc_tsc_b_ad, sbc_tsc_ad, STAT=ierr(4) )

         sbc_oce_dealloc_tam = MAXVAL( ierr )
         IF( lk_mpp                  )   CALL mpp_sum ( sbc_oce_dealloc_tam )
         IF( sbc_oce_dealloc_tam > 0 )   CALL ctl_warn('sbc_oce_dealloc_tam: allocation of adjoint arrays failed')
      !
         ll_allocad = .FALSE.
      END IF
   END FUNCTION sbc_oce_dealloc_tam
   !
   SUBROUTINE sbc_oce_tam_init( kmode )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE sbc_oce_tam_init  ***
      !!
      !! ** Purpose : Allocate and initialize the tangent linear and
      !!              adjoint arrays
      !!
      !! ** Method  : kindic = 0  allocate/initialize both tl and ad variables
      !!              kindic = 1  allocate/initialize only tl variables
      !!              kindic = 2  allocate/initialize only ad variables
      !!
      !! ** Action  :
      !!
      !! References :
      !!
      !! History :
      !!        ! 2009-03 (A. Weaver) Initial version (based on oce_tam_init)
      !!        ! 2010-04 (A. Vidard) Nemo3.2 update
      !!        ! 2012-09 (P.-A. Bouttier) Nemo3.4 update
      !!-----------------------------------------------------------------------
      INTEGER, INTENT(in) :: kmode
      INTEGER :: ierr
      IF ( ( kmode == 0 ) .OR. ( kmode == 1 ) ) THEN
         IF ( .NOT. ll_alloctl ) ierr = sbc_oce_alloc_tam ( 1 )
      qsr_hc_tl(:,:,:) = 0.0_wp
      qsr_hc_b_tl(:,:,:) = 0.0_wp
      sbc_tsc_tl(:,:,:) = 0.0_wp
      sbc_tsc_b_tl(:,:,:) = 0.0_wp
      utau_tl (:,:) = 0.0_wp
      vtau_tl (:,:) = 0.0_wp
      utau_b_tl (:,:) = 0.0_wp
      vtau_b_tl (:,:) = 0.0_wp
      taum_tl (:,:) = 0.0_wp
      wndm_tl (:,:) = 0.0_wp
      qns_tl  (:,:) = 0.0_wp
      qns_b_tl  (:,:) = 0.0_wp
      rnf_tl  (:,:) = 0.0_wp
      rnf_b_tl  (:,:) = 0.0_wp
      qsr_tl  (:,:) = 0.0_wp
      qns_tot_tl  (:,:) = 0.0_wp
      qsr_tot_tl  (:,:) = 0.0_wp
      emp_tl  (:,:) = 0.0_wp
      emps_tl (:,:) = 0.0_wp
      emp_b_tl  (:,:) = 0.0_wp
      emps_b_tl (:,:) = 0.0_wp
      fr_i_tl (:,:) = 0.0_wp
      ssu_m_tl(:,:) = 0.0_wp
      ssv_m_tl(:,:) = 0.0_wp
      sst_m_tl(:,:) = 0.0_wp
      sss_m_tl(:,:) = 0.0_wp
      ssh_m_tl(:,:) = 0.0_wp
      END IF
      IF ( ( kmode == 0 ) .OR. ( kmode == 2 ) ) THEN
         IF ( .NOT. ll_allocad ) ierr = sbc_oce_alloc_tam ( 2 )
      qsr_hc_ad(:,:,:) = 0.0_wp
      qsr_hc_b_ad(:,:,:) = 0.0_wp
      sbc_tsc_ad(:,:,:) = 0.0_wp
      sbc_tsc_b_ad(:,:,:) = 0.0_wp
      utau_ad (:,:) = 0.0_wp
      vtau_ad (:,:) = 0.0_wp
      utau_b_ad (:,:) = 0.0_wp
      vtau_b_ad (:,:) = 0.0_wp
      taum_ad (:,:) = 0.0_wp
      wndm_ad (:,:) = 0.0_wp
      qns_ad  (:,:) = 0.0_wp
      qns_b_ad  (:,:) = 0.0_wp
      rnf_ad  (:,:) = 0.0_wp
      rnf_b_ad  (:,:) = 0.0_wp
      qsr_ad  (:,:) = 0.0_wp
      qns_tot_ad  (:,:) = 0.0_wp
      qsr_tot_ad  (:,:) = 0.0_wp
      emp_ad  (:,:) = 0.0_wp
      emps_ad (:,:) = 0.0_wp
      emp_b_ad  (:,:) = 0.0_wp
      emps_b_ad (:,:) = 0.0_wp
      fr_i_ad (:,:) = 0.0_wp
      ssu_m_ad(:,:) = 0.0_wp
      ssv_m_ad(:,:) = 0.0_wp
      sst_m_ad(:,:) = 0.0_wp
      sss_m_ad(:,:) = 0.0_wp
      ssh_m_ad(:,:) = 0.0_wp
      END IF
   END SUBROUTINE sbc_oce_tam_init

END MODULE sbc_oce_tam
