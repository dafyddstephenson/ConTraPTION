MODULE istate_tam
   !!======================================================================
   !!                     ***  MODULE  istate_tam  ***
   !! Ocean state   :  initial state setting
   !!                  Tangent and Adjoint Module
   !!=====================================================================
   !! History of the direct module:
   !!             4.0  !  89-12  (P. Andrich)  Original code
   !!             5.0  !  91-11  (G. Madec)  rewritting
   !!             6.0  !  96-01  (G. Madec)  terrain following coordinates
   !!             8.0  !  01-09  (M. Levy, M. Ben Jelloul)  istate_eel
   !!             8.0  !  01-09  (M. Levy, M. Ben Jelloul)  istate_uvg
   !!             9.0  !  03-08  (G. Madec)  F90: Free form, modules
   !!             9.0  !  03-09  (G. Madec, C. Talandier)  add EEL R5
   !!             9.0  !  04-05  (A. Koch-Larrouy)  istate_gyre
   !!             9.0  !  06-07  (S. Masson)  distributed restart using iom
   !! History of the T&A module:
   !!             9.0  !  09-04  (F. Vigilant) TAM of the 06-07 version
   !!             9.0  !  10-05  (A. Vidard) TAM - NEMO 3.2
   !!        NEMO 3.4  ! 12-07   (P.-A. Bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   istate_init_tan : initial state setting for the tangent model
   !!----------------------------------------------------------------------
   USE par_oce        , ONLY: & ! Ocean space and time domain variables
      & jpi, jpj, jpk, jpiglo
   USE oce_tam ! ocean dynamics and active tracers
   USE oce
   USE dom_oce
   USE daymod
   USE c1d
   USE restart
   USE in_out_manager
   USE zpshde_tam
   USE eosbn2_tam
   USE divcur_tam
   USE tstool_tam
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE lbclnk
   USE lbclnk_tam

   IMPLICIT NONE
   PRIVATE

   PUBLIC   istate_init_tan       ! routine called by step.F90
   PUBLIC   istate_init_adj       ! routine called by step.F90
   PUBLIC   istate_init_adj_tst   ! routine called by tst.F90

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------
! reference for s- or zps-coordinate (3D no time dependency)
! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases




   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 2528 2010-12-27 17:33:53Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
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
   !!----------------------------------------------------------------------
   !!   OPA 9.0 , LOCEAN-IPSL (2006)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE istate_init_tan( keuler )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init Tangent ***
      !!
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!----------------------------------------------------------------------
      INTEGER, INTENT (IN), OPTIONAL :: keuler

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_ini_tan : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

      CALL lbc_lnk( tsn_tl( :,:,:,1), 'T',  1. )
      CALL lbc_lnk( tsn_tl(:,:,:,2 ), 'T',  1. )
      CALL lbc_lnk( un_tl(  :,:,:  ), 'U', -1. )
      CALL lbc_lnk( vn_tl(  :,:,:  ), 'V', -1. )
      CALL lbc_lnk( sshn_tl(:,:    ), 'T',  1. )
      !
      rhd_tl  (:,:,:) = 0.0_wp
      rhop_tl (:,:,:) = 0.0_wp
      rn2_tl  (:,:,:) = 0.0_wp
      rn2b_tl (:,:,:) = 0.0_wp
      !                                    ! Start from rest
      !                                    ! ---------------
      IF ( PRESENT( keuler ) ) THEN
         neuler = keuler
      ELSE
         neuler = 0                              ! Set time-step indicator at nit000 (euler forward)
      END IF
      numror = 0                              ! define numror = 0 -> no restart file to read
      !CALL day_init()                         ! model calendar (using both namelist and restart infos)
      !                                       ! Initialization of ocean to zero
      !     before fields       !       now fields
      IF ( neuler == 0 ) THEN
         ub_tl   (:,:,:) = un_tl   (:,:,:)
         vb_tl   (:,:,:) = vn_tl   (:,:,:)
         tsb_tl   (:,:,:,:) = tsn_tl   (:,:,:,:)
         sshb_tl (  :,:) = sshn_tl (  :,:)
         CALL div_cur_tan( nit000 - 1 )
      END IF
      !
      CALL eos_tan( tsb, tsb_tl, rhd_tl, rhop_tl )        ! before potential and in situ densities
      !
      IF( ln_zps .AND. .NOT. lk_c1d )   &
            &             CALL zps_hde_tan( nit000, jpts, tsb, tsb_tl, rhd_tl,  &  ! Partial steps: before Horizontal DErivative
            &                                  gtsu_tl, gru_tl,             &  ! of t, s, rd at the bottom ocean level
            &                                  gtsv_tl, grv_tl )
   END SUBROUTINE istate_init_tan

   SUBROUTINE istate_init_adj
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE istate_init  Adjoint Module ***
      !!
      !! ** Purpose :   Initialization of the dynamics and tracer fields.
      !!----------------------------------------------------------------------
      !! * Local declarations

      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'istate_ini_adj : Initialization of the dynamics and tracers'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~'

      neuler=0
      IF( ln_zps .AND. .NOT. lk_c1d )   &
            &             CALL zps_hde_adj( nit000, jpts, tsb, tsb_ad, rhd_ad, &  ! Partial steps: before Horizontal DErivative
            &                                  gtsu_ad, gru_ad,               &  ! of t, s, rd at the bottom ocean level
            &                                  gtsv_ad, grv_ad )


      CALL eos_adj( tsb, tsb_ad, rhd_ad, rhop_ad )        ! before potential and in situ densities

      !     before fields       !       now fields
      CALL div_cur_adj( nit000 - 1 )
      un_ad   (:,:,:) = un_ad(:,:,:) + ub_ad (:,:,:)
      ub_ad   (:,:,:) = 0.0_wp
      vn_ad   (:,:,:) = vn_ad(:,:,:) + vb_ad (:,:,:)
      vb_ad   (:,:,:) = 0.0_wp
      tsn_ad   (:,:,:,:) = tsn_ad(:,:,:,:) + tsb_ad (:,:,:,:)
      tsb_ad   (:,:,:,:) = 0.0_wp
      sshn_ad   (:,:  ) = sshn_ad(:,:  ) + sshb_ad (:,:  )
      sshb_ad   (:,:  ) = 0.0_wp
      !
      rhd_ad  (:,:,:) = 0.0_wp
      rhop_ad (:,:,:) = 0.0_wp
      rn2_ad  (:,:,:) = 0.0_wp
      rn2b_ad (:,:,:) = 0.0_wp
      CALL lbc_lnk_adj( tsn_ad( :,:,:,1), 'T',  1. )
      CALL lbc_lnk_adj( tsn_ad(:,:,:,2 ), 'T',  1. )
      CALL lbc_lnk_adj( un_ad(  :,:,:  ), 'U', -1. )
      CALL lbc_lnk_adj( vn_ad(  :,:,:  ), 'V', -1. )
      CALL lbc_lnk_adj( sshn_ad(:,:    ), 'T',  1. )
      !
   END SUBROUTINE istate_init_adj

   SUBROUTINE istate_init_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE istate_init_adj_tst ***
      !!
      !! ** Purpose : Test the adjoint routine.
      !!
      !! ** Method  : Verify the scalar product
      !!
      !!                 ( L dx )^T W dy  =  dx^T L^T W dy
      !!
      !!              where  L   = tangent routine
      !!                     L^T = adjoint routine
      !!                     W   = diagonal matrix of scale factors
      !!                     dx  = input perturbation (random field)
      !!                     dy  = L dx
      !!
      !!
      !! History :
      !!        ! 2009-05 (F. Vigilant)
      !!        ! 2010-05 (A. Vidard) Update for NEMO 3.2
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER ::         &
         & ji,           &        ! dummy loop indices
         & jj,           &
         & jk

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & ztn_tlin,     & ! Tangent input: temperature
         & zsn_tlin,     & ! Tangent input: salinity
         & zun_tlin,     & ! Tangent input: velocity
         & zvn_tlin,     & ! Tangent input: velocity
         & zrotn_tlin,   & ! Tangent input: rotational
         & zdivn_tlin,   & ! Tangent input: divergence
         & ztn_adout,    & ! Adjoint output: temperature
         & zsn_adout,    & ! Adjoint output: salinity
         & zun_adout,    & ! Adjoint output: velocity
         & zvn_adout,    & ! Adjoint output: velocity
         & zrotn_adout,  & ! adjoint output: rotational
         & zdivn_adout,  & ! adjoint output: divergence
         & ztb_tlout,    & ! Tangent output: temperature
         & zsb_tlout,    & ! Tangent output: salinity
         & zub_tlout,    & ! Tangent output: velocity
         & zvb_tlout,    & ! Tangent output: velocity
         & zrhd_tlout,   & ! Tangent output:
         & zrhop_tlout,  & ! Tangent output:
         & zrotb_tlout,  & ! Tangent output:
         & zdivb_tlout,  & ! Tangent output:
         & zrn2_tlout,   & ! Tangent output:
         & zrhd_adin,    & ! Adjoint input:
         & zrhop_adin,   & ! Adjoint input:
         & ztb_adin,     & ! Adjoint input: temperature
         & zsb_adin,     & ! Adjoint input: salinity
         & zub_adin,     & ! Adjoint input: velocity
         & zvb_adin,     & ! Adjoint input: velocity
         & zrotb_adin,   & ! Adjoint input:
         & zdivb_adin,   & ! Adjoint input:
         & zrn2_adin,    & ! Adjoint input:
         & z3r             ! 3D random field

      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zsshn_tlin,   & ! Tangent input : horizontal gradient
         & zsshb_tlout,  & ! Tangent output: horizontal gradient
         & zgtu_tlout,   & ! Tangent output: horizontal gradient
         & zgtv_tlout,   & ! Tangent output: horizontal gradient
         & zgsu_tlout,   & ! Tangent output: horizontal gradient
         & zgsv_tlout,   & ! Tangent output: horizontal gradient
         & zgru_tlout,   & ! Tangent output: horizontal gradient
         & zgrv_tlout,   & ! Tangent output: horizontal gradient
         & zgtu_adin,    & ! Adjoint input : horizontal gradient
         & zgtv_adin,    & ! Adjoint input : horizontal gradient
         & zgsu_adin,    & ! Adjoint input : horizontal gradient
         & zgsv_adin,    & ! Adjoint input : horizontal gradient
         & zgru_adin,    & ! Adjoint input : horizontal gradient
         & zgrv_adin,    & ! Adjoint input : horizontal gradient
         & zsshb_adin,   & ! Adjoint input : horizontal gradient
         & zsshn_adout,  & ! Adjoint output : horizontal gradient
         & z2r             ! 2D random field

      REAL(KIND=wp) ::   &
                           ! random field standard deviation for:
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp1_3,       & !   scalar product components
         & zsp1_4,       &
         & zsp1_5,       & !   scalar product components
         & zsp1_6,       &
         & zsp1_7,       & !   scalar product components
         & zsp1_8,       &
         & zsp1_9,       &
         & zsp1_10,      &
         & zsp1_11,      &
         & zsp1_12,      &
         & zsp1_13,      &
         & zsp1_14,      &
         & zsp1_15,      &
         & zsp1_16,      &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4,       &
         & zsp2_5,       &
         & zsp2_6,       &
         & zsp2_7

      CHARACTER (LEN=14) :: &
         & cl_name

      ! Allocate memory
      ALLOCATE( &
         & ztn_tlin(jpi,jpj,jpk),     &
         & zsn_tlin(jpi,jpj,jpk),     &
         & zun_tlin(jpi,jpj,jpk),     &
         & zvn_tlin(jpi,jpj,jpk),     &
         & zsshn_tlin(jpi,jpj),       &
         & zrotn_tlin(jpi,jpj,jpk),   &
         & zdivn_tlin(jpi,jpj,jpk),   &
         & ztn_adout(jpi,jpj,jpk),    &
         & zsn_adout(jpi,jpj,jpk),    &
         & zun_adout(jpi,jpj,jpk),    &
         & zvn_adout(jpi,jpj,jpk),    &
         & zrotn_adout(jpi,jpj,jpk),  &
         & zdivn_adout(jpi,jpj,jpk),  &
         & zsshn_adout(jpi, jpj),     &
         & z3r(jpi,jpj,jpk),          &
         & z2r(jpi,jpj),              &
         & zub_tlout(jpi,jpj,jpk),    &
         & zvb_tlout(jpi,jpj,jpk),    &
         & zsb_tlout(jpi,jpj,jpk),    &
         & ztb_tlout(jpi,jpj,jpk),    &
         & zrotb_tlout(jpi,jpj,jpk),  &
         & zdivb_tlout(jpi,jpj,jpk),  &
         & zrn2_tlout(jpi,jpj,jpk),   &
         & zsshb_tlout(jpi,jpj),      &
         & zgtu_tlout(jpi,jpj),       &
         & zgtv_tlout(jpi,jpj),       &
         & zgsu_tlout(jpi,jpj),       &
         & zgsv_tlout(jpi,jpj),       &
         & zgru_tlout(jpi,jpj),       &
         & zgrv_tlout(jpi,jpj),       &
         & zrhd_tlout(jpi,jpj,jpk),   &
         & zrhop_tlout(jpi,jpj,jpk),  &
         & zub_adin(jpi,jpj,jpk),     &
         & zvb_adin(jpi,jpj,jpk),     &
         & zsb_adin(jpi,jpj,jpk),     &
         & ztb_adin(jpi,jpj,jpk),     &
         & zrotb_adin(jpi,jpj,jpk),   &
         & zdivb_adin(jpi,jpj,jpk),   &
         & zrn2_adin(jpi,jpj,jpk),    &
         & zsshb_adin(jpi,jpj),       &
         & zgtu_adin(jpi,jpj),        &
         & zgtv_adin(jpi,jpj),        &
         & zgsu_adin(jpi,jpj),        &
         & zgsv_adin(jpi,jpj),        &
         & zgru_adin(jpi,jpj),        &
         & zgrv_adin(jpi,jpj),        &
         & zrhd_adin(jpi,jpj,jpk),    &
         & zrhop_adin(jpi,jpj,jpk)    &
         & )


      !=============================================================
      ! 1) dx = ( T ) and dy = ( T )
      !=============================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      ztn_tlin    = 0.0_wp
      zsn_tlin    = 0.0_wp
      zun_tlin    = 0.0_wp
      zvn_tlin    = 0.0_wp
      zsshn_tlin  = 0.0_wp
      zrotn_tlin  = 0.0_wp
      zdivn_tlin  = 0.0_wp
      ztn_adout   = 0.0_wp
      zsn_adout   = 0.0_wp
      zun_adout   = 0.0_wp
      zvn_adout   = 0.0_wp
      zrotn_adout = 0.0_wp
      zdivn_adout = 0.0_wp
      zsshn_adout = 0.0_wp
      zub_tlout   = 0.0_wp
      zvb_tlout   = 0.0_wp
      zsb_tlout   = 0.0_wp
      ztb_tlout   = 0.0_wp
      zrotb_tlout = 0.0_wp
      zdivb_tlout = 0.0_wp
      zrn2_tlout  = 0.0_wp
      zsshb_tlout = 0.0_wp
      zgtu_tlout  = 0.0_wp
      zgtv_tlout  = 0.0_wp
      zgsu_tlout  = 0.0_wp
      zgsv_tlout  = 0.0_wp
      zgru_tlout  = 0.0_wp
      zgrv_tlout  = 0.0_wp
      zrhd_tlout  = 0.0_wp
      zrhop_tlout = 0.0_wp
      zub_adin    = 0.0_wp
      zvb_adin    = 0.0_wp
      zsb_adin    = 0.0_wp
      ztb_adin    = 0.0_wp
      zrotb_adin  = 0.0_wp
      zdivb_adin  = 0.0_wp
      zrn2_adin   = 0.0_wp
      zsshb_adin  = 0.0_wp
      zgtu_adin   = 0.0_wp
      zgtv_adin   = 0.0_wp
      zgsu_adin   = 0.0_wp
      zgsv_adin   = 0.0_wp
      zgru_adin   = 0.0_wp
      zgrv_adin   = 0.0_wp
      zrhd_adin   = 0.0_wp
      zrhop_adin  = 0.0_wp

      tsn_tl      (:,:,:,:) = 0.0_wp
      un_tl      (:,:,:) = 0.0_wp
      vn_tl      (:,:,:) = 0.0_wp
      rotn_tl    (:,:,:) = 0.0_wp
      hdivn_tl   (:,:,:) = 0.0_wp
      sshn_tl    (:,:  ) = 0.0_wp
      tsb_tl      (:,:,:,:) = 0.0_wp
      ub_tl      (:,:,:) = 0.0_wp
      vb_tl      (:,:,:) = 0.0_wp
      sshb_tl    (:,:  ) = 0.0_wp
      rhd_tl     (:,:,:) = 0.0_wp
      rhop_tl    (:,:,:) = 0.0_wp
      gtsu_tl    (:,:,:) = 0.0_wp
      gru_tl     (:,:  ) = 0.0_wp
      gtsv_tl     (:,:,:  ) = 0.0_wp
      grv_tl     (:,:  ) = 0.0_wp
      tsb_ad      (:,:,:,:) = 0.0_wp
      ub_ad      (:,:,:) = 0.0_wp
      vb_ad      (:,:,:) = 0.0_wp
      sshb_ad    (:,:  ) = 0.0_wp
      tsn_ad      (:,:,:,:) = 0.0_wp
      un_ad      (:,:,:) = 0.0_wp
      vn_ad      (:,:,:) = 0.0_wp
      sshn_ad    (:,:  ) = 0.0_wp
      gtsu_ad     (:,:,:  ) = 0.0_wp
      gtsv_ad     (:,:,:  ) = 0.0_wp

      ! Warning, following variables used by istate
      hdivn_tl           = 0.0_wp
      hdivb_tl           = 0.0_wp
      rotn_tl            = 0.0_wp
      rotb_tl            = 0.0_wp
      hdivn_ad           = 0.0_wp
      hdivb_ad           = 0.0_wp
      rotn_ad            = 0.0_wp
      rotb_ad            = 0.0_wp

      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               ztn_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zsn_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'T', 0.0_wp, stdr )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zrotn_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'T', 0.0_wp, stdr )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zdivn_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z2r, 'T', 0.0_wp, stdssh )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsshn_tlin(ji,jj) = z2r(ji,jj)
         END DO
      END DO
      CALL grid_random(  z3r, 'U', 0.0_wp, stdu )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zun_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'V', 0.0_wp, stdv )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zvn_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      tsn_tl  (:,:,:,jp_tem) = ztn_tlin  (:,:,:)
      tsn_tl  (:,:,:,jp_sal) = zsn_tlin  (:,:,:)
      rotn_tl(:,:,:) = zrotn_tlin(:,:,:)
      hdivn_tl(:,:,:)= zdivn_tlin(:,:,:)
      sshn_tl(:,:  ) = zsshn_tlin(:,:  )
      un_tl  (:,:,:) = zun_tlin  (:,:,:)
      vn_tl  (:,:,:) = zvn_tlin  (:,:,:)

     !--------------------------------------------------------------------
     ! Call the tangent routine: dy = L dx
     !--------------------------------------------------------------------

     CALL istate_init_tan

     zrhd_tlout  (:,:,:) = rhd_tl  (:,:,:)
     zrhop_tlout (:,:,:) = rhop_tl (:,:,:)
     zgtu_tlout  (:,:  ) = gtsu_tl  (:,:,jp_tem  )
     zgtv_tlout  (:,:  ) = gtsv_tl  (:,:,jp_tem  )
     zgru_tlout  (:,:  ) = gru_tl  (:,:  )
     zgrv_tlout  (:,:  ) = grv_tl  (:,:  )
     zgsu_tlout  (:,:  ) = gtsu_tl  (:,:,jp_sal  )
     zgsv_tlout  (:,:  ) = gtsv_tl  (:,:,jp_sal  )
     zsshb_tlout (:,:  ) = sshb_tl (:,:  )
     ztb_tlout   (:,:,:) = tsb_tl   (:,:,:,jp_tem)
     zsb_tlout   (:,:,:) = tsb_tl   (:,:,:,jp_sal)
     zub_tlout   (:,:,:) = ub_tl   (:,:,:)
     zvb_tlout   (:,:,:) = vb_tl   (:,:,:)
     zsshb_tlout (:,:  ) = sshb_tl (:,:  )

     !--------------------------------------------------------------------
     ! Initialize the adjoint variables: dy^* = W dy
     !--------------------------------------------------------------------

     DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zrhd_adin(ji,jj,jk)  = zrhd_tlout(ji,jj,jk) &
                 &                 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)&
                 &                 * tmask(ji,jj,jk)
              zrhop_adin(ji,jj,jk) = zrhop_tlout(ji,jj,jk) &
                 &                 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)&
                 &                 * tmask(ji,jj,jk)
              zrn2_adin(ji,jj,jk)  = zrn2_tlout(ji,jj,jk) &
                 &                 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)&
                 &                 * tmask(ji,jj,jk)
              zrotb_adin(ji,jj,jk) = zrotb_tlout(ji,jj,jk) &
                 &                 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)&
                 &                 * tmask(ji,jj,jk)
              zdivb_adin(ji,jj,jk) = zdivb_tlout(ji,jj,jk) &
                 &                 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)&
                 &                 * tmask(ji,jj,jk)
           END DO
        END DO
     END DO
     DO jj = nldj, nlej
        DO ji = nldi, nlei
           zgtu_adin(ji,jj) = zgtu_tlout(ji,jj) &
              &             * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) &
              &             * umask(ji,jj,1)
           zgsu_adin(ji,jj) = zgsu_tlout(ji,jj) &
              &             * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) &
              &             * umask(ji,jj,1)
           zgru_adin(ji,jj) = zgru_tlout(ji,jj) &
              &             * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) &
              &             * umask(ji,jj,1)
           zgtv_adin(ji,jj) = zgtv_tlout(ji,jj) &
              &             * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,1) &
              &             * vmask(ji,jj,1)
           zgsv_adin(ji,jj) = zgsv_tlout(ji,jj) &
              &             * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,1) &
              &             * vmask(ji,jj,1)
           zgrv_adin(ji,jj) = zgrv_tlout(ji,jj) &
              &             * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,1) &
              &             * vmask(ji,jj,1)
        END DO
     END DO
     DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              ztb_adin(ji,jj,jk)   = ztb_tlout(ji,jj,jk) &
                 &                 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)&
                 &                 * tmask(ji,jj,jk)
           END DO
        END DO
     END DO
     DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsb_adin(ji,jj,jk)   = zsb_tlout(ji,jj,jk) &
                 &                 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)&
                 &                 * tmask(ji,jj,jk)
           END DO
        END DO
     END DO
     DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zub_adin(ji,jj,jk)   = zub_tlout(ji,jj,jk) &
                 &                 * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk)&
                 &                 * umask(ji,jj,jk)
           END DO
        END DO
     END DO
     DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zvb_adin(ji,jj,jk)   = zvb_tlout(ji,jj,jk) &
                 &                 * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk)&
                 &                 * vmask(ji,jj,jk)
           END DO
        END DO
     END DO
     DO jj = nldj, nlej
        DO ji = nldi, nlei
           zsshb_adin(ji,jj)    = zsshb_tlout(ji,jj) &
              &                 * e1t(ji,jj) * e2t(ji,jj) * wesp_ssh &
              &                 * tmask(ji,jj,1)
        END DO
     END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1_1    = DOT_PRODUCT( zrhd_tlout   , zrhd_adin   )
      zsp1_2    = DOT_PRODUCT( zrhop_tlout  , zrhop_adin  )
      zsp1_3    = DOT_PRODUCT( zgtu_tlout   , zgtu_adin   )
      zsp1_4    = DOT_PRODUCT( zgru_tlout   , zgru_adin   )
      zsp1_5    = DOT_PRODUCT( zgsu_tlout   , zgsu_adin   )
      zsp1_6    = DOT_PRODUCT( zgtv_tlout   , zgtv_adin   )
      zsp1_7    = DOT_PRODUCT( zgrv_tlout   , zgrv_adin   )
      zsp1_8    = DOT_PRODUCT( zgsv_tlout   , zgsv_adin   )
      zsp1_9    = DOT_PRODUCT( zub_tlout    , zub_adin    )
      zsp1_10   = DOT_PRODUCT( zvb_tlout    , zvb_adin    )
      zsp1_11   = DOT_PRODUCT( ztb_tlout    , ztb_adin    )
      zsp1_12   = DOT_PRODUCT( zsb_tlout    , zsb_adin    )
      zsp1_13   = DOT_PRODUCT( zsshb_tlout  , zsshb_adin  )
      zsp1_14   = DOT_PRODUCT( zrotb_tlout  , zrotb_adin  )
      zsp1_15   = DOT_PRODUCT( zdivb_tlout  , zdivb_adin  )
      zsp1_16   = DOT_PRODUCT( zrn2_tlout   , zrn2_adin   )

      zsp1      = zsp1_1 + zsp1_2  + zsp1_3  + zsp1_4  + &
                & zsp1_5 + zsp1_6  + zsp1_7  + zsp1_8  + &
                & zsp1_9 + zsp1_10 + zsp1_11 + zsp1_12 + &
                & zsp1_13+ zsp1_14 + zsp1_15 + zsp1_16

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      rhd_ad  (:,:,:) = zrhd_adin  (:,:,:)
      rhop_ad (:,:,:) = zrhop_adin (:,:,:)
      gtsu_ad  (:,:,jp_tem  ) = zgtu_adin  (:,:  )
      gtsv_ad  (:,:,jp_tem  ) = zgtv_adin  (:,:  )
      gru_ad  (:,:  ) = zgru_adin  (:,:  )
      grv_ad  (:,:  ) = zgrv_adin  (:,:  )
      gtsu_ad  (:,:,jp_sal  ) = zgsu_adin  (:,:  )
      gtsv_ad  (:,:,jp_sal  ) = zgsv_adin  (:,:  )
      ub_ad   (:,:,:) = zub_adin   (:,:,:)
      vb_ad   (:,:,:) = zvb_adin   (:,:,:)
      tsb_ad   (:,:,:,jp_tem) = ztb_adin   (:,:,:)
      tsb_ad   (:,:,:,jp_sal) = zsb_adin   (:,:,:)
      sshb_ad (:,:  ) = zsshb_adin (:,:  )
      rotb_ad (:,:,:) = zrotb_adin (:,:,:)
      hdivb_ad(:,:,:) = zdivb_adin (:,:,:)
      rn2_ad  (:,:,:) = zrn2_adin  (:,:,:)

      CALL istate_init_adj

      ztn_adout  (:,:,:) = tsn_ad   (:,:,:,jp_tem)
      zsn_adout  (:,:,:) = tsn_ad   (:,:,:,jp_sal)
      zrotn_adout(:,:,:) = rotn_ad (:,:,:)
      zdivn_adout(:,:,:) = hdivn_ad(:,:,:)
      zun_adout  (:,:,:) = un_ad   (:,:,:)
      zvn_adout  (:,:,:) = vn_ad   (:,:,:)
      zsshn_adout(:,:  ) = sshn_ad (:,:  )

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1    = DOT_PRODUCT( ztn_tlin  , ztn_adout    )
      zsp2_2    = DOT_PRODUCT( zsn_tlin  , zsn_adout    )
      zsp2_3    = DOT_PRODUCT( zrotn_tlin, zrotn_adout  )
      zsp2_4    = DOT_PRODUCT( zdivn_tlin, zdivn_adout  )
      zsp2_5    = DOT_PRODUCT( zun_tlin  , zun_adout    )
      zsp2_6    = DOT_PRODUCT( zvn_tlin  , zvn_adout    )
      zsp2_7    = DOT_PRODUCT( zsshn_tlin, zsshn_adout  )

      zsp2      = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5 + zsp2_6 + zsp2_7

      ! Compare the scalar products

      ! 14 char:'12345678901234'
      cl_name = 'istate_tst    '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      ! Deallocate memory
      DEALLOCATE(        &
         & ztn_tlin,     &
         & zsn_tlin,     &
         & ztn_adout,    &
         & zsn_adout,    &
         & zrotn_tlin,   &
         & zdivn_tlin,   &
         & zrotn_adout,  &
         & zdivn_adout,  &
         & z3r,          &
         & z2r,          &
         & zub_tlout,    &
         & zvb_tlout,    &
         & zsb_tlout,    &
         & ztb_tlout,    &
         & zsshb_tlout,  &
         & zgtu_tlout,   &
         & zgtv_tlout,   &
         & zgsu_tlout,   &
         & zgsv_tlout,   &
         & zgru_tlout,   &
         & zgrv_tlout,   &
         & zrhd_tlout,   &
         & zrhop_tlout,  &
         & zub_adin,     &
         & zvb_adin,     &
         & zsb_adin,     &
         & ztb_adin,     &
         & zsshb_adin,   &
         & zgtu_adin,    &
         & zgtv_adin,    &
         & zgsu_adin,    &
         & zgsv_adin,    &
         & zgru_adin,    &
         & zgrv_adin,    &
         & zrhd_adin,    &
         & zrhop_adin,   &
         & zrotb_adin,   &
         & zdivb_adin,   &
         & zrn2_adin,    &
         & zrotb_tlout,  &
         & zdivb_tlout,  &
         & zrn2_tlout    &
         & )

   END SUBROUTINE istate_init_adj_tst

   !!=====================================================================
END MODULE istate_tam
