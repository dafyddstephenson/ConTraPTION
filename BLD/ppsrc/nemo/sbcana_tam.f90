MODULE sbcana_tam
   !!======================================================================
   !!                       ***  MODULE  sbcana  ***
   !! Ocean forcing:  analytical momentum, heat and freshwater forcings
   !!=====================================================================
   !! History of the direct module :
   !!            3.0   ! 2006-06  (G. Madec)  Original code
   !!            3.2   ! 2009-07  (G. Madec)  Style only
   !! History of the T&A module :
   !!            3.0   ! 2009-10  (F. Vigilant) original verison
   !!            3.2   ! 2020-04  (A. Vidard) nemo 3.2 update
   !!            3.4   ! 2012-07  (P.-A. Bouttier) 3.4 update
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ana  : set an analytical ocean forcing
   !!   sbc_gyre : set the GYRE configuration analytical forcing
   !!   sbc_sqb  : set the SQB configuration analytical forcing
   !!----------------------------------------------------------------------
   USE par_kind
   USE oce_tam
   USE par_oce
   USE dom_oce
   USE in_out_manager
   USE lib_mpp
   USE sbc_oce_tam
   USE tstool_tam
   USE gridrandom
   USE dotprodfld
   USE lib_fortran

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_sqb_tan
   PUBLIC   sbc_sqb_adj
   PUBLIC   sbc_ana_tan        ! routine called sbcmod_tam
   PUBLIC   sbc_ana_adj        ! routine called sbcmod_tam
   PUBLIC   sbc_gyre_tan       ! routine called sbcmod_tam
   PUBLIC   sbc_gyre_adj       ! routine called sbcmod_tam
   PUBLIC   sbc_gyre_adj_tst   ! routine called by tst

   !! * Namelist namsbc_ana
   INTEGER  ::   nn_tau000 = 1       ! nb of time-step during which the surface stress
      !                              ! increase from 0 to its nominal value
   REAL(wp) ::   rn_utau0  = 0._wp   ! constant wind stress value in i-direction
   REAL(wp) ::   rn_vtau0  = 0._wp   ! constant wind stress value in j-direction
   REAL(wp) ::   rn_qns0   = 0._wp   ! non solar heat flux
   REAL(wp) ::   rn_qsr0   = 0._wp   !     solar heat flux
   REAL(wp) ::   rn_emp0   = 0._wp   ! net freshwater flux

   REAL(wp) ::   rhoa      = 1.22_wp   ! Air density kg/m3
   REAL(wp) ::   cdrag     = 1.5e-3_wp ! drag coefficient

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

   SUBROUTINE sbc_ana_tan( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ana_tan ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   Constant and uniform surface forcing specified from
      !!              namsbc_ana namelist parameters. All the fluxes are time
      !!              independant except the stresses which increase from zero
      !!              during the first nn_tau000 time-step
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, emps
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      !!
      NAMELIST/namsbc_ana/ nn_tau000, rn_utau0, rn_vtau0, rn_qns0, rn_qsr0, rn_emp0
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         !
         REWIND ( numnam )                   ! Read Namelist namsbc : surface fluxes
         READ   ( numnam, namsbc_ana )
         !
         IF(lwp) WRITE(numout,*)' '
         IF(lwp) WRITE(numout,*)' sbc_ana_tan : Constant surface fluxes read in namsbc_ana namelist'
         IF(lwp) WRITE(numout,*)' ~~~~~~~ '
         IF(lwp) WRITE(numout,*)'              spin up of the stress  nn_tau000 = ', nn_tau000, ' time-steps'
         IF(lwp) WRITE(numout,*)'              constant i-stress      rn_utau0  = ', rn_utau0 , ' N/m2'
         IF(lwp) WRITE(numout,*)'              constant j-stress      rn_vtau0  = ', rn_vtau0 , ' N/m2'
         IF(lwp) WRITE(numout,*)'              non solar heat flux    rn_qns0   = ', rn_qns0  , ' W/m2'
         IF(lwp) WRITE(numout,*)'              solar heat flux        rn_qsr0   = ', rn_qsr0  , ' W/m2'
         IF(lwp) WRITE(numout,*)'              net heat flux          rn_emp0   = ', rn_emp0  , ' Kg/m2/s'
         !
         nn_tau000 = MAX( nn_tau000, 1 )   ! must be >= 1
         qns_tl   (:,:) = 0.0_wp
         qsr_tl   (:,:) = 0.0_wp
         emp_tl   (:,:) = 0.0_wp
         emps_tl  (:,:) = 0.0_wp
         utau_tl  (:,:) = 0.0_wp
         vtau_tl  (:,:) = 0.0_wp
         taum_tl  (:,:) = 0.0_wp
         wndm_tl  (:,:) = 0.0_wp
         !
      ENDIF
      !
   END SUBROUTINE sbc_ana_tan

   SUBROUTINE sbc_ana_adj( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ana_adj ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   Constant and uniform surface forcing specified from
      !!              namsbc_ana namelist parameters. All the fluxes are time
      !!              independant except the stresses which increase from zero
      !!              during the first nn_tau000 time-step
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, emps
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      !!
      NAMELIST/namsbc_ana/ nn_tau000, rn_utau0, rn_vtau0, rn_qns0, rn_qsr0, rn_emp0
      !!---------------------------------------------------------------------
      !
      IF( kt == nitend ) THEN
         !
         IF(lwp) WRITE(numout,*)' '
         IF(lwp) WRITE(numout,*)' sbc_ana_adj : Constant surface fluxes read in namsbc_ana namelist'
         IF(lwp) WRITE(numout,*)' ~~~~~~~ '

         nn_tau000 = MAX( nn_tau000, 1 )   ! must be >= 1
         qns_ad   (:,:) = 0.0_wp
         qsr_ad   (:,:) = 0.0_wp
         emp_ad   (:,:) = 0.0_wp
         emps_ad  (:,:) = 0.0_wp
         utau_ad  (:,:) = 0.0_wp
         vtau_ad  (:,:) = 0.0_wp
         taum_ad  (:,:) = 0.0_wp
         wndm_ad  (:,:) = 0.0_wp
         !
      ENDIF
      !
   END SUBROUTINE sbc_ana_adj


   SUBROUTINE sbc_gyre_tan( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_gyre_tam ***
      !!
      !! ** Purpose :   provide at each time-step the GYRE surface boundary
      !!              condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical seasonal cycle for GYRE configuration.
      !!                CAUTION : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.
      !!                   utau, vtau, taum, wndm, qns, qsr, emp, emps
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step
      !!
      INTEGER  ::   ji, jj                 ! dummy loop indices
      REAL(wp) ::   ztstar
      REAL(wp) ::   ztrp
      REAL(wp) ::   zsumemp_tl, zsurf
      !!---------------------------------------------------------------------

      ! ---------------------------- !
      !  heat and freshwater fluxes  !
      ! ---------------------------- !
      !same temperature, E-P as in HAZELEGER 2000

      ztrp= - 40.e0        ! retroaction term on heat fluxes (W/m2/K)
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! domain from 15 deg to 50 deg between 27 and 28  degC at 15N, -3
            ! and 13 degC at 50N 53.5 + or - 11 = 1/4 period :
            ! 64.5 in summer, 42.5 in winter
            ! 23.5 deg : tropics
            qsr_tl (ji,jj) =  0.0_wp
            qns_tl (ji,jj) = ztrp *  tsb_tl(ji,jj,1,jp_tem)
            emp_tl (ji,jj) =  0.0_wp
         END DO
      END DO
      emps_tl(:,:) = emp_tl(:,:)

      ! Compute the emp flux such as its integration on the whole domain at each time is zero
      IF( nbench /= 1 ) THEN
         zsumemp_tl = GLOB_SUM( emp_tl(:,:) )
         zsurf      = GLOB_SUM( tmask (:,:,1) )
         ! Default GYRE configuration
         zsumemp_tl = zsumemp_tl / zsurf
      ELSE
         ! Benchmark GYRE configuration (to allow the bit to bit comparison between Mpp/Mono case)
         zsumemp_tl = 0.e0   ;    zsurf = 0.e0
      ENDIF

      !salinity terms
      emp_tl (:,:) = emp_tl(:,:) - zsumemp_tl * tmask(:,:,1)
      emps_tl(:,:) = emp_tl(:,:)

   END SUBROUTINE sbc_gyre_tan

   SUBROUTINE sbc_gyre_adj( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_gyre_adj ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!      condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   analytical seasonal cycle for GYRE configuration.
      !!      * C A U T I O N : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.
      !!                   utau, vtau, qns, qsr, emp, emps
      !!
      !! Reference : Hazeleger, W., and S. Drijfhout, JPO, 30, 677-695, 2000.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step

      INTEGER  ::   ji, jj                 ! dummy loop indices
      REAL(wp) ::   ztstar
      REAL(wp) ::   ztrp
      REAL(wp) ::   zsumemp_ad, zsurf
      !!---------------------------------------------------------------------

      zsumemp_ad = 0.0_wp
      zsurf = 0.0_wp
      DO jj = 1, jpj
         DO ji = 1, jpi
            emp_ad (ji,jj) = emp_ad(ji,jj) + emps_ad(ji,jj)
            emps_ad(ji,jj) = 0.0_wp
         END DO
      END DO
      DO jj = 1, jpj
         DO ji = 1, jpi
            zsumemp_ad = zsumemp_ad - emp_ad (ji,jj) * tmask(ji,jj,1)
         END DO
      END DO

      ! Compute the emp flux such as its integration on the whole domain at each time is zero
      IF( nbench /= 1 ) THEN
         zsurf = GLOB_SUM( tmask (:,:,1) )
         ! Default GYRE configuration
         zsumemp_ad = zsumemp_ad / zsurf

         IF( lk_mpp )   CALL mpp_sum( zsurf    )       ! sum over the global domain
         IF( lk_mpp )   CALL mpp_sum( zsumemp_ad  )       ! sum over the global domain

         DO jj = 1, jpj
            DO ji = 1, jpi
               emp_ad(ji,jj) = emp_ad(ji,jj) + zsumemp_ad * tmask(ji,jj,1) * tmask_i(ji,jj)
            END DO
         END DO
         zsumemp_ad = 0.0_wp   ;   zsurf = 0.0_wp
      ELSE
         ! Benchmark GYRE configuration (to allow the bit to bit comparison between Mpp/Mono case)
         zsumemp_ad = 0.0_wp   ;    zsurf = 0.0_wp
      ENDIF

      DO jj = 1, jpj
         DO ji = 1, jpi
            emp_ad (ji,jj) = emp_ad(ji,jj) + emps_ad(ji,jj)
            emps_ad(ji,jj) = 0.0_wp
         END DO
      END DO

      ! ---------------------------- !
      !  heat and freshwater fluxes  !
      ! ---------------------------- !
      !same temperature, E-P as in HAZELEGER 2000

      ztrp= - 40.e0        ! retroaction term on heat fluxes (W/m2/K)
      DO jj = 1, jpj
         DO ji = 1, jpi
            ! domain from 15 deg to 50 deg between 27 and 28  degC at 15N, -3
            ! and 13 degC at 50N 53.5 + or - 11 = 1/4 period :
            ! 64.5 in summer, 42.5 in winter
            ! 23.5 deg : tropics
            emp_ad (ji,jj)   =  0.0_wp
            tsb_ad  (ji,jj,1,jp_tem) = tsb_ad(ji,jj,1,jp_tem) + ztrp * qns_ad (ji,jj)
            qns_ad (ji,jj)   = 0.0_wp
            qsr_ad (ji,jj)   =  0.0_wp
         END DO
      END DO

   END SUBROUTINE sbc_gyre_adj

   SUBROUTINE sbc_gyre_adj_tst ( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE example_adj_tst ***
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
      !! History :
      !!        ! 09-10 (F. Vigilant)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      !! * Local declarations
      INTEGER ::  &
         & istp,  &
         & jstp,  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zemp_tlin ,     & ! Tangent input
         & zemp_tlout,     & ! Tangent output
         & zemps_tlout,    & ! Tangent output
         & zqns_tlout,     & ! Tangent output
         & zemp_adin ,     & ! Adjoint input
         & zemps_adin,     & ! Adjoint input
         & zqns_adin ,     & ! Adjoint input
         & zemp_adout,     & ! Adjoint output
         & zr                ! 2D random field
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & ztb_tlin  ,     & ! Tangent input
         & ztb_adout ,     & ! Adjoint output
         & z3r               ! 3D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & zemp_tlin(  jpi,jpj),     &
         & ztb_tlin(jpi,jpj,jpk),    &
         & zemp_tlout( jpi,jpj),     &
         & zemps_tlout( jpi,jpj),    &
         & zqns_tlout(  jpi,jpj),    &
         & ztb_adout(jpi,jpj,jpk),   &
         & zemp_adin(  jpi,jpj),     &
         & zemps_adin( jpi,jpj),     &
         & zqns_adin(  jpi,jpj),     &
         & zemp_adout( jpi,jpj),     &
         & zr(         jpi,jpj),     &
         & z3r(     jpi,jpj,jpk)     &
         & )
      !==================================================================
      ! 1) dx = ( emp_tl, emps_tl, ssh_tl ) and
      !    dy = ( emp_tl, emps_tl )
      !==================================================================
      ! Test for time steps nit000 and nit000 + 1 (the matrix changes)

      DO jstp = nit000, nit000 + 1

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          zemp_tlin  (:,:) = 0.0_wp
          ztb_tlin (:,:,:) = 0.0_wp
          zemp_tlout (:,:) = 0.0_wp
          zemps_tlout(:,:) = 0.0_wp
          zqns_tlout (:,:) = 0.0_wp
          zemp_adin  (:,:) = 0.0_wp
          zemps_adin (:,:) = 0.0_wp
          zqns_adin  (:,:) = 0.0_wp
          zemp_adout (:,:) = 0.0_wp
          ztb_adout(:,:,:) = 0.0_wp
          z3r(:,:,:)       = 0.0_wp
          zr(:,:)          = 0.0_wp

          qns_tl (:,:)     = 0.0_wp
          qsr_tl (:,:)     = 0.0_wp
          emps_tl(:,:)     = 0.0_wp
          qsr_ad (:,:)     = 0.0_wp
          tsb_ad(:,:,:,:)     = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

          CALL grid_random( zr, 'T', 0.0_wp, stdemp )
          DO jj = nldj, nlej
             DO ji = nldi, nlei
                zemp_tlin(ji,jj) = zr(ji,jj)
             END DO
          END DO
          CALL grid_random( z3r, 'T', 0.0_wp, stdt )
          DO jk = 1, jpk
             DO jj = nldj, nlej
                DO ji = nldi, nlei
                   ztb_tlin(ji,jj,jk) = z3r(ji,jj,jk)
                END DO
             END DO
          END DO

          tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
          emp_tl (:,:) = zemp_tlin (:,:)

          CALL sbc_gyre_tan( istp )

          zemps_tlout(:,:) = emps_tl(:,:)
          zemp_tlout (:,:) = emp_tl (:,:)
          zqns_tlout(:,:) = qns_tl(:,:)

         !-----------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !-----------------------------------------------------------------

          DO jj = nldj, nlej
             DO ji = nldi, nlei
                zemp_adin( ji,jj) = zemp_tlout( ji,jj) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
                     &               * tmask(ji,jj,1)
                zemps_adin(ji,jj) = zemps_tlout(ji,jj) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
                     &               * tmask(ji,jj,1)
                zqns_adin(ji,jj) = zqns_tlout(ji,jj) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
                     &               * tmask(ji,jj,1)
             END DO
          END DO

         !-----------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !-----------------------------------------------------------------

          zsp1 = DOT_PRODUCT( zemp_tlout,  zemp_adin  )   &
               & + DOT_PRODUCT( zemps_tlout, zemps_adin ) &
               & + DOT_PRODUCT( zqns_tlout, zqns_adin )

         !-----------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !-----------------------------------------------------------------

         emp_ad (:,:) = zemp_adin (:,:)
         emps_ad(:,:) = zemps_adin(:,:)
         qns_ad( :,:) = zqns_adin( :,:)

         CALL sbc_gyre_adj ( istp )

         ztb_adout(:,:,:) = tsb_ad(:,:,:,jp_tem)
         zemp_adout (:,:) = emp_ad (:,:)

         zsp2 = DOT_PRODUCT( zemp_tlin,  zemp_adout  ) &
            & + DOT_PRODUCT( ztb_tlin,  ztb_adout  )

         ! 14 char:'12345678901234'
         IF ( jstp == nit000 ) THEN
             WRITE (cl_name,"(A14)") 'sbc_gyre_adj 1'
         ELSEIF ( jstp == nit000 + 1 ) THEN
             WRITE (cl_name,"(A14)") 'sbc_gyre_adj 2'
         END IF
!         WRITE (cl_name,"(A11,2x,i1)") 'sbc_fwb_adj',jn_fwb
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      END DO

      DEALLOCATE(       &
         & zemp_tlin,   &
         & ztb_tlin,    &
         & zemp_tlout,  &
         & zemps_tlout, &
         & zqns_tlout,  &
         & zemp_adin,   &
         & zemps_adin,  &
         & zqns_adin,   &
         & zemp_adout,  &
         & ztb_adout,   &
         & z3r,         &
         & zr           &
         & )

   END SUBROUTINE sbc_gyre_adj_tst

   SUBROUTINE sbc_sqb_tan( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_sqb ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!      condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   Constant and uniform surface forcing specified from
      !!      namsbc_ana namelist parameters. All the fluxes are time inde-
      !!      pendant except the stresses which increase from zero during
      !!      the first nn_tau000 time-step
      !!      * C A U T I O N : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.
      !!                   utau, vtau, qns, qsr, emp, emps
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      !!
      NAMELIST/namsbc_ana/ nn_tau000, rn_utau0, rn_vtau0, rn_qns0, rn_qsr0, rn_emp0
      !!---------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         !
         REWIND ( numnam )                   ! Read Namelist namsbc : surface fluxes
         READ   ( numnam, namsbc_ana )
         !
         IF(lwp) WRITE(numout,*)' '
         IF(lwp) WRITE(numout,*)' sbc_ana_tan : Constant surface fluxes read in namsbc_ana namelist'
         IF(lwp) WRITE(numout,*)' ~~~~~~~ '
         IF(lwp) WRITE(numout,*)'              spin up of the stress  nn_tau000 = ', nn_tau000, ' time-steps'
         IF(lwp) WRITE(numout,*)'              constant i-stress      rn_utau0  = ', rn_utau0 , ' N/m2'
         IF(lwp) WRITE(numout,*)'              constant j-stress      rn_vtau0  = ', rn_vtau0 , ' N/m2'
         IF(lwp) WRITE(numout,*)'              non solar heat flux    rn_qns0   = ', rn_qns0  , ' W/m2'
         IF(lwp) WRITE(numout,*)'              solar heat flux        rn_qsr0   = ', rn_qsr0  , ' W/m2'
         IF(lwp) WRITE(numout,*)'              net heat flux          rn_emp0   = ', rn_emp0  , ' Kg/m2/s'
         !
         nn_tau000 = MAX( nn_tau000, 1 )   ! must be >= 1
         qns_tl   (:,:) = 0.0_wp
         qsr_tl   (:,:) = 0.0_wp
         emp_tl   (:,:) = 0.0_wp
         emps_tl  (:,:) = 0.0_wp
         !
      ENDIF
      !
   END SUBROUTINE sbc_sqb_tan

   SUBROUTINE sbc_sqb_adj( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_ana ***
      !!
      !! ** Purpose :   provide at each time-step the ocean surface boundary
      !!      condition, i.e. the momentum, heat and freshwater fluxes.
      !!
      !! ** Method  :   Constant and uniform surface forcing specified from
      !!      namsbc_ana namelist parameters. All the fluxes are time inde-
      !!      pendant except the stresses which increase from zero during
      !!      the first nn_tau000 time-step
      !!      * C A U T I O N : never mask the surface stress field !
      !!
      !! ** Action  : - set the ocean surface boundary condition, i.e.
      !!                   utau, vtau, qns, qsr, emp, emps
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt       ! ocean time step
      !!
      NAMELIST/namsbc_ana/ nn_tau000, rn_utau0, rn_vtau0, rn_qns0, rn_qsr0, rn_emp0
      !!---------------------------------------------------------------------
      !
      IF( kt == nitend ) THEN
         !
         IF(lwp) WRITE(numout,*)' '
         IF(lwp) WRITE(numout,*)' sbc_ana_adj : Constant surface fluxes read in namsbc_ana namelist'
         IF(lwp) WRITE(numout,*)' ~~~~~~~ '

         nn_tau000 = MAX( nn_tau000, 1 )   ! must be >= 1
         qns_ad   (:,:) = 0.0_wp
         qsr_ad   (:,:) = 0.0_wp
         emp_ad   (:,:) = 0.0_wp
         emps_ad  (:,:) = 0.0_wp
         !
      ENDIF
      !
   END SUBROUTINE sbc_sqb_adj



   !!======================================================================
END MODULE sbcana_tam
