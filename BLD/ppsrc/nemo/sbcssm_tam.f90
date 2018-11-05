MODULE sbcssm_tam
   !!======================================================================
   !!                       ***  MODULE  sbcssm_tam  ***
   !! Surface module :  provide time-mean ocean surface variables
   !!                   Tangent and adjoint module
   !!======================================================================
   !! History of the direct module:
   !!            9.0   !  06-07  (G. Madec)  Original code
   !! History of the TAM module:
   !!            9.0   !  08-11  (A. Vidard) Original code
   !!            9.0   !  10-04  (A. Vidard) Nemo3.2 update
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssm_[tan adj]: calculate sea surface mean currents, temperature,
   !!                    and salinity over nn_fsbc time-step
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_kind
   USE oce_tam
   USE dom_oce
   USE sbc_oce                ! Surface boundary condition: frequency of sbc computation (as well as sea-ice model)
   USE sbc_oce_tam
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE tstool_tam

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssm_tan     ! routine called by step_tam.F90
   PUBLIC   sbc_ssm_adj     ! routine called by step_tam.F90
   PUBLIC   sbc_ssm_adj_tst ! routine called by tst.F90

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
   !!   OPA 9.0 , LOCEAN-IPSL (2006)
   !! $Id: sbcssm.F90 1196 2008-09-19 07:07:00Z ctlod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE sbc_ssm_tan( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssm_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!                provide ocean surface variable to sea-surface boundary
      !!                condition computation
      !!
      !! ** Method of the direct routine:
      !!      compute mean surface velocity (2 components at U and
      !!      V-points) [m/s], temperature [Celcius] and salinity [psu] over
      !!      the periode (kt - nn_fsbc) to kt
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt        ! ocean time step
      !
      REAL(wp) ::   zcoef       ! temporary scalar
      REAL(wp) ::   zf_sbc      ! read sbc frequency
      !!---------------------------------------------------------------------
      !                                                   ! ---------------------------------------- !
      IF( nn_fsbc == 1 ) THEN                             !      Instantaneous surface fields        !
         !                                                ! ---------------------------------------- !
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_ssm_tan: sea surface mean fields, nn_fsbc=1 : instantaneous values'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         ENDIF
         !
         ssu_m_tl(:,:) = ub_tl(:,:,1)
         ssv_m_tl(:,:) = vb_tl(:,:,1)
         sst_m_tl(:,:) = tsn_tl(:,:,1,jp_tem)
         sss_m_tl(:,:) = tsn_tl(:,:,1,jp_sal)
         ssh_m_tl(:,:) = sshn_tl(:,:)
         !
      ELSE
         !                                                ! ---------------------------------------- !
         IF( kt == nit000) THEN                           !       Initialisation: 1st time-step      !
            !                                             ! ---------------------------------------- !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_ssm_tan : sea surface mean fields'
            !
            IF( ln_rstart ) THEN
                  ssu_m_tl(:,:) = 0.0_wp
                  ssv_m_tl(:,:) = 0.0_wp
                  sst_m_tl(:,:) = 0.0_wp
                  sss_m_tl(:,:) = 0.0_wp
                  ssh_m_tl(:,:) = 0.0_wp
             ELSE
               IF(lwp) WRITE(numout,*) '~~~~~~~   mean fields initialised to instantaneous values'
               zcoef = REAL( nn_fsbc - 1, wp )
               ssu_m_tl(:,:) = zcoef * ub_tl(:,:,1)
               ssv_m_tl(:,:) = zcoef * vb_tl(:,:,1)
               sst_m_tl(:,:) = zcoef * tsn_tl(:,:,1,jp_tem)
               sss_m_tl(:,:) = zcoef * tsn_tl(:,:,1,jp_sal)
               ssh_m_tl(:,:) = zcoef * sshn_tl(:,:)
            ENDIF
            !                                             ! ---------------------------------------- !
         ELSEIF( MOD( kt - 2 , nn_fsbc ) == 0 ) THEN      !   Initialisation: New mean computation   !
            !                                             ! ---------------------------------------- !
            ssu_m_tl(:,:) = 0.0_wp      ! reset to zero ocean mean sbc fields
            ssv_m_tl(:,:) = 0.0_wp
            sst_m_tl(:,:) = 0.0_wp
            sss_m_tl(:,:) = 0.0_wp
            ssh_m_tl(:,:) = 0.0_wp
         ENDIF
         !                                                ! ---------------------------------------- !
         !                                                !        Cumulate at each time step        !
         !                                                ! ---------------------------------------- !
         ssu_m_tl(:,:) = ssu_m_tl(:,:) + ub_tl(:,:,1)
         ssv_m_tl(:,:) = ssv_m_tl(:,:) + vb_tl(:,:,1)
         sst_m_tl(:,:) = sst_m_tl(:,:) + tsn_tl(:,:,1,jp_tem)
         sss_m_tl(:,:) = sss_m_tl(:,:) + tsn_tl(:,:,1,jp_sal)
         ssh_m_tl(:,:) = ssh_m_tl(:,:) + sshn_tl(:,:)
         !                                                ! ---------------------------------------- !
         IF( MOD( kt - 1 , nn_fsbc ) == 0 ) THEN          !   Mean value at each nn_fsbc time-step   !
            !                                             ! ---------------------------------------- !
            zcoef = 1. / REAL( nn_fsbc, wp )
            sst_m_tl(:,:) = sst_m_tl(:,:) * zcoef           ! mean SST             [Celcius]
            sss_m_tl(:,:) = sss_m_tl(:,:) * zcoef           ! mean SSS             [psu]
            ssu_m_tl(:,:) = ssu_m_tl(:,:) * zcoef           ! mean suface current  [m/s]
            ssv_m_tl(:,:) = ssv_m_tl(:,:) * zcoef           !
            ssh_m_tl(:,:) = ssh_m_tl(:,:) * zcoef           !
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE sbc_ssm_tan

   SUBROUTINE sbc_ssm_adj( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssm_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!                provide ocean surface variable to sea-surface boundary
      !!                condition computation
      !!
      !! ** Method of the direct routine:
      !!      compute mean surface velocity (2 components at U and
      !!      V-points) [m/s], temperature [Celcius] and salinity [psu] over
      !!      the periode (kt - nn_fsbc) to kt
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt        ! ocean time step
      !
      REAL(wp) ::   zcoef       ! temporary scalar
      REAL(wp) ::   zf_sbc      ! read sbc frequency
      !!---------------------------------------------------------------------
      !                                                   ! ---------------------------------------- !
      IF( nn_fsbc == 1 ) THEN                             !      Instantaneous surface fields        !
         !                                                ! ---------------------------------------- !
         IF( kt == nitend) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_ssm_adj: sea surface mean fields, nn_fsbc=1 : instantaneous values'
            IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         ENDIF
         !
         ub_ad(:,:,1)  = ub_ad(:,:,1) + ssu_m_ad(:,:)
         vb_ad(:,:,1)  = vb_ad(:,:,1) + ssv_m_ad(:,:)
         tsn_ad(:,:,1,jp_tem)  = tsn_ad(:,:,1,jp_tem) + sst_m_ad(:,:)
         tsn_ad(:,:,1,jp_sal)  = tsn_ad(:,:,1,jp_sal) + sss_m_ad(:,:)
         sshn_ad(:,:)  = sshn_ad(:,:) + ssh_m_ad(:,:)
         ssu_m_ad(:,:) = 0.0_wp
         ssv_m_ad(:,:) = 0.0_wp
         sst_m_ad(:,:) = 0.0_wp
         sss_m_ad(:,:) = 0.0_wp
         ssh_m_ad(:,:) = 0.0_wp
         !
      ELSE
         !                                                ! ---------------------------------------- !
         IF( MOD( kt - 1 , nn_fsbc ) == 0 ) THEN          !   Mean value at each nn_fsbc time-step   !
            !                                             ! ---------------------------------------- !
            zcoef = 1. / REAL( nn_fsbc, wp )
            sst_m_ad(:,:) = sst_m_ad(:,:) * zcoef           ! mean SST             [Celcius]
            sss_m_ad(:,:) = sss_m_ad(:,:) * zcoef           ! mean SSS             [psu]
            ssu_m_ad(:,:) = ssu_m_ad(:,:) * zcoef           ! mean suface current  [m/s]
            ssv_m_ad(:,:) = ssv_m_ad(:,:) * zcoef           !
            ssh_m_ad(:,:) = ssh_m_ad(:,:) * zcoef           !
            !
         ENDIF
         !                                                ! ---------------------------------------- !
         !                                                !        Cumulate at each time step        !
         !                                                ! ---------------------------------------- !
         ub_ad(:,:,1) = ssu_m_ad(:,:) + ub_ad(:,:,1)
         vb_ad(:,:,1) = ssv_m_ad(:,:) + vb_ad(:,:,1)
         tsn_ad(:,:,1,jp_tem) = sst_m_ad(:,:) + tsn_ad(:,:,1,jp_tem)
         tsn_ad(:,:,1,jp_sal) = sss_m_ad(:,:) + tsn_ad(:,:,1,jp_sal)
         sshn_ad(:,:) = ssh_m_ad(:,:) + sshn_ad(:,:)
         !                                                ! ---------------------------------------- !
         IF( kt == nitend) THEN                           !       Initialisation: 1st time-step      !
            !                                             ! ---------------------------------------- !
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'sbc_ssm_adj : sea surface mean fields'
            !
         ENDIF
         !                                                ! ---------------------------------------- !
         IF( kt == nit000) THEN                           !       Initialisation: 1st time-step      !
            !                                             ! ---------------------------------------- !

            IF( ln_rstart ) THEN
                  ssu_m_ad(:,:) = 0.0_wp
                  ssv_m_ad(:,:) = 0.0_wp
                  sst_m_ad(:,:) = 0.0_wp
                  sss_m_ad(:,:) = 0.0_wp
                  ssh_m_ad(:,:) = 0.0_wp
             ELSE
               IF(lwp) WRITE(numout,*) '~~~~~~~   mean fields initialised to instantaneous values'
               zcoef = REAL( nn_fsbc - 1, wp )
               ub_ad(:,:,1)  = ub_ad(:,:,1) + zcoef * ssu_m_ad(:,:)
               vb_ad(:,:,1)  = vb_ad(:,:,1) + zcoef * ssv_m_ad(:,:)
               tsn_ad(:,:,1,jp_tem)  = tsn_ad(:,:,1,jp_tem) + zcoef * sst_m_ad(:,:)
               tsn_ad(:,:,1,jp_sal)  = tsn_ad(:,:,1,jp_sal) + zcoef * sss_m_ad(:,:)
               sshn_ad(:,:)  = sshn_ad(:,:) + zcoef * ssh_m_ad(:,:)
               ssu_m_ad(:,:) = 0.0_wp
               ssv_m_ad(:,:) = 0.0_wp
               sst_m_ad(:,:) = 0.0_wp
               sss_m_ad(:,:) = 0.0_wp
               ssh_m_ad(:,:) = 0.0_wp
            ENDIF
            !                                             ! ---------------------------------------- !
         ELSEIF( MOD( kt - 2 , nn_fsbc ) == 0 ) THEN      !   Initialisation: New mean computation   !
            !                                             ! ---------------------------------------- !
            ssu_m_ad(:,:) = 0.0_wp      ! reset to zero ocean mean sbc fields
            ssv_m_ad(:,:) = 0.0_wp
            sst_m_ad(:,:) = 0.0_wp
            sss_m_ad(:,:) = 0.0_wp
            ssh_m_ad(:,:) = 0.0_wp
         ENDIF

         !
      ENDIF
      !
   END SUBROUTINE sbc_ssm_adj

   SUBROUTINE sbc_ssm_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE sbc_ssm_adj_tst ***
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
      !! ** Action
      !!              dx = ( un_tl, vn_tl, tn_tl, sn_tl )
      !!              dy = ( ssu_m_tl, ssv_m_tl, sst_m_tl, sss_m_tl )
      !!
      !! History :
      !!        ! 08-08 (A. Vidard)
      !!        ! 09-01 (A. Weaver) cleaning
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER ::  &
         & ji,    &        ! dummy loop indices
         & jj
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zub_tlin   ,     & ! Tangent input
         & zvb_tlin   ,     & ! Tangent input
         & ztn_tlin   ,     & ! Tangent input
         & zsn_tlin   ,     & ! Tangent input
         & zsshn_tlin ,     & ! Adjoint output
         & zssum_tlin ,     & ! Tangent input
         & zssvm_tlin ,     & ! Tangent input
         & zsstm_tlin ,     & ! Tangent input
         & zsssm_tlin ,     & ! Tangent input
         & zsshm_tlin ,     & ! Tangent input
         & zssum_tlout,     & ! Tangent output
         & zssvm_tlout,     & ! Tangent output
         & zsstm_tlout,     & ! Tangent output
         & zsssm_tlout,     & ! Tangent output
         & zsshm_tlout,     & ! Tangent output
         & zub_adout  ,     & ! Adjoint input
         & zvb_adout  ,     & ! Adjoint input
         & ztn_adout  ,     & ! Adjoint input
         & zsn_adout  ,     & ! Adjoint input
         & zsshn_adout,     & ! Adjoint output
         & zssum_adout,     & ! Adjoint input
         & zssvm_adout,     & ! Adjoint input
         & zsstm_adout,     & ! Adjoint input
         & zsssm_adout,     & ! Adjoint input
         & zsshm_adout,     & ! Adjoint input
         & zssum_adin ,     & ! Adjoint output
         & zssvm_adin ,     & ! Adjoint output
         & zsstm_adin ,     & ! Adjoint output
         & zsssm_adin ,     & ! Adjoint output
         & zsshm_adin ,     & ! Adjoint output
         & zr                 ! 2D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE(                      &
         & zub_tlin   (jpi,jpj),     &
         & zvb_tlin   (jpi,jpj),     &
         & ztn_tlin   (jpi,jpj),     &
         & zsn_tlin   (jpi,jpj),     &
         & zsshn_tlin (jpi,jpj),     &
         & zssum_tlin (jpi,jpj),     &
         & zssvm_tlin (jpi,jpj),     &
         & zsstm_tlin (jpi,jpj),     &
         & zsssm_tlin (jpi,jpj),     &
         & zsshm_tlin (jpi,jpj),     &
         & zssum_tlout(jpi,jpj),     &
         & zssvm_tlout(jpi,jpj),     &
         & zsstm_tlout(jpi,jpj),     &
         & zsssm_tlout(jpi,jpj),     &
         & zsshm_tlout(jpi,jpj),     &
         & zub_adout  (jpi,jpj),     &
         & zvb_adout  (jpi,jpj),     &
         & ztn_adout  (jpi,jpj),     &
         & zsn_adout  (jpi,jpj),     &
         & zsshn_adout(jpi,jpj),     &
         & zssum_adout(jpi,jpj),     &
         & zssvm_adout(jpi,jpj),     &
         & zsstm_adout(jpi,jpj),     &
         & zsssm_adout(jpi,jpj),     &
         & zsshm_adout(jpi,jpj),     &
         & zssum_adin (jpi,jpj),     &
         & zssvm_adin (jpi,jpj),     &
         & zsstm_adin (jpi,jpj),     &
         & zsssm_adin (jpi,jpj),     &
         & zsshm_adin (jpi,jpj),     &
         & zr         (jpi,jpj)      &
         & )
      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, tn_tl, sn_tl ) and
      !    dy = ( ssu_m_tl, ssv_m_tl, sst_m_tl, sss_m_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      zub_tlin   (:,:) = 0.0_wp
      zvb_tlin   (:,:) = 0.0_wp
      ztn_tlin   (:,:) = 0.0_wp
      zsn_tlin   (:,:) = 0.0_wp
      zssum_tlin (:,:) = 0.0_wp
      zssvm_tlin (:,:) = 0.0_wp
      zsstm_tlin (:,:) = 0.0_wp
      zsssm_tlin (:,:) = 0.0_wp
      zsshm_tlin (:,:) = 0.0_wp
      zssum_tlout(:,:) = 0.0_wp
      zssvm_tlout(:,:) = 0.0_wp
      zsstm_tlout(:,:) = 0.0_wp
      zsssm_tlout(:,:) = 0.0_wp
      zsshm_tlout(:,:) = 0.0_wp
      zub_adout  (:,:) = 0.0_wp
      zvb_adout  (:,:) = 0.0_wp
      ztn_adout  (:,:) = 0.0_wp
      zsn_adout  (:,:) = 0.0_wp
      zsshn_adout(:,:) = 0.0_wp
      zssum_adout(:,:) = 0.0_wp
      zssvm_adout(:,:) = 0.0_wp
      zsstm_adout(:,:) = 0.0_wp
      zsssm_adout(:,:) = 0.0_wp
      zsshm_adout(:,:) = 0.0_wp
      zssum_adin (:,:) = 0.0_wp
      zssvm_adin (:,:) = 0.0_wp
      zsstm_adin (:,:) = 0.0_wp
      zsssm_adin (:,:) = 0.0_wp
      zsshm_adin (:,:) = 0.0_wp
      zr         (:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random( zr, 'U', 0.0_wp, stdu )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zub_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      CALL grid_random(  zr, 'V', 0.0_wp, stdv )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zvb_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      CALL grid_random(  zr, 'T', 0.0_wp, stds )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsn_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            ztn_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stdssh )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsshn_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      CALL grid_random(  zr, 'U', 0.0_wp, stdu )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zssum_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO
      CALL grid_random(  zr, 'V', 0.0_wp, stdv )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zssvm_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsstm_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      CALL grid_random(  zr, 'T', 0.0_wp, stds )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsssm_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      CALL grid_random(  zr, 'T', 0.0_wp, stdssh )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsshm_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      ub_tl (:,:,1) = zub_tlin  (:,:)
      vb_tl (:,:,1) = zvb_tlin  (:,:)
      tsn_tl (:,:,1,jp_tem) = ztn_tlin  (:,:)
      tsn_tl (:,:,1,jp_sal) = zsn_tlin  (:,:)
      sshn_tl (:,:) = zsshn_tlin(:,:)
      ssu_m_tl(:,:) = zssum_tlin(:,:)
      ssv_m_tl(:,:) = zssvm_tlin(:,:)
      sst_m_tl(:,:) = zsstm_tlin(:,:)
      sss_m_tl(:,:) = zsssm_tlin(:,:)
      ssh_m_tl(:,:) = zsshm_tlin(:,:)

      CALL sbc_ssm_tan( nit000 + 1 )

      zssum_tlout (:,:) = ssu_m_tl(:,:)
      zssvm_tlout (:,:) = ssv_m_tl(:,:)
      zsstm_tlout (:,:) = sst_m_tl(:,:)
      zsssm_tlout (:,:) = sss_m_tl(:,:)
      zsshm_tlout (:,:) = ssh_m_tl(:,:)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zssum_adin(ji,jj) = zssum_tlout(ji,jj) &
               &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) &
               &               * umask(ji,jj,1)
            zssvm_adin(ji,jj) = zssvm_tlout(ji,jj) &
               &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) &
               &               * umask(ji,jj,1)
            zsstm_adin(ji,jj) = zsstm_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
            zsssm_adin(ji,jj) = zsssm_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_s(1)
            zsshm_adin(ji,jj) = zsshm_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_s(1)
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 = DOT_PRODUCT( zssum_tlout, zssum_adin ) &
         & + DOT_PRODUCT( zssvm_tlout, zssvm_adin ) &
         & + DOT_PRODUCT( zsstm_tlout, zsstm_adin ) &
         & + DOT_PRODUCT( zsssm_tlout, zsssm_adin ) &
         & + DOT_PRODUCT( zsshm_tlout, zsshm_adin )

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      ssu_m_ad(:,:) = zssum_adin(:,:)
      ssv_m_ad(:,:) = zssvm_adin(:,:)
      sst_m_ad(:,:) = zsstm_adin(:,:)
      sss_m_ad(:,:) = zsssm_adin(:,:)
      ssh_m_ad(:,:) = zsshm_adin(:,:)
      ub_ad(:,:,1)  = 0.0_wp
      vb_ad(:,:,1)  = 0.0_wp
      tsn_ad(:,:,1,jp_tem)  = 0.0_wp
      tsn_ad(:,:,1,jp_tem)  = 0.0_wp
      sshn_ad(:,:)  = 0.0_wp

      CALL sbc_ssm_adj( nit000 + 1 )

      zub_adout  (:,:) = ub_ad(:,:,1)
      zvb_adout  (:,:) = vb_ad(:,:,1)
      ztn_adout  (:,:) = tsn_ad(:,:,1,jp_tem)
      zsn_adout  (:,:) = tsn_ad(:,:,1,jp_sal)
      zsshn_adout(:,:) = sshn_ad(:,:)
      zssum_adout(:,:) = ssu_m_ad(:,:)
      zssvm_adout(:,:) = ssv_m_ad(:,:)
      zsstm_adout(:,:) = sst_m_ad(:,:)
      zsssm_adout(:,:) = sss_m_ad(:,:)
      zsshm_adout(:,:) = ssh_m_ad(:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T dx^*
      !--------------------------------------------------------------------

      zsp2 = DOT_PRODUCT( zub_tlin  , zub_adout   ) &
         & + DOT_PRODUCT( zvb_tlin  , zvb_adout   ) &
         & + DOT_PRODUCT( ztn_tlin  , ztn_adout   ) &
         & + DOT_PRODUCT( zsn_tlin  , zsn_adout   ) &
         & + DOT_PRODUCT( zsshn_tlin, zsshn_adout ) &
         & + DOT_PRODUCT( zssum_tlin, zssum_adout ) &
         & + DOT_PRODUCT( zssvm_tlin, zssvm_adout ) &
         & + DOT_PRODUCT( zsstm_tlin, zsstm_adout ) &
         & + DOT_PRODUCT( zsssm_tlin, zsssm_adout ) &
         & + DOT_PRODUCT( zsshm_tlin, zsshm_adout )

      ! 14 char:'12345678901234'
      cl_name = 'sbc_ssm_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(           &
         & zub_tlin   ,     & ! Tangent input
         & zvb_tlin   ,     & ! Tangent input
         & ztn_tlin   ,     & ! Tangent input
         & zsn_tlin   ,     & ! Tangent input
         & zssum_tlin ,     & ! Tangent input
         & zssvm_tlin ,     & ! Tangent input
         & zsstm_tlin ,     & ! Tangent input
         & zsssm_tlin ,     & ! Tangent input
         & zsshm_tlin ,     & ! Tangent input
         & zssum_tlout,     & ! Tangent output
         & zssvm_tlout,     & ! Tangent output
         & zsstm_tlout,     & ! Tangent output
         & zsssm_tlout,     & ! Tangent output
         & zsshm_tlout,     & ! Tangent output
         & zub_adout  ,     & ! Adjoint input
         & zvb_adout  ,     & ! Adjoint input
         & ztn_adout  ,     & ! Adjoint input
         & zsn_adout  ,     & ! Adjoint input
         & zssum_adout,     & ! Adjoint input
         & zssvm_adout,     & ! Adjoint input
         & zsstm_adout,     & ! Adjoint input
         & zsssm_adout,     & ! Adjoint input
         & zsshm_adout,     & ! Adjoint input
         & zssum_adin ,     & ! Adjoint output
         & zssvm_adin ,     & ! Adjoint output
         & zsstm_adin ,     & ! Adjoint output
         & zsssm_adin ,     & ! Adjoint output
         & zsshm_adin ,     & ! Adjoint output
         & zr               &
         & )

   END SUBROUTINE sbc_ssm_adj_tst
   !!======================================================================
END MODULE sbcssm_tam
