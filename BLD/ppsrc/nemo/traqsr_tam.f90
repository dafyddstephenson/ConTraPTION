MODULE traqsr_tam
   !!======================================================================
   !!                       ***  MODULE  traqsr_tam  ***
   !! Ocean physics: solar radiation penetration in the top ocean levels
   !!                Tangent and Adjoint Module
   !!======================================================================
   !! History :  OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1996-01  (G. Madec)  s-coordinates
   !!   NEMO     1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!             -   !  2005-11  (G. Madec) zco, zps, sco coordinate
   !!            3.2  !  2009-04  (G. Madec & NEMO team)
   !! History of the TAM:
   !!                 !  2008-05  (A. Vidard) Skeleton
   !!            3.0  !  2008-09  (A. Vidard)   TAM of the 2005-11 version
   !!            3.2  !  2010-03  (F. Vigilant) TAM of the 2009-11 version
   !!            3.4  !  2012-07  (P.-A. Bouttier) 3.4 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_qsr      : trend due to the solar radiation penetration
   !!   tra_qsr_init : solar radiation penetration initialization
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE in_out_manager
   USE fldread
   USE sbc_oce
   USE sbc_oce_tam
   USE phycst
   USE prtctl
   USE gridrandom
   USE dotprodfld
   USE traqsr
   USE trc_oce
   USE trc_oce_tam
   USE tstool_tam
   USE lib_mpp
   USE wrk_nemo
   USE timing
   USE restart
   USE fldread
   USE iom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_qsr_tan      ! routine called by step_tam.F90 (ln_traqsr=T)
   PUBLIC   tra_qsr_adj      ! routine called by step_tam.F90 (ln_traqsr=T)
   PUBLIC   tra_qsr_init_tam
   PUBLIC   tra_qsr_adj_tst  ! routine called by tst.F90

   REAL(wp) :: xsi0r
   REAL(wp) :: xsi1r
   REAL(wp), DIMENSION(3,61) :: rkrgb
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

CONTAINS

   SUBROUTINE tra_qsr_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_tan  ***
      !!
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!      penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 3 wavebands (RGB) and a ratio rn_abs
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration
      !!         is given by : zta = 1/e3t dk[ I ] / (rau0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         In z-coordinate case, the computation is only done down to the
      !!      level where I(k) < 1.e-15 W/m2. In addition, the coefficients
      !!      used for the computation are calculated one for once as they
      !!      depends on k only.
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!
      !! Reference  : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!              Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt        ! ocean time-step
      !
      !!
      INTEGER  ::    ji, jj, jk          ! dummy loop indexes
      INTEGER  ::   irgb                 ! temporary integers
      REAL(wp) ::   zchl, zcoef, zfact, z1_e3t   ! temporary scalars
      REAL(wp) ::   zc0, zc1, zc2, zc3   !    -         -
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zekb, zekg, zekr                      ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   ze0tl, ze1tl , ze2tl, ze3tl, zeatl    ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_qsr_tan')
      !
      CALL wrk_alloc( jpi, jpj,      zekb, zekg, zekr        )
      CALL wrk_alloc( jpi, jpj, jpk, ze0tl, ze1tl, ze2tl, ze3tl, zeatl )
      !

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_qsr_tan : penetration of the surface solar radiation'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF( .NOT.ln_traqsr )   RETURN
      ENDIF
      !                                        Set before qsr tracer content field
      !                                        ***********************************
      IF( kt == nit000 ) THEN                     ! Set the forcing field at nit000 - 1
         !                                        ! -----------------------------------
         IF( ln_rstart ) THEN !.AND.    &                    ! Restart: read in restart file
              !& iom_varid( numror, 'qsr_hc_b', ldstop = .FALSE. ) > 0 ) THEN
            !IF(lwp) WRITE(numout,*) '          nit000-1 qsr tracer content forcing field red in the restart file'
            zfact = 0.5e0
            !CALL iom_get( numror, jpdom_autoglo, 'qsr_hc_b', qsr_hc_b_tl )   ! before heat content trend due to Qsr flux
         ELSE                                           ! No restart or restart not found: Euler forward time stepping
            zfact = 1.e0
         ENDIF
         qsr_hc_b_tl(:,:,:) = 0.e0
      ELSE                                        ! Swap of forcing field
         !                                        ! ---------------------
         zfact = 0.5e0
         qsr_hc_b_tl(:,:,:) = qsr_hc_tl(:,:,:)
      ENDIF
      !                                        Compute now qsr tracer content field
      !
      !                                           ! ============================================== !
      IF( lk_qsr_bio .AND. ln_qsr_bio ) THEN      !  bio-model fluxes  : all vertical coordinates  !
         !                                        ! ============================================== !
         DO jk = 1, jpkm1
            qsr_hc_tl(:,:,jk) = ro0cpr * ( etot3_tl(:,:,jk) - etot3_tl(:,:,jk+1) )
         END DO

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z1_e3t = zfact / e3t(ji,jj,jk)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem) + ( qsr_hc_b_tl(ji,jj,jk) + qsr_hc_tl(ji,jj,jk) ) * z1_e3t
               END DO
            END DO
         END DO
         !                                        ! ============================================== !
      ELSE                                        !  Ocean alone :
         !                                        ! ============================================== !
         !
         !                                                ! ------------------------- !
         IF( ln_qsr_rgb) THEN                             !  R-G-B  light penetration !
            !                                             ! ------------------------- !
            IF( nn_chldta == 1 .OR. lk_vvl ) THEN            !*  Variable Chlorophyll or ocean volume
               !!! Set chlorophyl concentration
               !!IF( nn_chldta ==1 ) THEN                             !*  Variable Chlorophyll
                  !!!
                  !!CALL fld_read( kt, 1, sf_chl )                         ! Read Chl data and provides it at the current time step
                  !!!
   !!!CDIR COLLAPSE
   !!!CDIR NOVERRCHK
                  !!DO jj = 1, jpj                                         ! Separation in R-G-B depending of the surface Chl
   !!!CDIR NOVERRCHK
                     !!DO ji = 1, jpi
                        !!zchl = MIN( 10.0_wp , MAX( 0.03_wp, sf_chl(1)%fnow(ji,jj) ) )
                        !!irgb = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )
                        !!zekb(ji,jj) = rkrgb(1,irgb)
                        !!zekg(ji,jj) = rkrgb(2,irgb)
                        !!zekr(ji,jj) = rkrgb(3,irgb)
                     !!END DO
                  !!END DO
               !ELSE
                  !zchl = 0.05                                     ! constant chlorophyll
                  !irgb = NINT( 41 + 20.*LOG10( zchl ) + 1.e-15 )
                  !zekb(:,:) = rkrgb(1,irgb)                       ! Separation in R-G-B depending of the chlorophyll
                  !zekg(:,:) = rkrgb(2,irgb)
                  !zekr(:,:) = rkrgb(3,irgb)
               !ENDIF
               !
               !zcoef  = ( 1.0_wp - rn_abs ) / 3.0_wp                        ! equi-partition in R-G-B
               !ze0tl(:,:,1) = rn_abs  * qsr_tl(:,:)
               !ze1tl(:,:,1) = zcoef * qsr_tl(:,:)
               !ze2tl(:,:,1) = zcoef * qsr_tl(:,:)
               !ze3tl(:,:,1) = zcoef * qsr_tl(:,:)
               !zeatl(:,:,1) =         qsr_tl(:,:)
               !!
               !DO jk = 2, nksr+1
!!CDIR NOVERRCHK
                  !DO jj = 1, jpj
!!CDIR NOVERRCHK
                     !DO ji = 1, jpi
                        !zc0tl = ze0tl(ji,jj,jk-1) * EXP( - e3t(ji,jj,jk-1) * xsi0r     )
                        !zc1tl = ze1tl(ji,jj,jk-1) * EXP( - e3t(ji,jj,jk-1) * zekb(ji,jj) )
                        !zc2tl = ze2tl(ji,jj,jk-1) * EXP( - e3t(ji,jj,jk-1) * zekg(ji,jj) )
                        !zc3tl = ze3tl(ji,jj,jk-1) * EXP( - e3t(ji,jj,jk-1) * zekr(ji,jj) )
                        !ze0tl(ji,jj,jk) = zc0tl
                        !ze1tl(ji,jj,jk) = zc1tl
                        !ze2tl(ji,jj,jk) = zc2tl
                        !ze3tl(ji,jj,jk) = zc3tl
                        !zeatl(ji,jj,jk) = ( zc0tl + zc1tl + zc2tl + zc3tl ) * tmask(ji,jj,jk)
                     !END DO
                  !END DO
               !END DO
               !!
               !DO jk = 1, nksr                                        ! compute and add qsr trend to ta
                  !qsr_tl(:,:) = ro0cpr * ( zeatl(:,:,jk) - zeatl(:,:,jk+1) )
               !END DO
               !zeatl(:,:,nksr+1:jpk) = 0.0_wp     ! below 400m set to zero
               !!
               CALL ctl_stop('tra_qsr_tan: key_vvl or non-constant chlorophyll management(nn_chldta = 1) &
                          &   not implemented in TAM yet')
            ELSE                                                 !*  Constant Chlorophyll
               DO jk = 1, nksr
                  qsr_hc_tl(:,:,jk) = etot3_tl(:,:,jk) * qsr(:,:) + etot3(:,:,jk) * qsr_tl(:,:)
               END DO
            ENDIF

         ENDIF
         !                                                ! ------------------------- !
         IF( ln_qsr_2bd ) THEN                            !  2 band light penetration !
            !                                             ! ------------------------- !
            !
            IF( lk_vvl ) THEN                                  !* variable volume
               !zz0   =        rn_abs   * ro0cpr
               !zz1   = ( 1. - rn_abs ) * ro0cpr
               !DO jk = 1, nksr                    ! solar heat absorbed at T-point in the top 400m
                  !DO jj = 1, jpj
                     !DO ji = 1, jpi
                        !zc0 = zz0 * EXP( -gdepw(ji,jj,jk  )*xsi0r ) + zz1 * EXP( -gdepw(ji,jj,jk  )*xsi1r )
                        !zc1 = zz0 * EXP( -gdepw(ji,jj,jk+1)*xsi0r ) + zz1 * EXP( -gdepw(ji,jj,jk+1)*xsi1r )
                        !qsr_hc_tl(ji,jj,jk) = qsr_tl(ji,jj) * ( zc0*tmask(ji,jj,jk) - zc1*tmask(ji,jj,jk+1) )
                     !END DO
                  !END DO
               !END DO
               CALL ctl_stop('tra_qsr_tan: key_vvl or chlorophyll management not implemented in TAM yet')
            ELSE                                               !* constant volume: coef. computed one for all
               DO jk = 1, nksr
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        qsr_hc_tl(ji,jj,jk) =  etot3_tl(ji,jj,jk) * qsr(ji,jj) + etot3(ji,jj,jk) * qsr_tl(ji,jj)
                     END DO
                  END DO
               END DO
               !
            ENDIF
            !
         ENDIF
         DO jk = 1, nksr
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z1_e3t = zfact / e3t(ji,jj,jk)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem) + ( qsr_hc_b_tl(ji,jj,jk) + qsr_hc_tl(ji,jj,jk) ) * z1_e3t
               END DO
            END DO
         END DO
         !
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj,      zekb, zekg, zekr        )
      CALL wrk_dealloc( jpi, jpj, jpk, ze0tl, ze1tl, ze2tl, ze3tl, zeatl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_qsr_tan')
      !
   END SUBROUTINE tra_qsr_tan
   SUBROUTINE tra_qsr_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_adj  ***
      !!
      !! ** Purpose :   Compute the temperature trend due to the solar radiation
      !!      penetration and add it to the general temperature trend.
      !!
      !! ** Method  : The profile of the solar radiation within the ocean is defined
      !!      through 2 wavebands (rn_si0,rn_si1) or 3 wavebands (RGB) and a ratio rn_abs
      !!      Considering the 2 wavebands case:
      !!         I(k) = Qsr*( rn_abs*EXP(z(k)/rn_si0) + (1.-rn_abs)*EXP(z(k)/rn_si1) )
      !!         The temperature trend associated with the solar radiation penetration
      !!         is given by : zta = 1/e3t dk[ I ] / (rau0*Cp)
      !!         At the bottom, boudary condition for the radiation is no flux :
      !!      all heat which has not been absorbed in the above levels is put
      !!      in the last ocean level.
      !!         In z-coordinate case, the computation is only done down to the
      !!      level where I(k) < 1.e-15 W/m2. In addition, the coefficients
      !!      used for the computation are calculated one for once as they
      !!      depends on k only.
      !!
      !! ** Action  : - update ta with the penetrative solar radiation trend
      !!
      !! Reference  : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!              Lengaigne et al. 2007, Clim. Dyn., V28, 5, 503-516.
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt     ! ocean time-step
      !
      !!
      INTEGER  ::    ji, jj, jk       ! dummy loop indexes
      INTEGER  ::   irgb                 ! temporary integers
      REAL(wp) ::   zchl, zcoef, zfact, z1_e3t   ! temporary scalars
      REAL(wp) ::   zc0, zc1, zc2, zc3, zz0, zz1   !    -         -
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zekb, zekg, zekr                      ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ze0ad, ze1ad , ze2ad, ze3ad, zeaad    ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_qsr_adj')
      !
      CALL wrk_alloc( jpi, jpj,      zekb, zekg, zekr        )
      CALL wrk_alloc( jpi, jpj, jpk, ze0ad, ze1ad, ze2ad, ze3ad, zeaad )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_qsr_adj : penetration of the surface solar radiation'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF( .NOT.ln_traqsr )   RETURN
      ENDIF
      !                                        Set before qsr tracer content field
      !                                        ***********************************
      IF( kt == nit000 ) THEN                     ! Set the forcing field at nit000 - 1
         !                                        ! -----------------------------------
         IF( ln_rstart ) THEN !.AND.    &                    ! Restart: read in restart file
              !& iom_varid( numror, 'qsr_hc_b', ldstop = .FALSE. ) > 0 ) THEN
            !IF(lwp) WRITE(numout,*) '          nit000-1 qsr tracer content forcing field red in the restart file'
            zfact = 0.5e0
            !CALL iom_get( numror, jpdom_autoglo, 'qsr_hc_b', qsr_hc_b_ad )   ! before heat content trend due to Qsr flux
         ELSE                                           ! No restart or restart not found: Euler forward time stepping
            zfact = 1.e0
         ENDIF
      ELSE                                        ! Swap of forcing field
         !                                        ! ---------------------
         zfact = 0.5e0
      ENDIF
      !                                           ! ============================================== !
      IF( lk_qsr_bio .AND. ln_qsr_bio ) THEN      !  bio-model fluxes  : all vertical coordinates  !
         !                                        ! ============================================== !
         DO jk = jpkm1, 1, -1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z1_e3t = zfact / e3t(ji,jj,jk)
                  qsr_hc_b_ad(ji,jj,jk) = qsr_hc_b_ad(ji,jj,jk) + tsa_ad(ji,jj,jk,jp_tem) * z1_e3t
                  qsr_hc_ad(ji,jj,jk) = qsr_hc_ad(ji,jj,jk) + tsa_ad(ji,jj,jk,jp_tem) * z1_e3t
               END DO
            END DO
         END DO
         DO jk = jpkm1, 1, -1
            etot3_ad(:,:,jk)   = etot3_ad(:,:,jk) + ro0cpr * qsr_hc_ad(:,:,jk)
            etot3_ad(:,:,jk+1) = etot3_ad(:,:,jk+1) - ro0cpr * qsr_hc_ad(:,:,jk)
         END DO
         !                                        ! ============================================== !
      ELSE                                        !  Ocean alone :
         !                                        ! ============================================== !
         !
         DO jk = 1, nksr
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  z1_e3t = zfact / e3t(ji,jj,jk)
                  qsr_hc_b_ad(ji,jj,jk) = qsr_hc_b_ad(ji,jj,jk) + tsa_ad(ji,jj,jk,jp_tem) * z1_e3t
                  qsr_hc_ad(ji,jj,jk)   = qsr_hc_ad(ji,jj,jk)   + tsa_ad(ji,jj,jk,jp_tem) * z1_e3t
               END DO
            END DO
         END DO
         !                                                ! ------------------------- !
         IF( ln_qsr_2bd ) THEN                            !  2 band light penetration !
            !                                             ! ------------------------- !
            !
            IF( lk_vvl ) THEN                                  !* variable volume
               !zz0   =        rn_abs   * ro0cpr
               !zz1   = ( 1. - rn_abs ) * ro0cpr
               !DO jk = nksr, 1, -1                    ! solar heat absorbed at T-point in the top 400m
                  !DO jj = 1, jpj
                     !DO ji = 1, jpi
                        !zc0 = zz0 * EXP( -gdepw(ji,jj,jk  )*xsi0r ) + zz1 * EXP( -gdepw(ji,jj,jk  )*xsi1r )
                        !zc1 = zz0 * EXP( -gdepw(ji,jj,jk+1)*xsi0r ) + zz1 * EXP( -gdepw(ji,jj,jk+1)*xsi1r )
                        !qsr_ad(ji,jj) = qsr_hc_ad(ji,jj) * ( zc0*tmask(ji,jj,jk) - zc1*tmask(ji,jj,jk+1) )
                     !END DO
                  !END DO
               !END DO
               CALL ctl_stop('tra_qsr_adj: key_vvl or chlorophyll management not implemented in TAM yet')
            ELSE                                               !* constant volume: coef. computed one for all
               DO jk = 1, nksr
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        etot3_ad(ji,jj,jk) = etot3_ad(ji,jj,jk) + qsr(ji,jj) * qsr_hc_ad(ji,jj,jk)
                        qsr_ad(ji,jj)      = qsr_ad(ji,jj) + etot3(ji,jj,jk) * qsr_hc_ad(ji,jj,jk)
                        qsr_hc_ad(ji,jj,jk) = 0._wp
                     END DO
                  END DO
               END DO
               !
            ENDIF
            !
         ENDIF
         !
         !                                                ! ------------------------- !
         IF( ln_qsr_rgb) THEN                             !  R-G-B  light penetration !
            !                                             ! ------------------------- !
            ! Set chlorophyl concentration
            IF( nn_chldta == 1 .OR. lk_vvl ) THEN            !*  Variable Chlorophyll or ocean volume
               !!!
               !!IF( nn_chldta ==1 ) THEN                             !*  Variable Chlorophyll
                  !!zc0ad = 0.0_wp; zc1ad = 0.0_wp; zc2ad = 0.0_wp; zc3ad = 0.0_wp
                  !!ze0ad(:,:,:) = 0.0_wp; ze1ad(:,:,:) = 0.0_wp; ze2ad(:,:,:) = 0.0_wp; ze3ad(:,:,:) = 0.0_wp
                  !!zeaad(:,:,:) = 0.0_wp
                  !!!
                  !!CALL fld_read( kt, 1, sf_chl )                         ! Read Chl data and provides it at the current time step
                  !!!
   !!!CDIR COLLAPSE
   !!!CDIR NOVERRCHK
                  !!DO jj = 1, jpj                                         ! Separation in R-G-B depending of the surface Chl
   !!!CDIR NOVERRCHK
                     !!DO ji = 1, jpi
                        !!zchl = MIN( 10.0_wp , MAX( 0.03_wp, sf_chl(1)%fnow(ji,jj) ) )
                        !!irgb = NINT( 41 + 20.*LOG10(zchl) + 1.e-15 )
                        !!zekb(ji,jj) = rkrgb(1,irgb)
                        !!zekg(ji,jj) = rkrgb(2,irgb)
                        !!zekr(ji,jj) = rkrgb(3,irgb)
                     !!END DO
                  !!END DO
               !ELSE
                  !zchl = 0.05                                     ! constant chlorophyll
                  !irgb = NINT( 41 + 20.*LOG10( zchl ) + 1.e-15 )
                  !zekb(:,:) = rkrgb(1,irgb)                       ! Separation in R-G-B depending of the chlorophyll
                  !zekg(:,:) = rkrgb(2,irgb)
                  !zekr(:,:) = rkrgb(3,irgb)
               !ENDIF
               !!
               !zcoef  = ( 1.0_wp - rn_abs ) / 3.0_wp

               !zeaad(:,:,nksr+1:jpk) = 0.0_wp     ! below 400m set to zero
               !!
               !DO jk = 1, nksr                                        ! compute and add qsr trend to ta
                  !zeaad(:,:,jk  ) =   ro0cpr * qsr_hc_ad(:,:,jk)
                  !zeaad(:,:,jk+1) = - ro0cpr * qsr_hc_ad(:,:,jk)
               !END DO
               !!
               !DO jk = nksr+1, 2, -1
!!CDIR NOVERRCHK
                  !DO jj = 1, jpj
!!CDIR NOVERRCHK
                     !DO ji = 1, jpi
                        !zc0ad = zc0ad + zeaad(ji,jj,jk) * tmask(ji,jj,jk)
                        !zc1ad = zc1ad + zeaad(ji,jj,jk) * tmask(ji,jj,jk)
                        !zc2ad = zc2ad + zeaad(ji,jj,jk) * tmask(ji,jj,jk)
                        !zc3ad = zc3ad + zeaad(ji,jj,jk) * tmask(ji,jj,jk)
                        !zeaad(ji,jj,jk) = 0.0_wp
                        !zc0ad = zc0ad + ze0ad(ji,jj,jk)
                        !zc1ad = zc1ad + ze1ad(ji,jj,jk)
                        !zc2ad = zc2ad + ze2ad(ji,jj,jk)
                        !zc3ad = zc3ad + ze3ad(ji,jj,jk)
                        !ze0ad(ji,jj,jk) = 0.0_wp
                        !ze1ad(ji,jj,jk) = 0.0_wp
                        !ze2ad(ji,jj,jk) = 0.0_wp
                        !ze3ad(ji,jj,jk) = 0.0_wp
                        !ze0ad(ji,jj,jk-1) = ze0ad(ji,jj,jk-1) + zc0ad * EXP( - e3t(ji,jj,jk-1) * xsi0r     )
                        !ze1ad(ji,jj,jk-1) = ze1ad(ji,jj,jk-1) + zc1ad * EXP( - e3t(ji,jj,jk-1) * zekb(ji,jj) )
                        !ze2ad(ji,jj,jk-1) = ze2ad(ji,jj,jk-1) + zc2ad * EXP( - e3t(ji,jj,jk-1) * zekg(ji,jj) )
                        !ze3ad(ji,jj,jk-1) = ze3ad(ji,jj,jk-1) + zc3ad * EXP( - e3t(ji,jj,jk-1) * zekr(ji,jj) )
                        !zc0ad = 0.0_wp
                        !zc1ad = 0.0_wp
                        !zc2ad = 0.0_wp
                        !zc3ad = 0.0_wp
                     !END DO
                  !END DO
               !END DO
               !!
               !qsr_ad(:,:) = qsr_ad(:,:) + zeaad(:,:,1)
               !qsr_ad(:,:) = qsr_ad(:,:) + zcoef  * ze3ad(:,:,1)
               !qsr_ad(:,:) = qsr_ad(:,:) + zcoef  * ze2ad(:,:,1)
               !qsr_ad(:,:) = qsr_ad(:,:) + zcoef  * ze1ad(:,:,1)
               !qsr_ad(:,:) = qsr_ad(:,:) + rn_abs * ze0ad(:,:,1)
               !!
               CALL ctl_stop('tra_qsr_adj: key_vvl or chlorophyll management not implemented in TAM yet')
            ELSE                                                 !*  Constant Chlorophyll
               DO jk = 1, nksr
                  etot3_ad(:,:,jk) = etot3_ad(:,:,jk) + qsr_hc_ad(:,:,jk) * qsr(:,:)
                  qsr_ad(  :,:   ) = qsr_ad(:,:)      + qsr_hc_ad(:,:,jk) * etot3(:,:,jk)
                  qsr_hc_ad(:,:,jk) = 0._wp
               END DO
            ENDIF
         ENDIF
      ENDIF
      IF ( kt /= nit000 ) THEN
         qsr_hc_ad(:,:,:) = qsr_hc_ad(:,:,:) + qsr_hc_b_ad(:,:,:)
      ENDIF
      qsr_hc_b_ad(:,:,:) = 0._wp

      CALL wrk_dealloc( jpi, jpj,      zekb, zekg, zekr        )
      CALL wrk_dealloc( jpi, jpj, jpk, ze0ad, ze1ad, ze2ad, ze3ad, zeaad )

      IF( nn_timing == 1 )  CALL timing_stop('tra_qsr_adj')

   END SUBROUTINE tra_qsr_adj
   SUBROUTINE tra_qsr_adj_tst ( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!          ***  ROUTINE tra_sbc_adj_tst : TEST OF tra_sbc_adj  ***
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
      !!        ! 08-08 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER ::  &
         & jstp,  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zta_tlin,     &! Tangent input : after temperature
         & zta_tlout,    &! Tangent output: after temperature
         & zta_adout,    &! Adjoint output: after temperature
         & zta_adin,     &! Adjoint input : after temperature
         & zqsr_hc_tlin,     &! qsr_hcngent input : after temperature
         & zqsr_hc_tlout,    &! qsr_hcngent output: after temperature
         & zqsr_hc_adout,    &! Adjoint output: after temperature
         & zqsr_hc_adin,     &! Adjoint input : after temperature
         & zqsr_hc_b_tlout,    &! qsr_hc_bngent output: after temperature
         & zqsr_hc_b_adin,     &! Adjoint input : after temperature
         & zetot3_tlin,  &! Tangent input
         & zetot3_adout, &! Adjoint output
         & zta,          & ! temporary after temperature
         & zqsr_hc,          & ! temporary after temperature
         & zqsr_hc_b,          & ! temporary after temperature
         & zetot3          ! temporary
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zqsr_tlin,    &! Tangent input : solar radiation (w/m2)
         & zqsr_adout,   &! Adjoint output: solar radiation (w/m2)
         & zqsr           ! temporary solar radiation (w/m2)
      REAL(KIND=wp) ::       &
         & zsp1,             & ! scalar product involving the tangent routine
         & zsp2,             & ! scalar product involving the adjoint routine
         & zsp2_1,           & ! scalar product involving the adjoint routine
         & zsp2_2,           & ! scalar product involving the adjoint routine
         & zsp2_3              ! scalar product involving the adjoint routine
      CHARACTER(LEN=14) :: &
         & cl_name

      ALLOCATE( &
         & zta_tlin(jpi,jpj,jpk),     &
         & zta_tlout(jpi,jpj,jpk),    &
         & zta_adout(jpi,jpj,jpk),    &
         & zta_adin(jpi,jpj,jpk),     &
         & zta(jpi,jpj,jpk),          &
         & zqsr_hc_tlin(jpi,jpj,jpk),     &
         & zqsr_hc_tlout(jpi,jpj,jpk),    &
         & zqsr_hc_adout(jpi,jpj,jpk),    &
         & zqsr_hc_adin(jpi,jpj,jpk),     &
         & zqsr_hc(jpi,jpj,jpk),          &
         & zqsr_hc_b_tlout(jpi,jpj,jpk),    &
         & zqsr_hc_b_adin(jpi,jpj,jpk),     &
         & zqsr_hc_b(jpi,jpj,jpk),          &
         & zqsr_tlin(jpi,jpj),        &
         & zqsr_adout(jpi,jpj),       &
         & zetot3_tlin(jpi,jpj,jpk),  &
         & zetot3_adout(jpi,jpj,jpk), &
         & zqsr(jpi,jpj),             &
         & zetot3(jpi,jpj,jpk)        &
         & )
      ! Initialize the reference state
      qsr(:,:) = 1.0_wp ! ???
      !Initialize etot3 to non-zero value until traj(nit000-1) is fixed
      etot3(:,:,1) = 2.e-8   ; etot3(:,:,2) = 1.5e-9; etot3(:,:,3) = 8.5e-10
      etot3(:,:,4) = 5.4e-10 ; etot3(:,:,5) = 3.5e-10; etot3(:,:,6:jpk) = 0.0_wp
      ! Initialize random field standard deviations
      !=============================================================
      ! 1) dx = ( T ) and dy = ( T )
      !=============================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      zta_tlin(:,:,:)        = 0.0_wp
      zta_tlout(:,:,:)       = 0.0_wp
      zta_adout(:,:,:)       = 0.0_wp
      zta_adin(:,:,:)        = 0.0_wp
      zqsr_hc_tlin(:,:,:)    = 0.0_wp
      zqsr_hc_tlout(:,:,:)   = 0.0_wp
      zqsr_hc_adout(:,:,:)   = 0.0_wp
      zqsr_hc_adin(:,:,:)    = 0.0_wp
      zqsr_hc_b_tlout(:,:,:) = 0.0_wp
      zqsr_hc_b_adin(:,:,:)  = 0.0_wp
      zqsr_adout(:,:)        = 0.0_wp
      zqsr_tlin(:,:)         = 0.0_wp
      zetot3_tlin(:,:,:)     = 0.0_wp
      zetot3_adout(:,:,:)    = 0.0_wp
      tsa_ad(:,:,:,jp_tem)   = 0.0_wp
      qsr_ad(:,:)            = 0.0_wp
      qsr_hc_ad(:,:,:)       = 0.0_wp
      qsr_hc_b_ad(:,:,:)     = 0.0_wp
      etot3_ad(:,:,:)        = 0.0_wp

      CALL grid_random(  zqsr   , 'T', 0.0_wp, stdqsr )
      CALL grid_random(  zqsr_hc, 'T', 0.0_wp, stdqsr )
      CALL grid_random(  zta    , 'T', 0.0_wp, stdt )
      CALL grid_random(  zetot3 , 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zta_tlin(ji,jj,jk) = zta(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zqsr_hc_tlin(ji,jj,jk) = zqsr_hc(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zetot3_tlin(ji,jj,jk) = zetot3(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zqsr_tlin(ji,jj)  = zqsr(ji,jj)
         END DO
      END DO
      ! Test for time steps nit000 and nit000 + 1 (the matrix changes)
      DO jstp = nit000, nit000 + 1
         !--------------------------------------------------------------------
         ! Call the tangent routine: dy = L dx
         !--------------------------------------------------------------------

         tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
         etot3_tl(:,:,:)      = zetot3_tlin(:,:,:)
         qsr_tl(:,:)          = zqsr_tlin(:,:)
         qsr_hc_tl(:,:,:)     = zqsr_hc_tlin(:,:,:)

         CALL tra_qsr_tan( jstp )

         zta_tlout(:,:,:)       = tsa_tl(:,:,:,jp_tem)
         zqsr_hc_tlout(:,:,:)   = qsr_hc_tl(:,:,:)
         zqsr_hc_b_tlout(:,:,:) = qsr_hc_b_tl(:,:,:)

         !--------------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !--------------------------------------------------------------------

         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zta_adin(ji,jj,jk) = zta_tlout(ji,jj,jk) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                     &               * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zqsr_hc_adin(ji,jj,jk) = zqsr_hc_tlout(ji,jj,jk) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                     &               * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zqsr_hc_b_adin(ji,jj,jk) = zqsr_hc_b_tlout(ji,jj,jk) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                     &               * tmask(ji,jj,jk)
               END DO
            END DO
         END DO

         !--------------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !--------------------------------------------------------------------

         zsp1 = DOT_PRODUCT( zta_tlout, zta_adin ) &
            &  + DOT_PRODUCT( zqsr_hc_tlout, zqsr_hc_adin ) &
            &  + DOT_PRODUCT( zqsr_hc_b_tlout, zqsr_hc_b_adin )

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------

         etot3_ad(:,:,:)      = 0.0_wp
         qsr_ad(:,:)          = 0.0_wp
         tsa_ad(:,:,:,jp_tem) = zta_adin(:,:,:)
         qsr_hc_ad(:,:,:)     = zqsr_hc_adin(:,:,:)
         qsr_hc_b_ad(:,:,:)   = zqsr_hc_b_adin(:,:,:)

         CALL tra_qsr_adj( jstp )

         zta_adout(:,:,:)       = tsa_ad(:,:,:,jp_tem)
         zetot3_adout(:,:,:)    = etot3_ad(:,:,:)
         zqsr_adout(:,:)        = qsr_ad(:,:)
         zqsr_hc_adout(:,:,:)   = qsr_hc_ad(:,:,:)

         !--------------------------------------------------------------------
         ! Compute the scalar product: dx^T L^T W dy
         !--------------------------------------------------------------------

         zsp2_1 = DOT_PRODUCT( zta_tlin    , zta_adout     )
         zsp2_1 = zsp2_1 + DOT_PRODUCT( zqsr_hc_tlin    , zqsr_hc_adout   )
         zsp2_2 = DOT_PRODUCT( zqsr_tlin   , zqsr_adout    )
         zsp2_3 = DOT_PRODUCT( zetot3_tlin , zetot3_adout  )

         zsp2 = zsp2_1 + zsp2_2 + zsp2_3

         ! Compare the scalar products

         ! 14 char:   '12345678901234'
         IF (jstp == nit000) THEN
            cl_name = 'tra_qsr_adj  1'
         ELSE
            cl_name = 'tra_qsr_adj  2'
         END IF
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )
      END DO

      DEALLOCATE( &
         & zta_tlin,        &
         & zta_tlout,       &
         & zta_adout,       &
         & zta_adin,        &
         & zta,             &
         & zqsr_hc_tlin,    &
         & zqsr_hc_tlout,   &
         & zqsr_hc_adout,   &
         & zqsr_hc_adin,    &
         & zqsr_hc,         &
         & zqsr_hc_b_tlout, &
         & zqsr_hc_b_adin,  &
         & zqsr_hc_b,       &
         & zqsr_adout,      &
         & zqsr_tlin,       &
         & zqsr             &
         & )

      !
   END SUBROUTINE tra_qsr_adj_tst
   SUBROUTINE tra_qsr_init_tam
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_qsr_init_tan  ***
      !!
      !! ** Purpose :   Initialization for the penetrative solar radiation
      !!
      !! ** Method  :   The profile of solar radiation within the ocean is set
      !!      from two length scale of penetration (rn_si0,rn_si1) and a ratio
      !!      (rn_abs). These parameters are read in the namtra_qsr namelist. The
      !!      default values correspond to clear water (type I in Jerlov'
      !!      (1968) classification.
      !!         called by tra_qsr at the first timestep (nit000)
      !!
      !! ** Action  : - initialize rn_si0, rn_si1 and rn_abs
      !!
      !! Reference : Jerlov, N. G., 1968 Optical Oceanography, Elsevier, 194pp.
      !!----------------------------------------------------------------------

      IF( ln_traqsr  ) THEN      !  Initialisation of Light Penetration  !
         !                       ! ===================================== !
         !
         !                                ! ---------------------------------- !
         IF( ln_qsr_rgb ) THEN            !  Red-Green-Blue light penetration  !
            !                             ! ---------------------------------- !
            etot3_tl(:,:,:) = 0.0_wp
            etot3_ad(:,:,:) = 0.0_wp
            !
         ENDIF
            !                             ! ---------------------------------- !
         IF( ln_qsr_2bd ) THEN            !    2 bands    light penetration    !
            !                             ! ---------------------------------- !
            etot3_tl(:,:,:) = 0.0_wp
            etot3_ad(:,:,:) = 0.0_wp
            !
         ENDIF
         !
      ENDIF
      !
   END SUBROUTINE tra_qsr_init_tam

   !!======================================================================
END MODULE traqsr_tam
