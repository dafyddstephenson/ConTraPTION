MODULE dynadv_ubs_tam
   !!======================================================================
   !!                       ***  MODULE  dynadv_ubs_tam  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 trend using a 3rd order upstream biased scheme
   !!======================================================================
   !! History of the direct module :
   !!            2.0  ! 2006-08  (R. Benshila, L. Debreu)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !! History of the T&A module :
   !!            3.2  ! 2011-02  (A. Vidard)  Original
   !!            3.4  ! 2012-0è  (P.-A. Bouttier)  Phasing with 3.4
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv_ubs   : flux form momentum advection using    (ln_dynadv=T)
   !!                   an 3rd order Upstream Biased Scheme or Quick scheme
   !!                   combined with 2nd or 4th order finite differences
   !!----------------------------------------------------------------------
   USE oce_tam        ! ocean dynamics and tracers
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE prtctl         ! Print control
   USE lib_mpp        ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   REAL(wp), PARAMETER :: gamma1 = 1._wp/3._wp  ! =1/4 quick      ; =1/3  3rd order UBS
   REAL(wp), PARAMETER :: gamma2 = 0._wp  ! =0   2nd order  ; =1/8  4th order centred

   PUBLIC   dyn_adv_ubs_tan   ! routine called by step.F90

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
   !! NEMO/OPA 3.2 , LODYC-IPSL  (2009)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_adv_ubs_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_ubs_tan  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!              and the general trend of the momentum equation.
      !!
      !! ** Method  :   The scheme is the one implemeted in ROMS. It depends
      !!      on two parameter gamma1 and gamma2. The former control the
      !!      upstream baised part of the scheme and the later the centred
      !!      part:     gamma1 = 0    pure centered  (no diffusive part)
      !!                       = 1/4  Quick scheme
      !!                       = 1/3  3rd order Upstream biased scheme
      !!                gamma2 = 0    2nd order finite differencing
      !!                       = 1/8  4th order finite differencing
      !!      For stability reasons, the first term of the fluxes which cor-
      !!      responds to a second order centered scheme is evaluated using
      !!      the now velocity (centered in time) while the second term which
      !!      is the diffusive part of the scheme, is evaluated using the
      !!      before velocity (forward in time).
      !!      Default value (hard coded in the begining of the module) are
      !!      gamma1=1/4 and gamma2=1/8.
      !!
      !! ** Action : - (ua,va) updated with the 3D advective momentum trends
      !!
      !! Reference : Shchepetkin & McWilliams, 2005, Ocean Modelling.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt     ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      REAL(wp) ::   zbu, zbv    ! temporary scalars
      REAL(wp) ::   zui, zvj, zfuj, zfvi, zl_u, zl_v   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfu, zfv     ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfu_t, zfu_f     ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfv_t, zfv_f     !    "           "
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfw, zfu_uw, zfv_vw
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zlu_uu, zlu_uv   ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zlv_vv, zlv_vu   ! temporary workspace
      REAL(wp) ::   zuitl, zvjtl, zfujtl, zfvitl, zl_utl, zl_vtl   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfutl, zfvtl     ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfu_ttl, zfu_ftl     ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfv_ttl, zfv_ftl     !    "           "
      REAL(wp), POINTER, DIMENSION(:,:,:)   ::   zfwtl, zfu_uwtl, zfv_vwtl
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zlu_uutl, zlu_uvtl   ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:,:) ::   zlv_vvtl, zlv_vutl   ! temporary workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_adv_ubs_tan')
      !
      CALL wrk_alloc( jpi, jpj, jpk,       zfu_t , zfv_t , zfu_f , zfv_f, zfu_uw, zfv_vw, zfu, zfv, zfw )
      CALL wrk_alloc( jpi, jpj, jpk, jpts, zlu_uu, zlv_vv, zlu_uv, zlv_vu                               )
      CALL wrk_alloc( jpi, jpj, jpk,       zfu_ttl , zfv_ttl , zfu_ftl , zfv_ftl, zfu_uwtl, zfv_vwtl, zfutl, zfvtl, zfwtl )
      CALL wrk_alloc( jpi, jpj, jpk, jpts, zlu_uutl, zlv_vvtl, zlu_uvtl, zlv_vutl                                         )
      !
      IF( kt == nit000) THEN
      !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_adv_ubs_tan : UBS flux form momentum advection'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      zfu_t(:,:,:) = 0.e0
      zfv_t(:,:,:) = 0.e0
      zfu_f(:,:,:) = 0.e0
      zfv_f(:,:,:) = 0.e0
      !
      zlu_uu(:,:,:,:) = 0.e0
      zlv_vv(:,:,:,:) = 0.e0
      zlu_uv(:,:,:,:) = 0.e0
      zlv_vu(:,:,:,:) = 0.e0

      zfu_ttl(:,:,:) = 0.e0
      zfv_ttl(:,:,:) = 0.e0
      zfu_ftl(:,:,:) = 0.e0
      zfv_ftl(:,:,:) = 0.e0
      !
      zlu_uutl(:,:,:,:) = 0.e0
      zlv_vvtl(:,:,:,:) = 0.e0
      zlu_uvtl(:,:,:,:) = 0.e0
      zlv_vutl(:,:,:,:) = 0.e0

      !                                      ! =========================== !
      DO jk = 1, jpkm1                       !  Laplacian of the velocity  !
         !                                   ! =========================== !
         !                                         ! horizontal volume fluxes
         zfu(  :,:,jk) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
         zfv(  :,:,jk) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)
         zfutl(:,:,jk) = e2u(:,:) * e3u(:,:,jk) * un_tl(:,:,jk)
         zfvtl(:,:,jk) = e1v(:,:) * e3v(:,:,jk) * vn_tl(:,:,jk)
         !
         DO jj = 2, jpjm1                          ! laplacian
            DO ji = 2, jpim1   ! vector opt.
               zlu_uu(ji,jj,jk,1) = ( ub (ji+1,jj,jk)-2.*ub (ji,jj,jk)+ub (ji-1,jj,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,1) = ( vb (ji,jj+1,jk)-2.*vb (ji,jj,jk)+vb (ji,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,1) = ( ub (ji,jj+1,jk)-2.*ub (ji,jj,jk)+ub (ji,jj-1,jk) ) * umask(ji,jj,jk)
               zlv_vu(ji,jj,jk,1) = ( vb (ji+1,jj,jk)-2.*vb (ji,jj,jk)+vb (ji-1,jj,jk) ) * vmask(ji,jj,jk)

               zlu_uu(ji,jj,jk,2) = ( zfu(ji+1,jj,jk)-2.*zfu(ji,jj,jk)+zfu(ji-1,jj,jk) ) * umask(ji,jj,jk)
               zlv_vv(ji,jj,jk,2) = ( zfv(ji,jj+1,jk)-2.*zfv(ji,jj,jk)+zfv(ji,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uv(ji,jj,jk,2) = ( zfu(ji,jj+1,jk)-2.*zfu(ji,jj,jk)+zfu(ji,jj-1,jk) ) * umask(ji,jj,jk)
               zlv_vu(ji,jj,jk,2) = ( zfv(ji+1,jj,jk)-2.*zfv(ji,jj,jk)+zfv(ji-1,jj,jk) ) * vmask(ji,jj,jk)

               zlu_uutl(ji,jj,jk,1) = ( ub_tl (ji+1,jj,jk)-2.*ub_tl (ji,jj,jk)+ub_tl (ji-1,jj,jk) ) * umask(ji,jj,jk)
               zlv_vvtl(ji,jj,jk,1) = ( vb_tl (ji,jj+1,jk)-2.*vb_tl (ji,jj,jk)+vb_tl (ji,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uvtl(ji,jj,jk,1) = ( ub_tl (ji,jj+1,jk)-2.*ub_tl (ji,jj,jk)+ub_tl (ji,jj-1,jk) ) * umask(ji,jj,jk)
               zlv_vutl(ji,jj,jk,1) = ( vb_tl (ji+1,jj,jk)-2.*vb_tl (ji,jj,jk)+vb_tl (ji-1,jj,jk) ) * vmask(ji,jj,jk)

               zlu_uutl(ji,jj,jk,2) = ( zfutl(ji+1,jj,jk)-2.*zfutl(ji,jj,jk)+zfutl(ji-1,jj,jk) ) * umask(ji,jj,jk)
               zlv_vvtl(ji,jj,jk,2) = ( zfvtl(ji,jj+1,jk)-2.*zfvtl(ji,jj,jk)+zfvtl(ji,jj-1,jk) ) * vmask(ji,jj,jk)
               zlu_uvtl(ji,jj,jk,2) = ( zfutl(ji,jj+1,jk)-2.*zfutl(ji,jj,jk)+zfutl(ji,jj-1,jk) ) * umask(ji,jj,jk)
               zlv_vutl(ji,jj,jk,2) = ( zfvtl(ji+1,jj,jk)-2.*zfvtl(ji,jj,jk)+zfvtl(ji-1,jj,jk) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
!!!gm BUG !!!  just below this should be +1 in all the communications
      !CALL lbc_lnk( zlu_uu(:,:,:,1), 'U', -1.)   ;   CALL lbc_lnk( zlu_uv(:,:,:,1), 'U', -1.)
      !CALL lbc_lnk( zlu_uu(:,:,:,2), 'U', -1.)   ;   CALL lbc_lnk( zlu_uv(:,:,:,2), 'U', -1.)
      !CALL lbc_lnk( zlv_vv(:,:,:,1), 'V', -1.)   ;   CALL lbc_lnk( zlv_vu(:,:,:,1), 'V', -1.)
      !CALL lbc_lnk( zlv_vv(:,:,:,2), 'V', -1.)   ;   CALL lbc_lnk( zlv_vu(:,:,:,2), 'V', -1.)
!!gm corrected:
      CALL lbc_lnk( zlu_uu(:,:,:,1), 'U', 1. )   ;   CALL lbc_lnk( zlu_uv(:,:,:,1), 'U', 1. )
      CALL lbc_lnk( zlu_uu(:,:,:,2), 'U', 1. )   ;   CALL lbc_lnk( zlu_uv(:,:,:,2), 'U', 1. )
      CALL lbc_lnk( zlv_vv(:,:,:,1), 'V', 1. )   ;   CALL lbc_lnk( zlv_vu(:,:,:,1), 'V', 1. )
      CALL lbc_lnk( zlv_vv(:,:,:,2), 'V', 1. )   ;   CALL lbc_lnk( zlv_vu(:,:,:,2), 'V', 1. )
!!gm end
!!gm BUG !!!  just below this should be +1 in all the communications
      !CALL lbc_lnk( zlu_uutl(:,:,:,1), 'U', -1.)   ;   CALL lbc_lnk( zlu_uvtl(:,:,:,1), 'U', -1.)
      !CALL lbc_lnk( zlu_uutl(:,:,:,2), 'U', -1.)   ;   CALL lbc_lnk( zlu_uvtl(:,:,:,2), 'U', -1.)
      !CALL lbc_lnk( zlv_vvtl(:,:,:,1), 'V', -1.)   ;   CALL lbc_lnk( zlv_vutl(:,:,:,1), 'V', -1.)
      !CALL lbc_lnk( zlv_vvtl(:,:,:,2), 'V', -1.)   ;   CALL lbc_lnk( zlv_vutl(:,:,:,2), 'V', -1.)
!!gm corrected:
      CALL lbc_lnk( zlu_uutl(:,:,:,1), 'U', 1. )   ;   CALL lbc_lnk( zlu_uvtl(:,:,:,1), 'U', 1. )
      CALL lbc_lnk( zlu_uutl(:,:,:,2), 'U', 1. )   ;   CALL lbc_lnk( zlu_uvtl(:,:,:,2), 'U', 1. )
      CALL lbc_lnk( zlv_vvtl(:,:,:,1), 'V', 1. )   ;   CALL lbc_lnk( zlv_vutl(:,:,:,1), 'V', 1. )
      CALL lbc_lnk( zlv_vvtl(:,:,:,2), 'V', 1. )   ;   CALL lbc_lnk( zlv_vutl(:,:,:,2), 'V', 1. )
!!gm end

      !                                      ! ====================== !
      !                                      !  Horizontal advection  !
      DO jk = 1, jpkm1                       ! ====================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk) = 0.25 * e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk) = 0.25 * e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)
         zfutl(:,:,jk) = 0.25 * e2u(:,:) * e3u(:,:,jk) * un_tl(:,:,jk)
         zfvtl(:,:,jk) = 0.25 * e1v(:,:) * e3v(:,:,jk) * vn_tl(:,:,jk)
         !
         DO jj = 1, jpjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, jpim1   ! vector opt.
               zui = ( un(ji,jj,jk) + un(ji+1,jj  ,jk) )
               zvj = ( vn(ji,jj,jk) + vn(ji  ,jj+1,jk) )
               zuitl = ( un_tl(ji,jj,jk) + un_tl(ji+1,jj  ,jk) )
               zvjtl = ( vn_tl(ji,jj,jk) + vn_tl(ji  ,jj+1,jk) )
               !
               IF (zui > 0) THEN
                  zl_u = zlu_uu(ji  ,jj,jk,1)
                  zl_utl = zlu_uutl(ji  ,jj,jk,1)
               ELSE
                  zl_u = zlu_uu(ji+1,jj,jk,1)
                  zl_utl = zlu_uutl(ji+1,jj,jk,1)
               ENDIF
               IF (zvj > 0) THEN
                  zl_v = zlv_vv(ji,jj  ,jk,1)
                  zl_vtl = zlv_vvtl(ji,jj  ,jk,1)
               ELSE
                  zl_v = zlv_vv(ji,jj+1,jk,1)
                  zl_vtl = zlv_vvtl(ji,jj+1,jk,1)
               ENDIF
               zfu_t(ji+1,jj  ,jk) = ( zfu(ji,jj,jk) + zfu(ji+1,jj  ,jk)                               &
                  &                    - gamma2 * ( zlu_uu(ji,jj,jk,2) + zlu_uu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( zui - gamma1 * zl_u)
               zfv_t(ji  ,jj+1,jk) = ( zfv(ji,jj,jk) + zfv(ji  ,jj+1,jk)                               &
                  &                    - gamma2 * ( zlv_vv(ji,jj,jk,2) + zlv_vv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( zvj - gamma1 * zl_v)
               !
               zfu_ttl(ji+1,jj  ,jk) = ( zfutl(ji,jj,jk) + zfutl(ji+1,jj  ,jk)                               &
                  &                    - gamma2 * ( zlu_uutl(ji,jj,jk,2) + zlu_uutl(ji+1,jj  ,jk,2) )  )   &
                  &                * ( zui - gamma1 * zl_u) + ( zfu(ji,jj,jk) + zfu(ji+1,jj  ,jk)                               &
                  &                    - gamma2 * ( zlu_uu(ji,jj,jk,2) + zlu_uu(ji+1,jj  ,jk,2) )  )   &
                  &                * ( zuitl - gamma1 * zl_utl)
               zfv_ttl(ji  ,jj+1,jk) = ( zfvtl(ji,jj,jk) + zfvtl(ji  ,jj+1,jk)                               &
                  &                    - gamma2 * ( zlv_vvtl(ji,jj,jk,2) + zlv_vvtl(ji  ,jj+1,jk,2) )  )   &
                  &                * ( zvj - gamma1 * zl_v) + ( zfv(ji,jj,jk) + zfv(ji  ,jj+1,jk)                               &
                  &                    - gamma2 * ( zlv_vv(ji,jj,jk,2) + zlv_vv(ji  ,jj+1,jk,2) )  )   &
                  &                * ( zvjtl - gamma1 * zl_vtl)
               !
               zfuj = ( zfu(ji,jj,jk) + zfu(ji  ,jj+1,jk) )
               zfvi = ( zfv(ji,jj,jk) + zfv(ji+1,jj  ,jk) )
               !
               zfujtl = ( zfutl(ji,jj,jk) + zfutl(ji  ,jj+1,jk) )
               zfvitl = ( zfvtl(ji,jj,jk) + zfvtl(ji+1,jj  ,jk) )
               IF (zfuj > 0) THEN
                  zl_v = zlv_vu( ji  ,jj  ,jk,1)
                  zl_vtl = zlv_vutl( ji  ,jj  ,jk,1)
               ELSE
                  zl_v = zlv_vu( ji+1,jj,jk,1)
                  zl_vtl = zlv_vutl( ji+1,jj,jk,1)
               ENDIF
               IF (zfvi > 0) THEN
                  zl_u = zlu_uv( ji,jj  ,jk,1)
                  zl_utl = zlu_uvtl( ji,jj  ,jk,1)
               ELSE
                  zl_u = zlu_uv( ji,jj+1,jk,1)
                  zl_utl = zlu_uvtl( ji,jj+1,jk,1)
               ENDIF
               !
               zfv_ftl(ji  ,jj  ,jk) = ( zfvitl - gamma2 * ( zlv_vutl(ji,jj,jk,2) + zlv_vutl(ji+1,jj  ,jk,2) )  )   &
                  &                * ( un(ji,jj,jk) + un(ji  ,jj+1,jk) - gamma1 * zl_u )                            &
                  &                + ( zfvi - gamma2 * ( zlv_vu(ji,jj,jk,2) + zlv_vu(ji+1,jj  ,jk,2) )  )           &
                  &                * ( un_tl(ji,jj,jk) + un_tl(ji  ,jj+1,jk) - gamma1 * zl_utl )
               zfu_ftl(ji  ,jj  ,jk) = ( zfujtl - gamma2 * ( zlu_uvtl(ji,jj,jk,2) + zlu_uvtl(ji  ,jj+1,jk,2) )  )   &
                  &                * ( vn(ji,jj,jk) + vn(ji+1,jj  ,jk) - gamma1 * zl_v )                            &
                  &                + ( zfuj - gamma2 * ( zlu_uv(ji,jj,jk,2) + zlu_uv(ji  ,jj+1,jk,2) )  )           &
                  &                * ( vn_tl(ji,jj,jk) + vn_tl(ji+1,jj  ,jk) - gamma1 * zl_vtl )
            END DO
         END DO
         DO jj = 2, jpjm1                          ! divergence of horizontal momentum fluxes
            DO ji = 2, jpim1   ! vector opt.
               zbu = e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk)
               zbv = e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk)
               !
               ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) - (  zfu_ttl(ji+1,jj  ,jk) - zfu_ttl(ji  ,jj  ,jk)    &
                  &                           + zfv_ftl(ji  ,jj  ,jk) - zfv_ftl(ji  ,jj-1,jk)  ) / zbu
               va_tl(ji,jj,jk) = va_tl(ji,jj,jk) - (  zfu_ftl(ji  ,jj  ,jk) - zfu_ftl(ji-1,jj  ,jk)    &
                  &                           + zfv_ttl(ji  ,jj+1,jk) - zfv_ttl(ji  ,jj  ,jk)  ) / zbv
            END DO
         END DO
      END DO
      !                                      ! ==================== !
      !                                      !  Vertical advection  !
      DO jk = 1, jpkm1                       ! ==================== !
         !                                         ! Vertical volume fluxes 
         zfw(:,:,jk) = 0.25 * e1t(:,:) * e2t(:,:) * wn(:,:,jk)
         zfwtl(:,:,jk) = 0.25 * e1t(:,:) * e2t(:,:) * wn_tl(:,:,jk)
         !
         IF( jk == 1 ) THEN                        ! surface/bottom advective fluxes
            zfu_uwtl(:,:,jpk) = 0.e0                      ! Bottom  value : flux set to zero
            zfv_vwtl(:,:,jpk) = 0.e0
            !                                           ! Surface value :
            IF( lk_vvl ) THEN                                ! variable volume : flux set to zero
               zfu_uwtl(:,:, 1 ) = 0.e0
               zfv_vwtl(:,:, 1 ) = 0.e0
            ELSE                                             ! constant volume : advection through the surface
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zfu_uwtl(ji,jj, 1 ) = 2.e0 * ( zfwtl(ji,jj,1) + zfwtl(ji+1,jj  ,1) ) * un(ji,jj,1)   &
                                     &   + 2.e0 * ( zfw(ji,jj,1) + zfw(ji+1,jj  ,1) ) * un_tl(ji,jj,1)
                     zfv_vwtl(ji,jj, 1 ) = 2.e0 * ( zfwtl(ji,jj,1) + zfwtl(ji  ,jj+1,1) ) * vn(ji,jj,1)   &
                                     &   + 2.e0 * ( zfw(ji,jj,1) + zfw(ji  ,jj+1,1) ) * vn_tl(ji,jj,1)
                  END DO
               END DO
            ENDIF
         ELSE                                      ! interior fluxes
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zfu_uwtl(ji,jj,jk) = ( zfwtl(ji,jj,jk)+ zfwtl(ji+1,jj  ,jk) ) * ( un(ji,jj,jk) + un(ji,jj,jk-1) )   &
                                     &   + ( zfw(ji,jj,jk)+ zfw(ji+1,jj  ,jk) ) * ( un_tl(ji,jj,jk) + un_tl(ji,jj,jk-1) )
                  zfv_vwtl(ji,jj,jk) = ( zfwtl(ji,jj,jk)+ zfwtl(ji  ,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji,jj,jk-1) )   &
                                     &   + ( zfw(ji,jj,jk)+ zfw(ji  ,jj+1,jk) ) * ( vn_tl(ji,jj,jk) + vn_tl(ji,jj,jk-1) )
               END DO
            END DO
         ENDIF
      END DO
      DO jk = 1, jpkm1                             ! divergence of vertical momentum flux divergence
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ua_tl(ji,jj,jk) =  ua_tl(ji,jj,jk) - ( zfu_uwtl(ji,jj,jk) - zfu_uwtl(ji,jj,jk+1) )    &
                  &  / ( e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) )
               va_tl(ji,jj,jk) =  va_tl(ji,jj,jk) - ( zfv_vwtl(ji,jj,jk) - zfv_vwtl(ji,jj,jk+1) )    &
                  &  / ( e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) )
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk,       zfu_t , zfv_t , zfu_f , zfv_f, zfu_uw, zfv_vw, zfu, zfv, zfw )
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, zlu_uu, zlv_vv, zlu_uv, zlv_vu                               )
      CALL wrk_dealloc( jpi, jpj, jpk,       zfu_ttl , zfv_ttl , zfu_ftl , zfv_ftl, zfu_uwtl, zfv_vwtl, zfutl, zfvtl, zfwtl )
      CALL wrk_dealloc( jpi, jpj, jpk, jpts, zlu_uutl, zlv_vvtl, zlu_uvtl, zlv_vutl                               )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_adv_ubs_tan')
      !
   END SUBROUTINE dyn_adv_ubs_tan
   !!==============================================================================
END MODULE dynadv_ubs_tam
