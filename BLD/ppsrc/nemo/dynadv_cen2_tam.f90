MODULE dynadv_cen2_tam
   !!======================================================================
   !!                       ***  MODULE  dynadv_cen2_tam  ***
   !! Ocean dynamics: Update the momentum trend with the flux form advection
   !!                 using a 2nd order centred scheme
   !!======================================================================
   !! History of the direct module:
   !!            2.0  ! 2006-08  (G. Madec, S. Theetten)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !! History ot the T&A module
   !!            3.2  ! 2011-01  (A. Vidard) Original version
   !!            3.4  ! 2012-07  (P.-A. bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   dyn_adv_cen2       : flux form momentum advection (ln_dynadv_cen2=T)
   !!                        trends using a 2nd order centred scheme
   !!----------------------------------------------------------------------
   USE oce
   USE dom_oce
   USE oce_tam
   USE in_out_manager
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_adv_cen2_tan                 ! routine called by step_tam.F90

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
   SUBROUTINE dyn_adv_cen2_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_cen2_tan  ***
      !!
      !! ** Purpose :   Compute the now momentum advection trend in flux form
      !!      and the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zbu, zbv     ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zfu_ttl, zfu_ftl, zfu_uwtl   ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zfv_ttl, zfv_ftl, zfv_vwtl   !  -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zfw, zfu, zfv          !  -      -
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zfwtl, zfutl, zfvtl          !  -      -
      !!----------------------------------------------------------------------
      !
      IF ( nn_timing == 1 )  CALL timing_start('dyn_adv_cen2_tan')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zfu_ttl, zfv_ttl, zfu_ftl, zfv_ftl, zfu_uwtl, zfv_vwtl, zfutl, zfvtl, zfwtl )
      !
      IF ( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_adv_cen2_tan : 2nd order flux form momentum advection'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~'
      ENDIF
      !                                      ! ====================== !
      !                                      !  Horizontal advection  !
      DO jk = 1, jpkm1                       ! ====================== !
         !                                         ! horizontal volume fluxes
         zfu(:,:,jk)   = 0.25 * e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
         zfv(:,:,jk)   = 0.25 * e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)
         zfutl(:,:,jk) = 0.25 * e2u(:,:) * e3u(:,:,jk) * un_tl(:,:,jk)
         zfvtl(:,:,jk) = 0.25 * e1v(:,:) * e3v(:,:,jk) * vn_tl(:,:,jk)
         !
         DO jj = 1, jpjm1                          ! horizontal momentum fluxes at T- and F-point
            DO ji = 1, jpim1   ! vector opt.
               zfu_ttl(ji+1,jj  ,jk) = ( zfutl(ji,jj,jk) + zfutl(ji+1,jj  ,jk) ) * ( un(ji,jj,jk) + un(ji+1,jj  ,jk) ) &
                  &                  + ( zfu(ji,jj,jk) + zfu(ji+1,jj  ,jk) ) * ( un_tl(ji,jj,jk) + un_tl(ji+1,jj  ,jk) )
               zfv_ftl(ji  ,jj  ,jk) = ( zfvtl(ji,jj,jk) + zfvtl(ji+1,jj  ,jk) ) * ( un(ji,jj,jk) + un(ji  ,jj+1,jk) ) &
                  &                  + ( zfv(ji,jj,jk) + zfv(ji+1,jj  ,jk) ) * ( un_tl(ji,jj,jk) + un_tl(ji  ,jj+1,jk) )
               zfu_ftl(ji  ,jj  ,jk) = ( zfutl(ji,jj,jk) + zfutl(ji  ,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji+1,jj  ,jk) ) &
                  &                  + ( zfu(ji,jj,jk) + zfu(ji  ,jj+1,jk) ) * ( vn_tl(ji,jj,jk) + vn_tl(ji+1,jj  ,jk) )
               zfv_ttl(ji  ,jj+1,jk) = ( zfvtl(ji,jj,jk) + zfvtl(ji  ,jj+1,jk) ) * ( vn(ji,jj,jk) + vn(ji  ,jj+1,jk) ) &
                  &                  + ( zfv(ji,jj,jk) + zfv(ji  ,jj+1,jk) ) * ( vn_tl(ji,jj,jk) + vn_tl(ji  ,jj+1,jk) )
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
      !
      !                                      ! ==================== !
      !                                      !  Vertical advection  !
      DO jk = 1, jpkm1                       ! ==================== !
         !                                         ! Vertical volume fluxes 
         zfw(:,:,jk)   = 0.25 * e1t(:,:) * e2t(:,:) * wn(:,:,jk)
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
                     zfu_uwtl(ji,jj, 1 ) = 2.e0 * ( zfwtl(ji,jj,1) + zfwtl(ji+1,jj  ,1) ) * un(   ji,jj,1) &
                        &                + 2.e0 * ( zfw(  ji,jj,1) + zfw(  ji+1,jj  ,1) ) * un_tl(ji,jj,1)
                     zfv_vwtl(ji,jj, 1 ) = 2.e0 * ( zfwtl(ji,jj,1) + zfwtl(ji  ,jj+1,1) ) * vn(   ji,jj,1) &
                        &                + 2.e0 * ( zfw(  ji,jj,1) + zfw(  ji  ,jj+1,1) ) * vn_tl(ji,jj,1)
                  END DO
               END DO
            ENDIF
         ELSE                                      ! interior fluxes
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zfu_uwtl(ji,jj,jk) = ( zfwtl(ji,jj,jk)+ zfwtl(ji+1,jj  ,jk) ) * ( un(   ji,jj,jk) + un(   ji,jj,jk-1) ) &
                     &               + ( zfw(  ji,jj,jk)+ zfw(  ji+1,jj  ,jk) ) * ( un_tl(ji,jj,jk) + un_tl(ji,jj,jk-1) )
                  zfv_vwtl(ji,jj,jk) = ( zfwtl(ji,jj,jk)+ zfwtl(ji  ,jj+1,jk) ) * ( vn(   ji,jj,jk) + vn(   ji,jj,jk-1) ) &
                     &               + ( zfw(  ji,jj,jk)+ zfw(  ji  ,jj+1,jk) ) * ( vn_tl(ji,jj,jk) + vn_tl(ji,jj,jk-1) )
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
      CALL wrk_dealloc( jpi, jpj, jpk, zfu_ttl, zfv_ttl, zfu_ftl, zfv_ftl, zfu_uwtl, zfv_vwtl, zfutl, zfvtl, zfwtl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_adv_cen2_tan')
      !
   END SUBROUTINE dyn_adv_cen2_tan
   !!==============================================================================
END MODULE dynadv_cen2_tam
