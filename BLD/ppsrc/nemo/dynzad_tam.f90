MODULE dynzad_tam
   !!======================================================================
   !!                       ***  MODULE  dynzad_tam  ***
   !! Ocean dynamics : vertical advection trend
   !!                        Tangent and Adjoint module
   !!======================================================================
   !! History of the direct module:
   !!            6.0  !  91-01  (G. Madec) Original code
   !!            7.0  !  91-11  (G. Madec)
   !!            7.5  !  96-01  (G. Madec) statement function for e3
   !!            8.5  !  02-07  (G. Madec) j-k-i case: Original code
   !!            8.5  !  02-07  (G. Madec) Free form, F90
   !! History of the tam module:
   !!            9.0  !  08-08  (A. Vidard) first version
   !!    NEMO    3.4  !  12-07  (P.-A. Bouttier) phasing with 3.4
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zad_tan   : tangent of the vertical advection momentum trend
   !!   dyn_zad_adj   : adjoint of the vertical advection momentum trend
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE oce_tam
   USE dom_oce
   USE in_out_manager
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zad_tan   ! routine called by step_tam.F90
   PUBLIC   dyn_zad_adj   ! routine called by step_tam.F90

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
   !!   OPA 9.0 , LOCEAN-IPSL (2005)
   !! $Header$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   SUBROUTINE dyn_zad_tan ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynzad_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now vertical momentum advection trend and
      !!      add it to the general trend of momentum equation.
      !!
      !! ** Method of the direct routine:
      !!      The now vertical advection of momentum is given by:
      !!         w dz(u) = ua + 1/(e1u*e2u*e3u) mk+1[ mi(e1t*e2t*wn) dk(un) ]
      !!         w dz(v) = va + 1/(e1v*e2v*e3v) mk+1[ mj(e1t*e2t*wn) dk(vn) ]
      !!      Add this trend to the general trend (ua,va):
      !!         (ua,va) = (ua,va) + w dz(u,v)
      !!
      !! ** Action  : - Update (ua_tl,va_tl) with the vert. momentum advection trends
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step inedx
      !!
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   zuatl, zvatl        ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zww                ! 2D  workspace
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zwwtl              ! 2D  workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwuwtl, zwvwtl     ! 3D workspace
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('dyn_zad_tan')
      !
      CALL wrk_alloc( jpi,jpj, zww )
      CALL wrk_alloc( jpi,jpj, zwwtl )
      CALL wrk_alloc( jpi,jpj,jpk, zwuwtl , zwvwtl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'dyn_zad_tan : arakawa advection scheme'
         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~'
         CALL flush(numout)
      ENDIF

      DO jk = 2, jpkm1              ! Vertical momentum advection at level w and u- and v- vertical
         DO jj = 2, jpj                   ! vertical fluxes
            DO ji = 2, jpi             ! vector opt.
               zww(  ji,jj) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn(   ji,jj,jk)
               zwwtl(ji,jj) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn_tl(ji,jj,jk)
            END DO
         END DO
         DO jj = 2, jpjm1                 ! vertical momentum advection at w-point
            DO ji = 2, jpim1        ! vector opt.
               zwuwtl(ji,jj,jk) = ( zwwtl(ji+1,jj  ) + zwwtl(ji,jj) ) * ( un(   ji,jj,jk-1)-un(   ji,jj,jk) ) + &
                  &               ( zww(  ji+1,jj  ) + zww(  ji,jj) ) * ( un_tl(ji,jj,jk-1)-un_tl(ji,jj,jk) )
               zwvwtl(ji,jj,jk) = ( zwwtl(ji  ,jj+1) + zwwtl(ji,jj) ) * ( vn(   ji,jj,jk-1)-vn(   ji,jj,jk) ) + &
                  &               ( zww(  ji  ,jj+1) + zww(  ji,jj) ) * ( vn_tl(ji,jj,jk-1)-vn_tl(ji,jj,jk) )
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1              ! Surface and bottom values set to zero
         DO ji = 2, jpim1           ! vector opt.
            zwuwtl(ji,jj, 1 ) = 0.0_wp
            zwvwtl(ji,jj, 1 ) = 0.0_wp
            zwuwtl(ji,jj,jpk) = 0.0_wp
            zwvwtl(ji,jj,jpk) = 0.0_wp
         END DO
      END DO

      DO jk = 1, jpkm1              ! Vertical momentum advection at u- and v-points
         DO jj = 2, jpjm1
            DO ji = 2, jpim1       ! vector opt.
               !                         ! vertical momentum advective trends
               zuatl = - ( zwuwtl(ji,jj,jk) + zwuwtl(ji,jj,jk+1) ) / ( e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) )
               zvatl = - ( zwvwtl(ji,jj,jk) + zwvwtl(ji,jj,jk+1) ) / ( e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) )
               !                         ! add the trends to the general momentum trends
               ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) + zuatl
               va_tl(ji,jj,jk) = va_tl(ji,jj,jk) + zvatl
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj, zww )
      CALL wrk_dealloc( jpi,jpj, zwwtl )
      CALL wrk_dealloc( jpi,jpj,jpk, zwuwtl , zwvwtl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zad_tan')
      !
   END SUBROUTINE dyn_zad_tan

   SUBROUTINE dyn_zad_adj ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dynzad_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now vertical momentum advection trend and
      !!      add it to the general trend of momentum equation.
      !!
      !! ** Method of the direct routine:
      !!      The now vertical advection of momentum is given by:
      !!         w dz(u) = ua + 1/(e1u*e2u*e3u) mk+1[ mi(e1t*e2t*wn) dk(un) ]
      !!         w dz(v) = va + 1/(e1v*e2v*e3v) mk+1[ mj(e1t*e2t*wn) dk(vn) ]
      !!      Add this trend to the general trend (ua,va):
      !!         (ua,va) = (ua,va) + w dz(u,v)
      !!
      !! ** Action  : - Update (ua_tl,va_tl) with the vert. momentum advection trends
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step inedx
      !!
      INTEGER  ::   ji, jj, jk      ! dummy loop indices
      REAL(wp) ::   zuaad, zvaad        ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zww                ! 2D  workspace
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zwwad              ! 2D  workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwuwad, zwvwad     ! 3D workspace
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('dyn_zad_adj')
      !
      CALL wrk_alloc( jpi,jpj, zww )
      CALL wrk_alloc( jpi,jpj, zwwad )
      CALL wrk_alloc( jpi,jpj,jpk, zwuwad , zwvwad )
      !
      IF( kt == nitend ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'dyn_zad_adj : arakawa advection scheme'
         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      zuaad         = 0.0_wp ; zvaad         = 0.0_wp
      zwuwad(:,:,:) = 0.0_wp ; zwvwad(:,:,:) = 0.0_wp ; zwwad(:,:)    = 0.0_wp

      DO jk = jpkm1, 1, -1             ! Vertical momentum advection at u- and v-points
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1 ! vector opt.
               !                       ! add the trends to the general momentum trends
               zuaad = zuaad + ua_ad(ji,jj,jk)
               zvaad = zvaad + va_ad(ji,jj,jk)
               !                       ! vertical momentum advective trends
               zwuwad(ji,jj,jk  ) = zwuwad(ji,jj,jk  ) - zuaad / ( e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) )
               zwuwad(ji,jj,jk+1) = zwuwad(ji,jj,jk+1) - zuaad / ( e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) )
               zuaad = 0.0_wp

               zwvwad(ji,jj,jk  ) = zwvwad(ji,jj,jk  ) - zvaad / ( e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) )
               zwvwad(ji,jj,jk+1) = zwvwad(ji,jj,jk+1) - zvaad / ( e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) )
               zvaad = 0.0_wp
            END DO
         END DO
      END DO
      DO jj = 2, jpjm1              ! Surface and bottom values set to zero
         DO ji = 2, jpim1           ! vector opt.
            zwuwad(ji,jj, 1 ) = 0.0_wp
            zwvwad(ji,jj, 1 ) = 0.0_wp
            zwuwad(ji,jj,jpk) = 0.0_wp
            zwvwad(ji,jj,jpk) = 0.0_wp
         END DO
      END DO
      DO jk = jpkm1, 2, -1             ! Vertical momentum advection at level w and u- and v- vertical
         DO jj = 2, jpj                   ! vertical fluxes
            DO ji = 2, jpi             ! vector opt.
               zww(ji,jj) = 0.25 * e1t(ji,jj) * e2t(ji,jj) * wn(ji,jj,jk)
            END DO
         END DO
         DO jj = jpjm1, 2, -1          ! vertical momentum advection at w-point
            DO ji = jpim1, 2, -1 ! vector opt.
               zwwad(ji,jj+1) = zwwad(ji,jj+1) + zwvwad(ji,jj,jk) * ( vn(ji,jj,jk-1)-vn(ji,jj,jk) )
               zwwad(ji,jj  ) = zwwad(ji,jj  ) + zwvwad(ji,jj,jk) * ( vn(ji,jj,jk-1)-vn(ji,jj,jk) )
               vn_ad(ji,jj,jk-1) = vn_ad(ji,jj,jk-1) + zwvwad(ji,jj,jk) * ( zww(ji,jj+1) + zww(ji,jj) )
               vn_ad(ji,jj,jk  ) = vn_ad(ji,jj,jk  ) - zwvwad(ji,jj,jk) * ( zww(ji,jj+1) + zww(ji,jj) )
               zwvwad(ji,jj,jk) = 0.0_wp

               zwwad(ji+1,jj) = zwwad(ji+1,jj) + zwuwad(ji,jj,jk) * ( un(ji,jj,jk-1)-un(ji,jj,jk) )
               zwwad(ji  ,jj) = zwwad(ji  ,jj) + zwuwad(ji,jj,jk) * ( un(ji,jj,jk-1)-un(ji,jj,jk) )
               un_ad(ji,jj,jk-1) = un_ad(ji,jj,jk-1)  + zwuwad(ji,jj,jk) * ( zww(ji+1,jj) + zww(ji,jj) )
               un_ad(ji,jj,jk  ) = un_ad(ji,jj,jk  )  - zwuwad(ji,jj,jk) * ( zww(ji+1,jj) + zww(ji,jj) )
               zwuwad(ji,jj,jk) = 0.0_wp
            END DO
         END DO
         DO jj = jpj, 2, -1                   ! vertical fluxes
            DO ji = jpi, 2, -1             ! vector opt.
               wn_ad(ji,jj,jk) = wn_ad(ji,jj,jk) + zwwad(ji,jj) * 0.25 * e1t(ji,jj) * e2t(ji,jj)
               zwwad(ji,jj) = 0.0_wp
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj, zww )
      CALL wrk_dealloc( jpi,jpj, zwwad )
      CALL wrk_dealloc( jpi,jpj,jpk, zwuwad, zwvwad )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zad_adj')
      !
   END SUBROUTINE dyn_zad_adj
   SUBROUTINE dyn_zad_adj_tst( kumadt )
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit
      ! done in dynadv_tam
   END SUBROUTINE dyn_zad_adj_tst
END MODULE dynzad_tam
