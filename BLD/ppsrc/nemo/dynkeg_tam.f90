MODULE dynkeg_tam
   !!===========================================================================
   !!                       ***  MODULE  dynkeg_tam  ***
   !! Ocean dynamics:  kinetic energy gradient trend
   !!======================================================================
   !! History of the direct module:
   !!            1.0  !  87-09  (P. Andrich, m.-a. Foujols)  Original code
   !!            7.0  !  97-05  (G. Madec)  Split dynber into dynkeg and dynhpg
   !!            9.0  !  02-07  (G. Madec)  F90: Free form and module
   !! History of the TAM module:
   !!            9.0  !  08-08  (A. Vidard) first version
   !!      NEMO  3.4  !  12-07  (P.-A. Bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   dyn_keg_tan  : update the momentum trend with the horizontal tke
   !!   dyn_keg_adj  : update the momentum trend with the horizontal tke
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE dom_oce
   USE in_out_manager
   USE oce_tam
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_keg_tan ! routine called by step_tam module
   PUBLIC   dyn_keg_adj ! routine called by step_tam module

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

   SUBROUTINE dyn_keg_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the
      !!      general momentum trend.
      !!
      !! ** Method of the direct routine:
      !!      Compute the now horizontal kinetic energy
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend (ua,va).
      !!         ua = ua - 1/e1u di[ zhke ]
      !!         va = va - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (ua_tl, va_tl) with the hor. ke gradient trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zutl, zvtl   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zhketl   ! temporary 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_keg_tan')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zhketl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg_tan : kinetic energy gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpj         ! Horizontal kinetic energy at T-point
            DO ji = 2, jpi   ! vector opt.
               zutl = 0.5_wp * ( un_tl(ji-1,jj  ,jk) * un(ji-1,jj  ,jk)   &
                  &            + un_tl(ji  ,jj  ,jk) * un(ji  ,jj  ,jk)  )
               zvtl = 0.5_wp * ( vn_tl(ji  ,jj-1,jk) * vn(ji  ,jj-1,jk)   &
                  &            + vn_tl(ji  ,jj  ,jk) * vn(ji  ,jj  ,jk)  )
               zhketl(ji,jj,jk) = zvtl + zutl
            END DO
         END DO
         DO jj = 2, jpjm1       ! add the gradient of kinetic energy to the general momentum trends
            DO ji = 2, jpim1   ! vector opt.
               ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) - ( zhketl(ji+1,jj  ,jk) - zhketl(ji,jj,jk) ) / e1u(ji,jj)
               va_tl(ji,jj,jk) = va_tl(ji,jj,jk) - ( zhketl(ji  ,jj+1,jk) - zhketl(ji,jj,jk) ) / e2v(ji,jj)
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, jpk, zhketl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_keg_tan')
      !
   END SUBROUTINE dyn_keg_tan

   SUBROUTINE dyn_keg_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_keg_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now momentum trend due to the horizontal
      !!      gradient of the horizontal kinetic energy and add it to the
      !!      general momentum trend.
      !!
      !! ** Method of the direct routine:
      !!      Compute the now horizontal kinetic energy
      !!         zhke = 1/2 [ mi-1( un^2 ) + mj-1( vn^2 ) ]
      !!      Take its horizontal gradient and add it to the general momentum
      !!      trend (ua,va).
      !!         ua = ua - 1/e1u di[ zhke ]
      !!         va = va - 1/e2v dj[ zhke ]
      !!
      !! ** Action : - Update the (ua_ad, va_ad) with the hor. ke gradient trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zuad, zvad   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zhkead   ! temporary 3D workspace
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('dyn_keg_adj')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zhkead )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_keg_adj : kinetic energy gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      zhkead(:,:,:) = 0.0_wp ; zuad = 0.0_wp ; zvad = 0.0_wp
      !                                                ! ===============
      DO jk = jpkm1, 1, -1                             ! Horizontal slab
         !                                             ! ===============
         DO jj =  jpjm1, 2, -1 ! add the gradient of kinetic energy to the general momentum trends
            DO ji = jpim1, 2, -1 ! vector opt.
               zhkead(ji  ,jj+1,jk) = zhkead(ji  ,jj+1,jk)              &
                  &                  - va_ad(ji  ,jj  ,jk) / e2v(ji,jj)
               zhkead(ji  ,jj  ,jk) = zhkead(ji  ,jj  ,jk)              &
                  &                  + va_ad(ji  ,jj  ,jk) / e2v(ji,jj)
               zhkead(ji+1,jj  ,jk) = zhkead(ji+1,jj  ,jk) &
                  &                  - ua_ad(ji  ,jj  ,jk) / e1u(ji,jj)
               zhkead(ji  ,jj  ,jk) = zhkead(ji  ,jj  ,jk)              &
                  &                  + ua_ad(ji  ,jj  ,jk) / e1u(ji,jj)
            END DO
         END DO
         DO jj = jpj, 2, -1         ! Horizontal kinetic energy at T-point
            DO ji = jpi, 2, -1   ! vector opt.
               zuad = zhkead(ji,jj,jk)
               zvad = zhkead(ji,jj,jk)
               zhkead(ji,jj,jk) = 0.0_wp

               vn_ad(ji  ,jj-1,jk) = vn_ad(ji  ,jj-1,jk) + zvad * vn(ji  ,jj-1,jk) * 0.5_wp
               vn_ad(ji  ,jj  ,jk) = vn_ad(ji  ,jj  ,jk) + zvad * vn(ji  ,jj  ,jk) * 0.5_wp
               un_ad(ji-1,jj  ,jk) = un_ad(ji-1,jj  ,jk) + zuad * un(ji-1,jj  ,jk) * 0.5_wp
               un_ad(ji  ,jj  ,jk) = un_ad(ji  ,jj  ,jk) + zuad * un(ji  ,jj  ,jk) * 0.5_wp
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, jpk, zhkead )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_keg_adj')
      !
   END SUBROUTINE dyn_keg_adj
   SUBROUTINE dyn_keg_adj_tst( kumadt )
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit
      ! done in dynadv_tam
   END SUBROUTINE dyn_keg_adj_tst
END MODULE dynkeg_tam
