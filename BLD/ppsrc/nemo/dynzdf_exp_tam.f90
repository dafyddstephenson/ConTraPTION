MODULE dynzdf_exp_tam
   !!==============================================================================
   !!                     ***  MODULE  dynzdf_exp_tam  ***
   !! Ocean dynamics:  vertical component(s) of the momentum mixing trend
   !!                 Tangent and Adjoint Module
   !!==============================================================================
   !! History of the direct module:
   !!                !  90-10  (B. Blanke)  Original code
   !!                !  97-05  (G. Madec)  vertical component of isopycnal
   !!           8.5  !  02-08  (G. Madec)  F90: Free form and module
   !! History of the TAM module:
   !!           9.0  !  08-0!  (A. Vidard) Skeleton
   !!           3.4  !  12-07  (P.-A. Bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf_exp  : update the momentum trend with the vertical diffu-
   !!                  sion using an explicit time-stepping scheme.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce
   USE oce_tam
   USE zdf_oce
   USE dom_oce
   USE phycst
   USE in_out_manager
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_zdf_exp_tan     ! called by dynzdf_tam.F90
   PUBLIC dyn_zdf_exp_adj     ! called by dynzdf_tam.F90

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

CONTAINS

   SUBROUTINE dyn_zdf_exp_tan( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_exp_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the trend due to the vert. momentum diffusion
      !!
      !! ** Method of the direct routine:
      !!      Explicit forward time stepping with a time splitting
      !!      technique. The vertical diffusion of momentum is given by:
      !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
      !!      Surface boundary conditions: wind stress input
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F90)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) with the vertical diffusive trend
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER , INTENT( in ) ::   kt                           ! ocean time-step index
      REAL(wp), INTENT( in ) ::   p2dt                         ! time-step

      !! * Local declarations
      INTEGER ::   ji, jj, jk, jl                              ! dummy loop indices
      REAL(wp) ::   zrau0r, zlavmr, zuatl, zvatl                   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwxtl, zwytl, zwztl, zwwtl ! temporary workspace arrays
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_exp_tan')
      !
      CALL wrk_alloc( jpi,jpj,jpk, zwxtl, zwytl, zwztl, zwwtl )
      !
      IF( kt == nit000 .AND. lwp) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_exp_tan : vertical momentum diffusion explicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~ '
      ENDIF

      ! Local constant initialization
      ! -----------------------------
      zrau0r = 1. / rau0                                   ! inverse of the reference density
      zlavmr = 1. / REAL( nn_zdfexp )                      ! inverse of the number of sub time step
      !                                                ! ===============
                                                       !  Vertical slab
      !                                                ! ===============
      ! Surface boundary condition
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            zwytl(ji,jj,1) = 0.0_wp
            zwwtl(ji,jj,1) = 0.0_wp
         END DO
      END DO
      ! Initialization of x, z and contingently trends array
      DO jk = 1, jpk
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               zwxtl(ji,jj,jk) = ub_tl(ji,jj,jk)
               zwztl(ji,jj,jk) = vb_tl(ji,jj,jk)
            END DO
         END DO
      END DO
      ! Time splitting loop
      DO jl = 1, nn_zdfexp
         !
         ! First vertical derivative
         DO jk = 2, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zwytl(ji,jj,jk) = avmu(ji,jj,jk) * ( zwxtl(ji,jj,jk-1) - zwxtl(ji,jj,jk) ) / e3uw(ji,jj,jk)
                  zwwtl(ji,jj,jk) = avmv(ji,jj,jk) * ( zwztl(ji,jj,jk-1) - zwztl(ji,jj,jk) ) / e3vw(ji,jj,jk)
               END DO
            END DO
         END DO
         ! Second vertical derivative and trend estimation at kt+l*rdt/nn_zdfexp
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zuatl = zlavmr*( zwytl(ji,jj,jk) - zwytl(ji,jj,jk+1) ) / e3u(ji,jj,jk)
                  zvatl = zlavmr*( zwwtl(ji,jj,jk) - zwwtl(ji,jj,jk+1) ) / e3v(ji,jj,jk)
                  ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) + zuatl
                  va_tl(ji,jj,jk) = va_tl(ji,jj,jk) + zvatl
                  zwxtl(ji,jj,jk) = zwxtl(ji,jj,jk) + p2dt*zuatl*umask(ji,jj,jk)
                  zwztl(ji,jj,jk) = zwztl(ji,jj,jk) + p2dt*zvatl*vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      END DO
      !                                                ! ===============
      !                                                !   End of slab
      !                                                ! ===============
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zwxtl, zwytl, zwztl, zwwtl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf_exp_tan')
      !
   END SUBROUTINE dyn_zdf_exp_tan
   SUBROUTINE dyn_zdf_exp_adj( kt, p2dt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_exp_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the trend due to the vert. momentum diffusion
      !!
      !! ** Method of the direct routine:
      !!      Explicit forward time stepping with a time splitting
      !!      technique. The vertical diffusion of momentum is given by:
      !!         diffu = dz( avmu dz(u) ) = 1/e3u dk+1( avmu/e3uw dk(ub) )
      !!      Surface boundary conditions: wind stress input
      !!      Bottom boundary conditions : bottom stress (cf zdfbfr.F90)
      !!      Add this trend to the general trend ua :
      !!         ua = ua + dz( avmu dz(u) )
      !!
      !! ** Action : - Update (ua,va) with the vertical diffusive trend
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER , INTENT( in ) ::   kt                           ! ocean time-step index
      REAL(wp), INTENT( in ) ::   p2dt                         ! time-step

      !! * Local declarations
      INTEGER ::   ji, jj, jk, jl                                  ! dummy loop indices
      REAL(wp) ::   zrau0r, zlavmr, zuaad, zvaad                   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwxad, zwyad, zwzad, zwwad ! temporary workspace arrays
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_exp_adj')
      !
      CALL wrk_alloc( jpi,jpj,jpk, zwxad, zwyad, zwzad, zwwad )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_zdf_exp_adj : vertical momentum diffusion explicit operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~ '
      ENDIF
      ! Local constant initialization
      ! -----------------------------
      zrau0r = 1. / rau0                                   ! inverse of the reference density
      zlavmr = 1. / float( nn_zdfexp )                      ! inverse of the number of sub time step
      zwxad(:,:,:) = 0.0_wp ; zwyad(:,:,:) = 0.0_wp ; zwzad(:,:,:) = 0.0_wp ; zwwad(:,:,:) = 0.0_wp
      !                                                ! ===============
      !                                                !  Vertical slab
      !                                                ! ===============
      ! Time splitting loop
      DO jl = 1, nn_zdfexp
         !
         ! Second vertical derivative and trend estimation at kt+l*rdt/nn_zdfexp
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zuaad = p2dt * zwxad(ji,jj,jk) * umask(ji,jj,jk)
                  zvaad = p2dt * zwzad(ji,jj,jk) * vmask(ji,jj,jk)
                  zuaad = zuaad + ua_ad(ji,jj,jk)
                  zvaad = zvaad + va_ad(ji,jj,jk)
                  zwyad(ji,jj,jk  ) = zwyad(ji,jj,jk  ) + zlavmr * zuaad  / e3u(ji,jj,jk)
                  zwyad(ji,jj,jk+1) = zwyad(ji,jj,jk+1) - zlavmr * zuaad  / e3u(ji,jj,jk)
                  zwwad(ji,jj,jk  ) = zwwad(ji,jj,jk  ) + zlavmr * zvaad  / e3v(ji,jj,jk)
                  zwwad(ji,jj,jk+1) = zwwad(ji,jj,jk+1) - zlavmr * zvaad  / e3v(ji,jj,jk)
               END DO
            END DO
         END DO
         ! First vertical derivative
         DO jk = 2, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zwxad(ji,jj,jk-1) = zwxad(ji,jj,jk-1) + avmu(ji,jj,jk) * zwyad(ji,jj,jk) / e3uw(ji,jj,jk)
                  zwxad(ji,jj,jk  ) = zwxad(ji,jj,jk  ) - avmu(ji,jj,jk) * zwyad(ji,jj,jk) / e3uw(ji,jj,jk)
                  zwyad(ji,jj,jk  ) = 0.0_wp

                  zwzad(ji,jj,jk-1) = zwzad(ji,jj,jk-1) + avmv(ji,jj,jk) * zwwad(ji,jj,jk) / e3vw(ji,jj,jk)
                  zwzad(ji,jj,jk  ) = zwzad(ji,jj,jk  ) - avmv(ji,jj,jk) * zwwad(ji,jj,jk) / e3vw(ji,jj,jk)
                  zwwad(ji,jj,jk  ) = 0.0_wp
               END DO
            END DO
         END DO
         !
      END DO
      !
      ! Initialization of x, z and contingently trends array
      DO jk = 1, jpk
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ub_ad(ji,jj,jk) = ub_ad(ji,jj,jk) + zwxad(ji,jj,jk)
               vb_ad(ji,jj,jk) = vb_ad(ji,jj,jk) + zwzad(ji,jj,jk)
               zwxad(ji,jj,jk) = 0.0_wp
               zwzad(ji,jj,jk) = 0.0_wp
            END DO
         END DO
      END DO
      !                                                ! ===============
      !                                                !   End of slab
      !                                                ! ===============
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zwxad, zwyad, zwzad, zwwad )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf_exp_adj')
      !
   END SUBROUTINE dyn_zdf_exp_adj
   !!==============================================================================
END MODULE dynzdf_exp_tam
