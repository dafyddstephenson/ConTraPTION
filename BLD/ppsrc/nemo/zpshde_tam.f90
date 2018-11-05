MODULE zpshde_tam
   !!==============================================================================
   !!                       ***  MODULE zpshde_tam   ***
   !! z-coordinate - partial step : Horizontal Derivative
   !!                               Tangent and Adjoint Module
   !!==============================================================================

   !!----------------------------------------------------------------------
   !!   zps_hde      :  Horizontal DErivative of T, S and rd at the last
   !!                   ocean level (Z-coord. with Partial Steps)
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind
   USE par_oce
   USE oce
   USE oce_tam
   USE dom_oce
   USE in_out_manager
   USE eosbn2_tam
   USE lbclnk
   USE lbclnk_tam
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE lib_mpp
   USE lib_mpp_tam
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC zps_hde_tan      ! routine called by step_tam.F90
   PUBLIC zps_hde_adj      ! routine called by step_tam.F90
   PUBLIC zps_hde_adj_tst  ! routine called by tst.F90

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

   SUBROUTINE zps_hde_tan ( kt, kjpt, pta,            &
      &                      pta_tl,  prd_tl,    &
      &                      pgtu_tl, pgru_tl,   &
      &                      pgtv_tl, pgrv_tl     )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE zps_hde_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the horizontal derivative of T, S and rd
      !!      at u- and v-points with a linear interpolation for z-coordinate
      !!      with partial steps.
      !!
      !! ** Method of the direct routine:
      !!      In z-coord with partial steps, scale factors on last
      !!      levels are different for each grid point, so that T, S and rd
      !!      points are not at the same depth as in z-coord. To have horizontal
      !!      gradients again, we interpolate T and S at the good depth :
      !!      Linear interpolation of T, S
      !!         Computation of di(tb) and dj(tb) by vertical interpolation:
      !!          di(t) = t~ - t(i,j,k) or t(i+1,j,k) - t~
      !!          dj(t) = t~ - t(i,j,k) or t(i,j+1,k) - t~
      !!         This formulation computes the two cases:
      !!                 CASE 1                   CASE 2
      !!         k-1  ___ ___________   k-1   ___ ___________
      !!                    Ti  T~                  T~  Ti+1
      !!                  _____                        _____
      !!         k        |   |Ti+1     k           Ti |   |
      !!                  |   |____                ____|   |
      !!              ___ |   |   |           ___  |   |   |
      !!
      !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
      !!          t~ = t(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(Ti+1)/e3w(i+1)
      !!        ( t~ = t(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(Tj+1)/e3w(j+1)  )
      !!          or
      !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
      !!          t~ = t(i,j,k) + (e3w(i) - e3w(i+1)) * dk(Ti)/e3w(i )
      !!        ( t~ = t(i,j,k) + (e3w(j) - e3w(j+1)) * dk(Tj)/e3w(j ) )
      !!          Idem for di(s) and dj(s)
      !!
      !!      For rho, we call eos_insitu_2d which will compute rd~(t~,s~) at
      !!      the good depth zh from interpolated T and S for the different
      !!      formulation of the equation of state (eos).
      !!      Gradient formulation for rho :
      !!          di(rho) = rd~ - rd(i,j,k) or rd (i+1,j,k) - rd~
      !!
      !! ** Action  : - pgtu, pgsu, pgru: horizontal gradient of T, S
      !!                and rd at U-points
      !!              - pgtv, pgsv, pgrv: horizontal gradient of T, S
      !!                and rd at V-points
      !!
      !! History of the direct routine:
      !!   8.5  !  02-04  (A. Bozec)  Original code
      !!   8.5  !  02-08  (G. Madec E. Durand)  Optimization and Free form
      !! History of the TAM routine:
      !!   9.0  !  08-06 (A. Vidard) Skeleton
      !!        !  08-06 (A. Vidard) tangent of the 02-08 version
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt, kjpt          ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT( in ) ::   &
         pta, pta_tl            ! 3D T, S and rd direct fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( in ), OPTIONAL ::   &
         prd_tl
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( out ), OPTIONAL ::   &
         pgtu_tl, pgtv_tl                                    ! 3D T, S and rd tangent fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT( out ), OPTIONAL ::   &
         pgru_tl,                               &  ! horizontal grad. of T, S and rd at u-
         pgrv_tl                                   ! and v-points of the partial step level
      !! * Local declarations
      INTEGER ::   ji, jj,jk, jn,            &        ! Dummy loop indices
               &    iku,ikv, ikum1, ikvm1          ! partial step level at u- and v-points
      REAL(wp), POINTER, DIMENSION(:,:) ::   &
         zri, zrj,               &  ! and rd
         zritl, zrjtl,           &  ! and rdtl
         zhi, zhj                 ! depth of interpolation for eos2d
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zti, ztj, ztitl, ztjtl    ! interpolated value of tracer
      REAL(wp) ::   &
         ze3wu, ze3wv,           &  ! temporary scalars
         zmaxu, zmaxu2,         &  !    "         "
         zmaxv, zmaxv2             !    "         "
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'zps_hde_tan')
      !
      CALL wrk_alloc( jpi, jpj,       zri, zrj, zhi, zhj, zritl, zrjtl )
      CALL wrk_alloc( jpi, jpj, kjpt, zti, ztj, ztitl, ztjtl           )
      !
      DO jn = 1, kjpt
         ! Interpolation of T and S at the last ocean level
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)   ;   ikum1 = MAX( iku - 1 , 1 )    ! last and before last ocean level at u- & v-points
               ikv = mbkv(ji,jj)   ;   ikvm1 = MAX( ikv - 1 , 1 )    ! if level first is a p-step, ik.m1=1
               ze3wu = e3w(ji+1,jj  ,iku) - e3w(ji,jj,iku)
               ze3wv = e3w(ji  ,jj+1,ikv) - e3w(ji,jj,ikv)
               !
               ! i- direction
               IF( ze3wu >= 0. ) THEN      ! case 1
                  zmaxu =  ze3wu / e3w(ji+1,jj,iku)
                  ! interpolated values of T and S
                  zti(ji,jj,jn) = pta(ji+1,jj,iku,jn)            &
                     &       + zmaxu * ( pta(ji+1,jj,ikum1,jn) - pta(ji+1,jj,iku,jn) )
                  ztitl(ji,jj,jn) = pta_tl(ji+1,jj,iku,jn)       &
                     &         + zmaxu * ( pta_tl(ji+1,jj,ikum1,jn) - pta_tl(ji+1,jj,iku,jn) )
                  ! gradient of T and S
                  pgtu_tl(ji,jj,jn) = umask(ji,jj,1) * ( ztitl(ji,jj,jn) - pta_tl(ji,jj,iku,jn) )
               ELSE                        ! case 2
                  zmaxu = -ze3wu / e3w(ji,jj,iku)
                  ! interpolated values of T and S
                  zti(ji,jj,jn) = pta(ji,jj,iku,jn)              &
                     &       + zmaxu * ( pta(ji,jj,ikum1,jn) - pta(ji,jj,iku,jn) )
                  ! interpolated values of T and S
                  ztitl(ji,jj,jn) = pta_tl(ji,jj,iku,jn)         &
                     &         + zmaxu * ( pta_tl(ji,jj,ikum1,jn) - pta_tl(ji,jj,iku,jn) )
                  ! gradient of T and S
                  pgtu_tl(ji,jj,jn) = umask(ji,jj,1) * ( pta_tl(ji+1,jj,iku,jn) - ztitl (ji,jj,jn) )
               ENDIF
               ! j- direction
               IF( ze3wv >= 0. ) THEN      ! case 1
                  ! interpolated values of direct T and S
                  zmaxv =  ze3wv / e3w(ji,jj+1,ikv)
                  ztj(ji,jj,jn) = pta(ji,jj+1,ikv,jn)            &
                     &       + zmaxv * ( pta(ji,jj+1,ikvm1,jn) - pta(ji,jj+1,ikv,jn) )
                  ! interpolated values of tangent T and S
                  ztjtl(ji,jj,jn) = pta_tl(ji,jj+1,ikv,jn)       &
                     &         + zmaxv * ( pta_tl(ji,jj+1,ikvm1,jn) - pta_tl(ji,jj+1,ikv,jn) )
                  ! gradient of T and S
                  pgtv_tl(ji,jj,jn) = vmask(ji,jj,1) * ( ztjtl(ji,jj,jn) - pta_tl(ji,jj,ikv,jn) )
               ELSE                        ! case 2
                  zmaxv =  -ze3wv / e3w(ji,jj,ikv)
                  ! interpolated values of T and S
                  ztj(ji,jj,jn) = pta(ji,jj,ikv,jn)       &
                     &       + zmaxv * ( pta(ji,jj,ikvm1,jn) - pta(ji,jj,ikv,jn) )
                  ! interpolated values of T and S
                  ztjtl(ji,jj,jn) = pta_tl(ji,jj,ikv,jn)  &
                     &         + zmaxv * ( pta_tl(ji,jj,ikvm1,jn) - pta_tl(ji,jj,ikv,jn) )
                  ! gradient of T and S
                  pgtv_tl(ji,jj,jn) = vmask(ji,jj,1) * ( pta_tl(ji,jj+1,ikv,jn) - ztjtl(ji,jj,jn) )
               ENDIF
            END DO
         END DO
         CALL lbc_lnk( pgtu_tl(:,:,jn), 'U', -1. )   ;   CALL lbc_lnk( pgtv_tl(:,:,jn), 'V', -1. )   ! Lateral boundary cond.
         !
      END DO
      !
      ! horizontal derivative of density anomalies (rd)
      IF( PRESENT( prd_tl ) ) THEN         ! depth of the partial step level
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)
               ikv = mbkv(ji,jj)
               ze3wu  = e3w(ji+1,jj  ,iku) - e3w(ji,jj,iku)
               ze3wv  = e3w(ji  ,jj+1,ikv) - e3w(ji,jj,ikv)
               IF( ze3wu >= 0._wp ) THEN   ;   zhi(ji,jj) = gdept(ji  ,jj,iku)     ! i-direction: case 1
               ELSE                        ;   zhi(ji,jj) = gdept(ji+1,jj,iku)     ! -     -      case 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN   ;   zhj(ji,jj) = gdept(ji,jj  ,ikv)     ! j-direction: case 1
               ELSE                        ;   zhj(ji,jj) = gdept(ji,jj+1,ikv)     ! -     -      case 2
               ENDIF
            END DO
         END DO
         ! Compute interpolated rd from zti, zsi, ztj, zsj for the 2 cases at the depth of the partial
         ! step and store it in  zri, zrj for each  case
         CALL eos_tan( zti, zhi, ztitl, zritl )
         CALL eos_tan( ztj, zhj, ztjtl, zrjtl )


         ! Gradient of density at the last level
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)
               ikv = mbkv(ji,jj)
               ze3wu  = e3w(ji+1,jj  ,iku) - e3w(ji,jj,iku)
               ze3wv  = e3w(ji  ,jj+1,ikv) - e3w(ji,jj,ikv)
               IF( ze3wu >= 0. ) THEN    ! i-direction: case 1
                  pgru_tl(ji,jj) = umask(ji,jj,1) * ( zritl(ji,jj) - prd_tl(ji,jj,iku) )
               ELSE                      ! i-direction: case 2
                  pgru_tl(ji,jj) = umask(ji,jj,1) * ( prd_tl(ji+1,jj,iku) - zritl(ji,jj) )
               ENDIF
               IF( ze3wv >= 0. ) THEN    ! j-direction: case 1
                  pgrv_tl(ji,jj) = vmask(ji,jj,1) * ( zrjtl(ji,jj) - prd_tl(ji,jj,ikv) )
               ELSE                      ! j-direction: case 2
                  pgrv_tl(ji,jj) = vmask(ji,jj,1) * ( prd_tl(ji,jj+1,ikv) - zrjtl(ji,jj) )
               ENDIF
            END DO
         END DO
         ! Lateral boundary conditions on each gradient
         CALL lbc_lnk( pgru_tl , 'U', -1.0_wp )  ;  CALL lbc_lnk( pgrv_tl , 'V', -1.0_wp )
      END IF
      !
      CALL wrk_dealloc( jpi, jpj,       zri, zrj, zhi, zhj, zritl, zrjtl )
      CALL wrk_dealloc( jpi, jpj, kjpt, zti, ztj, ztitl, ztjtl           )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'zps_hde_tan')
      !
   END SUBROUTINE zps_hde_tan

   SUBROUTINE zps_hde_adj ( kt, kjpt, pta,       &
      &                      pta_ad,  prd_ad,    &
      &                      pgtu_ad, pgru_ad,   &
      &                      pgtv_ad, pgrv_ad     )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE zps_hde_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the horizontal derivative of T, S and rd
      !!      at u- and v-points with a linear interpolation for z-coordinate
      !!      with partial steps.
      !!
      !! ** Method of the direct routine:
      !!      In z-coord with partial steps, scale factors on last
      !!      levels are different for each grid point, so that T, S and rd
      !!      points are not at the same depth as in z-coord. To have horizontal
      !!      gradients again, we interpolate T and S at the good depth :
      !!      Linear interpolation of T, S
      !!         Computation of di(tb) and dj(tb) by vertical interpolation:
      !!          di(t) = t~ - t(i,j,k) or t(i+1,j,k) - t~
      !!          dj(t) = t~ - t(i,j,k) or t(i,j+1,k) - t~
      !!         This formulation computes the two cases:
      !!                 CASE 1                   CASE 2
      !!         k-1  ___ ___________   k-1   ___ ___________
      !!                    Ti  T~                  T~  Ti+1
      !!                  _____                        _____
      !!         k        |   |Ti+1     k           Ti |   |
      !!                  |   |____                ____|   |
      !!              ___ |   |   |           ___  |   |   |
      !!
      !!      case 1->   e3w(i+1) >= e3w(i) ( and e3w(j+1) >= e3w(j) ) then
      !!          t~ = t(i+1,j  ,k) + (e3w(i+1) - e3w(i)) * dk(Ti+1)/e3w(i+1)
      !!        ( t~ = t(i  ,j+1,k) + (e3w(j+1) - e3w(j)) * dk(Tj+1)/e3w(j+1)  )
      !!          or
      !!      case 2->   e3w(i+1) <= e3w(i) ( and e3w(j+1) <= e3w(j) ) then
      !!          t~ = t(i,j,k) + (e3w(i) - e3w(i+1)) * dk(Ti)/e3w(i )
      !!        ( t~ = t(i,j,k) + (e3w(j) - e3w(j+1)) * dk(Tj)/e3w(j ) )
      !!          Idem for di(s) and dj(s)
      !!
      !!      For rho, we call eos_insitu_2d which will compute rd~(t~,s~) at
      !!      the good depth zh from interpolated T and S for the different
      !!      formulation of the equation of state (eos).
      !!      Gradient formulation for rho :
      !!          di(rho) = rd~ - rd(i,j,k) or rd (i+1,j,k) - rd~
      !!
      !! ** Action  : - pgtu, pgsu, pgru: horizontal gradient of T, S
      !!                and rd at U-points
      !!              - pgtv, pgsv, pgrv: horizontal gradient of T, S
      !!                and rd at V-points
      !!
      !! History of the direct routine:
      !!   8.5  !  02-04  (A. Bozec)  Original code
      !!   8.5  !  02-08  (G. Madec E. Durand)  Optimization and Free form
      !! History of the TAM routine:
      !!   9.0  !  08-06 (A. Vidard) Skeleton
      !!        !  08-08 (A. Vidard) adjoint of the 02-08 version
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt, kjpt          ! ocean time-step index
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT( inout ) ::   &
         pta, pta_ad            ! 3D T, S and rd direct fields
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT( inout ), OPTIONAL ::   &
         prd_ad                                              ! 3D T, S and rd tangent fields
      REAL(wp), DIMENSION(jpi,jpj,kjpt), INTENT( inout ), OPTIONAL ::   &
         pgtu_ad, pgtv_ad                                     ! 3D T, S and rd tangent fields
      REAL(wp), DIMENSION(jpi,jpj), INTENT( inout ), OPTIONAL ::   &
         pgru_ad,                               &  ! horizontal grad. of T, S and rd at u-
         pgrv_ad                                   ! and v-points of the partial step level
      !! * Local declarations
      INTEGER ::   ji, jj,jk, jn,  &                    ! Dummy loop indices
                &   iku,ikv, ikum1, ikvm1          ! partial step level at u- and v-points
      REAL(wp), POINTER, DIMENSION(:,:) ::   &
         zri, zrj,               &  ! and rd
         zriad, zrjad,           &  ! and rdtl
         zhi, zhj                 ! depth of interpolation for eos2d
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  zti, ztj, ztiad, ztjad    ! interpolated value of tracer
      REAL(wp) ::   &
         ze3wu, ze3wv,           &  ! temporary scalars
         zmaxu, zmaxu2,         &  !    "         "
         zmaxv, zmaxv2             !    "         "
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'zps_hde_adj')
      !
      CALL wrk_alloc( jpi, jpj,       zri, zrj, zhi, zhj, zriad, zrjad )
      CALL wrk_alloc( jpi, jpj, kjpt, zti, ztj, ztiad, ztjad           )
      !
      ! 1. Direct model recomputation
      DO jn = 1, kjpt      !==   Interpolation of tracers at the last ocean level   ==!
         !
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)   ;   ikum1 = MAX( iku - 1 , 1 )    ! last and before last ocean level at u- & v-points
               ikv = mbkv(ji,jj)   ;   ikvm1 = MAX( ikv - 1 , 1 )    ! if level first is a p-step, ik.m1=1
               ze3wu = e3w(ji+1,jj  ,iku) - e3w(ji,jj,iku)
               ze3wv = e3w(ji  ,jj+1,ikv) - e3w(ji,jj,ikv)
               !
               ! i- direction
               IF( ze3wu >= 0._wp ) THEN      ! case 1
                  zmaxu =  ze3wu / e3w(ji+1,jj,iku)
                  ! interpolated values of tracers
                  zti(ji,jj,jn) = pta(ji+1,jj,iku,jn) + zmaxu * ( pta(ji+1,jj,ikum1,jn) - pta(ji+1,jj,iku,jn) )
               ELSE                           ! case 2
                  zmaxu = -ze3wu / e3w(ji,jj,iku)
                  ! interpolated values of tracers
                  zti(ji,jj,jn) = pta(ji,jj,iku,jn) + zmaxu * ( pta(ji,jj,ikum1,jn) - pta(ji,jj,iku,jn) )
               ENDIF
               !
               ! j- direction
               IF( ze3wv >= 0._wp ) THEN      ! case 1
                  zmaxv =  ze3wv / e3w(ji,jj+1,ikv)
                  ! interpolated values of tracers
                  ztj(ji,jj,jn) = pta(ji,jj+1,ikv,jn) + zmaxv * ( pta(ji,jj+1,ikvm1,jn) - pta(ji,jj+1,ikv,jn) )
               ELSE                           ! case 2
                  zmaxv =  -ze3wv / e3w(ji,jj,ikv)
                  ! interpolated values of tracers
                  ztj(ji,jj,jn) = pta(ji,jj,ikv,jn) + zmaxv * ( pta(ji,jj,ikvm1,jn) - pta(ji,jj,ikv,jn) )
               ENDIF
            END DO
         END DO
         !
      END DO

      !2. Adjoint model counterpart
      ztiad = 0.0_wp ; ztjad = 0.0_wp
      zriad = 0.0_wp ; zrjad = 0.0_wp
! horizontal derivative of density anomalies (rd)
      IF( PRESENT( prd_ad ) ) THEN         ! depth of the partial step level
         ! Lateral boundary conditions on each gradient
         CALL lbc_lnk_adj( pgru_ad , 'U', -1.0_wp )
         CALL lbc_lnk_adj( pgrv_ad , 'V', -1.0_wp )
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               iku = mbku(ji,jj)
               ikv = mbkv(ji,jj)
               ze3wu  = e3w(ji+1,jj  ,iku) - e3w(ji,jj,iku)
               ze3wv  = e3w(ji  ,jj+1,ikv) - e3w(ji,jj,ikv)
               IF( ze3wu >= 0._wp ) THEN   ;   zhi(ji,jj) = gdept(ji  ,jj,iku)     ! i-direction: case 1
               ELSE                        ;   zhi(ji,jj) = gdept(ji+1,jj,iku)     ! -     -      case 2
               ENDIF
               IF( ze3wv >= 0._wp ) THEN   ;   zhj(ji,jj) = gdept(ji,jj  ,ikv)     ! j-direction: case 1
               ELSE                        ;   zhj(ji,jj) = gdept(ji,jj+1,ikv)     ! -     -      case 2
               ENDIF
            END DO
         END DO
         ! Gradient of density at the last level
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1
               iku = mbku(ji,jj)
               ikv = mbkv(ji,jj)
               ze3wu  = e3w(ji+1,jj  ,iku) - e3w(ji,jj,iku)
               ze3wv  = e3w(ji  ,jj+1,ikv) - e3w(ji,jj,ikv)
               IF( ze3wv >= 0. ) THEN    ! j-direction: case 1
                  zrjad(ji,jj)        = zrjad(ji,jj)                     &
                     &                + pgrv_ad(ji,jj) * vmask(ji,jj,1)
                  prd_ad(ji,jj,ikv)   = prd_ad(ji,jj,ikv)                &
                     &                - pgrv_ad(ji,jj) * vmask(ji,jj,1)
                  pgrv_ad(ji,jj)      = 0.0_wp
               ELSE                      ! j-direction: case 2
                  prd_ad(ji,jj+1,ikv) = prd_ad(ji,jj+1,ikv)            &
                     &                + pgrv_ad(ji,jj) * vmask(ji,jj,1)
                  zrjad(ji,jj)        = zrjad(ji,jj)                     &
                     &                - pgrv_ad(ji,jj) * vmask(ji,jj,1)
                  pgrv_ad(ji,jj)      = 0.0_wp
               ENDIF
               IF( ze3wu >= 0. ) THEN    ! i-direction: case 1
                  zriad(ji,jj)        = zriad(ji,jj)                     &
                     &                + pgru_ad(ji,jj) * umask(ji,jj,1)
                  prd_ad(ji,jj,iku)   = prd_ad(ji,jj,iku)                &
                     &                - pgru_ad(ji,jj) * umask(ji,jj,1)
                  pgru_ad(ji,jj)      = 0.0_wp
               ELSE                      ! i-direction: case 2
                  prd_ad(ji+1,jj,iku) = prd_ad(ji+1,jj,iku)            &
                     &                + pgru_ad(ji,jj) * umask(ji,jj,1)
                  zriad(ji,jj)        = zriad(ji,jj)                   &
                     &                - pgru_ad(ji,jj) * umask(ji,jj,1)
                  pgru_ad(ji,jj)      = 0.0_wp
               ENDIF

            END DO
         END DO

         ! Compute interpolated rd from zti, zsi, ztj, zsj for the 2 cases at the depth of the partial
         ! step and store it in  zri, zrj for each  case
         CALL eos_adj( ztj, zhj, ztjad, zrjad )
         CALL eos_adj( zti, zhi, ztiad, zriad )
      END IF

      DO jn = 1, kjpt
         CALL lbc_lnk_adj( pgtu_ad(:,:,jn), 'U', -1. )   ;   CALL lbc_lnk_adj( pgtv_ad(:,:,jn), 'V', -1. )   ! Lateral boundary cond.
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1
               iku = mbku(ji,jj)   ;   ikum1 = MAX( iku - 1 , 1 )    ! last and before last ocean level at u- & v-points
               ikv = mbkv(ji,jj)   ;   ikvm1 = MAX( ikv - 1 , 1 )    ! if level first is a p-step, ik.m1=1
               ze3wu = e3w(ji+1,jj  ,iku) - e3w(ji,jj,iku)
               ze3wv = e3w(ji  ,jj+1,ikv) - e3w(ji,jj,ikv)
               !
               ! j- direction
               IF( ze3wv >= 0. ) THEN      ! case 1
                  zmaxv =  ze3wv / e3w(ji,jj+1,ikv)
                  ! gradient of T and S
                  ztjad(ji,jj,jn)           = ztjad(ji,jj,jn)                &
                     &                   + pgtv_ad(ji,jj,jn) * vmask(ji,jj,1)
                  pta_ad(ji,jj,ikv,jn)     = pta_ad(ji,jj,ikv,jn)          &
                     &                   - pgtv_ad(ji,jj,jn) * vmask(ji,jj,1)
                  pgtv_ad(ji,jj,jn)         = 0.0_wp
                  ! interpolated values of T and S
                  pta_ad(ji,jj+1,ikv,jn)   = pta_ad(ji,jj+1,ikv,jn)        &
                     &                   + ztjad(ji,jj,jn) * (1 - zmaxv)
                  pta_ad(ji,jj+1,ikvm1,jn) = pta_ad(ji,jj+1,ikvm1,jn)      &
                     &                   + ztjad(ji,jj,jn)* zmaxv
                  ztjad(ji,jj,jn)           = 0.0_wp
               ELSE                        ! case 2
                  zmaxv =  -ze3wv / e3w(ji,jj,ikv)
                  ! gradient of T and S
                  pta_ad(ji,jj+1,ikv,jn)   = pta_ad(ji,jj+1,ikv,jn)        &
                     &                   + pgtv_ad(ji,jj,jn) * vmask(ji,jj,1)
                  ztjad(ji,jj,jn)           = ztjad(ji,jj,jn)                &
                     &                   - pgtv_ad(ji,jj,jn) * vmask(ji,jj,1)
                  pgtv_ad(ji,jj,jn)         = 0.0_wp

                  ! interpolated values of T and S
                  pta_ad(ji,jj,ikv,jn)     = pta_ad(ji,jj,ikv,jn)          &
                     &                   + ztjad(ji,jj,jn) * (1 - zmaxv)
                  pta_ad(ji,jj,ikvm1,jn)   = pta_ad(ji,jj,ikvm1,jn)        &
                     &                   + ztjad(ji,jj,jn) * zmaxv
                  ztjad(ji,jj,jn)           = 0.0_wp
               ENDIF
               ! i- direction
               IF( ze3wu >= 0. ) THEN      ! case 1
                  zmaxu =  ze3wu / e3w(ji+1,jj,iku)
                  ! gradient of T and S
                  ztiad(ji,jj,jn)       = ztiad(ji,jj,jn)                    &
                     &               + pgtu_ad(ji,jj,jn) * umask(ji,jj,1)
                  pta_ad(ji,jj,iku,jn) = pta_ad(ji,jj,iku,jn)              &
                     &               - pgtu_ad(ji,jj,jn) * umask(ji,jj,1)
                  pgtu_ad(ji,jj,jn)     = 0.0_wp
                  ! interpolated values of T and S
                  pta_ad(ji+1,jj,iku,jn)   = pta_ad(ji+1,jj,iku,jn)        &
                     &                   + ztiad(ji,jj,jn) * (1 - zmaxu)
                  pta_ad(ji+1,jj,ikum1,jn) = pta_ad(ji+1,jj,ikum1,jn)      &
                     &                   + ztiad(ji,jj,jn) * zmaxu
                  ztiad(ji,jj,jn)           = 0.0_wp
               ELSE                        ! case 2
                  zmaxu = -ze3wu / e3w(ji,jj,iku)
                  ! gradient of T and S
                  pta_ad(ji+1,jj,iku,jn)   = pta_ad(ji+1,jj,iku,jn)        &
                     &                   + pgtu_ad(ji,jj,jn) * umask(ji,jj,1)
                  ztiad (ji,jj,jn)          = ztiad (ji,jj,jn)               &
                     &                   - pgtu_ad(ji,jj,jn) * umask(ji,jj,1)
                  pgtu_ad(ji,jj,jn)         = 0.0_wp
                  ! interpolated values of T and S
                  pta_ad(ji,jj,iku,jn)     = pta_ad(ji,jj,iku,jn)          &
                     &                   + ztiad(ji,jj,jn) * (1 - zmaxu)
                  pta_ad(ji,jj,ikum1,jn)   = pta_ad(ji,jj,ikum1,jn)        &
                     &                   + ztiad(ji,jj,jn) * zmaxu
                  ztiad(ji,jj,jn)           = 0.0_wp
               ENDIF
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj,       zri, zrj, zhi, zhj, zriad, zrjad )
      CALL wrk_dealloc( jpi, jpj, kjpt, zti, ztj, ztiad, ztjad           )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'zps_hde_adj')
      !
   END SUBROUTINE zps_hde_adj
      SUBROUTINE zps_hde_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE zps_hde_adj_tst ***
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
      !!        ! 08-08 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER :: &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & kt,    &
         & jn

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:,:), POINTER :: &
         & zts,         & ! Direct field : temperature
         & zts_tlin,    & ! Tangent input: temperature
         & zts_adout,   & ! Adjoint output: temperature
         & zats           ! 3D random field for t

      REAL(KIND=wp), DIMENSION(:,:,:), POINTER :: &
         & zgtu_tlout,   & ! Tangent output: horizontal gradient
         & zgtv_tlout,   & ! Tangent output: horizontal gradient
         & zrd_adout,    & ! Adjoint output:
         & zar,          & ! 3D random field for rd
         & zrd_tlin,     & ! Tangent input:
         & zgtu_adin,    & ! Adjoint input : horizontal gradient
         & zgtv_adin       ! Adjoint input : horizontal gradient

      REAL(KIND=wp), DIMENSION(:,:), POINTER :: &
         & zgru_tlout,   & ! Tangent output: horizontal gradient
         & zgrv_tlout,   & ! Tangent output: horizontal gradient
         & zgru_adin,    & ! Adjoint input : horizontal gradient
         & zgrv_adin       ! Adjoint input : horizontal gradient

      REAL(KIND=wp) :: &
                           ! random field standard deviation for:
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp1_3,       & !   scalar product components
         & zsp1_4,       &
         & zsp1_5,       & !   scalar product components
         & zsp1_6,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3
      CHARACTER (LEN=14) :: &
         & cl_name

      kt = nit000
      ! Allocate memory
      CALL wrk_alloc(jpi, jpj, jpk, jpts, zts, zts_tlin,zts_adout, zats )

      CALL wrk_alloc(jpi, jpj, jpts, zgtu_tlout, zgtv_tlout, zgtu_adin, zgtv_adin )

      CALL wrk_alloc(jpi, jpj, jpk , zrd_adout, zrd_tlin, zar )

      CALL wrk_alloc(jpi, jpj, zgru_tlout, zgrv_tlout, zgru_adin, zgrv_adin )

      ! Initialize random field standard deviationsthe reference state
      zts = tsn(:,:,:,:)

      !=============================================================
      ! 1) dx = ( T ) and dy = ( T )
      !=============================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      zts_tlin(:,:,:,:)  = 0.0_wp
      zrd_tlin(:,:,:)   = 0.0_wp
      zts_adout(:,:,:,:) = 0.0_wp
      zrd_adout(:,:,:)  = 0.0_wp
      zgtu_tlout(:,:,:)   = 0.0_wp
      zgtv_tlout(:,:,:)   = 0.0_wp
      zgru_tlout(:,:)   = 0.0_wp
      zgrv_tlout(:,:)   = 0.0_wp
      zgtu_adin(:,:,:)    = 0.0_wp
      zgtv_adin(:,:,:)    = 0.0_wp
      zgru_adin(:,:)    = 0.0_wp
      zgrv_adin(:,:)    = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------
      DO jn = 1, jpts
         CALL grid_random(  zats(:,:,:,jn), 'T', 0.0_wp, stdt )
      END DO
      CALL grid_random(  zar, 'T', 0.0_wp, stdr )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zts_tlin(ji,jj,:,:) = zats(ji,jj,:,:)
            zrd_tlin(ji,jj,:)   = zar (ji,jj,:)
         END DO
      END DO

      CALL zps_hde_tan ( nit000, jpts, zts,           &
         &                   zts_tlin , zrd_tlin  ,   &
         &                   zgtu_tlout, zgru_tlout,  &
         &                   zgtv_tlout, zgrv_tlout   )

      DO jn = 1, jpts
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               jk = mbku(ji,jj)
               zgtu_adin(ji,jj,jn) = zgtu_tlout(ji,jj,jn) &
                  &             * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
                  &             * umask(ji,jj,jk)
               jk = mbkv(ji,jj)
               zgtv_adin(ji,jj,jn) = zgtv_tlout(ji,jj,jn) &
                  &             * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) &
                  &             * vmask(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            jk = mbku(ji,jj)
            zgru_adin(ji,jj) = zgru_tlout(ji,jj) &
               &             * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
               &             * umask(ji,jj,jk)
            jk = mbkv(ji,jj)
            zgrv_adin(ji,jj) = zgrv_tlout(ji,jj) &
               &             * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) &
               &             * vmask(ji,jj,jk)
         END DO
      END DO
      !!--------------------------------------------------------------------
      !! Compute the scalar product: ( L dx )^T W dy
      !!--------------------------------------------------------------------

      zsp1_1 = DOT_PRODUCT( zgtu_adin(:,:,jp_tem), zgtu_tlout(:,:,jp_tem) )
      zsp1_2 = DOT_PRODUCT( zgtu_adin(:,:,jp_sal), zgtu_tlout(:,:,jp_sal) )
      zsp1_3 = DOT_PRODUCT( zgru_adin, zgru_tlout )
      zsp1_4 = DOT_PRODUCT( zgtv_adin(:,:,jp_tem), zgtv_tlout(:,:,jp_tem) )
      zsp1_5 = DOT_PRODUCT( zgtv_adin(:,:,jp_sal), zgtv_tlout(:,:,jp_sal) )
      zsp1_6 = DOT_PRODUCT( zgrv_adin, zgrv_tlout )
      zsp1 = zsp1_1 + zsp1_2 + zsp1_3 + zsp1_4 + zsp1_5 + zsp1_6


      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------
      CALL zps_hde_adj ( kt, jpts, zts ,   &
         &                   zts_adout , zrd_adout ,   &
         &                   zgtu_adin , zgru_adin ,   &
         &                   zgtv_adin , zgrv_adin     )

      zsp2_1 = DOT_PRODUCT( zts_tlin(:,:,:,jp_tem), zts_adout(:,:,:,jp_tem) )
      zsp2_2 = DOT_PRODUCT( zts_tlin(:,:,:,jp_sal), zts_adout(:,:,:,jp_sal) )
      zsp2_3 = DOT_PRODUCT( zrd_tlin , zrd_adout  )
      zsp2   = zsp2_1 + zsp2_2 + zsp2_3

      ! Compare the scalar products

      cl_name = 'zps_hde_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      ! Deallocate memory
      CALL wrk_dealloc(jpi, jpj, jpk, jpts, zts, zts_tlin,zts_adout, zats )

      CALL wrk_dealloc(jpi, jpj, jpts, zgtu_tlout, zgtv_tlout, zgtu_adin, zgtv_adin )

      CALL wrk_dealloc(jpi, jpj, jpk , zrd_adout, zrd_tlin, zar )

      CALL wrk_dealloc(jpi, jpj, zgru_tlout, zgrv_tlout, zgru_adin, zgrv_adin )
     

   END SUBROUTINE zps_hde_adj_tst
END MODULE zpshde_tam
