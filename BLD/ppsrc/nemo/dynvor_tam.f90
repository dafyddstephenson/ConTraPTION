MODULE dynvor_tam
   !!======================================================================
   !!                       ***  MODULE  dynvor_tam  ***
   !! Ocean dynamics: Update the momentum trend with the relative and
   !!                 planetary vorticity trends
   !!                 Tangent and Adjoint Module
   !!======================================================================
   !! History of the drect module:
   !! History :  OPA  !  1989-12  (P. Andrich)  vor_ens: Original code
   !!            5.0  !  1991-11  (G. Madec) vor_ene, vor_mix: Original code
   !!            6.0  !  1996-01  (G. Madec)  s-coord, suppress work arrays
   !!            8.5  !  2002-08  (G. Madec)  F90: Free form and module
   !!   NEMO     1.0  !  2004-02  (G. Madec)  vor_een: Original code
   !!             -   !  2003-08  (G. Madec)  add vor_ctl
   !!             -   !  2005-11  (G. Madec)  add dyn_vor (new step architecture)
   !!            2.0  !  2006-11  (G. Madec)  flux form advection: add metric term
   !!            3.2  !  2009-04  (R. Benshila)  vvl: correction of een scheme
   !! History of the TAM module:
   !!            9.0  ! 2008-06  (A. Vidard) Skeleton
   !!            9.0  ! 2009-01  (A. Vidard) TAM of the 06-11 version
   !!            9.0  ! 2010-01  (F. Vigilant) Add een TAM option
   !!  NEMO      3.2  ! 2010-04 (F. Vigilant) 3.2 version
   !!  NEMO      3.4  ! 2012-07 (P.-A. Bouttier) 3.4 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_vor     : Update the momentum trend with the vorticity trend
   !!       vor_ens : enstrophy conserving scheme       (ln_dynvor_ens=T)
   !!       vor_ene : energy conserving scheme          (ln_dynvor_ene=T)
   !!       vor_mix : mixed enstrophy/energy conserving (ln_dynvor_mix=T)
   !!       vor_een : energy and enstrophy conserving   (ln_dynvor_een=T)
   !!       vor_ctl : set and control of the different vorticity option
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE oce_tam
   USE divcur
   USE divcur_tam
   USE dom_oce
   USE dynadv
   USE dynvor
   USE lbclnk
   USE lbclnk_tam
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE dommsk
   USE lib_mpp
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_vor_init_tam
   PUBLIC   dyn_vor_tan     ! routine called by step_tam.F90
   PUBLIC   dyn_vor_adj     ! routine called by step_tam.F90
   PUBLIC   dyn_vor_adj_tst ! routine called by the tst.F90

   INTEGER ::   nvor = 0   ! type of vorticity trend used
   INTEGER ::   ncor = 1   ! coriolis
   INTEGER ::   nrvm = 2   ! =2 relative vorticity ; =3 metric term
   INTEGER ::   ntot = 4   ! =4 total vorticity (relative + planetary) ; =5 coriolis + metric term

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

   SUBROUTINE dyn_vor_tan( kt )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose of the direct routine:
      !!               compute the lateral ocean tracer physics.
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_vor_tan')
      !
      !!                                          ! vorticity term
      SELECT CASE ( nvor )                       ! compute the vorticity trend and add it to the general trend
      !
      CASE ( -1 )                                      ! esopa: test all possibility with control print
         CALL vor_ene_tan( kt, ntot, ua_tl, va_tl )
         CALL vor_ens_tan( kt, ntot, ua_tl, va_tl )
         CALL vor_mix_tan( kt )
         CALL vor_een_tan( kt, ntot, ua_tl, va_tl )
         !
      CASE ( 0 )                                       ! energy conserving scheme
         CALL vor_ene_tan( kt, ntot, ua_tl, va_tl )                ! total vorticity
         !
      CASE ( 1 )                                       ! enstrophy conserving scheme
         CALL vor_ens_tan( kt, ntot, ua_tl, va_tl )                ! total vorticity
         !
      CASE ( 2 )                                       ! mixed ene-ens scheme
         CALL vor_mix_tan( kt )                               ! total vorticity (mix=ens-ene)
         !
      CASE ( 3 )                                       ! energy and enstrophy conserving scheme
         CALL vor_een_tan( kt, ntot, ua_tl, va_tl )                ! total vorticity
         !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_vor_tan')
      !
   END SUBROUTINE dyn_vor_tan
   SUBROUTINE vor_ene_tan( kt, kvor, pua_tl, pva_tl )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE vor_ene  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Sadourny (1975) flux form formulation : conserves the
      !!      horizontal kinetic energy.
      !!      The trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivatives:
      !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f  mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f  mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1[ (rotn+f)  mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)  mj(e2u un) ]
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
      !                                                           ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua_tl ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva_tl ! total v-trend
      !!
      INTEGER  ::   ji, jj, jk         ! dummy loop indices
      REAL(wp) ::   zx1, zy1, zfact2   ! temporary scalars
      REAL(wp) ::   zx1tl, zy1tl   ! temporary scalars
      REAL(wp) ::   zx2, zy2           !    "         "
      REAL(wp) ::   zx2tl, zy2tl           !    "         "
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwx, zwy, zwz   ! temporary 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwxtl, zwytl, zwztl   ! temporary 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_ene_tan')
      !
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zwz )
      CALL wrk_alloc( jpi, jpj, zwxtl, zwytl, zwztl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_ene_tan : vorticity term: energy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      zfact2 = 0.5 * 0.5      ! Local constant initialization

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )   ;   zwz(:,:) =                  ff(:,:)  ; zwztl(:,:) = 0.
         CASE ( 2 )   ;    zwz(:,:) =   rotn(:,:,jk)  ; zwztl(:,:) =   rotn_tl(:,:,jk)                ! relative  vorticity
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
                  zwztl(ji,jj) = (   ( vn_tl(ji+1,jj  ,jk) + vn_tl (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un_tl(ji  ,jj+1,jk) + un_tl (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
               END DO
            END DO
         CASE ( 4 )    ;   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) ;   zwztl(:,:) = rotn_tl(:,:,jk)    ! total (relative + planetary vorticity)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                                                                       &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )                                               &
                       &       )
                  zwztl(ji,jj) = (   ( vn_tl(ji+1,jj  ,jk) + vn_tl (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un_tl(ji  ,jj+1,jk) + un_tl (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
               END DO
            END DO
         END SELECT

         IF( ln_sco ) THEN
            zwz(:,:) = zwz(:,:) / e3f(:,:,jk)
            zwx(:,:) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)
            zwztl(:,:) = zwztl(:,:) / e3f(:,:,jk)
            zwxtl(:,:) = e2u(:,:) * e3u(:,:,jk) * un_tl(:,:,jk)
            zwytl(:,:) = e1v(:,:) * e3v(:,:,jk) * vn_tl(:,:,jk)
         ELSE
            zwx(:,:) = e2u(:,:) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * vn(:,:,jk)
            zwxtl(:,:) = e2u(:,:) * un_tl(:,:,jk)
            zwytl(:,:) = e1v(:,:) * vn_tl(:,:,jk)
         ENDIF

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zy1 = zwy(ji,jj-1) + zwy(ji+1,jj-1)
               zy2 = zwy(ji,jj  ) + zwy(ji+1,jj  )
               zx1 = zwx(ji-1,jj) + zwx(ji-1,jj+1)
               zx2 = zwx(ji  ,jj) + zwx(ji  ,jj+1)
               zy1tl = zwytl(ji,jj-1) + zwytl(ji+1,jj-1)
               zy2tl = zwytl(ji,jj  ) + zwytl(ji+1,jj  )
               zx1tl = zwxtl(ji-1,jj) + zwxtl(ji-1,jj+1)
               zx2tl = zwxtl(ji  ,jj) + zwxtl(ji  ,jj+1)
               pua_tl(ji,jj,jk) = pua_tl(ji,jj,jk) + zfact2 / e1u(ji,jj) * ( zwztl(ji  ,jj-1) * zy1 + zwztl(ji,jj) * zy2 )   &
                                &                  + zfact2 / e1u(ji,jj) * ( zwz(ji  ,jj-1) * zy1tl + zwz(ji,jj) * zy2tl )
               pva_tl(ji,jj,jk) = pva_tl(ji,jj,jk) - zfact2 / e2v(ji,jj) * ( zwztl(ji-1,jj  ) * zx1 + zwztl(ji,jj) * zx2 )   &
                                &                  - zfact2 / e2v(ji,jj) * ( zwz(ji-1,jj  ) * zx1tl + zwz(ji,jj) * zx2tl )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zwz )
      CALL wrk_dealloc( jpi, jpj, zwxtl, zwytl, zwztl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_ene_tan')
      !
   END SUBROUTINE vor_ene_tan


   SUBROUTINE vor_mix_tan( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE vor_mix  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      Mixte formulation : conserves the potential enstrophy of a hori-
      !!      zontally non-divergent flow for (rotzu x uh), the relative vor-
      !!      ticity term and the horizontal kinetic energy for (f x uh), the
      !!      coriolis term. the now trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivatives:
      !!          voru = 1/e1u  mj-1(rotn/e3f) mj-1[ mi(e1v*e3v vn) ]
      !!              +1/e1u  mj-1[ f/e3f          mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1(rotn/e3f) mi-1[ mj(e2u*e3u un) ]
      !!              +1/e2v  mi-1[ f/e3f          mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1(rotn) mj-1[ mi(e1v vn) ]
      !!              +1/e1u  mj-1[ f          mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1(rotn) mi-1[ mj(e2u un) ]
      !!              +1/e2v  mi-1[ f          mj(e2u un) ]
      !!      Add this now trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
      !!             - Save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean timestep index
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zfact1, zua, zcua, zx1, zy1   ! temporary scalars
      REAL(wp) ::   zfact2, zva, zcva, zx2, zy2   !    "         "
      REAL(wp) ::   zuatl, zcuatl, zx1tl, zy1tl   ! temporary scalars
      REAL(wp) ::   zvatl, zcvatl, zx2tl, zy2tl   !    "         "
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwxtl, zwytl, zwztl, zwwtl   ! temporary 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwx, zwy, zwz, zww   ! temporary 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_mix_tan')
      !
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zwz, zww )
      CALL wrk_alloc( jpi, jpj, zwxtl, zwytl, zwztl, zwwtl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_mix : vorticity term: mixed energy/enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      zfact1 = 0.5 * 0.25      ! Local constant initialization
      zfact2 = 0.5 * 0.5

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz, zww )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         ! Relative and planetary potential vorticity and horizontal fluxes
         ! ----------------------------------------------------------------
         IF( ln_sco ) THEN
            IF( ln_dynadv_vec ) THEN
               zww(:,:) = rotn(:,:,jk) / e3f(:,:,jk)
               zwwtl(:,:) = rotn_tl(:,:,jk) / e3f(:,:,jk)
            ELSE
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1   ! vector opt.
                     zww(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                        &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                        &       * 0.5 / ( e1f(ji,jj) * e2f (ji,jj) * e3f(ji,jj,jk) )
                     zwwtl(ji,jj) = (   ( vn_tl(ji+1,jj  ,jk) + vn_tl (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                        &           - ( un_tl(ji  ,jj+1,jk) + un_tl (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                        &       * 0.5 / ( e1f(ji,jj) * e2f (ji,jj) * e3f(ji,jj,jk) )
                  END DO
               END DO
            ENDIF
            zwz(:,:) = ff  (:,:)    / e3f(:,:,jk)
            zwx(:,:) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)
            zwztl(:,:) = 0._wp
            zwxtl(:,:) = e2u(:,:) * e3u(:,:,jk) * un_tl(:,:,jk)
            zwytl(:,:) = e1v(:,:) * e3v(:,:,jk) * vn_tl(:,:,jk)
         ELSE
            IF( ln_dynadv_vec ) THEN
               zww(:,:) = rotn(:,:,jk)
               zwwtl(:,:) = rotn_tl(:,:,jk)
            ELSE
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1   ! vector opt.
                     zww(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                        &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                        &       * 0.5 / ( e1f(ji,jj) * e2f (ji,jj) )
                     zwwtl(ji,jj) = (   ( vn_tl(ji+1,jj  ,jk) + vn_tl (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                        &           - ( un_tl(ji  ,jj+1,jk) + un_tl (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                        &       * 0.5 / ( e1f(ji,jj) * e2f (ji,jj) )
                  END DO
               END DO
            ENDIF
            zwz(:,:) = ff (:,:)
            zwx(:,:) = e2u(:,:) * un(:,:,jk)
            zwy(:,:) = e1v(:,:) * vn(:,:,jk)
            zwztl(:,:) = 0.0_wp
            zwxtl(:,:) = e2u(:,:) * un_tl(:,:,jk)
            zwytl(:,:) = e1v(:,:) * vn_tl(:,:,jk)
         ENDIF

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zy1 = ( zwy(ji,jj-1) + zwy(ji+1,jj-1) ) / e1u(ji,jj)
               zy2 = ( zwy(ji,jj  ) + zwy(ji+1,jj  ) ) / e1u(ji,jj)
               zx1 = ( zwx(ji-1,jj) + zwx(ji-1,jj+1) ) / e2v(ji,jj)
               zx2 = ( zwx(ji  ,jj) + zwx(ji  ,jj+1) ) / e2v(ji,jj)
               zy1tl = ( zwytl(ji,jj-1) + zwytl(ji+1,jj-1) ) / e1u(ji,jj)
               zy2tl = ( zwytl(ji,jj  ) + zwytl(ji+1,jj  ) ) / e1u(ji,jj)
               zx1tl = ( zwxtl(ji-1,jj) + zwxtl(ji-1,jj+1) ) / e2v(ji,jj)
               zx2tl = ( zwxtl(ji  ,jj) + zwxtl(ji  ,jj+1) ) / e2v(ji,jj)
               ! enstrophy conserving formulation for relative vorticity term
               zua = zfact1 * ( zww(ji  ,jj-1) + zww(ji,jj) ) * ( zy1 + zy2 )
               zva =-zfact1 * ( zww(ji-1,jj  ) + zww(ji,jj) ) * ( zx1 + zx2 )
               zuatl = zfact1 * ( zwwtl(ji  ,jj-1) + zwwtl(ji,jj) ) * ( zy1 + zy2 )   &
                 &   + zfact1 * ( zww(ji  ,jj-1) + zww(ji,jj) ) * ( zy1tl + zy2tl )
               zvatl =-zfact1 * ( zwwtl(ji-1,jj  ) + zwwtl(ji,jj) ) * ( zx1 + zx2 )   &
                 &   - zfact1 * ( zww(ji-1,jj  ) + zww(ji,jj) ) * ( zx1tl + zx2tl )
               ! energy conserving formulation for planetary vorticity term
               zcua = zfact2 * ( zwz(ji  ,jj-1) * zy1 + zwz(ji,jj) * zy2 )
               zcva =-zfact2 * ( zwz(ji-1,jj  ) * zx1 + zwz(ji,jj) * zx2 )
               zcuatl = zfact2 * ( zwztl(ji  ,jj-1) * zy1 + zwztl(ji,jj) * zy2 )   &
                  &   + zfact2 * ( zwz(ji  ,jj-1) * zy1tl + zwz(ji,jj) * zy2tl )
               zcvatl =-zfact2 * ( zwztl(ji-1,jj  ) * zx1 + zwztl(ji,jj) * zx2 )   &
                  &    -zfact2 * ( zwz(ji-1,jj  ) * zx1tl + zwz(ji,jj) * zx2tl )
               ! mixed vorticity trend added to the momentum trends
               ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) + zcuatl + zuatl
               va_tl(ji,jj,jk) = va_tl(ji,jj,jk) + zcvatl + zvatl
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zwz, zww )
      CALL wrk_dealloc( jpi, jpj, zwxtl, zwytl, zwztl, zwwtl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_mix_tan')
      !
   END SUBROUTINE vor_mix_tan
   SUBROUTINE vor_ens_tan( kt, kvor, pua_tl, pva_tl )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_ens_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method of the direct routine:
      !!      Trend evaluated using now fields (centered in time)
      !!      and the Sadourny (1975) flux FORM formulation : conserves the
      !!      potential enstrophy of a horizontally non-divergent flow. the
      !!      trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivative:
      !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f ]  mj-1[ mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f ]  mi-1[ mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1[ rotn+f ]  mj-1[ mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1[ rotn+f ]  mi-1[ mj(e2u un) ]
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
      !!             - Save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
         !                                                        ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua_tl    ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva_tl    ! total v-trend
      !!
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zfact1, zuav, zvau   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwx, zwy, zwz         ! temporary 3D workspace
      REAL(wp) ::   zuavtl, zvautl   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwxtl, zwytl, zwztl   ! temporary 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_ens_tan')
      !
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zwz )
      CALL wrk_alloc( jpi, jpj, zwxtl, zwytl, zwztl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_vor_ens_tan : vorticity term: enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF
      ! Local constant initialization
      zfact1 = 0.5 * 0.25

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz  zwxtl, zwytl, zwztl)
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )   ;   zwz(:,:) =                  ff(:,:)      ! planetary vorticity (Coriolis)
         CASE ( 2 )   ;   zwz(:,:) =   rotn(:,:,jk)                ! relative  vorticity
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
               END DO
            END DO
         CASE ( 4 )   ;   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) )    ! total (relative + planetary vorticity)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                                                                       &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )                                               &
                       &       )
               END DO
            END DO
         END SELECT

         IF( ln_sco ) THEN
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwz(ji,jj) = zwz(ji,jj) / e3f(ji,jj,jk)
                  zwx(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk) * un(ji,jj,jk)
                  zwy(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk) * vn(ji,jj,jk)
               END DO
            END DO
         ELSE
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwx(ji,jj) = e2u(ji,jj) * un(ji,jj,jk)
                  zwy(ji,jj) = e1v(ji,jj) * vn(ji,jj,jk)
               END DO
            END DO
         ENDIF

         ! ===================
         ! Tangent counterpart
         ! ===================
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )     ;   zwztl(:,:) =  0.0_wp                      ! planetary vorticity (Coriolis)
         CASE ( 2 ,4)   ;   zwztl(:,:) =  rotn_tl(:,:,jk)             ! relative  vorticity
         CASE ( 3 ,5 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwztl(ji,jj) = (   ( vn_tl(ji+1,jj  ,jk) + vn_tl (ji,jj,jk) )  &
                     &               * ( e2v(ji+1,jj  ) - e2v(ji,jj) )           &
                     &             - ( un_tl(ji  ,jj+1,jk) + un_tl (ji,jj,jk) )  &
                     &               * ( e1u(ji  ,jj+1) - e1u(ji,jj) )           &
                     &           ) * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
               END DO
            END DO
         END SELECT

         IF( ln_sco ) THEN
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwztl(ji,jj) = zwztl(ji,jj) / e3f(ji,jj,jk)
                  zwxtl(ji,jj) = e2u(ji,jj)   * e3u(ji,jj,jk) * un_tl(ji,jj,jk)
                  zwytl(ji,jj) = e1v(ji,jj)   * e3v(ji,jj,jk) * vn_tl(ji,jj,jk)
               END DO
            END DO
         ELSE
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwxtl(ji,jj) = e2u(ji,jj) * un_tl(ji,jj,jk)
                  zwytl(ji,jj) = e1v(ji,jj) * vn_tl(ji,jj,jk)
               END DO
            END DO
         ENDIF

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zuav   = zfact1 / e1u(ji,jj) * ( zwy(  ji  ,jj-1) + zwy(  ji+1,jj-1)   &
                  &                           + zwy(  ji  ,jj  ) + zwy(  ji+1,jj  ) )
               zvau   =-zfact1 / e2v(ji,jj) * ( zwx(  ji-1,jj  ) + zwx(  ji-1,jj+1)   &
                  &                           + zwx(  ji  ,jj  ) + zwx(  ji  ,jj+1) )
               zuavtl = zfact1 / e1u(ji,jj) * ( zwytl(ji  ,jj-1) + zwytl(ji+1,jj-1)   &
                  &                           + zwytl(ji  ,jj  ) + zwytl(ji+1,jj  ) )
               zvautl =-zfact1 / e2v(ji,jj) * ( zwxtl(ji-1,jj  ) + zwxtl(ji-1,jj+1)   &
                  &                           + zwxtl(ji  ,jj  ) + zwxtl(ji  ,jj+1) )
               pua_tl(ji,jj,jk) = pua_tl(ji,jj,jk)                   &
                  &     + zuavtl * ( zwz(  ji,jj-1) + zwz(  ji,jj) ) &
                  &     + zuav   * ( zwztl(ji,jj-1) + zwztl(ji,jj) )
               pva_tl(ji,jj,jk) = pva_tl(ji,jj,jk)                   &
                  &     + zvautl * ( zwz(  ji-1,jj) + zwz(  ji,jj) ) &
                  &     + zvau   * ( zwztl(ji-1,jj) + zwztl(ji,jj) )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zwz )
      CALL wrk_dealloc( jpi, jpj, zwxtl, zwytl, zwztl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_ens_tan')
      !
   END SUBROUTINE vor_ens_tan

   SUBROUTINE vor_een_tan( kt, kvor, pua_tl, pva_tl )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_een_tan  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Arakawa and Lamb (1980) flux form formulation : conserves
      !!      both the horizontal kinetic energy and the potential enstrophy
      !!      when horizontal divergence is zero (see the NEMO documentation)
      !!      Add this trend to the general momentum trend (ua,va).
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References : Arakawa and Lamb 1980, Mon. Wea. Rev., 109, 18-36
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
         !                                                        ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua_tl ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva_tl ! total v-trend
      !!
      INTEGER ::   ji, jj, jk          ! dummy loop indices
      INTEGER :: ierr
      REAL(wp) ::   zfac12, zua, zva   ! temporary scalars
      REAL(wp) ::   zuatl, zvatl       ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwx, zwy, zwz                    ! temporary 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztnw, ztne, ztsw, ztse           ! temporary 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwxtl, zwytl, zwztl              ! temporary 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztnwtl, ztnetl, ztswtl, ztsetl   ! temporary 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:), SAVE ::   ze3f_tl
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('vor_een_tan')
      !
      CALL wrk_alloc( jpi, jpj,      zwx , zwy , zwz        )
      CALL wrk_alloc( jpi, jpj,      ztnw, ztne, ztsw, ztse )
      CALL wrk_alloc( jpi, jpj,      zwxtl , zwytl , zwztl        )
      CALL wrk_alloc( jpi, jpj,      ztnwtl, ztnetl, ztswtl, ztsetl )

      zuatl = 0.0_wp ; zvatl = 0.0_wp
      zwx (:,:) = 0.0_wp ; zwy (:,:) = 0.0_wp ; zwz (:,:) = 0.0_wp
      ztnw(:,:) = 0.0_wp ; ztne(:,:) = 0.0_wp ; ztsw(:,:) = 0.0_wp ; ztse(:,:) = 0.0_wp
      zwxtl (:,:) = 0.0_wp ; zwytl (:,:) = 0.0_wp ; zwztl (:,:) = 0.0_wp
      ztnwtl(:,:) = 0.0_wp ; ztnetl(:,:) = 0.0_wp ; ztswtl(:,:) = 0.0_wp ; ztsetl(:,:) = 0.0_wp
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_een_tam : vorticity term: energy and enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF( .NOT.lk_vvl ) THEN
            ALLOCATE( ze3f_tl(jpi,jpj,jpk) , STAT=ierr )
            IF( lk_mpp    )   CALL mpp_sum ( ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'dyn:vor_een : unable to allocate arrays' )
            ze3f_tl = 0._wp
         ENDIF
      ENDIF

      IF( kt == nit000 .OR. lk_vvl ) THEN      ! reciprocal of e3 at F-point (masked averaging of e3t)
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  ze3f_tl(ji,jj,jk) = ( e3t(ji,jj+1,jk)*tmask(ji,jj+1,jk) + e3t(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)   &
                     &             + e3t(ji,jj  ,jk)*tmask(ji,jj  ,jk) + e3t(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) ) * 0.25_wp
                  IF( ze3f_tl(ji,jj,jk) /= 0.0_wp )   ze3f_tl(ji,jj,jk) = 1.0_wp / ze3f_tl(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( ze3f_tl, 'F', 1._wp )
      ENDIF

      zfac12 = 1.0_wp / 12.0_wp      ! Local constant initialization


!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz, ztnw, ztne, ztsw, ztse )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )
            zwz(:,:) = ff(:,:)      * ze3f_tl(:,:,jk)   ! planetary vorticity (Coriolis)
            zwztl(:,:) = 0.0_wp
         CASE ( 2 )
            zwz(:,:) = rotn(:,:,jk) * ze3f_tl(:,:,jk)   ! relative  vorticity
            zwztl(:,:) = rotn_tl(:,:,jk) * ze3f_tl(:,:,jk)
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) ) * ze3f_tl(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1._wp )
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwztl(ji,jj) = (   ( vn_tl(ji+1,jj  ,jk) + vn_tl (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &         - ( un_tl(ji  ,jj+1,jk) + un_tl (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) ) * ze3f_tl(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwztl, 'F', 1._wp )
         CASE ( 4 )
            zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) * ze3f_tl(:,:,jk) ! total (relative + planetary vorticity)
            zwztl(:,:) = ( rotn_tl(:,:,jk) ) * ze3f_tl(:,:,jk)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                                                                       &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )                                               &
                       &       ) * ze3f_tl(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1._wp )
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwztl(ji,jj) = ( (   ( vn_tl(ji+1,jj  ,jk) + vn_tl (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un_tl(ji  ,jj+1,jk) + un_tl (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )                                               &
                       &       ) * ze3f_tl(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwztl, 'F', 1._wp )
         END SELECT

         zwx(:,:) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
         zwy(:,:) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)

         zwxtl(:,:) = e2u(:,:) * e3u(:,:,jk) * un_tl(:,:,jk)
         zwytl(:,:) = e1v(:,:) * e3v(:,:,jk) * vn_tl(:,:,jk)

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         jj=2
         ztne(1,:)   = 0.0_wp ; ztnw(1,:)   = 0.0_wp ; ztse(1,:)   = 0.0_wp ; ztsw(1,:)   = 0.0_wp
         ztnetl(1,:) = 0.0_wp ; ztnwtl(1,:) = 0.0_wp ; ztsetl(1,:) = 0.0_wp ; ztswtl(1,:) = 0.0_wp
         DO ji = 2, jpi
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )

               ztnetl(ji,jj) = zwztl(ji-1,jj  ) + zwztl(ji  ,jj  ) + zwztl(ji  ,jj-1)
               ztnwtl(ji,jj) = zwztl(ji-1,jj-1) + zwztl(ji-1,jj  ) + zwztl(ji  ,jj  )
               ztsetl(ji,jj) = zwztl(ji  ,jj  ) + zwztl(ji  ,jj-1) + zwztl(ji-1,jj-1)
               ztswtl(ji,jj) = zwztl(ji  ,jj-1) + zwztl(ji-1,jj-1) + zwztl(ji-1,jj  )
         END DO
         DO jj = 3, jpj
            DO ji = 2, jpi   ! vector opt.
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )

               ztnetl(ji,jj) = zwztl(ji-1,jj  ) + zwztl(ji  ,jj  ) + zwztl(ji  ,jj-1)
               ztnwtl(ji,jj) = zwztl(ji-1,jj-1) + zwztl(ji-1,jj  ) + zwztl(ji  ,jj  )
               ztsetl(ji,jj) = zwztl(ji  ,jj  ) + zwztl(ji  ,jj-1) + zwztl(ji-1,jj-1)
               ztswtl(ji,jj) = zwztl(ji  ,jj-1) + zwztl(ji-1,jj-1) + zwztl(ji-1,jj  )
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zuatl = + zfac12 / e1u(ji,jj) * (  ztnetl(ji,jj  ) * zwy(ji  ,jj  ) + ztne(ji,jj  ) * zwytl(ji  ,jj  )  &
                  &                             + ztnwtl(ji+1,jj) * zwy(ji+1,jj  ) + ztnw(ji+1,jj) * zwytl(ji+1,jj  )  &
                  &                             + ztsetl(ji,jj  ) * zwy(ji  ,jj-1) + ztse(ji,jj  ) * zwytl(ji  ,jj-1)    &
                  &                             + ztswtl(ji+1,jj) * zwy(ji+1,jj-1) + ztsw(ji+1,jj) * zwytl(ji+1,jj-1))

               zvatl = - zfac12 / e2v(ji,jj) * (  ztswtl(ji,jj+1) * zwx(ji-1,jj+1) + ztsw(ji,jj+1) * zwxtl(ji-1,jj+1)  &
                  &                           +   ztsetl(ji,jj+1) * zwx(ji  ,jj+1) + ztse(ji,jj+1) * zwxtl(ji  ,jj+1)   &
                  &                           +   ztnwtl(ji,jj  ) * zwx(ji-1,jj  ) + ztnw(ji,jj  ) * zwxtl(ji-1,jj  )   &
                  &                           +   ztnetl(ji,jj  ) * zwx(ji  ,jj  ) + ztne(ji,jj  ) * zwxtl(ji  ,jj  ) )
               pua_tl(ji,jj,jk) = pua_tl(ji,jj,jk) + zuatl
               pva_tl(ji,jj,jk) = pva_tl(ji,jj,jk) + zvatl
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj,      zwx , zwy , zwz                )
      CALL wrk_dealloc( jpi, jpj,      ztnw, ztne, ztsw, ztse         )
      CALL wrk_dealloc( jpi, jpj,      zwxtl , zwytl , zwztl          )
      CALL wrk_dealloc( jpi, jpj,      ztnwtl, ztnetl, ztswtl, ztsetl )
      IF ( kt == nitend .AND. .NOT. lk_vvl ) DEALLOCATE( ze3f_tl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_een_tan')
      !
   END SUBROUTINE vor_een_tan


   SUBROUTINE dyn_vor_adj( kt )
      !!----------------------------------------------------------------------
      !!
      !! ** Purpose of the direct routine:
      !!               compute the lateral ocean tracer physics.
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_vor_adj')
      !
      !                                          ! vorticity term
      SELECT CASE ( nvor )                       ! compute the vorticity trend and add it to the general trend
      !
      CASE ( -1 )                                      ! esopa: test all possibility with control print
         CALL vor_een_adj( kt, ntot, ua_ad, va_ad )
!         CALL vor_mix_adj( kt )
         CALL vor_ens_adj( kt, ntot, ua_ad, va_ad )
!         CALL vor_ene_adj( kt, ntot, ua_ad, va_ad )
         !
      CASE ( 0 )                                       ! energy conserving scheme
         CALL ctl_stop ('vor_ene_adj not available yet')
!         CALL vor_ene_adj( kt, ntot, ua_ad, va_ad )                ! total vorticity
         !
      CASE ( 1 )                                       ! enstrophy conserving scheme
         CALL vor_ens_adj( kt, ntot, ua_ad, va_ad )                ! total vorticity
         !
      CASE ( 2 )                                       ! mixed ene-ens scheme
         CALL ctl_stop ('vor_mix_adj not available yet')
!         CALL vor_mix_adj( kt )                               ! total vorticity (mix=ens-ene)
         !
      CASE ( 3 )                                       ! energy and enstrophy conserving scheme
         CALL vor_een_adj( kt, ntot, ua_ad, va_ad )                ! total vorticity
         !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_vor_adj')
      !
   END SUBROUTINE dyn_vor_adj
   SUBROUTINE vor_ens_adj( kt, kvor, pua_ad, pva_ad )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_ens_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method of the direct routine:
      !!      Trend evaluated using now fields (centered in time)
      !!      and the Sadourny (1975) flux FORM formulation : conserves the
      !!      potential enstrophy of a horizontally non-divergent flow. the
      !!      trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivative:
      !!          voru = 1/e1u  mj-1[ (rotn+f)/e3f ]  mj-1[ mi(e1v*e3v vn) ]
      !!          vorv = 1/e2v  mi-1[ (rotn+f)/e3f ]  mi-1[ mj(e2u*e3u un) ]
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!          voru = 1/e1u  mj-1[ rotn+f ]  mj-1[ mi(e1v vn) ]
      !!          vorv = 1/e2v  mi-1[ rotn+f ]  mi-1[ mj(e2u un) ]
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) arrays with the now vorticity term trend
      !!             - Save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References : Sadourny, r., 1975, j. atmos. sciences, 32, 680-689.
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
         !                                                        ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua_ad    ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva_ad    ! total v-trend
      !!
      INTEGER  ::   ji, jj, jk           ! dummy loop indices
      REAL(wp) ::   zfact1               ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwx, zwy, zwz         ! temporary 3D workspace
      REAL(wp) ::   zuav, zvau       ! temporary scalars
      REAL(wp) ::   zuavad, zvauad   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwxad, zwyad, zwzad   ! temporary 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_ens_adj')
      !
      CALL wrk_alloc( jpi, jpj, zwx, zwy, zwz )
      CALL wrk_alloc( jpi, jpj, zwxad, zwyad, zwzad )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_vor_ens_adj : vorticity term: enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      ! Local constant initialization
      zfact1 = 0.5 * 0.25

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz, zwxad, zwyad, zwzad )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )   ;   zwz(:,:) =                  ff(:,:)      ! planetary vorticity (Coriolis)
         CASE ( 2 )   ;   zwz(:,:) =   rotn(:,:,jk)                ! relative  vorticity
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )  &
                       &       - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) ) )&
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
               END DO
            END DO
         CASE ( 4 )   ;   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) )    ! total (relative + planetary vorticity)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                         &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )  &
                       &       )
               END DO
            END DO
         END SELECT

         IF( ln_sco ) THEN
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwz(ji,jj) = zwz(ji,jj) / e3f(ji,jj,jk)
                  zwx(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk) * un(ji,jj,jk)
                  zwy(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk) * vn(ji,jj,jk)
               END DO
            END DO
         ELSE
            DO jj = 1, jpj                      ! caution: don't use (:,:) for this loop
               DO ji = 1, jpi                   ! it causes optimization problems on NEC in auto-tasking
                  zwx(ji,jj) = e2u(ji,jj) * un(ji,jj,jk)
                  zwy(ji,jj) = e1v(ji,jj) * vn(ji,jj,jk)
               END DO
            END DO
         ENDIF
         ! ===================
         ! Adjoint counterpart
         ! ===================
         zuavad     = 0.0_wp
         zvauad     = 0.0_wp
         zwxad(:,:) = 0.0_wp
         zwyad(:,:) = 0.0_wp
         zwzad(:,:) = 0.0_wp
         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               ! Compute and add the vorticity term trend
               ! ----------------------------------------
               !- Direct counterpart
               zuav = zfact1 / e1u(ji,jj) * ( zwy(ji  ,jj-1) + zwy(ji+1,jj-1)   &
                  &                         + zwy(ji  ,jj  ) + zwy(ji+1,jj  ) )
               zvau =-zfact1 / e2v(ji,jj) * ( zwx(ji-1,jj  ) + zwx(ji-1,jj+1)   &
                  &                         + zwx(ji  ,jj  ) + zwx(ji  ,jj+1) )
               !-
               zvauad = zvauad + pva_ad(ji,jj,jk) * ( zwz(ji-1,jj) + zwz(ji,jj) )
               zwzad(ji-1,jj) = zwzad(ji-1,jj) + pva_ad(ji,jj,jk) * zvau
               zwzad(ji  ,jj) = zwzad(ji  ,jj) + pva_ad(ji,jj,jk) * zvau

               zuavad = zuavad + pua_ad(ji,jj,jk) * ( zwz(ji,jj-1) + zwz(ji,jj) )
               zwzad(ji,jj-1) = zwzad(ji,jj-1) + pua_ad(ji,jj,jk) * zuav
               zwzad(ji,jj  ) = zwzad(ji,jj  ) + pua_ad(ji,jj,jk) * zuav

               zwxad(ji-1,jj  ) = zwxad(ji-1,jj  ) - zvauad * zfact1 / e2v(ji,jj)
               zwxad(ji-1,jj+1) = zwxad(ji-1,jj+1) - zvauad * zfact1 / e2v(ji,jj)
               zwxad(ji  ,jj  ) = zwxad(ji  ,jj  ) - zvauad * zfact1 / e2v(ji,jj)
               zwxad(ji  ,jj+1) = zwxad(ji  ,jj+1) - zvauad * zfact1 / e2v(ji,jj)
               zvauad = 0.0_wp

               zwyad(ji  ,jj-1) = zwyad(ji  ,jj-1) + zuavad * zfact1 / e1u(ji,jj)
               zwyad(ji+1,jj-1) = zwyad(ji+1,jj-1) + zuavad * zfact1 / e1u(ji,jj)
               zwyad(ji  ,jj  ) = zwyad(ji  ,jj  ) + zuavad * zfact1 / e1u(ji,jj)
               zwyad(ji+1,jj  ) = zwyad(ji+1,jj  ) + zuavad * zfact1 / e1u(ji,jj)
               zuavad = 0.0_wp
            END DO
         END DO
         IF( ln_sco ) THEN
            DO jj = jpj, 1, -1     ! caution: don't use (:,:) for this loop
               DO ji = jpi, 1, -1  ! it causes optimization problems on NEC in auto-tasking
                  zwzad(ji,jj) = zwzad(ji,jj) / e3f(ji,jj,jk)
                  un_ad(ji,jj,jk) = un_ad(ji,jj,jk) + zwxad(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk)
                  vn_ad(ji,jj,jk) = vn_ad(ji,jj,jk) + zwyad(ji,jj) * e1v(ji,jj) * e3v(ji,jj,jk)
                  zwxad(ji,jj) = 0.0_wp
                  zwyad(ji,jj) = 0.0_wp
               END DO
            END DO
         ELSE
            DO jj = jpj, 1, -1    ! caution: don't use (:,:) for this loop
               DO ji = jpi, 1, -1 ! it causes optimization problems on NEC in auto-tasking
                  un_ad(ji,jj,jk) = un_ad(ji,jj,jk) + e2u(ji,jj) * zwxad(ji,jj)
                  vn_ad(ji,jj,jk) = vn_ad(ji,jj,jk) + e1v(ji,jj) * zwyad(ji,jj)
                  zwxad(ji,jj) = 0.0_wp
                  zwyad(ji,jj) = 0.0_wp
               END DO
            END DO
         ENDIF
         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )                     ! planetary vorticity (Coriolis)
            zwzad(:,:) =  0.0_wp
         CASE ( 2 ,4)                   ! relative  vorticity
            rotn_ad(:,:,jk) = rotn_ad(:,:,jk) + zwzad(:,:)
            zwzad(:,:) = 0.0_wp
         CASE ( 3 ,5 )                  ! metric term
            DO jj = jpjm1, 1, -1
               DO ji = jpim1, 1, -1  ! vector opt.
                  vn_ad(ji+1,jj,jk) = vn_ad(ji+1,jj,jk)                   &
                     &  + zwzad(ji,jj) * ( e2v(ji+1,jj) - e2v(ji,jj) )    &
                     &                 * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
                  vn_ad(ji  ,jj,jk) = vn_ad(ji  ,jj,jk)                   &
                     &  + zwzad(ji,jj) * ( e2v(ji+1,jj) - e2v(ji,jj) )    &
                     &                 * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
                  un_ad(ji,jj+1,jk) = un_ad(ji,jj+1,jk)                   &
                     &  - zwzad(ji,jj) * ( e1u(ji,jj+1) - e1u(ji,jj) )    &
                     &                 * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
                  un_ad(ji,jj  ,jk) = un_ad(ji,jj  ,jk)                   &
                     &  - zwzad(ji,jj) * ( e1u(ji,jj+1) - e1u(ji,jj) )    &
                     &                 * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )
                  zwzad(ji,jj) = 0.0_wp
               END DO
            END DO
         END SELECT
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj, zwx, zwy, zwz )
      CALL wrk_dealloc( jpi, jpj, zwxad, zwyad, zwzad )
      !
      IF( nn_timing == 1 )  CALL timing_stop('vor_ens_adj')
      !
   END SUBROUTINE vor_ens_adj


   SUBROUTINE vor_een_adj( kt, kvor, pua_ad, pva_ad )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE vor_een_adj  ***
      !!
      !! ** Purpose :   Compute the now total vorticity trend and add it to
      !!      the general trend of the momentum equation.
      !!
      !! ** Method  :   Trend evaluated using now fields (centered in time)
      !!      and the Arakawa and Lamb (19XX) flux form formulation : conserves
      !!      both the horizontal kinetic energy and the potential enstrophy
      !!      when horizontal divergence is zero.
      !!      The trend of the vorticity term is given by:
      !!       * s-coordinate (ln_sco=T), the e3. are inside the derivatives:
      !!       * z-coordinate (default key), e3t=e3u=e3v, the trend becomes:
      !!      Add this trend to the general momentum trend (ua,va):
      !!          (ua,va) = (ua,va) + ( voru , vorv )
      !!
      !! ** Action : - Update (ua,va) with the now vorticity term trend
      !!             - save the trends in (ztrdu,ztrdv) in 2 parts (relative
      !!               and planetary vorticity trends) ('key_trddyn')
      !!
      !! References : Arakawa and Lamb 1980, Mon. Wea. Rev., 109, 18-36
      !!----------------------------------------------------------------------
      INTEGER , INTENT(in   )                         ::   kt     ! ocean time-step index
      INTEGER , INTENT(in   )                         ::   kvor   ! =ncor (planetary) ; =ntot (total) ;
         !                                                        ! =nrvm (relative vorticity or metric)
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pua_ad ! total u-trend
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pva_ad ! total v-trend
      !!
      INTEGER ::   ji, jj, jk          ! dummy loop indices
      INTEGER :: ierr
      REAL(wp) ::   zfac12             ! temporary scalars
      REAL(wp) ::   zuaad, zvaad       ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwx, zwy, zwz                    ! temporary 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztnw, ztne, ztsw, ztse           ! temporary 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   zwxad, zwyad, zwzad              ! temporary 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztnwad, ztnead, ztswad, ztsead   ! temporary 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:), SAVE ::   ze3f_ad
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('vor_een_adj')
      !
      CALL wrk_alloc( jpi, jpj,      zwx , zwy , zwz        )
      CALL wrk_alloc( jpi, jpj,      ztnw, ztne, ztsw, ztse )
      CALL wrk_alloc( jpi, jpj,      zwxad , zwyad , zwzad        )
      CALL wrk_alloc( jpi, jpj,      ztnwad, ztnead, ztswad, ztsead )

      ! local adjoint initailization
      zuaad = 0.0_wp ; zvaad = 0.0_wp
      zwx (:,:) = 0.0_wp ; zwy (:,:) = 0.0_wp ; zwz (:,:) = 0.0_wp
      ztnw(:,:) = 0.0_wp ; ztne(:,:) = 0.0_wp ; ztsw(:,:) = 0.0_wp ; ztse(:,:) = 0.0_wp
      zwxad (:,:) = 0.0_wp ; zwyad (:,:) = 0.0_wp ; zwzad (:,:) = 0.0_wp
      ztnwad(:,:) = 0.0_wp ; ztnead(:,:) = 0.0_wp ; ztswad(:,:) = 0.0_wp ; ztsead(:,:) = 0.0_wp


      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:vor_een_adj : vorticity term: energy and enstrophy conserving scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
         IF( .NOT.lk_vvl ) THEN
            ALLOCATE( ze3f_ad(jpi,jpj,jpk) , STAT=ierr )
            IF( lk_mpp    )   CALL mpp_sum ( ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'dyn:vor_een : unable to allocate arrays' )
            ze3f_ad = 0._wp
         ENDIF
      ENDIF

      IF( kt == nitend .OR. lk_vvl ) THEN      ! reciprocal of e3 at F-point (masked averaging of e3t)
         DO jk = 1, jpk
            DO jj = 1, jpjm1
               DO ji = 1, jpim1
                  ze3f_ad(ji,jj,jk) = ( e3t(ji,jj+1,jk)*tmask(ji,jj+1,jk) + e3t(ji+1,jj+1,jk)*tmask(ji+1,jj+1,jk)   &
                     &             + e3t(ji,jj  ,jk)*tmask(ji,jj  ,jk) + e3t(ji+1,jj  ,jk)*tmask(ji+1,jj  ,jk) ) * 0.25_wp
                  IF( ze3f_ad(ji,jj,jk) /= 0.0_wp )   ze3f_ad(ji,jj,jk) = 1.0_wp / ze3f_ad(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL lbc_lnk( ze3f_ad, 'F', 1._wp )
      ENDIF

      ! Local constant initialization
      zfac12 = 1.0_wp / 12.0_wp

!CDIR PARALLEL DO PRIVATE( zwx, zwy, zwz, ztnw, ztne, ztsw, ztse )
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============

         ! Potential vorticity and horizontal fluxes  (Direct local variables init)
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )
            zwz(:,:) =  ff(:,:) * ze3f_ad(:,:,jk)      ! planetary vorticity (Coriolis)
         CASE ( 2 )   ;   zwz(:,:) =   rotn(:,:,jk) * ze3f_ad(:,:,jk)                ! relative  vorticity
         CASE ( 3 )                                                ! metric term
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )  &
                       &       - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) ) )&
                       &     * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) ) * ze3f_ad(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1._wp )
         CASE ( 4 )   ;   zwz(:,:) = ( rotn(:,:,jk) + ff(:,:) ) * ze3f_ad(:,:,jk)    ! total (relative + planetary vorticity)
         CASE ( 5 )                                                ! total (coriolis + metric)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zwz(ji,jj) = ( ff (ji,jj)                         &
                       &       + (   ( vn(ji+1,jj  ,jk) + vn (ji,jj,jk) ) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )       &
                       &           - ( un(ji  ,jj+1,jk) + un (ji,jj,jk) ) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )   )   &
                       &       * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) )  &
                       &       ) * ze3f_ad(ji,jj,jk)
               END DO
            END DO
            CALL lbc_lnk( zwz, 'F', 1._wp )
         END SELECT

         zwx(:,:) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)
         zwy(:,:) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)

         ! Compute and add the vorticity term trend
         ! ----------------------------------------
         jj=2
         ztne(1,:)   = 0.0_wp ; ztnw(1,:)   = 0.0_wp ; ztse(1,:)   = 0.0_wp ; ztsw(1,:)   = 0.0_wp
         DO ji = 2, jpi
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
         END DO
         DO jj = 3, jpj
            DO ji = 2, jpi   ! vector opt.
               ztne(ji,jj) = zwz(ji-1,jj  ) + zwz(ji  ,jj  ) + zwz(ji  ,jj-1)
               ztnw(ji,jj) = zwz(ji-1,jj-1) + zwz(ji-1,jj  ) + zwz(ji  ,jj  )
               ztse(ji,jj) = zwz(ji  ,jj  ) + zwz(ji  ,jj-1) + zwz(ji-1,jj-1)
               ztsw(ji,jj) = zwz(ji  ,jj-1) + zwz(ji-1,jj-1) + zwz(ji-1,jj  )
            END DO
         END DO

      ! ===================
      ! Adjoint counterpart
      ! ===================

         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               zuaad = zuaad + pua_ad(ji,jj,jk)
               zvaad = zvaad + pva_ad(ji,jj,jk)

               zvaad = - zvaad * zfac12 / e2v(ji,jj)
               ztswad(ji  ,jj+1) = ztswad(ji  ,jj+1) + zvaad * zwx (ji-1,jj+1)
               zwxad (ji-1,jj+1) = zwxad (ji-1,jj+1) + zvaad * ztsw(ji  ,jj+1)
               ztsead(ji  ,jj+1) = ztsead(ji  ,jj+1) + zvaad * zwx (ji  ,jj+1)
               zwxad (ji  ,jj+1) = zwxad (ji  ,jj+1) + zvaad * ztse(ji  ,jj+1)
               ztnwad(ji  ,jj  ) = ztnwad(ji  ,jj  ) + zvaad * zwx (ji-1,jj  )
               zwxad (ji-1,jj  ) = zwxad (ji-1,jj  ) + zvaad * ztnw(ji  ,jj  )
               ztnead(ji  ,jj  ) = ztnead(ji  ,jj  ) + zvaad * zwx (ji  ,jj  )
               zwxad (ji  ,jj  ) = zwxad (ji  ,jj  ) + zvaad * ztne(ji  ,jj  )
               zvaad = 0.0_wp

               zuaad = zuaad * zfac12 / e1u(ji,jj)
               ztnead(ji  ,jj  ) = ztnead(ji  ,jj  ) + zuaad * zwy (ji  ,jj  )
               zwyad (ji  ,jj  ) = zwyad (ji  ,jj  ) + zuaad * ztne(ji  ,jj  )
               ztnwad(ji+1,jj  ) = ztnwad(ji+1,jj  ) + zuaad * zwy (ji+1,jj  )
               zwyad (ji+1,jj  ) = zwyad (ji+1,jj  ) + zuaad * ztnw(ji+1,jj  )
               ztsead(ji  ,jj  ) = ztsead(ji  ,jj  ) + zuaad * zwy (ji  ,jj-1)
               zwyad (ji  ,jj-1) = zwyad (ji  ,jj-1) + zuaad * ztse(ji  ,jj  )
               ztswad(ji+1,jj  ) = ztswad(ji+1,jj  ) + zuaad * zwy (ji+1,jj-1)
               zwyad (ji+1,jj-1) = zwyad (ji+1,jj-1) + zuaad * ztsw(ji+1,jj  )
               zuaad = 0.0_wp
            END DO
         END DO
         DO jj = jpj, 3, -1
            DO ji = jpi, 2, -1   ! vector opt.
               zwzad (ji  ,jj-1) = zwzad(ji  ,jj-1) + ztswad(ji,jj)
               zwzad (ji-1,jj-1) = zwzad(ji-1,jj-1) + ztswad(ji,jj)
               zwzad (ji-1,jj  ) = zwzad(ji-1,jj  ) + ztswad(ji,jj)
               ztswad(ji  ,jj  ) = 0.0_wp
               zwzad (ji  ,jj  ) = zwzad(ji  ,jj  ) + ztsead(ji,jj)
               zwzad (ji  ,jj-1) = zwzad(ji  ,jj-1) + ztsead(ji,jj)
               zwzad (ji-1,jj-1) = zwzad(ji-1,jj-1) + ztsead(ji,jj)
               ztsead(ji,jj)     = 0.0_wp
               zwzad (ji-1,jj-1) = zwzad(ji-1,jj-1) + ztnwad(ji,jj)
               zwzad (ji-1,jj  ) = zwzad(ji-1,jj  ) + ztnwad(ji,jj)
               zwzad (ji  ,jj  ) = zwzad(ji  ,jj  ) + ztnwad(ji,jj)
               ztnwad(ji  ,jj  ) = 0.0_wp
               zwzad (ji-1,jj  ) = zwzad(ji-1,jj  ) + ztnead(ji,jj)
               zwzad (ji  ,jj  ) = zwzad(ji  ,jj  ) + ztnead(ji,jj)
               zwzad (ji  ,jj-1) = zwzad(ji  ,jj-1) + ztnead(ji,jj)
               ztnead(ji,jj)     = 0.0_wp
            END DO
         END DO
         jj=2
         DO ji = jpi, 2, -1
               zwzad (ji  ,jj-1) = zwzad(ji  ,jj-1) + ztswad(ji,jj)
               zwzad (ji-1,jj-1) = zwzad(ji-1,jj-1) + ztswad(ji,jj)
               zwzad (ji-1,jj  ) = zwzad(ji-1,jj  ) + ztswad(ji,jj)
               ztswad(ji,jj) = 0.0_wp
               zwzad (ji  ,jj  ) = zwzad(ji  ,jj  ) + ztsead(ji,jj)
               zwzad (ji  ,jj-1) = zwzad(ji  ,jj-1) + ztsead(ji,jj)
               zwzad (ji-1,jj-1) = zwzad(ji-1,jj-1) + ztsead(ji,jj)
               ztsead(ji  ,jj  ) = 0.0_wp
               zwzad (ji-1,jj-1) = zwzad(ji-1,jj-1) + ztnwad(ji,jj)
               zwzad (ji-1,jj  ) = zwzad(ji-1,jj  ) + ztnwad(ji,jj)
               zwzad (ji  ,jj  ) = zwzad(ji  ,jj  ) + ztnwad(ji,jj)
               ztnwad(ji  ,jj ) = 0.0_wp
               zwzad (ji-1,jj  ) = zwzad(ji-1,jj  ) + ztnead(ji,jj)
               zwzad (ji  ,jj  ) = zwzad(ji  ,jj  ) + ztnead(ji,jj)
               zwzad (ji  ,jj-1) = zwzad(ji  ,jj-1) + ztnead(ji,jj)
               ztnead(ji  ,jj  ) = 0.0_wp
         END DO
         ztnead(1,:)   = 0.0_wp ; ztnwad(1,:)   = 0.0_wp
         ztsead(1,:)   = 0.0_wp ; ztswad(1,:)   = 0.0_wp

         vn_ad(:,:,jk) = vn_ad(:,:,jk) + zwyad(:,:) * e1v(:,:) * e3v(:,:,jk)
         un_ad(:,:,jk) = un_ad(:,:,jk) + zwxad(:,:) * e2u(:,:) * e3u(:,:,jk)
         zwyad(:,:)    = 0.0_wp
         zwxad(:,:)    = 0.0_wp

         ! Potential vorticity and horizontal fluxes
         ! -----------------------------------------
         SELECT CASE( kvor )      ! vorticity considered
         CASE ( 1 )
            zwzad(:,:) = 0.0_wp
         CASE ( 2 )
            rotn_ad(:,:,jk) = rotn_ad(:,:,jk) +  zwzad(:,:) * ze3f_ad(:,:,jk)
            zwzad(:,:)      = 0.0_wp
         CASE ( 3 )
            CALL lbc_lnk_adj( zwzad, 'F', 1._wp )                                                ! metric term
            DO jj = jpjm1, 1, -1
               DO ji = jpim1, 1, -1   ! vector opt.
                  zwzad(ji  ,jj     ) = zwzad(ji,jj) * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) ) * ze3f_ad(ji,jj,jk)
                  vn_ad(ji+1,jj  ,jk) =  vn_ad(ji+1,jj  ,jk) + zwzad(ji,jj) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )
                  vn_ad(ji  ,jj  ,jk) =  vn_ad(ji  ,jj  ,jk) + zwzad(ji,jj) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )
                  un_ad(ji  ,jj+1,jk) =  un_ad(ji  ,jj+1,jk) - zwzad(ji,jj) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )
                  un_ad(ji  ,jj  ,jk) =  un_ad(ji  ,jj  ,jk) - zwzad(ji,jj) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )
                  zwzad(ji  ,jj     ) = 0.0_wp
               END DO
            END DO
         CASE ( 4 )
            rotn_ad(:,:,jk) = rotn_ad(:,:,jk) + zwzad(:,:) * ze3f_ad(:,:,jk)
            zwzad(:,:) = 0.0_wp
         CASE ( 5 )
            CALL lbc_lnk_adj( zwzad, 'F', 1._wp )                                               ! total (coriolis + metric)
            DO jj = jpjm1, 1, -1
               DO ji = jpim1, 1, -1   ! vector opt.
                  zwzad(ji  ,jj     ) = zwzad(ji,jj) * 0.5 / ( e1f(ji,jj) * e2f(ji,jj) ) * ze3f_ad(ji,jj,jk)
                  vn_ad(ji+1,jj  ,jk) = vn_ad(ji+1,jj  ,jk) + zwzad(ji,jj) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )
                  vn_tl(ji  ,jj  ,jk) = vn_tl(ji  ,jj  ,jk) + zwzad(ji,jj) * ( e2v(ji+1,jj  ) - e2v(ji,jj) )
                  un_ad(ji  ,jj+1,jk) = un_ad(ji  ,jj+1,jk) - zwzad(ji,jj) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )
                  un_ad(ji  ,jj  ,jk) = un_ad(ji  ,jj  ,jk) - zwzad(ji,jj) * ( e1u(ji  ,jj+1) - e1u(ji,jj) )
                  zwzad(ji  ,jj     ) = 0.0_wp
               END DO
            END DO
         END SELECT
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj,      zwx , zwy , zwz        )
      CALL wrk_dealloc( jpi, jpj,      ztnw, ztne, ztsw, ztse )
      CALL wrk_dealloc( jpi, jpj,      zwxad , zwyad , zwzad        )
      CALL wrk_dealloc( jpi, jpj,      ztnwad, ztnead, ztswad, ztsead )

      IF ( kt == nit000 .AND. .NOT. lk_vvl ) DEALLOCATE( ze3f_ad )

      IF( nn_timing == 1 )  CALL timing_stop('vor_een_adj')
      !
   END SUBROUTINE vor_een_adj

   SUBROUTINE dyn_vor_init_tam
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE vor_ctl_tam  ***
      !!
      !! ** Purpose :   Control the consistency between cpp options for
      !!      tracer advection schemes
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio          ! temporary integer
      INTEGER ::   ji, jj, jk      ! dummy loop indices
      NAMELIST/namdyn_vor/ ln_dynvor_ens, ln_dynvor_ene, ln_dynvor_mix, ln_dynvor_een
      !!----------------------------------------------------------------------

      REWIND ( numnam )               ! Read Namelist namdyn_vor : Vorticity scheme options
      READ   ( numnam, namdyn_vor )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn:vor_init_tam : vorticity term : read namelist and control the consistency'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '        Namelist namdyn_vor : choice of the vorticity term scheme'
         WRITE(numout,*) '           energy    conserving scheme                ln_dynvor_ene = ', ln_dynvor_ene
         WRITE(numout,*) '           enstrophy conserving scheme                ln_dynvor_ens = ', ln_dynvor_ens
         WRITE(numout,*) '           mixed enstrophy/energy conserving scheme   ln_dynvor_mix = ', ln_dynvor_mix
         WRITE(numout,*) '           enstrophy and energy conserving scheme     ln_dynvor_een = ', ln_dynvor_een
      ENDIF
      ! If energy, enstrophy or mixed advection of momentum in vector form change the value for masks
      ! at angles with three ocean points and one land point
      IF( ln_vorlat .AND. ( ln_dynvor_ene .OR. ln_dynvor_ens .OR. ln_dynvor_mix ) ) THEN
         DO jk = 1, jpk
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  IF( tmask(ji,jj,jk)+tmask(ji+1,jj,jk)+tmask(ji,jj+1,jk)+tmask(ji+1,jj+1,jk) == 3._wp ) &
                      fmask(ji,jj,jk) = 1._wp
               END DO
            END DO
         END DO
          !
          CALL lbc_lnk( fmask, 'F', 1._wp )      ! Lateral boundary conditions on fmask
          !
      ENDIF
      !
      ioptio = 0                     ! Control of vorticity scheme options
      IF( ln_dynvor_ene )   ioptio = ioptio + 1
      IF( ln_dynvor_ens )   ioptio = ioptio + 1
      IF( ln_dynvor_mix )   ioptio = ioptio + 1
      IF( ln_dynvor_een )   ioptio = ioptio + 1
      IF( lk_esopa      )   ioptio =          1

      IF( ioptio /= 1 ) CALL ctl_stop( ' use ONE and ONLY one vorticity scheme' )

      !                              ! Set nvor (type of scheme for vorticity)
      IF( ln_dynvor_ene )   nvor =  0
      IF( ln_dynvor_ens )   nvor =  1
      IF( ln_dynvor_mix )   nvor =  2
      IF( ln_dynvor_een )   nvor =  3
      IF( lk_esopa      )   nvor = -1

      !                              ! Set ncor, nrvm, ntot (type of vorticity)
      IF(lwp) WRITE(numout,*)
      ncor = 1
      IF( ln_dynadv_vec ) THEN
         IF(lwp) WRITE(numout,*) '         Vector form advection : vorticity = Coriolis + relative vorticity'
         nrvm = 2
         ntot = 4
      ELSE
         IF(lwp) WRITE(numout,*) '         Flux form advection   : vorticity = Coriolis + metric term'
         nrvm = 3
         ntot = 5
      ENDIF

      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( nvor ==  0 )   WRITE(numout,*) '         vorticity scheme : energy conserving scheme'
         IF( nvor ==  1 )   WRITE(numout,*) '         vorticity scheme : enstrophy conserving scheme'
         IF( nvor ==  2 )   WRITE(numout,*) '         vorticity scheme : mixed enstrophy/energy conserving scheme'
         IF( nvor ==  3 )   WRITE(numout,*) '         vorticity scheme : energy and enstrophy conserving scheme'
         IF( nvor == -1 )   WRITE(numout,*) '         esopa test: use all lateral physics options'
      ENDIF
      !
   END SUBROUTINE dyn_vor_init_tam

   SUBROUTINE dyn_vor_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_adv_adj_tst ***
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
      !! ** Action  : Separate tests are applied for the following dx and dy:
      !!
      !!              1) dx = ( SSH ) and dy = ( SSH )
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
         & jt
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zun_tlin,     & ! Tangent input:  now   u-velocity
         & zvn_tlin,     & ! Tangent input:  now   v-velocity
         & zrotn_tlin,   & ! Tangent input:  now   rot
         & zun_adout,    & ! Adjoint output: now   u-velocity
         & zvn_adout,    & ! Adjoint output: now   v-velocity
         & zrotn_adout,  & ! Adjoint output: now   rot
         & zua_adout,    & ! Tangent output: after u-velocity
         & zva_adout,    & ! Tangent output: after v-velocity
         & zua_tlin,     & ! Tangent output: after u-velocity
         & zva_tlin,     & ! Tangent output: after v-velocity
         & zua_tlout,    & ! Tangent output: after u-velocity
         & zva_tlout,    & ! Tangent output: after v-velocity
         & zua_adin,     & ! Tangent output: after u-velocity
         & zva_adin,     & ! Tangent output: after v-velocity
         & zau,          & ! 3D random field for rotn
         & zav,          & ! 3D random field for rotn
         & znu,          & ! 3D random field for u
         & znv             ! 3D random field for v
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4,       &
         & zsp2_5
      CHARACTER(LEN=14) :: cl_name

      ! Allocate memory

      ALLOCATE( &
         & zun_tlin(jpi,jpj,jpk),     &
         & zvn_tlin(jpi,jpj,jpk),     &
         & zrotn_tlin(jpi,jpj,jpk),   &
         & zun_adout(jpi,jpj,jpk),    &
         & zvn_adout(jpi,jpj,jpk),    &
         & zrotn_adout(jpi,jpj,jpk),  &
         & zua_adout(jpi,jpj,jpk),    &
         & zva_adout(jpi,jpj,jpk),    &
         & zua_tlin(jpi,jpj,jpk),     &
         & zva_tlin(jpi,jpj,jpk),     &
         & zua_tlout(jpi,jpj,jpk),    &
         & zva_tlout(jpi,jpj,jpk),    &
         & zua_adin(jpi,jpj,jpk),     &
         & zva_adin(jpi,jpj,jpk),     &
         & zau(jpi,jpj,jpk),          &
         & zav(jpi,jpj,jpk),          &
         & znu(jpi,jpj,jpk),          &
         & znv(jpi,jpj,jpk)           &
         & )

      ! init ntot parameter
      CALL dyn_vor_init_tam          ! initialisation & control of options

      DO jt = 1, 2
         IF (jt == 1) nvor=1 ! enstrophy conserving scheme
         IF (jt == 2) nvor=3 ! energy and enstrophy conserving scheme

         ! Initialize rotn
!AV: it calls cla_init the tries to reallocate already allocated arrays...nasty
!AV         CALL div_cur ( nit000 )

         !==================================================================
         ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
         !    dy = ( hdivb_tl, hdivn_tl )
         !==================================================================

         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------

         zun_tlin(:,:,:) = 0.0_wp
         zvn_tlin(:,:,:) = 0.0_wp
         zrotn_tlin(:,:,:) = 0.0_wp
         zun_adout(:,:,:) = 0.0_wp
         zvn_adout(:,:,:) = 0.0_wp
         zrotn_adout(:,:,:) = 0.0_wp
         zua_tlout(:,:,:) = 0.0_wp
         zva_tlout(:,:,:) = 0.0_wp
         zua_adin(:,:,:) = 0.0_wp
         zva_adin(:,:,:) = 0.0_wp
         zua_adout(:,:,:) = 0.0_wp
         zva_adout(:,:,:) = 0.0_wp
         zua_tlin(:,:,:) = 0.0_wp
         zva_tlin(:,:,:) = 0.0_wp
         znu(:,:,:) = 0.0_wp
         znv(:,:,:) = 0.0_wp
         zau(:,:,:) = 0.0_wp
         zav(:,:,:) = 0.0_wp


         un_tl(:,:,:) = 0.0_wp
         vn_tl(:,:,:) = 0.0_wp
         ua_tl(:,:,:) = 0.0_wp
         va_tl(:,:,:) = 0.0_wp
         un_ad(:,:,:) = 0.0_wp
         vn_ad(:,:,:) = 0.0_wp
         ua_ad(:,:,:) = 0.0_wp
         va_ad(:,:,:) = 0.0_wp
         rotn_tl(:,:,:) = 0.0_wp
         rotn_ad(:,:,:) = 0.0_wp

         !--------------------------------------------------------------------
         ! Initialize the tangent input with random noise: dx
         !--------------------------------------------------------------------

         CALL grid_random(  znu, 'U', 0.0_wp, stdu )
         CALL grid_random(  znv, 'V', 0.0_wp, stdv )
         CALL grid_random(  zau, 'U', 0.0_wp, stdu )
         CALL grid_random(  zav, 'V', 0.0_wp, stdv )

         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zun_tlin(ji,jj,jk) = znu(ji,jj,jk)
                  zvn_tlin(ji,jj,jk) = znv(ji,jj,jk)
                  zua_tlin(ji,jj,jk) = zau(ji,jj,jk)
                  zva_tlin(ji,jj,jk) = zav(ji,jj,jk)
               END DO
            END DO
         END DO
         un_tl(:,:,:) = zun_tlin(:,:,:)
         vn_tl(:,:,:) = zvn_tlin(:,:,:)
         ua_tl(:,:,:) = zua_tlin(:,:,:)
         va_tl(:,:,:) = zva_tlin(:,:,:)

         ! initialize rotn_tl with noise
         CALL div_cur_tan ( nit000 )

         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zrotn_tlin(ji,jj,jk) = rotn_tl(ji,jj,jk)
               END DO
            END DO
         END DO
         rotn_tl(:,:,:) = zrotn_tlin(:,:,:)


         IF (nvor == 1 )  CALL vor_ens_tan( nit000, ntot, ua_tl, va_tl )
         IF (nvor == 3 )  CALL vor_een_tan( nit000, ntot, ua_tl, va_tl )
         zua_tlout(:,:,:) = ua_tl(:,:,:)
         zva_tlout(:,:,:) = va_tl(:,:,:)

         !--------------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !--------------------------------------------------------------------

         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zua_adin(ji,jj,jk) = zua_tlout(ji,jj,jk) &
                       &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
                       &               * umask(ji,jj,jk)
                  zva_adin(ji,jj,jk) = zva_tlout(ji,jj,jk) &
                       &               * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) &
                       &               * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !--------------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !--------------------------------------------------------------------

         zsp1_1 = DOT_PRODUCT( zua_tlout, zua_adin )
         zsp1_2 = DOT_PRODUCT( zva_tlout, zva_adin )
         zsp1   = zsp1_1 + zsp1_2

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------

         ua_ad(:,:,:) = zua_adin(:,:,:)
         va_ad(:,:,:) = zva_adin(:,:,:)


         IF (nvor == 1 )  CALL vor_ens_adj( nitend, ntot, ua_ad, va_ad )
         IF (nvor == 3 )  CALL vor_een_adj( nitend, ntot, ua_ad, va_ad )
         zun_adout(:,:,:)   = un_ad(:,:,:)
         zvn_adout(:,:,:)   = vn_ad(:,:,:)
         zrotn_adout(:,:,:) = rotn_ad(:,:,:)
         zua_adout(:,:,:)   = ua_ad(:,:,:)
         zva_adout(:,:,:)   = va_ad(:,:,:)

         zsp2_1 = DOT_PRODUCT( zun_tlin, zun_adout )
         zsp2_2 = DOT_PRODUCT( zvn_tlin, zvn_adout )
         zsp2_3 = DOT_PRODUCT( zrotn_tlin, zrotn_adout )
         zsp2_4 = DOT_PRODUCT( zua_tlin, zua_adout )
         zsp2_5 = DOT_PRODUCT( zva_tlin, zva_adout )
         zsp2   = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5

         ! Compare the scalar products

         ! 14 char:'12345678901234'
         IF (nvor == 1 )  cl_name = 'dynvor_adj ens'
         IF (nvor == 3 )  cl_name = 'dynvor_adj een'

         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )
      END DO

      DEALLOCATE( &
         & zun_tlin,     &
         & zvn_tlin,     &
         & zrotn_tlin,   &
         & zun_adout,    &
         & zvn_adout,    &
         & zrotn_adout,  &
         & zua_adout,    &
         & zva_adout,    &
         & zua_tlin,     &
         & zva_tlin,     &
         & zua_tlout,    &
         & zva_tlout,    &
         & zua_adin,     &
         & zva_adin,     &
         & zau,          &
         & zav,          &
         & znu,          &
         & znv           &
         & )
   END SUBROUTINE dyn_vor_adj_tst
   !!=============================================================================
END MODULE dynvor_tam
