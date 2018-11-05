MODULE dynspg_flt_tam
   !!----------------------------------------------------------------------
   !!    This software is governed by the CeCILL licence (Version 2)
   !!----------------------------------------------------------------------
   !!======================================================================
   !!  ***  MODULE  dynspg_flt_tam : TANGENT/ADJOINT OF MODULE dynspg_flt  ***
   !!
   !! Ocean dynamics:  surface pressure gradient trend
   !!
   !!======================================================================
   !! History of the direct module:
   !! History    OPA  !  1998-05  (G. Roullet)  free surface
   !!                 !  1998-10  (G. Madec, M. Imbard)  release 8.2
   !!   NEMO     O.1  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!            1.0  !  2004-08  (C. Talandier) New trends organization
   !!             -   !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.0  !  2006-07  (S. Masson)  distributed restart using iom
   !!             -   !  2006-08  (J.Chanut, A.Sellar) Calls to BDY routines.
   !!            3.2  !  2009-03  (G. Madec, M. Leclair, R. Benshila) introduce sshwzv module
   !! History of the TAM module:
   !!            9.0  ! 2008-08  (A. Vidard) skeleton
   !!             -   ! 2008-08  (A. Weaver) original version
   !!   NEMO     3.2  ! 2010-04  (F. Vigilant) converison to 3.2
   !!   NEMO     3.4  ! 2012-07  (P.-A. Bouttier) converison to 3.4
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_dynspg_flt'                              filtered free surface
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   dyn_spg_flt_tan  : update the momentum trend with the surface pressure
   !!                      gradient in the filtered free surface case
   !!                      (tangent routine)
   !!   dyn_spg_flt_adj  : update the momentum trend with the surface pressure
   !!                      gradient in the filtered free surface case
   !!                      (adjoint routine)
   !!   dyn_spg_flt_adj_tst : Test of the adjoint routine
   !!----------------------------------------------------------------------
   USE prtctl                  ! Print control
   USE par_oce
   USE in_out_manager
   USE phycst
   USE lib_mpp
   USE dom_oce
   USE solver
   USE dynspg_flt
   USE sol_oce
   USE oce_tam
   USE sbc_oce_tam
   USE sol_oce
   USE sol_oce_tam
   USE solsor_tam
   USE lbclnk
   USE lbclnk_tam
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE dynadv
   USE cla_tam
   USE tstool_tam
   USE wrk_nemo        ! Memory Allocation
   USE lib_fortran
   USE timing
   USE iom



   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_flt_tan,     & ! routine called by step_tan.F90
      &   dyn_spg_flt_adj,     & ! routine called by step_adj.F90
      &   dyn_spg_flt_adj_tst    ! routine called by the tst.F90
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

   SUBROUTINE dyn_spg_flt_tan( kt, kindic )
      !!---------------------------------------------------------------------
      !!                  ***  routine dyn_spg_flt_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now trend due to the surface pressure
      !!      gradient in case of filtered free surface formulation  and add
      !!      it to the general trend of momentum equation.
      !!
      !! ** Method of the direct routine:
      !!      Filtered free surface formulation. The surface
      !!      pressure gradient is given by:
      !!         spgu = 1/rau0 d/dx(ps) =  1/e1u di( sshn + btda )
      !!         spgv = 1/rau0 d/dy(ps) =  1/e2v dj( sshn + btda )
      !!      where sshn is the free surface elevation and btda is the after
      !!      time derivative of the free surface elevation
      !!       -1- evaluate the surface presure trend (including the addi-
      !!      tional force) in three steps:
      !!        a- compute the right hand side of the elliptic equation:
      !!            gcb = 1/(e1t e2t) [ di(e2u spgu) + dj(e1v spgv) ]
      !!         where (spgu,spgv) are given by:
      !!            spgu = vertical sum[ e3u (ub+ 2 rdt ua ) ]
      !!                 - grav 2 rdt hu /e1u di[sshn + emp]
      !!            spgv = vertical sum[ e3v (vb+ 2 rdt va) ]
      !!                 - grav 2 rdt hv /e2v dj[sshn + emp]
      !!         and define the first guess from previous computation :
      !!            zbtd = btda
      !!            btda = 2 zbtd - btdb
      !!            btdb = zbtd
      !!        b- compute the relative accuracy to be reached by the
      !!         iterative solver
      !!        c- apply the solver by a call to sol... routine
      !!       -2- compute and add the free surface pressure gradient inclu-
      !!      ding the additional force used to stabilize the equation.
      !!
      !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
      !!
      !! References : Roullet and Madec 1999, JGR.
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( IN  ) :: &
         &  kt         ! ocean time-step index
      INTEGER, INTENT( OUT ) :: &
         &  kindic     ! solver convergence flag (<0 if not converge)

      !! * Local declarations
      INTEGER  :: &
         & ji,    &     ! dummy loop indices
         & jj,    &
         & jk
      REAL(wp) ::  &
         & z2dt,   & ! temporary scalars
         & z2dtg,  &
         & zgcb,   &
         & zbtd,   &
         & ztdgu,  &
         & ztdgv
      !!----------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_flt_tan')
      !
      IF( kt == nit000 ) THEN

         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_flt_tan : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   (free surface constant volume case)'
         ! set to zero free surface specific arrays
         spgu_tl(:,:) = 0.0_wp                    ! surface pressure gradient (i-direction)
         spgv_tl(:,:) = 0.0_wp                    ! surface pressure gradient (j-direction)
         ! Reinitialize the solver arrays
         gcxb_tl(:,:) = 0.e0
         gcx_tl (:,:) = 0.e0
         CALL sol_mat( nit000 )
      ENDIF
      ! Local constant initialization
      z2dt = 2. * rdt                                             ! time step: leap-frog
      IF( neuler == 0 .AND. kt == nit000   )   z2dt = rdt         ! time step: Euler if restart from rest
      IF( neuler == 0 .AND. kt == nit000+1 )   CALL sol_mat( kt )
      z2dtg  = grav * z2dt
      ! Evaluate the masked next velocity (effect of the additional force not included)
      IF( lk_vvl ) THEN          ! variable volume  (surface pressure gradient already included in dyn_hpg)
         CALL ctl_stop('key_vvl is not implemented in TAM yet')
         !
      ELSE                       ! fixed volume  (add the surface pressure gradient + unweighted time stepping)
         !
         DO jj = 2, jpjm1              ! Surface pressure gradient (now)
            DO ji = 2, jpim1   ! vector opt.
               spgu_tl(ji,jj) = - grav * ( sshn_tl(ji+1,jj) - sshn_tl(ji,jj) ) / e1u(ji,jj)
               spgv_tl(ji,jj) = - grav * ( sshn_tl(ji,jj+1) - sshn_tl(ji,jj) ) / e2v(ji,jj)
            END DO
         END DO

         DO jk = 1, jpkm1              ! unweighted time stepping
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ua_tl(ji,jj,jk) = (  ub_tl(ji,jj,jk) + z2dt * ( ua_tl(ji,jj,jk) + spgu_tl(ji,jj) )  ) * umask(ji,jj,jk)
                  va_tl(ji,jj,jk) = (  vb_tl(ji,jj,jk) + z2dt * ( va_tl(ji,jj,jk) + spgv_tl(ji,jj) )  ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !
      ENDIF
      IF( nn_cla == 1 )   CALL cla_dynspg_tan( kt )      ! Cross Land Advection (update (ua,va))

      ! compute the next vertically averaged velocity (effect of the additional force not included)
      ! ---------------------------------------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            spgu_tl(ji,jj) = 0.0_wp
            spgv_tl(ji,jj) = 0.0_wp
         END DO
      END DO

      ! vertical sum
!CDIR NOLOOPCHG
      IF( lk_vopt_loop ) THEN          ! vector opt., forced unroll
         DO jk = 1, jpkm1
            DO ji = 1, jpij
               spgu_tl(ji,1) = spgu_tl(ji,1) + e3u(ji,1,jk) * ua_tl(ji,1,jk)
               spgv_tl(ji,1) = spgv_tl(ji,1) + e3v(ji,1,jk) * va_tl(ji,1,jk)
            END DO
         END DO
      ELSE                        ! No  vector opt.
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  spgu_tl(ji,jj) = spgu_tl(ji,jj) + e3u(ji,jj,jk) * ua_tl(ji,jj,jk)
                  spgv_tl(ji,jj) = spgv_tl(ji,jj) + e3v(ji,jj,jk) * va_tl(ji,jj,jk)
               END DO
            END DO
         END DO
      ENDIF

      ! transport: multiplied by the horizontal scale factor
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            spgu_tl(ji,jj) = spgu_tl(ji,jj) * e2u(ji,jj)
            spgv_tl(ji,jj) = spgv_tl(ji,jj) * e1v(ji,jj)
         END DO
      END DO

      CALL lbc_lnk( spgu_tl, 'U', -1.0_wp )       ! lateral boundary conditions
      CALL lbc_lnk( spgv_tl, 'V', -1.0_wp )

      IF( lk_vvl ) CALL ctl_stop( 'dyn_spg_flt_tan: lk_vvl is not available' )

      ! Right hand side of the elliptic equation and first guess
      ! -----------------------------------------------------------
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            ! Divergence of the after vertically averaged velocity
            zgcb =  spgu_tl(ji,jj) - spgu_tl(ji-1,jj)   &
               &  + spgv_tl(ji,jj) - spgv_tl(ji,jj-1)
            gcb_tl(ji,jj) = gcdprc(ji,jj) * zgcb
            ! First guess of the after barotropic transport divergence
            zbtd = gcx_tl(ji,jj)
            gcx_tl (ji,jj) = 2.0_wp * zbtd - gcxb_tl(ji,jj)
            gcxb_tl(ji,jj) =          zbtd
         END DO
      END DO

      ! apply the lateral boundary conditions
      IF( nn_solv == 2 .AND. MAX( jpr2di, jpr2dj ) > 0 ) CALL lbc_lnk_e( gcb_tl, c_solver_pt, 1.0_wp )

      ! Relative precision
      ! ------------------

      rnorme = GLOB_SUM( gcb_tl(1:jpi,1:jpj) * gcdmat(1:jpi,1:jpj) * gcb_tl(1:jpi,1:jpj) * bmask(:,:) )

      epsr = eps * eps * rnorme
      ncut = 0
      ! if rnorme is 0, the solution is 0, the solver is not called
      IF( rnorme == 0.0_wp ) THEN
         gcx_tl(:,:) = 0.0_wp
         res   = 0.0_wp
         niter = 0
         ncut  = 999
      ENDIF

      ! Evaluate the next transport divergence
      ! --------------------------------------
      !    Iterarive solver for the elliptic equation (except IF sol.=0)
      !    (output in gcx with boundary conditions applied)
      kindic = 0
      IF( ncut == 0 ) THEN
         IF    ( nn_solv == 1 ) THEN  ;  CALL ctl_stop('sol_pcg_tan not available in TAM yet')   ! diagonal preconditioned conjuguate gradient
         ELSEIF( nn_solv == 2 ) THEN  ;  CALL sol_sor_tan( kt, kindic )   ! successive-over-relaxation
         ENDIF
      ENDIF

      ! Transport divergence gradient multiplied by z2dt
      ! --------------------------------------------====
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            ! trend of Transport divergence gradient
            ztdgu = z2dtg * ( gcx_tl(ji+1,jj  ) - gcx_tl(ji,jj) ) / e1u(ji,jj)
            ztdgv = z2dtg * ( gcx_tl(ji  ,jj+1) - gcx_tl(ji,jj) ) / e2v(ji,jj)
            ! multiplied by z2dt
            spgu_tl(ji,jj) = z2dt * ztdgu
            spgv_tl(ji,jj) = z2dt * ztdgv
         END DO
      END DO

      ! Add the trends multiplied by z2dt to the after velocity
      ! -------------------------------------------------------
      !     ( c a u t i o n : (ua_tl,va_tl) here are the after velocity not the
      !                       trend, the leap-frog time stepping will not
      !                       be done in dynnxt_tan.F90 routine)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ua_tl(ji,jj,jk) = ( ua_tl(ji,jj,jk) + spgu_tl(ji,jj) ) * umask(ji,jj,jk)
               va_tl(ji,jj,jk) = ( va_tl(ji,jj,jk) + spgv_tl(ji,jj) ) * vmask(ji,jj,jk)
            END DO
         END DO
      END DO

      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_flt_tan')
      !
   END SUBROUTINE dyn_spg_flt_tan

   SUBROUTINE dyn_spg_flt_adj( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  routine dyn_spg_flt_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now trend due to the surface pressure
      !!      gradient in case of filtered free surface formulation  and add
      !!      it to the general trend of momentum equation.
      !!
      !! ** Method of the direct routine:
      !!      Filtered free surface formulation. The surface
      !!      pressure gradient is given by:
      !!         spgu = 1/rau0 d/dx(ps) =  1/e1u di( sshn + btda )
      !!         spgv = 1/rau0 d/dy(ps) =  1/e2v dj( sshn + btda )
      !!      where sshn is the free surface elevation and btda is the after
      !!      time derivative of the free surface elevation
      !!       -1- evaluate the surface presure trend (including the addi-
      !!      tional force) in three steps:
      !!        a- compute the right hand side of the elliptic equation:
      !!            gcb = 1/(e1t e2t) [ di(e2u spgu) + dj(e1v spgv) ]
      !!         where (spgu,spgv) are given by:
      !!            spgu = vertical sum[ e3u (ub+ 2 rdt ua ) ]
      !!                 - grav 2 rdt hu /e1u di[sshn + emp]
      !!            spgv = vertical sum[ e3v (vb+ 2 rdt va) ]
      !!                 - grav 2 rdt hv /e2v dj[sshn + emp]
      !!         and define the first guess from previous computation :
      !!            zbtd = btda
      !!            btda = 2 zbtd - btdb
      !!            btdb = zbtd
      !!        b- compute the relative accuracy to be reached by the
      !!         iterative solver
      !!        c- apply the solver by a call to sol... routine
      !!       -2- compute and add the free surface pressure gradient inclu-
      !!      ding the additional force used to stabilize the equation.
      !!
      !! ** Action : - Update (ua,va) with the surf. pressure gradient trend
      !!
      !! References : Roullet and Madec 1999, JGR.
      !!---------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( IN  ) :: &
         &  kt         ! ocean time-step index
      INTEGER, INTENT( OUT ) :: &
         &  kindic     ! solver convergence flag (<0 if not converge)

      !! * Local declarations
      INTEGER  :: &
         & ji, &     ! dummy loop indices
         & jj, &
         & jk
      REAL(wp) :: &
         & z2dt,   & ! temporary scalars
         & z2dtg,  &
         & zgcb,   &
         & zbtd,   &
         & ztdgu,  &
         & ztdgv
      !!----------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_flt_adj')
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_spg_flt_adj : surface pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   (free surface constant volume case)'
      ENDIF

      ! Local constant initialization
      IF     ( neuler == 0 .AND. kt == nit000     ) THEN

         z2dt = rdt             ! time step: Euler if restart from rest
         CALL sol_mat(kt)       ! initialize matrix

      ELSEIF (                   kt == nitend     ) THEN

         z2dt = 2.0_wp * rdt    ! time step: leap-frog
         CALL sol_mat(kt)       ! reinitialize matrix

      ELSE

         z2dt = 2.0_wp * rdt    ! time step: leap-frog

      ENDIF

      z2dtg  = grav * z2dt

      ! set to zero free surface specific arrays (they are actually local variables)
      spgu_ad(:,:) = 0.0_wp    ;      spgv_ad(:,:) = 0.0_wp

      ! Add the trends multiplied by z2dt to the after velocity
      ! -----------------------------------------------------------
      !     ( c a u t i o n : (ua_ad,va_ad) here are the after velocity not the
      !                       trend, the leap-frog time stepping will not
      !                       be done in dynnxt_adj.F90 routine)
      DO jk = 1, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ua_ad(ji,jj,jk) = ua_ad(ji,jj,jk) * umask(ji,jj,jk)
               va_ad(ji,jj,jk) = va_ad(ji,jj,jk) * vmask(ji,jj,jk)
               spgu_ad(ji,jj)  = spgu_ad(ji,jj)  + ua_ad(ji,jj,jk)
               spgv_ad(ji,jj)  = spgv_ad(ji,jj)  + va_ad(ji,jj,jk)
            END DO
         END DO
      END DO

      ! Transport divergence gradient multiplied by z2dt
      ! --------------------------------------------====
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            ! multiplied by z2dt
            ztdgu = z2dt * spgu_ad(ji,jj)
            ztdgv = z2dt * spgv_ad(ji,jj)
            spgu_ad(ji,jj) = 0.0_wp
            spgv_ad(ji,jj) = 0.0_wp
            ! trend of Transport divergence gradient
            ztdgu = ztdgu * z2dtg / e1u(ji,jj)
            ztdgv = ztdgv * z2dtg / e2v(ji,jj)
            gcx_ad(ji  ,jj  ) = gcx_ad(ji  ,jj  ) - ztdgu - ztdgv
            gcx_ad(ji  ,jj+1) = gcx_ad(ji  ,jj+1) + ztdgv
            gcx_ad(ji+1,jj  ) = gcx_ad(ji+1,jj  ) + ztdgu
         END DO
      END DO

      ! Evaluate the next transport divergence
      ! --------------------------------------
      !    Iterative solver for the elliptic equation (except IF sol.=0)
      !    (output in gcx_ad with boundary conditions applied)

      kindic = 0
      ncut = 0    !  Force solver
      IF( ncut == 0 ) THEN
         IF    ( nn_solv == 1 ) THEN ;  CALL ctl_stop('sol_pcg_adj not available in TMA yet')   ! diagonal preconditioned conjuguate gradient
         ELSEIF( nn_solv == 2 ) THEN ;  CALL sol_sor_adj( kt, kindic )   ! successive-over-relaxation
         ENDIF
      ENDIF

      ! Right hand side of the elliptic equation and first guess
      ! --------------------------------------------------------
      ! apply the lateral boundary conditions
      IF( nn_solv == 2 .AND. MAX( jpr2di, jpr2dj ) > 0 ) &
         &    CALL lbc_lnk_e_adj( gcb_ad, c_solver_pt, 1.0_wp )

      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            ! First guess of the after barotropic transport divergence
            zbtd = gcxb_ad(ji,jj) + 2.0_wp * gcx_ad(ji,jj)
            gcxb_ad(ji,jj) = - gcx_ad(ji,jj)
            gcx_ad (ji,jj) = zbtd
            ! Divergence of the after vertically averaged velocity
            zgcb = gcb_ad(ji,jj) * gcdprc(ji,jj)
            gcb_ad(ji,jj) = 0.0_wp
            spgu_ad(ji-1,jj  ) = spgu_ad(ji-1,jj  ) - zgcb
            spgu_ad(ji  ,jj  ) = spgu_ad(ji  ,jj  ) + zgcb
            spgv_ad(ji  ,jj-1) = spgv_ad(ji  ,jj-1) - zgcb
            spgv_ad(ji  ,jj  ) = spgv_ad(ji  ,jj  ) + zgcb
         END DO
      END DO

      IF( lk_vvl ) CALL ctl_stop( 'dyn_spg_flt_adj: lk_vvl is not available' )

      ! Boundary conditions on (spgu_ad,spgv_ad)
      CALL lbc_lnk_adj( spgu_ad, 'U', -1.0_wp )
      CALL lbc_lnk_adj( spgv_ad, 'V', -1.0_wp )

      ! transport: multiplied by the horizontal scale factor
      DO jj = 2,jpjm1
         DO ji = 2,jpim1   ! vector opt.
            spgu_ad(ji,jj) = spgu_ad(ji,jj) * e2u(ji,jj)
            spgv_ad(ji,jj) = spgv_ad(ji,jj) * e1v(ji,jj)
         END DO
      END DO

      ! compute the next vertically averaged velocity (effect of the additional force not included)
      ! ---------------------------------------------

      ! vertical sum
!CDIR NOLOOPCHG
      IF( lk_vopt_loop ) THEN     ! vector opt., forced unroll
         DO jk = 1, jpkm1
            DO ji = 1, jpij
               ua_ad(ji,1,jk) = ua_ad(ji,1,jk) + e3u(ji,1,jk) * spgu_ad(ji,1)
               va_ad(ji,1,jk) = va_ad(ji,1,jk) + e3v(ji,1,jk) * spgv_ad(ji,1)
            END DO
         END DO
      ELSE                        ! No  vector opt.
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ua_ad(ji,jj,jk) = ua_ad(ji,jj,jk) + e3u(ji,jj,jk) * spgu_ad(ji,jj)
                  va_ad(ji,jj,jk) = va_ad(ji,jj,jk) + e3v(ji,jj,jk) * spgv_ad(ji,jj)
               END DO
            END DO
         END DO
      ENDIF

      DO jj = 2, jpjm1
         DO ji = 2, jpim1 ! vector opt.
            spgu_ad(ji,jj) = 0.0_wp
            spgv_ad(ji,jj) = 0.0_wp
         END DO
      END DO

      IF( nn_cla == 1 )   CALL cla_dynspg_adj( kt )      ! Cross Land Advection (update (ua,va))

      ! Evaluate the masked next velocity (effect of the additional force not included)
      IF( lk_vvl ) THEN          ! variable volume  (surface pressure gradient already included in dyn_hpg)
         !
         IF( ln_dynadv_vec ) THEN      ! vector form : applied on velocity
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ub_ad(ji,jj,jk) = ub_ad(ji,jj,jk) + z2dt * ua_ad(ji,jj,jk) * umask(ji,jj,jk)
                     ua_ad(ji,jj,jk) = z2dt * ua_ad(ji,jj,jk) * umask(ji,jj,jk)
                     vb_ad(ji,jj,jk) = vb_ad(ji,jj,jk) + z2dt * va_ad(ji,jj,jk) * vmask(ji,jj,jk)
                     va_ad(ji,jj,jk) = z2dt * va_ad(ji,jj,jk) * vmask(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
         ELSE                          ! flux form : applied on thickness weighted velocity
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ua_ad(ji,jj,jk) = ua_ad(ji,jj,jk) / e3u(ji,jj,jk) * umask(ji,jj,jk)
                     ub_ad(ji,jj,jk) = ub_ad(ji,jj,jk) + ua_ad(ji,jj,jk) * e3u(ji,jj,jk)
                     ua_ad(ji,jj,jk) = ua_ad(ji,jj,jk) * z2dt * e3u(ji,jj,jk)
                     va_ad(ji,jj,jk) = va_ad(ji,jj,jk) / e3v(ji,jj,jk) * vmask(ji,jj,jk)
                     vb_ad(ji,jj,jk) = vb_ad(ji,jj,jk) + va_ad(ji,jj,jk) * e3v(ji,jj,jk)
                     va_ad(ji,jj,jk) = va_ad(ji,jj,jk) * z2dt * e3v(ji,jj,jk)
                 END DO
               END DO
            END DO
            !
         ENDIF
         !
      ELSE                       ! fixed volume  (add the surface pressure gradient + unweighted time stepping)
         !
         DO jk = 1, jpkm1               ! unweighted time stepping
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ua_ad(  ji,jj,jk) = ua_ad(ji,jj,jk) * umask(ji,jj,jk)
                  ub_ad(  ji,jj,jk) = ub_ad(ji,jj,jk) + ua_ad(ji,jj,jk)
                  ua_ad(  ji,jj,jk) = ua_ad(ji,jj,jk) * z2dt
                  spgu_ad(ji,jj   ) = spgu_ad(ji,jj)  + ua_ad(ji,jj,jk)
                  va_ad(  ji,jj,jk) = va_ad(ji,jj,jk) * vmask(ji,jj,jk)
                  vb_ad(  ji,jj,jk) = vb_ad(ji,jj,jk) + va_ad(ji,jj,jk)
                  va_ad(  ji,jj,jk) = va_ad(ji,jj,jk) * z2dt
                  spgv_ad(ji,jj   ) = spgv_ad(ji,jj)  + va_ad(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jj = 2, jpjm1              ! Surface pressure gradient (now)
            DO ji = 2, jpim1     ! vector opt.
               spgu_ad(ji  ,jj  ) = spgu_ad(ji  ,jj  ) * grav / e1u(ji,jj)
               spgv_ad(ji  ,jj  ) = spgv_ad(ji  ,jj  ) * grav / e2v(ji,jj)
               sshn_ad(ji  ,jj  ) = sshn_ad(ji  ,jj  ) + spgv_ad(ji,jj)
               sshn_ad(ji  ,jj+1) = sshn_ad(ji  ,jj+1) - spgv_ad(ji,jj)
               sshn_ad(ji  ,jj  ) = sshn_ad(ji  ,jj  ) + spgu_ad(ji,jj)
               sshn_ad(ji+1,jj  ) = sshn_ad(ji+1,jj  ) - spgu_ad(ji,jj)
               spgu_ad(ji  ,jj  ) = 0.0_wp
               spgv_ad(ji  ,jj  ) = 0.0_wp
            END DO
         END DO
      ENDIF

      IF( kt == nit000 ) THEN
         ! set to zero free surface specific arrays
         spgu_ad(:,:) = 0.0_wp                    ! surface pressure gradient (i-direction)
         spgv_ad(:,:) = 0.0_wp                    ! surface pressure gradient (j-direction)
         ! Reinitialize the solver arrays
         gcxb_ad(:,:) = 0.e0
         gcx_ad (:,:) = 0.e0
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_flt_adj')
      !
   END SUBROUTINE dyn_spg_flt_adj

   SUBROUTINE dyn_spg_flt_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_spg_flt_adj_tst ***
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
      !! ** Action  :
      !!
      !! History :
      !!        ! 09-01 (A. Weaver)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt        ! Output unit

      !! * Local declarations
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zua_tlin,    & ! Tangent input: ua_tl
         & zva_tlin,    & ! Tangent input: va_tl
         & zub_tlin,    & ! Tangent input: ub_tl
         & zvb_tlin,    & ! Tangent input: vb_tl
         & zua_tlout,   & ! Tangent output: ua_tl
         & zva_tlout,   & ! Tangent output: va_tl
         & zua_adin,    & ! Adjoint input: ua_ad
         & zva_adin,    & ! Adjoint input: va_ad
         & zua_adout,   & ! Adjoint output: ua_ad
         & zva_adout,   & ! Adjoint output: va_ad
         & zub_adout,   & ! Adjoint oputput: ub_ad
         & zvb_adout,   & ! Adjoint output: vb_ad
         & znu            ! 3D random field for u

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zgcx_tlin, zgcxb_tlin, zgcx_tlout, zgcxb_tlout,  &
         & zgcx_adin, zgcxb_adin, zgcx_adout, zgcxb_adout

      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zsshn_tlin, & ! Tangent input: sshn_tl
         & zsshn_adout,& ! Adjoint output: sshn_ad
         & znssh         ! 2D random field for SSH
      REAL(wp) :: &
         & zsp1,    &   ! scalar product involving the tangent routine
         & zsp2         ! scalar product involving the adjoint routine
      INTEGER :: &
         & indic, &
         & istp
      INTEGER :: &
         & ji, &
         & jj, &
         & jk, &
         & kmod, &
         & jstp
      CHARACTER (LEN=14) :: &
         & cl_name
      INTEGER ::             &
         & jpert
      INTEGER, PARAMETER :: jpertmax = 6

      ! Allocate memory

      ALLOCATE( &
         & zua_tlin(jpi,jpj,jpk),  &
         & zva_tlin(jpi,jpj,jpk),  &
         & zub_tlin(jpi,jpj,jpk),  &
         & zvb_tlin(jpi,jpj,jpk),  &
         & zua_tlout(jpi,jpj,jpk), &
         & zva_tlout(jpi,jpj,jpk), &
         & zua_adin(jpi,jpj,jpk),  &
         & zva_adin(jpi,jpj,jpk),  &
         & zua_adout(jpi,jpj,jpk), &
         & zva_adout(jpi,jpj,jpk), &
         & zub_adout(jpi,jpj,jpk), &
         & zvb_adout(jpi,jpj,jpk), &
         & znu(jpi,jpj,jpk)        &
         & )
      ALLOCATE( &
         & zsshn_tlin(jpi,jpj), &
         & zsshn_adout(jpi,jpj),&
         & znssh(jpi,jpj)       &
         & )

      ALLOCATE( zgcx_tlin (jpi,jpj), zgcx_tlout (jpi,jpj), zgcx_adin (jpi,jpj), zgcx_adout (jpi,jpj),  &
                zgcxb_tlin(jpi,jpj), zgcxb_tlout(jpi,jpj), zgcxb_adin(jpi,jpj), zgcxb_adout(jpi,jpj)   )

      !=========================================================================
      !     dx = ( sshb_tl, sshn_tl, ub_tl, ua_tl, vb_tl, va_tl, wn_tl, emp_tl )
      ! and dy = ( sshb_tl, sshn_tl, ua_tl, va_tl )
      !=========================================================================

      ! Test for time steps nit000 and nit000 + 1 (the matrix changes)

      DO jstp = nit000, nitend, nitend-nit000
         DO jpert = jpertmax, jpertmax
            istp = jstp

            !--------------------------------------------------------------------
            ! Reset the tangent and adjoint variables
            !--------------------------------------------------------------------

            zub_tlin (:,:,:) = 0.0_wp
            zvb_tlin (:,:,:) = 0.0_wp
            zua_tlin (:,:,:) = 0.0_wp
            zva_tlin (:,:,:) = 0.0_wp
            zua_tlout(:,:,:) = 0.0_wp
            zva_tlout(:,:,:) = 0.0_wp
            zua_adin (:,:,:) = 0.0_wp
            zva_adin (:,:,:) = 0.0_wp
            zub_adout(:,:,:) = 0.0_wp
            zvb_adout(:,:,:) = 0.0_wp
            zua_adout(:,:,:) = 0.0_wp
            zva_adout(:,:,:) = 0.0_wp

            zsshn_tlin (:,:) = 0.0_wp
            zsshn_adout(:,:) = 0.0_wp

            zgcx_tlout (:,:) = 0.0_wp ; zgcx_adin (:,:) = 0.0_wp ; zgcx_adout (:,:) = 0.0_wp
            zgcxb_tlout(:,:) = 0.0_wp ; zgcxb_adin(:,:) = 0.0_wp ; zgcxb_adout(:,:) = 0.0_wp

            ub_tl(:,:,:) = 0.0_wp
            vb_tl(:,:,:) = 0.0_wp
            ua_tl(:,:,:) = 0.0_wp
            va_tl(:,:,:) = 0.0_wp
            sshn_tl(:,:) = 0.0_wp
            gcx_tl(:,:)  = 0.0_wp
            gcxb_tl(:,:) = 0.0_wp
            spgu_tl(:,:) = 0.0_wp
            spgv_tl(:,:) = 0.0_wp
            ub_ad(:,:,:) = 0.0_wp
            vb_ad(:,:,:) = 0.0_wp
            ua_ad(:,:,:) = 0.0_wp
            va_ad(:,:,:) = 0.0_wp
            sshn_ad(:,:) = 0.0_wp
            gcb_ad(:,:)  = 0.0_wp
            gcx_ad(:,:)  = 0.0_wp
            gcxb_ad(:,:) = 0.0_wp
            spgu_ad(:,:) = 0.0_wp
            spgv_ad(:,:) = 0.0_wp
            !--------------------------------------------------------------------
            ! Initialize the tangent input with random noise: dx
            !--------------------------------------------------------------------
            IF ( (jpert == 1) .OR. (jpert == jpertmax) ) THEN

               CALL grid_random(  znu, 'U', 0.0_wp, stdu )

               DO jk = 1, jpk
                  DO jj = nldj, nlej
                     DO ji = nldi, nlei
                        zua_tlin(ji,jj,jk) = znu(ji,jj,jk)
                     END DO
                  END DO
               END DO

            ENDIF
            IF ( (jpert == 2) .OR. (jpert == jpertmax) ) THEN
               CALL grid_random(  znu, 'V', 0.0_wp, stdv )

               DO jk = 1, jpk
                  DO jj = nldj, nlej
                     DO ji = nldi, nlei
                        zva_tlin(ji,jj,jk) = znu(ji,jj,jk)
                     END DO
                  END DO
               END DO

            ENDIF
            IF ( (jpert == 3) .OR. (jpert == jpertmax) ) THEN
               CALL grid_random(  znu, 'U', 0.0_wp, stdu )

               DO jk = 1, jpk
                  DO jj = nldj, nlej
                     DO ji = nldi, nlei
                        zub_tlin(ji,jj,jk) = znu(ji,jj,jk)
                     END DO
                  END DO
               END DO

            ENDIF
            IF ( (jpert == 4) .OR. (jpert == jpertmax) ) THEN
               CALL grid_random(  znu, 'V', 0.0_wp, stdv )

               DO jk = 1, jpk
                  DO jj = nldj, nlej
                     DO ji = nldi, nlei
                        zvb_tlin(ji,jj,jk) = znu(ji,jj,jk)
                     END DO
                  END DO
               END DO

            ENDIF
            IF ( (jpert == 5) .OR. (jpert == jpertmax) ) THEN

               CALL grid_random(  znssh, 'T', 0.0_wp, stdssh )
               DO jj = nldj, nlej
                  DO ji = nldi, nlei
                     zsshn_tlin(ji,jj) = znssh(ji,jj)
                  END DO
               END DO
            END IF
            zgcx_tlin  (:,:) = ( zua_tlin(:,:,1) + zub_tlin(:,:,1) ) / 10.
            zgcxb_tlin (:,:) = ( zua_tlin(:,:,2) + zub_tlin(:,:,2) ) / 10.
            !--------------------------------------------------------------------
            ! Call the tangent routine: dy = L dx
            !--------------------------------------------------------------------

            ua_tl(:,:,:) = zua_tlin(:,:,:)
            va_tl(:,:,:) = zva_tlin(:,:,:)
            ub_tl(:,:,:) = zub_tlin(:,:,:)
            vb_tl(:,:,:) = zvb_tlin(:,:,:)
            sshn_tl(:,:) = zsshn_tlin(:,:)

            gcb_tl (:,:) = 0.e0
            gcx_tl (:,:) = zgcx_tlin (:,:)   ;   gcxb_tl(:,:) = zgcxb_tlin(:,:)

            CALL sol_mat( istp ) ! for nitend, it is not called in _tan so it is still set to the nit000 case
            CALL dyn_spg_flt_tan( istp, indic )

            zua_tlout(:,:,:) = ua_tl(:,:,:)   ;   zva_tlout(:,:,:) = va_tl(:,:,:)
            zgcxb_tlout(:,:) = gcxb_tl(:,:)   ;   zgcx_tlout (:,:) = gcx_tl (:,:)

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
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zgcx_adin (ji,jj) = zgcx_tlout (ji,jj) &
                     &              * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) * tmask(ji,jj,1)
                  zgcxb_adin(ji,jj) = zgcxb_tlout(ji,jj) &
                     &              * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) * tmask(ji,jj,1)
               END DO
            END DO

            !--------------------------------------------------------------------
            ! Compute the scalar product: ( L dx )^T W dy
            !--------------------------------------------------------------------

            zsp1 =   DOT_PRODUCT( zua_tlout   , zua_adin   ) &
               &   + DOT_PRODUCT( zgcx_tlout  , zgcx_adin  ) &
               &   + DOT_PRODUCT( zgcxb_tlout , zgcxb_adin ) &
               &   + DOT_PRODUCT( zva_tlout   , zva_adin   )


            !--------------------------------------------------------------------
            ! Call the adjoint routine: dx^* = L^T dy^*
            !--------------------------------------------------------------------

            ua_ad(:,:,:) = zua_adin(:,:,:)
            va_ad(:,:,:) = zva_adin(:,:,:)

            gcx_ad (:,:)   = zgcx_adin (:,:)   ;   gcxb_ad(:,:) = zgcxb_adin (:,:)
            ub_ad  (:,:,:) = 0.0_wp ;  vb_ad  (:,:,:) = 0.0_wp

            CALL dyn_spg_flt_adj( istp, indic )

            zua_adout(:,:,:) = ua_ad(:,:,:)
            zva_adout(:,:,:) = va_ad(:,:,:)
            zub_adout(:,:,:) = ub_ad(:,:,:)
            zvb_adout(:,:,:) = vb_ad(:,:,:)
            zsshn_adout(:,:) = sshn_ad(:,:)
            zgcx_adout (:,:) = gcx_ad (:,:)
            zgcxb_adout(:,:) = gcxb_ad(:,:)

            !--------------------------------------------------------------------
            ! Compute the scalar product: dx^T L^T W dy
            !--------------------------------------------------------------------

            zsp2 =   DOT_PRODUCT( zua_tlin  , zua_adout   ) &
               &   + DOT_PRODUCT( zva_tlin  , zva_adout   ) &
               &   + DOT_PRODUCT( zub_tlin  , zub_adout   ) &
               &   + DOT_PRODUCT( zvb_tlin  , zvb_adout   ) &
               &   + DOT_PRODUCT( zgcx_tlin , zgcx_adout  ) &
               &   + DOT_PRODUCT( zgcxb_tlin, zgcxb_adout ) &
               &   + DOT_PRODUCT( zsshn_tlin, zsshn_adout )

            ! Compare the scalar products

            !    14 char:'12345678901234'
            IF ( istp == nit000 ) THEN
               SELECT CASE (jpert)
               CASE(1)
                  cl_name = 'spg_flt  Ua T1'
               CASE(2)
                  cl_name = 'spg_flt  Va T1'
               CASE(3)
                  cl_name = 'spg_flt  Ub T1'
               CASE(4)
                  cl_name = 'spg_flt  Vb T1'
               CASE(5)
                  cl_name = 'spg_flt ssh T1'
               CASE(jpertmax)
                  cl_name = 'dyn_spg_flt T1'
               END SELECT
            ELSEIF ( istp == nitend ) THEN
               SELECT CASE (jpert)
               CASE(1)
                  cl_name = 'spg_flt  Ua T2'
               CASE(2)
                  cl_name = 'spg_flt  Va T2'
               CASE(3)
                  cl_name = 'spg_flt  Ub T2'
               CASE(4)
                  cl_name = 'spg_flt  Vb T2'
               CASE(5)
                  cl_name = 'spg_flt ssh T2'
               CASE(jpertmax)
                  cl_name = 'dyn_spg_flt T2'
               END SELECT
            END IF
            CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

         END DO
      END DO

      nitsor(:) = jp_it0adj  ! restore nitsor to avoid non reproducible results with or without the tests

      ! Deallocate memory

      DEALLOCATE( &
         & zua_tlin,  &
         & zva_tlin,  &
         & zub_tlin,  &
         & zvb_tlin,  &
         & zua_tlout, &
         & zva_tlout, &
         & zua_adin,  &
         & zva_adin,  &
         & zua_adout, &
         & zva_adout, &
         & zub_adout, &
         & zvb_adout, &
         & znu        &
         & )
      DEALLOCATE( &
         & zsshn_tlin, &
         & zsshn_adout,&
         & znssh       &
         & )
      DEALLOCATE( zgcx_tlin , zgcx_tlout , zgcx_adin , zgcx_adout,  &
         & zgcxb_tlin, zgcxb_tlout, zgcxb_adin, zgcxb_adout         )
   END SUBROUTINE dyn_spg_flt_adj_tst

END MODULE dynspg_flt_tam
