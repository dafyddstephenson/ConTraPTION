MODULE dynnxt_tam
   !!======================================================================
   !!                       ***  MODULE  dynnxt_tam  ***
   !! Ocean dynamics: time stepping
   !!                 Tangent and Adjoint Module
   !!======================================================================
   !! History of the direct module:
   !!            OPA  !  1987-02  (P. Andrich, D. L Hostis)  Original code
   !!                 !  1990-10  (C. Levy, G. Madec)
   !!            7.0  !  1993-03  (M. Guyon)  symetrical conditions
   !!            8.0  !  1997-02  (G. Madec & M. Imbard)  opa, release 8.0
   !!            8.2  !  1997-04  (A. Weaver)  Euler forward step
   !!             -   !  1997-06  (G. Madec)  lateral boudary cond., lbc routine
   !!    NEMO    1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-10  (C. Talandier, A-M. Treguier) Open boundary cond.
   !!            2.0  !  2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.3  !  2007-07  (D. Storkey) Calls to BDY routines.
   !!            3.2  !  2009-06  (G. Madec, R.Benshila)  re-introduce the vvl option
   !! History of the TAM routine:
   !!            9.0  !  2008-06  (A. Vidard) Skeleton
   !!                 !  2008-08  (A. Vidard) tangent of the 05-11 version
   !!                 !  2008-08  (A. Vidard) tangent of the 07-07 version
   !!            3.2  !  2010-04  (F. Vigilant) 3.2 conversion
   !!            3.4  !  2012-07  (P.-A. Bouttier) 3.4 conversion
   !!-------------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   dyn_nxt_tan  : update the horizontal velocity from the momentum trend
   !!   dyn_nxt_adj  : update the horizontal velocity from the momentum trend
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE in_out_manager
   USE dynspg_oce
   USE dynadv
   USE lbclnk
   USE lbclnk_tam
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE lib_mpp
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_nxt_tan            ! routine called by step.F90
   PUBLIC dyn_nxt_adj            ! routine called by step.F90
   PUBLIC dyn_nxt_adj_tst        ! routine called by step.F90
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

CONTAINS

   SUBROUTINE dyn_nxt_tan ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_nxt_tan  ***
      !!
      !! ** Purpose :   Compute the after horizontal velocity. Apply the boundary
      !!             condition on the after velocity, achieved the time stepping
      !!             by applying the Asselin filter on now fields and swapping
      !!             the fields.
      !!
      !! ** Method  : * After velocity is compute using a leap-frog scheme:
      !!                       (ua,va) = (ub,vb) + 2 rdt (ua,va)
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the leap-frog is applied on thickness weighted
      !!             velocity.
      !!             Note also that in filtered free surface (lk_dynspg_flt=T),
      !!             the time stepping has already been done in dynspg module
      !!
      !!              * Apply lateral boundary conditions on after velocity
      !!             at the local domain boundaries through lbc_lnk call,
      !!             at the radiative open boundaries (lk_obc=T),
      !!             at the relaxed   open boundaries (lk_bdy=T), and
      !!             at the AGRIF zoom     boundaries (lk_agrif=T)
      !!
      !!              * Apply the time filter applied and swap of the dynamics
      !!             arrays to start the next time step:
      !!                (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
      !!                (un,vn) = (ua,va).
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the time filter is applied on thickness weighted
      !!             velocity.
      !!
      !! ** Action :   ub,vb   filtered before horizontal velocity of next time-step
      !!               un,vn   now horizontal velocity of next time-step
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !! * Local declarations
      INTEGER  ::   ji, jj, jk, iku, ikv         ! dummy loop indices
      REAL(wp) ::   zue3atl , zue3ntl , zue3btl  ! temporary scalar
      REAL(wp) ::   zve3atl , zve3ntl , zve3btl  !    -         -
      REAL(wp) ::   zuftl   , zvftl              !    -         -
      !!----------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_nxt_tan')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_nxt_tan : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      !
      ! Next velocity :   Leap-frog time stepping already done in dynspg_flt.F routine
      ! -------------

      ! Update after velocity on domain lateral boundaries      (only local domain required)
      ! --------------------------------------------------
      CALL lbc_lnk( ua_tl, 'U', -1.0_wp )         ! local domain boundaries
      CALL lbc_lnk( va_tl, 'V', -1.0_wp )
      !
      ! Time filter and swap of dynamics arrays
      ! ------------------------------------------
      IF( neuler == 0 .AND. kt == nit000 ) THEN        !* Euler at first time-step: only swap
         DO jk = 1, jpkm1
            un_tl(:,:,jk) = ua_tl(:,:,jk)                          ! un <-- ua
            vn_tl(:,:,jk) = va_tl(:,:,jk)
         END DO
      ELSE                                             !* Leap-Frog : Asselin filter and swap
         IF( .NOT. lk_vvl ) THEN          ! applied on velocity
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zuftl = un_tl(ji,jj,jk) + atfp * ( ub_tl(ji,jj,jk) - 2._wp * un_tl(ji,jj,jk) + ua_tl(ji,jj,jk) )
                     zvftl = vn_tl(ji,jj,jk) + atfp * ( vb_tl(ji,jj,jk) - 2._wp * vn_tl(ji,jj,jk) + va_tl(ji,jj,jk) )
                     !
                     ub_tl(ji,jj,jk) = zuftl                      ! ub <-- filtered velocity
                     vb_tl(ji,jj,jk) = zvftl
                     un_tl(ji,jj,jk) = ua_tl(ji,jj,jk)             ! un <-- ua
                     vn_tl(ji,jj,jk) = va_tl(ji,jj,jk)
                  END DO
               END DO
            END DO
         ELSE                                                ! applied on thickness weighted velocity
            CALL ctl_stop ( 'dyn_nxt_tan: key_vvl not available yet in TAM' )
         ENDIF
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_nxt_tan')
      !
   END SUBROUTINE dyn_nxt_tan

   SUBROUTINE dyn_nxt_adj ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_nxt_tan  ***
      !!
      !! ** Purpose :   Compute the after horizontal velocity. Apply the boundary
      !!             condition on the after velocity, achieved the time stepping
      !!             by applying the Asselin filter on now fields and swapping
      !!             the fields.
      !!
      !! ** Method  : * After velocity is compute using a leap-frog scheme:
      !!                       (ua,va) = (ub,vb) + 2 rdt (ua,va)
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the leap-frog is applied on thickness weighted
      !!             velocity.
      !!             Note also that in filtered free surface (lk_dynspg_flt=T),
      !!             the time stepping has already been done in dynspg module
      !!
      !!              * Apply lateral boundary conditions on after velocity
      !!             at the local domain boundaries through lbc_lnk call,
      !!             at the radiative open boundaries (lk_obc=T),
      !!             at the relaxed   open boundaries (lk_bdy=T), and
      !!             at the AGRIF zoom     boundaries (lk_agrif=T)
      !!
      !!              * Apply the time filter applied and swap of the dynamics
      !!             arrays to start the next time step:
      !!                (ub,vb) = (un,vn) + atfp [ (ub,vb) + (ua,va) - 2 (un,vn) ]
      !!                (un,vn) = (ua,va).
      !!             Note that with flux form advection and variable volume layer
      !!             (lk_vvl=T), the time filter is applied on thickness weighted
      !!             velocity.
      !!
      !! ** Action :   ub,vb   filtered before horizontal velocity of next time-step
      !!               un,vn   now horizontal velocity of next time-step
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      INTEGER  ::   ji, jj, jk, iku, ikv   ! dummy loop indices
      REAL(wp) ::   zue3aad , zue3nad , zue3bad  ! temporary scalar
      REAL(wp) ::   zve3aad , zve3nad , zve3bad  !    -         -
      REAL(wp) ::   zufad   , zvfad              !    -         -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_nxt_adj')
      !
      ! adjoint local variables initialization
      zue3aad = 0.0_wp ;  zue3nad = 0.0_wp ; zue3bad = 0.0_wp
      zve3aad = 0.0_wp ;  zve3nad = 0.0_wp ; zve3bad = 0.0_wp
      zufad   = 0.0_wp ;  zvfad   = 0.0_wp

      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_nxt_adj : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      ! Time filter and swap of dynamics arrays
      ! ------------------------------------------
      IF( neuler == 0 .AND. kt == nit000 ) THEN        !* Euler at first time-step: only swap
         DO jk = 1, jpkm1
            ua_ad(:,:,jk) = ua_ad(:,:,jk) + un_ad(:,:,jk)          ! un_ad <-- ua_ad
            va_ad(:,:,jk) = va_ad(:,:,jk) + vn_ad(:,:,jk)
            un_ad(:,:,jk) = 0.0_wp
            vn_ad(:,:,jk) = 0.0_wp
         END DO
      ELSE                                             !* Leap-Frog : Asselin filter and swap
         IF( .NOT. lk_vvl ) THEN          ! applied on velocity
            DO jk = 1, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     va_ad(ji,jj,jk) = va_ad(ji,jj,jk) + vn_ad(ji,jj,jk)
                     ua_ad(ji,jj,jk) = ua_ad(ji,jj,jk) + un_ad(ji,jj,jk)
                     un_ad(ji,jj,jk) = 0.0_wp
                     vn_ad(ji,jj,jk) = 0.0_wp
                     zvfad           = zvfad + vb_ad(ji,jj,jk)
                     zufad           = zufad + ub_ad(ji,jj,jk)
                     ub_ad(ji,jj,jk) = 0.0_wp
                     vb_ad(ji,jj,jk) = 0.0_wp

                     ub_ad(ji,jj,jk) = ub_ad(ji,jj,jk) + atfp  * zufad
                     ua_ad(ji,jj,jk) = ua_ad(ji,jj,jk) + atfp  * zufad
                     un_ad(ji,jj,jk) = un_ad(ji,jj,jk) + ( 1 - 2._wp * atfp ) * zufad
                     vb_ad(ji,jj,jk) = vb_ad(ji,jj,jk) + atfp  * zvfad
                     va_ad(ji,jj,jk) = va_ad(ji,jj,jk) + atfp  * zvfad
                     vn_ad(ji,jj,jk) = vn_ad(ji,jj,jk) + ( 1 - 2._wp * atfp ) * zvfad
                     zufad           = 0.0_wp
                     zvfad           = 0.0_wp
                  END DO
               END DO
            END DO
         ELSE                                                ! applied on thickness weighted velocity
            CALL ctl_stop ( 'dyn_nxt_adj: key_vvl not available yet in TAM' )
         ENDIF
      ENDIF

      !
      ! Next velocity :   Leap-frog time stepping already done in dynspg_flt.F routine
      ! -------------

      ! Update after velocity on domain lateral boundaries      (only local domain required)
      ! --------------------------------------------------
      CALL lbc_lnk_adj( ua_ad, 'U', -1.0_wp )         ! local domain boundaries
      CALL lbc_lnk_adj( va_ad, 'V', -1.0_wp )
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_nxt_adj')
      !
   END SUBROUTINE dyn_nxt_adj

   SUBROUTINE dyn_nxt_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_nxt_adj_tst ***
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
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zun_tlin,     & ! Tangent input:   now   u-velocity
         & zvn_tlin,     & ! Tangent input:   now   v-velocity
         & zua_tlin,     & ! Tangent input:  after  u-velocity
         & zva_tlin,     & ! Tangent input:  after  u-velocity
         & zub_tlin,     & ! Tangent input:  before u-velocity
         & zvb_tlin,     & ! Tangent input:  before u-velocity
         & zun_adin,     & ! Adjoint input:   now   u-velocity
         & zvn_adin,     & ! Adjoint input:   now   v-velocity
         & zua_adin,     & ! Adjoint input:  after  u-velocity
         & zva_adin,     & ! Adjoint input:  after  u-velocity
         & zub_adin,     & ! Adjoint input:  before u-velocity
         & zvb_adin,     & ! Adjoint input:  before u-velocity
         & zun_tlout,    & ! Tangent output:  now   u-velocity
         & zvn_tlout,    & ! Tangent output:  now   v-velocity
         & zua_tlout,    & ! Tangent output: after  u-velocity
         & zva_tlout,    & ! Tangent output: after  u-velocity
         & zub_tlout,    & ! Tangent output: before u-velocity
         & zvb_tlout,    & ! Tangent output: before u-velocity
         & zun_adout,    & ! Adjoint output:  now   u-velocity
         & zvn_adout,    & ! Adjoint output:  now   v-velocity
         & zua_adout,    & ! Adjoint output: after  u-velocity
         & zva_adout,    & ! Adjoint output: after  u-velocity
         & zub_adout,    & ! Adjoint output: before u-velocity
         & zvb_adout,    & ! Adjoint output: before u-velocity
         & znu,          & ! 3D random field for u
         & znv,          & ! 3D random field for v
         & zbu,          & ! 3D random field for u
         & zbv,          & ! 3D random field for v
         & zau,          & ! 3D random field for u
         & zav             ! 3D random field for v

      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp1_3,       &
         & zsp1_4,       &
         & zsp1_5,       &
         & zsp1_6,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4,       &
         & zsp2_5,       &
         & zsp2_6
      CHARACTER(LEN=14) :: cl_name

      ! Allocate memory

      ALLOCATE( &
         & zun_tlin(jpi,jpj,jpk),     &
         & zvn_tlin(jpi,jpj,jpk),     &
         & zua_tlin(jpi,jpj,jpk),     &
         & zva_tlin(jpi,jpj,jpk),     &
         & zub_tlin(jpi,jpj,jpk),     &
         & zvb_tlin(jpi,jpj,jpk),     &
         & zun_adin(jpi,jpj,jpk),     &
         & zvn_adin(jpi,jpj,jpk),     &
         & zua_adin(jpi,jpj,jpk),     &
         & zva_adin(jpi,jpj,jpk),     &
         & zub_adin(jpi,jpj,jpk),     &
         & zvb_adin(jpi,jpj,jpk),     &
         & zun_tlout(jpi,jpj,jpk),    &
         & zvn_tlout(jpi,jpj,jpk),    &
         & zua_tlout(jpi,jpj,jpk),    &
         & zva_tlout(jpi,jpj,jpk),    &
         & zub_tlout(jpi,jpj,jpk),    &
         & zvb_tlout(jpi,jpj,jpk),    &
         & zun_adout(jpi,jpj,jpk),    &
         & zvn_adout(jpi,jpj,jpk),    &
         & zua_adout(jpi,jpj,jpk),    &
         & zva_adout(jpi,jpj,jpk),    &
         & zub_adout(jpi,jpj,jpk),    &
         & zvb_adout(jpi,jpj,jpk),    &
         & znu(jpi,jpj,jpk),          &
         & znv(jpi,jpj,jpk),          &
         & zbu(jpi,jpj,jpk),          &
         & zbv(jpi,jpj,jpk),          &
         & zau(jpi,jpj,jpk),          &
         & zav(jpi,jpj,jpk)           &
         & )


      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

          zun_tlin(:,:,:) = 0.0_wp
          zvn_tlin(:,:,:) = 0.0_wp
          zua_tlin(:,:,:) = 0.0_wp
          zva_tlin(:,:,:) = 0.0_wp
          zub_tlin(:,:,:) = 0.0_wp
          zvb_tlin(:,:,:) = 0.0_wp
          zun_adin(:,:,:) = 0.0_wp
          zvn_adin(:,:,:) = 0.0_wp
          zua_adin(:,:,:) = 0.0_wp
          zva_adin(:,:,:) = 0.0_wp
          zub_adin(:,:,:) = 0.0_wp
          zvb_adin(:,:,:) = 0.0_wp
          zun_tlout(:,:,:) = 0.0_wp
          zvn_tlout(:,:,:) = 0.0_wp
          zua_tlout(:,:,:) = 0.0_wp
          zva_tlout(:,:,:) = 0.0_wp
          zub_tlout(:,:,:) = 0.0_wp
          zvb_tlout(:,:,:) = 0.0_wp
          zun_adout(:,:,:) = 0.0_wp
          zvn_adout(:,:,:) = 0.0_wp
          zua_adout(:,:,:) = 0.0_wp
          zva_adout(:,:,:) = 0.0_wp
          zub_adout(:,:,:) = 0.0_wp
          zvb_adout(:,:,:) = 0.0_wp
          znu(:,:,:) = 0.0_wp
          znv(:,:,:) = 0.0_wp
          zbu(:,:,:) = 0.0_wp
          zbv(:,:,:) = 0.0_wp
          zau(:,:,:) = 0.0_wp
          zav(:,:,:) = 0.0_wp

          un_tl(:,:,:) = 0.0_wp
          vn_tl(:,:,:) = 0.0_wp
          ua_tl(:,:,:) = 0.0_wp
          va_tl(:,:,:) = 0.0_wp
          ub_tl(:,:,:) = 0.0_wp
          vb_tl(:,:,:) = 0.0_wp
          un_ad(:,:,:) = 0.0_wp
          vn_ad(:,:,:) = 0.0_wp
          ua_ad(:,:,:) = 0.0_wp
          va_ad(:,:,:) = 0.0_wp
          ub_ad(:,:,:) = 0.0_wp
          vb_ad(:,:,:) = 0.0_wp


      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  znu, 'U', 0.0_wp, stdu )
      CALL grid_random(  znv, 'V', 0.0_wp, stdv )
      CALL grid_random(  zbu, 'U', 0.0_wp, stdu )
      CALL grid_random(  zbv, 'V', 0.0_wp, stdv )
      CALL grid_random(  zau, 'U', 0.0_wp, stdu )
      CALL grid_random(  zav, 'V', 0.0_wp, stdv )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zun_tlin(ji,jj,jk) = znu(ji,jj,jk)
               zvn_tlin(ji,jj,jk) = znv(ji,jj,jk)
               zub_tlin(ji,jj,jk) = zbu(ji,jj,jk)
               zvb_tlin(ji,jj,jk) = zbv(ji,jj,jk)
               zua_tlin(ji,jj,jk) = zau(ji,jj,jk)
               zva_tlin(ji,jj,jk) = zav(ji,jj,jk)
            END DO
         END DO
      END DO

      un_tl(:,:,:) = zun_tlin(:,:,:)
      vn_tl(:,:,:) = zvn_tlin(:,:,:)
      ub_tl(:,:,:) = zub_tlin(:,:,:)
      vb_tl(:,:,:) = zvb_tlin(:,:,:)
      ua_tl(:,:,:) = zua_tlin(:,:,:)
      va_tl(:,:,:) = zva_tlin(:,:,:)

      call dyn_nxt_tan ( nit000 )

      zun_tlout(:,:,:) = un_tl(:,:,:)
      zvn_tlout(:,:,:) = vn_tl(:,:,:)
      zub_tlout(:,:,:) = ub_tl(:,:,:)
      zvb_tlout(:,:,:) = vb_tl(:,:,:)
      zua_tlout(:,:,:) = ua_tl(:,:,:)
      zva_tlout(:,:,:) = va_tl(:,:,:)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zun_adin(ji,jj,jk) = zun_tlout(ji,jj,jk) &
                 &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
                 &               * umask(ji,jj,jk)
              zvn_adin(ji,jj,jk) = zvn_tlout(ji,jj,jk) &
                 &               * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) &
                 &               * vmask(ji,jj,jk)
              zub_adin(ji,jj,jk) = zub_tlout(ji,jj,jk) &
                 &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
                 &               * umask(ji,jj,jk)
              zvb_adin(ji,jj,jk) = zvb_tlout(ji,jj,jk) &
                 &               * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) &
                 &               * vmask(ji,jj,jk)
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

      zsp1_1 = DOT_PRODUCT( zun_tlout, zun_adin )
      zsp1_2 = DOT_PRODUCT( zvn_tlout, zvn_adin )
      zsp1_3 = DOT_PRODUCT( zub_tlout, zub_adin )
      zsp1_4 = DOT_PRODUCT( zvb_tlout, zvb_adin )
      zsp1_5 = DOT_PRODUCT( zua_tlout, zua_adin )
      zsp1_6 = DOT_PRODUCT( zva_tlout, zva_adin )
      zsp1   = zsp1_1 + zsp1_2 + zsp1_3 + zsp1_4 + zsp1_5 + zsp1_6

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      un_ad(:,:,:) = zun_adin(:,:,:)
      vn_ad(:,:,:) = zvn_adin(:,:,:)
      ub_ad(:,:,:) = zub_adin(:,:,:)
      vb_ad(:,:,:) = zvb_adin(:,:,:)
      ua_ad(:,:,:) = zua_adin(:,:,:)
      va_ad(:,:,:) = zva_adin(:,:,:)

      CALL dyn_nxt_adj ( nit000 )

      zun_adout(:,:,:) = un_ad(:,:,:)
      zvn_adout(:,:,:) = vn_ad(:,:,:)
      zub_adout(:,:,:) = ub_ad(:,:,:)
      zvb_adout(:,:,:) = vb_ad(:,:,:)
      zua_adout(:,:,:) = ua_ad(:,:,:)
      zva_adout(:,:,:) = va_ad(:,:,:)

      zsp2_1 = DOT_PRODUCT( zun_tlin, zun_adout )
      zsp2_2 = DOT_PRODUCT( zvn_tlin, zvn_adout )
      zsp2_3 = DOT_PRODUCT( zub_tlin, zub_adout )
      zsp2_4 = DOT_PRODUCT( zvb_tlin, zvb_adout )
      zsp2_5 = DOT_PRODUCT( zua_tlin, zua_adout )
      zsp2_6 = DOT_PRODUCT( zva_tlin, zva_adout )
      zsp2   = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5 + zsp2_6

      ! Compare the scalar products
      ! 14 char:'12345678901234'
      cl_name = 'dyn_nxt_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE( &
         & zun_tlin,     &
         & zvn_tlin,     &
         & zua_tlin,     &
         & zva_tlin,     &
         & zub_tlin,     &
         & zvb_tlin,     &
         & zun_adin,     &
         & zvn_adin,     &
         & zua_adin,     &
         & zva_adin,     &
         & zub_adin,     &
         & zvb_adin,     &
         & zun_tlout,    &
         & zvn_tlout,    &
         & zua_tlout,    &
         & zva_tlout,    &
         & zub_tlout,    &
         & zvb_tlout,    &
         & zun_adout,    &
         & zvn_adout,    &
         & zua_adout,    &
         & zva_adout,    &
         & zub_adout,    &
         & zvb_adout,    &
         & znu,          &
         & znv,          &
         & zbu,          &
         & zbv,          &
         & zau,          &
         & zav           &
         & )


   END SUBROUTINE dyn_nxt_adj_tst
   !!======================================================================
END MODULE dynnxt_tam
