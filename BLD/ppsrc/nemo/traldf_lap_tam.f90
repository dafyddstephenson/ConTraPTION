MODULE traldf_lap_tam
   !!==============================================================================
   !!                       ***  MODULE  traldf_lap  ***
   !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
   !!                        Tangent and adjoint module
   !!==============================================================================
   !! History :   9.0  !  ?????
   !! History of the T&A module:
   !!             9.0  !  09-03  (F. Vigilant) original version
   !!       NEMO  3.4  ! 12-07   (P.-A. Bouttier)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf_lap  : update the tracer trend with the horizontal diffusion
   !!                 using a iso-level harmonic (laplacien) operator.
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE ldftra_oce
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE paresp
   USE trc_oce
   USE lib_mpp
   USE timing
   USE wrk_nemo

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_ldf_lap_tan      ! routine called by tradldf_tam.F90
   PUBLIC tra_ldf_lap_adj      ! routine called by tradldf_tam.F90
   PUBLIC tra_ldf_lap_adj_tst  ! routine called by tradldf_tam.F90

   REAL(wp), SAVE, ALLOCATABLE, DIMENSION(:,:) ::   e1ur, e2vr   ! scale factor coefficients

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
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c2d' :                 aht: 2D coefficient
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
   !! $Id: traldf_lap.F90 1152 2008-06-26 14:11:13Z rblod $
   !! This software is governed by the CeCILL licence see modipsl/doc/NEMO_CeCILL.txt
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_ldf_lap_tan( kt, kit000, cdtype, pgu_tl, pgv_tl,      &
      &                                ptb_tl, pta_tl, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_lap_tan  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   Second order diffusive operator evaluated using before
      !!      fields (forward time scheme). The horizontal diffusive trends of
      !!      temperature (idem for salinity) is given by:
      !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(tb) ]
      !!                                   + dj-1[ aht e1v*e3v/e2v dj(tb) ] }
      !!     Note: key_zco defined, the e3t=e3u=e3v, the trend becomes:
      !!          difft = 1/(e1t*e2t) {  di-1[ aht e2u/e1u di(tb) ]
      !!                               + dj-1[ aht e1v/e2v dj(tb) ] }
      !!      Add this trend to the general tracer trend (ta,sa):
      !!          (ta,sa) = (ta,sa) + ( difft , diffs )
      !!
      !! ** Action  : - Update (ta,sa) arrays with the before iso-level
      !!                harmonic mixing trend.
      !!
      !! History :
      !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!        !  05-11  (G. Madec)  add zps case
      !! History of the tangent routine
      !!   9.0  !  03-09 (F. Vigilant) tangent of 9.0
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000           ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype           ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt             ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu_tl, pgv_tl   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb_tl           ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_tl           ! tracer trend
      !! * Local declarations
      INTEGER  ::   ji, jj, jk, jn, ierr         ! dummy loop indices
      INTEGER  ::   iku, ikv               ! temporary integers
      REAL(wp) ::  zabe1, zabe2, zbtr      ! temporary scalars
      REAL(wp) ::  ztatl, zsatl            ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztutl, ztvtl     ! 3D workspace
      !!----------------------------------------------------------------------
      !
      CALL timing_start('tra_ldf_lap_tan')
      !
      CALL wrk_alloc(jpi,jpj,jpk, ztutl, ztvtl)
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_lap_tan : iso-level laplacian diffusion on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         IF( .NOT. ALLOCATED( e1ur ) ) THEN
            ! This routine may be called for both active and passive tracers.
            ! Allocate and set saved arrays on first call only.
            ALLOCATE( e1ur(jpi,jpj), e2vr(jpi,jpj), STAT=ierr )
            IF( lk_mpp    )   CALL mpp_sum( ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'tra_ldf_lap : unable to allocate arrays' )
            !
            e1ur(:,:) = e2u(:,:) / e1u(:,:)
            e2vr(:,:) = e1v(:,:) / e2v(:,:)
         ENDIF
      ENDIF

      !
      DO jn = 1, kjpt
         !                                                  ! =============
         DO jk = 1, jpkm1                                   ! Vertical slab
            !                                               ! =============
            ! 1. First derivative (gradient)
            ! -------------------
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zabe1 = rldf * ahtu(ji,jj) * umask(ji,jj,jk) * e1ur(ji,jj) * e3u(ji,jj,jk)
                  zabe2 = rldf * ahtv(ji,jj) * vmask(ji,jj,jk) * e2vr(ji,jj) * e3v(ji,jj,jk)
                  ztutl(ji,jj,jk) = zabe1 * ( ptb_tl(ji+1,jj  ,jk,jn) - ptb_tl(ji,jj,jk,jn) )
                  ztvtl(ji,jj,jk) = zabe2 * ( ptb_tl(ji  ,jj+1,jk,jn) - ptb_tl(ji,jj,jk,jn) )
               END DO
            END DO
            IF( ln_zps ) THEN      ! set gradient at partial step level
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1   ! vector opt.
                     ! last level
                     iku = mbku(ji,jj)
                     ikv = mbku(ji,jj)
                     IF( iku == jk ) THEN
                        zabe1 = rldf * ahtu(ji,jj) * umask(ji,jj,iku) * e1ur(ji,jj) * e3u(ji,jj,iku)
                        ztutl(ji,jj,jk) = zabe1 * pgu_tl(ji,jj,jn)
                     ENDIF
                     IF( ikv == jk ) THEN
                        zabe2 = rldf * ahtv(ji,jj) * vmask(ji,jj,ikv) * e2vr(ji,jj) * e3v(ji,jj,ikv)
                        ztvtl(ji,jj,jk) = zabe2 * pgv_tl(ji,jj,jn)
                     ENDIF
                  END DO
               END DO
            ENDIF


            ! 2. Second derivative (divergence)
            ! --------------------
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zbtr = 1._wp / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk))
                  ! horizontal diffusive trends
                  ztatl = zbtr * (  ztutl(ji,jj,jk) - ztutl(ji-1,jj  ,jk)   &
                     &            + ztvtl(ji,jj,jk) - ztvtl(ji  ,jj-1,jk)  )
                  ! add it to the general tracer trends
                  pta_tl(ji,jj,jk,jn) = pta_tl(ji,jj,jk,jn) + ztatl
               END DO
            END DO
            !                                               ! =============
         END DO                                             !  End of slab
         !                                                  ! =============
      END DO
      !
      CALL timing_stop('tra_ldf_lap_tan')
      !
      CALL wrk_dealloc(jpi,jpj,jpk, ztutl, ztvtl)
      !

   END SUBROUTINE tra_ldf_lap_tan


   SUBROUTINE tra_ldf_lap_adj( kt, kit000, cdtype, pgu_ad, pgv_ad,      &
      &                                ptb_ad, pta_ad, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_lap_adj  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   Second order diffusive operator evaluated using before
      !!      fields (forward time scheme). The horizontal diffusive trends of
      !!      temperature (idem for salinity) is given by:
      !!          difft = 1/(e1t*e2t*e3t) {  di-1[ aht e2u*e3u/e1u di(tb) ]
      !!                                   + dj-1[ aht e1v*e3v/e2v dj(tb) ] }
      !!     Note: key_zco defined, the e3t=e3u=e3v, the trend becomes:
      !!          difft = 1/(e1t*e2t) {  di-1[ aht e2u/e1u di(tb) ]
      !!                               + dj-1[ aht e1v/e2v dj(tb) ] }
      !!      Add this trend to the general tracer trend (ta,sa):
      !!          (ta,sa) = (ta,sa) + ( difft , diffs )
      !!
      !! ** Action  : - Update (ta,sa) arrays with the before iso-level
      !!                harmonic mixing trend.
      !!
      !! History :
      !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !!        !  05-11  (G. Madec)  add zps case
      !! History of the tangent routine
      !!   9.0  !  03-09 (F. Vigilant) adjoint of 9.0
      !!----------------------------------------------------------------------

      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
            INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(inout) ::   pgu_ad, pgv_ad   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   ptb_ad        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_ad        ! tracer trend
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn             ! dummy loop indices
      INTEGER ::   iku, ikv, ierr               ! temporary integers
      REAL(wp) ::  zabe1, zabe2, zbtr     ! temporary scalars
      REAL(wp) ::  ztaad, zsaad           ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) :: ztuad, ztvad    ! 3D workspace
      !!----------------------------------------------------------------------
      !
      CALL timing_start('tra_ldf_lap_adj')
      !
      CALL wrk_alloc(jpi,jpj,jpk, ztuad, ztvad)
      !
      ztvad(:,:,:) = 0.0_wp
      ztuad(:,:,:) = 0.0_wp
      ztaad        = 0.0_wp
      zsaad        = 0.0_wp

      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_lap_adj : iso-level laplacian diffusion on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         IF( .NOT. ALLOCATED( e1ur ) ) THEN
            ! This routine may be called for both active and passive tracers.
            ! Allocate and set saved arrays on first call only.
            ALLOCATE( e1ur(jpi,jpj), e2vr(jpi,jpj), STAT=ierr )
            IF( lk_mpp    )   CALL mpp_sum( ierr )
            IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'tra_ldf_lap : unable to allocate arrays' )
            !
            e1ur(:,:) = e2u(:,:) / e1u(:,:)
            e2vr(:,:) = e1v(:,:) / e2v(:,:)
         ENDIF
      ENDIF


      DO jn = 1, kjpt
         !                                                  ! =============
         DO jk = 1, jpkm1                                   ! Vertical slab
            !                                               ! =============

            ! 2. Second derivative (divergence)
            ! --------------------
            DO jj = jpjm1, 2, -1
               DO ji = jpim1, 2, -1   ! vector opt.
                  zbtr = 1._wp / ( e1t(ji,jj) *e2t(ji,jj) * e3t(ji,jj,jk) )
                  ! horizontal diffusive trends
                  ztaad             = zbtr * pta_ad(ji,jj,jk,jn)
                  ztuad(ji  ,jj  ,jk) = ztuad(ji  ,jj  ,jk) + ztaad
                  ztuad(ji-1,jj  ,jk) = ztuad(ji-1,jj  ,jk) - ztaad
                  ztvad(ji  ,jj  ,jk) = ztvad(ji  ,jj  ,jk) + ztaad
                  ztvad(ji  ,jj-1,jk) = ztvad(ji  ,jj-1,jk) - ztaad
                  !ztaad   = 0.0_wp
                  !zsaad   = 0.0_wp
               END DO
            END DO

            ! 1. First derivative (gradient)
            ! -------------------

            IF( ln_zps ) THEN      ! set gradient at partial step level
               DO jj = jpjm1, 1, -1
                  DO ji = jpim1, 1, -1   ! vector opt.?
                     ! last level
                     iku = mbku(ji,jj)
                     ikv = mbkv(ji,jj)
                     IF( iku == jk ) THEN
                        zabe1 = rldf * ahtu(ji,jj) * umask(ji,jj,iku) * e1ur(ji,jj) * e3u(ji,jj,iku)
                        pgu_ad(ji,jj,jn) = zabe1 * ztuad(ji,jj,jk) + pgu_ad(ji,jj,jn)
                        ztuad(ji,jj,jk)  = 0.0_wp
                     ENDIF
                     IF( ikv == jk ) THEN
                        zabe2 = rldf * ahtv(ji,jj) * vmask(ji,jj,ikv) * e2vr(ji,jj) * e3v(ji,jj,ikv)
                        pgv_ad(ji,jj,jn) = zabe2 * ztvad(ji,jj,jk) + pgv_ad(ji,jj,jn)
                        ztvad(ji,jj,jk)  = 0.0_wp
                     ENDIF
                  END DO
               END DO
            ENDIF
            DO jj = jpjm1, 1, -1
               DO ji = jpim1, 1, -1   ! vector opt. ?
                  zabe1 = rldf * ahtu(ji,jj) * umask(ji,jj,jk) * e1ur(ji,jj) * e3u(ji,jj,jk)
                  zabe2 = rldf * ahtv(ji,jj) * vmask(ji,jj,jk) * e2vr(ji,jj) * e3v(ji,jj,jk)
                  ptb_ad(ji  ,jj  ,jk,jn) = ptb_ad(ji  ,jj  ,jk,jn) - (zabe1 * ztuad(ji,jj,jk) + zabe2 * ztvad(ji,jj,jk))
                  ptb_ad(ji+1,jj  ,jk,jn) = ptb_ad(ji+1,jj  ,jk,jn) +  zabe1 * ztuad(ji,jj,jk)
                  ptb_ad(ji  ,jj+1,jk,jn) = ptb_ad(ji  ,jj+1,jk,jn) +  zabe2 * ztvad(ji,jj,jk)
                  ztuad(ji,jj,jk)   = 0.0_wp
                  ztvad(ji,jj,jk)   = 0.0_wp
               END DO
            END DO
             !                                               ! =============
         END DO                                              !  End of slab
      !                                                      ! =============
      END DO
      !
      CALL timing_stop('tra_ldf_lap_adj')
      !
      CALL wrk_dealloc(jpi,jpj,jpk, ztuad, ztvad)
      !
   END SUBROUTINE tra_ldf_lap_adj


   SUBROUTINE tra_ldf_lap_adj_tst ( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE example_adj_tst ***
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

      !! * Local declarations
      INTEGER ::  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_T,       &
         & zsp1_S,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       &
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4,       &
         & zsp2_5,       &
         & zsp2_6,       &
         & zsp2_7,       &
         & zsp2_8,       &
         & zsp2_T,       &
         & zsp2_S
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & ztb_tlin ,      & ! Tangent input
         & zsb_tlin ,      & ! Tangent input
         & zta_tlin ,      & ! Tangent input
         & zsa_tlin ,      & ! Tangent input
         & zta_tlout,      & ! Tangent output
         & zsa_tlout,      & ! Tangent output
         & zta_adin,       & ! Adjoint input
         & zsa_adin,       & ! Adjoint input
         & ztb_adout ,     & ! Adjoint output
         & zsb_adout ,     & ! Adjoint output
         & zta_adout ,     & ! Adjoint output
         & zsa_adout ,     & ! Adjoint output
         & z3r               ! 3D random field
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zgtu_tlin ,     & ! Tangent input
         & zgsu_tlin ,     & ! Tangent input
         & zgtv_tlin ,     & ! Tangent input
         & zgsv_tlin ,     & ! Tangent input
         & zgtu_adout ,    & ! Adjoint output
         & zgsu_adout ,    & ! Adjoint output
         & zgtv_adout ,    & ! Adjoint output
         & zgsv_adout ,    & ! Adjoint output
         & z2r               ! 2D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & ztb_tlin(jpi,jpj,jpk),      &
         & zsb_tlin(jpi,jpj,jpk),      &
         & zta_tlin(jpi,jpj,jpk),      &
         & zsa_tlin(jpi,jpj,jpk),      &
         & zgtu_tlin(jpi,jpj),         &
         & zgsu_tlin(jpi,jpj),         &
         & zgtv_tlin(jpi,jpj),         &
         & zgsv_tlin(jpi,jpj),         &
         & zta_tlout(jpi,jpj,jpk),     &
         & zsa_tlout(jpi,jpj,jpk),     &
         & zta_adin(jpi,jpj,jpk),      &
         & zsa_adin(jpi,jpj,jpk),      &
         & ztb_adout(jpi,jpj,jpk),     &
         & zsb_adout(jpi,jpj,jpk),     &
         & zta_adout(jpi,jpj,jpk),     &
         & zsa_adout(jpi,jpj,jpk),     &
         & zgtu_adout(jpi,jpj),        &
         & zgsu_adout(jpi,jpj),        &
         & zgtv_adout(jpi,jpj),        &
         & zgsv_adout(jpi,jpj),        &
         & z3r(jpi,jpj,jpk),           &
         & z2r(jpi,jpj)                &
         & )


      !=======================================================================
      ! 1) dx = ( tb_tl, ta_tl, sb_tl, sa_tl, gtu_tl, gtv_tl, gsu_tl, gsv_tl )
      !    dy = ( ta_tl, sa_tl )
      !=======================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      ztb_tlin(:,:,:)  = 0.0_wp
      zsb_tlin(:,:,:)  = 0.0_wp
      zta_tlin(:,:,:)  = 0.0_wp
      zsa_tlin(:,:,:)  = 0.0_wp
      zgtu_tlin(:,:)   = 0.0_wp
      zgsu_tlin(:,:)   = 0.0_wp
      zgtv_tlin(:,:)   = 0.0_wp
      zgsv_tlin(:,:)   = 0.0_wp
      zta_tlout(:,:,:) = 0.0_wp
      zsa_tlout(:,:,:) = 0.0_wp
      zta_adin(:,:,:)  = 0.0_wp
      zsa_adin(:,:,:)  = 0.0_wp
      ztb_adout(:,:,:) = 0.0_wp
      zsb_adout(:,:,:) = 0.0_wp
      zta_adout(:,:,:) = 0.0_wp
      zsa_adout(:,:,:) = 0.0_wp
      zgtu_adout(:,:)  = 0.0_wp
      zgsu_adout(:,:)  = 0.0_wp
      zgtv_adout(:,:)  = 0.0_wp
      zgsv_adout(:,:)  = 0.0_wp

      tsb_tl(:,:,:,:) = 0.0_wp
      tsa_tl(:,:,:,:) = 0.0_wp
      gtsu_tl(:,:,:)  = 0.0_wp
      gtsv_tl(:,:,:)  = 0.0_wp
      tsb_ad(:,:,:,:) = 0.0_wp
      tsa_ad(:,:,:,:) = 0.0_wp
      gtsu_ad(:,:,:)  = 0.0_wp
      gtsv_ad(:,:,:)  = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------
      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      ztb_tlin(:,:,:) = z3r(:,:,:)
      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      zsb_tlin(:,:,:) = z3r(:,:,:)
      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      zta_tlin(:,:,:) = z3r(:,:,:)
      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      zsa_tlin(:,:,:) = z3r(:,:,:)
      CALL grid_random(  z2r, 'U', 0.0_wp, stds )
      zgtu_tlin(:,:) = z2r(:,:)
      CALL grid_random(  z2r, 'U', 0.0_wp, stds )
      zgsu_tlin(:,:) = z2r(:,:)
      CALL grid_random(  z2r, 'V', 0.0_wp, stds )
      zgtv_tlin(:,:) = z2r(:,:)
      CALL grid_random(  z2r, 'V', 0.0_wp, stds )
      zgsv_tlin(:,:) = z2r(:,:)

      tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
      tsb_tl(:,:,:,jp_sal) = zsb_tlin(:,:,:)
      tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
      tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)
      gtsu_tl(:,:,jp_tem)  = zgtu_tlin(:,:)
      gtsu_tl(:,:,jp_sal)  = zgsu_tlin(:,:)
      gtsv_tl(:,:,jp_tem)  = zgtv_tlin(:,:)
      gtsv_tl(:,:,jp_sal)  = zgsv_tlin(:,:)

      CALL tra_ldf_lap_tan( nit000, nit000, 'TRA', gtsu_tl, gtsv_tl, tsb_tl, tsa_tl, jpts )

      zta_tlout(:,:,:) = tsa_tl(:,:,:,jp_tem)
      zsa_tlout(:,:,:) = tsa_tl(:,:,:,jp_sal)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsa_adin(ji,jj,jk) = zsa_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_s(jk)
              zta_adin(ji,jj,jk) = zta_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_t(jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !-------------------------------------------------------------------

      zsp1_T = DOT_PRODUCT( zta_tlout, zta_adin )
      zsp1_S = DOT_PRODUCT( zsa_tlout, zsa_adin )
      zsp1 = zsp1_T + zsp1_S

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      tsa_ad(:,:,:,jp_tem) = zta_adin(:,:,:)
      tsa_ad(:,:,:,jp_sal) = zsa_adin(:,:,:)

      CALL tra_ldf_lap_adj( nit000, nit000, 'TRA', gtsu_ad, gtsv_ad, tsb_ad, tsa_ad, jpts )

      ztb_adout(:,:,:) = tsb_ad(:,:,:,jp_tem)
      zsb_adout(:,:,:) = tsb_ad(:,:,:,jp_sal)
      zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)
      zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)
      zgtu_adout(:,:)  = gtsu_ad(:,:,jp_tem)
      zgsu_adout(:,:)  = gtsu_ad(:,:,jp_sal)
      zgtv_adout(:,:)  = gtsv_ad(:,:,jp_tem)
      zgsv_adout(:,:)  = gtsv_ad(:,:,jp_sal)

      zsp2_1 = DOT_PRODUCT( ztb_tlin , ztb_adout  )
      zsp2_2 = DOT_PRODUCT( zta_tlin , zta_adout  )
      zsp2_3 = DOT_PRODUCT( zgtu_tlin, zgtu_adout )
      zsp2_4 = DOT_PRODUCT( zgtv_tlin, zgtv_adout )
      zsp2_5 = DOT_PRODUCT( zsb_tlin , zsb_adout  )
      zsp2_6 = DOT_PRODUCT( zsa_tlin , zsa_adout  )
      zsp2_7 = DOT_PRODUCT( zgsu_tlin, zgsu_adout )
      zsp2_8 = DOT_PRODUCT( zgsv_tlin, zgsv_adout )

      zsp2_T = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4
      zsp2_S = zsp2_5 + zsp2_6 + zsp2_7 + zsp2_8
      zsp2   = zsp2_T + zsp2_S

      cl_name = 'tra_ldf_lap_ad'
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(         &
         & ztb_tlin,      & ! Tangent input
         & zsb_tlin,      & ! Tangent input
         & zta_tlin,      & ! Tangent input
         & zsa_tlin,      & ! Tangent input
         & zgtu_tlin,     & ! Tangent input
         & zgsu_tlin,     & ! Tangent input
         & zgtv_tlin,     & ! Tangent input
         & zgsv_tlin,     & ! Tangent input
         & zta_tlout,     & ! Tangent output
         & zsa_tlout,     & ! Tangent output
         & zta_adin,      & ! Adjoint input
         & zsa_adin,      & ! Adjoint input
         & ztb_adout,     & ! Adjoint output
         & zsb_adout,     & ! Adjoint output
         & zta_adout,     & ! Adjoint output
         & zsa_adout,     & ! Adjoint output
         & zgtu_adout,    & ! Adjoint output
         & zgsu_adout,    & ! Adjoint output
         & zgtv_adout,    & ! Adjoint output
         & zgsv_adout,    & ! Adjoint output
         & z3r,           & ! 3D random field
         & z2r            &
         & )


   END SUBROUTINE tra_ldf_lap_adj_tst
   !!==============================================================================
END MODULE traldf_lap_tam
