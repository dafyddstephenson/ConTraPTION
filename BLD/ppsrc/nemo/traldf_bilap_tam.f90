MODULE traldf_bilap_tam
   !!===========================================================================
   !!                       ***  MODULE  dynldf_bilap_tam  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!                  tangent and Adjoint Module
   !!===========================================================================

   !!---------------------------------------------------------------------------
   !!   dyn_ldf_bilap_tan : update the momentum trend with the lateral diffusion
   !!                       using an iso-level bilaplacian operator (tangent)
   !!   dyn_ldf_bilap_adj : update the momentum trend with the lateral diffusion
   !!                       using an iso-level bilaplacian operator (adjoint)
   !!---------------------------------------------------------------------------

	!! * Modules used
   USE lbclnk
   USE lbclnk_tam
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
   USE wrk_nemo
   USE timing

IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_ldf_bilap_tan  ! called by dynldf_tam.F90
   PUBLIC tra_ldf_bilap_adj  ! called by dynldf_tam.F90
   PUBLIC tra_ldf_bilap_adj_tst  ! routine called by tradldf_tam.F90
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
   !!                   ***  ldfeiv_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaei. the eddy induced velocity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfeiv_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'traldf_c2d' :                           eiv: 2D coefficient
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
   SUBROUTINE tra_ldf_bilap_tan( kt, kit000, cdtype, pgu_tl, pgv_tl,      &
      &                                  ptb_tl, pta_tl, kjpt  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   4th order diffusive operator along model level surfaces
      !!      evaluated using before fields (forward time scheme). The hor.
      !!      diffusive trends of temperature (idem for salinity) is given by:
      !!      Laplacian of tb:
      !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(tb) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(tb) ]  }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
      !!      Note: if key_zco defined, e3t=e3u=e3v, they are simplified.
      !!
      !!      Add this trend to the general trend (ta,sa):
      !!         (ta,sa) = (ta,sa) + ( difft , diffs )
      !!
      !! ** Action : - Update (ta,sa) arrays with the before iso-level
      !!               biharmonic mixing trend.
      !!
      !! History :
      !!        !  91-11  (G. Madec)  Original code
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  96-01  (M. Imbard)  mpp exchange
      !!        !  97-07  (G. Madec)  optimization, and ahtt
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. talandier) New trends organization
      !!        !  05-11  (G. Madec)  zps or sco as default option
      !!----------------------------------------------------------------------
      !! History of the tangent routine
      !!   9.0  !  10-05 (P.A. Bouttier) tangent of 9.0
      !!----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000           ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype           ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt             ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(in   ) ::   pgu_tl, pgv_tl   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb_tl           ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_tl           ! tracer trend
      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn      ! dummy loop indices
      INTEGER ::   iku, ikv               ! temporary integers
      REAL(wp) ::   ztatl, zsatl, zbtr              ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   &
         & zeeu, zeev,               & ! 2D workspace
         & zlttl
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   &
         & ztutl, ztvtl       ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_ldf_bilap_tan')
      !
      CALL wrk_alloc( jpi, jpj, zeeu, zeev, zlttl )
      CALL wrk_alloc( jpi, jpj, jpk, ztutl, ztvtl )
      !
      IF( kt == kit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_bilap_tan : iso-level biharmonic operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF

      DO jn = 1, kjpt
      !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============

            ! 0. Initialization of metric arrays (for z- or s-coordinates)
            ! ----------------------------------

            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zeeu(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk) / e1u(ji,jj) * umask(ji,jj,jk)
                  zeev(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk) / e2v(ji,jj) * vmask(ji,jj,jk)
               END DO
            END DO

            ! 1. Laplacian
            ! ------------

            ! First derivative (gradient)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  ztutl(ji,jj,jk) = zeeu(ji,jj) * ( ptb_tl(ji+1,jj  ,jk,jn) - ptb_tl(ji,jj,jk,jn) )
                  ztvtl(ji,jj,jk) = zeev(ji,jj) * ( ptb_tl(ji  ,jj+1,jk,jn) - ptb_tl(ji,jj,jk,jn) )
               END DO
            END DO
            IF( ln_zps ) THEN      ! set gradient at partial step level
               DO jj = 1, jpjm1
                  DO ji = 1, jpim1
                     IF( mbku(ji,jj) == jk )  ztutl(ji,jj,jk) = zeeu(ji,jj) * pgu_tl(ji,jj,jn)
                     IF( mbkv(ji,jj) == jk )  ztvtl(ji,jj,jk) = zeev(ji,jj) * pgv_tl(ji,jj,jn)
                  END DO
               END DO
            ENDIF

            ! Second derivative (divergence)
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zbtr = 1.0 / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                  zlttl(ji,jj) = rldf * ahtt(ji,jj) * zbtr * (  ztutl(ji,jj,jk) - ztutl(ji-1,jj,jk)   &
                               &                            + ztvtl(ji,jj,jk) - ztvtl(ji,jj-1,jk)  )
               END DO
            END DO

            ! Lateral boundary conditions on the laplacian (zlttl,zlstl)   (unchanged sgn)
            CALL lbc_lnk( zlttl, 'T', 1.0_wp )

            ! 2. Bilaplacian
            ! --------------

            ! third derivative (gradient)
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  ztutl(ji,jj,jk) = zeeu(ji,jj) * ( zlttl(ji+1,jj  ) - zlttl(ji,jj) )
                  ztvtl(ji,jj,jk) = zeev(ji,jj) * ( zlttl(ji  ,jj+1) - zlttl(ji,jj) )
               END DO
            END DO

            ! fourth derivative (divergence) and add to the general tracer trend
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ! horizontal diffusive trends
                  zbtr = 1.0 / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                  ztatl = zbtr * (  ztutl(ji,jj,jk) - ztutl(ji-1,jj,jk) + ztvtl(ji,jj,jk) - ztvtl(ji,jj-1,jk)  )
                  ! add it to the general tracer trends
                  pta_tl(ji,jj,jk,jn) = pta_tl(ji,jj,jk,jn) + ztatl
               END DO
            END DO
            !                                             ! ===============
         END DO                                           ! Horizontal slab
         !                                                ! ===============
      END DO
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_ldf_bilap_tan')
      !
      CALL wrk_dealloc( jpi, jpj, zeeu, zeev, zlttl )
      CALL wrk_dealloc( jpi, jpj, jpk, ztutl, ztvtl )
      !
   END SUBROUTINE tra_ldf_bilap_tan


   SUBROUTINE tra_ldf_bilap_adj( kt, kit000, cdtype, pgu_ad, pgv_ad,      &
      &                                  ptb_ad, pta_ad, kjpt   )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_bilap  ***
      !!
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   4th order diffusive operator along model level surfaces
      !!      evaluated using before fields (forward time scheme). The hor.
      !!      diffusive trends of temperature (idem for salinity) is given by:
      !!      Laplacian of tb:
      !!         zlt   = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(tb) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(tb) ]  }
      !!      Multiply by the eddy diffusivity coef. and insure lateral bc:
      !!        zlt   = ahtt * zlt
      !!        call to lbc_lnk
      !!      Bilaplacian (laplacian of zlt):
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ e2u*e3u/e1u di(zlt) ]
      !!                                  + dj-1[ e1v*e3v/e2v dj(zlt) ]  }
      !!      Note: if key_zco defined, e3t=e3u=e3v, they are simplified.
      !!
      !!      Add this trend to the general trend (ta,sa):
      !!         (ta,sa) = (ta,sa) + ( difft , diffs )
      !!
      !! ** Action : - Update (ta,sa) arrays with the before iso-level
      !!               biharmonic mixing trend.
      !!
      !! History :
      !!        !  91-11  (G. Madec)  Original code
      !!        !  93-03  (M. Guyon)  symetrical conditions
      !!        !  95-11  (G. Madec)  suppress volumetric scale factors
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  96-01  (M. Imbard)  mpp exchange
      !!        !  97-07  (G. Madec)  optimization, and ahtt
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. talandier) New trends organization
      !!        !  05-11  (G. Madec)  zps or sco as default option
      !!----------------------------------------------------------------------
      !! History of the tangent routine
      !!   9.0  !  10-05 (P.A. Bouttier) tangent of 9.0
      !!----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000           ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype           ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt             ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,    kjpt), INTENT(inout) ::   pgu_ad, pgv_ad   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   ptb_ad           ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_ad           ! tracer trend

      !! * Local declarations
      INTEGER ::   ji, jj, jk, jn             ! dummy loop indices
      INTEGER ::   iku, ikv               ! temporary integers
      REAL(wp) ::  ztaad, zsaad, ztmp, zbtr     ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   &
         zeeu, zeev,              & ! 2D workspace
         zltad
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   &
         ztuad, ztvad                         ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_ldf_bilap_adj')
      !
      CALL wrk_alloc( jpi, jpj, zeeu, zeev, zltad )
      CALL wrk_alloc( jpi, jpj, jpk, ztuad, ztvad )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_bilap_adj : iso-level biharmonic operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~'
      ENDIF
      ztaad = 0.e0_wp
      zsaad = 0.e0_wp
      zltad(:,:)   = 0.e0_wp
      ztuad(:,:,:) = 0.e0_wp
      ztvad(:,:,:) = 0.e0_wp

      DO jn = 1, kjpt
         DO jk = jpkm1, 1, -1
               ! 0. Initialization of metric arrays (for z- or s-coordinates)
               ! ----------------------------------
                  DO jj = jpjm1, 1, -1
                     DO ji = jpim1, 1 ,-1   ! vector opt.
                        zeeu(ji,jj) = e2u(ji,jj) * e3u(ji,jj,jk) / e1u(ji,jj) * umask(ji,jj,jk)
                        zeev(ji,jj) = e1v(ji,jj) * e3v(ji,jj,jk) / e2v(ji,jj) * vmask(ji,jj,jk)
                     END DO
                  END DO
               !
               ! 2. Bilaplacian
               ! --------------
               ! fourth derivative (divergence) and add to the general tracer trend
               DO jj = jpjm1, 2, -1
                  DO ji = jpim1, 2, -1   ! vector opt.
                     zbtr = 1.0 / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                     ! add it to the general tracer trends
                     ztaad = pta_ad(ji,jj,jk,jn) * zbtr
                     ! horizontal diffusive trends
                     ztuad(ji  ,jj  ,jk) = ztuad(ji  ,jj  ,jk) + ztaad
                     ztuad(ji-1,jj  ,jk) = ztuad(ji-1,jj  ,jk) - ztaad
                     ztvad(ji  ,jj  ,jk) = ztvad(ji  ,jj  ,jk) + ztaad
                     ztvad(ji  ,jj-1,jk) = ztvad(ji  ,jj-1,jk) - ztaad
                  END DO
               END DO

               ! third derivative (gradient)
               DO jj = jpjm1, 1, -1
                  DO ji = jpim1, 1 ,-1   ! vector opt.
                     ztmp = zeev(ji,jj) * ztvad(ji,jj,jk)
                     zltad(ji  ,jj+1) = zltad(ji  ,jj+1) + ztmp
                     zltad(ji  ,jj  ) = zltad(ji  ,jj  ) - ztmp
                     ztmp = zeeu(ji,jj) * ztuad(ji,jj,jk)
                     zltad(ji+1,jj  ) = zltad(ji+1,jj  ) + ztmp
                     zltad(ji  ,jj  ) = zltad(ji  ,jj  ) - ztmp
                     ztuad(ji  ,jj  ,jk) = 0.0_wp
                     ztvad(ji  ,jj  ,jk) = 0.0_wp
                  END DO
               END DO
               ! Lateral boundary conditions on the laplacian (zltad,zlsad)   (unchanged sgn)
               CALL lbc_lnk_adj( zltad, 'T', 1.0_wp )
               ! Multiply by the eddy diffusivity coefficient
               DO jj = jpjm1, 2, -1
                  DO ji = jpim1, 2, -1   ! vector opt.
                     zltad(ji,jj) = rldf * ahtt(ji,jj) * zltad(ji,jj)
                  END DO
               END DO
               ! Second derivative (divergence)
               DO jj = jpjm1, 2, -1
                  DO ji = jpim1, 2, -1   ! vector opt.
                     zbtr = 1.0 / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                     ztmp = zbtr * zltad(ji,jj)
                     ztvad(ji  ,jj-1,jk) = ztvad(ji  ,jj-1,jk) - ztmp
                     ztvad(ji  ,jj  ,jk) = ztvad(ji  ,jj  ,jk) + ztmp
                     ztuad(ji-1,jj  ,jk) = ztuad(ji-1,jj  ,jk) - ztmp
                     ztuad(ji  ,jj  ,jk) = ztuad(ji  ,jj  ,jk) + ztmp
                     zltad(ji,jj) = 0.0_wp
                  END DO
               END DO
               IF( ln_zps ) THEN      ! set gradient at partial step level
               DO jj = jpjm1, 1, -1
                  DO ji = jpim1, 1 ,-1   ! vector opt.
                        ! last level
                        IF( mbku(ji,jj) == jk ) THEN
                           pgu_ad(ji,jj,jn)   = pgu_ad(ji,jj,jn) + zeeu(ji,jj) * ztuad(ji,jj,jk)
                           ztuad(ji,jj,jk) = 0.0_wp
                        ENDIF
                        IF( mbkv(ji,jj) == jk ) THEN
                           pgv_ad(ji,jj,jn)   = pgv_ad(ji,jj,jn) + zeev(ji,jj) * ztvad(ji,jj,jk)
                           ztvad(ji,jj,jk) = 0.0_wp
                        ENDIF
                     END DO
                  END DO
               ENDIF
               ! 1. Laplacian
               ! ------------

               ! First derivative (gradient)
               DO jj = jpjm1, 1, -1
                  DO ji = jpim1, 1 ,-1   ! vector opt.
                     ztmp = zeev(ji,jj) * ztvad(ji,jj,jk)
                     ptb_ad(ji  ,jj  ,jk,jn) = ptb_ad(ji  ,jj  ,jk,jn) - ztmp
                     ptb_ad(ji  ,jj+1,jk,jn) = ptb_ad(ji  ,jj+1,jk,jn) + ztmp
                     ztmp = zeeu(ji,jj) * ztuad(ji,jj,jk)
                     ptb_ad(ji  ,jj  ,jk,jn) = ptb_ad(ji  ,jj  ,jk,jn) - ztmp
                     ptb_ad(ji+1,jj  ,jk,jn) = ptb_ad(ji+1,jj  ,jk,jn) + ztmp
                     ztuad(ji,jj,jk) = 0.0_wp
                     ztvad(ji,jj,jk) = 0.0_wp
                  END DO
               END DO
          END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, zeeu, zeev, zltad )
      CALL wrk_dealloc( jpi, jpj, jpk, ztuad, ztvad )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_ldf_bilap_adj')
      !
   END SUBROUTINE tra_ldf_bilap_adj

   SUBROUTINE tra_ldf_bilap_adj_tst ( kumadt )
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
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              ztb_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsb_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsa_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  z2r, 'U', 0.0_wp, stds )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgtu_tlin(ji,jj) = z2r(ji,jj)
         END DO
      END DO

      CALL grid_random(  z2r, 'U', 0.0_wp, stds )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgsu_tlin(ji,jj) = z2r(ji,jj)
         END DO
      END DO

      CALL grid_random(  z2r, 'V', 0.0_wp, stds )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgtv_tlin(ji,jj) = z2r(ji,jj)
         END DO
      END DO

      CALL grid_random(  z2r, 'V', 0.0_wp, stds )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgsv_tlin(ji,jj) = z2r(ji,jj)
         END DO
      END DO

      tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
      tsb_tl(:,:,:,jp_sal) = zsb_tlin(:,:,:)
      tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
      tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)
      gtsu_tl(:,:,jp_tem)  = zgtu_tlin(:,:)
      gtsu_tl(:,:,jp_sal)  = zgsu_tlin(:,:)
      gtsv_tl(:,:,jp_tem)  = zgtv_tlin(:,:)
      gtsv_tl(:,:,jp_sal)  = zgsv_tlin(:,:)

      CALL tra_ldf_bilap_tan( nit000, nit000, 'TRA', gtsu_tl, gtsv_tl,tsb_tl, tsa_tl, jpts  )

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

      CALL tra_ldf_bilap_adj( nit000 , nit000, 'TRA', gtsu_ad, gtsv_ad,tsb_ad, tsa_ad, jpts)

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

      cl_name = 'tra_ldf_bilap'
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


   END SUBROUTINE tra_ldf_bilap_adj_tst

END MODULE traldf_bilap_tam
