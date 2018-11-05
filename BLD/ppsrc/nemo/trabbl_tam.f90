MODULE trabbl_tam
   !!==============================================================================
   !!                       ***  MODULE  trabbl  ***
   !! Ocean physics :  advective and/or diffusive bottom boundary layer scheme
   !!==============================================================================
   !! History :  OPA  ! 1996-06  (L. Mortier)  Original code
   !!            8.0  ! 1997-11  (G. Madec)    Optimization
   !!   NEMO     1.0  ! 2002-08  (G. Madec)  free form + modules
   !!             -   ! 2004-01  (A. de Miranda, G. Madec, J.M. Molines ) add advective bbl
   !!            3.3  ! 2009-11  (G. Madec)  merge trabbl and trabbl_adv + style + optimization
   !!             -   ! 2010-04  (G. Madec)  Campin & Goosse advective bbl
   !!             -   ! 2010-06  (C. Ethe, G. Madec)  merge TRA-TRC
   !!             -   ! 2010-11  (G. Madec) add mbk. arrays associated to the deepest ocean level
   !! History of the T&A module
   !!   NEMO     3.2  ! 2011-02  (A. Vidard) Original version
   !!            3.4  ! 2012-09  (A. Vidard) Update to 3.4
   !!
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_trabbl'   or                             bottom boundary layer
   !!----------------------------------------------------------------------
   !!   tra_bbl_alloc : allocate trabbl arrays
   !!   tra_bbl       : update the tracer trends due to the bottom boundary layer (advective and/or diffusive)
   !!   tra_bbl_dif   : generic routine to compute bbl diffusive trend
   !!   tra_bbl_adv   : generic routine to compute bbl advective trend
   !!   bbl           : computation of bbl diffu. flux coef. & transport in bottom boundary layer
   !!   tra_bbl_init  : initialization, namelist read, parameters control
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and active tracers
   USE oce_tam
   USE dom_oce        ! ocean space and time domain
   USE phycst         ! physical constant
   USE eosbn2         ! equation of state
   USE iom            ! IOM server
   USE in_out_manager ! I/O manager
   USE lbclnk         ! ocean lateral boundary conditions
   USE prtctl         ! Print control
   USE wrk_nemo       ! Memory Allocation
   USE timing         ! Timing
   USE trabbl
   USE gridrandom
   USE dotprodfld
   USE tstool_tam

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_bbl_tan       !  routine called by step.F90
   PUBLIC   tra_bbl_init_tam  !  routine called by opa.F90
   PUBLIC   tra_bbl_dif_tan   !  routine called by trcbbl.F90
   PUBLIC   tra_bbl_adv_tan   !  -          -          -              -
   PUBLIC   bbl_tan           !  routine called by trcbbl.F90 and dtadyn.F90
   PUBLIC   tra_bbl_adj       !  routine called by step.F90
   PUBLIC   tra_bbl_dif_adj   !  routine called by trcbbl.F90
   PUBLIC   tra_bbl_adv_adj   !  -          -          -              -
   PUBLIC   bbl_adj           !  routine called by trcbbl.F90 and dtadyn.F90
   PUBLIC   bbl_adj_tst       !  routine called by tamtst
   PUBLIC   tra_bbl_adj_tst   !  routine called by tamtst

   REAL(WP), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: utr_bbl_tl, vtr_bbl_tl
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: ahu_bbl_tl, ahv_bbl_tl
   REAL(WP), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: utr_bbl_ad, vtr_bbl_ad
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:), PUBLIC :: ahu_bbl_ad, ahv_bbl_ad

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ahu_bbl_0_tl, ahv_bbl_0_tl
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION(:,:) ::   ahu_bbl_0_ad, ahv_bbl_0_ad

   LOGICAL, PRIVATE :: ll_alloctl = .FALSE., ll_allocad = .FALSE.

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION tra_bbl_alloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_bbl_alloc_tam  ***
      !!----------------------------------------------------------------------
      INTEGER, OPTIONAL :: kmode
      INTEGER, DIMENSION(2) :: ierr
      INTEGER :: jmode

      IF ( PRESENT( kmode ) ) THEN
         jmode = kmode
      ELSE
         jmode = 0
      END IF
      ierr(:) = 0.0_wp

      IF ( ( jmode == 0 ) .OR. ( jmode == 1 ) .AND. ( .NOT. ll_alloctl ) ) THEN
         ALLOCATE( utr_bbl_tl  (jpi,jpj), vtr_bbl_tl  (jpi,jpj), &
            &      ahu_bbl_tl  (jpi,jpj), ahv_bbl_tl  (jpi,jpj), &
            &      ahu_bbl_0_tl  (jpi,jpj), ahv_bbl_0_tl  (jpi,jpj), &
            &      STAT= ierr(1) )
         ll_alloctl = .TRUE.
      END IF
         !
      IF ( ( jmode == 0 ) .OR. ( jmode == 2 ) .AND. ( .NOT. ll_allocad ) ) THEN
         ALLOCATE( utr_bbl_ad  (jpi,jpj), vtr_bbl_ad  (jpi,jpj), &
            &      ahu_bbl_ad  (jpi,jpj), ahv_bbl_ad  (jpi,jpj), &
            &      ahu_bbl_0_ad  (jpi,jpj), ahv_bbl_0_ad  (jpi,jpj), &
            &      STAT= ierr(2) )
         ll_allocad = .TRUE.
      END IF
      tra_bbl_alloc_tam = SUM(ierr)
         !
      IF( lk_mpp            )   CALL mpp_sum ( tra_bbl_alloc_tam )
      IF( tra_bbl_alloc_tam > 0 )   CALL ctl_warn('tra_bbl_alloc_tam: allocation of arrays failed.')
   END FUNCTION tra_bbl_alloc_tam


   INTEGER FUNCTION tra_bbl_dealloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_bbl_dealloc  ***
      !!----------------------------------------------------------------------
      INTEGER, OPTIONAL :: kmode
      INTEGER, DIMENSION(2) :: ierr

      IF ( .NOT. PRESENT( kmode ) ) kmode=0
      ierr(:) = 0.0_wp

      IF ( ( kmode == 0 ) .OR. ( kmode == 1 ) .AND. ( ll_alloctl ) ) THEN
         DEALLOCATE( utr_bbl_tl, vtr_bbl_tl, &
            &      ahu_bbl_0_tl, ahv_bbl_0_tl, &
            &      STAT= ierr(1) )
         ll_alloctl = .FALSE.
      END IF
         !
      IF ( ( kmode == 0 ) .OR. ( kmode == 1 ) .AND. ( ll_allocad ) ) THEN
         DEALLOCATE( utr_bbl_ad, vtr_bbl_ad, &
            &      ahu_bbl_ad, ahv_bbl_ad, &
            &      ahu_bbl_0_ad, ahv_bbl_0_ad, &
            &      STAT= ierr(2) )
         ll_allocad = .FALSE.
      END IF
      tra_bbl_dealloc_tam = SUM(ierr)
         !
      IF( lk_mpp            )   CALL mpp_sum ( tra_bbl_dealloc_tam )
      IF( tra_bbl_dealloc_tam > 0 )   CALL ctl_warn('tra_bbl_dealloc_tam: allocation of arrays failed.')
   END FUNCTION tra_bbl_dealloc_tam


   SUBROUTINE tra_bbl_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl_tan  ***
      !!
      !! ** Purpose :   Compute the before tracer (t & s) trend associated
      !!              with the bottom boundary layer and add it to the general
      !!              trend of tracer equations.
      !!
      !! ** Method  :   Depending on namtra_bbl namelist parameters the bbl
      !!              diffusive and/or advective contribution to the tracer trend
      !!              is added to the general tracer trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step
      !!
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdttl, ztrdstl
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl_tan')
      !
      IF( l_bbl )  THEN 
         CALL bbl( kt, nit000, 'TRA' )
         CALL bbl_tan( kt, nit000, 'TRA' )   !* bbl coef. and transport (only if not already done in trcbbl)
      END IF

      IF( nn_bbl_ldf == 1 ) THEN                   !* Diffusive bbl
         !
         CALL tra_bbl_dif_tan( tsb, tsb_tl, tsa_tl, jpts )
         !
      END IF

      IF( nn_bbl_adv /= 0 ) THEN                !* Advective bbl
         !
         CALL tra_bbl_adv_tan( tsb, tsb_tl, tsa_tl, jpts )
         !
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl_tan')
      !
   END SUBROUTINE tra_bbl_tan


   SUBROUTINE tra_bbl_dif_tan( ptb, ptb_tl, pta_tl, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_dif_tan  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  :
      !!        * diffusive bbl (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!
      !! ** Action  :   pta   increased by the bbl diffusive trend
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      !
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb    ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb_tl    ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_tl    ! tracer trend
      !
      INTEGER  ::   ji, jj, jn   ! dummy loop indices
      INTEGER  ::   ik           ! local integers
      REAL(wp) ::   zbtr         ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) :: zptb, zptbtl
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_bbl_dif_tan')
      !
      CALL wrk_alloc( jpi, jpj, zptb, zptbtl )
      !
      DO jn = 1, kjpt                                     ! tracer loop
         !                                                ! ===========
         DO jj = 1, jpj
            DO ji = 1, jpi
               ik = mbkt(ji,jj)                        ! bottom T-level index
               zptbtl(ji,jj) = ptb_tl(ji,jj,ik,jn)              ! bottom before T and S
               zptb  (ji,jj) = ptb   (ji,jj,ik,jn)     ! bottom before T and S
            END DO
         END DO
         !                                                ! Compute the trend
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ik = mbkt(ji,jj)                            ! bottom T-level index
               zbtr = r1_e1e2t(ji,jj)  / e3t(ji,jj,ik)
               pta_tl(ji,jj,ik,jn) = pta_tl(ji,jj,ik,jn)                                                         &
                  &               + (   ahu_bbl(ji  ,jj  ) * ( zptbtl(ji+1,jj  ) - zptbtl(ji  ,jj  ) )   &
                  &                   - ahu_bbl(ji-1,jj  ) * ( zptbtl(ji  ,jj  ) - zptbtl(ji-1,jj  ) )   &
                  &                   + ahv_bbl(ji  ,jj  ) * ( zptbtl(ji  ,jj+1) - zptbtl(ji  ,jj  ) )   &
                  &                   - ahv_bbl(ji  ,jj-1) * ( zptbtl(ji  ,jj  ) - zptbtl(ji  ,jj-1) )   ) * zbtr
            END DO
         END DO
         !                                                  ! ===========
      END DO                                                ! end tracer
      !                                                     ! ===========
      CALL wrk_dealloc( jpi, jpj, zptbtl, zptb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_bbl_dif_tan')
      !
   END SUBROUTINE tra_bbl_dif_tan


   SUBROUTINE tra_bbl_adv_tan( ptb, ptb_tl, pta_tl, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_bbl  ***
      !!
      !! ** Purpose :   Compute the before passive tracer trend associated
      !!     with the bottom boundary layer and add it to the general trend
      !!     of tracer equations.
      !! ** Method  :   advective bbl (nn_bbl_adv = 1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean near bottom velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation i.e.
      !!                       transport proportional to the along-slope density gradient
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb    ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb_tl ! before and now tangent tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_tl    ! tracer trend
      !
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      INTEGER  ::   iis , iid , ijs , ijd    ! local integers
      INTEGER  ::   ikus, ikud, ikvs, ikvd   !   -       -
      REAL(wp) ::   zbtr, ztra               ! local scalars
      REAL(wp) ::   ztratl                   !   -      -
      REAL(wp) ::   zu_bbl, zv_bbl           !   -      -
      REAL(wp) ::   zu_bbltl, zv_bbltl       !   -      -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl_adv_tan')
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         DO jj = 1, jpjm1
            DO ji = 1, jpim1            ! CAUTION start from i=1 to update i=2 when cyclic east-west
               IF( utr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero i-direction bbl advection
                  ! down-slope i/k-indices (deep)      &   up-slope i/k indices (shelf)
                  iid  = ji + MAX( 0, mgrhu(ji,jj) )   ;   iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
                  ikud = mbku_d(ji,jj)                 ;   ikus = mbku(ji,jj)
                  zu_bbl = ABS( utr_bbl(ji,jj) )
                  zu_bbltl = SIGN( utr_bbl_tl(ji,jj), utr_bbl(ji,jj) )
                  !
                  !                                               ! up  -slope T-point (shelf bottom point)
                  zbtr = r1_e1e2t(iis,jj) / e3t(iis,jj,ikus)
                  ztratl = ( zu_bbltl * ( ptb   (iid,jj,ikus,jn) - ptb   (iis,jj,ikus,jn) )   &
                     &     + zu_bbl   * ( ptb_tl(iid,jj,ikus,jn) - ptb_tl(iis,jj,ikus,jn) ) ) * zbtr
                  pta_tl(iis,jj,ikus,jn) = pta_tl(iis,jj,ikus,jn) + ztratl
                  !
                  DO jk = ikus, ikud-1                            ! down-slope upper to down T-point (deep column)
                     zbtr = r1_e1e2t(iid,jj) / e3t(iid,jj,jk)
                     ztratl = ( zu_bbltl * ( ptb   (iid,jj,jk+1,jn) - ptb   (iid,jj,jk,jn) ) &
                        &     + zu_bbl   * ( ptb_tl(iid,jj,jk+1,jn) - ptb_tl(iid,jj,jk,jn) ) ) * zbtr
                     pta_tl(iid,jj,jk,jn) = pta_tl(iid,jj,jk,jn) + ztratl
                  END DO
                  !
                  zbtr = r1_e1e2t(iid,jj) / e3t(iid,jj,ikud)
                  ztratl = ( zu_bbltl * ( ptb   (iis,jj,ikus,jn) - ptb   (iid,jj,ikud,jn) ) &
                     &     + zu_bbl   * ( ptb_tl(iis,jj,ikus,jn) - ptb_tl(iid,jj,ikud,jn) ) ) * zbtr
                  pta_tl(iid,jj,ikud,jn) = pta_tl(iid,jj,ikud,jn) + ztratl
               ENDIF
               !
               IF( vtr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero j-direction bbl advection
                  ! down-slope j/k-indices (deep)        &   up-slope j/k indices (shelf)
                  ijd  = jj + MAX( 0, mgrhv(ji,jj) )     ;   ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
                  ikvd = mbkv_d(ji,jj)                   ;   ikvs = mbkv(ji,jj)
                  zv_bbl = ABS( vtr_bbl(ji,jj) )
                  zv_bbltl = SIGN( vtr_bbl_tl(ji,jj), vtr_bbl(ji,jj) )
                  !
                  ! up  -slope T-point (shelf bottom point)
                  zbtr = r1_e1e2t(ji,ijs) / e3t(ji,ijs,ikvs)
                  ztratl = ( zv_bbltl * ( ptb   (ji,ijd,ikvs,jn) - ptb   (ji,ijs,ikvs,jn) ) &
                     &     + zv_bbl   * ( ptb_tl(ji,ijd,ikvs,jn) - ptb_tl(ji,ijs,ikvs,jn) ) ) * zbtr
                  pta_tl(ji,ijs,ikvs,jn) = pta_tl(ji,ijs,ikvs,jn) + ztratl
                  !
                  DO jk = ikvs, ikvd-1                            ! down-slope upper to down T-point (deep column)
                     zbtr = r1_e1e2t(ji,ijd) / e3t(ji,ijd,jk)
                     ztratl = ( zv_bbltl * ( ptb   (ji,ijd,jk+1,jn) - ptb   (ji,ijd,jk,jn) ) &
                        &     + zv_bbl   * ( ptb_tl(ji,ijd,jk+1,jn) - ptb_tl(ji,ijd,jk,jn) ) ) * zbtr
                     pta_tl(ji,ijd,jk,jn) = pta_tl(ji,ijd,jk,jn)  + ztratl
                  END DO
                  !                                               ! down-slope T-point (deep bottom point)
                  zbtr = r1_e1e2t(ji,ijd) / e3t(ji,ijd,ikvd)
                  ztratl = ( zv_bbltl * ( ptb   (ji,ijs,ikvs,jn) - ptb   (ji,ijd,ikvd,jn) ) &
                     &     + zv_bbl   * ( ptb_tl(ji,ijs,ikvs,jn) - ptb_tl(ji,ijd,ikvd,jn) ) ) * zbtr
                  pta_tl(ji,ijd,ikvd,jn) = pta_tl(ji,ijd,ikvd,jn) + ztratl
               ENDIF
            END DO
            !
         END DO
         !                                                  ! ===========
      END DO                                                ! end tracer
      !                                                     ! ===========
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl_adv_tan')
      !
   END SUBROUTINE tra_bbl_adv_tan


   SUBROUTINE bbl_tan( kt, kit000, cdtype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  :
      !!        * diffusive bbl (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!        * advective bbl (nn_bbl_adv=1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation
      !!        i.e. transport proportional to the along-slope density gradient
      !!
      !!      NB: the along slope density gradient is evaluated using the
      !!      local density (i.e. referenced at a common local depth).
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      !
      INTEGER         , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER         , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3), INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      !!
      INTEGER  ::   ji, jj                    ! dummy loop indices
      INTEGER  ::   ik                        ! local integers
      INTEGER  ::   iis , iid , ijs , ijd     !   -       -
      INTEGER  ::   ikus, ikud, ikvs, ikvd    !   -       -
      REAL(wp) ::   zsign, zsigna, zgbbl      ! local scalars
      REAL(wp) ::   zgdrho, zt, zs, zh        !   -      -
      REAL(wp) ::   zgdrhotl, zttl, zstl, zhtl!   -      -
      !!
      REAL(wp) ::   fsalbt, fsbeta, pft, pfs, pfh   ! statement function
      REAL(wp) ::   fsalbt_tan, fsbeta_tan, pfttl, pfstl, pfhtl   ! statement function
      REAL(wp), POINTER, DIMENSION(:,:) :: zub  , zvb  , ztb  , zsb  , zdep
      REAL(wp), POINTER, DIMENSION(:,:) :: zubtl, zvbtl, ztbtl, zsbtl
      !!----------------------- zv_bbl-----------------------------------------------
      ! ratio alpha/beta = fsalbt : ratio of thermal over saline expension coefficients
      ! ================            pft :  potential temperature in degrees celcius
      !                             pfs :  salinity anomaly (s-35) in psu
      !                             pfh :  depth in meters
      ! nn_eos = 0  (Jackett and McDougall 1994 formulation)
      fsalbt( pft, pfs, pfh ) =                                              &   ! alpha/beta
         ( ( ( -0.255019e-07 * pft + 0.298357e-05 ) * pft                    &
                                   - 0.203814e-03 ) * pft                    &
                                   + 0.170907e-01 ) * pft                    &
                                   + 0.665157e-01                            &
         +(-0.678662e-05 * pfs - 0.846960e-04 * pft + 0.378110e-02 ) * pfs   &
         +  ( ( - 0.302285e-13 * pfh                                         &
                - 0.251520e-11 * pfs                                         &
                + 0.512857e-12 * pft * pft          ) * pfh                  &
                                     - 0.164759e-06   * pfs                  &
             +(   0.791325e-08 * pft - 0.933746e-06 ) * pft                  &
                                     + 0.380374e-04 ) * pfh
      fsbeta( pft, pfs, pfh ) =                                              &   ! beta
         ( ( -0.415613e-09 * pft + 0.555579e-07 ) * pft                      &
                                 - 0.301985e-05 ) * pft                      &
                                 + 0.785567e-03                              &
         + (     0.515032e-08 * pfs                                          &
               + 0.788212e-08 * pft - 0.356603e-06 ) * pfs                   &
               +(  (   0.121551e-17 * pfh                                    &
                     - 0.602281e-15 * pfs                                    &
                     - 0.175379e-14 * pft + 0.176621e-12 ) * pfh             &
                                          + 0.408195e-10   * pfs             &
                 + ( - 0.213127e-11 * pft + 0.192867e-09 ) * pft             &
                                          - 0.121555e-07 ) * pfh

      fsalbt_tan( pft, pfs, pfh, pfttl, pfstl, pfhtl ) =  &   ! alpha/beta
         &         ( - 0.255019e-07 * 4 * pft * pft * pft            &
         &           + 0.298357e-05 * 3 * pft * pft                  &
         &           - 0.203814e-03 * 2 * pft                        &
         &           - 0.846960e-04 * pfs                            &
         &           + 0.512857e-12 * 2 * pft * pfh * pfh            &
         &           + 0.791325e-08 * pft * pfh                      &
         &           - 0.933746e-06 * pfh                            &
         &           + 0.170907e-01                        ) * pfttl &
         &       + ( - 0.678662e-05 * 2 * pfs                        &
         &           - 0.846960e-04 * pft                            &
         &           - 0.251520e-11 * pfh  * pfh                     &
         &           - 0.164759e-06 * pfh                            &
         &           + 0.378110e-02                        ) * pfstl &
         &       + ( - 0.302285e-13 * 3 * pfh  * pfh                 &
         &           - 0.251520e-11 * pfs * pfh                      &
         &           + 0.512857e-12 * pft * pft * pfh                &
         &           - 0.164759e-06 * pfs                            &
         &           + 0.791325e-08 * pft * pft                      &
         &           - 0.933746e-06 * pft                            &
         &           + 0.380374e-04                        ) * pfhtl


      fsbeta_tan( pft, pfs, pfh, pfttl, pfstl, pfhtl ) =  &   ! beta
         &         ( - 0.415613e-09 * 3 * pft * pft                  &
         &           + 0.555579e-07 * 2 * pft                        &
         &           - 0.301985e-05                                  &
         &           + 0.788212e-08 * pfs                            &
         &           - 0.213127e-11 * 2 * pfh * pft                  &
         &           - 0.175379e-14 * pfh * pfh            ) * pfttl &
         &       + (   0.788212e-08 * pft                            &
         &           + 0.515032e-08 * 2 * pfs                        &
         &           - 0.356603e-06                                  &
         &           + 0.408195e-10 * pfh                            &
         &           - 0.602281e-15 * pfh * pfh            ) * pfstl &
         &       + (   0.121551e-17 * 3 * pfh * pfh                  &
         &           - 0.602281e-15 * 2 * pfs * pfh                  &
         &           - 0.175379e-14 * 2 * pft * pfh                  &
         &           + 0.176621e-12 * 2 * pfh                        &
         &           + 0.408195e-10 * pfs                            &
         &           + 0.192867e-09 * pfh                            &
         &           - 0.213127e-11 * pft * pft                      &
         &           + 0.192867e-09 * pft                            &
         &           - 0.121555e-07                        ) * pfhtl


      !!----------------------------------------------------------------------

      !
      IF( nn_timing == 1 )  CALL timing_start( 'bbl_tan')
      !
      CALL wrk_alloc( jpi, jpj, zub  , zvb  , ztb  , zsb  , zdep, &
         &                      zubtl, zvbtl, ztbtl, zsbtl        )
      !

      IF( kt == kit000 )  THEN
         IF(lwp)  WRITE(numout,*)
         IF(lwp)  WRITE(numout,*) 'trabbl_tam:bbl_tan : Compute bbl velocities and diffusive coefficients in ', cdtype
         IF(lwp)  WRITE(numout,*) '~~~~~~~~~~'
      ENDIF

      !                                        !* bottom temperature, salinity, velocity and depth
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = mbkt(ji,jj)                        ! bottom T-level index
            ztb  (ji,jj) = tsb(ji,jj,ik,jp_tem) * tmask(ji,jj,1)      ! bottom before T and S
            zsb  (ji,jj) = tsb(ji,jj,ik,jp_sal) * tmask(ji,jj,1)
            ztbtl(ji,jj) = tsb_tl(ji,jj,ik,jp_tem) * tmask(ji,jj,1)      ! bottom before T and S
            zsbtl(ji,jj) = tsb_tl(ji,jj,ik,jp_sal) * tmask(ji,jj,1)
            zdep(ji,jj) = gdept(ji,jj,ik)        ! bottom T-level reference depth
            !
            zub(ji,jj) = un(ji,jj,mbku(ji,jj))      ! bottom velocity
            zvb(ji,jj) = vn(ji,jj,mbkv(ji,jj))

            zubtl(ji,jj) = un_tl(ji,jj,mbku(ji,jj))      ! bottom velocity
            zvbtl(ji,jj) = vn_tl(ji,jj,mbkv(ji,jj))
         END DO
      END DO

      !                                   !-------------------!
      IF( nn_bbl_ldf == 1 ) THEN          !   diffusive bbl   !
         !                                !-------------------!
         ! AV NOTE : while rn_ahtbbl remains a passive variable, the code below will only yield ah_bbl_tl=0, so i put it under key
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               ahu_bbl_tl(ji,jj)=0.0_wp
               ahv_bbl_tl(ji,jj)=0.0_wp
            END DO
         END DO
         !
      ENDIF

      !                                   !-------------------!
      IF( nn_bbl_adv /= 0 ) THEN          !   advective bbl   !
         !                                !-------------------!
         SELECT CASE ( nn_bbl_adv )             !* bbl transport type
         !
         CASE( 1 )                                   != use of upper velocity
            ! AV NOTE: not much needed for deriving, almost all the computations are for the SIGN, which is kept identical as in the NL
            DO jj = 1, jpjm1                                 ! criteria: grad(rho).grad(h)<0  and grad(rho).grad(h)<0
               DO ji = 1, jpim1   ! vector opt.
                  !                                               ! i-direction
                  zt = 0.5 * ( ztb (ji,jj) + ztb (ji+1,jj) )                  ! T, S anomalie, and depth
                  zs = 0.5 * ( zsb (ji,jj) + zsb (ji+1,jj) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )
                  !                                                           ! masked bbl i-gradient of density
                  zgdrho = (  fsalbt( zt, zs, zh ) * ( ztb(ji+1,jj) - ztb(ji,jj) )    &
                     &                             - ( zsb(ji+1,jj) - zsb(ji,jj) )  ) * umask(ji,jj,1)
                  !
                  zsign = SIGN(  0.5, - zgdrho   * REAL( mgrhu(ji,jj) )  )    ! sign of i-gradient * i-slope
                  zsigna= SIGN(  0.5, zub(ji,jj) * REAL( mgrhu(ji,jj) )  )    ! sign of u * i-slope
                  !
                  !                                                           ! bbl velocity
                  utr_bbl_tl(ji,jj) = ( 0.5 + zsigna ) * ( 0.5 - zsign ) * e2u(ji,jj) * e3u_bbl_0(ji,jj) * zubtl(ji,jj)
                  !
                  !                                               ! j-direction
                  zt = 0.5 * ( ztb (ji,jj+1) + ztb (ji,jj) )                  ! T, S anomalie, and depth
                  zs = 0.5 * ( zsb (ji,jj+1) + zsb (ji,jj) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj+1) + zdep(ji,jj) )
                  !                                                           ! masked bbl j-gradient of density
                  zgdrho = (  fsalbt( zt, zs, zh ) * ( ztb(ji,jj+1) - ztb(ji,jj) )    &
                     &                             - ( zsb(ji,jj+1) - zsb(ji,jj) )  ) * vmask(ji,jj,1)
                  zsign = SIGN(  0.5, - zgdrho   * REAL( mgrhv(ji,jj) )  )    ! sign of j-gradient * j-slope
                  zsigna= SIGN(  0.5, zvb(ji,jj) * REAL( mgrhv(ji,jj) )  )    ! sign of u * i-slope
                  !
                  !                                                           ! bbl velocity
                  vtr_bbl_tl(ji,jj) = ( 0.5 + zsigna ) * ( 0.5 - zsign ) * e1v(ji,jj) * e3v_bbl_0(ji,jj) * zvbtl(ji,jj)
               END DO
            END DO
            !
         CASE( 2 )                                 != bbl velocity = F( delta rho )
            ! AV NOTE: this one is nastier
            zgbbl = grav * rn_gambbl
            DO jj = 1, jpjm1                            ! criteria: rho_up > rho_down
               DO ji = 1, jpim1   ! vector opt.
                  !                                         ! i-direction
                  ! down-slope T-point i/k-index (deep)  &   up-slope T-point i/k-index (shelf)
                  iid  = ji + MAX( 0, mgrhu(ji,jj) )     ;    iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
                  ikud = mbku_d(ji,jj)                   ;    ikus = mbku(ji,jj)
                  !
                  !                                             ! mid-depth density anomalie (up-slope minus down-slope)
                  zt = 0.5 * ( ztb (ji,jj) + ztb (ji+1,jj) )           ! mid slope depth of T, S, and depth
                  zs = 0.5 * ( zsb (ji,jj) + zsb (ji+1,jj) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )
                  zgdrho =    fsbeta( zt, zs, zh )                                    &
                     &   * (  fsalbt( zt, zs, zh ) * ( ztb(iid,jj) - ztb(iis,jj) )    &
                     &                             - ( zsb(iid,jj) - zsb(iis,jj) )  ) * umask(ji,jj,1)
                  zttl = 0.5 * ( ztbtl (ji,jj) + ztbtl (ji+1,jj) )           ! mid slope depth of T, S, and depth
                  zstl = 0.5 * ( zsbtl (ji,jj) + zsbtl (ji+1,jj) )
                  zhtl = 0.0_wp
                  zgdrhotl =  ( fsbeta_tan( zt, zs, zh, zttl, zstl, zhtl )                     &
                     &        * ( fsalbt( zt, zs, zh ) * ( ztb(iid,jj) - ztb(iis,jj) )         &
                     &                                 - ( zsb(iid,jj) - zsb(iis,jj) ) )       &
                     &        + fsbeta( zt, zs, zh )                                           &
                     &        * ( fsalbt_tan( zt, zs, zh, zttl, zstl, zhtl )                   &
                     &                                     * ( ztb  (iid,jj) - ztb  (iis,jj) ) &
                     &          + fsalbt    ( zt, zs, zh ) * ( ztbtl(iid,jj) - ztbtl(iis,jj) ) &
                     &                                  - ( zsbtl(iid,jj) - zsbtl(iis,jj) )  ) ) * umask(ji,jj,1)

                  zsign  = SIGN( 0.5_wp, zgdrho ) ! tangent of zgdrho = MAX( 0.e0, zgdrho )
                  !                                             ! bbl transport (down-slope direction)
                  utr_bbl_tl(ji,jj) = zsign * e2u(ji,jj) * e3u_bbl_0(ji,jj) * zgbbl * zgdrhotl * REAL( mgrhu(ji,jj) )
                  !
                  !                                         ! j-direction
                  !  down-slope T-point j/k-index (deep)  &   of the up  -slope T-point j/k-index (shelf)
                  ijd  = jj + MAX( 0, mgrhv(ji,jj) )      ;    ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
                  ikvd = mbkv_d(ji,jj)                    ;    ikvs = mbkv(ji,jj)
                  !
                  !                                             ! mid-depth density anomalie (up-slope minus down-slope)
                  zt = 0.5 * ( ztb (ji,jj) + ztb (ji,jj+1) )           ! mid slope depth of T, S, and depth
                  zs = 0.5 * ( zsb (ji,jj) + zsb (ji,jj+1) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj) + zdep(ji,jj+1) )
                  zgdrho =    fsbeta( zt, zs, zh )                                    &
                     &   * (  fsalbt( zt, zs, zh ) * ( ztb(ji,ijd) - ztb(ji,ijs) )    &
                     &                             - ( zsb(ji,ijd) - zsb(ji,ijs) )  ) * vmask(ji,jj,1)
                  zttl = 0.5 * ( ztbtl (ji,jj) + ztbtl (ji,jj+1) )           ! mid slope depth of T, S, and depth
                  zstl = 0.5 * ( zsbtl (ji,jj) + zsbtl (ji,jj+1) )
                  zhtl = 0.0_wp
                  zgdrhotl =  ( fsbeta_tan( zt, zs, zh, zttl, zstl, zhtl )                     &
                     &        * ( fsalbt( zt, zs, zh ) * ( ztb(ji,ijd) - ztb(ji,ijs) )         &
                     &                                 - ( zsb(ji,ijd) - zsb(ji,ijs) ) )       &
                     &        + fsbeta( zt, zs, zh )                                           &
                     &        * ( fsalbt_tan( zt, zs, zh, zttl, zstl, zhtl )                   &
                     &                                     * ( ztb  (ji,ijd) - ztb  (ji,ijs) ) &
                     &          + fsalbt    ( zt, zs, zh ) * ( ztbtl(ji,ijd) - ztbtl(ji,ijs) ) &
                     &                                     - ( zsbtl(ji,ijd) - zsbtl(ji,ijs) ) ) ) * vmask(ji,jj,1)
                  !
                  zsign  = SIGN( 0.5_wp, zgdrho ) ! tangent of zgdrho = MAX( 0.e0, zgdrho )
                  !                                             ! bbl transport (down-slope direction)
                  vtr_bbl_tl(ji,jj) = zsign * e1v(ji,jj) * e3v_bbl_0(ji,jj) * zgbbl * zgdrhotl * REAL( mgrhv(ji,jj) )
               END DO
            END DO
         END SELECT
         !
      ENDIF
      !
      CALL wrk_dealloc( jpi, jpj, zub  , zvb  , ztb  , zsb  , zdep,   &
         &                        zubtl, zvbtl, ztbtl, zsbtl        )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'bbl_tan')
      !
   END SUBROUTINE bbl_tan


   SUBROUTINE tra_bbl_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl_adj  ***
      !!
      !! ** Purpose :   Compute the before tracer (t & s) trend associated
      !!              with the bottom boundary layer and add it to the general
      !!              trend of tracer equations.
      !!
      !! ** Method  :   Depending on namtra_bbl namelist parameters the bbl
      !!              diffusive and/or advective contribution to the tracer trend
      !!              is added to the general tracer trend
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl_adj')
      !
      IF( l_bbl )  CALL bbl( kt, nitend, 'TRA' )   !* bbl coef. and transport (only if not already done in trcbbl)

      IF( nn_bbl_adv /= 0 ) THEN                !* Advective bbl
         !
         CALL tra_bbl_adv_adj( tsb, tsb_ad, tsa_ad, jpts )
         !
      END IF

      IF( nn_bbl_ldf == 1 ) THEN                   !* Diffusive bbl
         !
         CALL tra_bbl_dif_adj( tsb, tsb_ad, tsa_ad, jpts )
         !
      END IF

      IF( l_bbl )  CALL bbl_adj( kt, nitend, 'TRA' )   !* bbl coef. and transport (only if not already done in trcbbl)
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl_adj')
      !
   END SUBROUTINE tra_bbl_adj


   SUBROUTINE tra_bbl_dif_adj( ptb, ptb_ad, pta_ad, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_dif_adj  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  :
      !!        * diffusive bbl (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!
      !! ** Action  :   pta   increased by the bbl diffusive trend
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      !
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb    ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   ptb_ad ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_ad    ! tracer trend
      !
      INTEGER  ::   ji, jj, jn   ! dummy loop indices
      INTEGER  ::   ik           ! local integers
      REAL(wp) ::   zbtr         ! local scalars
      REAL(wp), POINTER, DIMENSION(:,:) :: zptb, zptbad
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_bbl_dif_adj')
      !
      CALL wrk_alloc( jpi, jpj, zptb, zptbad )
      zptbad(:,:) = 0.0_wp
      !
      DO jn = kjpt, 1, -1                                 ! tracer loop
         !                                                ! ===========
         !                                                  ! ===========
         !                                             ! Compute the trend
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1
               ik = mbkt(ji,jj)                            ! bottom T-level index
               zbtr = r1_e1e2t(ji,jj)  / e3t(ji,jj,ik)
               zptbad(ji  ,jj  ) = zptbad(ji  ,jj  ) - pta_ad(ji,jj,ik,jn) * ( ahu_bbl(ji  ,jj  ) + ahu_bbl(ji-1,jj  ) &
                  &                                                          + ahv_bbl(ji  ,jj  ) + ahv_bbl(ji  ,jj-1) ) * zbtr
               zptbad(ji+1,jj  ) = zptbad(ji+1,jj  ) + pta_ad(ji,jj,ik,jn) * ahu_bbl(ji  ,jj  ) * zbtr
               zptbad(ji-1,jj  ) = zptbad(ji-1,jj  ) + pta_ad(ji,jj,ik,jn) * ahu_bbl(ji-1,jj  ) * zbtr
               zptbad(ji  ,jj+1) = zptbad(ji  ,jj+1) + pta_ad(ji,jj,ik,jn) * ahv_bbl(ji  ,jj  ) * zbtr
               zptbad(ji  ,jj-1) = zptbad(ji  ,jj-1) + pta_ad(ji,jj,ik,jn) * ahv_bbl(ji  ,jj-1) * zbtr

               pta_ad(ji,jj,ik,jn) = pta_ad(ji,jj,ik,jn)                                                       &
                  &                 + ( ahu_bbl(ji  ,jj  ) * ( zptbad(ji+1,jj  ) - zptbad(ji  ,jj  ) )   &
                  &                   - ahu_bbl(ji-1,jj  ) * ( zptbad(ji  ,jj  ) - zptbad(ji-1,jj  ) )   &
                  &                   + ahv_bbl(ji  ,jj  ) * ( zptbad(ji  ,jj+1) - zptbad(ji  ,jj  ) )   &
                  &                   - ahv_bbl(ji  ,jj-1) * ( zptbad(ji  ,jj  ) - zptbad(ji  ,jj-1) )   ) * zbtr
            END DO
         END DO
         DO jj = jpj, 1, -1
            DO ji = jpi, 1, -1
               ik = mbkt(ji,jj)                        ! bottom T-level index
               ptb_ad(ji,jj,ik,jn) = ptb_ad(ji,jj,ik,jn) + zptbad(ji,jj)
               zptbad(ji,jj) = 0.0_wp                  ! bottom before T and S
            END DO
         END DO
         !                                                  ! ===========
      END DO                                                ! end tracer
      !                                                     ! ===========
      CALL wrk_dealloc( jpi, jpj, zptbad, zptb )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_bbl_dif_adj')
      !
   END SUBROUTINE tra_bbl_dif_adj


   SUBROUTINE tra_bbl_adv_adj( ptb, ptb_ad, pta_ad, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trc_bbl  ***
      !!
      !! ** Purpose :   Compute the before passive tracer trend associated
      !!     with the bottom boundary layer and add it to the general trend
      !!     of tracer equations.
      !! ** Method  :   advective bbl (nn_bbl_adv = 1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean near bottom velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation i.e.
      !!                       transport proportional to the along-slope density gradient
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kjpt   ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb    ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   ptb_ad ! before and now adjoint tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_ad    ! tracer trend
      !
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      INTEGER  ::   iis , iid , ijs , ijd    ! local integers
      INTEGER  ::   ikus, ikud, ikvs, ikvd   !   -       -
      REAL(wp) ::   zbtr, ztra               ! local scalars
      REAL(wp) ::   ztraad                   !   -      -
      REAL(wp) ::   zu_bbl, zv_bbl           !   -      -
      REAL(wp) ::   zu_bblad, zv_bblad       !   -      -
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl_adv_adj')
      !
      zu_bblad = 0.0_wp ; zv_bblad = 0.0_wp
      !                                                          ! ===========
      DO jn = 1, kjpt                                            ! tracer loop
         !                                                       ! ===========
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1            ! CAUTION start from i=1 to update i=2 when cyclic east-west
               IF( vtr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero j-direction bbl advection
                  ! down-slope j/k-indices (deep)        &   up-slope j/k indices (shelf)
                  ijd  = jj + MAX( 0, mgrhv(ji,jj) )     ;   ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
                  ikvd = mbkv_d(ji,jj)                   ;   ikvs = mbkv(ji,jj)
                  zv_bbl   = ABS ( vtr_bbl(ji,jj) )
                  !                                               ! down-slope T-point (deep bottom point)
                  zbtr = r1_e1e2t(ji,ijd) / e3t(ji,ijd,ikvd)
                  ztraad = pta_ad(ji,ijd,ikvd,jn)
                  zv_bblad = zv_bblad + ztraad * ( ptb(ji,ijs,ikvs,jn) - ptb(ji,ijd,ikvd,jn) ) * zbtr
                  ptb_ad(ji,ijs,ikvs,jn) = ptb_ad(ji,ijs,ikvs,jn) + ztraad * zv_bbl * zbtr
                  ptb_ad(ji,ijd,ikvd,jn) = ptb_ad(ji,ijd,ikvd,jn) - ztraad * zv_bbl * zbtr
                  !
                  DO jk = ikvd-1, ikvs, -1                        ! down-slope upper to down T-point (deep column)
                     zbtr = r1_e1e2t(ji,ijd) / e3t(ji,ijd,jk)
                     ztraad = pta_ad(ji,ijd,jk,jn)
                     zv_bblad = zv_bblad + ztraad * ( ptb(ji,ijd,jk+1,jn) - ptb(ji,ijd,jk,jn) ) * zbtr
                     ptb_ad(ji,ijd,jk+1,jn) = ptb_ad(ji,ijd,jk+1,jn) + ztraad * zv_bbl * zbtr
                     ptb_ad(ji,ijd,jk  ,jn) = ptb_ad(ji,ijd,jk  ,jn) - ztraad * zv_bbl * zbtr
                  END DO
                  ! up  -slope T-point (shelf bottom point)
                  zbtr = r1_e1e2t(ji,ijs) / e3t(ji,ijs,ikvs)
                  ztraad = pta_ad(ji,ijs,ikvs,jn)
                  zv_bblad = zv_bblad + ztraad * ( ptb(ji,ijd,ikvs,jn) - ptb(ji,ijs,ikvs,jn) ) * zbtr
                  ptb_ad(ji,ijd,ikvs,jn) = ptb_ad(ji,ijd,ikvs,jn) + ztraad * zv_bbl * zbtr
                  ptb_ad(ji,ijs,ikvs,jn) = ptb_ad(ji,ijs,ikvs,jn) - ztraad * zv_bbl * zbtr

                  !
                  vtr_bbl_ad(ji,jj) = vtr_bbl_ad(ji,jj) + SIGN( zv_bblad, vtr_bbl(ji,jj) )
                  zv_bblad = 0.0_wp
                  !
               ENDIF



               IF( utr_bbl(ji,jj) /= 0.e0 ) THEN            ! non-zero i-direction bbl advection
                  ! down-slope i/k-indices (deep)      &   up-slope i/k indices (shelf)
                  iid  = ji + MAX( 0, mgrhu(ji,jj) )   ;   iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
                  ikud = mbku_d(ji,jj)                 ;   ikus = mbku(ji,jj)
                  zu_bbl = ABS( utr_bbl(ji,jj) )
                  !
                  zbtr   = r1_e1e2t(iid,jj) / e3t(iid,jj,ikud)
                  ztraad = pta_ad(iid,jj,ikud,jn)
                  zu_bblad = zu_bblad + ztraad * ( ptb(iis,jj,ikus,jn) - ptb(iid,jj,ikud,jn) ) * zbtr
                  ptb_ad(iis,jj,ikus,jn) = ptb_ad(iis,jj,ikus,jn) + ztraad * zu_bbl * zbtr
                  ptb_ad(iid,jj,ikud,jn) = ptb_ad(iid,jj,ikud,jn) - ztraad * zu_bbl * zbtr
                  !
                  DO jk = ikud-1, ikus, -1                            ! down-slope upper to down T-point (deep column)
                     zbtr = r1_e1e2t(iid,jj) / e3t(iid,jj,jk)
                     ztraad = pta_ad(iid,jj,jk,jn)
                     zu_bblad = zu_bblad + ztraad * ( ptb(iid,jj,jk+1,jn) - ptb(iid,jj,jk,jn) ) * zbtr
                     ptb_ad(iid,jj,jk+1,jn) = ptb_ad(iid,jj,jk+1,jn) + ztraad * zu_bbl * zbtr
                     ptb_ad(iid,jj,jk  ,jn) = ptb_ad(iid,jj,jk  ,jn) - ztraad * zu_bbl * zbtr
                  END DO
                  !                                               ! up  -slope T-point (shelf bottom point)
                  zbtr   = r1_e1e2t(iis,jj) / e3t(iis,jj,ikus)
                  ztraad = pta_ad(iis,jj,ikus,jn)
                  zu_bblad = zu_bblad + ztraad * ( ptb(iid,jj,ikus,jn) - ptb(iis,jj,ikus,jn) ) * zbtr
                  ptb_ad(iid,jj,ikus,jn) = ptb_ad(iid,jj,ikus,jn) + ztraad * zu_bbl * zbtr
                  ptb_ad(iis,jj,ikus,jn) = ptb_ad(iis,jj,ikus,jn) - ztraad * zu_bbl * zbtr
                  !
                  utr_bbl_ad(ji,jj) = utr_bbl_ad(ji,jj) + SIGN( zu_bblad, utr_bbl(ji,jj) )
                  zu_bblad = 0.0_wp
                  !
               ENDIF
               !
            END DO
            !
         END DO
         !                                                  ! ===========
      END DO                                                ! end tracer
      !                                                     ! ===========
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl_adv_adj')
      !
   END SUBROUTINE tra_bbl_adv_adj


   SUBROUTINE bbl_adj( kt, kit000, cdtype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE bbl  ***
      !!
      !! ** Purpose :   Computes the bottom boundary horizontal and vertical
      !!                advection terms.
      !!
      !! ** Method  :
      !!        * diffusive bbl (nn_bbl_ldf=1) :
      !!        When the product grad( rho) * grad(h) < 0 (where grad is an
      !!      along bottom slope gradient) an additional lateral 2nd order
      !!      diffusion along the bottom slope is added to the general
      !!      tracer trend, otherwise the additional trend is set to 0.
      !!      A typical value of ahbt is 2000 m2/s (equivalent to
      !!      a downslope velocity of 20 cm/s if the condition for slope
      !!      convection is satified)
      !!        * advective bbl (nn_bbl_adv=1 or 2) :
      !!      nn_bbl_adv = 1   use of the ocean velocity as bbl velocity
      !!      nn_bbl_adv = 2   follow Campin and Goosse (1999) implentation
      !!        i.e. transport proportional to the along-slope density gradient
      !!
      !!      NB: the along slope density gradient is evaluated using the
      !!      local density (i.e. referenced at a common local depth).
      !!
      !! References : Beckmann, A., and R. Doscher, 1997, J. Phys.Oceanogr., 581-591.
      !!              Campin, J.-M., and H. Goosse, 1999, Tellus, 412-430.
      !!----------------------------------------------------------------------
      !
      INTEGER         , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER         , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3), INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      !!
      INTEGER  ::   ji, jj                    ! dummy loop indices
      INTEGER  ::   ik                        ! local integers
      INTEGER  ::   iis , iid , ijs , ijd     !   -       -
      INTEGER  ::   ikus, ikud, ikvs, ikvd    !   -       -
      REAL(wp) ::   zsign, zsigna, zgbbl      ! local scalars
      REAL(wp) ::   zgdrho, zt, zs, zh        !   -      -
      REAL(wp) ::   zgdrhoad, ztad, zsad, zhad!   -      -
      !!
      REAL(wp) ::   fsalbt, fsbeta, pft, pfs, pfh   ! statement function
      REAL(wp) ::   pftad, pfsad, pfhad
      REAL(wp) ::   fsalbt_adj_t, fsbeta_adj_t
      REAL(wp) ::   fsalbt_adj_s, fsbeta_adj_s
      REAL(wp) ::   fsalbt_adj_h, fsbeta_adj_h
      REAL(wp), POINTER, DIMENSION(:,:) :: zub, zvb, ztb, zsb, zdep
      REAL(wp), POINTER, DIMENSION(:,:) :: zubad, zvbad, ztbad, zsbad
      !!----------------------- zv_bbl-----------------------------------------------
      ! ratio alpha/beta = fsalbt : ratio of thermal over saline expension coefficients
      ! ================            pft :  potential temperature in degrees celcius
      !                             pfs :  salinity anomaly (s-35) in psu
      !                             pfh :  depth in meters
      ! nn_eos = 0  (Jackett and McDougall 1994 formulation)
      fsalbt( pft, pfs, pfh ) =                                              &   ! alpha/beta
         ( ( ( -0.255019e-07 * pft + 0.298357e-05 ) * pft                    &
                                   - 0.203814e-03 ) * pft                    &
                                   + 0.170907e-01 ) * pft                    &
                                   + 0.665157e-01                            &
         +(-0.678662e-05 * pfs - 0.846960e-04 * pft + 0.378110e-02 ) * pfs   &
         +  ( ( - 0.302285e-13 * pfh                                         &
                - 0.251520e-11 * pfs                                         &
                + 0.512857e-12 * pft * pft          ) * pfh                  &
                                     - 0.164759e-06   * pfs                  &
             +(   0.791325e-08 * pft - 0.933746e-06 ) * pft                  &
                                     + 0.380374e-04 ) * pfh
      fsbeta( pft, pfs, pfh ) =                                              &   ! beta
         ( ( -0.415613e-09 * pft + 0.555579e-07 ) * pft                      &
                                 - 0.301985e-05 ) * pft                      &
                                 + 0.785567e-03                              &
         + (     0.515032e-08 * pfs                                          &
               + 0.788212e-08 * pft - 0.356603e-06 ) * pfs                   &
               +(  (   0.121551e-17 * pfh                                    &
                     - 0.602281e-15 * pfs                                    &
                     - 0.175379e-14 * pft + 0.176621e-12 ) * pfh             &
                                          + 0.408195e-10   * pfs             &
                 + ( - 0.213127e-11 * pft + 0.192867e-09 ) * pft             &
                                          - 0.121555e-07 ) * pfh

      fsalbt_adj_t( pft, pfs, pfh, pftad ) =  &   ! alpha/beta
         &         ( - 0.255019e-07 * 4 * pft * pft * pft            &
         &           + 0.298357e-05 * 3 * pft * pft                  &
         &           - 0.203814e-03 * 2 * pft                        &
         &           - 0.846960e-04 * pfs                            &
         &           + 0.512857e-12 * 2 * pft * pfh * pfh            &
         &           + 0.791325e-08 * pft * pfh                      &
         &           - 0.933746e-06 * pfh                            &
         &           + 0.170907e-01                        ) * pftad

      fsalbt_adj_s( pft, pfs, pfh, pfsad ) =  &   ! alpha/beta
         &       + ( - 0.678662e-05 * 2 * pfs                        &
         &           - 0.846960e-04 * pft                            &
         &           - 0.251520e-11 * pfh  * pfh                     &
         &           - 0.164759e-06 * pfh                            &
         &           + 0.378110e-02                        ) * pfsad

      fsalbt_adj_h( pft, pfs, pfh, pfhad ) =  &   ! alpha/beta
         &       + ( - 0.302285e-13 * 3 * pfh  * pfh                 &
         &           - 0.251520e-11 * pfs * pfh                      &
         &           + 0.512857e-12 * pft * pft * pfh                &
         &           - 0.164759e-06 * pfs                            &
         &           + 0.791325e-08 * pft * pft                      &
         &           - 0.933746e-06 * pft                            &
         &           + 0.380374e-04                        ) * pfhad


      fsbeta_adj_t( pft, pfs, pfh, pftad ) =  &   ! beta
         &         ( - 0.415613e-09 * 3 * pft * pft                  &
         &           + 0.555579e-07 * 2 * pft                        &
         &           - 0.301985e-05                                  &
         &           + 0.788212e-08 * pfs                            &
         &           - 0.213127e-11 * 2 * pfh * pft                  &
         &           - 0.175379e-14 * pfh * pfh            ) * pftad
      fsbeta_adj_s( pft, pfs, pfh, pfsad ) =  &   ! beta
         &         (   0.788212e-08 * pft                            &
         &           + 0.515032e-08 * 2 * pfs                        &
         &           - 0.356603e-06                                  &
         &           + 0.408195e-10 * pfh                            &
         &           - 0.602281e-15 * pfh * pfh            ) * pfsad
      fsbeta_adj_h( pft, pfs, pfh, pfhad ) =  &   ! beta
         &         (   0.121551e-17 * 3 * pfh * pfh                  &
         &           - 0.602281e-15 * 2 * pfs * pfh                  &
         &           - 0.175379e-14 * 2 * pft * pfh                  &
         &           + 0.176621e-12 * 2 * pfh                        &
         &           + 0.408195e-10 * pfs                            &
         &           + 0.192867e-09 * pfh                            &
         &           - 0.213127e-11 * pft * pft                      &
         &           + 0.192867e-09 * pft                            &
         &           - 0.121555e-07                        ) * pfhad
      !!----------------------------------------------------------------------

      !
      IF( nn_timing == 1 )  CALL timing_start( 'bbl_adj')
      !
      CALL wrk_alloc( jpi, jpj, zub  , zvb  , ztb  , zsb, zdep, &
         &                      zubad, zvbad, ztbad, zsbad      )
      !
      zubad(:,:) = 0.0_wp ; zvbad(:,:) = 0.0_wp ; ztbad(:,:) = 0.0_wp ; zsbad(:,:) = 0.0_wp

      IF( kt == kit000 )  THEN
         IF(lwp)  WRITE(numout,*)
         IF(lwp)  WRITE(numout,*) 'trabbl_tam:bbl_adj : Compute bbl velocities and diffusive coefficients in ', cdtype
         IF(lwp)  WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !                                        !* bottom temperature, salinity, velocity and depth
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = mbkt(ji,jj)                        ! bottom T-level index
            ztb (ji,jj) = tsb(ji,jj,ik,jp_tem) * tmask(ji,jj,1)      ! bottom before T and S
            zsb (ji,jj) = tsb(ji,jj,ik,jp_sal) * tmask(ji,jj,1)
            zdep(ji,jj) = gdept(ji,jj,ik)        ! bottom T-level reference depth
            !
            zub(ji,jj) = un(ji,jj,mbku(ji,jj))      ! bottom velocity
            zvb(ji,jj) = vn(ji,jj,mbkv(ji,jj))
         END DO
      END DO
      !                                   !-------------------!
      IF( nn_bbl_adv /= 0 ) THEN          !   advective bbl   !
         !                                !-------------------!
         SELECT CASE ( nn_bbl_adv )             !* bbl transport type
               !
         CASE( 1 )                                   != use of upper velocity
            ! NOTE: not much needed for deriving, almost all the computations are for the SIGN, which is kept as in the NL
            DO jj = 1, jpjm1                                 ! criteria: grad(rho).grad(h)<0  and grad(rho).grad(h)<0
               DO ji = 1, jpim1   ! vector opt.
                  !                                                ! j-direction
                  zt = 0.5 * ( ztb (ji,jj+1) + ztb (ji,jj) )                ! T, S anomalie, and depth
                  zs = 0.5 * ( zsb (ji,jj+1) + zsb (ji,jj) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj+1) + zdep(ji,jj) )
                  !                                                         ! masked bbl j-gradient of density
                  zgdrho = (  fsalbt( zt, zs, zh ) * ( ztb(ji,jj+1) - ztb(ji,jj) )    &
                  &                             - ( zsb(ji,jj+1) - zsb(ji,jj) )  ) * vmask(ji,jj,1)
                  zsign = SIGN(  0.5, - zgdrho   * REAL( mgrhv(ji,jj) )  )    ! sign of j-gradient * j-slope
                  zsigna= SIGN(  0.5, zvb(ji,jj) * REAL( mgrhv(ji,jj) )  )    ! sign of u * i-slope
                  !
                  !                                                           ! bbl velocity
                  zvbad(ji,jj) = zvbad(ji,jj) + vtr_bbl_ad(ji,jj) * ( 0.5 + zsigna ) * ( 0.5 - zsign )   &
                               &                                  *     e1v(ji,jj) * e3v_bbl_0(ji,jj)
                  vtr_bbl_ad(ji,jj) = 0.0_wp
                  !                                               ! i-direction
                  zt = 0.5 * ( ztb (ji,jj) + ztb (ji+1,jj) )                  ! T, S anomalie, and depth
                  zs = 0.5 * ( zsb (ji,jj) + zsb (ji+1,jj) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )
                  !                                                           ! masked bbl i-gradient of density
                  zgdrho = (  fsalbt( zt, zs, zh ) * ( ztb(ji+1,jj) - ztb(ji,jj) )    &
                     &                             - ( zsb(ji+1,jj) - zsb(ji,jj) )  ) * umask(ji,jj,1)
                  !
                  zsign = SIGN(  0.5, - zgdrho   * REAL( mgrhu(ji,jj) )  )    ! sign of i-gradient * i-slope
                  zsigna= SIGN(  0.5, zub(ji,jj) * REAL( mgrhu(ji,jj) )  )    ! sign of u * i-slope
                  !
                  !                                                           ! bbl velocity
                  zubad(ji,jj) = zubad(ji,jj) + utr_bbl_ad(ji,jj) * ( 0.5 + zsigna ) * ( 0.5 - zsign )   &
                               &                                  * e2u(ji,jj) * e3u_bbl_0(ji,jj)
                  utr_bbl_ad(ji,jj) = 0.0_wp
               !
            END DO
         END DO
            !
         CASE( 2 )                                 != bbl velocity = F( delta rho )
            ! NOTE: this one is nastier
            zgbbl = grav * rn_gambbl
            DO jj = jpjm1, 1, -1                            ! criteria: rho_up > rho_down
               DO ji = jpim1, 1, -1   ! vector opt.
                  !                                               ! j-direction
                  !  down-slope T-point j/k-index (deep)  &   of the up  -slope T-point j/k-index (shelf)
                  ijd  = jj + MAX( 0, mgrhv(ji,jj) )      ;    ijs  = jj + 1 - MAX( 0, mgrhv(ji,jj) )
                  ikvd = mbkv_d(ji,jj)                    ;    ikvs = mbkv(ji,jj)
                  !
                  !                                             ! mid-depth density anomalie (up-slope minus down-slope)
                  zt = 0.5 * ( ztb (ji,jj) + ztb (ji,jj+1) )           ! mid slope depth of T, S, and depth
                  zs = 0.5 * ( zsb (ji,jj) + zsb (ji,jj+1) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj) + zdep(ji,jj+1) )
                  zgdrho =    fsbeta( zt, zs, zh )                                    &
                     &   * (  fsalbt( zt, zs, zh ) * ( ztb(ji,ijd) - ztb(ji,ijs) )    &
                     &                             - ( zsb(ji,ijd) - zsb(ji,ijs) )  ) * vmask(ji,jj,1)
                  !
                  zsign  = SIGN( 0.5_wp, zgdrho ) ! adjoint of zgdrho = MAX( 0.e0, zgdrho )
                  !                                             ! bbl transport (down-slope direction)
                  zgdrhoad = zsign * e1v(ji,jj) * e3v_bbl_0(ji,jj) * zgbbl * vtr_bbl_ad(ji,jj) * REAL( mgrhv(ji,jj) )
                  vtr_bbl_ad(ji,jj) = 0.0_wp

                  ztad = ( fsbeta_adj_t( zt, zs, zh, zgdrhoad )                        &
                     & * ( fsalbt( zt, zs, zh ) * ( ztb(ji,ijd) - ztb(ji,ijs) )        &
                     &                          - ( zsb(ji,ijd) - zsb(ji,ijs) ) )      &
                     & +   fsbeta( zt, zs, zh ) * fsalbt_adj_t( zt, zs, zh, zgdrhoad ) &
                     &                          * ( ztb(ji,ijd) - ztb(ji,ijs) )        &
                     &   ) * vmask(ji,jj,1)
                  zsad = ( fsbeta_adj_s( zt, zs, zh, zgdrhoad )                        &
                     & * ( fsalbt( zt, zs, zh ) * ( ztb(ji,ijd) - ztb(ji,ijs) )        &
                     &                          - ( zsb(ji,ijd) - zsb(ji,ijs) ) )      &
                     & +   fsbeta( zt, zs, zh ) * fsalbt_adj_s( zt, zs, zh, zgdrhoad ) &
                     &                          * ( ztb(ji,ijd) - ztb(ji,ijs) )        &
                     &   ) * vmask(ji,jj,1)
                  zhad = ( fsbeta_adj_h( zt, zs, zh, zgdrhoad )                        &
                     & * ( fsalbt( zt, zs, zh ) * ( ztb(ji,ijd) - ztb(ji,ijs) )        &
                     &                          - ( zsb(ji,ijd) - zsb(ji,ijs) ) )      &
                     & +   fsbeta( zt, zs, zh ) * fsalbt_adj_h( zt, zs, zh, zgdrhoad ) &
                     &                          * ( ztb(ji,ijd) - ztb(ji,ijs) )        &
                     &   ) * vmask(ji,jj,1)

                  ztbad(ji,ijd) = ztbad(ji,ijd) + zgdrhoad * fsbeta( zt, zs, zh ) * fsalbt( zt, zs, zh ) * vmask(ji,jj,1)
                  ztbad(ji,ijs) = ztbad(ji,ijs) - zgdrhoad * fsbeta( zt, zs, zh ) * fsalbt( zt, zs, zh ) * vmask(ji,jj,1)
                  zsbad(ji,ijd) = zsbad(ji,ijd) - zgdrhoad * fsbeta( zt, zs, zh ) * vmask(ji,jj,1)
                  zsbad(ji,ijs) = zsbad(ji,ijs) + zgdrhoad * fsbeta( zt, zs, zh ) * vmask(ji,jj,1)

                  ztbad (ji,jj  ) = ztbad (ji,jj  ) + 0.5 * ztad
                  ztbad (ji,jj+1) = ztbad (ji,jj+1) + 0.5 * ztad
                  zsbad (ji,jj  ) = zsbad (ji,jj  ) + 0.5 * zsad
                  zsbad (ji,jj+1) = zsbad (ji,jj+1) + 0.5 * zsad
                  ztad = 0.0_wp ; zsad = 0.0_wp ; zhad = 0.0_wp

                  !                                         ! i-direction
                  ! down-slope T-point i/k-index (deep)  &   up-slope T-point i/k-index (shelf)
                  iid  = ji + MAX( 0, mgrhu(ji,jj) )     ;    iis  = ji + 1 - MAX( 0, mgrhu(ji,jj) )
                  ikud = mbku_d(ji,jj)                   ;    ikus = mbku(ji,jj)
                  !
                  !                                             ! mid-depth density anomalie (up-slope minus down-slope)
                  zt = 0.5 * ( ztb (ji,jj) + ztb (ji+1,jj) )           ! mid slope depth of T, S, and depth
                  zs = 0.5 * ( zsb (ji,jj) + zsb (ji+1,jj) ) - 35.0
                  zh = 0.5 * ( zdep(ji,jj) + zdep(ji+1,jj) )
                  zgdrho =    fsbeta( zt, zs, zh )                                    &
                     &   * (  fsalbt( zt, zs, zh ) * ( ztb(iid,jj) - ztb(iis,jj) )    &
                     &                             - ( zsb(iid,jj) - zsb(iis,jj) )  ) * umask(ji,jj,1)
                  zsign  = SIGN( 0.5_wp, zgdrho ) ! adjoint of zgdrho = MAX( 0.e0, zgdrho )
                  !                                             ! bbl transport (down-slope direction)
                  zgdrhoad = zsign * e2u(ji,jj) * e3u_bbl_0(ji,jj) * zgbbl * utr_bbl_ad(ji,jj) * REAL( mgrhu(ji,jj) )
                  utr_bbl_ad(ji,jj) = 0.0_wp
                  !
                  ztad = ( fsbeta_adj_t( zt, zs, zh, zgdrhoad )                        &
                     & * ( fsalbt( zt, zs, zh ) * ( ztb(iid,jj) - ztb(iis, jj) )       &
                     &                          - ( zsb(iid,jj) - zsb(iis, jj) ) )     &
                     & +   fsbeta( zt, zs, zh ) * fsalbt_adj_t( zt, zs, zh, zgdrhoad ) &
                     &                          * ( ztb(iid,jj) - ztb(iis, jj) )       &
                     &   ) * umask(ji,jj,1)
                  zsad = ( fsbeta_adj_s( zt, zs, zh, zgdrhoad )                        &
                     & * ( fsalbt( zt, zs, zh ) * ( ztb(iid,jj) - ztb(iis, jj) )       &
                     &                          - ( zsb(iid,jj) - zsb(iis, jj) ) )     &
                     & +   fsbeta( zt, zs, zh ) * fsalbt_adj_s( zt, zs, zh, zgdrhoad ) &
                     &                          * ( ztb(iid,jj) - ztb(iis, jj) )       &
                     &   ) * umask(ji,jj,1)
                  zhad = ( fsbeta_adj_h( zt, zs, zh, zgdrhoad )                        &
                     & * ( fsalbt( zt, zs, zh ) * ( ztb(iid,jj) - ztb(iis, jj) )       &
                     &                          - ( zsb(iid,jj) - zsb(iis, jj) ) )     &
                     & +   fsbeta( zt, zs, zh ) * fsalbt_adj_h( zt, zs, zh, zgdrhoad ) &
                     &                          * ( ztb(iid,jj) - ztb(iis, jj) )       &
                     &   ) * umask(ji,jj,1)

                  ztbad(iid,jj) = ztbad(iid,jj) + zgdrhoad * fsbeta( zt, zs, zh ) * fsalbt( zt, zs, zh ) * umask(ji,jj,1)
                  ztbad(iis,jj) = ztbad(iis,jj) - zgdrhoad * fsbeta( zt, zs, zh ) * fsalbt( zt, zs, zh ) * umask(ji,jj,1)
                  zsbad(iid,jj) = zsbad(iid,jj) - zgdrhoad * fsbeta( zt, zs, zh ) * umask(ji,jj,1)
                  zsbad(iis,jj) = zsbad(iis,jj) + zgdrhoad * fsbeta( zt, zs, zh ) * umask(ji,jj,1)
                  zgdrhoad = 0.0_wp

                  ztbad (ji,jj  ) = ztbad (ji,jj  ) + 0.5 * ztad
                  ztbad (ji+1,jj) = ztbad (ji+1,jj) + 0.5 * ztad
                  zsbad (ji,jj  ) = zsbad (ji,jj  ) + 0.5 * zsad
                  zsbad (ji+1,jj) = zsbad (ji+1,jj) + 0.5 * zsad
                  ztad = 0.0_wp ; zsad = 0.0_wp ; zhad = 0.0_wp
                  !
               END DO
            END DO
         END SELECT
         !
      ENDIF
      IF( nn_bbl_ldf == 1 ) THEN          !   diffusive bbl   !
         !                                !-------------------!
         ! NOTE : while rn_ahtbbl remains a passive variable, the code below will only yield ah_bbl_ad=0
         DO jj = 1, jpjm1
            DO ji = 1, jpim1
               ahu_bbl_ad(ji,jj)=0.0_wp
               ahv_bbl_ad(ji,jj)=0.0_wp
            END DO
         END DO
            !
      ENDIF
      !                                        !* bottom temperature, salinity, velocity and depth
      DO jj = 1, jpj
         DO ji = 1, jpi
            ik = mbkt(ji,jj)                        ! bottom T-level index
            tsb_ad(ji,jj,ik,jp_tem) = tsb_ad(ji,jj,ik,jp_tem) + ztbad(ji,jj) * tmask(ji,jj,1)
            tsb_ad(ji,jj,ik,jp_sal) = tsb_ad(ji,jj,ik,jp_sal) + zsbad(ji,jj) * tmask(ji,jj,1)
            ztbad (ji,jj) = 0.0_wp
            zsbad (ji,jj) = 0.0_wp

            un_ad(ji,jj,mbku(ji,jj)) = un_ad(ji,jj,mbku(ji,jj)) + zubad(ji,jj)
            vn_ad(ji,jj,mbkv(ji,jj)) = vn_ad(ji,jj,mbkv(ji,jj)) + zvbad(ji,jj)
            zvbad(ji,jj) = 0.0_wp
            zubad(ji,jj) = 0.0_wp
         END DO
      END DO
      !                                   !-------------------!
      !
      CALL wrk_dealloc( jpi, jpj, zub  , zvb  , ztb  , zsb, zdep, &
         &                        zubad, zvbad, ztbad, zsbad      )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'bbl_adj')
      !
   END SUBROUTINE bbl_adj


   SUBROUTINE tra_bbl_init_tam
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbl_init  ***
      !!
      !! ** Purpose :   Initialization for the bottom boundary layer scheme.
      !!
      !! ** Method  :
      !!----------------------------------------------------------------------
      !
      integer :: ierr
      IF( nn_timing == 1 )  CALL timing_start( 'tra_bbl_init_tam')
      !
      ierr = tra_bbl_alloc_tam( 0 )

      ahu_bbl_0_tl = 0.0_wp
      ahv_bbl_0_tl = 0.0_wp
      ahu_bbl_0_ad = 0.0_wp
      ahv_bbl_0_ad = 0.0_wp
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_bbl_init_tam')
      !
   END SUBROUTINE tra_bbl_init_tam

   SUBROUTINE bbl_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE tra_bbl_adj_tst ***
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
      !! * Local declarations
      INTEGER ::  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & jtst
      INTEGER ::  &
         & jsav1, &
         & jsav2
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), POINTER, DIMENSION(:,:,:,:) :: &
         & ztsb_tlin ,     & ! Tangent input
         & ztsb_adout,     & ! Adjoint input
         & zrts            ! 2*3D random field
      REAL(KIND=wp), POINTER, DIMENSION(:,:,:) :: &
         & zun_tlin,       &
         & zvn_tlin,       &
         & zun_adout,      &
         & zvn_adout,      &
         & zr3             ! 3D random field
      REAL(KIND=wp), POINTER, DIMENSION(:,:) :: &
         & zahu_tlout,      &
         & zahv_tlout,      &
         & zahu_adin ,      &
         & zahv_adin ,      &
         & zutr_tlout,      &
         & zvtr_tlout,      &
         & zutr_adin,       &
         & zvtr_adin,       &
         & zr2
      CHARACTER(LEN=14) :: &
         & cl_name
      ! Allocate memory
      CALL wrk_alloc( jpi, jpj, jpk, jpts, ztsb_tlin , ztsb_adout, zrts )
      CALL wrk_alloc( jpi, jpj, jpk, zun_tlin, zvn_tlin, zun_adout, zvn_adout, zr3 )
      CALL wrk_alloc( jpi, jpj, zahu_tlout, zahv_tlout, zahu_adin, zahv_adin, &
         &                      zutr_tlout, zvtr_tlout, zutr_adin, zvtr_adin, zr2 )

      CALL grid_random( utr_bbl(:,:), 'U', 0.0_wp, stdu )
      CALL grid_random( vtr_bbl(:,:), 'V', 0.0_wp, stdv )

      jsav1 = nn_bbl_ldf
      jsav2 = nn_bbl_adv

      DO jtst = 1, 3
         !==================================================================
         ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
         !    dy = ( hdivb_tl, hdivn_tl )
         !==================================================================

         SELECT CASE( jtst)
         CASE ( 1 )
            nn_bbl_ldf = 1
            nn_bbl_adv = 0
         CASE ( 2 )
            nn_bbl_ldf = 0
            nn_bbl_adv = 1
         CASE ( 3 )
            nn_bbl_ldf = 0
            nn_bbl_adv = 2
         END SELECT
         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------
         ztsb_tlin = 0.0_wp
         zun_tlin  = 0.0_wp
         zvn_tlin  = 0.0_wp
         zahu_adin = 0.0_wp
         zahv_adin = 0.0_wp
         zutr_adin = 0.0_wp
         zvtr_adin = 0.0_wp

         ahu_bbl_tl = 0.0_wp
         ahv_bbl_tl = 0.0_wp
         utr_bbl_tl = 0.0_wp
         vtr_bbl_tl = 0.0_wp

         un_ad  = 0.0_wp
         vn_ad  = 0.0_wp
         tsb_ad = 0.0_wp

         !--------------------------------------------------------------------
         ! Initialize the tangent input with random noise: dx
         !--------------------------------------------------------------------

         CALL grid_random( zrts(:,:,:,jp_tem), 'T', 0.0_wp, stdt )
         CALL grid_random( zrts(:,:,:,jp_sal), 'T', 0.0_wp, stds )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  ztsb_tlin(ji,jj,jk,:) = zrts(ji,jj,jk,:)
               END DO
            END DO
         END DO

         CALL grid_random( zr3(:,:,:), 'U', 0.0_wp, stdu )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zun_tlin(ji,jj,jk) = zr3(ji,jj,jk)
               END DO
            END DO
         END DO

         CALL grid_random( zr3(:,:,:), 'V', 0.0_wp, stdv )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zvn_tlin(ji,jj,jk) = zr3(ji,jj,jk)
               END DO
            END DO
         END DO

         tsb_tl(:,:,:,:) = ztsb_tlin(:,:,:,:)
         un_tl(:,:,:)    = zun_tlin(:,:,:)
         vn_tl(:,:,:)    = zvn_tlin(:,:,:)

         CALL bbl_tan (0, 1, 'TRA')

         zahu_tlout(:,:) = ahu_bbl_tl(:,:)
         zahv_tlout(:,:) = ahv_bbl_tl(:,:)
         zutr_tlout(:,:) = utr_bbl_tl(:,:)
         zvtr_tlout(:,:) = vtr_bbl_tl(:,:)
         
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zahu_adin(ji,jj) = zahu_tlout(ji,jj) &
                  &             * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) &
                  &             * umask(ji,jj,1)
               zahv_adin(ji,jj) = zahv_tlout(ji,jj) &
                  &             * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,1) &
                  &             * vmask(ji,jj,1)
               zutr_adin(ji,jj) = zutr_tlout(ji,jj) &
                  &             * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) &
                  &             * umask(ji,jj,1)
               zvtr_adin(ji,jj) = zvtr_tlout(ji,jj) &
                  &             * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,1) &
                  &             * vmask(ji,jj,1)
            END DO
         END DO
         !--------------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !--------------------------------------------------------------------

         zsp1 = DOT_PRODUCT( zahu_tlout, zahu_adin ) &
            & + DOT_PRODUCT( zutr_tlout, zutr_adin ) &
            & + DOT_PRODUCT( zahv_tlout, zahv_adin ) &
            & + DOT_PRODUCT( zvtr_tlout, zvtr_adin )

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------

         ahu_bbl_ad(:,:) = zahu_adin(:,:)
         ahv_bbl_ad(:,:) = zahv_adin(:,:)
         utr_bbl_ad(:,:) = zutr_adin(:,:)
         vtr_bbl_ad(:,:) = zvtr_adin(:,:)

         CALL bbl_adj (0, 1, 'TRA')

         ztsb_adout = tsb_ad
         zun_adout  = un_ad
         zvn_adout  = vn_ad

         zsp2 = DOT_PRODUCT( ztsb_tlin(:,:,:,jp_tem), ztsb_adout(:,:,:,jp_tem) ) &
            & + DOT_PRODUCT( ztsb_tlin(:,:,:,jp_sal), ztsb_adout(:,:,:,jp_sal) ) &
            & + DOT_PRODUCT( zun_tlin (:,:,:       ), zun_adout (:,:,:       ) ) &
            & + DOT_PRODUCT( zvn_tlin (:,:,:       ), zvn_adout (:,:,:       ) ) 

         SELECT CASE ( jtst )
         CASE ( 1 )
            ! 14 char:'12345678901234'
            cl_name = 'bbl_adj_dif   '
         CASE ( 2 )
            ! 14 char:'12345678901234'
            cl_name = 'bbl_adj_adv  1'
         CASE ( 3 )
            ! 14 char:'12345678901234'
            cl_name = 'bbl_adj_adv  2'
         END SELECT
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      END DO

      nn_bbl_ldf = jsav1
      nn_bbl_adv = jsav2


      CALL wrk_dealloc( jpi, jpj, jpk, jpts, ztsb_tlin , ztsb_adout, zrts )
      CALL wrk_dealloc( jpi, jpj, jpk, zun_tlin, zvn_tlin, zun_adout, zvn_adout, zr3 )
      CALL wrk_dealloc( jpi, jpj, zahu_tlout, zahv_tlout, zahu_adin, zahv_adin, &
         &                      zutr_tlout, zvtr_tlout, zutr_adin, zvtr_adin, zr2 )

   END SUBROUTINE bbl_adj_tst

   SUBROUTINE tra_bbl_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE tra_bbl_adj_tst ***
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

      !! * Local declarations
      INTEGER ::  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & jtst
      INTEGER ::  &
         & jsav1, &
         & jsav2
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), POINTER, DIMENSION(:,:,:,:) :: &
         & ztsa_tlin ,     & ! Tangent input
         & ztsa_tlout,     &
         & ztsb_tlin ,     &
         & ztsa_adout,     & ! Adjoint input
         & ztsa_adin ,     &
         & ztsb_adout,     &
         & zrts            ! 2*3D random field
      REAL(KIND=wp), POINTER, DIMENSION(:,:,:) :: &
         & zun_tlin,       &
         & zvn_tlin,       &
         & zun_adout,      &
         & zvn_adout,      &
         & zr3             ! 3D random field
      CHARACTER(LEN=14) :: &
         & cl_name
      ! Allocate memory

      CALL wrk_alloc( jpi, jpj, jpk, jpts, ztsa_tlin , ztsa_tlout, ztsb_tlin , &
         &                                  ztsa_adout, ztsa_adin , ztsb_adout, &
         &                                  zrts )
      CALL wrk_alloc( jpi, jpj, jpk, zun_tlin, zvn_tlin, zun_adout, zvn_adout, zr3 )

      CALL grid_random( utr_bbl(:,:), 'U', 0.0_wp, stdu )
      CALL grid_random( vtr_bbl(:,:), 'V', 0.0_wp, stdv )

      jsav1 = nn_bbl_ldf
      jsav2 = nn_bbl_adv

      DO jtst = 1, 3
         !==================================================================
         ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
         !    dy = ( hdivb_tl, hdivn_tl )
         !==================================================================

         SELECT CASE( jtst)
         CASE ( 1 )
            nn_bbl_ldf = 1
            nn_bbl_adv = 0
         CASE ( 2 )
            nn_bbl_ldf = 0
            nn_bbl_adv = 1
         CASE ( 3 )
            nn_bbl_ldf = 0
            nn_bbl_adv = 2
         END SELECT
         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------
         ztsa_tlin (:,:,:,:) = 0.0_wp
         ztsa_tlout(:,:,:,:) = 0.0_wp
         ztsb_tlin (:,:,:,:) = 0.0_wp
         ztsa_adout(:,:,:,:) = 0.0_wp
         ztsa_adin (:,:,:,:) = 0.0_wp
         ztsb_adout(:,:,:,:) = 0.0_wp

         zun_tlin (:,:,:) = 0.0_wp
         zun_adout(:,:,:) = 0.0_wp
         zvn_tlin (:,:,:) = 0.0_wp
         zvn_adout(:,:,:) = 0.0_wp

         tsb_tl(:,:,:,:) = 0.0_wp
         tsa_tl(:,:,:,:) = 0.0_wp
         tsb_ad(:,:,:,:) = 0.0_wp
         tsa_ad(:,:,:,:) = 0.0_wp

         un_tl(:,:,:) = 0.0_wp
         vn_tl(:,:,:) = 0.0_wp
         un_ad(:,:,:) = 0.0_wp
         vn_ad(:,:,:) = 0.0_wp
         !--------------------------------------------------------------------
         ! Initialize the tangent input with random noise: dx
         !--------------------------------------------------------------------

         CALL grid_random( zrts(:,:,:,jp_tem), 'T', 0.0_wp, stdt )
         CALL grid_random( zrts(:,:,:,jp_sal), 'T', 0.0_wp, stds )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  ztsa_tlin(ji,jj,jk,:) = zrts(ji,jj,jk,:)
               END DO
            END DO
         END DO

         CALL grid_random( zrts(:,:,:,jp_tem), 'T', 0.0_wp, stdt )
         CALL grid_random( zrts(:,:,:,jp_sal), 'T', 0.0_wp, stds )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  ztsb_tlin(ji,jj,jk,:) = zrts(ji,jj,jk,:)
               END DO
            END DO
         END DO

         CALL grid_random( zr3(:,:,:), 'U', 0.0_wp, stdu )
         DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
                  zun_tlin(ji,jj,jk) = zr3(ji,jj,jk)
            END DO
         END DO
         END DO

         CALL grid_random( zr3(:,:,:), 'V', 0.0_wp, stdv )
         DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
                  zvn_tlin(ji,jj,jk) = zr3(ji,jj,jk)
            END DO
         END DO
         END DO

         tsa_tl(:,:,:,:) = ztsa_tlin(:,:,:,:)
         tsb_tl(:,:,:,:) = ztsb_tlin(:,:,:,:)
         un_tl(:,:,:)    = zun_tlin(:,:,:)
         vn_tl(:,:,:)    = zvn_tlin(:,:,:)

         CALL tra_bbl_tan ( nit000 )
         ztsa_tlout(:,:,:,:) = tsa_tl(:,:,:,:)
         !--------------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !--------------------------------------------------------------------

         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  ztsa_adin(ji,jj,jk,jp_tem) = ztsa_tlout(ji,jj,jk,jp_tem) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                     &               * tmask(ji,jj,jk)
                  ztsa_adin(ji,jj,jk,jp_sal) = ztsa_tlout(ji,jj,jk,jp_sal) &
                     &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                     &               * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !--------------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !--------------------------------------------------------------------

         zsp1 = DOT_PRODUCT( ztsa_tlout(:,:,:,jp_tem), ztsa_adin(:,:,:,jp_tem) ) &
            & + DOT_PRODUCT( ztsa_tlout(:,:,:,jp_sal), ztsa_adin(:,:,:,jp_sal) ) 

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------

         tsa_ad(:,:,:,:) = ztsa_adin(:,:,:,:)
         CALL tra_bbl_adj ( nit000 )
         ztsa_adout(:,:,:,:) = tsa_ad(:,:,:,:)
         ztsb_adout(:,:,:,:) = tsb_ad(:,:,:,:)
         zun_adout(:,:,:)    = un_ad (:,:,:)
         zvn_adout(:,:,:)    = vn_ad (:,:,:)

         zsp2 = DOT_PRODUCT( ztsa_tlin(:,:,:,jp_tem), ztsa_adout(:,:,:,jp_tem) ) &
            & + DOT_PRODUCT( ztsa_tlin(:,:,:,jp_sal), ztsa_adout(:,:,:,jp_sal) ) &
            & + DOT_PRODUCT( ztsb_tlin(:,:,:,jp_tem), ztsb_adout(:,:,:,jp_tem) ) &
            & + DOT_PRODUCT( ztsb_tlin(:,:,:,jp_sal), ztsb_adout(:,:,:,jp_sal) ) &
            & + DOT_PRODUCT( zun_tlin (:,:,:       ), zun_adout (:,:,:       ) ) &
            & + DOT_PRODUCT( zvn_tlin (:,:,:       ), zvn_adout (:,:,:       ) ) 

         SELECT CASE ( jtst )
         CASE ( 1 )
            ! 14 char:'12345678901234'
            cl_name = 'trabbl_adj_dif'
         CASE ( 2 )
            ! 14 char:'12345678901234'
            cl_name = 'trabbl_ad_adv1'
         CASE ( 3 )
            ! 14 char:'12345678901234'
            cl_name = 'trabbl_ad_adv2'
         END SELECT
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      END DO

      CALL wrk_dealloc( jpi, jpj, jpk, jpts, ztsa_tlin , ztsa_tlout, ztsb_tlin , &
         &                                    ztsa_adout, ztsa_adin , ztsb_adout, &
         &                                    zrts )
      CALL wrk_dealloc( jpi, jpj, jpk, zun_tlin, zvn_tlin, zun_adout, zvn_adout, zr3 )

      nn_bbl_ldf = jsav1
      nn_bbl_adv = jsav2


   END SUBROUTINE tra_bbl_adj_tst

   !!======================================================================
END MODULE trabbl_tam
