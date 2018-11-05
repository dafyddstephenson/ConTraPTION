MODULE sbcssr_tam
   !!======================================================================
   !!                       ***  MODULE  sbcssr_tam  ***
   !! Surface module :  add to heat and fresh water fluxes a restoring term
   !!                   toward observed SST/SSS
   !!                   Tangent and Adjoint Module
   !!======================================================================
   !! History of the direct routine:
   !!            3.0  !  2006-06  (G. Madec)   Original code
   !!            3.2  !  2009-04  (B. Lemaire) Introduce iom_put
   !! History of the T&A routine:
   !!            3.0  !  2008-11  (A. Vidard)  Original code (simplification: no linear salinity damping)
   !!            3.2  !  2010-04  (A. Vidard)  Nemo3.2 update
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_ssr        : add to sbc a restoring term toward SST/SSS climatology
   !!----------------------------------------------------------------------
   USE par_oce
   USE par_kind
   USE dom_oce
   USE sbc_oce
   USE sbc_oce_tam
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE timing
   USE prtctl
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_ssr_tan    ! routine called in sbcmod_tam
   PUBLIC   sbc_ssr_adj    ! routine called in sbcmod_tam
   PUBLIC   sbc_ssr_ini_tam    ! routine called in sbcmod_tam
   PUBLIC   sbc_ssr_adj_tst! routine called in tst

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:) ::  &
     & erp_tl,    &    !: evaporation damping                          [kg/m2/s]
     & qrp_tl,    &    !: heat flux damping                            [w/m2]
     & erp_ad,    &    !: evaporation damping                          [kg/m2/s]
     & qrp_ad          !: heat flux damping                            [w/m2]

   !                                           !!* Namelist namsbc_ssr *
   INTEGER, PUBLIC ::   nn_sstr     =   0       ! SST/SSS restoring indicator
   INTEGER, PUBLIC ::   nn_sssr     =   0       ! SST/SSS restoring indicator
   REAL(wp)        ::   rn_dqdt     = -40.e0    ! restoring factor on SST and SSS
   REAL(wp)        ::   rn_deds     = -27.70    ! restoring factor on SST and SSS
   LOGICAL         ::   ln_sssr_bnd = .false.   ! flag to bound erp term
   REAL(wp)        ::   rn_sssr_bnd =   0.e0    ! ABS(Max./Min.) value of erp term [mm/day]

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
   !! NEMO/OPA 3.2 , LOCEAN-IPSL (2009)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE sbc_ssr_tan( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssr_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!                Add to heat and/or freshwater fluxes a damping term
      !!                toward observed SST and/or SSS.
      !!
      !! ** Method of the direct routine:
      !!            : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS
      !!              - at each nscb time step
      !!                   add a retroaction term on qns    (nn_sstr = 1)
      !!                   add a damping term on emps       (nn_sssr = 1)
      !!                   add a damping term on emp & emps (nn_sssr = 2)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj, ierror   ! dummy loop indices
      REAL(wp) ::   zqrptl   ! local scalar for heat flux damping
      !!
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_ssr_tan')
      !
      !                                               ! -------------------- !
      !IF( kt == nit000 ) THEN                         ! First call kt=nit000 !
         !!                                            ! -------------------- !
         !! Allocate erp and qrp array
         !ALLOCATE( qrp_tl(jpi,jpj), erp_tl(jpi,jpj), STAT=ierror )
         !IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate erp and qrp array' )
         !CALL sbc_ssr_ini_tam ( 0 )
      !ENDIF

      IF( nn_sstr + nn_sssr /= 0 ) THEN

         !                                         ! ========================= !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN      !    Add restoring term     !
            !                                      ! ========================= !
            !
            IF( nn_sstr == 1 ) THEN                   !* Temperature restoring term
!CDIR COLLAPSE
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zqrptl        = rn_dqdt * sst_m_tl(ji,jj)
                     qns_tl(ji,jj) = qns_tl(ji,jj) + zqrptl
                     qrp_tl(ji,jj) = zqrptl
                  END DO
               END DO
            ENDIF
            !
            ! No linear Salinity damping term  (simplification)
            !
         ENDIF
         !
      ENDIF
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_ssr_tan')
      !
   END SUBROUTINE sbc_ssr_tan
   SUBROUTINE sbc_ssr_adj( kt )
      !!---------------------------------------------------------------------
      !!                     ***  ROUTINE sbc_ssr_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!                Add to heat and/or freshwater fluxes a damping term
      !!                toward observed SST and/or SSS.
      !!
      !! ** Method of the direct routine:
      !!            : - Read namelist namsbc_ssr
      !!              - Read observed SST and/or SSS
      !!              - at each nscb time step
      !!                   add a retroaction term on qns    (nn_sstr = 1)
      !!                   add a damping term on emps       (nn_sssr = 1)
      !!                   add a damping term on emp & emps (nn_sssr = 2)
      !!---------------------------------------------------------------------
      INTEGER, INTENT(in   ) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj, ierror   ! dummy loop indices
      REAL(wp) ::   zqrpad   ! local scalar for heat flux damping
      !!
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_ssr_adj')
      !
      zqrpad = 0.0
      !                                               ! -------------------- !
      !IF( kt == nit000 ) THEN                         ! First call kt=nit000 !
         !!                                            ! -------------------- !
         !! Allocate erp and qrp array
         !ALLOCATE( qrp_ad(jpi,jpj), erp_ad(jpi,jpj), STAT=ierror )
         !IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate erp and qrp array' )
         !CALL sbc_ssr_ini_tam ( 1 )
      !ENDIF

      IF( nn_sstr + nn_sssr /= 0 ) THEN
         !                                         ! ========================= !
         IF( MOD( kt-1, nn_fsbc ) == 0 ) THEN      !    Add restoring term     !
            !                                      ! ========================= !
            !
            IF( nn_sstr == 1 ) THEN                   ! Temperature restoring term
!CDIR COLLAPSE
               ! use zqrp scalar to optimize memory access (speedup the loop)
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     zqrpad = qrp_ad(ji,jj)
                     qrp_ad(ji,jj) = 0.0_wp

                     zqrpad = zqrpad + qns_ad(ji,jj)
                     sst_m_ad(ji,jj) = sst_m_ad(ji,jj) + rn_dqdt * zqrpad
                  END DO
               END DO
            ENDIF
            !
            ! No linear Salinity damping term  (simplification)
            !
         ENDIF
         !
      ENDIF
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_ssr_adj')
      !
   END SUBROUTINE sbc_ssr_adj
   SUBROUTINE sbc_ssr_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE sbc_ssr_adj_tst ***
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
      !!        ! 08-11 (A. Vidard)
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
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zsst_m_tlin ,   & ! Tangent input
         & zqns_tlin ,     & ! Tangent input
         & zqns_tlout,     & ! Tangent output
         & zqrp_tlout,     & ! Tangent output
         & zqns_adin ,     & ! Adjoint input
         & zqrp_adin ,     & ! Adjoint input
         & zsst_m_adout,   & ! Adjoint output
         & zqns_adout,     & ! Adjoint output
         & zr                ! 3D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & zqns_tlin(   jpi,jpj),   &
         & zsst_m_tlin( jpi,jpj),   &
         & zqns_tlout(  jpi,jpj),   &
         & zqrp_tlout(  jpi,jpj),   &
         & zqns_adin(   jpi,jpj),   &
         & zqrp_adin(   jpi,jpj),   &
         & zqns_adout(  jpi,jpj),   &
         & zsst_m_adout(jpi,jpj),   &
         & zr(          jpi,jpj)    &
         & )
      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          zqns_tlin(   :,:) = 0.0_wp
          zsst_m_tlin( :,:) = 0.0_wp
          zqns_tlout(  :,:) = 0.0_wp
          zqrp_tlout(  :,:) = 0.0_wp
          zqns_adin(   :,:) = 0.0_wp
          zqrp_adin(   :,:) = 0.0_wp
          zqns_adout(  :,:) = 0.0_wp
          zsst_m_adout(:,:) = 0.0_wp
          zr(          :,:) = 0.0_wp
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random( zr, 'T', 0.0_wp, stdqns )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zqns_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      CALL grid_random( zr, 'T', 0.0_wp, stdt )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsst_m_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      sst_m_tl(:,:) = zsst_m_tlin(:,:)
      qns_tl(  :,:) = zqns_tlin(  :,:)

      CALL sbc_ssr_tan (nit000)
      zqns_tlout(:,:) = qns_tl(:,:)
      zqrp_tlout(:,:) = qrp_tl(:,:)
      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zqns_adin(ji,jj)   = zqns_tlout(ji,jj) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
                 &               * tmask(ji,jj,1)
            END DO
         END DO
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zqrp_adin(ji,jj)   = zqrp_tlout(ji,jj) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
                 &               * tmask(ji,jj,1)
            END DO
         END DO
      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 = DOT_PRODUCT( zqns_tlout, zqns_adin ) &
         & + DOT_PRODUCT( zqrp_tlout, zqrp_adin )

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------
      qns_ad(:,:) = zqns_adin(:,:)
      qrp_ad(:,:) = zqrp_adin(:,:)
      CALL sbc_ssr_adj (nit000)
      zqns_adout(  :,:) = qns_ad(  :,:)
      zsst_m_adout(:,:) = sst_m_ad(:,:)

      zsp2 = DOT_PRODUCT( zqns_tlin,   zqns_adout   ) &
         & + DOT_PRODUCT( zsst_m_tlin, zsst_m_adout )

      ! 14 char:'12345678901234'
      cl_name = 'sbc_ssr_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(        &
         & zqns_tlin,    &
         & zsst_m_tlin,  &
         & zqns_tlout,   &
         & zqrp_tlout,   &
         & zqns_adin,    &
         & zqrp_adin,    &
         & zqns_adout,   &
         & zsst_m_adout, &
         & zr            &
         & )



   END SUBROUTINE sbc_ssr_adj_tst

   SUBROUTINE sbc_ssr_ini_tam
      USE fldread
      INTEGER :: ierror
      CHARACTER(len=100) ::  cn_dir          ! Root directory for location of ssr files
      TYPE(FLD_N) ::   sn_sst, sn_sss        ! informations about the fields to be read
      !!----------------------------------------------------------------------
      NAMELIST/namsbc_ssr/ cn_dir, nn_sstr, nn_sssr, rn_dqdt, rn_deds, sn_sst, &
         &                 sn_sss, ln_sssr_bnd, rn_sssr_bnd

      REWIND ( numnam )         ! ... read in namlist namflx
      READ( numnam, namsbc_ssr )

      IF(lwp) THEN              ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_ssr_tam : SST and/or SSS damping term '
         WRITE(numout,*) '~~~~~~~~~~~ '
         WRITE(numout,*) '          SST restoring term (Yes=1)             nn_sstr = ', nn_sstr
         WRITE(numout,*) '          SSS damping term (Yes=1, salt flux)    nn_sssr = ', nn_sssr
         WRITE(numout,*) '                           (Yes=2, volume flux) '
         WRITE(numout,*) '          dQ/dT (restoring magnitude on SST)     dqdt    = ', rn_dqdt, ' W/m2/K'
         WRITE(numout,*) '          dE/dS (restoring magnitude on SST)     deds    = ', rn_deds, ' mm/day'
      ENDIF
         ALLOCATE( qrp_ad(jpi,jpj), erp_ad(jpi,jpj), &
            &      qrp_tl(jpi,jpj), erp_tl(jpi,jpj), STAT=ierror )
         IF( ierror > 0 )   CALL ctl_stop( 'STOP', 'sbc_ssr: unable to allocate erp and qrp array' )

      !
      ! Initialize qrp and erp if no restoring
      qrp_tl(:,:) = 0.e0
      erp_tl(:,:) = 0.e0
      qrp_ad(:,:) = 0.e0
      erp_ad(:,:) = 0.e0
   END SUBROUTINE sbc_ssr_ini_tam

   !!======================================================================
END MODULE sbcssr_tam
