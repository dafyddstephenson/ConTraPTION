MODULE tstool_tam
   !!==========================================================================
   !!     ***  MODULE  tstool_tam : TAM testing utilities  ***
   !!==========================================================================
   !! History of the NEMOTAM module:
   !!             3.0  !  08-11  (A. Vidard) initial version
   USE par_oce,        ONLY: & ! Ocean space and time domain variables
      & jpi,                 &
      & jpj,                 &
      & jpk,                 &
      & jpiglo,              &
      & jpjglo
   USE dom_oce       , ONLY: & ! Ocean space and time domain
      & e1u,                 &
      & e2u,                 &
      & e3t_0,               &
      & e3u,                 &
      & umask,               &
      & mig,                 &
      & mjg,                 &
      & nldi,                &
      & nldj,                &
      & nlei,                &
      & nlej
   USE wrk_nemo       ! Memory Allocation
   USE par_kind      , ONLY: & ! Precision variables
      & wp
   USE in_out_manager, ONLY: & ! I/O manager
      & lwp
   USE gridrandom    , ONLY: & ! Random Gaussian noise on grids
      & grid_random
   USE dotprodfld    , ONLY: & ! Computes dot product for 3D and 2D fields
      & dot_product
   IMPLICIT NONE
   PRIVATE
   REAL(KIND=wp), PUBLIC ::  & ! random field standard deviation for:
      & stdu   =  0.01_wp,    & !   u-velocity
      & stdv   =  0.01_wp,    & !   v-velocity
      & stdw   = 0.01_wp,    & !   w-velocity
      & stds   =  0.01_wp,    & !   salinity
      & stdt   =  0.1_wp,    & !   temperature
      & stdssh = 0.01_wp,    & !   sea surface height
      & stdemp = 0.001_wp,    & !   evaporation minus precip 0.1_wp / SQRT( wesp_emp )
      & stdqns =  0.1_wp,    & !   non solar heat flux
      & stdqsr =  0.1_wp,    & !   solar heat flux
      & stdgc  =  0.01_wp,    & !   gcx, gcb
      & stdr   =  0.01_wp,    & !   rotb, rhd
      & stdh   =  0.01_wp       !   hdivb

   PUBLIC &
      & prntst_adj,          &
      & prntst_tlm

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

CONTAINS

   SUBROUTINE prntst_adj( cd_name, kumadt, psp1, psp2 )
      CHARACTER(LEN=14), INTENT(in) :: cd_name
      REAL(wp), INTENT(in) :: psp1, psp2
      INTEGER, INTENT(in) :: kumadt
      REAL(KIND=wp) :: &
         & zspdif,       & ! scalar product difference
         & ztol            ! accepted tolerance
      CHARACTER (LEN=47) :: &
         & FMT
      CHARACTER (LEN=9)  :: &
         & cl_stat         ! Accuracy status of adjoint routine (ok or warning)

      ! Compare the scalar products

      zspdif = ABS( psp1 - psp2 )
      IF ( psp1 /= 0.0_wp ) zspdif = zspdif / ABS( psp1 )

      ztol = EPSILON( zspdif ) * 10.

      IF ( zspdif < ztol ) THEN
         cl_stat = '   ok    '
      ELSEIF ( zspdif < ztol*1000._wp ) THEN
         cl_stat = ' warning '
      ELSE
         cl_stat = 'RED ALERT'
      ENDIF

      IF (lwp) THEN
         FMT = "(A14,1X,E20.15,2X,E20.15,2X,E6.1,1X,E6.1,1x,A9)"
         WRITE(kumadt,FMT) cd_name, psp1, psp2, zspdif, ztol, cl_stat
         CALL FLUSH( kumadt )
      ENDIF
   END SUBROUTINE prntst_adj

   SUBROUTINE prntst_tlm( cd_name, kumadt, psp1, psp2 )
      CHARACTER(LEN=14), INTENT(in) :: cd_name
      REAL(wp), INTENT(in) :: psp1, psp2
      INTEGER, INTENT(in) :: kumadt
      REAL(KIND=wp) :: &
         & zspratio       ! scalar product difference
      CHARACTER (LEN=47) :: &
         & FMT
      CHARACTER (LEN=9)  :: &
         & cl_stat         ! Accuracy status of adjoint routine (ok or warning)

      ! Compare the scalar products

      IF ( psp1 /= 0.0_wp ) zspratio = 100 * psp1 /  psp2

      IF (lwp) THEN
         FMT = "(A14,1X,E20.13,2X,E20.15,2X,E6.1,1X)"
         WRITE(kumadt,FMT) cd_name, psp1, psp2, zspratio
         CALL FLUSH( kumadt )
      ENDIF
   END SUBROUTINE prntst_tlm

   SUBROUTINE example_adj_tst( kumadt )
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
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), POINTER, DIMENSION(:,:,:) :: &
         & z_tlin ,     & ! Tangent input
         & z_tlout,     & ! Tangent output
         & z_adin ,     & ! Adjoint input
         & z_adout,     & ! Adjoint output
         & zr             ! 3D random field
      CHARACTER(LEN=14) :: &
         & cl_name
      ! Allocate memory

      CALL wrk_alloc( jpi, jpj, jpk, z_tlin , z_tlout, z_adin, z_adout, zr )

      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          z_tlin( :,:,:) = 0.0_wp
          z_tlout(:,:,:) = 0.0_wp
          z_adin( :,:,:) = 0.0_wp
          z_adout(:,:,:) = 0.0_wp
          zr(     :,:,:) = 0.0_wp
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random( zr, 'U', 0.0_wp, stdr )
      z_tlin(:,:,:) = zr(:,:,:)

      CALL example_tan
      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              z_adin(ji,jj,jk) = z_tlout(ji,jj,jk) &
                 &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
                 &               * umask(ji,jj,jk)
            END DO
         END DO
      END DO
      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 = DOT_PRODUCT( z_tlout, z_adin )

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      CALL example_adj

      zsp2 = DOT_PRODUCT( z_tlin, z_adout )

      ! 14 char:'12345678901234'
      cl_name = 'example_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      CALL wrk_dealloc( jpi, jpj, jpk, z_tlin , z_tlout, z_adin, z_adout, zr )

   END SUBROUTINE example_adj_tst

   SUBROUTINE example_tan
   END SUBROUTINE example_tan
   SUBROUTINE example_adj
   END SUBROUTINE example_adj

END MODULE tstool_tam
