MODULE dynspg_tam
   !!----------------------------------------------------------------------
   !!    This software is governed by the CeCILL licence (Version 2)
   !!----------------------------------------------------------------------
   !!======================================================================
   !!                       ***  MODULE  dynspg_tam  ***
   !! Ocean dynamics:  surface pressure gradient control
   !!                  Tangent and Adjoint Module
   !!======================================================================
   !! History of the direct module:
   !!            1.0  ! 2005-12  (C. Talandier, G. Madec, V. Garnier)  Original code
   !!            3.2  ! 2009-07  (R. Benshila)  Suppression of rigid-lid option
   !! History of the T&A module:
   !!            9.0  ! 2008-06  (A. Vidard) Skeleton
   !!                 ! 2008-11  (A. Vidard) nemo v3
   !!                 ! 2009-03  (A. Weaver) dynspg_flt_tam
   !!            3.2  ! 2010-04  (F. Vigilant) modification for 3.2
   !!            3.4  ! 2012-07  (P.-A. Bouttier) phasing with 3.2
   !!----------------------------------------------------------------------
   !!   dyn_spg_tan     : update the dynamics trend with the surface pressure
   !!                     gradient (tangent routine)
   !!   dyn_spg_adj     : update the dynamics trend with the surface pressure
   !!                     gradient (adjoint routine)
   !!   dyn_spg_adj_tst : Test of the adjoint routine
   !!----------------------------------------------------------------------
   USE par_oce
   USE phycst
   USE sbc_oce
   USE dom_oce
   USE oce_tam
   USE dynspg_oce
   USE in_out_manager
   USE dynspg_exp_tam ! surface pressure gradient     (dyn_spg_exp routine)
!   USE dynspg_ts_tam  ! surface pressure gradient     (dyn_spg_ts  routine)
   USE dynspg_flt_tam  ! surface pressure gradient     (dyn_spg_flt routine)
   USE lib_mpp        ! MPP library
   USE solver          ! solver initialization
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC dyn_spg_tan,      &   ! routine called by steptan module
      &   dyn_spg_adj,      &   ! routine called by stepadj module
      &   dyn_spg_adj_tst,  &   ! routine controlling adjoint tests
      &   dyn_spg_init_tam

   !! * module variables
   INTEGER ::   nspg = 0   ! type of surface pressure gradient scheme defined from lk_dynspg_...

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

   SUBROUTINE dyn_spg_tan( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!              achieve the momentum time stepping by computing the
      !!              last trend, the surface pressure gradient, and performing
      !!              the Leap-Frog integration.
      !!gm              In the current version only the filtered solution provide
      !!gm            the after velocity, in the 2 other (ua,va) are still the trends
      !!
      !! ** Method  :   Three schemes:
      !!              - explicit computation      : the spg is evaluated at now
      !!              - filtered computation      : the Roulet & madec (2000) technique is used
      !!              - split-explicit computation: a time splitting technique is used
      !!
      !! N.B. : When key_esopa is used all the scheme are tested, regardless
      !!        of the physical meaning of the results.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( IN  ) :: &
         & kt      ! ocean time-step index
      INTEGER, INTENT( OUT ) :: &
         & kindic  ! solver flag
      INTEGER  ::   ji, jj, jk                             ! dummy loop indices
      REAL(wp) ::   z2dt, zg_2                             ! temporary scalar
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_tan')
      !
      SELECT CASE ( nspg )                       ! compute surf. pressure gradient
                                                 ! trend and add it to the general trend
      CASE (  0 )
         CALL dyn_spg_exp_tan( kt )              ! explicit
      CASE (  1 )
         CALL ctl_stop ( 'dyn_spg_ts_tan not available yet' )
      CASE (  2 )
         CALL dyn_spg_flt_tan( kt, kindic )      ! filtered
      !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_tan')
      !
   END SUBROUTINE dyn_spg_tan

   SUBROUTINE dyn_spg_adj( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the lateral ocean dynamics physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( IN  ) :: &
         & kt      ! ocean time-step index
      INTEGER, INTENT( OUT ) :: &
         & kindic  ! solver flag
      !
      INTEGER  ::   ji, jj, jk                             ! dummy loop indices
      REAL(wp) ::   z2dt, zg_2                             ! temporary scalar
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_adj')
      !
      kindic = 0
      spgu_ad(:,:) = 0._wp
      spgv_ad(:,:) = 0._wp

      SELECT CASE ( nspg )                       ! compute surf. pressure gradient
                                                 ! trend and add it to the general trend
      CASE (  0 )
         CALL dyn_spg_exp_adj( kt )              ! explicit
      CASE (  1 )
         CALL ctl_stop ( 'dyn_spg_ts_adj not available yet' )
!!!      CALL dyn_spg_ts_adj ( kt )              ! time-splitting
      CASE (  2 )
         CALL dyn_spg_flt_adj( kt, kindic )      ! filtered
      !
      END SELECT
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_adj')
      !
   END SUBROUTINE dyn_spg_adj

   SUBROUTINE dyn_spg_adj_tst( kumadt )
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
      !! ** Action  : Call the appropriate test routine depending on the
      !!              choice of free surface.
      !!
      !! History :
      !!        ! 09-01 (A. Weaver)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt        ! Output unit

      SELECT CASE ( nspg )
      CASE (  0 )
         CALL dyn_spg_exp_adj_tst( kumadt )      ! explicit
      CASE (  1 )
         CALL ctl_stop ( 'dyn_spg_ts_adj_tst not available yet' )
!!!      CALL dyn_spg_ts_adj_tst ( kumadt )      ! time-splitting
      CASE (  2 )
         CALL dyn_spg_flt_adj_tst( kumadt )      ! filtered
      !
      END SELECT
      !
   END SUBROUTINE dyn_spg_adj_tst

   SUBROUTINE dyn_spg_init_tam
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_spg_ctl_tam  ***
      !!
      !! ** Purpose :  Control the consistency between cpp options for
      !!               surface pressure gradient schemes
      !!----------------------------------------------------------------------
      !! * Local declarations
      INTEGER :: &
        & ioptio

      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_spg_init_tam')
      !
      IF(lwp) THEN             ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_spg_init_tam : choice of the surface pressure gradient scheme'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '     Explicit free surface                  lk_dynspg_exp = ', lk_dynspg_exp
         WRITE(numout,*) '     Free surface with time splitting       lk_dynspg_ts  = ', lk_dynspg_ts
         WRITE(numout,*) '     Filtered free surface cst volume       lk_dynspg_flt = ', lk_dynspg_flt
      ENDIF

      ! Control of surface pressure gradient scheme options
      ! ---------------------------------------------------
      ioptio = 0
      IF(lk_dynspg_exp)   ioptio = ioptio + 1
      IF(lk_dynspg_ts )   ioptio = ioptio + 1
      IF(lk_dynspg_flt)   ioptio = ioptio + 1

      IF( ( ioptio > 1 .AND. .NOT. lk_esopa ) .OR. ioptio == 0 )   &
           &   CALL ctl_stop( ' Choose only one surface pressure gradient scheme with a key cpp' )

      IF( lk_esopa     )   nspg = -1
      IF( lk_dynspg_exp)   nspg =  0
      IF( lk_dynspg_ts )   nspg =  1
      IF( lk_dynspg_flt)   nspg =  2

      IF( lk_esopa     )   nspg = -1

     IF(lwp) THEN
         WRITE(numout,*)
         IF( nspg == -1 )   WRITE(numout,*) '     ESOPA test All scheme used'
         IF( nspg ==  0 )   WRITE(numout,*) '     explicit free surface'
         IF( nspg ==  1 )   WRITE(numout,*) '     free surface with time splitting scheme'
         IF( nspg ==  2 )   WRITE(numout,*) '     filtered free surface'
      ENDIF
      !CALL solver_init( nit000 )   ! Elliptic solver initialisation
      ! Control of timestep choice
      ! --------------------------
      IF( lk_dynspg_ts .OR. lk_dynspg_exp) THEN
         IF( nn_cla == 1 )   &
           &   CALL ctl_stop( ' Crossland advection not implemented for this free surface formulation ' )
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_spg_init_tam')
      !
   END SUBROUTINE dyn_spg_init_tam
  !!======================================================================
END MODULE dynspg_tam
