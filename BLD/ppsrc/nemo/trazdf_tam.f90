MODULE trazdf_tam
   !!==============================================================================
   !!                 ***  MODULE  trazdf_zdf  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!                 Tangent and Adjoint Module
   !!==============================================================================
   !! History of the direct module:
   !!            9.0  !  05-11  (G. Madec)  Original code
   !! History of the TAM module:
   !!            9.0  !  08-06  (A. Vidard) Skeleton
   !!            9.0  !  09-01  (A. Vidard) TAM of the 05-11 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_zdf_tan  : Update the tracer trend with the vertical diffusion (tan)
   !!   tra_zdf_adj  : Update the tracer trend with the vertical diffusion (adj)
   !!       zdf_ctl  : ???
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE ldftra_oce
   USE zdf_oce
   USE trazdf_exp_tam
   USE trazdf_imp_tam
   USE in_out_manager
   USE prtctl
   USE lib_mpp
   USE wrk_nemo
   USE timing
   USE phycst

   IMPLICIT NONE
   PRIVATE

   PUBLIC  &
      & tra_zdf_tan, &
      & tra_zdf_adj         ! routines called by step_tam.F90
   PUBLIC  tra_zdf_adj_tst  ! routine called by tst.F90
   PUBLIC  tra_zdf_init_tam
   INTEGER ::   nzdf = 0               ! type vertical diffusion algorithm used
      !                                ! defined from ln_zdf...  namlist logicals)

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
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   Defautl option :                     avs = avt
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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

   SUBROUTINE tra_zdf_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the vertical ocean tracer physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_tan')
      !
      !                                             ! set time step
      IF( neuler == 0 .AND. kt == nit000 ) THEN     ! at nit000
         r2dtra =  rdttra(:)                          ! = rdtra (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1) THEN                ! at nit000 or nit000+1
         r2dtra = 2. * rdttra(:)                      ! = 2 rdttra (leapfrog)
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( -1 )                                       ! esopa: test all possibility with control print
         CALL tra_zdf_exp_tan    ( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb_tl, tsa_tl, jpts )
         CALL tra_zdf_imp_tan    ( kt, nit000, 'TRA', r2dtra,            tsb_tl, tsa_tl, jpts  )

      CASE ( 0 )                                       ! explicit scheme
         CALL tra_zdf_exp_tan    ( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb_tl, tsa_tl, jpts  )

      CASE ( 1 )                                       ! implicit scheme
         CALL tra_zdf_imp_tan    ( kt, nit000, 'TRA', r2dtra,            tsb_tl, tsa_tl, jpts  )

      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_tan')
      !
   END SUBROUTINE tra_zdf_tan
   SUBROUTINE tra_zdf_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the vertical ocean tracer physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index

      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_adj')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000 ) THEN     ! at nit000
         r2dtra =  rdttra(:)                          ! = rdtra (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1) THEN                ! at nit000 or nit000+1
         r2dtra = 2. * rdttra(:)                      ! = 2 rdttra (leapfrog)
      ENDIF

      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( -1 )                                       ! esopa: test all possibility with control print
         CALL tra_zdf_exp_adj    ( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb_ad, tsa_ad, jpts )
         CALL tra_zdf_imp_adj    ( kt, nit000, 'TRA', r2dtra,            tsb_ad, tsa_ad, jpts )

      CASE ( 0 )                                       ! explicit scheme
         CALL tra_zdf_exp_adj    ( kt, nit000, 'TRA', r2dtra, nn_zdfexp, tsb_ad, tsa_ad, jpts )

      CASE ( 1 )                                       ! implicit scheme
         CALL tra_zdf_imp_adj    ( kt, nit000, 'TRA', r2dtra,            tsb_ad, tsa_ad, jpts )

      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_adj')
      !
   END SUBROUTINE tra_zdf_adj
   SUBROUTINE tra_zdf_adj_tst( kumadt )
      !! ** Purpose : Test the adjoint routines.
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
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      !! * Local declarations
      ! init
      CALL tra_zdf_init_tam
      ! Test the explicit formulation
      CALL tra_zdf_exp_adj_tst    ( kumadt )
      ! Test the implicit formulation
      CALL tra_zdf_imp_adj_tst    ( kumadt )
   END SUBROUTINE tra_zdf_adj_tst
   !!==============================================================================
   SUBROUTINE tra_zdf_init_tam
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE zdf_ctl_tam  ***
      !!
      !! ** Purpose :   Choose the vertical mixing scheme
      !!
      !! ** Method  :   Set nzdf from ln_zdfexp
      !!      nzdf = 0   explicit (time-splitting) scheme (ln_zdfexp=T)
      !!           = 1   implicit (euler backward) scheme (ln_zdfexp=F)
      !!      NB: rotation of lateral mixing operator or TKE or KPP scheme,
      !!      the implicit scheme is required.
      !!----------------------------------------------------------------------
      USE zdftke
      USE zdfkpp
      USE zdfgls
      !!----------------------------------------------------------------------

      !  Define the vertical tracer physics scheme
      ! ==========================================

      ! Choice from ln_zdfexp already read in namelist in zdfini module
      IF( ln_zdfexp ) THEN               ! use explicit scheme
         nzdf = 0
      ELSE                               ! use implicit scheme
         nzdf = 1
      ENDIF

      ! Force implicit schemes
      IF( lk_zdfgls .OR. lk_zdftke .OR. lk_zdfkpp       )   nzdf = 1      ! TKE or KPP physics
      IF( ln_traldf_iso                                 )   nzdf = 1      ! iso-neutral lateral physics
      IF( ln_traldf_hor .AND. ln_sco                    )   nzdf = 1      ! horizontal lateral physics in s-coordinate

      IF( ln_zdfexp .AND. nzdf == 1 )   THEN
         CALL ctl_stop( 'tra_zdf_tam : If using the rotation of lateral mixing operator or TKE ', &
            &           '            or KPP scheme, the implicit scheme is required, set ln_zdfexp = .false.' )
      ENDIF

      ! Test: esopa
      IF( lk_esopa )    nzdf = -1                      ! All schemes used

      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra_zdf_init_tam : vertical tracer physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         IF( nzdf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF
   END SUBROUTINE tra_zdf_init_tam
END MODULE trazdf_tam
