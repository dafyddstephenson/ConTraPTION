MODULE traldf_tam
   !!======================================================================
   !!                       ***  MODULE  traldf_tam  ***
   !! Ocean Active tracers : lateral diffusive trends
   !!                       Tangent and Adjoint Module
   !!=====================================================================
   !! History of the direct module:
   !!          9.0  ! 05-11 (G. Madec)  Original code
   !! History of the TAM module:
   !!          9.0  ! 08-06 (A. Vidard) Skeleton
   !!          9.0  ! 09-03 (F. Vigilant) adding tra_ldf_lap option
   !!          9.0  ! 10-06 (P.A. Bouttier) adding tra_ldf_bilap option
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf     : update the tracer trend with the lateral diffusion
   !!       ldf_ctl : initialization, namelist read, and parameters control
   !!       ldf_ano : compute lateral diffusion for constant T-S profiles
   !!----------------------------------------------------------------------
   USE traldf_iso_tam
   USE traldf_lap_tam
   USE traldf_bilap_tam
   USE in_out_manager
   USE ldftra_oce
   USE dom_oce
   USE oce
   USE oce_tam
   USE ldfslp
   USE wrk_nemo
   USE timing
   USE prtctl
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_tan     ! called by step_tam.F90
   PUBLIC   tra_ldf_adj     ! called by step_tam.F90
   PUBLIC   tra_ldf_init_tam
   PUBLIC   tra_ldf_adj_tst ! called by tamtst.F90

   INTEGER ::  nldf = 0

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

   SUBROUTINE tra_ldf_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the lateral ocean tracer physics.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_ldf_tan')
      !
      rldf = 1.0
      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 )   ;   CALL tra_ldf_lap_tan   ( kt, nit000, 'TRA', gtsu_tl, gtsv_tl, tsb_tl, tsa_tl, jpts )      ! iso-level laplacian
      CASE ( 1 )   ;   CALL tra_ldf_iso_tan   ( kt, nit000, 'TRA', gtsu_tl, gtsv_tl, tsb_tl, tsa_tl, jpts, ahtb0 )      ! rotated laplacian (except dk[ dk[.] ] part)
      CASE ( 2 )   ;   CALL tra_ldf_bilap_tan ( kt, nit000, 'TRA', gtsu_tl, gtsv_tl, tsb_tl, tsa_tl, jpts )      ! iso-level bilaplacian
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_ldf_tan')
      !
   END SUBROUTINE tra_ldf_tan

   SUBROUTINE tra_ldf_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the lateral ocean tracer physics.
      !!
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      !
      rldf = 1.0
      IF( nn_timing == 1 )  CALL timing_start('tra_ldf_adj')
      !
      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 )   ;   CALL tra_ldf_lap_adj   ( kt, nit000, 'TRA', gtsu_ad, gtsv_ad, tsb_ad, tsa_ad, jpts )      ! rotated laplacian (except dk[ dk[.] ] part)
      CASE ( 1 )   ;   CALL tra_ldf_iso_adj   ( kt, nit000, 'TRA', gtsu_ad, gtsv_ad, tsb_ad, tsa_ad, jpts, ahtb0 )      ! rotated laplacian (except dk[ dk[.] ] part)
      CASE ( 2 )   ;   CALL tra_ldf_bilap_adj ( kt, nit000, 'TRA', gtsu_ad, gtsv_ad, tsb_ad, tsa_ad, jpts )      ! iso-level bilaplacian
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_ldf_adj')
      !
   END SUBROUTINE tra_ldf_adj

   SUBROUTINE tra_ldf_adj_tst( kumadt )
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

      rldf = 1.0
      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
       CASE ( 0 )   ;   CALL tra_ldf_lap_adj_tst( kumadt )
       CASE ( 1 )   ;   CALL tra_ldf_iso_adj_tst( kumadt )
       CASE ( 2 )   ;   CALL tra_ldf_bilap_adj_tst( kumadt )
      END SELECT

   END SUBROUTINE tra_ldf_adj_tst

   SUBROUTINE tra_ldf_init_tam
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_ctl  ***
      !!
      !! ** Purpose :   Choice of the operator for the lateral tracer diffusion
      !!
      !! ** Method  :   set nldf from the namtra_ldf logicals
      !!      nldf == -1   ESOPA test: ALL operators are used
      !!      nldf ==  0   laplacian operator
      !!      nldf ==  1   Rotated laplacian operator
      !!      nldf ==  2   bilaplacian operator
      !!      nldf ==  3   Rotated bilaplacian
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr         ! temporary integers
!
!     NAMELIST/nam_traldf/ ln_traldf_lap  , ln_traldf_bilap,                &
!        &                 ln_traldf_level, ln_traldf_hor, ln_traldf_iso,   &
!        &                 aht0, ahtb0, aeiv0
      !!----------------------------------------------------------------------

      !  Define the lateral mixing oparator for tracers
      ! ===============================================

      ! Namelist nam_traldf already read in ldftra module
!     ! Read Namelist nam_traldf : Lateral physics on tracers
!     REWIND( numnam )
!     READ  ( numnam, namtra_ldf )

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_ldf_init : lateral tracer diffusive operator'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_ldf already read in ldftra module'
         WRITE(numout,*) '   see ldf_tra_init report for lateral mixing parameters'
         WRITE(numout,*)
      ENDIF

      !                               ! control the input
      ioptio = 0
      IF( ln_traldf_lap   )   ioptio = ioptio + 1
      IF( ln_traldf_bilap )   ioptio = ioptio + 1
      IF( ioptio >  1 )   CALL ctl_stop( '          use ONE or NONE of the 2 lap/bilap operator type on tracer' )
      IF( ioptio == 0 )   nldf = -2   ! No lateral diffusion
      ioptio = 0
      IF( ln_traldf_level )   ioptio = ioptio + 1
      IF( ln_traldf_hor   )   ioptio = ioptio + 1
      IF( ln_traldf_iso   )   ioptio = ioptio + 1
      IF( ioptio > 1 )   CALL ctl_stop( '          use only ONE direction (level/hor/iso)' )

      ! defined the type of lateral diffusion from ln_traldf_... logicals
      ! CAUTION : nldf = 1 is used in trazdf_imp, change it carefully
      ierr = 0
      IF( ln_traldf_lap ) THEN       ! laplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_traldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_traldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_traldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_traldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! z-coordinate
           CALL ctl_stop( '          You shouldn t have seen this error message, ln_sco option not impemented yet for tam' )
         ENDIF
      ENDIF

      IF( ln_traldf_bilap ) THEN      ! bilaplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_traldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_traldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   ierr = 2      ! isoneutral (   rotation)
      ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_traldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_traldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ierr == 1 )   CALL ctl_stop( ' iso-level in z-coordinate - partial step, not allowed' )
      IF( ierr == 2 )   CALL ctl_stop( ' isoneutral bilaplacian operator does not exist' )
      IF( lk_traldf_eiv .AND. .NOT.ln_traldf_iso )   &
           CALL ctl_stop( '          eddy induced velocity on tracers',   &
           &              ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion' )
      IF( nldf == 1 .OR. nldf == 3 ) THEN      ! rotation
         IF( .NOT.lk_ldfslp )   CALL ctl_stop( '          the rotation of the diffusive tensor require key_ldfslp' )
         l_traldf_rot = .TRUE.                 ! needed for trazdf_imp
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         IF( nldf == -2 )   WRITE(numout,*) '          NO lateral diffusion'
         IF( nldf ==  0 )   WRITE(numout,*) '          laplacian operator'
         IF( nldf ==  1 )   WRITE(numout,*) '          Rotated laplacian operator'
      ENDIF

   END SUBROUTINE tra_ldf_init_tam
END MODULE traldf_tam
