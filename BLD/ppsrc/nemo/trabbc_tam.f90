MODULE trabbc_tam
   !!==============================================================================
   !!                       ***  MODULE  trabbc  ***
   !! Ocean active tracers:  bottom boundary condition (geothermal heat flux)
   !!==============================================================================
   !! History :  OPA  ! 1999-10 (G. Madec)  original code
   !!   NEMO     1.0  ! 2002-08 (G. Madec)  free form + modules
   !!             -   ! 2002-11 (A. Bozec)  tra_bbc_init: original code
   !!            3.3  ! 2010-10 (G. Madec)  dynamical allocation + suppression of key_trabbc
   !!             -   ! 2010-11 (G. Madec)  use mbkt array (deepest ocean t-level)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_bbc      : update the tracer trend at ocean bottom
   !!   tra_bbc_init : initialization of geothermal heat flux trend
   !!----------------------------------------------------------------------
   USE oce             ! ocean variables
   USE oce_tam
   USE dom_oce         ! domain: ocean
   USE phycst          ! physical constants
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE trabbc

   IMPLICIT NONE
   PRIVATE

   PUBLIC tra_bbc_tan          ! routine called by step.F90
   PUBLIC tra_bbc_adj          ! routine called by step.F90

   !                                                !!* Namelist nambbc: bottom boundary condition *
   INTEGER         ::   nn_geoflx     = 1            !  Geothermal flux (=1:constant flux, =2:read in file )
   REAL(wp)        ::   rn_geoflx_cst = 86.4e-3_wp   !  Constant value of geothermal heat flux

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_bbc_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbc_tan  ***
      !!
      !! ** Purpose :   Compute the bottom boundary contition on temperature
      !!              associated with geothermal heating and add it to the
      !!              general trend of temperature equations.
      !!
      !! ** Method  :   The geothermal heat flux set to its constant value of
      !!              86.4 mW/m2 (Stein and Stein 1992, Huang 1999).
      !!       The temperature trend associated to this heat flux through the
      !!       ocean bottom can be computed once and is added to the temperature
      !!       trend juste above the bottom at each time step:
      !!            ta = ta + Qsf / (rau0 rcp e3T) for k= mbkt
      !!       Where Qsf is the geothermal heat flux.
      !!
      !! ** Action  : - update the temperature trends (ta) with the trend of
      !!                the ocean bottom boundary condition
      !!
      !! References : Stein, C. A., and S. Stein, 1992, Nature, 359, 123-129.
      !!              Emile-Geay and Madec, 2009, Ocean Science.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, ik    ! dummy loop indices
      REAL(wp) ::   zqgh_trdtl      ! geothermal heat flux trend
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrdt
      !!----------------------------------------------------------------------
      !
      !...Nothing to do...
      !
   END SUBROUTINE tra_bbc_tan


   SUBROUTINE tra_bbc_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_bbc_adj  ***
      !!
      !! ** Purpose :   Compute the bottom boundary contition on temperature
      !!              associated with geothermal heating and add it to the
      !!              general trend of temperature equations.
      !!
      !! ** Method  :   The geothermal heat flux set to its constant value of
      !!              86.4 mW/m2 (Stein and Stein 1992, Huang 1999).
      !!       The temperature trend associated to this heat flux through the
      !!       ocean bottom can be computed once and is added to the temperature
      !!       trend juste above the bottom at each time step:
      !!            ta = ta + Qsf / (rau0 rcp e3T) for k= mbkt
      !!       Where Qsf is the geothermal heat flux.
      !!
      !! ** Action  : - update the temperature trends (ta) with the trend of
      !!                the ocean bottom boundary condition
      !!
      !! References : Stein, C. A., and S. Stein, 1992, Nature, 359, 123-129.
      !!              Emile-Geay and Madec, 2009, Ocean Science.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, ik    ! dummy loop indices
      REAL(wp) ::   zqgh_trdtl      ! geothermal heat flux trend
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   ztrdt
      !!----------------------------------------------------------------------
      !
      !...Nothing to do...
      !
   END SUBROUTINE tra_bbc_adj
   !!======================================================================
END MODULE trabbc_tam
