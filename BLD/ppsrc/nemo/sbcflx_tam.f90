MODULE sbcflx_tam
   !!======================================================================
   !!                       ***  MODULE  sbcflx_tam  ***
   !! Ocean forcing:  momentum, heat and freshwater flux formulation
   !!                 Angent and Andjoint Module
   !!=====================================================================
   !! History of the direct module:
   !!      9.0   !  06-06  (G. Madec)  Original code
   !! History of the T&A module:
   !!      9.0   !  08-11  (A. Vidard)  Original code
   !!      3.2   !  10-04  (A. Vidard)  Update to  new reference
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   namflx   : flux formulation namlist
   !!   sbc_flx  : flux formulation as ocean surface boundary condition
   !!              (forced mode, fluxes read in NetCDF files)
   !!----------------------------------------------------------------------
   !! question diverses
   !!  *   ajouter un test sur la division entier de freqh et rdttra ???
   !!  **  ajoute dans namelist: 1 year forcing files
   !!                         or forcing file starts at the begining of the run
   !!  *** we assume that the forcing file start and end with the previous
   !!      year last record and the next year first record (useful for
   !!      time interpolation, required even if no time interp???)
   !!  *   ajouter un test sur la division de la frequence en pas de temps
   !!  ==> daymod ajout de nsec_year = number of second since the begining of the year
   !!      assumed to be 0 at 0h january the 1st (i.e. 24h december the 31)
   !!
   !!  *** regrouper dtatem et dtasal
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE sbc_oce_tam

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_flx_tan       ! routine called by sbcmod_tam.F90
   PUBLIC sbc_flx_adj       ! routine called by sbcmod_tam.F90
   PUBLIC sbc_flx_adj_tst   ! routine called by tst.F90

   REAL(wp) ::   rhoa  = 1.22         ! Air density kg/m3
   REAL(wp) ::   cdrag = 1.5e-3       ! drag coefficient

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
   !!   OPA 9.0 , LOCEAN-IPSL (2006)
   !! $Id: sbcflx.F90 1200 2008-09-24 13:05:20Z rblod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE sbc_flx_tan( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!                provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  : - READ each fluxes in NetCDF files:
      !!                   i-component of the stress              utau  (N/m2)
      !!                   j-component of the stress              vtau  (N/m2)
      !!                   net downward heat flux                 qtot  (watt/m2)
      !!                   net downward radiative flux            qsr   (watt/m2)
      !!                   net upward freshwater (evapo - precip) emp   (kg/m2/s)
      !!
      !!      CAUTION :  - never mask the surface stress fields
      !!                 - the stress is assumed to be in the mesh referential
      !!                   i.e. the (i,j) referential
      !!
      !! ** Action  :   update at each time-step
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point
      !!              - qns, qsr    non-slor and solar heat flux
      !!              - emp, emps   evaporation minus precipitation
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj          ! dummy indices
     !!
      !!---------------------------------------------------------------------
      ! set the ocean fluxes from read fields
!CDIR COLLAPSE
      DO jj = 1, jpj
         DO ji = 1, jpi
            utau_tl(ji,jj) = 0.0_wp
            vtau_tl(ji,jj) = 0.0_wp
            qns_tl (ji,jj) = 0.0_wp
            qsr_tl (ji,jj) = 0.0_wp
            emp_tl (ji,jj) = 0.0_wp
         END DO
      END DO

      ! Initialization of emps (when no ice model)
      emps_tl(:,:) = emp_tl (:,:)

      ! Estimation of wind speed as a function of wind stress ( |tau|=rhoa*Cd*|U|^2 )
      ! WARNING !!! TANGENT NOT VALID if TAU is part of the state vector.
      taum_tl(:,:) = 0.0_wp
      wndm_tl(:,:) = 0.0_wp
      !
   END SUBROUTINE sbc_flx_tan
   SUBROUTINE sbc_flx_adj( kt )
      !!---------------------------------------------------------------------
      !!                    ***  ROUTINE sbc_flx_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!                provide at each time step the surface ocean fluxes
      !!                (momentum, heat, freshwater and runoff)
      !!
      !! ** Method  : - READ each fluxes in NetCDF files:
      !!                   i-component of the stress              utau  (N/m2)
      !!                   j-component of the stress              vtau  (N/m2)
      !!                   net downward heat flux                 qtot  (watt/m2)
      !!                   net downward radiative flux            qsr   (watt/m2)
      !!                   net upward freshwater (evapo - precip) emp   (kg/m2/s)
      !!
      !!      CAUTION :  - never mask the surface stress fields
      !!                 - the stress is assumed to be in the mesh referential
      !!                   i.e. the (i,j) referential
      !!
      !! ** Action  :   update at each time-step
      !!              - utau, vtau  i- and j-component of the wind stress
      !!              - taum        wind stress module at T-point
      !!              - wndm        10m wind module at T-point
      !!              - qns, qsr    non-slor and solar heat flux
      !!              - emp, emps   evaporation minus precipitation
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time step
      !!
      INTEGER  ::   ji, jj   ! dummy indices
      !!
      !!---------------------------------------------------------------------
      ! Estimation of wind speed as a function of wind stress ( |tau|=rhoa*Cd*|U|^2 )
      ! WARNING !!! ADJGENT NOT VALID if TAU is par of the state vector.
      taum_ad(:,:) = 0.0_wp
      wndm_ad(:,:) = 0.0_wp
      ! Initialization of emps (when no ice model)
      emp_ad (:,:) =  emp_ad (:,:) + emps_ad(:,:)
      emps_ad(:,:) = 0.0_wp
       ! set the ocean fluxes from read fields
!CDIR COLLAPSE
      DO jj = 1, jpj
         DO ji = 1, jpi
            utau_ad(ji,jj) = 0.0_wp
            vtau_ad(ji,jj) = 0.0_wp
            qns_ad (ji,jj) = 0.0_wp
            qsr_ad (ji,jj) = 0.0_wp
            emp_ad (ji,jj) = 0.0_wp
         END DO
      END DO
     !
   END SUBROUTINE sbc_flx_adj
   SUBROUTINE sbc_flx_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE sbc_flx_adj_tst ***
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

      ! adjoint is obvious, no need to test...
   END SUBROUTINE sbc_flx_adj_tst
   !!======================================================================
END MODULE sbcflx_tam
