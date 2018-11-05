MODULE sbctide
  !!=================================================================================
  !!                       ***  MODULE  sbctide  ***
  !! Initialization of tidal forcing
  !! History :  9.0  !  07  (O. Le Galloudec)  Original code
  !!=================================================================================
  !! * Modules used
  USE oce             ! ocean dynamics and tracers variables
  USE dom_oce         ! ocean space and time domain
  USE in_out_manager  ! I/O units
  USE ioipsl          ! NetCDF IPSL library
  USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
  USE phycst
  USE daymod
  USE dynspg_oce
  USE tide_mod
  USE iom

  IMPLICIT NONE
  PUBLIC

  REAL(wp), PUBLIC, ALLOCATABLE, DIMENSION(:,:) :: pot_astro
  LOGICAL, PUBLIC :: ln_tide_pot = .false.
  !!----------------------------------------------------------------------
  !!   Default case :   Empty module
  !!----------------------------------------------------------------------
  LOGICAL, PUBLIC, PARAMETER ::   lk_tide = .FALSE.
CONTAINS
  SUBROUTINE sbc_tide( kt )      ! Empty routine
    INTEGER         , INTENT(in) ::   kt         ! ocean time-step
    WRITE(*,*) 'sbc_tide: You should not have seen this print! error?', kt
  END SUBROUTINE sbc_tide
  !!======================================================================

END MODULE sbctide
