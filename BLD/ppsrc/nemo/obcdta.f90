MODULE obcdta
   !!==============================================================================
   !!                            ***  MODULE obcdta  ***
   !! Open boundary data : read the data for the open boundaries.
   !!==============================================================================
   !! History :  OPA  ! 1998-05 (J.M. Molines) Original code
   !!            8.5  ! 2002-10 (C. Talandier, A-M. Treguier) Free surface, F90
   !!   NEMO     1.0  ! 2004-06 (F. Durand, A-M. Treguier) Netcdf BC files on input
   !!            3.0  ! 2007-2008 (C. Langlais, P. Mathiot, J.M. Molines) high frequency boundaries data
   !!------------------------------------------------------------------------------
      !!------------------------------------------------------------------------------
      !!   default option:           Dummy module          NO Open Boundary Conditions
      !!------------------------------------------------------------------------------
   CONTAINS
      SUBROUTINE obc_dta( kt )             ! Dummy routine
         INTEGER, INTENT (in) :: kt
         WRITE(*,*) 'obc_dta: You should not have seen this print! error?', kt
      END SUBROUTINE obc_dta
      !!-----------------------------------------------------------------------------
      !!   Default option
      !!-----------------------------------------------------------------------------
      SUBROUTINE obc_dta_bt ( kt, kbt )     ! Empty routine
         INTEGER,INTENT(in) :: kt
         INTEGER, INTENT( in ) ::   kbt     ! barotropic ocean time-step index
         WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kt
         WRITE(*,*) 'obc_dta_bt: You should not have seen this print! error?', kbt
      END SUBROUTINE obc_dta_bt
   !!==============================================================================
   END MODULE obcdta
