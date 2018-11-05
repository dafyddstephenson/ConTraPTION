MODULE sbcice_lim
   !!======================================================================
   !!                       ***  MODULE  sbcice_lim  ***
   !! Surface module :  update the ocean surface boundary condition over ice
   !!       &           covered area using LIM sea-ice model
   !! Sea-Ice model  :  LIM-3 Sea ice model time-stepping
   !!=====================================================================
   !! History :  2.0  ! 2006-12  (M. Vancoppenolle) Original code
   !!            3.0  ! 2008-02  (C. Talandier)  Surface module from icestp.F90
   !!             -   ! 2008-04  (G. Madec)  sltyle and lim_ctl routine
   !!            3.3  ! 2010-11  (G. Madec) ice-ocean stress always computed at each ocean time-step
   !!            4.0  ! 2011-01  (A Porter)  dynamical allocation
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option           Dummy module      NO LIM 3.0 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE sbc_ice_lim ( kt, kblk )     ! Dummy routine
      WRITE(*,*) 'sbc_ice_lim: You should not have seen this print! error?', kt, kblk
   END SUBROUTINE sbc_ice_lim

   !!======================================================================
END MODULE sbcice_lim
