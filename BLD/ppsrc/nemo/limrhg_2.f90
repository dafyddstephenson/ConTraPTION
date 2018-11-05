MODULE limrhg_2
   !!======================================================================
   !!                     ***  MODULE  limrhg_2  ***
   !!   Ice rheology :  performs sea ice rheology
   !!======================================================================
   !! History :  0.0  !  1993-12  (M.A. Morales Maqueda.)  Original code
   !!            1.0  !  1994-12  (H. Goosse) 
   !!            2.0  !  2003-12  (C. Ethe, G. Madec)  F90, mpp
   !!             -   !  2006-08  (G. Madec)  surface module, ice-stress at I-point
   !!             -   !  2009-09  (G. Madec)  Huge verctor optimisation
   !!            3.3  !  2009-05  (G.Garric, C. Bricaud) addition of the lim2_evp case
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option        Dummy module      NO VP & LIM-2 sea-ice model
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE lim_rhg_2( k1 , k2 )       ! Dummy routine
      WRITE(*,*) 'lim_rhg_2: You should not have seen this print! error?', k1, k2
   END SUBROUTINE lim_rhg_2

   !!==============================================================================
END MODULE limrhg_2
