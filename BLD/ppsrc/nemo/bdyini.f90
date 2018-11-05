MODULE bdyini
   !!======================================================================
   !!                       ***  MODULE  bdyini  ***
   !! Unstructured open boundaries : initialisation
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-01  (D. Storkey) Tidal forcing
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) updates for Shelf configurations
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!---------------------------------------------------------------------------------
   !!   Dummy module                                   NO open boundaries
   !!---------------------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_init      ! Dummy routine
   END SUBROUTINE bdy_init

   !!=================================================================================
END MODULE bdyini
