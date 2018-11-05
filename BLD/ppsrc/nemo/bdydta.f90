MODULE bdydta
   !!======================================================================
   !!                       ***  MODULE bdydta  ***
   !! Open boundary data : read the data for the unstructured open boundaries.
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-01  (D. Storkey) Update to use IOM module
   !!             -   !  2007-07  (D. Storkey) add bdy_dta_fla
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dta( kt, jit, time_offset ) ! Empty routine
      INTEGER, INTENT( in )           ::   kt    
      INTEGER, INTENT( in ), OPTIONAL ::   jit   
      INTEGER, INTENT( in ), OPTIONAL ::   time_offset
      WRITE(*,*) 'bdy_dta: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dta
   SUBROUTINE bdy_dta_init()                  ! Empty routine
      WRITE(*,*) 'bdy_dta_init: You should not have seen this print! error?'
   END SUBROUTINE bdy_dta_init

   !!==============================================================================
END MODULE bdydta
