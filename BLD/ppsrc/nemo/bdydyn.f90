MODULE bdydyn
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Apply boundary conditions to velocities
   !!======================================================================
   !! History :  1.0  !  2005-02  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2007-07  (D. Storkey) Move Flather implementation to separate routine.
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.2  !  2008-04  (R. Benshila) consider velocity instead of transport 
   !!            3.3  !  2010-09  (E.O'Dea) modifications for Shelf configurations 
   !!            3.3  !  2010-09  (D.Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn

   !!======================================================================
END MODULE bdydyn
