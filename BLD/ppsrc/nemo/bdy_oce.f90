MODULE bdy_oce
   !!======================================================================
   !!                       ***  MODULE bdy_oce   ***
   !! Unstructured Open Boundary Cond. :   define related variables
   !!======================================================================
   !! History :  1.0  !  2001-05  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version     
   !!            3.3  !  2010-09  (D. Storkey) add ice boundary conditions
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                NO Unstructured Open Boundary Condition
   !!----------------------------------------------------------------------
   LOGICAL ::   ln_tides = .false.  !: =T apply tidal harmonic forcing along open boundaries

   !!======================================================================
END MODULE bdy_oce

