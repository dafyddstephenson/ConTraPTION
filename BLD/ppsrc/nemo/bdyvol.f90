MODULE bdyvol
   !!======================================================================
   !!                       ***  MODULE  bdyvol  ***
   !! Ocean dynamic :  Volume constraint when unstructured boundary 
   !!                  and filtered free surface are used
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!             -   !  2006-01  (J. Chanut) Bug correction
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_vol( kt )        ! Empty routine
      WRITE(*,*) 'bdy_vol: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_vol

   !!======================================================================
END MODULE bdyvol
