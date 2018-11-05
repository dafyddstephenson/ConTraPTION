MODULE bdydyn3d
   !!======================================================================
   !!                       ***  MODULE  bdydyn3d  ***
   !! Unstructured Open Boundary Cond. :   Flow relaxation scheme on baroclinic velocities
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn3d( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn_frs: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn3d

   !!======================================================================
END MODULE bdydyn3d
