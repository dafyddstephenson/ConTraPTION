MODULE bdydyn2d
   !!======================================================================
   !!                       ***  MODULE  bdydyn  ***
   !! Unstructured Open Boundary Cond. :   Apply boundary conditions to barotropic solution
   !!======================================================================
   !! History :  3.4  !  2011     (D. Storkey) new module as part of BDY rewrite
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                   NO Unstruct Open Boundary Conditions
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE bdy_dyn2d( kt )      ! Empty routine
      WRITE(*,*) 'bdy_dyn_frs: You should not have seen this print! error?', kt
   END SUBROUTINE bdy_dyn2d

   !!======================================================================
END MODULE bdydyn2d
