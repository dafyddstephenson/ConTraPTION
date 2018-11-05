MODULE bdytides
   !!======================================================================
   !!                       ***  MODULE  bdytides  ***
   !! Ocean dynamics:   Tidal forcing at open boundaries
   !!======================================================================
   !! History :  2.0  !  2007-01  (D.Storkey)  Original code
   !!            2.3  !  2008-01  (J.Holt)  Add date correction. Origins POLCOMS v6.3 2007
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D.Storkey and E.O'Dea)  bug fixes
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module         NO Unstruct Open Boundary Conditions for tides
   !!----------------------------------------------------------------------
!!gm  are you sure we need to define filtide and tide_cpt ?
   CHARACTER(len=80), PUBLIC               ::   filtide                !: Filename root for tidal input files
   CHARACTER(len=4 ), PUBLIC, DIMENSION(1) ::   tide_cpt               !: Names of tidal components used.

CONTAINS
   SUBROUTINE tide_init                ! Empty routine
   END SUBROUTINE tide_init
   SUBROUTINE tide_data                ! Empty routine
   END SUBROUTINE tide_data
   SUBROUTINE tide_update( kt, kit )   ! Empty routine
      WRITE(*,*) 'tide_update: You should not have seen this print! error?', kt, kit
   END SUBROUTINE tide_update

   !!======================================================================
END MODULE bdytides
