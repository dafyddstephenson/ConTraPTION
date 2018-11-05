MODULE obcrad 
   !!=================================================================================
   !!                       ***  MODULE  obcrad  ***
   !! Ocean dynamic :   Phase velocities for each open boundary
   !!=================================================================================
   !!=================================================================================
   !!                       ***  MODULE  obcrad  ***
   !! Ocean dynamic :   Phase velocities for each open boundary
   !!=================================================================================
CONTAINS
   SUBROUTINE obc_rad( kt )            ! No open boundaries ==> empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'obc_rad: You should not have seen this print! error?', kt
   END SUBROUTINE obc_rad

END MODULE obcrad
