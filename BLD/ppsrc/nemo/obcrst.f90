MODULE obcrst
   !!=================================================================================
   !!                       ***  MODULE  obcrst  ***
   !! Ocean dynamic :  Input/Output files for restart on OBC
   !!=================================================================================
CONTAINS
   SUBROUTINE obc_rst_write( kt )           !  No Open boundary ==> empty routine
      INTEGER,INTENT(in) :: kt
      WRITE(*,*) 'obc_rst_write: You should not have seen this print! error?', kt
   END SUBROUTINE obc_rst_write
   SUBROUTINE obc_rst_read                 !  No Open boundary ==> empty routine
   END SUBROUTINE obc_rst_read

   !!=================================================================================
END MODULE obcrst
