MODULE trdmld_rst
   !!=================================================================================
   !!                       ***  MODULE  trdmld_rst  ***
   !! Ocean dynamic :  Input/Output files for restart on mixed-layer diagnostics
   !!=================================================================================
   !! History :  9.0  ! 05-05 (C. Deltel) Original code
   !!---------------------------------------------------------------------------------
  !!=================================================================================
  !!                       ***  MODULE  trdmld_rst  ***
  !! Ocean dynamic :  Input/Output files for restart on mixed-layer diagnostics
  !!=================================================================================
CONTAINS
  SUBROUTINE trd_mld_rst_write( kt )           !  No ML diags ==> empty routine
    WRITE(*,*) 'trd_mld_rst_wri: You should not have seen this print! error?', kt
  END SUBROUTINE trd_mld_rst_write
  SUBROUTINE trd_mld_rst_read                  !  No ML Diags ==> empty routine
  END SUBROUTINE trd_mld_rst_read

  !!=================================================================================
END MODULE trdmld_rst
