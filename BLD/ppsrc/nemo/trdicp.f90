MODULE trdicp
   !!======================================================================
   !!                       ***  MODULE  trdicp  ***
   !! Ocean diagnostics:  ocean tracers and dynamic trends
   !!=====================================================================
   !! History :  1.0  !  2004-08 (C. Talandier) New trends organization
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :                                         Empty module
   !!----------------------------------------------------------------------
   INTERFACE trd_icp
      MODULE PROCEDURE trd_2d, trd_3d
   END INTERFACE

CONTAINS
   SUBROUTINE trd_2d( ptrd2dx, ptrd2dy, ktrd , ctype )       ! Empty routine
      REAL, DIMENSION(:,:) ::   ptrd2dx, ptrd2dy
      INTEGER                     , INTENT(in   ) ::   ktrd         ! tracer trend index
      CHARACTER(len=3)            , INTENT(in   ) ::   ctype        ! momentum ('DYN') or tracers ('TRA') trends
      WRITE(*,*) 'trd_2d: You should not have seen this print! error ?', &
          &       ptrd2dx(1,1), ptrd2dy(1,1), ktrd, ctype
   END SUBROUTINE trd_2d
   SUBROUTINE trd_3d( ptrd3dx, ptrd3dy, ktrd , ctype )       ! Empty routine
      REAL, DIMENSION(:,:,:) ::   ptrd3dx, ptrd3dy
      INTEGER                     , INTENT(in   ) ::   ktrd         ! tracer trend index
      CHARACTER(len=3)            , INTENT(in   ) ::   ctype        ! momentum ('DYN') or tracers ('TRA') trends
      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', &
          &       ptrd3dx(1,1,1), ptrd3dy(1,1,1), ktrd, ctype
   END SUBROUTINE trd_3d
   SUBROUTINE trd_icp_init               ! Empty routine
   END SUBROUTINE trd_icp_init
   SUBROUTINE trd_dwr( kt )          ! Empty routine
      WRITE(*,*) 'trd_dwr: You should not have seen this print! error ?', kt
   END SUBROUTINE trd_dwr
   SUBROUTINE trd_twr( kt )          ! Empty routine
      WRITE(*,*) 'trd_twr: You should not have seen this print! error ?', kt
   END SUBROUTINE trd_twr

   !!======================================================================
END MODULE trdicp
