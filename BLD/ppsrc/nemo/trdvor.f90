MODULE trdvor
   !!======================================================================
   !!                       ***  MODULE  trdvor  ***
   !! Ocean diagnostics:  momentum trends
   !!=====================================================================
   !! History :  1.0  !  04-2006  (L. Brunier, A-M. Treguier) Original code 
   !!            2.0  !  04-2008  (C. Talandier) New trends organization
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
   INTERFACE trd_vor_zint
      MODULE PROCEDURE trd_vor_zint_2d, trd_vor_zint_3d
   END INTERFACE
CONTAINS
   SUBROUTINE trd_vor( kt )        ! Empty routine
      WRITE(*,*) 'trd_vor: You should not have seen this print! error?', kt
   END SUBROUTINE trd_vor
   SUBROUTINE trd_vor_zint_2d( putrdvor, pvtrdvor, ktrd )
      REAL, DIMENSION(:,:), INTENT( inout ) ::   putrdvor, pvtrdvor
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
      WRITE(*,*) 'trd_vor_zint_2d: You should not have seen this print! error?', putrdvor(1,1), pvtrdvor(1,1), ktrd
   END SUBROUTINE trd_vor_zint_2d
   SUBROUTINE trd_vor_zint_3d( putrdvor, pvtrdvor, ktrd )
      REAL, DIMENSION(:,:,:), INTENT( inout ) ::   putrdvor, pvtrdvor
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
      WRITE(*,*) 'trd_vor_zint_3d: You should not have seen this print! error?', putrdvor(1,1,1), pvtrdvor(1,1,1), ktrd
   END SUBROUTINE trd_vor_zint_3d
   SUBROUTINE trd_vor_init              ! Empty routine
      WRITE(*,*) 'trd_vor_init: You should not have seen this print! error?'
   END SUBROUTINE trd_vor_init
   !!======================================================================
END MODULE trdvor
