MODULE trdmld
   !!======================================================================
   !!                       ***  MODULE  trdmld  ***
   !! Ocean diagnostics:  mixed layer T-S trends 
   !!=====================================================================
   !! History :       !  95-04  (J. Vialard)    Original code
   !!                 !  97-02  (E. Guilyardi)  Adaptation global + base cmo
   !!                 !  99-09  (E. Guilyardi)  Re-writing + netCDF output
   !!            8.5  !  02-06  (G. Madec)      F90: Free form and module
   !!            9.0  !  04-08  (C. Talandier)  New trends organization
   !!                 !  05-05  (C. Deltel)     Diagnose trends of time averaged ML T & S
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trd_mld( kt )             ! Empty routine
      INTEGER, INTENT( in) ::   kt
      WRITE(*,*) 'trd_mld: You should not have seen this print! error?', kt
   END SUBROUTINE trd_mld
   SUBROUTINE trd_mld_zint( pttrdmld, pstrdmld, ktrd, ctype )
      REAL, DIMENSION(:,:,:), INTENT( in ) ::   &
         pttrdmld, pstrdmld                   ! Temperature and Salinity trends
      INTEGER, INTENT( in ) ::   ktrd         ! ocean trend index
      CHARACTER(len=2), INTENT( in ) ::   &  
         ctype                                ! surface/bottom (2D arrays) or
         !                                    ! interior (3D arrays) physics
      WRITE(*,*) 'trd_mld_zint: You should not have seen this print! error?', pttrdmld(1,1,1)
      WRITE(*,*) '  "      "  : You should not have seen this print! error?', pstrdmld(1,1,1)
      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ctype
      WRITE(*,*) '  "      "  : You should not have seen this print! error?', ktrd
   END SUBROUTINE trd_mld_zint
   SUBROUTINE trd_mld_init              ! Empty routine
      WRITE(*,*) 'trd_mld_init: You should not have seen this print! error?'
   END SUBROUTINE trd_mld_init

   !!======================================================================
END MODULE trdmld
