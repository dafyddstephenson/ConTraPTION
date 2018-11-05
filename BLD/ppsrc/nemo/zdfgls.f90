MODULE zdfgls
   !!======================================================================
   !!                       ***  MODULE  zdfgls  ***
   !! Ocean physics:  vertical mixing coefficient computed from the gls 
   !!                 turbulent closure parameterization
   !!======================================================================
   !! History :   3.0  !  2009-09  (G. Reffray)  Original code
   !!             3.3  !  2010-10  (C. Bricaud)  Add in the reference
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        NO TKE scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfgls = .FALSE.   !: TKE flag
CONTAINS
   SUBROUTINE zdf_gls_init           ! Empty routine
      WRITE(*,*) 'zdf_gls_init: You should not have seen this print! error?'
   END SUBROUTINE zdf_gls_init
   SUBROUTINE zdf_gls( kt )          ! Empty routine
      WRITE(*,*) 'zdf_gls: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_gls
   SUBROUTINE gls_rst( kt, cdrw )          ! Empty routine
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      WRITE(*,*) 'gls_rst: You should not have seen this print! error?', kt, cdrw
   END SUBROUTINE gls_rst

   !!======================================================================
END MODULE zdfgls

