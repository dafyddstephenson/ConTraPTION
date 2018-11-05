MODULE domvvl
   !!======================================================================
   !!                       ***  MODULE domvvl   ***
   !! Ocean : 
   !!======================================================================
   !! History :  2.0  !  2006-06  (B. Levier, L. Marie)  original code
   !!            3.1  !  2009-02  (G. Madec, M. Leclair, R. Benshila)  pure z* coordinate
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                                      Empty routine
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dom_vvl
   END SUBROUTINE dom_vvl
   SUBROUTINE dom_vvl_2(kdum, pudum, pvdum )
      USE par_kind
      INTEGER                   , INTENT(in   ) ::   kdum       
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pudum, pvdum
   END SUBROUTINE dom_vvl_2

   !!======================================================================
END MODULE domvvl
