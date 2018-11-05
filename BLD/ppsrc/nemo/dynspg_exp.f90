MODULE dynspg_exp
   !!======================================================================
   !!                   ***  MODULE  dynspg_exp  ***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
   !! History :  2.0  !  2005-11  (V. Garnier, G. Madec, L. Bessieres) Original code
   !!            3.2  !  2009-06  (G. Madec, M. Leclair, R. Benshila) introduce sshwzv module
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart explicit free surface 
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dyn_spg_exp( kt )       ! Empty routine
      WRITE(*,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp

   !!======================================================================
END MODULE dynspg_exp
