MODULE dynspg_exp_tam
   !!======================================================================
   !!                   ***  MODULE  dynspg_exp_tam  TANGENT/ADJOINT OF MODULE dynspg_exp***
   !! Ocean dynamics:  surface pressure gradient trend
   !!======================================================================
   !! History of the direct module:
   !!            2.0  !  2005-11  (V. Garnier, G. Madec, L. Bessieres) Original code
   !!            3.2  !  2009-06  (G. Madec, M. Leclair, R. Benshila) introduce sshwzv module
   !! History of the tam module:
   !!            3.2  !  2010-06  (A. Vidard) tam of the 2009-06 version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart explicit free surface
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE dyn_spg_exp_tan( kt )       ! Empty routine
      WRITE(*,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp_tan
   SUBROUTINE dyn_spg_exp_adj( kt )       ! Empty routine
      WRITE(*,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp_adj
   SUBROUTINE dyn_spg_exp_adj_tst( kt )       ! Empty routine
      WRITE(*,*) 'dyn_spg_exp: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_exp_adj_tst

   !!======================================================================
END MODULE dynspg_exp_tam
