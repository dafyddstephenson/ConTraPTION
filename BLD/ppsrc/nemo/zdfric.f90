MODULE zdfric
   !!======================================================================
   !!                       ***  MODULE  zdfric  ***
   !! Ocean physics:  vertical mixing coefficient compute from the local
   !!                 Richardson number dependent formulation
   !!======================================================================
   !! History :  OPA  ! 1987-09  (P. Andrich)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!            7.0  ! 1996-01  (G. Madec)  complete rewriting of multitasking suppression of common work arrays
   !!            8.0  ! 1997-06  (G. Madec)  complete rewriting of zdfmix
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!            3.3.1! 2011-09  (P. Oddo) Mixed layer depth parameterization
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :              NO Richardson dependent vertical mixing
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfric = .FALSE.   !: Richardson mixing flag
CONTAINS
   SUBROUTINE zdf_ric_init         ! Dummy routine
   END SUBROUTINE zdf_ric_init
   SUBROUTINE zdf_ric( kt )        ! Dummy routine
      WRITE(*,*) 'zdf_ric: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_ric

   !!======================================================================
END MODULE zdfric
