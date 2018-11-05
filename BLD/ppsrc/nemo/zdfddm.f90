MODULE zdfddm
   !!======================================================================
   !!                       ***  MODULE  zdfddm  ***
   !! Ocean physics : double diffusion mixing parameterization
   !!======================================================================
   !! History :  OPA  ! 2000-08  (G. Madec)  double diffusive mixing
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :          Dummy module          No double diffusion
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfddm = .FALSE.   !: double diffusion flag
CONTAINS
   SUBROUTINE zdf_ddm( kt )           ! Dummy routine
      WRITE(*,*) 'zdf_ddm: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_ddm
   SUBROUTINE zdf_ddm_init            ! Dummy routine
   END SUBROUTINE zdf_ddm_init

   !!======================================================================
END MODULE zdfddm
