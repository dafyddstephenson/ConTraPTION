MODULE zdfkpp
   !!======================================================================
   !!                       ***  MODULE  zdfkpp  ***
   !! Ocean physics:  vertical mixing coefficient compute from the KPP 
   !!                 turbulent closure parameterization
   !!=====================================================================
   !! History :  OPA  ! 2000-03 (W.G. Large, J. Chanut) Original code
   !!            8.1  ! 2002-06 (J.M. Molines) for real case CLIPPER  
   !!            8.2  ! 2003-10 (Chanut J.) re-writting
   !!   NEMO     1.0  ! 2005-01 (C. Ethe, G. Madec) Free form, F90 + creation of tra_kpp routine
   !!            3.3  ! 2010-10 (C. Ethe, G. Madec) reorganisation of initialisation phase + merge TRC-TRA
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                                        NO KPP scheme
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_zdfkpp = .FALSE.   !: KPP flag
CONTAINS
   SUBROUTINE zdf_kpp_init           ! Dummy routine
      WRITE(*,*) 'zdf_kpp_init: You should not have seen this print! error?'
   END SUBROUTINE zdf_kpp_init
   SUBROUTINE zdf_kpp( kt )          ! Dummy routine
      WRITE(*,*) 'zdf_kpp: You should not have seen this print! error?', kt
   END SUBROUTINE zdf_kpp
   SUBROUTINE tra_kpp( kt )          ! Dummy routine
      WRITE(*,*) 'tra_kpp: You should not have seen this print! error?', kt
   END SUBROUTINE tra_kpp
   SUBROUTINE trc_kpp( kt )          ! Dummy routine
      WRITE(*,*) 'trc_kpp: You should not have seen this print! error?', kt
   END SUBROUTINE trc_kpp

   !!======================================================================
END MODULE zdfkpp
