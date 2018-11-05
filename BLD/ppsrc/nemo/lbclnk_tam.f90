MODULE lbclnk_tam
   !!======================================================================
   !!                       ***  MODULE  lbclnk_tam  ***
   !! Ocean        : TAM of lateral boundary conditions
   !!=====================================================================
   !! History :  OPA  ! 1997-06  (G. Madec)     Original code
   !!   NEMO     1.0  ! 2002-09  (G. Madec)     F90: Free form and module
   !!            3.2  ! 2009-03  (R. Benshila)  External north fold treatment
   !! History of TAM : 3.2 ! ???
   !!                  3.4 ! 2012-03 (P.-A. Bouttier) Update
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!   lbc_lnk_adj      : generic interface for mpp_lnk_3d_adj and mpp_lnk_2d_adj routines defined in lib_mpp
   !!   lbc_lnk_e_adj    : generic interface for mpp_lnk_2d_e_adj routine defined in lib_mpp
   !!----------------------------------------------------------------------
   USE lib_mpp_tam          ! distributed memory computing library

   INTERFACE lbc_lnk_adj
      MODULE PROCEDURE mpp_lnk_3d_gather_adj, mpp_lnk_3d_adj, mpp_lnk_2d_adj
   END INTERFACE

   INTERFACE lbc_lnk_e_adj
      MODULE PROCEDURE mpp_lnk_2d_e_adj
   END INTERFACE

   PUBLIC lbc_lnk_adj       ! ocean lateral boundary conditions
   PUBLIC lbc_lnk_e_adj

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------


   !!======================================================================
END MODULE lbclnk_tam
