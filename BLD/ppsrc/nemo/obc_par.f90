MODULE obc_par
   !!==============================================================================
   !!                  ***  MODULE obc_par   ***
   !! Open Boundary Cond. :   define related parameters
   !!==============================================================================
   !! history :  OPA  ! 1991-01 (CLIPPER)  Original code 
   !!   NEMO     1.0  ! 2002-04   (C. Talandier)  modules
   !!             -   ! 2004/06   (F. Durand) jptobc is defined as a parameter
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :                         NO open boundary condition
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_obc = .FALSE.     !: Ocean Boundary Condition flag

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obc_par.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE obc_par
