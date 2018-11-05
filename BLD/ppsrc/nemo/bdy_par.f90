MODULE bdy_par
   !!======================================================================
   !!                      ***  MODULE bdy_par   ***
   !! Unstructured Open Boundary Cond. :   define related parameters
   !!======================================================================
   !! History :  1.0  !  2005-01  (J. Chanut, A. Sellar)  Original code
   !!            3.0  !  2008-04  (NEMO team)  add in the reference version
   !!            3.3  !  2010-09  (D. Storkey and E. O'Dea) update for Shelf configurations
   !!            3.4  !  2011     (D. Storkey) rewrite in preparation for OBC-BDY merge
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default option :            NO Unstructured open boundary condition
   !!----------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_bdy  = .FALSE.   !: Unstructured Ocean Boundary Condition flag

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: bdy_par.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE bdy_par
