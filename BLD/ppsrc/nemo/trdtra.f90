MODULE trdtra
   !!======================================================================
   !!                       ***  MODULE  trdtra  ***
   !! Ocean diagnostics:  ocean tracers trends
   !!=====================================================================
   !! History :  1.0  !  2004-08  (C. Talandier) Original code
   !!            2.0  !  2005-04  (C. Deltel)    Add Asselin trend in the ML budget
   !!            3.3  !  2010-06  (C. Ethe) merge TRA-TRC 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :          Dummy module           No trend diagnostics
   !!----------------------------------------------------------------------
   USE par_oce      ! ocean variables trends
CONTAINS
   SUBROUTINE trd_tra( kt, ctype, ktra, ktrd, ptrd, pu, ptra )
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in)           ::  kt      ! time step
      CHARACTER(len=3)                , INTENT(in)           ::  ctype   ! tracers trends type 'TRA'/'TRC'
      INTEGER                         , INTENT(in)           ::  ktra    ! tracer index
      INTEGER                         , INTENT(in)           ::  ktrd    ! tracer trend index
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in)           ::  ptrd    ! tracer trend 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  pu      ! velocity 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(in), OPTIONAL ::  ptra    ! Tracer variable 
      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd(1,1,1), ptra(1,1,1), pu(1,1,1),   &
         &                                                               ktrd, ktra, ctype, kt
   END SUBROUTINE trd_tra
   !!======================================================================
END MODULE trdtra
