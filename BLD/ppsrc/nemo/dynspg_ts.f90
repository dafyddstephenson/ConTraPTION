MODULE dynspg_ts
   !!======================================================================
   !! History :   1.0  ! 2004-12  (L. Bessieres, G. Madec)  Original code
   !!              -   ! 2005-11  (V. Garnier, G. Madec)  optimization
   !!              -   ! 2006-08  (S. Masson)  distributed restart using iom
   !!             2.0  ! 2007-07  (D. Storkey) calls to BDY routines
   !!              -   ! 2008-01  (R. Benshila)  change averaging method
   !!             3.2  ! 2009-07  (R. Benshila, G. Madec) Complete revisit associated to vvl reactivation
   !!             3.3  ! 2010-09  (D. Storkey, E. O'Dea) update for BDY for Shelf configurations
   !!             3.3  ! 2011-03  (R. Benshila, R. Hordoir, P. Oddo) update calculation of ub_b
   !!---------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :   Empty module   No standart free surface cst volume
   !!----------------------------------------------------------------------
CONTAINS
   INTEGER FUNCTION dyn_spg_ts_alloc()    ! Dummy function
      dyn_spg_ts_alloc = 0
   END FUNCTION dyn_spg_ts_alloc
   SUBROUTINE dyn_spg_ts( kt )            ! Empty routine
      INTEGER, INTENT(in) :: kt
      WRITE(*,*) 'dyn_spg_ts: You should not have seen this print! error?', kt
   END SUBROUTINE dyn_spg_ts
   SUBROUTINE ts_rst( kt, cdrw )          ! Empty routine
      INTEGER         , INTENT(in) ::   kt         ! ocean time-step
      CHARACTER(len=*), INTENT(in) ::   cdrw       ! "READ"/"WRITE" flag
      WRITE(*,*) 'ts_rst    : You should not have seen this print! error?', kt, cdrw
   END SUBROUTINE ts_rst    
   
   !!======================================================================
END MODULE dynspg_ts
