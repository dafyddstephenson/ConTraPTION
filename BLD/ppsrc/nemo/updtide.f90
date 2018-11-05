MODULE updtide
  !!=================================================================================
  !!                       ***  MODULE  updtide  ***
  !! Initialization of tidal forcing
  !! History :  9.0  !  07  (O. Le Galloudec)  Original code
  !!=================================================================================
  !!----------------------------------------------------------------------
  !!   Dummy module :                                        NO TIDE
  !!----------------------------------------------------------------------
CONTAINS
  SUBROUTINE upd_tide( kt,kit )          ! Empty routine
    INTEGER,INTENT (IN) :: kt, kit
    WRITE(*,*) 'upd_tide: You should not have seen this print! error?', kt
  END SUBROUTINE upd_tide


  !!======================================================================

END MODULE updtide
