MODULE obcfla
   !!======================================================================
   !!                       ***  MODULE  obcfla  ***
   !! Ocean dynamics:   Flather's algorithm at open boundaries for the time-splitting
   !!======================================================================
   !! History :  2.0  ! 2005-12  (V. Garnier) original code
   !!            3.3  ! 2010-11  (G. Madec) 
   !!            4.0  ! 2011-02  (G. Madec) velocity & ssh passed in argument
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module :                             No OBC or time-splitting
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE obc_fla_ts( pua, pva, p_sshn, p_ssha )
      REAL, DIMENSION(:,:)::   pua, pva, p_sshn, p_ssha
      WRITE(*,*) 'obc_fla_ts: You should not have seen this print! error?', pua(1,1), pva(1,1), p_sshn(1,1), p_ssha(1,1)
   END SUBROUTINE obc_fla_ts
   !!======================================================================
END MODULE obcfla
