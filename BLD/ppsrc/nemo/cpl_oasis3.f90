MODULE cpl_oasis3
   !!======================================================================
   !!                    ***  MODULE cpl_oasis  ***
   !! Coupled O/A : coupled ocean-atmosphere case using OASIS3 V. prism_2_4
   !!               special case: NEMO OPA/LIM coupled to ECHAM5
   !!=====================================================================
   !! History :   
   !!   9.0  !  04-06  (R. Redler, NEC Laboratories Europe, Germany) Original code
   !!   " "  !  04-11  (R. Redler, NEC Laboratories Europe; N. Keenlyside, W. Park, IFM-GEOMAR, Germany) revision
   !!   " "  !  04-11  (V. Gayler, MPI M&D) Grid writing
   !!   " "  !  05-08  (R. Redler, W. Park) frld initialization, paral(2) revision
   !!   " "  !  05-09  (R. Redler) extended to allow for communication over root only
   !!   " "  !  06-01  (W. Park) modification of physical part
   !!   " "  !  06-02  (R. Redler, W. Park) buffer array fix for root exchange
   !!   3.4  !  11-11  (C. Harris) Changes to allow mutiple category fields
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case          Dummy module          Forced Ocean/Atmosphere
   !!----------------------------------------------------------------------
   USE in_out_manager               ! I/O manager
   LOGICAL, PUBLIC, PARAMETER :: lk_cpl = .FALSE.   !: coupled flag
   PUBLIC cpl_prism_init
   PUBLIC cpl_prism_finalize
CONTAINS
   SUBROUTINE cpl_prism_init (kl_comm) 
      INTEGER, INTENT(out)   :: kl_comm       ! local communicator of the model
      kl_comm = -1
      WRITE(numout,*) 'cpl_prism_init: Error you sould not be there...'
   END SUBROUTINE cpl_prism_init
   SUBROUTINE cpl_prism_finalize
      WRITE(numout,*) 'cpl_prism_finalize: Error you sould not be there...'
   END SUBROUTINE cpl_prism_finalize

   !!=====================================================================
END MODULE cpl_oasis3
