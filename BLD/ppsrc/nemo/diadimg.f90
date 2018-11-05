MODULE diadimg
   !!======================================================================
   !!                     ***  MODULE  diadimg  ***
   !! Ocean diagnostics :  write ocean output files in dimg direct access format (mpp)
   !!=====================================================================
   !!----------------------------------------------------------------------
   !!   Default option :                                       Empty module
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dia_wri_dimg( cd_name, cd_exper, ptab, klev, cd_type )
      REAL, DIMENSION(:,:,:) :: ptab
      INTEGER :: klev
      CHARACTER(LEN=80) :: cd_name, cd_exper,cd_type
      WRITE(*,*) ' This print must never occur ', cd_name, cd_exper,ptab, klev, cd_type
      WRITE(*,*) ' this routine is here just for compilation '
   END SUBROUTINE dia_wri_dimg
   !!======================================================================
END MODULE diadimg
