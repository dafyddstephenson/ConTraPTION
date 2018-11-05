MODULE zdf_oce_tam
   !!======================================================================
   !!              ***  MODULE  zdf_oce  ***
   !! Ocean physics : define vertical mixing variables
   !!=====================================================================
   !! history :  1.0  !  2002-06  (G. Madec) Original code
   !!            3.2  !  2009-07  (G.Madec) addition of avm
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! MPP library

   IMPLICIT NONE
   PRIVATE

   PUBLIC  zdf_oce_alloc_tam    ! Called in nemogcm.F90

   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   bfrua_ad, bfrva_ad   !: Bottom friction coefficients set in zdfbfr
   REAL(wp), PUBLIC, SAVE, ALLOCATABLE, DIMENSION(:,:)   ::   bfrua_tl, bfrva_tl   !: Bottom friction coefficients set in zdfbfr
   LOGICAL, PRIVATE, SAVE :: ll_alloctl = .FALSE.
   LOGICAL, PRIVATE, SAVE :: ll_allocad = .FALSE.
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION zdf_oce_alloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!            *** FUNCTION zdf_oce_alloc ***
      !!----------------------------------------------------------------------
      !
      INTEGER, OPTIONAL :: kmode
      INTEGER :: jmode
      INTEGER, DIMENSION(2) :: ierr

      IF ( PRESENT( kmode ) ) THEN
         jmode = kmode
      ELSE
         jmode = 0
      END IF


      IF ( ( jmode == 0 ) .OR. ( jmode == 1 ) .AND. ( .NOT. ll_alloctl ) ) THEN
         ALLOCATE(bfrua_tl(jpi,jpj) ,    &
            &     bfrva_tl(jpi,jpj) , STAT =  ierr(1) )
         ll_alloctl = .TRUE.
      END IF
      IF ( ( jmode == 0 ) .OR. ( jmode == 2 ) .AND. ( .NOT. ll_allocad ) ) THEN
         ALLOCATE(bfrua_ad(jpi,jpj) ,    &
            &     bfrva_ad(jpi,jpj) , STAT =  ierr(2) )
         ll_allocad = .TRUE.
      END IF
      zdf_oce_alloc_tam = SUM( ierr )
         !
      IF( lk_mpp                 )   CALL mpp_sum( zdf_oce_alloc_tam ) 
      IF( zdf_oce_alloc_tam /= 0 )   CALL ctl_warn('zdf_oce_alloc_tam: failed to allocate arrays')
      !
   END FUNCTION zdf_oce_alloc_tam

   INTEGER FUNCTION zdf_oce_dealloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!            *** FUNCTION zdf_oce_alloc ***
      !!----------------------------------------------------------------------
      !
      INTEGER, OPTIONAL :: kmode
      INTEGER :: jmode
      INTEGER, DIMENSION(2) :: ierr

      IF ( PRESENT( kmode ) ) THEN
         jmode = kmode
      ELSE
         jmode = 0
      END IF


      IF ( ( jmode == 0 ) .OR. ( jmode == 1 ) .AND. ( ll_alloctl) ) THEN
         DEALLOCATE(bfrua_tl, bfrva_tl, STAT =  ierr(1) )
         ll_alloctl = .FALSE.
      END IF
      IF ( ( jmode == 0 ) .OR. ( jmode == 2 ) .AND. ( ll_allocad) ) THEN
         DEALLOCATE(bfrua_ad, bfrva_ad, STAT =  ierr(2) )
         ll_allocad = .FALSE.
      END IF
      zdf_oce_dealloc_tam = SUM( ierr )
         !
      IF( lk_mpp                 )   CALL mpp_sum( zdf_oce_dealloc_tam ) 
      IF( zdf_oce_dealloc_tam /= 0 )   CALL ctl_warn('zdf_oce_dealloc_tam: failed to deallocate arrays')
      !
   END FUNCTION zdf_oce_dealloc_tam

   !!======================================================================
END MODULE zdf_oce_tam
