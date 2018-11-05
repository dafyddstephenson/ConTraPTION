MODULE trc_oce_tam
   !!======================================================================
   !!                      ***  MODULE  trc_oce_tam  ***
   !! Ocean passive tracer  :  share SMS/Ocean variables
   !!                          Tangent and Adjoint Module
   !!======================================================================
   !! History of the direct module:
   !!   9.0  !  04-03  (C. Ethe)  F90: Free form and module
   !! History of the T&A module:
   !!   9.0  !  08-11  (A. Vidard)  original version
   !!       NEMO 3.4  ! 2012-09 (A. Vidard) Deallocating and initialising options
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! Default option                         No Biological fluxes for light
   !!----------------------------------------------------------------------
   USE par_oce
   USE dom_oce
   USE lib_mpp
   IMPLICIT NONE
   !! * Module variables
   REAL(wp), PUBLIC, SAVE, DIMENSION(:,:,:), ALLOCATABLE :: &
      & etot3_tl, &
      & etot3_ad
   LOGICAL, PRIVATE, SAVE :: ll_alloctl = .FALSE., ll_allocad = .FALSE.
   PUBLIC &
      & trc_oce_alloc_tam,   &       !: Initialize the trend TAM fields
      & trc_oce_dealloc_tam, &       !: Deallocate the trend TAM fields
      & trc_oce_tam_init
CONTAINS

   INTEGER FUNCTION trc_oce_alloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                  ***  trc_oce_alloc_tam  ***
      !!----------------------------------------------------------------------
      INTEGER, optional :: kmode
      INTEGER ::   ierr(2)        ! Local variables
      INTEGER :: imode
      !!----------------------------------------------------------------------
      IF ( PRESENT(kmode) ) THEN
         imode = kmode
      ELSE
         imode = 0
      END IF
      ierr(:) = 0
      IF ( ( imode == 0 ) .OR. ( imode == 1 ) .AND. ( .NOT. ll_alloctl ) ) THEN
      ALLOCATE( etot3_tl (jpi,jpj,jpk), STAT=ierr(1) )
      ll_alloctl = .TRUE.
      END IF
      IF ( ( imode == 0 ) .OR. ( imode == 2 ) .AND. ( .NOT. ll_allocad ) ) THEN
      ALLOCATE( etot3_ad (jpi,jpj,jpk), STAT=ierr(2) )
      ll_allocad = .TRUE.
      END IF
      trc_oce_alloc_tam  = MAXVAL( ierr )
      !
      IF( trc_oce_alloc_tam /= 0 )   CALL ctl_warn('trc_oce_alloc: failed to allocate etot3 array')
   END FUNCTION trc_oce_alloc_tam
   !
   INTEGER FUNCTION trc_oce_dealloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                  ***  trc_oce_dealloc_tam  ***
      !!----------------------------------------------------------------------
      INTEGER, optional :: kmode
      INTEGER :: imode
      INTEGER ::   ierr(2)        ! Local variables
      !!----------------------------------------------------------------------
      IF ( PRESENT(kmode) ) THEN
         imode = kmode
      ELSE
         imode = 0
      END IF
      ierr(:) = 0
      IF ( ( imode == 0 ) .OR. ( imode == 1 ) .AND. ( ll_alloctl ) ) THEN
         DEALLOCATE( etot3_tl, STAT=ierr(1) )
         ll_alloctl = .FALSE.
      END IF
      IF ( ( imode == 0 ) .OR. ( imode == 2 ) .AND. ( ll_allocad ) ) THEN
         DEALLOCATE( etot3_ad, STAT=ierr(2) )
         ll_allocad = .FALSE.
      END IF
      trc_oce_dealloc_tam  = MAXVAL( ierr )
      !
      IF( trc_oce_dealloc_tam /= 0 )   CALL ctl_warn('trc_oce_dealloc: failed to deallocate etot3 array')
   END FUNCTION trc_oce_dealloc_tam
   !
   SUBROUTINE trc_oce_tam_init( kmode )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE trc_oce_tam_init  ***
      !!
      !! ** Purpose : Allocate and initialize the tangent linear and
      !!              adjoint arrays
      !!
      !! ** Method  : kindic = 0  allocate/initialize both tl and ad variables
      !!              kindic = 1  allocate/initialize only tl variables
      !!              kindic = 2  allocate/initialize only ad variables
      !!
      !! ** Action  :
      !!
      !! References :
      !!
      !! History :
      !!        ! 2009-03 (A. Weaver) Initial version (based on oce_tam_init)
      !!        ! 2010-04 (A. Vidard) Nemo3.2 update
      !!        ! 2012-09 (A. Vidard) Nemo3.4 update
      !!-----------------------------------------------------------------------
      !! * Arguments
      INTEGER :: kmode
      INTEGER :: ierr
      IF ( ( kmode == 0 ) .OR. ( kmode == 1 ) ) THEN
         IF ( .NOT. ll_alloctl ) ierr = trc_oce_alloc_tam ( 1 )
         etot3_tl(:,:,:) = 0.0_wp
      END IF
      IF ( ( kmode == 0 ) .OR. ( kmode == 2 ) ) THEN
         IF ( .NOT. ll_allocad ) ierr = trc_oce_alloc_tam ( 2 )
         etot3_ad(:,:,:) = 0.0_wp
      END IF
   END SUBROUTINE trc_oce_tam_init
END MODULE trc_oce_tam

