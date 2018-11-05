MODULE sol_oce_tam
   !!======================================================================
   !!                    ***  MODULE  sol_oce_tam  ***
   !! NEMOVAR : variables controlling the solver for tangent and adjoint models
   !!======================================================================
   !! History of the direct module:
   !!            1.0  ! 2002-11  (G. Madec)  F90: Free form and module
   !!            4.0  ! 2011-01  (A. R. Porter, STFC Daresbury) dynamical allocation
   !! History of the TAM module:
   !!            9.0  !  09-03  (A. Weaver) original version
   !!       NEMO 3.4  ! 2012-04 (P.-A. Bouttier) update
   !!       NEMO 3.4  ! 2012-09 (A. Vidard) Deallocating and initialising options
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sol_oce_alloc : allocate the solver arrays
   !!----------------------------------------------------------------------
   USE par_oce        ! ocean parameters
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distributed memory computing
   USE sol_oce
   USE solver

   IMPLICIT NONE
   PRIVATE

   PUBLIC                  &
      & sol_oce_alloc_tam,   &
      & sol_oce_dealloc_tam, &
      & sol_oce_tam_init

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gcx_tl, gcx_ad   !: TA of now solution of the elliptic eq.
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gcxb_tl, gcxb_ad !: TA of before solution of the elliptic eq.
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gcb_tl, gcb_ad   !: TA of second member of the elliptic eq.
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:)   ::   gcr_tl           !: Tangent of residu =b-a.x
   INTEGER,  PUBLIC, DIMENSION(:), SAVE, ALLOCATABLE :: nitsor          ! Number of SOR iterations
   LOGICAL, PRIVATE, SAVE :: ll_alloctl = .FALSE., ll_allocad = .FALSE.

   INTEGER, PARAMETER,  PUBLIC :: &
      & jp_it0adj = 300 ! initial value of the number of solver iteration (adj)

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence    (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS
   INTEGER FUNCTION sol_oce_alloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION sol_oce_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER, OPTIONAL :: kmode
      INTEGER  :: ierr(3)
      INTEGER :: imode
      !!----------------------------------------------------------------------
      IF ( PRESENT(kmode) ) THEN 
         imode = kmode
      ELSE 
         imode = 0
      END IF
      ierr(:) = 0
      !
      IF ( .NOT. ALLOCATED(nitsor) ) THEN
      ALLOCATE( nitsor(nitend - nit000 + 1) , STAT=ierr(1) )
      nitsor(:) = jp_it0adj
      END IF
      IF ( ( imode == 0 ) .OR. ( imode == 1 ) .AND. ( .NOT. ll_alloctl ) ) THEN 
      ALLOCATE( gcx_tl (1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)   ,     &
         &      gcxb_tl(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)   ,     &
         &      gcb_tl (1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)   ,     &
         &      gcr_tl (1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)   , STAT=ierr(2) )
         ll_alloctl = .TRUE.
      END IF
      IF ( ( imode == 0 ) .OR. ( imode == 2 ) .AND. ( .NOT. ll_allocad ) ) THEN 
      ALLOCATE( gcx_ad (1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)   ,     &
         &      gcxb_ad(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)   ,     &
         &      gcb_ad (1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj)   , STAT=ierr(3) )
         ll_allocad = .TRUE.
      END IF
         !
      sol_oce_alloc_tam = MAXVAL(ierr)
      !
      IF( lk_mpp            )   CALL mpp_sum ( sol_oce_alloc_tam )
      IF( sol_oce_alloc_tam > 0 )   CALL ctl_warn('sol_oce_alloc_tam: allocation of arrays failed')
      !
   END FUNCTION sol_oce_alloc_tam
   !
   INTEGER FUNCTION sol_oce_dealloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION sol_oce_dealloc  ***
      !!----------------------------------------------------------------------
      INTEGER, OPTIONAL :: kmode
      INTEGER  :: ierr(3)
      INTEGER :: imode
      !!----------------------------------------------------------------------
      IF ( PRESENT(kmode) ) THEN 
         imode = kmode
      ELSE 
         imode = 0
      END IF

      ierr(:) = 0
      !

      IF ( ( imode == 0 ) .OR. ( imode == 1 ) .AND. ( ll_alloctl ) ) THEN 
         DEALLOCATE( gcx_tl ,     &
              &      gcxb_tl,     &
              &      gcb_tl ,     &
              &      gcr_tl , STAT=ierr(2) )
         ll_alloctl = .FALSE.
      END IF
      IF ( ( imode == 0 ) .OR. ( imode == 2 ) .AND. ( ll_allocad ) ) THEN 
         DEALLOCATE( gcx_ad  ,     &
            &        gcxb_ad ,     &
            &        gcb_ad  , STAT=ierr(3) )
         ll_allocad = .FALSE.
      END IF
         !
      sol_oce_dealloc_tam = MAXVAL(ierr)
      !
      IF( lk_mpp                  )   CALL mpp_sum ( sol_oce_dealloc_tam )
      IF( sol_oce_dealloc_tam > 0 )   CALL ctl_warn('sol_oce_dealloc_tam: allocation of arrays failed')
      !
   END FUNCTION sol_oce_dealloc_tam
   !
   SUBROUTINE sol_oce_tam_init( kmode )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE sol_oce_tam_init  ***
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
      INTEGER, INTENT(IN) :: kmode
      INTEGER :: ierr
      IF ( ( kmode == 0 ) .OR. ( kmode == 1 ) ) THEN
         IF ( .NOT. ll_alloctl ) ierr = sol_oce_alloc_tam ( 1 )
         gcx_tl (:,:) = 0.0_wp
         gcxb_tl(:,:) = 0.0_wp
         gcb_tl (:,:) = 0.0_wp
         gcr_tl (:,:) = 0.0_wp
      END IF
      IF ( ( kmode == 0 ) .OR. ( kmode == 2 ) ) THEN
         IF ( .NOT. ll_allocad ) ierr = sol_oce_alloc_tam ( 2 )
         gcx_ad (:,:) = 0.0_wp
         gcxb_ad(:,:) = 0.0_wp
         gcb_ad (:,:) = 0.0_wp
      END IF
   END SUBROUTINE sol_oce_tam_init
   !!======================================================================

END MODULE sol_oce_tam
