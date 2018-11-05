MODULE ran_num
   !!======================================================================
   !!                       *** MODULE ran_num ***
   !! NEMOVAR: Random number routines
   !!======================================================================
   !!----------------------------------------------------------------------
   !! gaustb    : returns a gaussian randum number with mean and std.
   !! gaustb_2d : returns a gaussian randum number with mean and std for a
   !!             given horizontal grid point
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind       ! Kind variables
   USE dom_oce        ! Domain variables
   USE in_out_manager ! I/O stuff
   USE mt19937ar, ONLY : &
    & init_mtrand,       &
    & mtrand_real1,      &
    & mtrand_seeddump,   &
    & mtrand_seedread
 
   IMPLICIT NONE

   !! * Routine accessibility

   PRIVATE

   PUBLIC &
      & gaustb, &
      & gaustb_2d, &
      & random_seeddump, &
      & random_seedread

CONTAINS

   FUNCTION gaustb( pamp, pmean  )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE gaustb ***
      !!          
      !! ** Purpose : Returns a gaussian randum number with mean and std.
      !!
      !! ** Method  : Generate Gaussian random variables.
      !!              The standard deviation and mean of the variables are 
      !!              specified by the variables pamp and pmean.
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  07-07  (K. Mogensen)  Original code based on gaustb.F
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp) :: &
         & gaustb
      !! * Arguments
      REAL(wp), INTENT(IN) :: &
         & pamp, &      ! Amplitude
         & pmean        ! Mean value
      !! * Local declarations
      
      gaustb = pamp * gausva( ) + pmean

   END FUNCTION gaustb


   FUNCTION gausva( )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE gausva ***
      !!          
      !! ** Purpose : Returns a normally distributed deviate with 0 mean 
      !!              and unit variance using the unifva(kdum) as the 
      !!              source of uniform deviates.
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  07-07  (K. Mogensen)  Original code based on gaustb.F
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp) :: &
         & gausva
      !! * Local declarations
      REAL(wp), SAVE :: &
         & gset
      INTEGER, SAVE :: &
         & niset = 0
      REAL(wp) :: &
         & zfac, &
         & zrsq, &
         & zv1,  &
         & zv2

      ! Begin main

      IF ( niset == 0 ) THEN

         zv1   = 2.0_wp * psrandom( ) - 1.0_wp
         zv2   = 2.0_wp * psrandom( ) - 1.0_wp
         zrsq  = zv1**2 + zv2**2
         
         DO WHILE ( ( zrsq >= 1.0_wp ) .OR. ( zrsq == 0.0_wp ) ) 

            zv1   = 2.0_wp * psrandom( ) - 1.0_wp
            zv2   = 2.0_wp * psrandom( ) - 1.0_wp
            zrsq  = zv1**2 + zv2**2

         END DO

         zfac   = SQRT( -2.0_wp * LOG( zrsq ) / zrsq )
         gset   = zv1 * zfac
         gausva = zv2 * zfac
         niset  = 1

      ELSE

         gausva = gset
         niset  = 0

      ENDIF
       
   END FUNCTION gausva

   FUNCTION gaustb_2d( ki, kj, pamp, pmean  )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE gaustb_2d ***
      !!          
      !! ** Purpose : Returns a Gaussian randum number with mean and std
      !!              for a given horizontal grid point.
      !!
      !! ** Method  : Generate Gaussian random variables.
      !!              The standard deviation and mean of the variables are 
      !!              specified by the variables pamp and pmean.
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  07-07  (K. Mogensen)  Original code based on gaustb.F
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp) :: &
         & gaustb_2d
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & ki, &        ! Indices in seed array
         & kj
      REAL(wp), INTENT(IN) :: &
         & pamp, &      ! Amplitude
         & pmean        ! Mean value
      !! * Local declarations
      
      gaustb_2d = pamp * gausva_2d( ki, kj ) + pmean

   END FUNCTION gaustb_2d

   FUNCTION gausva_2d( ki, kj )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE gausva_2d ***
      !!          
      !! ** Purpose : Returns a normally distributed deviate with 0 mean 
      !!              and unit variance.
      !!
      !! ** Method  : 
      !!
      !! ** Action  : 
      !!
      !! History :
      !!        !  07-07  (K. Mogensen)  Original code based on gaustb.F
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp) :: &
         & gausva_2d
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & ki, &         ! Indices in seed array
         & kj
      !! * Local declarations
      REAL(wp), DIMENSION(jpi,jpj) :: &
         & gset
      INTEGER, DIMENSION(jpi,jpj) :: &
         & niset
      LOGICAL, SAVE :: &
         & llinit = .FALSE.
      REAL(wp) :: &
         & zfac, &
         & zrsq, &
         & zv1,  &
         & zv2

      ! Initialization

      IF ( .NOT. llinit ) THEN

         niset(:,:) = 0
         llinit     = .TRUE.

      ENDIF

      ! Begin main

      IF ( niset(ki,kj) == 0 ) THEN

         zv1   = 2.0_wp * psrandom( ) - 1.0_wp
         zv2   = 2.0_wp * psrandom( ) - 1.0_wp
         zrsq  = zv1**2 + zv2**2

         DO WHILE ( ( zrsq >= 1.0_wp ) .OR. ( zrsq == 0.0_wp ) ) 

            zv1   = 2.0_wp * psrandom( ) - 1.0_wp
            zv2   = 2.0_wp * psrandom( ) - 1.0_wp
            zrsq  = zv1**2 + zv2**2

         END DO

         zfac   = SQRT( -2.0_wp * LOG( zrsq ) / zrsq )
         gset(ki,kj) = zv1 * zfac
         gausva_2d   = zv2 * zfac
         niset(ki,kj) = 1

      ELSE

         gausva_2d = gset(ki,kj)
         niset(ki,kj)  = 0

      ENDIF
       
   END FUNCTION gausva_2d

   FUNCTION psrandom()
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE psrandom ***
      !!
      !! ** Purpose : Pseudo-Random number generator.
      !!
      !! ** Method  : Returns a pseudo-random number from a uniform distribution
      !!              between 0 and 1
      !!              Call with kdum a negative integer to initialize.
      !!              Thereafter, do not alter kdum between successive deviates
      !!              in sequence.
      !!
      !! ** Action  :
      !!
      !! History :
      !!        !  10-02  (F. Vigilant)  Original code
      !!        ! 2011-08 (D. Lea)  initialize with kdum1
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp) ::  &
         & psrandom
      !! * Arguments
      INTEGER  :: &
         & kdum          ! Seed
      LOGICAL, SAVE :: &
         & llinit = .TRUE.    ! Initialize on first call

      ! Initialization
      IF ( llinit ) THEN
         kdum = 456953 + nproc * 596035
         CALL init_mtrand(kdum)
         llinit     = .FALSE.    ! On subsequent calls do not reinitialize 
      ENDIF

      psrandom = mtrand_real1()

   END FUNCTION psrandom

   
   SUBROUTINE random_seeddump(cdfile)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE psrandom ***
      !!
      !! ** Purpose : Pseudo-Random number generator.
      !!
      !! ** Method  : Dump the seed
      !!
      !! ** Action  :
      !!
      !! History :
      !!        ! 2013-03 (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=128) :: cdfile
      CALL mtrand_seeddump(cdfile)

   END SUBROUTINE random_seeddump

   SUBROUTINE random_seedread(cdfile)
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE psrandom ***
      !!
      !! ** Purpose : Pseudo-Random number generator.
      !!
      !! ** Method  : Read the seed
      !!
      !! ** Action  :
      !!
      !! History :
      !!        ! 2013-03 (K. Mogensen) Original code
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(len=128) :: cdfile
      CALL mtrand_seedread(cdfile)

   END SUBROUTINE random_seedread

END MODULE ran_num
