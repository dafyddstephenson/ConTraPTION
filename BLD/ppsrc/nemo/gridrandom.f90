MODULE gridrandom
   !!======================================================================
   !!                       *** MODULE gridrandom ***
   !!
   !! NEMOVAR: Construct gridded random noise fields
   !!======================================================================
   !!----------------------------------------------------------------------
   !! grid_2d_ran     : Fill a 2D array with uncorrelated Gaussian noise
   !!                   using a constant seed
   !! grid_3d_ran     : Fill a 3D array with uncorrelated Gaussian noise
   !!                   using a constant seed
   !! grid_2d_ran_2d  : Fill a 2D array with uncorrelated Gaussian noise
   !!                   using a 2D seed array (for MPP)
   !! grid_3d_ran_2d  : Fill a 3D array with uncorrelated Gaussian noise
   !!                   using a 2D seed array (for MPP)
   !! grid_write_seed : Write out the 2D seed array for the random number 
   !!                   generator
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind       ! Kind variables
   USE dom_oce        ! Domain variables
   USE in_out_manager ! I/O stuff
   USE iom            ! I/O manager
   USE ran_num        ! Random number routines
   USE lbclnk         ! Boundary conditions and halos

   IMPLICIT NONE

   INTERFACE grid_random
      MODULE PROCEDURE grid_2d_ran, grid_3d_ran
   END INTERFACE

   INTERFACE grid_rd_sd
      MODULE PROCEDURE grid_2d_rd_sd_loc,  &
         &             grid_3d_rd_sd_loc
   END INTERFACE

   !! * Routine accessibility
   PRIVATE

   PUBLIC &
      & grid_random,      &
      & grid_rd_sd

CONTAINS

   SUBROUTINE grid_2d_ran( pt2d, cd_type, pmean, pstd )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE grid_2d_ran ***
      !!          
      !! ** Purpose : Fill a 2D (surface) array with uncorrelated Gaussian 
      !!              noise.
      !!
      !! ** Method  : The value of kseed is the seed for the random number
      !!              generator. On the first call to "grid_2d_ran" it should 
      !!              be set to a large negative number.
      !!
      !!              Apply the appropriate grid-point mask and lateral 
      !!              boundary conditions before exiting.
      !!
      !! ** Action  : 
      !!
      !! References : 
      !!
      !! History :
      !!        !  07-11  (A. Weaver)  
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      REAL, INTENT(INOUT), DIMENSION(jpi,jpj) :: &
         & pt2d        ! 2D field
      REAL(wp), INTENT(IN) :: &
         & pmean, &    ! Mean of noise
         & pstd        ! Standard deviation of noise

      !! * Local declarations
      CHARACTER(LEN=1), INTENT(IN) ::   &
         & cd_type     ! Nature of pt2d grid-point
                       !   = T , U  or V  grid-point
      INTEGER :: &
         & ji, &
         & jj

      !--------------------------------------------------------------------
      ! Fill in the 2D field with Gaussian noise
      !--------------------------------------------------------------------

      DO jj = 1, jpj
         DO ji = 1, jpi
            pt2d(ji,jj) = gaustb( pstd, pmean )
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Apply masks and lateral boundary conditions
      !--------------------------------------------------------------------

      SELECT CASE ( cd_type )
      CASE( 'T' )
 
         pt2d(:,:) = pt2d(:,:) * tmask(:,:,1)
         CALL lbc_lnk( pt2d, 'T',  1.0 )

      CASE( 'U' ) 

         pt2d(:,:) = pt2d(:,:) * umask(:,:,1)
         CALL lbc_lnk( pt2d, 'U', -1.0 )

      CASE ( 'V' )

         pt2d(:,:) = pt2d(:,:) * vmask(:,:,1)
         CALL lbc_lnk( pt2d, 'V', -1.0 )

      CASE ( 'S' ) !: AV: S ???
         CALL lbc_lnk( pt2d, 'S', 1.0 )

      CASE ( 'F' )         
      
         pt2d(:,:) = pt2d(:,:) * fmask(:,:,1)
         CALL lbc_lnk( pt2d, 'F',  1.0 )

      CASE ( 'W' )         
      
         pt2d(:,:) = pt2d(:,:) * tmask(:,:,1)
         CALL lbc_lnk( pt2d, 'W',  1.0 )

      END SELECT
            
   END SUBROUTINE grid_2d_ran
         
   SUBROUTINE grid_3d_ran( pt3d, cd_type, pmean, pstd )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE grid_3d_ran ***
      !!          
      !! ** Purpose : Fill a 3D array with uncorrelated Gaussian noise.
      !!
      !! ** Method  : The value of kseed is the seed for the random number
      !!              generator. On the first call to "grid_3d_ran" it should 
      !!              be set to a large negative number.
      !!
      !!              Apply the appropriate grid-point mask and lateral 
      !!              boundary conditions before exiting.
      !!
      !! ** Action  : 
      !!
      !! References : 
      !!
      !! History :
      !!        !  07-11  (A. Weaver)  
      !!----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      REAL, INTENT(INOUT), DIMENSION(jpi,jpj,jpk) :: &
         & pt3d        ! 3D field
      REAL(wp), INTENT(IN) :: &
         & pmean, &    ! Mean of noise
         & pstd        ! Standard deviation of noise
 
      !! * Local declarations
      CHARACTER(LEN=1), INTENT(IN) ::   &
         & cd_type     ! Nature of pt3d grid-point
                       !   = T , U  or V  grid-point
      INTEGER  :: &
         & ji, &
         & jj, &
         & jk

      !--------------------------------------------------------------------
      ! Fill in the 3D field with Gaussian noise
      !--------------------------------------------------------------------

      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               pt3d(ji,jj,jk) = gaustb( pstd, pmean )
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Apply masks and lateral boundary conditions
      !--------------------------------------------------------------------

      SELECT CASE ( cd_type )
      CASE( 'T' )

         pt3d(:,:,:) = pt3d(:,:,:) * tmask(:,:,:)
         CALL lbc_lnk( pt3d, 'T',  1.0 )
         
      CASE( 'U' )
            
         pt3d(:,:,:) = pt3d(:,:,:) * umask(:,:,:)
         CALL lbc_lnk( pt3d, 'U', -1.0 )
            
      CASE( 'V' )
         
         pt3d(:,:,:) = pt3d(:,:,:) * vmask(:,:,:)
         CALL lbc_lnk( pt3d, 'V', -1.0 )

      CASE( 'S' ) !: AV: S ???

         CALL lbc_lnk( pt3d, 'S', 1.0 )

      CASE( 'W' )

         pt3d(:,:,:) = pt3d(:,:,:) * tmask(:,:,:)
         CALL lbc_lnk( pt3d, 'W',  1.0 )

      CASE( 'F' )

         pt3d(:,:,:) = pt3d(:,:,:) * fmask(:,:,:)
         CALL lbc_lnk( pt3d, 'F',  1.0 )

      END SELECT

   END SUBROUTINE grid_3d_ran
      
   SUBROUTINE grid_2d_rd_sd_loc( pt2d, cd_type, pmean, pstd )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE grid_2d_rd_sd ***
      !!          
      !! ** Purpose : Fill a 2D (surface) array with uncorrelated Gaussian 
      !!              noise.
      !!
      !! ** Method  : The value of kseed is an integer from which a seed is 
      !!              generated for the random number
      !!              and then call grid_random routine
      !!
      !!              Apply the appropriate grid-point mask and lateral 
      !!              boundary conditions before exiting.
      !!
      !! ** Action  : 
      !!
      !! References : 
      !!
      !! History :
      !!        !  09-07  (F. Vigilant)  
      !!----------------------------------------------------------------------
      !! * Modules used
      USE par_oce       , ONLY: & ! Ocean space and time domain variables
         & jpiglo
      !! * Arguments
      REAL, INTENT(INOUT), DIMENSION(jpi,jpj) :: &
         & pt2d        ! 2D field
      REAL(wp), INTENT(IN) :: &
         & pmean, &    ! Mean of noise
         & pstd        ! Standard deviation of noise

      !! * Local declarations
      CHARACTER(LEN=1), INTENT(IN) ::   &
         & cd_type     ! Nature of pt2d grid-point
                       !   = T , U  or V  grid-point
      INTEGER :: &
         & ji, &
         & jj

      !--------------------------------------------------------------------
      ! Generate the noise
      !--------------------------------------------------------------------
      CALL grid_random( pt2d, cd_type, pmean, pstd )
            
   END SUBROUTINE grid_2d_rd_sd_loc

   SUBROUTINE grid_3d_rd_sd_loc( pt3d, cd_type, pmean, pstd )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE grid_3d_rd_sd ***
      !!          
      !! ** Purpose : Fill a 3D array with uncorrelated Gaussian 
      !!              noise.
      !!
      !! ** Method  : The value of kseed is an integer from which a seed is 
      !!              generated for the random number
      !!              and then call grid_random routine
      !!
      !!              Apply the appropriate grid-point mask and lateral 
      !!              boundary conditions before exiting.
      !!
      !! ** Action  : 
      !!
      !! References : 
      !!
      !! History :
      !!        !  09-07  (F. Vigilant)  
      !!----------------------------------------------------------------------
      !! * Modules used
      USE par_oce       , ONLY: & ! Ocean space and time domain variables
         & jpiglo
      !! * Arguments
      REAL, INTENT(INOUT), DIMENSION(jpi,jpj,jpk) :: &
         & pt3d        ! 3D field
      REAL(wp), INTENT(IN) :: &
         & pmean, &    ! Mean of noise
         & pstd        ! Standard deviation of noise

      !! * Local declarations
      CHARACTER(LEN=1), INTENT(IN) ::   &
         & cd_type     ! Nature of pt2d grid-point
                       !   = T , U  or V  grid-point
      INTEGER :: &
         & ji, &
         & jj


      !--------------------------------------------------------------------
      ! Generate the noise
      !--------------------------------------------------------------------
      CALL grid_random( pt3d, cd_type, pmean, pstd )
            
   END SUBROUTINE grid_3d_rd_sd_loc


END MODULE gridrandom
