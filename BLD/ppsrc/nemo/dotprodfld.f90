MODULE dotprodfld
   !!======================================================================
   !!                       ***  MODULE dotprodfld ***
   !! NEMOVAR dotprodfld : Computes dot prodoct for 3D and 2D fields
   !!======================================================================
   !!
   !!----------------------------------------------------------------------
   !!     dot_product     : Computes the dot_product for two 3D/2D fields
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind
   USE dom_oce, ONLY :       &
      & nldi,                &
      & nldj,                &
      & nlei,                &
      & nlej
   USE par_oce       , ONLY: & ! Ocean space and time domain variables
      & jpi,                 &
      & jpj,                 &
      & jpk

   USE lib_fortran

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC &
      & dot_product

   !! * Interfaces

   INTERFACE dot_product
      MODULE PROCEDURE dot_product_3d
      MODULE PROCEDURE dot_product_2d
   END INTERFACE

CONTAINS

   FUNCTION dot_product_3d( pvec1, pvec2 )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE dot_product_3d  ***
      !!
      !! ** Purpose : Computes the dot_product for two 3D fields
      !!
      !! ** Method  : Use the mppsum module
      !!
      !! ** Action  :
      !!
      !! References :
      !!
      !! History :
      !!        !  07-08  (K. Mogensen)  Original code
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp) dot_product_3d
      !! * Arguments
      REAL(wp), INTENT(IN), DIMENSION(jpi,jpj,jpk) :: &
         & pvec1, &     ! 3D fields to compute dot_product of
         & pvec2
      !! * Local declarations

      dot_product_3d = glob_sum( &
         &                       PACK( pvec1(nldi:nlei,nldj:nlej,:),.TRUE.) * &
         &                       PACK( pvec2(nldi:nlei,nldj:nlej,:),.TRUE.),  &
         &                       (nlei-nldi+1) * (nlej-nldj+1) * jpk )

   END FUNCTION dot_product_3d

   FUNCTION dot_product_2d( pvec1, pvec2 )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE dot_product_2d  ***
      !!
      !! ** Purpose : Computes the dot_product for two 2D fields
      !!
      !! ** Method  : Use the mppsum module
      !!
      !! ** Action  :
      !!
      !! References :
      !!
      !! History :
      !!        !  07-08  (K. Mogensen)  Original code
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp) dot_product_2d
      !! * Arguments
      REAL(wp), INTENT(IN), DIMENSION(jpi,jpj) :: &
         & pvec1, &     ! 2D fields to compute dot_product of
         & pvec2

      dot_product_2d = glob_sum( &
         &                       PACK( pvec1(nldi:nlei,nldj:nlej),.TRUE.) * &
         &                       PACK( pvec2(nldi:nlei,nldj:nlej),.TRUE.),  &
         &                       (nlei-nldi+1) * (nlej-nldj+1) )

   END FUNCTION dot_product_2d

END MODULE dotprodfld
