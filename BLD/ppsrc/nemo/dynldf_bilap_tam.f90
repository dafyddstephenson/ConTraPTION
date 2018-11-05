MODULE dynldf_bilap_tam
   !!===========================================================================
   !!                       ***  MODULE  dynldf_bilap_tam  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!                  Tangent and Adjoint Module
   !!===========================================================================

   !!---------------------------------------------------------------------------
   !!   dyn_ldf_bilap_tan : update the momentum trend with the lateral diffusion
   !!                       using an iso-level bilaplacian operator (tangent)
   !!   dyn_ldf_bilap_adj : update the momentum trend with the lateral diffusion
   !!                       using an iso-level bilaplacian operator (adjoint)
   !!---------------------------------------------------------------------------
   !! * Modules used
   USE lbclnk
   USE lbclnk_tam
   USE par_oce
   USE oce_tam
   USE ldfdyn_oce
   USE dom_oce
   USE in_out_manager
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC dyn_ldf_bilap_tan  ! called by dynldf_tam.F90
   PUBLIC dyn_ldf_bilap_adj  ! called by dynldf_tam.F90

   !! * Substitutions
   !!----------------------------------------------------------------------
   !!                    ***  domzgr_substitute.h90   ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsdep. and fse.., the vert. depth and scale
   !!      factors depending on the vertical coord. used, using CPP macro.
   !!----------------------------------------------------------------------
   !! History :  1.0  !  2005-10  (A. Beckmann, G. Madec) generalisation to all coord.
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  pure z* coordinate
   !!----------------------------------------------------------------------
! reference for s- or zps-coordinate (3D no time dependency)
! z- or s-coordinate (1D or 3D + no time dependency) use reference in all cases




   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: domzgr_substitute.h90 2528 2010-12-27 17:33:53Z rblod $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!                    ***  ldfdyn_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsahm., the lateral eddy viscosity coeff. 
   !!      with a constant, or 1D, or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfdyn_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!
   !! fsahmt, fsahmf - used for laplaian operator only
   !! fsahmu, fsahmv - used for bilaplacian operator only
   !!
!   ' key_dynldf_c3d' :                  3D coefficient
   !!----------------------------------------------------------------------
   !!                   ***  vectopt_loop_substitute  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute the inner loop starting and inding indices 
   !!      to allow unrolling of do-loop using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: vectopt_loop_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_ldf_bilap_tan( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_bilap_tan  ***
      !!
      !! ** Purpose :   Compute the before trend of the lateral momentum
      !!      diffusion and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The before horizontal momentum diffusion trend is a
      !!      bi-harmonic operator (bilaplacian type) which separates the
      !!      divergent and rotational parts of the flow.
      !!      Its horizontal components are computed as follow:
      !!      laplacian:
      !!          zlu = 1/e1u di[ hdivb ] - 1/(e2u*e3u) dj-1[ e3f rotb ]
      !!          zlv = 1/e2v dj[ hdivb ] + 1/(e1v*e3v) di-1[ e3f rotb ]
      !!      third derivative:
      !!       * multiply by the eddy viscosity coef. at u-, v-point, resp.
      !!          zlu = ahmu * zlu
      !!          zlv = ahmv * zlv
      !!       * curl and divergence of the laplacian
      !!          zuf = 1/(e1f*e2f) ( di[e2v zlv] - dj[e1u zlu] )
      !!          zut = 1/(e1t*e2t*e3t) ( di[e2u*e3u zlu] + dj[e1v*e3v zlv] )
      !!      bilaplacian:
      !!              diffu = 1/e1u di[ zut ] - 1/(e2u*e3u) dj-1[ e3f zuf ]
      !!              diffv = 1/e2v dj[ zut ] + 1/(e1v*e3v) di-1[ e3f zuf ]
      !!      If ln_sco=F and ln_zps=F, the vertical scale factors in the
      !!      rotational part of the diffusion are simplified
      !!      Add this before trend to the general trend (ua,va):
      !!            (ua,va) = (ua,va) + (diffu,diffv)
      !!      'key_trddyn' defined: the two components of the horizontal
      !!                               diffusion trend are saved.
      !!
      !! ** Action : - Update (ua,va) with the before iso-level biharmonic
      !!               mixing trend.
      !!
      !! History :
      !!        !  90-09  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-07  (G. Madec)  lbc calls
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !! History of the tangent routine
      !!   9.0  !  09-12 (F. Vigilant) tangent of 9.0
      !!   3.4  !  11-07 (P.-A. Bouttier) 3.4 version
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt           ! ocean time-step index
      !! * Local declarations
      INTEGER  ::   ji, jj, jk                ! dummy loop indices
      REAL(wp) ::   &
         zuatl, zvatl, zbt, ze2u, ze2v        ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::       &
         zcutl, zcvtl                         ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   &
         zuftl, zuttl, zlutl, zlvtl           ! temporary workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_bilap_tan')
      !
      CALL wrk_alloc( jpi, jpj,      zcutl, zcvtl           )
      CALL wrk_alloc( jpi, jpj, jpk, zuftl, zuttl, zlutl, zlvtl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_bilap_tan: iso-level bilaplacien operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~ '
      ENDIF

      zuftl(:,:,:) = 0.0_wp
      zuttl(:,:,:) = 0.0_wp
      zlutl(:,:,:) = 0.0_wp
      zlvtl(:,:,:) = 0.0_wp
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Laplacian
         ! ---------

         IF( ln_sco .OR. ln_zps ) THEN   ! s-coordinate or z-coordinate with partial steps
            zuftl(:,:,jk) = rotb_tl(:,:,jk) * e3f(:,:,jk)
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zlutl(ji,jj,jk) = - ( zuftl(ji,jj,jk) - zuftl(ji,jj-1,jk) ) / ( e2u(ji,jj) * e3u(ji,jj,jk) )   &
                     &         + ( hdivb_tl(ji+1,jj,jk) - hdivb_tl(ji,jj,jk) ) / e1u(ji,jj)

                  zlvtl(ji,jj,jk) = + ( zuftl(ji,jj,jk) - zuftl(ji-1,jj,jk) ) / ( e1v(ji,jj) * e3v(ji,jj,jk) )   &
                     &         + ( hdivb_tl(ji,jj+1,jk) - hdivb_tl(ji,jj,jk) ) / e2v(ji,jj)
               END DO
            END DO
         ELSE                            ! z-coordinate - full step
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zlutl(ji,jj,jk) = - ( rotb_tl (ji  ,jj,jk) - rotb_tl (ji,jj-1,jk) ) / e2u(ji,jj)   &
                     &         + ( hdivb_tl(ji+1,jj,jk) - hdivb_tl(ji,jj  ,jk) ) / e1u(ji,jj)

                  zlvtl(ji,jj,jk) = + ( rotb_tl (ji,jj  ,jk) - rotb_tl (ji-1,jj,jk) ) / e1v(ji,jj)   &
                     &         + ( hdivb_tl(ji,jj+1,jk) - hdivb_tl(ji  ,jj,jk) ) / e2v(ji,jj)
               END DO
            END DO
         ENDIF
      ENDDO

      ! Boundary conditions on the laplacian  (zlu,zlv)
      CALL lbc_lnk( zlutl, 'U', -1.0_wp )
      CALL lbc_lnk( zlvtl, 'V', -1.0_wp )

      DO jk = 1, jpkm1

         ! Third derivative
         ! ----------------

         ! Multiply by the eddy viscosity coef. (at u- and v-points)
         zlutl(:,:,jk) = zlutl(:,:,jk) * ahm3(:,:,jk)
         zlvtl(:,:,jk) = zlvtl(:,:,jk) * ahm4(:,:,jk)

         ! Contravariant "laplacian"
         zcutl(:,:) = e1u(:,:) * zlutl(:,:,jk)
         zcvtl(:,:) = e2v(:,:) * zlvtl(:,:,jk)

         ! Laplacian curl ( * e3f if s-coordinates or z-coordinate with partial steps)
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zuftl(ji,jj,jk) = fmask(ji,jj,jk) * (  zcvtl(ji+1,jj  ) - zcvtl(ji,jj)      &
                  &                            - zcutl(ji  ,jj+1) + zcutl(ji,jj)  )   &

                  &       * e3f(ji,jj,jk) / ( e1f(ji,jj)*e2f(ji,jj) )
            END DO
         END DO

         ! Laplacian Horizontal fluxes
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zlutl(ji,jj,jk) = e2u(ji,jj) * e3u(ji,jj,jk) * zlutl(ji,jj,jk)
               zlvtl(ji,jj,jk) = e1v(ji,jj) * e3v(ji,jj,jk) * zlvtl(ji,jj,jk)
            END DO
         END DO

         ! Laplacian divergence
         DO jj = 2, jpj
            DO ji = 2, jpi   ! vector opt.
               zbt = e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)
               zuttl(ji,jj,jk) = (  zlutl(ji,jj,jk) - zlutl(ji-1,jj  ,jk)   &
                  &             + zlvtl(ji,jj,jk) - zlvtl(ji  ,jj-1,jk) ) / zbt
            END DO
         END DO
      END DO

      ! boundary conditions on the laplacian curl and div (zuf,zut)
!!bug gm no need to do this 2 following lbc...
      CALL lbc_lnk( zuftl, 'F', 1.0_wp )
      CALL lbc_lnk( zuttl, 'T', 1.0_wp )

      DO jk = 1, jpkm1

         ! Bilaplacian
         ! -----------

         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ze2u = e2u(ji,jj) * e3u(ji,jj,jk)
               ze2v = e1v(ji,jj) * e3v(ji,jj,jk)
               ! horizontal biharmonic diffusive trends
               zuatl = - ( zuftl(ji  ,jj,jk) - zuftl(ji,jj-1,jk) ) / ze2u   &
                  &  + ( zuttl(ji+1,jj,jk) - zuttl(ji,jj  ,jk) ) / e1u(ji,jj)
               zvatl = + ( zuftl(ji,jj  ,jk) - zuftl(ji-1,jj,jk) ) / ze2v   &
                  &  + ( zuttl(ji,jj+1,jk) - zuttl(ji  ,jj,jk) ) / e2v(ji,jj)
               ! add it to the general momentum trends
               ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) + zuatl
               va_tl(ji,jj,jk) = va_tl(ji,jj,jk) + zvatl
            END DO
         END DO

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      CALL wrk_dealloc( jpi, jpj,      zcutl, zcvtl               )
      CALL wrk_dealloc( jpi, jpj, jpk, zuftl, zuttl, zlutl, zlvtl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_bilap_tan')
      !
   END SUBROUTINE dyn_ldf_bilap_tan


   SUBROUTINE dyn_ldf_bilap_adj( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_bilap_adj  ***
      !!
      !! ** Purpose :   Compute the before trend of the lateral momentum
      !!      diffusion and add it to the general trend of momentum equation.
      !!
      !! ** Method  :   The before horizontal momentum diffusion trend is a
      !!      bi-harmonic operator (bilaplacian type) which separates the
      !!      divergent and rotational parts of the flow.
      !!      Its horizontal components are computed as follow:
      !!      laplacian:
      !!          zlu = 1/e1u di[ hdivb ] - 1/(e2u*e3u) dj-1[ e3f rotb ]
      !!          zlv = 1/e2v dj[ hdivb ] + 1/(e1v*e3v) di-1[ e3f rotb ]
      !!      third derivative:
      !!       * multiply by the eddy viscosity coef. at u-, v-point, resp.
      !!          zlu = ahmu * zlu
      !!          zlv = ahmv * zlv
      !!       * curl and divergence of the laplacian
      !!          zuf = 1/(e1f*e2f) ( di[e2v zlv] - dj[e1u zlu] )
      !!          zut = 1/(e1t*e2t*e3t) ( di[e2u*e3u zlu] + dj[e1v*e3v zlv] )
      !!      bilaplacian:
      !!              diffu = 1/e1u di[ zut ] - 1/(e2u*e3u) dj-1[ e3f zuf ]
      !!              diffv = 1/e2v dj[ zut ] + 1/(e1v*e3v) di-1[ e3f zuf ]
      !!      If ln_sco=F and ln_zps=F, the vertical scale factors in the
      !!      rotational part of the diffusion are simplified
      !!      Add this before trend to the general trend (ua,va):
      !!            (ua,va) = (ua,va) + (diffu,diffv)
      !!      'key_trddyn' defined: the two components of the horizontal
      !!                               diffusion trend are saved.
      !!
      !! ** Action : - Update (ua,va) with the before iso-level biharmonic
      !!               mixing trend.
      !!
      !! History :
      !!        !  90-09  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!        !  93-03  (M. Guyon)  symetrical conditions (M. Guyon)
      !!        !  96-01  (G. Madec)  statement function for e3
      !!        !  97-07  (G. Madec)  lbc calls
      !!   8.5  !  02-08  (G. Madec)  F90: Free form and module
      !!   9.0  !  04-08  (C. Talandier) New trends organization
      !! History of the adjoint routine
      !!   9.0  !  09-12 (F. Vigilant) adjoint of 9.0
      !!   3.4  !  12-07  (P.-A. Bouttier) 3.4 version
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      !! * Local declarations
      INTEGER  ::   ji, jj, jk            ! dummy loop indices
      REAL(wp) ::   &
         zuaad, zvaad, zbt, ze2u, ze2v        ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::       &
         zcuad, zcvad                         ! temporary workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   &
         zufad, zutad, zluad, zlvad           ! temporary workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_bilap_adj')
      !
      CALL wrk_alloc( jpi, jpj,      zcuad, zcvad           )
      CALL wrk_alloc( jpi, jpj, jpk, zufad, zutad, zluad, zlvad )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf_bilap_adj: bilaplacien operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~~~ '
      ENDIF

      zuaad = 0.0_wp
      zvaad = 0.0_wp

      zufad(:,:,:) = 0.0_wp
      zutad(:,:,:) = 0.0_wp
      zluad(:,:,:) = 0.0_wp
      zlvad(:,:,:) = 0.0_wp

      zcvad(:,:) = 0.0_wp
      zcuad(:,:) = 0.0_wp

      DO jk = 1, jpkm1

         ! Bilaplacian
         ! -----------

         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               ze2u = e2u(ji,jj) * e3u(ji,jj,jk)
               ze2v = e1v(ji,jj) * e3v(ji,jj,jk)
               ! add it to the general momentum trends
               zvaad = zvaad + va_ad(ji,jj,jk)
               zuaad = zuaad + ua_ad(ji,jj,jk)
               ! horizontal biharmonic diffusive trends
               zufad(ji  ,jj  ,jk) = zufad(ji  ,jj  ,jk) + zvaad / ze2v
               zufad(ji-1,jj  ,jk) = zufad(ji-1,jj  ,jk) - zvaad / ze2v
               zutad(ji  ,jj  ,jk) = zutad(ji  ,jj  ,jk) - zvaad / e2v(ji,jj)
               zutad(ji  ,jj+1,jk) = zutad(ji  ,jj+1,jk) + zvaad / e2v(ji,jj)

               zufad(ji  ,jj  ,jk) = zufad(ji  ,jj  ,jk) - zuaad / ze2u
               zufad(ji  ,jj-1,jk) = zufad(ji  ,jj-1,jk) + zuaad / ze2u
               zutad(ji  ,jj  ,jk) = zutad(ji  ,jj  ,jk) - zuaad / e1u(ji,jj)
               zutad(ji+1,jj  ,jk) = zutad(ji+1,jj  ,jk) + zuaad / e1u(ji,jj)
               
               zuaad = 0.0_wp
               zvaad = 0.0_wp
            END DO
         END DO

         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      ! boundary conditions on the laplacian curl and div (zuf,zut)
!!bug gm no need to do this 2 following lbc...
      CALL lbc_lnk_adj( zutad, 'T', 1.0_wp )
      CALL lbc_lnk_adj( zufad, 'F', 1.0_wp )

      DO jk = 1, jpkm1

         ! Third derivative
         ! ----------------

         ! Laplacian divergence
         DO jj = jpj, 2, -1
            DO ji = jpi, 2, -1   ! vector opt.
               zbt = e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk)
               zluad(ji  ,jj  ,jk) = zluad(ji  ,jj  ,jk) + zutad(ji,jj,jk) / zbt
               zluad(ji-1,jj  ,jk) = zluad(ji-1,jj  ,jk) - zutad(ji,jj,jk) / zbt
               zlvad(ji  ,jj  ,jk) = zlvad(ji  ,jj  ,jk) + zutad(ji,jj,jk) / zbt
               zlvad(ji  ,jj-1,jk) = zlvad(ji  ,jj-1,jk) - zutad(ji,jj,jk) / zbt

               zutad(ji,jj,jk) = 0.0_wp
            END DO
         END DO

         ! Laplacian Horizontal fluxes
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1   ! vector opt.
               zluad(ji,jj,jk) = e2u(ji,jj) * e3u(ji,jj,jk) * zluad(ji,jj,jk)
               zlvad(ji,jj,jk) = e1v(ji,jj) * e3v(ji,jj,jk) * zlvad(ji,jj,jk)
            END DO
         END DO

         ! Laplacian curl ( * e3f if s-coordinates or z-coordinate with partial steps)
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1   ! vector opt.
               zufad(ji,jj,jk) = fmask(ji,jj,jk) * zufad(ji,jj,jk) * e3f(ji,jj,jk) / ( e1f(ji,jj)*e2f(ji,jj) )
               zcvad(ji  ,jj  ) = zcvad(ji  ,jj  ) - zufad(ji,jj,jk)
               zcvad(ji+1,jj  ) = zcvad(ji+1,jj  ) + zufad(ji,jj,jk)
               zcuad(ji  ,jj  ) = zcuad(ji  ,jj  ) + zufad(ji,jj,jk)
               zcuad(ji  ,jj+1) = zcuad(ji  ,jj+1) - zufad(ji,jj,jk)

               zufad(ji,jj,jk) = 0.0_wp
            END DO
         END DO

         ! Contravariant "laplacian"
         DO jj = 1, jpj
            DO ji = 1, jpi
               zlvad(ji,jj,jk) = zlvad(ji,jj,jk) + e2v(ji,jj) * zcvad(ji,jj)
               zluad(ji,jj,jk) = zluad(ji,jj,jk) + e1u(ji,jj) * zcuad(ji,jj)
               zcvad(ji,jj) = 0.0_wp
               zcuad(ji,jj) = 0.0_wp
            END DO
         END DO

         ! Multiply by the eddy viscosity coef. (at u- and v-points)
         zluad(:,:,jk) = zluad(:,:,jk) * ahm3(:,:,jk)
         zlvad(:,:,jk) = zlvad(:,:,jk) * ahm4(:,:,jk)

      END DO

      ! Boundary conditions on the laplacian  (zlu,zlv)
      CALL lbc_lnk_adj( zlvad, 'V', -1.0_wp )
      CALL lbc_lnk_adj( zluad, 'U', -1.0_wp )

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         ! Laplacian
         ! ---------

         IF( ln_sco .OR. ln_zps ) THEN   ! s-coordinate or z-coordinate with partial steps
            DO jj = jpjm1, 2, -1
               DO ji = jpim1, 2, -1   ! vector opt.
                  zufad   (ji  ,jj  ,jk) = zufad   (ji  ,jj  ,jk) + zlvad(ji  ,jj  ,jk) / ( e1v(ji,jj) * e3v(ji,jj,jk) )
                  zufad   (ji-1,jj  ,jk) = zufad   (ji-1,jj  ,jk) - zlvad(ji  ,jj  ,jk) / ( e1v(ji,jj) * e3v(ji,jj,jk) )
                  hdivb_ad(ji  ,jj  ,jk) = hdivb_ad(ji  ,jj  ,jk) - zlvad(ji,jj,jk) / e2v(ji,jj)
                  hdivb_ad(ji  ,jj+1,jk) = hdivb_ad(ji  ,jj+1,jk) + zlvad(ji,jj,jk) / e2v(ji,jj)
                  zlvad(ji,jj,jk) = 0.0_wp

                  zufad   (ji  ,jj  ,jk) = zufad   (ji  ,jj  ,jk) - zluad(ji  ,jj  ,jk) / ( e2u(ji,jj) * e3u(ji,jj,jk) )
                  zufad   (ji  ,jj-1,jk) = zufad   (ji  ,jj-1,jk) + zluad(ji  ,jj  ,jk) / ( e2u(ji,jj) * e3u(ji,jj,jk) )
                  hdivb_ad(ji  ,jj  ,jk) = hdivb_ad(ji  ,jj  ,jk) - zluad(ji,jj,jk) / e1u(ji,jj)
                  hdivb_ad(ji+1,jj  ,jk) = hdivb_ad(ji+1,jj  ,jk) + zluad(ji,jj,jk) / e1u(ji,jj)
                  zluad(ji,jj,jk) = 0.0_wp
               END DO
            END DO
            rotb_ad(:,:,jk) = rotb_ad(:,:,jk) + zufad(:,:,jk) * e3f(:,:,jk)
         ELSE                            ! z-coordinate - full step
            DO jj = jpjm1, 2, -1
               DO ji = jpim1, 2, -1   ! vector opt.
                  rotb_ad (ji  ,jj  ,jk) = rotb_ad (ji  ,jj  ,jk) + zlvad(ji,jj,jk) / e1v(ji,jj)
                  rotb_ad (ji-1,jj  ,jk) = rotb_ad (ji-1,jj  ,jk) - zlvad(ji,jj,jk) / e1v(ji,jj)
                  hdivb_ad(ji  ,jj  ,jk) = hdivb_ad(ji  ,jj  ,jk) - zlvad(ji,jj,jk) / e2v(ji,jj)
                  hdivb_ad(ji  ,jj+1,jk) = hdivb_ad(ji  ,jj+1,jk) + zlvad(ji,jj,jk) / e2v(ji,jj)
                  zlvad(ji,jj,jk) = 0.0_wp
                  
                  rotb_ad (ji  ,jj  ,jk) = rotb_ad (ji  ,jj  ,jk) - zluad(ji,jj,jk) / e2u(ji,jj)
                  rotb_ad (ji  ,jj-1,jk) = rotb_ad (ji  ,jj-1,jk) + zluad(ji,jj,jk) / e2u(ji,jj)
                  hdivb_ad(ji  ,jj  ,jk) = hdivb_ad(ji  ,jj  ,jk) - zluad(ji,jj,jk) / e1u(ji,jj)
                  hdivb_ad(ji+1,jj  ,jk) = hdivb_ad(ji+1,jj  ,jk) + zluad(ji,jj,jk) / e1u(ji,jj)
                  zlvad(ji,jj,jk) = 0.0_wp
               END DO
            END DO
         ENDIF
      ENDDO
      CALL wrk_dealloc( jpi, jpj,      zcuad, zcvad           )
      CALL wrk_dealloc( jpi, jpj, jpk, zufad, zutad, zluad, zlvad )
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_bilap_adj')
      !
   END SUBROUTINE dyn_ldf_bilap_adj

   !!======================================================================
END MODULE dynldf_bilap_tam
