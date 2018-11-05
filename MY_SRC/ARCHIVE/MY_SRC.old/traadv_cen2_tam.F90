MODULE traadv_cen2_tam
#if defined key_tam
   !!======================================================================
   !!                     ***  MODULE  traadv_cen2_tam  ***
   !! Ocean active tracers:  horizontal & vertical advective trend
   !!                         Tangent and Adjoint module
   !!======================================================================
   !! History :  8.2  ! 2001-08  (G. Madec, E. Durand)  trahad+trazad=traadv
   !!            1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            9.0  ! 2004-08  (C. Talandier) New trends organization
   !!             -   ! 2005-11  (V. Garnier) Surface pressure gradient organization
   !!            2.0  ! 2006-04  (R. Benshila, G. Madec) Step reorganization
   !!             -   ! 2006-07  (G. madec)  add ups_orca_set routine
   !!            3.2  ! 2009-07  (G. Madec) add avmb, avtb in restart for cen2 advection
   !! History of the T&A module
   !!            9.0  ! 2008-12  (A. Vidard) original version
   !!            3.2  ! 2010-04  (F. Vigilant)
   !!            3.4  ! 2012-07  (P.-A. Bouttier)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv_cen2 : update the tracer trend with the horizontal and
   !!                  vertical advection trends using a seconder order
   !!   ups_orca_set : allow mixed upstream/centered scheme in specific
   !!                  area (set for orca 2 and 4 only)
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce
   USE oce_tam
   USE dom_oce
   USE trc_oce
   USE gridrandom
   USE dotprodfld
   USE in_out_manager
   USE zdf_oce
   USE tstool_tam
   USE paresp
   USE lib_mpp
   USE wrk_nemo
   use timing
! 2016-11-16: making diffusion co-efficients available
   USE ldftra_oce
   USE zdf_oce
! 2016-11-18: making diffusion co-efficients available (vertical)

! 2016-11-22
   USE iom
   USE trj_tam, ONLY: nn_ittrjfrq_tan

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_cen2_tan    ! routine called by traadv_tam.F90
   PUBLIC   tra_adv_cen2_adj    ! routine called by traadv_tam.F90
   PUBLIC   tra_adv_cen2_adj_tst! routine called by tst.F90

   !! * Substitutions
#  include "domzgr_substitute.h90"
#  include "vectopt_loop_substitute.h90"
! 2016-11-11: making diffusion coefficients available
#  include "ldftra_substitute.h90"
#  include "zdfddm_substitute.h90"
! 2016-11-18: making diffusion coefficients available (vertical)


!#  include "ldftra_oce.F90"
! 2016-11-16: added inclusion of ldftra_oce.F90 as values subbed by ldftra_substitute.h90
   !          not recognised / even more errors / commented out
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.2 , LOCEAN-IPSL (2009)
   !! $Id: traadv_cen2.F90 1201 2008-09-24 13:24:21Z rblod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_adv_cen2_tan( kt, kit000, pun, pvn, pwn, ptn, &
      &                          pun_tl, pvn_tl, pwn_tl, ptn_tl, pta_tl, kjpt, pr_ups )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_cen2_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method  :   The advection is evaluated by a second order centered
      !!      scheme using now fields (leap-frog scheme). In specific areas
      !!      (vicinity of major river mouths, some straits, or where tn is
      !!      approaching the freezing point) it is mixed with an upstream
      !!      scheme for stability reasons.
      !!         Part 0 : compute the upstream / centered flag
      !!                  (3D array, zind, defined at T-point (0<zind<1))
      !!         Part I : horizontal advection
      !!       * centered flux:
      !!               zcenu = e2u*e3u  un  mi(tn)
      !!               zcenv = e1v*e3v  vn  mj(tn)
      !!       * upstream flux:
      !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
      !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
      !!       * mixed upstream / centered horizontal advection scheme
      !!               zcofi = max(zind(i+1), zind(i))
      !!               zcofj = max(zind(j+1), zind(j))
      !!               zwx = zcofi * zupsu + (1-zcofi) * zcenu
      !!               zwy = zcofj * zupsv + (1-zcofj) * zcenv
      !!       * horizontal advective trend (divergence of the fluxes)
      !!               zta = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
      !!       * Add this trend now to the general trend of tracer (ta,sa):
      !!              (ta,sa) = (ta,sa) + ( zta , zsa )
      !!       * trend diagnostic ('key_trdtra' defined): the trend is
      !!      saved for diagnostics. The trends saved is expressed as
      !!      Uh.gradh(T), i.e.
      !!                     save trend = zta + tn divn
      !!         In addition, the advective trend in the two horizontal direc-
      !!      tion is also re-computed as Uh gradh(T). Indeed hadt+tn divn is
      !!      equal to (in s-coordinates, and similarly in z-coord.):
      !!         zta+tn*divn=1/(e1t*e2t*e3t) { mi-1( e2u*e3u  un  di[tn] )
      !!                                      +mj-1( e1v*e3v  vn  mj[tn] )  }
      !!         NB:in z-coordinate - full step (ln_zco=T) e3u=e3v=e3t, so
      !!      they vanish from the expression of the flux and divergence.
      !!
      !!         Part II : vertical advection
      !!      For temperature (idem for salinity) the advective trend is com-
      !!      puted as follows :
      !!            zta = 1/e3t dk+1[ zwz ]
      !!      where the vertical advective flux, zwz, is given by :
      !!            zwz = zcofk * zupst + (1-zcofk) * zcent
      !!      with
      !!        zupsv = upstream flux = wn * (tb(k) or tb(k-1) ) [wn>0 or <0]
      !!        zcenu = centered flux = wn * mk(tn)
      !!         The surface boundary condition is :
      !!      variable volume (lk_vvl = T) : zero advective flux
      !!      lin. free-surf  (lk_vvl = F) : wn(:,:,1) * tn(:,:,1)
      !!         Add this trend now to the general trend of tracer (ta,sa):
      !!            (ta,sa) = (ta,sa) + ( zta , zsa )
      !!         Trend diagnostic ('key_trdtra' defined): the trend is
      !!      saved for diagnostics. The trends saved is expressed as :
      !!             save trend =  w.gradz(T) = zta - tn divn.
      !!
      !! ** Action :  - update (ta,sa) with the now advective tracer trends
      !!              - save trends in (ztrdt,ztrds) ('key_trdtra')
      !!----------------------------------------------------------------------
      !!
      INTEGER , INTENT(in)                         ::   kt       ! ocean time-step index
      INTEGER , INTENT(in   )                      ::   kit000   ! first time step index
      INTEGER , INTENT(in   )                      ::   kjpt   ! first time step index
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pun_tl   ! ocean velocity u-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pvn_tl   ! ocean velocity v-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pwn_tl   ! ocean velocity w-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pun      ! ocean velocity u-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pvn      ! ocean velocity v-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pwn      ! ocean velocity w-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk,kjpt) ::   ptn     ! ocean velocity v-component
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt) ::   ptn_tl, pta_tl     ! ocean velocity w-component
      ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
      REAL(wp), OPTIONAL, INTENT(in) :: pr_ups ! lin.-interpolation ratio between fully centred sheme
                                     ! (pr_ups=0) and fully trajectory-upstream scheme (pr_ups=1)
      !!
      INTEGER  ::   ji, jj, jk, jn                       ! dummy loop indices
      REAL(wp) ::   zbtr, zhw, zhwtl,                 &  ! temporary scalars
         &          ze3tr, zfui  , zfuitl  ,          &  !    "         "
         &          zfvj  , zfvjtl, ztra                 !    "         "

      ! 2016-11-22
      REAL(wp) :: zsgmu   , zsgmv , zsgmw                   !    "         "
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zthtu , zthtv , zthtw  !    "         "
      CHARACTER(LEN=32) :: zctheta
      INTEGER :: itheta
      LOGICAL :: ithetaIO
      !2016-11-16 added sigma (sgm) and theta (tht) corresponding to parameters in Fiadeiro&Veronis1977
    
     ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
      REAL(wp) ::   zr_ups, zpu, zpv, zpw
      REAL(wp) ::   zice                                 !    -         -
      REAL(wp), POINTER, DIMENSION(:,:)     ::   ztfreez      ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwztl ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwxtl, zwytl ! 3D workspace
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_cen2_tan')
      ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
      zr_ups = 0.0_wp
      IF (PRESENT(pr_ups)) zr_ups = pr_ups
      !
      CALL wrk_alloc( jpi, jpj, ztfreez )
      CALL wrk_alloc( jpi, jpj, jpk, zwztl, zwxtl, zwytl )
      ! 2016-11-22
      IF (zr_ups < 0.0_wp) THEN
         CALL wrk_alloc(jpi, jpj, jpk, zthtu, zthtv, zthtw)
         IF (MOD(kt-kit000+1, nn_ittrjfrq_tan) == 0) THEN
            ithetaIO = .TRUE.
         ELSE
            ithetaIO = .FALSE.
         END IF
      END IF
           ithetaIO = .FALSE. !2016-11-30 ignoring theta output trigger
      !
      zwztl = 0._wp ; zwxtl = 0._wp ; zwytl = 0._wp
      !
      IF( kt == kit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_cen2_tan : 2nd order centered advection scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   Vector optimization case'
         IF(lwp) WRITE(numout,*)
         !
      ENDIF
      
      ! 2016-11-22
      IF (zr_ups < 0.0_wp) THEN
         IF (ithetaIO) THEN
            WRITE(zctheta, FMT='(A,I6.6,".nc")' ) , 'theta_', kt
            zctheta = TRIM(zctheta)
            CALL iom_open(zctheta, itheta, ldwrt=.TRUE., kiolib = jprstlib)
         END IF
         
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1
                  zthtu(ji,jj,jk) = ( pun(ji,jj,jk) / (e2t(ji,jj)*fse3t(ji,jj,jk)) ) / ( 2.0_wp * fsahtu(ji,jj,jk))
                  zthtv(ji,jj,jk) = ( pvn(ji,jj,jk) / (e1t(ji,jj)*fse3t(ji,jj,jk)) ) / ( 2.0_wp * fsahtv(ji,jj,jk))
               END DO
            END DO
         END DO
         IF (ithetaIO) CALL iom_rstput(kt, kt, itheta, 'thtu', zthtu)
         IF (ithetaIO) CALL iom_rstput(kt, kt, itheta, 'thtv', zthtv)
      END IF

      !
      ! I. Horizontal advection
      !    ====================
      !
    DO jn = 1, kjpt 
         DO jk = 1, jpkm1
            !                        ! Second order centered tracer flux at u- and v-points
            DO jj = 1, jpjm1
               DO ji = 1, fs_jpim1   ! vector opt.
               ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
                  ! 2016-11-15
                  IF (zr_ups >= 0.0_wp) THEN
 !                    zpu = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pun(ji,jj,jk))))
 !                    zpv = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pvn(ji,jj,jk))))

                      zpu = 0.5_wp
                      zpv = 0.5_wp !!! 2016-12-05 temporarily forced horizontal advection to be always centred

                  ELSE

!!! 2016-11-22
!                    zthtu = ( pun(ji,jj,jk) / (e2t(ji,jj)*fse3t(ji,jj,jk)) ) / ( 2.0_wp * fsahtu(ji,jj,jk))
!                    zthtv = ( pvn(ji,jj,jk) / (e1t(ji,jj)*fse3t(ji,jj,jk)) ) / ( 2.0_wp * fsahtv(ji,jj,jk))
! 2016-11-28 2 previous lines were commented, uncommented them. 2016-11-30 (recommented)
                     ! \theta = [ {u h} / {2 K} ]
                     
                     !           sgmu = (cosh(thtu)/sinh(thtu)) - (1/thtu)
                     !           sgmv = (cosh(thtv)/sinh(thtv)) - (1/thtv)
                     
!                     zsgmu = zthtu(ji,jj,jk) / (1.0_wp + abs(zthtu(ji,jj,jk)))
!                     zsgmv = zthtv(ji,jj,jk) / (1.0_wp + abs(zthtv(ji,jj,jk)))
                     zsgmu = 0.0_wp 
                     zsgmv = 0.0_wp  !2016-11-29 forced horizontal advection to be centred - TEMPORARY

                     ! \sigma = [ coth(\theta) - (1 / \theta) 
                     !zr_ups = sgmu
                     !zr_ups=\sigma
                     ! zpu = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pun(ji,jj,jk))))
                       zpu = 0.5_wp * (1.0_wp + zsgmu)
                     ! zr_ups = sgmv

                     !zpv = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pvn(ji,jj,jk))))
                       zpv = 0.5_wp * (1.0_wp + zsgmv)

                     !2016-11-16: modified scheme so negative value of rn_tl_adv_ups in namelist means value is
                     !ignored and NEMO reverts to weighted-mean scheme (c.f. Fiadeiro and Veronis, 1977)
                     !added sgmu,sgmv (sigma) and thtu,thtv (theta) parameters (as in page 6)

                     !2016-11-18: used expansion of coth(x)-1(x) = x/(1+|x|) for sigma, changed zpu,zpv so full range
                     !            of sigma (-1,1) is realised. Previous zr_ups only in (0,1) 
                  ENDIF
                  ! volume fluxes * 1/2
                  zfuitl = zpu * pun_tl(ji,jj,jk)
                  zfvjtl = zpv * pvn_tl(ji,jj,jk)
                  zfui   = zpu * pun(   ji,jj,jk)
                  zfvj   = zpv * pvn(   ji,jj,jk)
                  !
                  ! centered scheme
                  
                  zwxtl(ji,jj,jk) = zfuitl * ptn(   ji,jj,jk, jn) + (pun_tl(ji,jj,jk) - zfuitl) * ptn(   ji+1,jj  ,jk, jn) &
                     &            + zfui   * ptn_tl(ji,jj,jk, jn) + (pun   (ji,jj,jk) - zfui  ) * ptn_tl(ji+1,jj  ,jk, jn)
                  zwytl(ji,jj,jk) = zfvjtl * ptn(   ji,jj,jk, jn) + (pvn_tl(ji,jj,jk) - zfvjtl) * ptn(   ji  ,jj+1,jk, jn) &
                     &            + zfvj   * ptn_tl(ji,jj,jk, jn) + (pvn   (ji,jj,jk) - zfvj  ) * ptn_tl(ji  ,jj+1,jk, jn)
             END DO
            END DO
         END DO

         ! "zonal" mean advective heat and salt transport
         ! ----------------------------------------------

         ! II. Vertical advection
         !     ==================
         !
         zwztl(:,:,jpk) = 0.0_wp                              ! Bottom value : flux set to zero

         IF( lk_vvl ) THEN
            zwztl(:,:, 1 ) = 0.e0                         ! volume variable
         ELSE
            zwztl(:,:, 1 ) = pwn(:,:,1) * ptn_tl(:,:,1,jn) &   ! linear free surface
               &           + pwn_tl(:,:,1) * ptn(:,:,1,jn)
         ENDIF
         !

        ! ##########################################################################

         ! 2016-11-22
         IF (zr_ups < 0.0_wp) THEN
            DO jk = 2, jpk
               DO jj = 2, jpjm1
                  DO ji = fs_2, fs_jpim1
                     IF (jn==1) THEN
                        IF (avt(ji,jj,jk)>0.0_wp) THEN !2016-11-29 added conditional statement setting theta=0 everywhere vertical diffusion is 0 (e.g. land) to avoid singularity
                        zthtw(ji,jj,jk) = ( pwn(ji,jj,jk) /( e1t(ji,jj)*e2t(ji,jj)) ) / ( 2.0_wp * avt(ji,jj,jk))          
                        ELSE
                        zthtw(ji,jj,jk) = 0.0_wp
                        END IF 
                     ELSE !2016-11-29 applied above change to salinity field
                        IF (fsavs(ji,jj,jk)>0.0_wp) THEN !2016-11-29 added conditional statement setting theta=0 everywhere vertical diffusion is 0 (e.g. land) to avoid singularity
                        zthtw(ji,jj,jk) = ( pwn(ji,jj,jk) /( e1t(ji,jj)*e2t(ji,jj)) ) / ( 2.0_wp * fsavs(ji,jj,jk))          
                        ELSE
                        zthtw(ji,jj,jk) = 0.0_wp
                        END IF
                     ENDIF
                  END DO
               END DO
            END DO
            IF (ithetaIO) THEN
               IF (jn==jp_tem) THEN
                  CALL iom_rstput(kt, kt, itheta, 'thtw_1', zthtw)
               ELSEIF (jn==jp_sal) THEN
                  CALL iom_rstput(kt, kt, itheta, 'thtw_2', zthtw)
               END IF
            END IF
         END IF

         DO jk = 2, jpk              ! Second order centered tracer flux at w-point
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
               ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
                  
                  IF (zr_ups >= 0.0_wp) THEN
                  
                  zpw = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pwn(ji,jj,jk))))

!!! 2016-11-22                  
                  ELSE !2016-11-23: this was commented out. Uncommented.
!!!                     IF (jn==1) THEN 
!!!                        thtw = ( pwn(ji,jj,jk) /( e1t(ji,jj)*e2t(ji,jj)) ) / ( 2.0_wp * avt(ji,jj,jk))          
!!!                     ELSE
!!!                        thtw = ( pwn(ji,jj,jk) /( e1t(ji,jj)*e2t(ji,jj)) ) / ( 2.0_wp * fsavs(ji,jj,jk))          
!!!                     ENDIF
                  zsgmw = zthtw(ji,jj,jk) / (1.0_wp + abs(zthtw(ji,jj,jk)))
                  zpw = 0.5_wp * ( 1.0_wp + zsgmw )
                  
                  !2016-11-18: repeated above adaptations for vertical advection

                  ENDIF

                  ! velocity * 1/2
                  zhwtl = zpw * pwn_tl(ji,jj,jk)
                  zhw   = zpw * pwn(   ji,jj,jk)
                  ! centered scheme
                  zwztl(ji,jj,jk) = zhwtl * ptn(   ji,jj,jk,jn) + (pwn_tl(ji,jj,jk) - zhwtl) * ptn(   ji,jj,jk-1,jn) &
                     &            + zhw   * ptn_tl(ji,jj,jk,jn) + (pwn   (ji,jj,jk) - zhw  ) * ptn_tl(ji,jj,jk-1,jn)
               END DO
            END DO
         END DO



         ! #################################################

         !
         DO jk = 1, jpkm1            ! divergence of Tracer flux added to the general trend
            DO jj = 2, jpjm1
               DO ji = fs_2, fs_jpim1   ! vector opt.
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) *  fse3t(ji,jj,jk) )
                  ! vertical advective trends
                  ztra = - zbtr * (  zwxtl(ji,jj,jk) - zwxtl(ji-1,jj  ,jk  )   &
                  &                + zwytl(ji,jj,jk) - zwytl(ji  ,jj-1,jk  )   &
                  &                + zwztl(ji,jj,jk) - zwztl(ji  ,jj  ,jk+1)  )
                  ! advective trends added to the general tracer trends
                  pta_tl(ji,jj,jk,jn) = pta_tl(ji,jj,jk,jn) + ztra
               END DO
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, ztfreez )
      CALL wrk_dealloc( jpi, jpj, jpk, zwztl, zwxtl, zwytl )

      ! 2016-11-22
      IF (zr_ups < 0.0_wp) THEN
         CALL wrk_dealloc(jpi, jpj, jpk, zthtu, zthtv, zthtw)
         IF (ithetaIO) CALL iom_close(itheta)
      END IF

      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_cen2_tan')
      !
      !
   END SUBROUTINE tra_adv_cen2_tan

   SUBROUTINE tra_adv_cen2_adj( kt, kit000, pun, pvn, pwn, ptn, &
      &                          pun_ad, pvn_ad, pwn_ad, ptn_ad, pta_ad, kjpt, pr_ups )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_cen2_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations.
      !!
      !! ** Method  :   The advection is evaluated by a second order centered
      !!      scheme using now fields (leap-frog scheme). In specific areas
      !!      (vicinity of major river mouths, some straits, or where tn is
      !!      approaching the freezing point) it is mixed with an upstream
      !!      scheme for stability reasons.
      !!         Part 0 : compute the upstream / centered flag
      !!                  (3D array, zind, defined at T-point (0<zind<1))
      !!         Part I : horizontal advection
      !!       * centered flux:
      !!               zcenu = e2u*e3u  un  mi(tn)
      !!               zcenv = e1v*e3v  vn  mj(tn)
      !!       * upstream flux:
      !!               zupsu = e2u*e3u  un  (tb(i) or tb(i-1) ) [un>0 or <0]
      !!               zupsv = e1v*e3v  vn  (tb(j) or tb(j-1) ) [vn>0 or <0]
      !!       * mixed upstream / centered horizontal advection scheme
      !!               zcofi = max(zind(i+1), zind(i))
      !!               zcofj = max(zind(j+1), zind(j))
      !!               zwx = zcofi * zupsu + (1-zcofi) * zcenu
      !!               zwy = zcofj * zupsv + (1-zcofj) * zcenv
      !!       * horizontal advective trend (divergence of the fluxes)
      !!               zta = 1/(e1t*e2t*e3t) { di-1[zwx] + dj-1[zwy] }
      !!       * Add this trend now to the general trend of tracer (ta,sa):
      !!              (ta,sa) = (ta,sa) + ( zta , zsa )
      !!       * trend diagnostic ('key_trdtra' defined): the trend is
      !!      saved for diagnostics. The trends saved is expressed as
      !!      Uh.gradh(T), i.e.
      !!                     save trend = zta + tn divn
      !!         In addition, the advective trend in the two horizontal direc-
      !!      tion is also re-computed as Uh gradh(T). Indeed hadt+tn divn is
      !!      equal to (in s-coordinates, and similarly in z-coord.):
      !!         zta+tn*divn=1/(e1t*e2t*e3t) { mi-1( e2u*e3u  un  di[tn] )
      !!                                      +mj-1( e1v*e3v  vn  mj[tn] )  }
      !!         NB:in z-coordinate - full step (ln_zco=T) e3u=e3v=e3t, so
      !!      they vanish from the expression of the flux and divergence.
      !!
      !!         Part II : vertical advection
      !!      For temperature (idem for salinity) the advective trend is com-
      !!      puted as follows :
      !!            zta = 1/e3t dk+1[ zwz ]
      !!      where the vertical advective flux, zwz, is given by :
      !!            zwz = zcofk * zupst + (1-zcofk) * zcent
      !!      with
      !!        zupsv = upstream flux = wn * (tb(k) or tb(k-1) ) [wn>0 or <0]
      !!        zcenu = centered flux = wn * mk(tn)
      !!         The surface boundary condition is :
      !!      rigid-lid (lk_dynspg_frd = T) : zero advective flux
      !!      free-surf (lk_dynspg_fsc = T) : wn(:,:,1) * tn(:,:,1)
      !!         Add this trend now to the general trend of tracer (ta,sa):
      !!            (ta,sa) = (ta,sa) + ( zta , zsa )
      !!         Trend diagnostic ('key_trdtra' defined): the trend is
      !!      saved for diagnostics. The trends saved is expressed as :
      !!             save trend =  w.gradz(T) = zta - tn divn.
      !!
      !! ** Action :  - update (ta,sa) with the now advective tracer trends
      !!              - save trends in (ztrdt,ztrds) ('key_trdtra')
      !!----------------------------------------------------------------------
      !!
      INTEGER , INTENT(in)                         ::   kt          ! ocean time-step index
      INTEGER , INTENT(in)                         ::   kit000
      INTEGER , INTENT(in)                         ::   kjpt
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pun_ad   ! ocean velocity u-component
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pvn_ad   ! ocean velocity v-component
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk) ::   pwn_ad   ! ocean velocity w-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pun      ! ocean velocity u-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pvn      ! ocean velocity v-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk) ::   pwn      ! ocean velocity w-component
      REAL(wp), INTENT(in), DIMENSION(jpi,jpj,jpk,kjpt) ::   ptn
      REAL(wp), INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt) ::   ptn_ad, pta_ad
      ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
      REAL(wp), OPTIONAL, INTENT(in) :: pr_ups ! lin.-interpolation ratio between fully centred sheme
                                     ! (pr_ups=0) and fully trajectory-upstream scheme (pr_ups=1)
      !!
      INTEGER  ::   ji, jj, jk, jn                       ! dummy loop indices
      REAL(wp) ::   zbtr, zhw, zhwad,                 &  ! temporary scalars
         &          ze3tr, zfui  , zfuiad  ,          &  !    "         "
         &          zfvj , zfvjad                        !    "         "
      ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
      REAL(wp) ::   zr_ups, zpu, zpv, zpw
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwzad ! 3D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwxad, zwyad ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_cen2_adj')
      ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
      zr_ups = 0.0_wp
      IF (PRESENT(pr_ups)) zr_ups = pr_ups
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwzad, zwyad, zwxad )
      !
      zhwad = 0.0_wp  ;  zfuiad = 0.0_wp  ;  zfvjad = 0.0_wp
      zwxad(:,:,:) = 0.0_wp ; zwyad(:,:,:) = 0.0_wp ; zwzad(:,:,:) = 0.0_wp

      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_cen2_adj : 2nd order centered advection scheme'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~   Vector optimization case'
         IF(lwp) WRITE(numout,*)
         !
      ENDIF
      ! II. Vertical advection
      !     ==================
      !
      DO jn = 1, kjpt
         DO jk = jpkm1, 1, -1            ! divergence of Tracer flux added to the general trend
            DO jj = jpjm1, 2, -1
               DO ji = fs_jpim1, fs_2, -1   ! vector opt.
                  zbtr = 1. / ( e1t(ji,jj) * e2t(ji,jj) *  fse3t(ji,jj,jk) )
                  zwxad(ji,jj,jk) = zwxad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk,jn)
                  zwyad(ji,jj,jk) = zwyad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk,jn)
                  zwzad(ji,jj,jk) = zwzad(ji,jj,jk) - zbtr * pta_ad(ji,jj,jk,jn)
                  zwxad(ji-1,jj,jk) = zwxad(ji-1,jj,jk) + zbtr * pta_ad(ji,jj,jk,jn)
                  zwyad(ji,jj-1,jk) = zwyad(ji,jj-1,jk) + zbtr * pta_ad(ji,jj,jk,jn)
                  zwzad(ji,jj,jk+1) = zwzad(ji,jj,jk+1) + zbtr * pta_ad(ji,jj,jk,jn)
               END DO
            END DO
         END DO
         !
         DO jk = jpk, 2, -1              ! Second order centered tracer flux at w-point
            DO jj = jpjm1, 2, -1
               DO ji = fs_jpim1, fs_2, -1   ! vector opt.
               ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
                  zpw   = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pwn(ji,jj,jk))))
                  zhw   = zpw * pwn(   ji,jj,jk)
                  ! centered scheme
                  zhwad = zwzad(ji,jj,jk) * ( zpw * ptn(   ji,jj,jk,jn) + (1.0_wp-zpw) * ptn(   ji,jj,jk-1,jn) )
                  ptn_ad(ji,jj,jk,jn) = ptn_ad(ji,jj,jk,jn) + zhw * zwzad(ji,jj,jk)
                  ptn_ad(ji,jj,jk-1,jn) = ptn_ad(ji,jj,jk-1,jn) + (pwn(ji,jj,jk) - zhw) * zwzad(ji,jj,jk)
                  zwzad(ji,jj,jk) = 0.0_wp
                  pwn_ad(ji,jj,jk) = pwn_ad(ji,jj,jk) + zhwad
                  zhwad = 0._wp
              END DO
            END DO
         END DO

         IF( lk_vvl ) THEN                                         ! Surface value : zero in variable volume
            zwzad(:,:, 1 ) = 0.0_wp
         ELSE                                                      !               : linear free surface case
            pwn_ad(:,:,1) = pwn_ad(:,:,1) + zwzad(:,:, 1 ) * ptn(:,:,1,jn)
            ptn_ad(:,:,1,jn)  = ptn_ad(:,:,1,jn)  + zwzad(:,:, 1 ) * pwn(:,:,1)
            zwzad(:,:, 1 ) = 0.0_wp
         ENDIF
         !
         zwzad(:,:,jpk) = 0.0_wp                                   ! Bottom value  : flux set to zero
         ! "zonal" mean advective heat and salt transport
         ! ----------------------------------------------
         ! I. Horizontal advective fluxes
         !    ====================
         !
         DO jk = 1, jpkm1
            !                        ! Second order centered tracer flux at u- and v-points
            DO jj = jpjm1, 1, -1
               DO ji = fs_jpim1, 1, -1   ! vector opt.
               ! 2016-05-20 - SAM: added optional trajectory-upstream advection scheme
                  zpu = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pun(ji,jj,jk))))
                  zpv = 0.5_wp * (1.0_wp-zr_ups+zr_ups*(1.0_wp+sign(1.0_wp,pvn(ji,jj,jk))))
                  ! volume fluxes * 1/2
                  zfui = zpu * pun(ji,jj,jk)
                  zfvj = zpv * pvn(ji,jj,jk)
                  ! centered scheme
                  zfvjad = zwyad(ji,jj,jk) * ( zpv * ptn(ji,jj,jk,jn) + (1.0_wp-zpv) * ptn(ji,jj+1,jk,jn) )
                  zfuiad = zwxad(ji,jj,jk) * ( zpu * ptn(ji,jj,jk,jn) + (1.0_wp-zpu) * ptn(ji+1,jj,jk,jn) )
                  ptn_ad(ji  ,jj  ,jk,jn) = ptn_ad(ji  ,jj  ,jk,jn) + zfvj * zwyad(ji,jj,jk)
                  ptn_ad(ji  ,jj+1,jk,jn) = ptn_ad(ji  ,jj+1,jk,jn) + (pvn(ji,jj,jk) - zfvj) * zwyad(ji,jj,jk)
                  ptn_ad(ji  ,jj  ,jk,jn) = ptn_ad(ji  ,jj  ,jk,jn) + zfui * zwxad(ji,jj,jk)
                  ptn_ad(ji+1,jj  ,jk,jn) = ptn_ad(ji+1,jj  ,jk,jn) + (pun(ji,jj,jk) - zfui) * zwxad(ji,jj,jk)
                  zwyad(ji  ,jj  ,jk) = 0.0_wp
                  zwxad(ji  ,jj  ,jk) = 0.0_wp
                  ! volume fluxes * 1/2
                  pun_ad(ji,jj,jk) = pun_ad(ji,jj,jk) + zfuiad
                  pvn_ad(ji,jj,jk) = pvn_ad(ji,jj,jk) + zfvjad
                  zfuiad = 0._wp
                  zfvjad = 0._wp
                END DO
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwzad, zwyad, zwxad )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_cen2_adj')
      !
   END SUBROUTINE tra_adv_cen2_adj
SUBROUTINE tra_adv_cen2_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE tra_adv_cen2_adj_tst ***
      !!
      !! ** Purpose : Test the adjoint routine.
      !!
      !! ** Method  : Verify the scalar product
      !!
      !!                 ( L dx )^T W dy  =  dx^T L^T W dy
      !!
      !!              where  L   = tangent routine
      !!                     L^T = adjoint routine
      !!                     W   = diagonal matrix of scale factors
      !!                     dx  = input perturbation (random field)
      !!                     dy  = L dx
      !!
      !!
      !! History :
      !!        ! 08-08 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      !! * Local declarations
      INTEGER ::  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zun_tlin ,     zvn_tlin ,     zwn_tlin ,     ztn_tlin ,     & ! Tangent input
         & zsn_tlin ,     zta_tlin ,     zsa_tlin ,     & ! Tangent input
         & zun_adout,     zvn_adout,     zwn_adout,     ztn_adout,     & ! Adjoint output
         & zsn_adout,     zta_adout,     zsa_adout,     & ! Adjoint output
         & zta_tlout,     zsa_tlout,     & ! Tangent output
         & zta_adin ,     zsa_adin ,     & ! Adjoint input
         & zr             ! 3D random field
      CHARACTER(LEN=14) ::&
         & cl_name
      ! Allocate memory

      ALLOCATE( &
         & zun_tlin( jpi,jpj,jpk),     zvn_tlin( jpi,jpj,jpk),     zwn_tlin( jpi,jpj,jpk),     &
         & ztn_tlin( jpi,jpj,jpk),     zsn_tlin( jpi,jpj,jpk),     zta_tlin( jpi,jpj,jpk),     &
         & zsa_tlin( jpi,jpj,jpk),     zta_tlout(jpi,jpj,jpk),     zsa_tlout(jpi,jpj,jpk),     &
         & zta_adin( jpi,jpj,jpk),     zsa_adin( jpi,jpj,jpk),     zun_adout(jpi,jpj,jpk),     &
         & zvn_adout(jpi,jpj,jpk),     zwn_adout(jpi,jpj,jpk),     ztn_adout(jpi,jpj,jpk),     &
         & zsn_adout(jpi,jpj,jpk),     zta_adout(jpi,jpj,jpk),     zsa_adout(jpi,jpj,jpk),     &
         & zr(       jpi,jpj,jpk)      &
         & )
      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      zun_tlin( :,:,:) = 0.0_wp
      zvn_tlin( :,:,:) = 0.0_wp
      zwn_tlin( :,:,:) = 0.0_wp
      ztn_tlin( :,:,:) = 0.0_wp
      zsn_tlin( :,:,:) = 0.0_wp
      zta_tlin( :,:,:) = 0.0_wp
      zsa_tlin( :,:,:) = 0.0_wp
      zta_tlout(:,:,:) = 0.0_wp
      zsa_tlout(:,:,:) = 0.0_wp
      zta_adin( :,:,:) = 0.0_wp
      zsa_adin( :,:,:) = 0.0_wp
      zun_adout(:,:,:) = 0.0_wp
      zvn_adout(:,:,:) = 0.0_wp
      zwn_adout(:,:,:) = 0.0_wp
      ztn_adout(:,:,:) = 0.0_wp
      zsn_adout(:,:,:) = 0.0_wp
      zta_adout(:,:,:) = 0.0_wp
      zsa_adout(:,:,:) = 0.0_wp
      zr(       :,:,:) = 0.0_wp

      tsn_ad(:,:,:,jp_tem) = 0.0_wp
      tsn_ad(:,:,:,jp_sal) = 0.0_wp


      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  zr, 'U', 0.0_wp, stdu )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zun_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'V', 0.0_wp, stdv )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zvn_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'W', 0.0_wp, stdw )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zwn_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              ztn_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsn_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsa_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO

      tsn_tl(:,:,:,jp_tem) = ztn_tlin(:,:,:)
      tsn_tl(:,:,:,jp_sal) = zsn_tlin(:,:,:)
      tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
      tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)

      CALL tra_adv_cen2_tan(nit000, nit000, un, vn, wn, tsn, zun_tlin, zvn_tlin, zwn_tlin, tsn_tl, tsa_tl, 2)

      zta_tlout(:,:,:) = tsa_tl(:,:,:,jp_tem)
      zsa_tlout(:,:,:) = tsa_tl(:,:,:,jp_sal)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_adin(ji,jj,jk) = zta_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_t(jk)
              zsa_adin(ji,jj,jk) = zsa_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * fse3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_s(jk)
            END DO
         END DO
      END DO
      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 = DOT_PRODUCT( zsa_tlout, zsa_adin ) &
         & + DOT_PRODUCT( zta_tlout, zta_adin )

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------
      tsa_ad(:,:,:,jp_tem) = zta_adin(:,:,:)
      tsa_ad(:,:,:,jp_sal) = zsa_adin(:,:,:)

      CALL tra_adv_cen2_adj(nit000, nit000, un, vn, wn, tsn, zun_adout, zvn_adout, zwn_adout, tsn_ad, tsa_ad, 2)

      ztn_adout(:,:,:) = tsn_ad(:,:,:,jp_tem)
      zsn_adout(:,:,:) = tsn_ad(:,:,:,jp_sal)
      zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)
      zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)

      zsp2 = DOT_PRODUCT( zun_tlin, zun_adout ) &
         & + DOT_PRODUCT( zvn_tlin, zvn_adout ) &
         & + DOT_PRODUCT( zwn_tlin, zwn_adout ) &
         & + DOT_PRODUCT( ztn_tlin, ztn_adout ) &
         & + DOT_PRODUCT( zsn_tlin, zsn_adout ) &
         & + DOT_PRODUCT( zta_tlin, zta_adout ) &
         & + DOT_PRODUCT( zsa_tlin, zsa_adout )

      ! 14 char:'12345678901234'
      cl_name = 'tra_adv_cen2  '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(         &
         & zun_tlin ,     zvn_tlin ,     zwn_tlin ,     ztn_tlin ,     zsn_tlin ,     &
         & zta_tlin ,     zsa_tlin ,     zta_tlout,     zsa_tlout,     zta_adin ,     &
         & zsa_adin ,     zun_adout,     zvn_adout,     zwn_adout,     ztn_adout,     &
         & zsn_adout,     zta_adout,     zsa_adout,     zr                            &
         & )

   END SUBROUTINE tra_adv_cen2_adj_tst
#endif
   !!======================================================================
END MODULE traadv_cen2_tam
