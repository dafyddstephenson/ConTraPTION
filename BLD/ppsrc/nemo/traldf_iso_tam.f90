MODULE traldf_iso_tam
   !!======================================================================
   !!                   ***  MODULE  traldf_iso_tam  ***
   !! Ocean active tracers:  horizontal component of the lateral tracer mixing trend
   !!                        Tangent and adjoint module
   !!======================================================================
   !! History of the direct module:
   !!                  !  94-08  (G. Madec, M. Imbard)
   !!                  !  97-05  (G. Madec)  split into traldf and trazdf
   !!             8.5  !  02-08  (G. Madec)  Free form, F90
   !!             9.0  !  05-11  (G. Madec)  merge traldf and trazdf :-)
   !! History of the T&A module:
   !!             9.0  !  2008-12 (A. Vidard) original version
   !!              -   !  2009-01 (A. Weaver) misc. bug fixes
   !!      NEMO   3.2  !  2010-04 (F. Vigilant) 3.2 version
   !!      NEMO   3.4  !  2012-07 (P.-A. Bouttier) 3.4 version
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_ldfslp'               slope of the lateral diffusive direction
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   tra_ldf_iso  : update the tracer trend with the horizontal
   !!                  component of a iso-neutral laplacian operator
   !!                  and with the vertical part of
   !!                  the isopycnal or geopotential s-coord. operator
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE ldftra_oce
!   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager
   USE ldfslp
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE paresp
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE trc_oce

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf_iso_tan     ! routine called by tralfd_tam.F90
   PUBLIC   tra_ldf_iso_adj     ! routine called by traldf_tam.F90
   PUBLIC   tra_ldf_iso_adj_tst ! routine called by traldf_tam.F90

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
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c2d' :                 aht: 2D coefficient
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
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_ldf_iso_tan( kt, kit000, cdtype, pgu_tl, pgv_tl,              &
      &                                ptb_tl, pta_tl, kjpt, pahtb0 )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_iso_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the before horizontal tracer (t & s) diffusive
      !!      trend for a laplacian tensor (ezxcept the dz[ dz[.] ] term) and
      !!      add it to the general trend of tracer equation.
      !!
      !! ** Method of the direct routine:
      !!      The horizontal component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      1st part :  masked horizontal derivative of T & S ( di[ t ] )
      !!      ========    with partial cell update if ln_zps=T.
      !!
      !!      2nd part :  horizontal fluxes of the lateral mixing operator
      !!      ========
      !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
      !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
      !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
      !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !!      3rd part: vertical trends of the lateral mixing operator
      !!      ========  (excluding the vertical flux proportional to dk[t] )
      !!      vertical fluxes associated with the rotated lateral mixing:
      !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
      !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ]
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !! ** Action :   Update (ta,sa) arrays with the before rotated diffusion
      !!            trend (except the dk[ dk[.] ] term)
      !!----------------------------------------------------------------------

      INTEGER, INTENT( in ) ::   kt    ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000           ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype           ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt             ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(in   ) ::   pgu_tl, pgv_tl   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb_tl           ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_tl           ! tracer trend
      REAL(wp)                             , INTENT(in   ) ::   pahtb0           ! background diffusion coef
      !!
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   iku, ikv     ! temporary integer
      REAL(wp) ::   zmsku, zabe1, zcof1, zcoef3, ztatl   ! temporary scalars
      REAL(wp) ::   zmskv, zabe2, zcof2, zcoef4, zsatl   !    "         "
      REAL(wp) ::   zcoef0, zbtr                       !    "         "
      REAL(wp), POINTER, DIMENSION(:,:)   ::   zdkttl , zdk1ttl, zftutl, zftvtl           ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zdittl, zdjttl, ztfwtl     ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_ldf_iso_tan')
      !
      CALL wrk_alloc( jpi, jpj,      zdkttl, zdk1ttl, zftutl, zftvtl )
      CALL wrk_alloc( jpi, jpj, jpk, zdittl, zdjttl, ztfwtl  )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_iso_tan : rotated laplacian diffusion operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      DO jn = 1, kjpt
         !!----------------------------------------------------------------------
         !!   I - masked horizontal derivative of T & S
         !!----------------------------------------------------------------------
   !!bug ajout.... why?   ( 1,jpj,:) and (jpi,1,:) should be sufficient....
         zdittl (1,:,:) = 0.e0     ;     zdittl (jpi,:,:) = 0.e0
         zdjttl (1,:,:) = 0.e0     ;     zdjttl (jpi,:,:) = 0.e0
   !!end

         ! Horizontal temperature and salinity gradient
         DO jk = 1, jpkm1
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  zdittl(ji,jj,jk) = ( ptb_tl(ji+1,jj  ,jk, jn) - ptb_tl(ji,jj,jk, jn) ) * umask(ji,jj,jk)
                  zdjttl(ji,jj,jk) = ( ptb_tl(ji  ,jj+1,jk, jn) - ptb_tl(ji,jj,jk, jn) ) * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         IF( ln_zps ) THEN      ! partial steps correction at the last level
            DO jj = 1, jpjm1
               DO ji = 1, jpim1   ! vector opt.
                  ! last level
                  zdittl(ji,jj,mbku(ji,jj)) = pgu_tl(ji,jj,jn)
                  zdjttl(ji,jj,mbkv(ji,jj)) = pgv_tl(ji,jj,jn)
               END DO
            END DO
         ENDIF

         !!----------------------------------------------------------------------
         !!   II - horizontal trend of T & S (full)
         !!----------------------------------------------------------------------

         !                                                ! ===============
         DO jk = 1, jpkm1                                 ! Horizontal slab
            !                                             ! ===============
            ! 1. Vertical tracer gradient at level jk and jk+1
            ! ------------------------------------------------
            ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)

            zdk1ttl(:,:) = ( ptb_tl(:,:,jk,jn) - ptb_tl(:,:,jk+1,jn) ) * tmask(:,:,jk+1)

            IF( jk == 1 ) THEN
               zdkttl(:,:) = zdk1ttl(:,:)
            ELSE
               zdkttl(:,:) = ( ptb_tl(:,:,jk-1,jn) - ptb_tl(:,:,jk,jn) ) * tmask(:,:,jk)
            ENDIF

            ! 2. Horizontal fluxes
            ! --------------------

            DO jj = 1 , jpjm1
               DO ji = 1, jpim1   ! vector opt.

                  zabe1 = ( rldf * ahtu(ji,jj) + pahtb0 ) * e2u(ji,jj) * e3u(ji,jj,jk) / e1u(ji,jj)
                  zabe2 = ( rldf * ahtv(ji,jj) + pahtb0 ) * e1v(ji,jj) * e3v(ji,jj,jk) / e2v(ji,jj)

                  zmsku = 1.0_wp / MAX(  tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
                     &                 + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1.0_wp )

                  zmskv = 1.0_wp / MAX(  tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
                     &                 + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1.0_wp )

                  ! *** NOTE ***  uslp() and vslp() are not linearized.

                  zcof1 = -rldf * ahtu(ji,jj) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku
                  zcof2 = -rldf * ahtv(ji,jj) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv

                  zftutl(ji,jj) = (  zabe1 * zdittl(ji,jj,jk)   &
                     &              + zcof1 * (  zdkttl (ji+1,jj) + zdk1ttl(ji,jj)      &
                     &                         + zdk1ttl(ji+1,jj) + zdkttl (ji,jj)  )  ) * umask(ji,jj,jk)
                  zftvtl(ji,jj) = (  zabe2 * zdjttl(ji,jj,jk)   &
                     &              + zcof2 * (  zdkttl (ji,jj+1) + zdk1ttl(ji,jj)      &
                     &                         + zdk1ttl(ji,jj+1) + zdkttl (ji,jj)  )  ) * vmask(ji,jj,jk)
               END DO
            END DO


            ! II.4 Second derivative (divergence) and add to the general trend
            ! ----------------------------------------------------------------
            DO jj = 2 , jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zbtr= 1.0_wp / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                  ztatl = zbtr * (   zftutl(ji,jj) - zftutl(ji-1,jj  ) &
                     &             + zftvtl(ji,jj) - zftvtl(ji  ,jj-1)  )
                  pta_tl(ji,jj,jk,jn) = pta_tl(ji,jj,jk,jn) + ztatl
               END DO
            END DO
            !                                          ! ===============
         END DO                                        !   End of slab
         !                                             ! ===============

         !!----------------------------------------------------------------------
         !!   III - vertical trend of T & S (extra diagonal terms only)
         !!----------------------------------------------------------------------

         ! Local constant initialization
         ! -----------------------------
         ztfwtl(1,:,:) = 0.0_wp     ;     ztfwtl(jpi,:,:) = 0.0_wp


         ! Vertical fluxes
         ! ---------------

         ! Surface and bottom vertical fluxes set to zero
         ztfwtl(:,:, 1 ) = 0.0_wp   ;     ztfwtl(:,:,jpk) = 0.0_wp

         ! interior (2=<jk=<jpk-1)
         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zcoef0 = - rldf * ahtw(ji,jj) * tmask(ji,jj,jk)

                  zmsku = 1.0_wp / MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
                     &                  + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.0_wp  )

                  zmskv = 1.0_wp / MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
                     &                  + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.0_wp  )

                  ! *** NOTE ***  wslpi() and wslpj() are not linearized.

                  zcoef3 = zcoef0 * e2t(ji,jj) * zmsku * wslpi(ji,jj,jk)
                  zcoef4 = zcoef0 * e1t(ji,jj) * zmskv * wslpj(ji,jj,jk)

                  ztfwtl(ji,jj,jk) = zcoef3 * ( zdittl(ji  ,jj  ,jk-1) + zdittl(ji-1,jj  ,jk)      &
                     &                        + zdittl(ji-1,jj  ,jk-1) + zdittl(ji  ,jj  ,jk)  )   &
                     &             + zcoef4 * ( zdjttl(ji  ,jj  ,jk-1) + zdjttl(ji  ,jj-1,jk)      &
                     &                        + zdjttl(ji  ,jj-1,jk-1) + zdjttl(ji  ,jj  ,jk)  )
               END DO
            END DO
         END DO


         ! I.5 Divergence of vertical fluxes added to the general tracer trend
         ! -------------------------------------------------------------------

         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zbtr =  1.0_wp / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )

                  ztatl = ( ztfwtl(ji,jj,jk) - ztfwtl(ji,jj,jk+1) ) * zbtr
                  pta_tl(ji,jj,jk,jn) = pta_tl(ji,jj,jk,jn) + ztatl
               END DO
            END DO
         END DO
      END DO
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_ldf_iso_tan')
      !
      CALL wrk_dealloc( jpi, jpj,      zdkttl, zdk1ttl, zftutl, zftvtl )
      CALL wrk_dealloc( jpi, jpj, jpk, zdittl, zdjttl, ztfwtl  )
      !
   END SUBROUTINE tra_ldf_iso_tan

   SUBROUTINE tra_ldf_iso_adj( kt, kit000, cdtype, pgu_ad, pgv_ad,              &
      &                                ptb_ad, pta_ad, kjpt, pahtb0  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_iso_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the before horizontal tracer (t & s) diffusive
      !!      trend for a laplacian tensor (ezxcept the dz[ dz[.] ] term) and
      !!      add it to the general trend of tracer equation.
      !!
      !! ** Method of the direct routine:
      !!      The horizontal component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neural or geopo-
      !!      tential surfaces to which an eddy induced advection can be added
      !!      It is computed using before fields (forward in time) and isopyc-
      !!      nal or geopotential slopes computed in routine ldfslp.
      !!
      !!      1st part :  masked horizontal derivative of T & S ( di[ t ] )
      !!      ========    with partial cell update if ln_zps=T.
      !!
      !!      2nd part :  horizontal fluxes of the lateral mixing operator
      !!      ========
      !!         zftu = (aht+ahtb0) e2u*e3u/e1u di[ tb ]
      !!               - aht       e2u*uslp    dk[ mi(mk(tb)) ]
      !!         zftv = (aht+ahtb0) e1v*e3v/e2v dj[ tb ]
      !!               - aht       e2u*vslp    dk[ mj(mk(tb)) ]
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) {  di-1[ zftu ] +  dj-1[ zftv ]  }
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !!      3rd part: vertical trends of the lateral mixing operator
      !!      ========  (excluding the vertical flux proportional to dk[t] )
      !!      vertical fluxes associated with the rotated lateral mixing:
      !!         zftw =-aht {  e2t*wslpi di[ mi(mk(tb)) ]
      !!                     + e1t*wslpj dj[ mj(mk(tb)) ]  }
      !!      take the horizontal divergence of the fluxes:
      !!         difft = 1/(e1t*e2t*e3t) dk[ zftw ]
      !!      Add this trend to the general trend (ta,sa):
      !!         ta = ta + difft
      !!
      !! ** Action :   Update (ta,sa) arrays with the before rotated diffusion
      !!            trend (except the dk[ dk[.] ] term)
      !!----------------------------------------------------------------------

      INTEGER                              , INTENT(in   ) ::   kt         ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype     ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt       ! number of tracers
      REAL(wp), DIMENSION(jpi,jpj    ,kjpt), INTENT(inout   ) ::   pgu_ad, pgv_ad   ! tracer gradient at pstep levels
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout   ) ::   ptb_ad        ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout   ) ::   pta_ad        ! tracer trend
      REAL(wp)                             , INTENT(in   ) ::   pahtb0     ! background diffusion coef
      !!
      INTEGER  ::   ji, jj, jk, jn   ! dummy loop indices
      INTEGER  ::   iku, ikv     ! temporary integer
      REAL(wp) ::   zmsku, zabe1, zcof1, zcoef3, ztaad   ! temporary scalars
      REAL(wp) ::   zmskv, zabe2, zcof2, zcoef4, zsaad   !    "         "
      REAL(wp) ::   ztf3, ztf4, zsf3, zsf4               !
      REAL(wp) ::   zcoef0, zbtr                       !    "         "
      REAL(wp), POINTER, DIMENSION(:,:)     ::   zdktad , zdk1tad, zftuad, zftvad   ! 2D workspace
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zditad, zdjtad, zdisad, zdjsad, ztfwad     ! 3D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_ldf_iso_adj')
      !
      CALL wrk_alloc( jpi, jpj,      zdktad, zdk1tad, zftuad, zftvad )
      CALL wrk_alloc( jpi, jpj, jpk, zditad, zdjtad, ztfwad, zdjsad, zdisad  )
      !
      zditad(:,:,:) = 0.0_wp ;  zdjtad(:,:,:) = 0.0_wp ;  ztfwad(:,:,:) = 0.0_wp
      zdktad(:,:) = 0.0_wp ;  zdk1tad(:,:) = 0.0_wp
      zftuad(:,:) = 0.0_wp ;  zftvad (:,:) = 0.0_wp

      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_ldf_iso_adj : rotated laplacian diffusion operator on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      DO jn = 1, kjpt
         !!----------------------------------------------------------------------
         !!   III - vertical trend of T & S (extra diagonal terms only)
         !!----------------------------------------------------------------------
         ! I.5 Divergence of vertical fluxes added to the general tracer trend
         ! -------------------------------------------------------------------

         DO jk = jpkm1, 1, -1
            DO jj = jpjm1, 2, -1
               DO ji = jpim1, 2, -1   ! vector opt.
                  zbtr =  1.0_wp / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                  ztaad = pta_ad(ji,jj,jk,jn) * zbtr
                  ztfwad(ji,jj,jk  ) = ztfwad(ji,jj,jk  ) + ztaad
                  ztfwad(ji,jj,jk+1) = ztfwad(ji,jj,jk+1) - ztaad
               END DO
            END DO
         END DO
         ! interior (2=<jk=<jpk-1)
         DO jk = jpkm1, 2, -1
            DO jj = jpjm1, 2, -1
               DO ji = jpim1, 2, -1   ! vector opt.
                  zcoef0 = - rldf * ahtw(ji,jj) * tmask(ji,jj,jk)

                  zmsku = 1.0_wp / MAX(   umask(ji  ,jj,jk-1) + umask(ji-1,jj,jk)      &
                     &                  + umask(ji-1,jj,jk-1) + umask(ji  ,jj,jk), 1.0_wp  )

                  zmskv = 1.0_wp / MAX(   vmask(ji,jj  ,jk-1) + vmask(ji,jj-1,jk)      &
                     &                  + vmask(ji,jj-1,jk-1) + vmask(ji,jj  ,jk), 1.0_wp  )

                  ! *** NOTE ***  wslpi() and wslpj() are not linearized.

                  zcoef3 = zcoef0 * e2t(ji,jj) * zmsku * wslpi(ji,jj,jk)
                  zcoef4 = zcoef0 * e1t(ji,jj) * zmskv * wslpj(ji,jj,jk)

                  ztf3 = ztfwad(ji,jj,jk) * zcoef3
                  ztf4 = ztfwad(ji,jj,jk) * zcoef4

                  zditad(ji  ,jj  ,jk-1) = zditad(ji  ,jj  ,jk-1) + ztf3
                  zditad(ji-1,jj  ,jk  ) = zditad(ji-1,jj  ,jk  ) + ztf3
                  zditad(ji-1,jj  ,jk-1) = zditad(ji-1,jj  ,jk-1) + ztf3
                  zditad(ji  ,jj  ,jk  ) = zditad(ji  ,jj  ,jk  ) + ztf3

                  zdjtad(ji  ,jj  ,jk-1) = zdjtad(ji  ,jj  ,jk-1) + ztf4
                  zdjtad(ji  ,jj-1,jk  ) = zdjtad(ji  ,jj-1,jk  ) + ztf4
                  zdjtad(ji  ,jj-1,jk-1) = zdjtad(ji  ,jj-1,jk-1) + ztf4
                  zdjtad(ji  ,jj  ,jk  ) = zdjtad(ji  ,jj  ,jk  ) + ztf4

                  ztfwad(ji,jj,jk) = 0.0_wp
               END DO
            END DO
         END DO

         ! Local constant initialization
         ! -----------------------------
         ztfwad(1,:,:) = 0.0_wp     ;     ztfwad(jpi,:,:) = 0.0_wp

         ! Vertical fluxes
         ! ---------------

         ! Surface and bottom vertical fluxes set to zero
         ztfwad(:,:, 1 ) = 0.0_wp   ;     ztfwad(:,:,jpk) = 0.0_wp

         !!----------------------------------------------------------------------
         !!   II - horizontal trend of T & S (full)
         !!----------------------------------------------------------------------

         !                                                ! ===============
         DO jk = jpkm1, 1, -1                             ! Horizontal slab
            !                                             ! ===============
            ! II.4 Second derivative (divergence) and add to the general trend
            ! ----------------------------------------------------------------
            DO jj = jpjm1, 2, -1
               DO ji = jpim1, 2, -1   ! vector opt.

                  zbtr= 1.0_wp / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                  ztaad = pta_ad(ji,jj,jk,jn) * zbtr

                  zftuad(ji  ,jj  ) = zftuad(ji  ,jj  ) + ztaad
                  zftuad(ji-1,jj  ) = zftuad(ji-1,jj  ) - ztaad
                  zftvad(ji  ,jj  ) = zftvad(ji  ,jj  ) + ztaad
                  zftvad(ji  ,jj-1) = zftvad(ji  ,jj-1) - ztaad
               END DO
            END DO

            ! 2. Horizontal fluxes
            ! --------------------
            DO jj = jpjm1, 1, -1
               DO ji = jpim1, 1, -1   ! vector opt.
                  zabe1 = umask(ji,jj,jk) * ( rldf * ahtu(ji,jj) + ahtb0 ) &
                     &                    * e2u(ji,jj) * e3u(ji,jj,jk) / e1u(ji,jj)
                  zabe2 = vmask(ji,jj,jk) * ( rldf * ahtv(ji,jj) + ahtb0 ) &
                     &                    * e1v(ji,jj) * e3v(ji,jj,jk) / e2v(ji,jj)

                  zmsku = 1.0_wp / MAX(  tmask(ji+1,jj,jk  ) + tmask(ji,jj,jk+1)   &
                     &                 + tmask(ji+1,jj,jk+1) + tmask(ji,jj,jk  ), 1.0_wp )

                  zmskv = 1.0_wp / MAX(  tmask(ji,jj+1,jk  ) + tmask(ji,jj,jk+1)   &
                     &                 + tmask(ji,jj+1,jk+1) + tmask(ji,jj,jk  ), 1.0_wp )

                  ! *** NOTE ***  uslp() and vslp() are not linearized.

                  zcof1 = -rldf * ahtu(ji,jj) * e2u(ji,jj) * uslp(ji,jj,jk) * zmsku * umask(ji,jj,jk)
                  zcof2 = -rldf * ahtv(ji,jj) * e1v(ji,jj) * vslp(ji,jj,jk) * zmskv * vmask(ji,jj,jk)

                  zditad(ji,jj,jk) = zditad(ji,jj,jk) + zftuad(ji,jj) * zabe1

                  zdktad (ji+1,jj) = zdktad (ji+1,jj) + zftuad(ji,jj) * zcof1
                  zdktad (ji  ,jj) = zdktad (ji  ,jj) + zftuad(ji,jj) * zcof1
                  zdk1tad(ji  ,jj) = zdk1tad(ji  ,jj) + zftuad(ji,jj) * zcof1
                  zdk1tad(ji+1,jj) = zdk1tad(ji+1,jj) + zftuad(ji,jj) * zcof1
                  zftuad (ji  ,jj) = 0.0_wp
                  !
                  zdjtad(ji,jj,jk) = zdjtad(ji,jj,jk) + zftvad(ji,jj) * zabe2

                  zdktad (ji,jj+1) = zdktad (ji,jj+1) + zftvad(ji,jj) * zcof2
                  zdktad (ji,jj  ) = zdktad (ji,jj  ) + zftvad(ji,jj) * zcof2
                  zdk1tad(ji,jj  ) = zdk1tad(ji,jj  ) + zftvad(ji,jj) * zcof2
                  zdk1tad(ji,jj+1) = zdk1tad(ji,jj+1) + zftvad(ji,jj) * zcof2
                  zftvad (ji,jj  ) = 0.0_wp
                  !
               END DO
            END DO

            ! 1. Vertical tracer gradient at level jk and jk+1
            ! ------------------------------------------------
            ! surface boundary condition: zdkt(jk=1)=zdkt(jk=2)

            IF( jk == 1 ) THEN
               zdk1tad(:,:) = zdk1tad(:,:) + zdktad(:,:)
               zdktad(:,:)  = 0.0_wp
            ELSE
               ptb_ad(:,:,jk-1,jn) = ptb_ad(:,:,jk-1,jn) + zdktad(:,:) * tmask(:,:,jk)
               ptb_ad(:,:,jk,jn  ) = ptb_ad(:,:,jk, jn ) - zdktad(:,:) * tmask(:,:,jk)
               zdktad(:,:) = 0.0_wp
            ENDIF
            ptb_ad(:,:,jk,jn  ) = ptb_ad(:,:,jk,jn  ) + zdk1tad(:,:) * tmask(:,:,jk+1)
            ptb_ad(:,:,jk+1,jn) = ptb_ad(:,:,jk+1,jn) - zdk1tad(:,:) * tmask(:,:,jk+1)
            zdk1tad(:,:) = 0.0_wp
            !                                          ! ===============
         END DO                                        !   End of slab
         !                                             ! ===============
         !!----------------------------------------------------------------------
         !!   I - masked horizontal derivative of T & S
         !!----------------------------------------------------------------------
         IF( ln_zps ) THEN      ! partial steps correction at the last level
            DO jj = jpjm1, 1, -1
               DO ji = jpim1, 1, -1   ! vector opt.
                  ! last level
                  pgu_ad(ji,jj,jn) = pgu_ad(ji,jj,jn) + zditad(ji,jj,mbku(ji,jj))
                  pgv_ad(ji,jj,jn) = pgv_ad(ji,jj,jn) + zdjtad(ji,jj,mbkv(ji,jj))

                  zditad(ji,jj,mbku(ji,jj)) = 0.0_wp
                  zdjtad(ji,jj,mbkv(ji,jj)) = 0.0_wp

                  zdisad(ji,jj,mbku(ji,jj)) = 0.0_wp
                  zdjsad(ji,jj,mbkv(ji,jj)) = 0.0_wp
              END DO
            END DO
         ENDIF

         ! Horizontal temperature and salinity gradient
         DO jk = jpkm1, 1, -1
            DO jj = jpjm1, 1, -1
               DO ji = jpim1, 1, -1   ! vector opt.
                  zditad(ji,jj,jk) = zditad(ji,jj,jk) * umask(ji,jj,jk)
                  zdjtad(ji,jj,jk) = zdjtad(ji,jj,jk) * vmask(ji,jj,jk)
                  ptb_ad(ji+1,jj  ,jk,jn) = ptb_ad(ji+1,jj  ,jk,jn) + zditad(ji,jj,jk)
                  ptb_ad(ji  ,jj  ,jk,jn) = ptb_ad(ji  ,jj  ,jk,jn) - zditad(ji,jj,jk)
                  ptb_ad(ji  ,jj+1,jk,jn) = ptb_ad(ji  ,jj+1,jk,jn) + zdjtad(ji,jj,jk)
                  ptb_ad(ji  ,jj  ,jk,jn) = ptb_ad(ji  ,jj  ,jk,jn) - zdjtad(ji,jj,jk)
                  zditad(ji,jj,jk) = 0._wp
                  zdjtad(ji,jj,jk) = 0._wp
               END DO
            END DO
         END DO
         zditad (1,:,:) = 0.e0     ;     zditad (jpi,:,:) = 0.e0
         zdjtad (1,:,:) = 0.e0     ;     zdjtad (jpi,:,:) = 0.e0
      END DO
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_ldf_iso_adj')
      !
      CALL wrk_dealloc( jpi, jpj,      zdktad, zdk1tad, zftuad, zftvad )
      CALL wrk_dealloc( jpi, jpj, jpk, zditad, zdjtad, ztfwad, zdjsad, zdisad  )
      !
   END SUBROUTINE tra_ldf_iso_adj

   SUBROUTINE tra_ldf_iso_adj_tst ( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE example_adj_tst ***
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
         & zsp1_T,       &
         & zsp1_S,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       &
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4,       &
         & zsp2_5,       &
         & zsp2_6,       &
         & zsp2_7,       &
         & zsp2_8,       &
         & zsp2_T,       &
         & zsp2_S
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & ztb_tlin ,      & ! Tangent input
         & zsb_tlin ,      & ! Tangent input
         & zta_tlin ,      & ! Tangent input
         & zsa_tlin ,      & ! Tangent input
         & zta_tlout,      & ! Tangent output
         & zsa_tlout,      & ! Tangent output
         & zta_adin,       & ! Adjoint input
         & zsa_adin,       & ! Adjoint input
         & ztb_adout ,     & ! Adjoint output
         & zsb_adout ,     & ! Adjoint output
         & zta_adout ,     & ! Adjoint output
         & zsa_adout ,     & ! Adjoint output
         & z3r               ! 3D random field
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zgtu_tlin ,     & ! Tangent input
         & zgsu_tlin ,     & ! Tangent input
         & zgtv_tlin ,     & ! Tangent input
         & zgsv_tlin ,     & ! Tangent input
         & zgtu_adout ,    & ! Adjoint output
         & zgsu_adout ,    & ! Adjoint output
         & zgtv_adout ,    & ! Adjoint output
         & zgsv_adout ,    & ! Adjoint output
         & z2r               ! 2D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & ztb_tlin(jpi,jpj,jpk),      &
         & zsb_tlin(jpi,jpj,jpk),      &
         & zta_tlin(jpi,jpj,jpk),      &
         & zsa_tlin(jpi,jpj,jpk),      &
         & zgtu_tlin(jpi,jpj),         &
         & zgsu_tlin(jpi,jpj),         &
         & zgtv_tlin(jpi,jpj),         &
         & zgsv_tlin(jpi,jpj),         &
         & zta_tlout(jpi,jpj,jpk),     &
         & zsa_tlout(jpi,jpj,jpk),     &
         & zta_adin(jpi,jpj,jpk),      &
         & zsa_adin(jpi,jpj,jpk),      &
         & ztb_adout(jpi,jpj,jpk),     &
         & zsb_adout(jpi,jpj,jpk),     &
         & zta_adout(jpi,jpj,jpk),     &
         & zsa_adout(jpi,jpj,jpk),     &
         & zgtu_adout(jpi,jpj),        &
         & zgsu_adout(jpi,jpj),        &
         & zgtv_adout(jpi,jpj),        &
         & zgsv_adout(jpi,jpj),        &
         & z3r(jpi,jpj,jpk),           &
         & z2r(jpi,jpj)                &
         & )

      ! Initialize the reference state
      uslp (:,:,:) = 2.0_wp
      vslp (:,:,:) = 3.0_wp
      wslpi(:,:,:) = 4.0_wp
      wslpj(:,:,:) = 5.0_wp

      !=======================================================================
      ! 1) dx = ( tb_tl, ta_tl, sb_tl, sa_tl, gtu_tl, gtv_tl, gsu_tl, gsv_tl )
      !    dy = ( ta_tl, sa_tl )
      !=======================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      ztb_tlin(:,:,:)  = 0.0_wp
      zsb_tlin(:,:,:)  = 0.0_wp
      zta_tlin(:,:,:)  = 0.0_wp
      zsa_tlin(:,:,:)  = 0.0_wp
      zgtu_tlin(:,:)   = 0.0_wp
      zgsu_tlin(:,:)   = 0.0_wp
      zgtv_tlin(:,:)   = 0.0_wp
      zgsv_tlin(:,:)   = 0.0_wp
      zta_tlout(:,:,:) = 0.0_wp
      zsa_tlout(:,:,:) = 0.0_wp
      zta_adin(:,:,:)  = 0.0_wp
      zsa_adin(:,:,:)  = 0.0_wp
      ztb_adout(:,:,:) = 0.0_wp
      zsb_adout(:,:,:) = 0.0_wp
      zta_adout(:,:,:) = 0.0_wp
      zsa_adout(:,:,:) = 0.0_wp
      zgtu_adout(:,:)  = 0.0_wp
      zgsu_adout(:,:)  = 0.0_wp
      zgtv_adout(:,:)  = 0.0_wp
      zgsv_adout(:,:)  = 0.0_wp

      tsb_tl(:,:,:,:) = 0.0_wp
      tsa_tl(:,:,:,:) = 0.0_wp
      gtsu_tl(:,:,:)  = 0.0_wp
      gtsv_tl(:,:,:)  = 0.0_wp
      tsb_ad(:,:,:,:) = 0.0_wp
      tsa_ad(:,:,:,:) = 0.0_wp
      gtsu_ad(:,:,:)  = 0.0_wp
      gtsv_ad(:,:,:)  = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              ztb_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsb_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsa_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z2r, 'U', 0.0_wp, stds )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgtu_tlin(ji,jj) = z2r(ji,jj)
         END DO
      END DO
      CALL grid_random(  z2r, 'U', 0.0_wp, stds )
      DO jj = nldj, nlej
        DO ji = nldi, nlei
           zgsu_tlin(ji,jj) = z2r(ji,jj)
        END DO
      END DO
      CALL grid_random(  z2r, 'V', 0.0_wp, stds )
      DO jj = nldj, nlej
        DO ji = nldi, nlei
           zgtv_tlin(ji,jj) = z2r(ji,jj)
        END DO
      END DO
      CALL grid_random(  z2r, 'V', 0.0_wp, stds )
      DO jj = nldj, nlej
        DO ji = nldi, nlei
           zgsv_tlin(ji,jj) = z2r(ji,jj)
        END DO
      END DO

      tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
      tsb_tl(:,:,:,jp_sal) = zsb_tlin(:,:,:)
      tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
      tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)
      gtsu_tl(:,:,jp_tem)  = zgtu_tlin(:,:)
      gtsu_tl(:,:,jp_sal)  = zgsu_tlin(:,:)
      gtsv_tl(:,:,jp_tem)  = zgtv_tlin(:,:)
      gtsv_tl(:,:,jp_sal)  = zgsv_tlin(:,:)

      CALL tra_ldf_iso_tan( nit000, nit000, 'TRA', gtsu_tl, gtsv_tl, tsb_tl, tsa_tl, jpts, ahtb0 )

      zta_tlout(:,:,:) = tsa_tl(:,:,:,jp_tem)
      zsa_tlout(:,:,:) = tsa_tl(:,:,:,jp_sal)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsa_adin(ji,jj,jk) = zsa_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_s(jk)
              zta_adin(ji,jj,jk) = zta_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_t(jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1_T = DOT_PRODUCT( zta_tlout, zta_adin )
      zsp1_S = DOT_PRODUCT( zsa_tlout, zsa_adin )
      zsp1 = zsp1_T + zsp1_S

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      tsa_ad(:,:,:,jp_tem) = zta_adin(:,:,:)
      tsa_ad(:,:,:,jp_sal) = zsa_adin(:,:,:)

      CALL tra_ldf_iso_adj( nit000, nit000, 'TRA', gtsu_ad, gtsv_ad, tsb_ad, tsa_ad, jpts, ahtb0 )

      ztb_adout(:,:,:) = tsb_ad(:,:,:,jp_tem)
      zsb_adout(:,:,:) = tsb_ad(:,:,:,jp_sal)
      zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)
      zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)
      zgtu_adout(:,:)  = gtsu_ad(:,:,jp_tem)
      zgsu_adout(:,:)  = gtsu_ad(:,:,jp_sal)
      zgtv_adout(:,:)  = gtsv_ad(:,:,jp_tem)
      zgsv_adout(:,:)  = gtsv_ad(:,:,jp_sal)

      zsp2_1 = DOT_PRODUCT( ztb_tlin , ztb_adout  )
      zsp2_2 = DOT_PRODUCT( zta_tlin , zta_adout  )
      zsp2_3 = DOT_PRODUCT( zgtu_tlin, zgtu_adout )
      zsp2_4 = DOT_PRODUCT( zgtv_tlin, zgtv_adout )
      zsp2_5 = DOT_PRODUCT( zsb_tlin , zsb_adout  )
      zsp2_6 = DOT_PRODUCT( zsa_tlin , zsa_adout  )
      zsp2_7 = DOT_PRODUCT( zgsu_tlin, zgsu_adout )
      zsp2_8 = DOT_PRODUCT( zgsv_tlin, zgsv_adout )

      zsp2_T = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4
      zsp2_S = zsp2_5 + zsp2_6 + zsp2_7 + zsp2_8
      zsp2   = zsp2_T + zsp2_S

      cl_name = 'tra_ldf_iso_ad'
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(         &
         & ztb_tlin,      & ! Tangent input
         & zsb_tlin,      & ! Tangent input
         & zta_tlin,      & ! Tangent input
         & zsa_tlin,      & ! Tangent input
         & zgtu_tlin,     & ! Tangent input
         & zgsu_tlin,     & ! Tangent input
         & zgtv_tlin,     & ! Tangent input
         & zgsv_tlin,     & ! Tangent input
         & zta_tlout,     & ! Tangent output
         & zsa_tlout,     & ! Tangent output
         & zta_adin,      & ! Adjoint input
         & zsa_adin,      & ! Adjoint input
         & ztb_adout,     & ! Adjoint output
         & zsb_adout,     & ! Adjoint output
         & zta_adout,     & ! Adjoint output
         & zsa_adout,     & ! Adjoint output
         & zgtu_adout,    & ! Adjoint output
         & zgsu_adout,    & ! Adjoint output
         & zgtv_adout,    & ! Adjoint output
         & zgsv_adout,    & ! Adjoint output
         & z3r,           & ! 3D random field
         & z2r            &
         & )

   END SUBROUTINE tra_ldf_iso_adj_tst

   !!==============================================================================
END MODULE traldf_iso_tam
