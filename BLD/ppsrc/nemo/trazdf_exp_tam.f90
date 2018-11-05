MODULE trazdf_exp_tam
   !!==============================================================================
   !!                    ***  MODULE  trazdf_exp_tam  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend using
   !!                        a split-explicit time-stepping
   !!                        Tangent and Adjoint module
   !!==============================================================================
   !! History of the direct module :
   !!   OPA           !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1992-06  (M. Imbard)  correction on tracer trend loops
   !!                 !  1996-01  (G. Madec)  statement function for e3
   !!                 !  1997-05  (G. Madec)  vertical component of isopycnal
   !!                 !  1997-07  (G. Madec)  geopotential diffusion in s-coord
   !!                 !  2000-08  (G. Madec)  double diffusive mixing
   !!   NEMO     1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2004-08  (C. Talandier) New trends organisation
   !!             -   !  2005-11  (G. Madec)  New organisation
   !!            3.0  !  2008-04  (G. Madec)  leap-frog time stepping done in trazdf
   !! History of the T&A module :
   !!                 !  2009-01  (A. Vidard) tam of the 2008-04 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_zdf_exp_tan  : compute the tracer the vertical diffusion trend using a
   !!                  split-explicit time stepping and provide the after tracer (tangent)
   !!   tra_zdf_exp_adj  : compute the tracer the vertical diffusion trend using a
   !!                  split-explicit time stepping and provide the after tracer (adjoint)
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE zdf_oce
   USE zdfddm
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE tstool_tam
   USE trc_oce
   USE trc_oce_tam
   USE lib_mpp
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_zdf_exp_tan       ! routine called by tra_zdf_tan.F90
   PUBLIC   tra_zdf_exp_adj       ! routine called by tra_zdf_adj.F90
   PUBLIC   tra_zdf_exp_adj_tst   ! routine called by tst.F90

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
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   Defautl option :                     avs = avt
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
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
   !! NEMO/OPA  3.0 , LOCEAN-IPSL (2008)
   !! $Id: trazdf_exp.F90 1146 2008-06-25 11:42:56Z rblod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE tra_zdf_exp_tan( kt, kit000, cdtype, p2dt, kn_zdfexp,   &
      &                                ptb_tl , pta_tl      , kjpt  )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_exp_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the after tracer fields due to the vertical
      !!      tracer mixing alone, and then due to the whole tracer trend.
      !!
      !! ** Method of the direct routine :
      !!               - The after tracer fields due to the vertical diffusion
      !!      of tracers alone is given by:
      !!                zwx = tb + p2dt difft
      !!      where difft = dz( avt dz(tb) ) = 1/e3t dk+1( avt/e3w dk(tb) )
      !!           (if lk_zdfddm=T use avs on salinity instead of avt)
      !!      difft is evaluated with an Euler split-explit scheme using a
      !!      no flux boundary condition at both surface and bottomi boundaries.
      !!      (N.B. bottom condition is applied through the masked field avt).
      !!              - the after tracer fields due to the whole trend is
      !!      obtained in leap-frog environment by :
      !!          ta = zwx + p2dt ta
      !!              - in case of variable level thickness (lk_vvl=T) the
      !!     the leap-frog is applied on thickness weighted tracer. That is:
      !!          ta = [ tb*e3tb + e3tn*( zwx - tb + p2dt ta ) ] / e3tn
      !!
      !! ** Action : - after tracer fields (ta,sa)
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in)                 ::   kt     ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000      ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype      ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt        ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kn_zdfexp   ! number of sub-time step
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt        ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb_tl      ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_tl      ! tracer trend
      !!
      INTEGER  ::   ji, jj, jk, jl, jn            ! dummy loop indices
      REAL(wp) ::   zlavmr, zave3r, ze3tr     ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwxtl, zwytl   ! 3D workspace
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_exp_tan')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwxtl, zwytl )
      !
      IF( kt == kit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_zdf_exp_tan : explicit vertical mixing on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      ! Initializations
      ! ---------------
      zlavmr = 1. / float( kn_zdfexp )                           ! Local constant
      !
      Do jn = 1, kjpt
         zwytl(:,:, 1 ) = 0.0_wp                ! surface boundary conditions: no flux
         zwytl(:,:,jpk) = 0.0_wp                ! bottom  boundary conditions: no flux
         !
         zwxtl(:,:,:)   = ptb_tl(:,:,:,jn)      ! zwx and zwz arrays set to before tracer values

         ! Split-explicit loop  (after tracer due to the vertical diffusion alone)
         ! -------------------
         !
         DO jl = 1, kn_zdfexp
            !                     ! first vertical derivative
            DO jk = 2, jpk
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     zave3r = 1.e0 / e3w(ji,jj,jk)
                     IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN  ! temperature : use of avt
                        zwytl(ji,jj,jk) =   avt(ji,jj,jk) * ( zwxtl(ji,jj,jk-1) - zwxtl(ji,jj,jk) ) * zave3r
                     ELSE
                        zwytl(ji,jj,jk) = avt(ji,jj,jk) * ( zwxtl(ji,jj,jk-1) - zwxtl(ji,jj,jk) ) * zave3r
                     END IF
                  END DO
               END DO
            END DO
            !
            DO jk = 1, jpkm1      ! second vertical derivative   ==> tracer at kt+l*2*rdt/n_zdfexp
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ze3tr = zlavmr / e3t(ji,jj,jk)
                     zwxtl(ji,jj,jk) = zwxtl(ji,jj,jk) + p2dt(jk) * ( zwytl(ji,jj,jk) - zwytl(ji,jj,jk+1) ) * ze3tr
                  END DO
               END DO
            END DO
            !
         END DO

         ! After tracer due to all trends
         ! ------------------------------
         IF( lk_vvl ) THEN          ! variable level thickness : leap-frog on tracer*e3t
            IF(lwp) WRITE(numout,*) "key_vvl net available in tangent yet"
            CALL abort
         ELSE                       ! fixed level thickness : leap-frog on tracers
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     pta_tl(ji,jj,jk,jn) = ( zwxtl(ji,jj,jk) + p2dt(jk) * pta_tl(ji,jj,jk,jn) ) *tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF
         !
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwxtl, zwytl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_exp_tan')
      !
   END SUBROUTINE tra_zdf_exp_tan

   SUBROUTINE tra_zdf_exp_adj( kt, kit000, cdtype, p2dt, kn_zdfexp,   &
      &                                ptb_ad , pta_ad      , kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_exp_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!      Compute the after tracer fields due to the vertical
      !!      tracer mixing alone, and then due to the whole tracer trend.
      !!
      !! ** Method of the direct routine :
      !!               - The after tracer fields due to the vertical diffusion
      !!      of tracers alone is given by:
      !!                zwx = tb + p2dt difft
      !!      where difft = dz( avt dz(tb) ) = 1/e3t dk+1( avt/e3w dk(tb) )
      !!           (if lk_zdfddm=T use avs on salinity instead of avt)
      !!      difft is evaluated with an Euler split-explit scheme using a
      !!      no flux boundary condition at both surface and bottomi boundaries.
      !!      (N.B. bottom condition is applied through the masked field avt).
      !!              - the after tracer fields due to the whole trend is
      !!      obtained in leap-frog environment by :
      !!          ta = zwx + p2dt ta
      !!              - in case of variable level thickness (lk_vvl=T) the
      !!     the leap-frog is applied on thickness weighted tracer. That is:
      !!          ta = [ tb*e3tb + e3tn*( zwx - tb + p2dt ta ) ] / e3tn
      !!
      !! ** Action : - after tracer fields (ta,sa)
      !!---------------------------------------------------------------------
      INTEGER , INTENT(in)                 ::   kt     ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000      ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype      ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt        ! number of tracers
      INTEGER                              , INTENT(in   ) ::   kn_zdfexp   ! number of sub-time step
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt        ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   ptb_ad      ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_ad      ! tracer trend
      !!
      INTEGER  ::   ji, jj, jk, jl, jn            ! dummy loop indices
      REAL(wp) ::   zlavmr, zave3r, ze3tr     ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zwxad, zwyad                 ! 3D workspace
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_exp_adj')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwxad, zwyad )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_zdf_exp_adj : explicit vertical mixing on ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      ENDIF

      ! Initializations
      ! ---------------
      zlavmr = 1. / float( kn_zdfexp )                           ! Local constant
      DO jn = 1, kjpt
         !
         zwxad(:,:,:) = 0.0_wp
         zwyad(:,:,:) = 0.0_wp
         ! After tracer due to all trends
         ! ------------------------------
         IF( lk_vvl ) THEN          ! variable level thickness : leap-frog on tracer*e3t
            IF(lwp) WRITE(numout,*) "key_vvl net available in adjoint yet"
            CALL abort
         ELSE                       ! fixed level thickness : leap-frog on tracers
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     zwxad(ji,jj,jk) = zwxad(ji,jj,jk) + pta_ad(ji,jj,jk,jn) * tmask(ji,jj,jk)
                     pta_ad(ji,jj,jk,jn) = p2dt(jk) * pta_ad(ji,jj,jk,jn) * tmask(ji,jj,jk)
                  END DO
               END DO
            END DO
         ENDIF
         !

         ! Split-explicit loop  (after tracer due to the vertical diffusion alone)
         ! -------------------
         !
         DO jl = 1, kn_zdfexp
            DO jk =  jpkm1, 1, -1      ! second vertical derivative   ==> tracer at kt+l*2*rdt/n_zdfexp
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ze3tr = zlavmr / e3t(ji,jj,jk)
                     zwyad(ji,jj,jk  ) = zwyad(ji,jj,jk  ) + p2dt(jk) * zwxad(ji,jj,jk) * ze3tr
                     zwyad(ji,jj,jk+1) = zwyad(ji,jj,jk+1) - p2dt(jk) * zwxad(ji,jj,jk) * ze3tr
                  END DO
               END DO
            END DO
            !                     ! first vertical derivative
            DO jk = jpk, 2, -1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     zave3r = 1.e0 / e3w(ji,jj,jk)
                     IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN  ! temperature : use of avt
                        zwxad(ji,jj,jk-1) = zwxad(ji,jj,jk-1) + avt(ji,jj,jk) * zwyad(ji,jj,jk) * zave3r
                        zwxad(ji,jj,jk  ) = zwxad(ji,jj,jk  ) - avt(ji,jj,jk) * zwyad(ji,jj,jk) * zave3r
                        zwyad(ji,jj,jk  ) = 0.0_wp
                     ELSE
                        zwxad(ji,jj,jk-1) = zwxad(ji,jj,jk-1) + avt(ji,jj,jk) * zwyad(ji,jj,jk) * zave3r
                        zwxad(ji,jj,jk  ) = zwxad(ji,jj,jk  ) - avt(ji,jj,jk) * zwyad(ji,jj,jk) * zave3r
                        zwyad(ji,jj,jk  ) = 0.0_wp
                     ENDIF
                  END DO
               END DO
            END DO
            !
            !
         END DO
         !
         ptb_ad(:,:,:,jn) = ptb_ad(:,:,:,jn) + zwxad(:,:,:)
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwxad, zwyad )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_exp_adj')
      !
   END SUBROUTINE tra_zdf_exp_adj

   SUBROUTINE tra_zdf_exp_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE tra_zdf_exp_adj_tst ***
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
      REAL(KIND=wp) ::   &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zta_tlin ,     & ! Tangent input
         & ztb_tlin ,     & ! Tangent input
         & zsa_tlin ,     & ! Tangent input
         & zsb_tlin ,     & ! Tangent input
         & zta_tlout,     & ! Tangent output
         & zsa_tlout,     & ! Tangent output
         & zta_adin ,     & ! Adjoint input
         & zsa_adin ,     & ! Adjoint input
         & zta_adout,     & ! Adjoint output
         & ztb_adout,     & ! Adjoint output
         & zsa_adout,     & ! Adjoint output
         & zsb_adout,     & ! Adjoint output
         & zr             ! 3D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & zta_tlin( jpi,jpj,jpk),     &
         & zsa_tlin( jpi,jpj,jpk),     &
         & ztb_tlin( jpi,jpj,jpk),     &
         & zsb_tlin( jpi,jpj,jpk),     &
         & zta_tlout(jpi,jpj,jpk),     &
         & zsa_tlout(jpi,jpj,jpk),     &
         & zta_adin( jpi,jpj,jpk),     &
         & zsa_adin( jpi,jpj,jpk),     &
         & zta_adout(jpi,jpj,jpk),     &
         & zsa_adout(jpi,jpj,jpk),     &
         & ztb_adout(jpi,jpj,jpk),     &
         & zsb_adout(jpi,jpj,jpk),     &
         & zr(       jpi,jpj,jpk)      &
         & )


      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          zta_tlin( :,:,:) = 0.0_wp
          ztb_tlin( :,:,:) = 0.0_wp
          zsa_tlin( :,:,:) = 0.0_wp
          zsb_tlin( :,:,:) = 0.0_wp
          zta_tlout(:,:,:) = 0.0_wp
          zsa_tlout(:,:,:) = 0.0_wp
          zta_adin( :,:,:) = 0.0_wp
          zsa_adin( :,:,:) = 0.0_wp
          zta_adout(:,:,:) = 0.0_wp
          zsa_adout(:,:,:) = 0.0_wp
          ztb_adout(:,:,:) = 0.0_wp
          zsb_adout(:,:,:) = 0.0_wp
          zr(       :,:,:) = 0.0_wp

          tsa_tl(:,:,:,:)     = 0.0_wp
          tsb_tl(:,:,:,:)     = 0.0_wp
          tsa_ad(:,:,:,:)     = 0.0_wp
          tsb_ad(:,:,:,:)     = 0.0_wp

          r2dtra(:) =  2.* rdttra(:)
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------
      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_tlin(ji,jj,jk) = zr(ji,jj,jk)
           END DO
        END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              ztb_tlin(ji,jj,jk) = zr(ji,jj,jk)
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
      CALL grid_random(  zr, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsb_tlin(ji,jj,jk) = zr(ji,jj,jk)
           END DO
        END DO
      END DO
      tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
      tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)
      tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
      tsb_tl(:,:,:,jp_sal) = zsb_tlin(:,:,:)
      CALL tra_zdf_exp_tan ( nit000, nit000, 'TRA', r2dtra, 1, tsb_tl, tsa_tl, jpts )
      zta_tlout(:,:,:) = tsa_tl(:,:,:,jp_tem)
      zsa_tlout(:,:,:) = tsa_tl(:,:,:,jp_sal)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_adin(ji,jj,jk) = zta_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_t(jk)
              zsa_adin(ji,jj,jk) = zsa_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_s(jk)
            END DO
         END DO
      END DO
      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 = DOT_PRODUCT( zta_tlout, zta_adin ) &
         & + DOT_PRODUCT( zsa_tlout, zsa_adin )


      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      tsa_ad(:,:,:,jp_tem) = zta_adin(:,:,:)
      tsa_ad(:,:,:,jp_sal) = zsa_adin(:,:,:)


      CALL tra_zdf_exp_adj ( nit000, nit000, 'TRA', r2dtra, 1, tsb_ad, tsa_ad, jpts )

      zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)
      zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)
      ztb_adout(:,:,:) = tsb_ad(:,:,:,jp_tem)
      zsb_adout(:,:,:) = tsb_ad(:,:,:,jp_sal)

      zsp2 = DOT_PRODUCT( zta_tlin, zta_adout ) &
         & + DOT_PRODUCT( zsa_tlin, zsa_adout ) &
         & + DOT_PRODUCT( ztb_tlin, ztb_adout ) &
         & + DOT_PRODUCT( zsb_tlin, zsb_adout )

      ! 14 char:'12345678901234'
      cl_name = 'trazdf_exp_adj'
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(   &
         & zta_tlin,  &
         & ztb_tlin,  &
         & zsa_tlin,  &
         & zsb_tlin,  &
         & zta_tlout, &
         & zsa_tlout, &
         & zta_adin,  &
         & zsa_adin,  &
         & zta_adout, &
         & ztb_adout, &
         & zsa_adout, &
         & zsb_adout, &
         & zr       &
         & )



   END SUBROUTINE tra_zdf_exp_adj_tst

   !!==============================================================================
END MODULE trazdf_exp_tam
