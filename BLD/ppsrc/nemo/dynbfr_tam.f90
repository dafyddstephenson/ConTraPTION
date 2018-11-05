MODULE dynbfr_tam
   !!==============================================================================
   !!                 ***  MODULE  dynbfr  ***
   !! Ocean dynamics :  bottom friction component of the momentum mixing trend
   !!==============================================================================
   !! History of the drect module:
   !!            9.0  !  2008-11  (A. C. Coward)  Original code
   !! History of the TAM module:
   !!  NEMO      3.2  ! 2010-04 (F. Vigilant) Original code
   !!  NEMO      3.4  ! 2012-07 (P.-A. bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_bfr      : Update the momentum trend with the bottom friction contribution
   !!----------------------------------------------------------------------
   USE oce
   USE oce_tam
   USE par_oce
   USE dom_oce
   USE zdf_oce         ! ocean vertical physics variables
   USE zdfbfr
   USE zdf_oce_tam
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE timing
   USE wrk_nemo

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_bfr_tan     ! routine called by step_tam.F90
   PUBLIC   dyn_bfr_adj     ! routine called by step_tam.F90
   PUBLIC   dyn_bfr_adj_tst ! routine called by the tst.F90

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

CONTAINS

   SUBROUTINE dyn_bfr_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_bfr_tan  ***
      !!
      !! ** Purpose of direct routine:
      !!                compute the bottom friction ocean dynamics physics.
      !!
      !! ** Action  :   (ua,va)   momentum trend increased by bottom friction trend
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj          ! dummy loop indexes
      INTEGER  ::   ikbu  , ikbv    ! temporary integers
      REAL(wp) ::   zm1_2dt, zbfru, zbfrv   ! temporary scalar
      REAL(wp) ::   zbfrutl, zbfrvtl     ! temporary scalar
      !!---------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_bfr_tan')
      !
      IF ( .NOT. ln_bfrimp ) THEN
         zm1_2dt = -1._wp / ( 2._wp * rdt )
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ikbu = mbku(ji,jj)
               ikbv = mbkv(ji,jj)
               !
               ! Apply stability criteria on absolute value  : Min abs(bfr) => Max (bfr)
               IF ( bfrua(ji,jj) >= e3u(ji,jj,ikbu)*zm1_2dt ) THEN
                  zbfru   = bfrua(   ji,jj)
                  zbfrutl  = bfrua_tl(ji,jj)
               ELSE
                  zbfru   = e3u(ji,jj,ikbu)*zm1_2dt
                  zbfrutl  = 0.0_wp
               END IF

               IF ( bfrva(ji,jj) >= e3v(ji,jj,ikbv)*zm1_2dt ) THEN
                  zbfrv   = bfrva(   ji,jj)
                  zbfrvtl = bfrva_tl(ji,jj)
               ELSE
                  zbfrv   = e3v(ji,jj,ikbv)*zm1_2dt
                  zbfrvtl = 0.0_wp
               END IF
               !
               ua_tl(ji,jj,ikbu) = ua_tl(ji,jj,ikbu) &
                  &                + ( zbfru * ub_tl(ji,jj,ikbu) + zbfrutl * ub(ji,jj,ikbu) ) &
                  &                / e3u(ji,jj,ikbu)
               va_tl(ji,jj,ikbv) = va_tl(ji,jj,ikbv) &
                  &                + ( zbfrv * vb_tl(ji,jj,ikbv) + zbfrvtl * vb(ji,jj,ikbv) ) &
                  &                / e3v(ji,jj,ikbv)
               !
            END DO
         END DO
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_bfr_tan')
      !
   END SUBROUTINE dyn_bfr_tan

   SUBROUTINE dyn_bfr_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_bfr_adj  ***
      !!
      !! ** Purpose of direct routine:
      !!                compute the bottom friction ocean dynamics physics.
      !!
      !! ** Action  :   (ua,va)   momentum trend increased by bottom friction trend
      !!---------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj          ! dummy loop indexes
      INTEGER  ::   ikbu  , ikbv    ! temporary integers
      REAL(wp) ::   zm1_2dt, zbfru, zbfrv   ! temporary scalar
      REAL(wp) ::   zbfruad, zbfrvad     ! temporary scalar
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_bfr_adj')
      !
      IF( .NOT.ln_bfrimp) THEN
         zm1_2dt = -1._wp / ( 2._wp * rdt )
         zbfruad = 0.0_wp  ;  zbfrvad = 0.0_wp
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1
               ikbu = mbku(ji,jj)
               ikbv = mbkv(ji,jj)
               !
               ! Apply stability criteria on absolute value  : Min abs(bfr) => Max (bfr)
               zbfru = MAX( bfrua(ji,jj), e3u(ji,jj,ikbu)*zm1_2dt )
               zbfrv = MAX( bfrva(ji,jj), e3v(ji,jj,ikbv)*zm1_2dt )
               !
               ub_ad(ji,jj,ikbu) = ub_ad(ji,jj,ikbu) + zbfru * ua_ad(ji,jj,ikbu) / e3u(ji,jj,ikbu)
               zbfruad = zbfruad + ub(ji,jj,ikbu) * ua_ad(ji,jj,ikbu) / e3u(ji,jj,ikbu)
               vb_ad(ji,jj,ikbv) = vb_ad(ji,jj,ikbv) + zbfrv * va_ad(ji,jj,ikbv) / e3v(ji,jj,ikbv)
               zbfrvad = zbfrvad + vb(ji,jj,ikbv) * va_ad(ji,jj,ikbv) / e3v(ji,jj,ikbv)
               !
               ! Apply stability criteria on absolute value  : Min abs(bfr) => Max (bfr)
               IF ( bfrua(ji,jj) >= e3u(ji,jj,ikbu)*zm1_2dt ) THEN
                  bfrua_ad(ji,jj) = bfrua_ad(ji,jj) + zbfruad
                  zbfruad  = 0.0_wp
               ELSE
                  zbfruad  = 0.0_wp
               END IF
               !
               IF ( bfrva(ji,jj) >= e3v(ji,jj,ikbv)*zm1_2dt ) THEN
                  bfrva_ad(ji,jj) = bfrva_ad(ji,jj) + zbfrvad
                  zbfrvad = 0.0_wp
               ELSE
                  zbfrvad = 0.0_wp
               END IF
               !
            END DO
         END DO
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_bfr_adj')
      !
   END SUBROUTINE dyn_bfr_adj

   SUBROUTINE dyn_bfr_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_bfr_adj_tst ***
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
      !! ** Action  :
      !!
      !! History :
      !!        ! 2010-04 (F. Vigilant)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt        ! Output unit

      !! * Local declarations
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zua_tlin,    & ! Tangent input: ua_tl
         & zva_tlin,    & ! Tangent input: va_tl
         & zub_tlin,    & ! Tangent input: ub_tl
         & zvb_tlin,    & ! Tangent input: vb_tl
         & zua_tlout,   & ! Tangent output: ua_tl
         & zva_tlout,   & ! Tangent output: va_tl
         & zua_adin,    & ! Adjoint input: ua_ad
         & zva_adin,    & ! Adjoint input: va_ad
         & zua_adout,   & ! Adjoint output: ua_ad
         & zva_adout,   & ! Adjoint output: va_ad
         & zub_adout,   & ! Adjoint oputput: ub_ad
         & zvb_adout,   & ! Adjoint output: vb_ad
         & znu            ! 3D random field for u
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zbu_tlin,     & ! Tangent input: bfrua_tl
         & zbv_tlin,     & ! Tangent input: bfrva_tl
         & zbu_adout,    & ! Adjoint output: bfrua_ad
         & zbv_adout       ! Adjoint output: bfrva_ad
      REAL(wp) :: &
         & zsp1,    &   ! scalar product involving the tangent routine
         & zsp2         ! scalar product involving the adjoint routine
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d    ! 2D seed for the random number generator
      INTEGER :: &
         & ji, &
         & jj, &
         & jk
      CHARACTER (LEN=14) :: &
         & cl_name

      ! Allocate memory

      ALLOCATE( &
         & zua_tlin(jpi,jpj,jpk),  &
         & zva_tlin(jpi,jpj,jpk),  &
         & zub_tlin(jpi,jpj,jpk),  &
         & zvb_tlin(jpi,jpj,jpk),  &
         & zua_tlout(jpi,jpj,jpk), &
         & zva_tlout(jpi,jpj,jpk), &
         & zua_adin(jpi,jpj,jpk),  &
         & zva_adin(jpi,jpj,jpk),  &
         & zua_adout(jpi,jpj,jpk), &
         & zva_adout(jpi,jpj,jpk), &
         & zub_adout(jpi,jpj,jpk), &
         & zvb_adout(jpi,jpj,jpk), &
         & znu(jpi,jpj,jpk)        &
         & )
      ALLOCATE( &
         & zbu_tlin(jpi,jpj),      &
         & zbv_tlin(jpi,jpj),      &
         & zbu_adout(jpi,jpj),     &
         & zbv_adout(jpi,jpj)      &
         & )

      !=========================================================================
      !     dx = ( ub_tl, ua_tl, vb_tl, va_tl )
      ! and dy = ( ua_tl, va_tl )
      !=========================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      zub_tlin (:,:,:) = 0.0_wp
      zvb_tlin (:,:,:) = 0.0_wp
      zua_tlin (:,:,:) = 0.0_wp
      zva_tlin (:,:,:) = 0.0_wp
      zua_tlout(:,:,:) = 0.0_wp
      zva_tlout(:,:,:) = 0.0_wp
      zua_adin (:,:,:) = 0.0_wp
      zva_adin (:,:,:) = 0.0_wp
      zub_adout(:,:,:) = 0.0_wp
      zvb_adout(:,:,:) = 0.0_wp
      zua_adout(:,:,:) = 0.0_wp
      zva_adout(:,:,:) = 0.0_wp

      ub_tl(:,:,:) = 0.0_wp
      vb_tl(:,:,:) = 0.0_wp
      ua_tl(:,:,:) = 0.0_wp
      va_tl(:,:,:) = 0.0_wp
      ub_ad(:,:,:) = 0.0_wp
      vb_ad(:,:,:) = 0.0_wp
      ua_ad(:,:,:) = 0.0_wp
      va_ad(:,:,:) = 0.0_wp
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  znu, 'U', 0.0_wp, stdu )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zua_tlin(ji,jj,jk) = znu(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  znu, 'V', 0.0_wp, stdv )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
              zva_tlin(ji,jj,jk) = znu(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  znu, 'U', 0.0_wp, stdu )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
             zub_tlin(ji,jj,jk) = znu(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  znu, 'V', 0.0_wp, stdv )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
             zvb_tlin(ji,jj,jk) = znu(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  znu, 'U', 0.0_wp, stdv )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zbu_tlin(ji,jj) = znu(ji,jj,1)
         END DO
      END DO
      CALL grid_random(  znu, 'V', 0.0_wp, stdv )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zbv_tlin(ji,jj) = znu(ji,jj,1)
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      ua_tl(:,:,:) = zua_tlin(:,:,:)
      va_tl(:,:,:) = zva_tlin(:,:,:)
      ub_tl(:,:,:) = zub_tlin(:,:,:)
      vb_tl(:,:,:) = zvb_tlin(:,:,:)
      bfrua_tl(:,:) = zbu_tlin(:,:)
      bfrva_tl(:,:) = zbv_tlin(:,:)

      CALL dyn_bfr_tan( nit000 )

      zua_tlout(:,:,:) = ua_tl(:,:,:)
      zva_tlout(:,:,:) = va_tl(:,:,:)
      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zua_adin(ji,jj,jk) = zua_tlout(ji,jj,jk) &
                 &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
                 &               * umask(ji,jj,jk)
              zva_adin(ji,jj,jk) = zva_tlout(ji,jj,jk) &
                 &               * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) &
                 &               * vmask(ji,jj,jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 =   DOT_PRODUCT( zua_tlout  , zua_adin   ) + DOT_PRODUCT( zva_tlout  , zva_adin   )


      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      ua_ad(:,:,:) = zua_adin(:,:,:)
      va_ad(:,:,:) = zva_adin(:,:,:)

      CALL dyn_bfr_adj( nit000 )

      zua_adout(:,:,:) = ua_ad(:,:,:)
      zva_adout(:,:,:) = va_ad(:,:,:)
      zub_adout(:,:,:) = ub_ad(:,:,:)
      zvb_adout(:,:,:) = vb_ad(:,:,:)
      zbu_adout(:,:) = bfrua_ad(:,:)
      zbv_adout(:,:) = bfrva_ad(:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2 =   DOT_PRODUCT( zua_tlin  , zua_adout   ) &
            &   + DOT_PRODUCT( zva_tlin  , zva_adout   ) &
            &   + DOT_PRODUCT( zub_tlin  , zub_adout   ) &
            &   + DOT_PRODUCT( zvb_tlin  , zvb_adout   ) &
            &   + DOT_PRODUCT( zbu_tlin  , zbu_adout   ) &
            &   + DOT_PRODUCT( zbv_tlin  , zbv_adout   )

      ! Compare the scalar products

      ! 14 char:'12345678901234'
      cl_name = 'dyn_bfr_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      ! Deallocate memory

      DEALLOCATE( &
         & zua_tlin,  &
         & zva_tlin,  &
         & zub_tlin,  &
         & zvb_tlin,  &
         & zua_tlout, &
         & zva_tlout, &
         & zua_adin,  &
         & zva_adin,  &
         & zua_adout, &
         & zva_adout, &
         & zub_adout, &
         & zvb_adout, &
         & znu        &
         & )
      DEALLOCATE( &
         & zbu_tlin,      &
         & zbv_tlin,      &
         & zbu_adout,     &
         & zbv_adout      &
         & )

   END SUBROUTINE dyn_bfr_adj_tst
   !!==============================================================================
END MODULE dynbfr_tam
