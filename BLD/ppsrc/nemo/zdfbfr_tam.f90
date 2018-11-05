MODULE zdfbfr_tam
   !!======================================================================
   !!                       ***  MODULE  zdfbfr_tam  ***
   !! Ocean physics: Bottom friction
   !!======================================================================
   !! History of the direct module:
   !!            OPA  ! 1997-06  (G. Madec, A.-M. Treguier)  Original code
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            3.2  ! 2009-09  (A.C.Coward)  Correction to include barotropic contribution
   !! History of the T&A module:
   !!   NEMO     3.2.2! 2011-02  (A. Vidard) Original version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   zdf_bfr_tan  : update momentum Kz at the ocean bottom due to the type of bottom friction chosen
   !!   zdf_bfr_adj  : update momentum Kz at the ocean bottom due to the type of bottom friction chosen
   !!                  parameters.
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE zdf_oce         ! ocean vertical physics variables
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! distributed memory computing
   USE zdfbfr
   USE oce_tam
   USE lbclnk_tam
   USE timing
   USE zdf_oce_tam
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE tstool_tam

   IMPLICIT NONE
   PRIVATE

   PUBLIC   zdf_bfr_tan    ! called by step_tam.F90
   PUBLIC   zdf_bfr_adj    ! called by step_tam.F90
   PUBLIC   zdf_bfr_adj_tst
   PUBLIC   zdf_bfr_init_tam

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
   !! NEMO/OPA 3.2 , LOCEAN-IPSL (2009)
   !! $Id$
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   SUBROUTINE zdf_bfr_init_tam
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE zdf_bfr_init  ***
      !!
      !! ** Purpose :   Initialization of the bottom friction
      !!
      !! ** Method  :   Read the nammbf namelist and check their consistency
      !!              called at the first timestep (nit000)
      !!----------------------------------------------------------------------
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_bfr_init_tam')
      !
      bfrua_tl = 0._wp
      bfrva_tl = 0._wp
      bfrua_ad = 0._wp
      bfrva_ad = 0._wp

      IF( nn_timing == 1 )  CALL timing_stop('zdf_bfr_init_tam')
      !
   END SUBROUTINE zdf_bfr_init_tam
   SUBROUTINE zdf_bfr_tan( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_bfr_tan  ***
      !!
      !! ** Purpose :   tangent of the computation of the bottom friction coefficient.
      !!
      !! ** Method  :   Calculate and store part of the momentum trend due
      !!              to bottom friction following the chosen friction type
      !!              (free-slip, linear, or quadratic). The component
      !!              calculated here is multiplied by the bottom velocity in
      !!              dyn_bfr to provide the trend term.
      !!                The coefficients are updated at each time step only
      !!              in the quadratic case.
      !!
      !! ** Action  :   bfrua , bfrva   bottom friction coefficients
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj         ! dummy loop indices
      INTEGER  ::   ikbu           ! temporary integers
      INTEGER  ::   ikbv           !    -          -
      REAL(wp) ::   zvu, zuv, zecu, zecv   ! temporary scalars
      REAL(wp) ::   zvutl, zuvtl, zecutl, zecvtl   ! temporary scalars
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_bfr_tan')
      !
      IF( kt == nit000 )   THEN
         bfrua_tl = 0.0_wp
         bfrva_tl = 0.0_wp
      END IF
      IF( nn_bfr == 2 ) THEN                 ! quadratic botton friction
         ! Calculate and store the quadratic bottom friction coefficient bfrua and bfrva
         ! where bfrUa = C_d*SQRT(u_bot^2 + v_bot^2 + e_b) {U=[u,v]}
         ! from these the trend due to bottom friction:  -F_h/e3U  can be calculated
         ! where -F_h/e3U_bot = bfrUa*Ub/e3U_bot {U=[u,v]}
         !
!CDIR NOVERRCHK
         DO jj = 2, jpjm1
!CDIR NOVERRCHK
            DO ji = 2, jpim1
               ikbu   = mbku(ji,jj)
               ikbv   = mbkv(ji,jj)
               !
               zvu   = 0.25 * (  vn(ji,jj  ,ikbu) + vn(ji+1,jj  ,ikbu)     &
                  &            + vn(ji,jj-1,ikbu) + vn(ji+1,jj-1,ikbu)  )
               zuv   = 0.25 * (  un(ji,jj  ,ikbv) + un(ji-1,jj  ,ikbv)     &
                  &            + un(ji,jj+1,ikbv) + un(ji-1,jj+1,ikbv)  )
               zvutl = 0.25 * (  vn_tl(ji,jj  ,ikbu) + vn_tl(ji+1,jj  ,ikbu)     &
                  &           + vn_tl(ji,jj-1,ikbu) + vn_tl(ji+1,jj-1,ikbu)  )
               zuvtl = 0.25 * (  un_tl(ji,jj  ,ikbv) + un_tl(ji-1,jj  ,ikbv)     &
                  &           + un_tl(ji,jj+1,ikbv) + un_tl(ji-1,jj+1,ikbv)  )
               !
               zecu = SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_bfeb2  )
               zecv = SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_bfeb2  )
               zecutl = ( un(ji,jj,ikbu) * un_tl(ji,jj,ikbu) + zvu*zvutl )  &
                  &     / SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_bfeb2  )
               zecvtl = ( vn(ji,jj,ikbv) * vn_tl(ji,jj,ikbv) + zuv*zuvtl )  &
                  &     / SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_bfeb2  )
               !
               bfrua_tl(ji,jj) = - 0.5_wp * ( bfrcoef2d(ji,jj) + bfrcoef2d(ji+1,jj  ) ) * zecutl
               bfrva_tl(ji,jj) = - 0.5_wp * ( bfrcoef2d(ji,jj) + bfrcoef2d(ji  ,jj+1) ) * zecvtl
            END DO
         END DO
         !
         CALL lbc_lnk( bfrua_tl, 'U', 1. )   ;   CALL lbc_lnk( bfrva_tl, 'V', 1. ) ! Lateral boundary condition
         !
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_bfr_tan')
      !
   END SUBROUTINE zdf_bfr_tan
   SUBROUTINE zdf_bfr_adj( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE zdf_bfr_adj  ***
      !!
      !! ** Purpose :   adjoint of the computation of the bottom friction coefficient.
      !!
      !! ** Method  :   Calculate and store part of the momentum trend due
      !!              to bottom friction following the chosen friction type
      !!              (free-slip, linear, or quadratic). The component
      !!              calculated here is multiplied by the bottom velocity in
      !!              dyn_bfr to provide the trend term.
      !!                The coefficients are updated at each time step only
      !!              in the quadratic case.
      !!
      !! ** Action  :   bfrua , bfrva   bottom friction coefficients
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj         ! dummy loop indices
      INTEGER  ::   ikbu, ikbum1   ! temporary integers
      INTEGER  ::   ikbv, ikbvm1   !    -          -
      REAL(wp) ::   zvu, zuv, zecu, zecv   ! temporary scalars
      REAL(wp) ::   zvuad, zuvad, zecuad, zecvad   ! temporary scalars
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('zdf_bfr_adj')
      !
      zvuad = 0.0_wp   ;   zuvad = 0.0_wp   ;   zecuad = 0.0_wp   ;   zecvad = 0.0_wp

      IF( kt == nitend )   THEN
         bfrua_ad = 0.0_wp
         bfrva_ad = 0.0_wp
      END IF
      IF( nn_bfr == 2 ) THEN                 ! quadratic botton friction
         ! Calculate and store the quadratic bottom friction coefficient bfrua and bfrva
         ! where bfrUa = C_d*SQRT(u_bot^2 + v_bot^2 + e_b) {U=[u,v]}
         ! from these the trend due to bottom friction:  -F_h/e3U  can be calculated
         ! where -F_h/e3U_bot = bfrUa*Ub/e3U_bot {U=[u,v]}
         !
         CALL lbc_lnk_adj( bfrua_ad, 'U', 1. )   ;   CALL lbc_lnk_adj( bfrva_ad, 'V', 1. ) ! Lateral boundary condition
         !
!CDIR NOVERRCHK
         DO jj = 2, jpjm1
!CDIR NOVERRCHK
            DO ji = 2, jpim1
               ikbu   = mbku(ji,jj)
               ikbv   = mbkv(ji,jj)
               ! direct computation
               zvu   = 0.25 * (  vn(ji,jj  ,ikbu) + vn(ji+1,jj  ,ikbu)     &
                  &            + vn(ji,jj-1,ikbu) + vn(ji+1,jj-1,ikbu)  )
               zuv   = 0.25 * (  un(ji,jj  ,ikbv) + un(ji-1,jj  ,ikbv)     &
                  &            + un(ji,jj+1,ikbv) + un(ji-1,jj+1,ikbv)  )
               zecu = SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_bfeb2  )
               zecv = SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_bfeb2  )
               ! Adjoint counterpart

               zecuad = - 0.5_wp * ( bfrcoef2d(ji,jj) + bfrcoef2d(ji+1,jj  ) ) * bfrua_ad(ji,jj)
               bfrua_ad(ji,jj) = 0.0_wp
               zecvad = - 0.5_wp * ( bfrcoef2d(ji,jj) + bfrcoef2d(ji  ,jj+1) ) * bfrva_ad(ji,jj)
               bfrva_ad(ji,jj) = 0.0_wp
               !
               un_ad(ji,jj,ikbu) = un_ad(ji,jj,ikbu) + zecuad * un(ji,jj,ikbu) &
                  &                  / SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_bfeb2  )
               zvuad  = zecuad * zvu / SQRT(  un(ji,jj,ikbu) * un(ji,jj,ikbu) + zvu*zvu + rn_bfeb2  )

               vn_ad(ji,jj,ikbv) = vn_ad(ji,jj,ikbv) + zecvad * vn(ji,jj,ikbv) &
                  &                  / SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_bfeb2  )
               zuvad  = zecvad * zuv / SQRT(  vn(ji,jj,ikbv) * vn(ji,jj,ikbv) + zuv*zuv + rn_bfeb2  )
               !
               vn_ad(ji  ,jj  ,ikbu) = vn_ad(ji,jj    ,ikbu) + zvuad * 0.25
               vn_ad(ji+1,jj  ,ikbu) = vn_ad(ji+1,jj  ,ikbu) + zvuad * 0.25
               vn_ad(ji  ,jj-1,ikbu) = vn_ad(ji,jj-1  ,ikbu) + zvuad * 0.25
               vn_ad(ji+1,jj-1,ikbu) = vn_ad(ji+1,jj-1,ikbu) + zvuad * 0.25

               un_ad(ji  ,jj  ,ikbv) = un_ad(ji  ,jj  ,ikbv) + zuvad * 0.25
               un_ad(ji-1,jj  ,ikbv) = un_ad(ji-1,jj  ,ikbv) + zuvad * 0.25
               un_ad(ji  ,jj+1,ikbv) = un_ad(ji  ,jj+1,ikbv) + zuvad * 0.25
               un_ad(ji-1,jj+1,ikbv) = un_ad(ji-1,jj+1,ikbv) + zuvad * 0.25
               !
            END DO
         END DO
         !
      ENDIF
      !
      IF ( kt == nit000 ) THEN
         bfrua_ad = 0.0_wp
         bfrva_ad = 0.0_wp
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop('zdf_bfr_adj')
      !
   END SUBROUTINE zdf_bfr_adj
   SUBROUTINE zdf_bfr_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE zdf_bfr_adj_tst ***
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
      !!        ! 09-01 (A. Weaver)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt        ! Output unit

      !! * Local declarations
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zua_tlin,    & ! Tangent input: ua_tl
         & zva_tlin,    & ! Tangent input: va_tl
         & zua_tlout,   & ! Tangent output: ua_tl
         & zva_tlout,   & ! Tangent output: va_tl
         & zua_adin,    & ! Adjoint input: ua_ad
         & zva_adin,    & ! Adjoint input: va_ad
         & zua_adout,   & ! Adjoint output: ua_ad
         & zva_adout,   & ! Adjoint output: va_ad
         & znu            ! 3D random field for u

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zspgu_tlout, zspgv_tlout, zspgu_adin, zspgv_adin

      REAL(wp) :: &
         & zsp1,    &   ! scalar product involving the tangent routine
         & zsp2         ! scalar product involving the adjoint routine
      INTEGER :: &
         & ji, &
         & jj, &
         & jk
      CHARACTER (LEN=14) :: &
         & cl_name

      ALLOCATE( &
         & zua_tlin(jpi,jpj,jpk),  &
         & zva_tlin(jpi,jpj,jpk),  &
         & zua_tlout(jpi,jpj,jpk), &
         & zva_tlout(jpi,jpj,jpk), &
         & zua_adin(jpi,jpj,jpk),  &
         & zva_adin(jpi,jpj,jpk),  &
         & zua_adout(jpi,jpj,jpk), &
         & zva_adout(jpi,jpj,jpk), &
         & znu(jpi,jpj,jpk)        &
         & )

      ALLOCATE( zspgu_tlout (jpi,jpj), zspgv_tlout (jpi,jpj), zspgu_adin (jpi,jpj), zspgv_adin (jpi,jpj))

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      zua_tlin (:,:,:) = 0.0_wp
      zva_tlin (:,:,:) = 0.0_wp
      zua_tlout(:,:,:) = 0.0_wp
      zva_tlout(:,:,:) = 0.0_wp
      zua_adin (:,:,:) = 0.0_wp
      zva_adin (:,:,:) = 0.0_wp
      zua_adout(:,:,:) = 0.0_wp
      zva_adout(:,:,:) = 0.0_wp

      zspgu_adin (:,:) = 0.0_wp
      zspgv_adin (:,:) = 0.0_wp
      zspgu_tlout(:,:) = 0.0_wp
      zspgv_tlout(:,:) = 0.0_wp

      ua_tl(:,:,:) = 0.0_wp
      va_tl(:,:,:) = 0.0_wp
      spgu_tl(:,:) = 0.0_wp
      spgv_tl(:,:) = 0.0_wp
      ua_ad(:,:,:) = 0.0_wp
      va_ad(:,:,:) = 0.0_wp
      spgu_ad(:,:) = 0.0_wp
      spgv_ad(:,:) = 0.0_wp
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
      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      ua_tl(:,:,:) = zua_tlin(:,:,:)
      va_tl(:,:,:) = zva_tlin(:,:,:)

      CALL zdf_bfr_tan( nit000 )

      zua_tlout(:,:,:) = ua_tl(:,:,:)   ;   zva_tlout(:,:,:) = va_tl(:,:,:)
      zspgu_tlout(:,:) = spgu_tl(:,:)   ;   zspgv_tlout(:,:) = spgv_tl(:,:)

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
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zspgu_adin (ji,jj) = zspgu_tlout (ji,jj) &
               &              * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,1) * umask(ji,jj,1)
            zspgv_adin(ji,jj) = zspgv_tlout(ji,jj) &
               &              * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,1) * vmask(ji,jj,1)
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 =   DOT_PRODUCT( zua_tlout  , zua_adin   ) &
         &   + DOT_PRODUCT( zspgu_tlout , zspgu_adin  ) &
         &   + DOT_PRODUCT( zspgv_tlout , zspgv_adin  ) &
         &   + DOT_PRODUCT( zva_tlout  , zva_adin   )


      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      ua_ad(:,:,:) = zua_adin(:,:,:)
      va_ad(:,:,:) = zva_adin(:,:,:)

      spgu_ad(:,:)   = zspgu_adin(:,:)
      spgv_ad(:,:)   = zspgv_adin(:,:)

      CALL zdf_bfr_adj( nitend )

      zua_adout(:,:,:) = ua_ad(:,:,:)
      zva_adout(:,:,:) = va_ad(:,:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2 =   DOT_PRODUCT( zua_tlin  , zua_adout   ) &
         &   + DOT_PRODUCT( zva_tlin  , zva_adout   )

      cl_name = 'zdf_bfr_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE( &
         & zua_tlin,  &
         & zva_tlin,  &
         & zua_tlout, &
         & zva_tlout, &
         & zua_adin,  &
         & zva_adin,  &
         & zua_adout, &
         & zva_adout, &
         & znu        &
         & )
   END SUBROUTINE zdf_bfr_adj_tst

   !!======================================================================
END MODULE zdfbfr_tam
