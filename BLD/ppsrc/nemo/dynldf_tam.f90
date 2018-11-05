MODULE dynldf_tam
   !!======================================================================
   !!                       ***  MODULE  dynldf_tam  ***
   !! Ocean physics:  lateral diffusivity trends
   !!                 Tangent and Adjoint module
   !!=====================================================================
   !! History of the direct module:
   !!          9.0  !  05-11  (G. Madec)  Original code (new step architecture)
   !! History of the TAM module
   !!          9.0  !  08-06  (A. Vidard) Skeleton
   !!               !  08-08  (A. Vidard) TAM of 9.0
   !!   NEMO   3.4  !  12-07  (P.-A. Bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   dyn_ldf     : update the dynamics trend with the lateral diffusion
   !!   dyn_ldf_init_tam : initialization, namelist read, and parameters control
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE ldfdyn_oce
   USE ldfslp
!   USE dynldf_bilapg_tam    ! lateral mixing       (dyn_ldf_bilapg routine)
   USE dynldf_bilap_tam
!   USE dynldf_iso_tam       ! lateral mixing       (dyn_ldf_iso    routine)
   USE dynldf_lap_tam
   USE in_out_manager
!   USE lib_mpp        , ONLY: & ! distribued memory computing library
!   USE lbclnk         , ONLY: & ! ocean lateral boundary conditions (or mpp link)
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE timing
   USE wrk_nemo

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_ldf_tan        ! called by step_tam module
   PUBLIC   dyn_ldf_adj        ! called by step_tam module
   PUBLIC   dyn_ldf_adj_tst    ! called by the tst  module
   PUBLIC   dyn_ldf_init_tam

   INTEGER ::   nldf = -2   ! type of lateral diffusion used defined from ln_dynldf_... namlist logicals)

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
   !!---------------------------------------------------------------------------------

CONTAINS

   SUBROUTINE dyn_ldf_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the lateral ocean dynamics physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_tan')
      !
      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )
         CALL dyn_ldf_lap_tan    ( kt )      ! iso-level laplacian
      CASE ( 1 )
         CALL ctl_stop('dyn_ldf_iso_tan not available yet')
!         CALL dyn_ldf_iso_tan    ( kt )      ! rotated laplacian (except dk[ dk[.] ] part)
      CASE ( 2 )
         CALL dyn_ldf_bilap_tan  ( kt )      ! iso-level bilaplacian
      CASE ( 3 )
         CALL ctl_stop('dyn_ldf_bilapg_tan not available yet')
!         CALL dyn_ldf_bilapg_tan ( kt )      ! s-coord. horizontal bilaplacian
      CASE ( 4 )                                        ! iso-level laplacian + bilaplacian
         CALL dyn_ldf_lap_tan    ( kt )
         CALL dyn_ldf_bilap_tan  ( kt )
      CASE ( 5 )                                        ! rotated laplacian + bilaplacian (s-coord)
         CALL ctl_stop('dyn_ldf_bilapg_tan not available yet')
         !CALL dyn_ldf_iso    ( kt )
         !CALL dyn_ldf_bilapg ( kt )
      !
      CASE ( -2 )                                       ! neither laplacian nor bilaplacian schemes used
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_ldf_tan : no lateral diffusion on momentum setup'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
         ENDIF
      END SELECT
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_tan')
      !
   END SUBROUTINE dyn_ldf_tan

   SUBROUTINE dyn_ldf_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the lateral ocean dynamics physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_adj')
      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )
         CALL dyn_ldf_lap_adj    ( kt )      ! iso-level laplacian
      CASE ( 1 )
         CALL ctl_stop('dyn_ldf_iso_adj not available yet')
!         CALL dyn_ldf_iso_adj    ( kt )      ! rotated laplacian (except dk[ dk[.] ] part)
      CASE ( 2 )
         CALL dyn_ldf_bilap_adj  ( kt )      ! iso-level bilaplacian
      CASE ( 3 )
         CALL ctl_stop('dyn_ldf_bilapg_adj not available yet')
!         CALL dyn_ldf_bilapg_adj ( kt )      ! s-coord. horizontal bilaplacian
      CASE ( 4 )                                        ! iso-level laplacian + bilaplacian
         CALL dyn_ldf_lap_adj    ( kt )
         CALL dyn_ldf_bilap_adj  ( kt )
      CASE ( 5 )                                        ! rotated laplacian + bilaplacian (s-coord)
         CALL ctl_stop('dyn_ldf_bilapg_tan not available yet')
         !CALL dyn_ldf_iso    ( kt )
         !CALL dyn_ldf_bilapg ( kt )
      !
      CASE ( -2 )                                       ! neither laplacian nor bilaplacian schemes used
         IF( kt == nit000 ) THEN
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) 'dyn_ldf_adj : no lateral diffusion on momentum setup'
            IF(lwp) WRITE(numout,*) '~~~~~~~ '
         ENDIF
      !
      END SELECT      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_adj')
   END SUBROUTINE dyn_ldf_adj

   SUBROUTINE dyn_ldf_init_tam
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_ldf_ctl_tam  ***
      !!
      !! ** Purpose of the direct routine:
      !!            initializations of the horizontal ocean dynamics physics
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr         ! temporary integers
      !!----------------------------------------------------------------------

      !                                   ! Namelist nam_dynldf: already read in ldfdyn module

      IF(lwp) THEN                        ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ldf_init_tam : Choice of the lateral diffusive operator on dynamics'
         WRITE(numout,*) '~~~~~~~~~~~~~~'
         WRITE(numout,*) '       Namelist nam_dynldf : set lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '          laplacian operator          ln_dynldf_lap   = ', ln_dynldf_lap
         WRITE(numout,*) '          bilaplacian operator        ln_dynldf_bilap = ', ln_dynldf_bilap
         WRITE(numout,*) '          iso-level                   ln_dynldf_level = ', ln_dynldf_level
         WRITE(numout,*) '          horizontal (geopotential)   ln_dynldf_hor   = ', ln_dynldf_hor
         WRITE(numout,*) '          iso-neutral                 ln_dynldf_iso   = ', ln_dynldf_iso
      ENDIF
      !                                   ! control the consistency
      ioptio = 0
      IF( ln_dynldf_lap   )   ioptio = ioptio + 1
      IF( ln_dynldf_bilap )   ioptio = ioptio + 1
      IF( ioptio < 1 ) CALL ctl_stop( '          use ONE of the 2 lap/bilap operator type on dynamics' )
      ioptio = 0
      IF( ln_dynldf_level )   ioptio = ioptio + 1
      IF( ln_dynldf_hor   )   ioptio = ioptio + 1
      IF( ln_dynldf_iso   )   ioptio = ioptio + 1
      IF( ioptio > 1 ) CALL ctl_stop( '          use only ONE direction (level/hor/iso)' )
      !                                   ! Set nldf, the type of lateral diffusion, from ln_dynldf_... logicals
      ierr = 0
      IF ( ln_dynldf_lap ) THEN      ! laplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_dynldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_dynldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_dynldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_dynldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 1      ! horizontal (   rotation)
            IF ( ln_dynldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ln_dynldf_bilap ) THEN      ! bilaplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_dynldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_dynldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_dynldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_dynldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 3      ! horizontal (   rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ln_dynldf_lap .AND. ln_dynldf_bilap ) THEN  ! mixed laplacian and bilaplacian operators
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_dynldf_level )   nldf = 4      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 4      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_dynldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_dynldf_hor   )   nldf = 4      ! horizontal (no rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! s-coordinate
            IF ( ln_dynldf_level )   nldf = 4      ! iso-level  (no rotation)
            IF ( ln_dynldf_hor   )   nldf = 5      ! horizontal (   rotation)
            IF ( ln_dynldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF


      IF( ierr == 1 )   CALL ctl_stop( 'iso-level in z-coordinate - partial step, not allowed' )
      IF( ierr == 2 )   CALL ctl_stop( 'isoneutral bilaplacian operator does not exist' )
      IF( nldf == 1 .OR. nldf == 3 ) THEN      ! rotation
         IF( .NOT.lk_ldfslp )   CALL ctl_stop( 'the rotation of the diffusive tensor require key_ldfslp' )
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         IF( nldf == -2 )   WRITE(numout,*) '              neither laplacian nor bilaplacian schemes used'
         IF( nldf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nldf ==  0 )   WRITE(numout,*) '              laplacian operator'
         IF( nldf ==  1 )   WRITE(numout,*) '              rotated laplacian operator'
         IF( nldf ==  2 )   WRITE(numout,*) '              bilaplacian operator'
         IF( nldf ==  3 )   WRITE(numout,*) '              rotated bilaplacian'
         IF( nldf ==  4 )   WRITE(numout,*) '              laplacian and bilaplacian operators'
         IF( nldf ==  5 )   WRITE(numout,*) '              rotated laplacian and bilaplacian operators'
      ENDIF
      !
   END SUBROUTINE dyn_ldf_init_tam

   SUBROUTINE dyn_ldf_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_ldf_adj_tst ***
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
      !! ** Action  : Separate tests are applied for the following dx and dy:
      !!
      !!              1) dx = ( SSH ) and dy = ( SSH )
      !!
      !! History :
      !!        ! 08-08 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER :: &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & jt
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zua_tlin,     & ! Tangent input: after u-velocity
         & zva_tlin,     & ! Tangent input: after u-velocity
         & zua_tlout,    & ! Tangent output:after u-velocity
         & zva_tlout,    & ! Tangent output:after v-velocity
         & zua_adin,     & ! adjoint input: after u-velocity
         & zva_adin,     & ! adjoint input: after v-velocity
         & zua_adout,    & ! adjoint output:after v-velocity
         & zva_adout,    & ! adjoint output:after u-velocity
         & zrotb_tlin,   &
         & zhdivb_tlin,  &
         & zrotb_adout,  &
         & zhdivb_adout, &
         & zrotb,        & ! 3D random field for rotb
         & zhdivb,       & ! 3D random field for hdivb
         & zau,          & ! 3D random field for u
         & zav             ! 3D random field for v
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4
      CHARACTER(LEN=14) :: cl_name

      INTEGER :: ildf !: store the nldf flag

      ! Allocate memory

      ALLOCATE( &
         & zua_tlin(jpi,jpj,jpk),     &
         & zva_tlin(jpi,jpj,jpk),     &
         & zua_tlout(jpi,jpj,jpk),    &
         & zva_tlout(jpi,jpj,jpk),    &
         & zua_adin(jpi,jpj,jpk),     &
         & zva_adin(jpi,jpj,jpk),     &
         & zua_adout(jpi,jpj,jpk),    &
         & zva_adout(jpi,jpj,jpk),    &
         & zrotb_tlin(jpi,jpj,jpk),   &
         & zhdivb_tlin(jpi,jpj,jpk),  &
         & zrotb_adout(jpi,jpj,jpk),  &
         & zhdivb_adout(jpi,jpj,jpk), &
         & zrotb(jpi,jpj,jpk),        &
         & zhdivb(jpi,jpj,jpk),       &
         & zau(jpi,jpj,jpk),          &
         & zav(jpi,jpj,jpk)           &
         & )

      ildf = nldf

      DO jt = 1, 2

         IF (jt == 1) nldf=0  ! iso-level laplacian
         IF (jt == 2) nldf=2  ! iso-level bilaplacian

         !==================================================================
         ! 1)      dx = ( ua_tl, va_tl, rotb_tl, hdivb_tl )
         !    and  dy = ( ua_tl, va_tl )
         !==================================================================

         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------
         zua_tlin(:,:,:)     = 0.0_wp
         zva_tlin(:,:,:)     = 0.0_wp
         zrotb_tlin(:,:,:)   = 0.0_wp
         zhdivb_tlin(:,:,:)  = 0.0_wp
         zua_tlout(:,:,:)    = 0.0_wp
         zva_tlout(:,:,:)    = 0.0_wp
         zua_adin(:,:,:)     = 0.0_wp
         zva_adin(:,:,:)     = 0.0_wp
         zrotb_adout(:,:,:)  = 0.0_wp
         zhdivb_adout(:,:,:) = 0.0_wp
         zua_adout(:,:,:)    = 0.0_wp
         zva_adout(:,:,:)    = 0.0_wp
         zrotb(:,:,:)        = 0.0_wp
         zhdivb(:,:,:)       = 0.0_wp
         zau(:,:,:)          = 0.0_wp
         zav(:,:,:)          = 0.0_wp

         ua_tl(:,:,:)    = 0.0_wp
         va_tl(:,:,:)    = 0.0_wp
         ua_ad(:,:,:)    = 0.0_wp
         va_ad(:,:,:)    = 0.0_wp
         rotb_tl(:,:,:)  = 0.0_wp
         hdivb_tl(:,:,:) = 0.0_wp
         rotb_ad(:,:,:)  = 0.0_wp
         hdivb_ad(:,:,:) = 0.0_wp

         !--------------------------------------------------------------------
         ! Initialize the tangent input with random noise: dx
         !--------------------------------------------------------------------

         CALL grid_random(  zau, 'U', 0.0_wp, stdu )
         CALL grid_random(  zav, 'V', 0.0_wp, stdv )
         CALL grid_random(  zrotb, 'F', 0.0_wp, stdr )
         CALL grid_random(  zhdivb, 'T', 0.0_wp, stdh )

         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zua_tlin   (ji,jj,jk) = zau   (ji,jj,jk)
                  zva_tlin   (ji,jj,jk) = zav   (ji,jj,jk)
                  zhdivb_tlin(ji,jj,jk) = zhdivb(ji,jj,jk)
                  zrotb_tlin (ji,jj,jk) = zrotb (ji,jj,jk)
               END DO
            END DO
         END DO
         hdivb_tl(:,:,:) = zhdivb_tlin(:,:,:)
         rotb_tl (:,:,:) = zrotb_tlin (:,:,:)
         ua_tl   (:,:,:) = zua_tlin   (:,:,:)
         va_tl   (:,:,:) = zva_tlin   (:,:,:)

         IF (nldf == 0 )  CALL dyn_ldf_lap_tan(   nit000 )
         IF (nldf == 2 )  CALL dyn_ldf_bilap_tan( nit000 )

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

         zsp1_1 = DOT_PRODUCT( zua_tlout, zua_adin )
         zsp1_2 = DOT_PRODUCT( zva_tlout, zva_adin )
         zsp1   = zsp1_1 + zsp1_2

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------

         ua_ad(:,:,:) = zua_adin(:,:,:)
         va_ad(:,:,:) = zva_adin(:,:,:)

         IF (nldf == 0 )  CALL dyn_ldf_lap_adj(   nit000 )
         IF (nldf == 2 )  CALL dyn_ldf_bilap_adj( nit000 )

         zua_adout   (:,:,:) = ua_ad   (:,:,:)
         zva_adout   (:,:,:) = va_ad   (:,:,:)
         zrotb_adout (:,:,:) = rotb_ad (:,:,:)
         zhdivb_adout(:,:,:) = hdivb_ad(:,:,:)

         !--------------------------------------------------------------------
         ! Compute the scalar product: dx^T dx^*
         !--------------------------------------------------------------------

         zsp2_1 = DOT_PRODUCT( zua_tlin,    zua_adout    )
         zsp2_2 = DOT_PRODUCT( zva_tlin,    zva_adout    )
         zsp2_3 = DOT_PRODUCT( zrotb_tlin,  zrotb_adout  )
         zsp2_4 = DOT_PRODUCT( zhdivb_tlin, zhdivb_adout )
         zsp2   = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4

      ! Compare the scalar products
      ! 14 char:'12345678901234'
         IF (nldf == 0 )  cl_name = 'dynldf_adj lap'
         IF (nldf == 2 )  cl_name = 'dynldf_adj blp'
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      END DO

      nldf = ildf ! restore nldf
      
      DEALLOCATE( &
         & zua_tlin,     &
         & zva_tlin,     &
         & zua_tlout,    &
         & zva_tlout,    &
         & zua_adin,     &
         & zva_adin,     &
         & zua_adout,    &
         & zva_adout,    &
         & zrotb_tlin,   &
         & zhdivb_tlin,  &
         & zrotb_adout,  &
         & zhdivb_adout, &
         & zrotb,        &
         & zhdivb,       &
         & zau,          &
         & zav           &
         & )
   END SUBROUTINE dyn_ldf_adj_tst
   !!======================================================================
END MODULE dynldf_tam
