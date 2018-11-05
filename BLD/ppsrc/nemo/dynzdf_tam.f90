MODULE dynzdf_tam
   !!==============================================================================
   !!                 ***  MODULE  dynzdf_tam  ***
   !! Ocean dynamics :  vertical component of the momentum mixing trend
   !!                   Tangent and adjoint module
   !!==============================================================================
   !! History of the direct module:
   !!         9.0  !  05-11  (G. Madec)  Original code
   !! History of the T&A module:
   !!         9.0  !  08-06  (A. Vidard) Skeleton
   !!         9.0  !  08-08  (A. Vidard) tam of the 05-11 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_zdf      : Update the momentum trend with the vertical diffusion
   !!       zdf_ctl  : initializations of the vertical diffusion scheme
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE zdf_oce
   USE dynzdf_exp_tam
   USE dynzdf_imp_tam
   USE ldfdyn_oce
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE lib_mpp
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_zdf_tan    !  routine called by step_tam.F90
   PUBLIC   dyn_zdf_adj    !  routine called by step_tam.F90
   PUBLIC   dyn_zdf_adj_tst!  routine called by tst.F90
   PUBLIC   dyn_zdf_init_tam

   INTEGER  ::   nzdf = 0              ! type vertical diffusion algorithm used
      !                                ! defined from ln_zdf...  namlist logicals)

   REAL(wp) ::   r2dt                  ! time-step, = 2 rdttra
      !                                ! except at nit000 (=rdttra) if neuler=0
   LOGICAL :: lfirst = .TRUE.

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

   SUBROUTINE dyn_zdf_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the vertical ocean dynamics physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_tan')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000    ) THEN   ;   r2dt =      rdt      ! = rdtra (restarting with Euler time stepping)
      ELSEIF(               kt <= nit000 + 1) THEN   ;   r2dt = 2. * rdt      ! = 2 rdttra (leapfrog)
      ENDIF
      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )   ;   CALL dyn_zdf_exp_tan    ( kt, r2dt )      ! explicit scheme
      CASE ( 1 )   ;   CALL dyn_zdf_imp_tan    ( kt, r2dt )      ! implicit scheme (k-j-i loop)
      !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf_tan')
      !
   END SUBROUTINE dyn_zdf_tan

   SUBROUTINE dyn_zdf_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_zdf_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the vertical ocean dynamics physics.
      !!---------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_zdf_adj')
      !
      !                                          ! set time step
      IF( neuler == 0 .AND. kt == nit000    ) THEN   ;   r2dt =      rdt      ! = rdtra (restarting with Euler time stepping)
      ELSEIF(               kt <= nit000 + 1) THEN   ;   r2dt = 2. * rdt      ! = 2 rdttra (leapfrog)
      ELSEIF( kt == nitend ) THEN ; r2dt = 2. * rdt
      ENDIF
      SELECT CASE ( nzdf )                       ! compute lateral mixing trend and add it to the general trend
      !
      CASE ( 0 )   ;   CALL dyn_zdf_exp_adj    ( kt, r2dt )      ! explicit scheme
      CASE ( 1 )   ;   CALL dyn_zdf_imp_adj    ( kt, r2dt )      ! implicit scheme (k-j-i loop)
      !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_zdf_adj')
      !
   END SUBROUTINE dyn_zdf_adj
   SUBROUTINE dyn_zdf_init_tam
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE zdf_ctl_tam  ***
      !!
      !! ** Purpose :   initializations of the vertical diffusion scheme
      !!
      !! ** Method  :   implicit (euler backward) scheme (default)
      !!                explicit (time-splitting) scheme if ln_zdfexp=T
      !!----------------------------------------------------------------------
      USE zdfgls
      USE zdftke, ONLY : lk_zdftke
      USE zdfkpp, ONLY : lk_zdfkpp
      !!----------------------------------------------------------------------

      ! Choice from ln_zdfexp read in namelist in zdfini
      IF( ln_zdfexp ) THEN   ;   nzdf = 0           ! use explicit scheme
      ELSE                   ;   nzdf = 1           ! use implicit scheme
      ENDIF

      ! Force implicit schemes
      IF( lk_zdfgls .OR. lk_zdftke .OR. lk_zdfkpp )   nzdf = 1   ! TKE or KPP physics
      IF( ln_dynldf_iso                               )   nzdf = 1   ! iso-neutral lateral physics
      IF( ln_dynldf_hor .AND. ln_sco                  )   nzdf = 1   ! horizontal lateral physics in s-coordinate

      IF( lk_esopa )    nzdf = -1                   ! Esopa key: All schemes used

      IF(lwp) THEN                                  ! Print the choice
         WRITE(numout,*)
         WRITE(numout,*) 'dyn:zdf_init_tam : vertical dynamics physics scheme'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         IF( nzdf == -1 )   WRITE(numout,*) '              ESOPA test All scheme used'
         IF( nzdf ==  0 )   WRITE(numout,*) '              Explicit time-splitting scheme'
         IF( nzdf ==  1 )   WRITE(numout,*) '              Implicit (euler backward) scheme'
      ENDIF
      !
      lfirst = .FALSE.
   END SUBROUTINE dyn_zdf_init_tam
   SUBROUTINE dyn_zdf_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_zdf_adj_tst ***
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
         & kumadt          ! Output unit
      INTEGER :: &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk

      !! * Local declarations
       REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zub_tlin,     & ! Tangent input: before u-velocity
         & zvb_tlin,     & ! Tangent input: before v-velocity
         & zua_tlin,     & ! Tangent input: after u-velocity
         & zva_tlin,     & ! Tangent input: after v-velocity
         & zua_tlout,    & ! Tangent output: after u-velocity
         & zva_tlout,    & ! Tangent output: after v-velocity
         & zua_adin,     & ! Adjoint input: after u-velocity
         & zva_adin,     & ! Adjoint input: after v-velocity
         & zub_adout,    & ! Adjoint output: before u-velocity
         & zvb_adout,    & ! Adjoint output: before v-velocity
         & zua_adout,    & ! Adjoint output: after u-velocity
         & zva_adout,    & ! Adjoint output: after v-velocity
         & zau,          & ! 3D random field for ua
         & zav,          & ! 3D random field for va
         & zbu,          & ! 3D random field for ub
         & zbv             ! 3D random field for vb
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
      ! Allocate memory

      ALLOCATE( &
         & zua_tlin(jpi,jpj,jpk),     &
         & zva_tlin(jpi,jpj,jpk),     &
         & zub_tlin(jpi,jpj,jpk),     &
         & zvb_tlin(jpi,jpj,jpk),     &
         & zua_tlout(jpi,jpj,jpk),    &
         & zva_tlout(jpi,jpj,jpk),    &
         & zua_adin(jpi,jpj,jpk),     &
         & zva_adin(jpi,jpj,jpk),     &
         & zua_adout(jpi,jpj,jpk),    &
         & zva_adout(jpi,jpj,jpk),    &
         & zub_adout(jpi,jpj,jpk),    &
         & zvb_adout(jpi,jpj,jpk),    &
         & zau(jpi,jpj,jpk),          &
         & zav(jpi,jpj,jpk),          &
         & zbu(jpi,jpj,jpk),          &
         & zbv(jpi,jpj,jpk)           &
         & )

      ! Initialize the direct trajectory
      avmu(:,:,:) = rn_avm0 * umask(:,:,:)
      avmv(:,:,:) = rn_avm0 * vmask(:,:,:)

      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      zua_tlin(:,:,:)  = 0.0_wp
      zva_tlin(:,:,:)  = 0.0_wp
      zub_tlin(:,:,:)  = 0.0_wp
      zvb_tlin(:,:,:)  = 0.0_wp
      zua_tlout(:,:,:) = 0.0_wp
      zva_tlout(:,:,:) = 0.0_wp
      zua_adin(:,:,:)  = 0.0_wp
      zva_adin(:,:,:)  = 0.0_wp
      zua_adout(:,:,:) = 0.0_wp
      zva_adout(:,:,:) = 0.0_wp
      zub_adout(:,:,:) = 0.0_wp
      zvb_adout(:,:,:) = 0.0_wp
      zau(:,:,:)       = 0.0_wp
      zav(:,:,:)       = 0.0_wp
      zbu(:,:,:)       = 0.0_wp
      zbv(:,:,:)       = 0.0_wp

      ub_tl(:,:,:)     = 0.0_wp
      vb_tl(:,:,:)     = 0.0_wp
      ua_tl(:,:,:)     = 0.0_wp
      va_tl(:,:,:)     = 0.0_wp
      ub_ad(:,:,:)     = 0.0_wp
      vb_ad(:,:,:)     = 0.0_wp
      ua_ad(:,:,:)     = 0.0_wp
      va_ad(:,:,:)     = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  zbu, 'U', 0.0_wp, stdu )
      CALL grid_random(  zbv, 'V', 0.0_wp, stdv )
      CALL grid_random(  zau, 'U', 0.0_wp, stdu )
      CALL grid_random(  zav, 'V', 0.0_wp, stdv )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zub_tlin(ji,jj,jk) = zbu(ji,jj,jk)
               zvb_tlin(ji,jj,jk) = zbv(ji,jj,jk)
               zua_tlin(ji,jj,jk) = zau(ji,jj,jk)
               zva_tlin(ji,jj,jk) = zav(ji,jj,jk)
            END DO
         END DO
      END DO

      ub_tl(:,:,:) = zub_tlin(:,:,:)
      vb_tl(:,:,:) = zvb_tlin(:,:,:)
      ua_tl(:,:,:) = zua_tlin(:,:,:)
      va_tl(:,:,:) = zva_tlin(:,:,:)

      CALL dyn_zdf_tan( nit000 )

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

      CALL dyn_zdf_adj ( nit000 )
      zub_adout(:,:,:)   = ub_ad(:,:,:)
      zvb_adout(:,:,:)   = vb_ad(:,:,:)
      zua_adout(:,:,:)   = ua_ad(:,:,:)
      zva_adout(:,:,:)   = va_ad(:,:,:)

      zsp2_1 = DOT_PRODUCT( zub_tlin, zub_adout )
      zsp2_2 = DOT_PRODUCT( zvb_tlin, zvb_adout )
      zsp2_3 = DOT_PRODUCT( zua_tlin, zua_adout )
      zsp2_4 = DOT_PRODUCT( zva_tlin, zva_adout )
      zsp2   = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4

      ! Compare the scalar products

      ! 14 char:'12345678901234'
      cl_name = 'dyn_zdf_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE( &
         & zua_tlin,     &
         & zva_tlin,     &
         & zub_tlin,     &
         & zvb_tlin,     &
         & zua_tlout,    &
         & zva_tlout,    &
         & zua_adin,     &
         & zva_adin,     &
         & zua_adout,    &
         & zva_adout,    &
         & zub_adout,    &
         & zvb_adout,    &
         & zau,          &
         & zav,          &
         & zbu,          &
         & zbv           &
         & )

  END SUBROUTINE dyn_zdf_adj_tst
   !!==============================================================================
END MODULE dynzdf_tam
