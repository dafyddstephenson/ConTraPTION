MODULE sshwzv_tam
   !!==============================================================================
   !!                       ***  MODULE  sshwzv  ***
   !! Ocean dynamics : sea surface height and vertical velocity
   !!==============================================================================
   !! History of the direct module:
   !!            3.1  !  2009-02  (G. Madec, M. Leclair)  Original code
   !! History of the TAM module:
   !!            3.2  !  2010-04  (F. Vigilant) Original code
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ssh_wzv        : after ssh & now vertical velocity
   !!   ssh_nxt        : filter ans swap the ssh arrays
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce
   USE in_out_manager
   USE dom_oce
   USE prtctl
   USE phycst
   USE lbclnk
   USE lbclnk_tam
   USE divcur_tam
   USE cla_tam
   USE oce_tam
   USE sbc_oce_tam
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE tstool_tam
   USE lib_mpp
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ssh_wzv_tan       ! called by step.F90
   PUBLIC   ssh_nxt_tan       ! called by step.F90
   PUBLIC   ssh_wzv_adj       ! called by step.F90
   PUBLIC   ssh_nxt_adj       ! called by step.F90
   PUBLIC   ssh_wzv_adj_tst   ! called by tamtst.F90
   PUBLIC   ssh_nxt_adj_tst   ! called by tamtst.F90

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

CONTAINS

   SUBROUTINE ssh_wzv_tan( kt , kdum )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE ssh_wzv_tan  ***
      !!
      !! ** Purpose of direct routine :
      !!              compute the after ssh (ssha), the now vertical velocity
      !!              and update the now vertical coordinate (lk_vvl=T).
      !!
      !! ** Method  : -
      !!
      !!              - Using the incompressibility hypothesis, the vertical
      !!      velocity is computed by integrating the horizontal divergence
      !!      from the bottom to the surface minus the scale factor evolution.
      !!        The boundary conditions are w=0 at the bottom (no flux) and.
      !!
      !! ** action  :   ssha    : after sea surface height
      !!                wn      : now vertical velocity
      !! if lk_vvl=T:   sshu_a, sshv_a, sshf_a  : after sea surface height
      !!                          at u-, v-, f-point s
      !!                hu, hv, hur, hvr : ocean depth and its inverse at u-,v-points
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt   ! time step
      !!
      INTEGER  ::   jk              ! dummy loop indices
      REAL(wp) ::   z2dt, z1_rau0     ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::   z2d, zhdivtl     ! 2D workspace
      INTEGER, OPTIONAL            ::   kdum        ! dummy argument to compute only vertical velocity
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ssh_wzv_tan')
      !
      CALL wrk_alloc( jpi, jpj, z2d, zhdivtl )
      z2d = 0._wp
      zhdivtl = 0._wp
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_wzv_tan : after sea surface height and now vertical velocity '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
         !
         wn_tl(:,:,jpk) = 0.0_wp                   ! bottom boundary condition: w=0 (set once for all)
         !
         IF( lk_vvl ) THEN                    ! before and now Sea SSH at u-, v-, f-points (vvl case only)
            IF (lwp) WRITE(numout,*) 'lk_vvl not available yet'
            CALL abort
         ENDIF
         !
      ENDIF
      !                                           !------------------------------!
      IF( lk_vvl ) THEN                           !  Update Now Vertical coord.  !   (only in vvl case)
                                                  !------------------------------!
         IF (lwp) WRITE(numout,*) 'lk_vvl not available yet'
         CALL abort
         !
      ENDIF

      CALL div_cur_tan( kt )            ! Horizontal divergence & Relative vorticity

      ! set time step size (Euler/Leapfrog)
      z2dt = 2. * rdt
      IF( neuler == 0 .AND. kt == nit000 )   z2dt =rdt

      z1_rau0 = 0.5_wp / rau0

      IF ( .NOT. PRESENT(kdum) ) THEN             ! jump ssh computing
         !                                        !------------------------------!
         !                                        !   After Sea Surface Height   !
         !                                        !------------------------------!
         zhdivtl(:,:) = 0.0_wp
         DO jk = 1, jpkm1                                 ! Horizontal divergence of barotropic transports
            zhdivtl(:,:) = zhdivtl(:,:) + e3t(:,:,jk) * hdivn_tl(:,:,jk)
         END DO

         !                                                ! Sea surface elevation time stepping
         ssha_tl(:,:) = (  sshb_tl(:,:) - z2dt * ( z1_rau0 * ( emp_b_tl(:,:) + emp_tl(:,:) ) + zhdivtl(:,:) )  ) * tmask(:,:,1)

         !                                                ! Sea Surface Height at u-,v- and f-points (vvl case only)
         IF( lk_vvl ) THEN                                ! (required only in key_vvl case)
            IF (lwp) WRITE(numout,*) 'lk_vvl not available yet'
            CALL abort
         ENDIF

      ENDIF
      !                                           !------------------------------!
      !                                           !     Now Vertical Velocity    !
      !                                           !------------------------------!
      !                                                ! integrate from the bottom the hor. divergence
      DO jk = jpkm1, 1, -1
         wn_tl(:,:,jk) = wn_tl(:,:,jk+1) - e3t(:,:,jk) * hdivn_tl(:,:,jk)
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, z2d, zhdivtl )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ssh_wzv_tan')
      !
   END SUBROUTINE ssh_wzv_tan

   SUBROUTINE ssh_nxt_tan( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ssh_nxt_tan  ***
      !!
      !! ** Purpose of the direct :
      !!              achieve the sea surface  height time stepping by
      !!              applying Asselin time filter and swapping the arrays
      !!              ssha  already computed in ssh_wzv
      !!
      !! ** Method  : - apply Asselin time fiter to now ssh and swap :
      !!             sshn = ssha + atfp * ( sshb -2 sshn + ssha )
      !!             sshn = ssha
      !!
      !! ** action  : - sshb, sshn   : before & now sea surface height
      !!                               ready for the next time step
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER  ::   ji, jj               ! dummy loop indices
      REAL(wp) :: zec
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ssh_nxt_tan')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_nxt_tan : next sea surface height (Asselin time filter + swap)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      ! Time filter and swap of the ssh
      ! -------------------------------

      IF( lk_vvl ) THEN      ! Variable volume levels :   ssh at t-, u-, v, f-points
         !                   ! ---------------------- !
         IF (lwp) WRITE(numout,*) 'lk_vvl not available yet'
         CALL abort
         !
      ELSE                   ! fixed levels :   ssh at t-point only
         !                   ! ------------ !
         IF( neuler == 0 .AND. kt == nit000 ) THEN      ! Euler time-stepping at first time-step : no filter
            sshn_tl(:,:) = ssha_tl(:,:)                            ! now <-- after  (before already = now)
            !
         ELSE                                           ! Leap-Frog time-stepping: Asselin filter + swap
            DO jj = 1, jpj
               DO ji = 1, jpi                                ! before <-- now filtered
                  sshb_tl(ji,jj) = sshn_tl(ji,jj) + atfp * ( sshb_tl(ji,jj) - 2 * sshn_tl(ji,jj) + ssha_tl(ji,jj) )
                  sshn_tl(ji,jj) = ssha_tl(ji,jj)                  ! now <-- after
               END DO
            END DO
         ENDIF
      ENDIF
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('ssh_nxt_tan')
      !
   END SUBROUTINE ssh_nxt_tan

   SUBROUTINE ssh_wzv_adj( kt , kdum )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE ssh_wzv_adj  ***
      !!
      !! ** Purpose of direct routine :
      !!              compute the after ssh (ssha), the now vertical velocity
      !!              and update the now vertical coordinate (lk_vvl=T).
      !!
      !! ** Method  : -
      !!
      !!              - Using the incompressibility hypothesis, the vertical
      !!      velocity is computed by integrating the horizontal divergence
      !!      from the bottom to the surface minus the scale factor evolution.
      !!        The boundary conditions are w=0 at the bottom (no flux) and.
      !!
      !! ** action  :   ssha    : after sea surface height
      !!                wn      : now vertical velocity
      !! if lk_vvl=T:   sshu_a, sshv_a, sshf_a  : after sea surface height
      !!                          at u-, v-, f-point s
      !!                hu, hv, hur, hvr : ocean depth and its inverse at u-,v-points
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt   ! time step
      !!
      INTEGER  ::   jk              ! dummy loop indices
      REAL(wp) ::   z2dt, z1_rau0     ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:) ::  z2d, zhdivad     ! 2D workspace
      INTEGER, OPTIONAL            ::   kdum        ! dummy argument to compute only vertical velocity
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ssh_wzv_adj')
      !
      CALL wrk_alloc( jpi, jpj, z2d, zhdivad )
      !
      ! adjoint variable initialization
      zhdivad = 0._wp
      z2d = 0._wp

      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_wzv_adj : after sea surface height and now vertical velocity '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      !
      !                                           !------------------------------!
      IF( lk_vvl ) THEN                           !  Update Now Vertical coord.  !   (only in vvl case)
                                                  !------------------------------!
         IF (lwp) WRITE(numout,*) 'lk_vvl not available yet'
         CALL abort
         !
      ENDIF
      !                                           !------------------------------!
      !                                           !     Now Vertical Velocity    !
      !                                           !------------------------------!
      !                                                ! integrate from the bottom the hor. divergence
      DO jk = 1, jpkm1
         hdivn_ad(:,:,jk  ) = hdivn_ad(:,:,jk  ) - wn_ad(:,:,jk) * e3t(:,:,jk)
         wn_ad(   :,:,jk+1) = wn_ad(   :,:,jk+1) + wn_ad(:,:,jk)
         wn_ad(   :,:,jk  ) = 0.0_wp
      END DO
      !
      ! set time step size (Euler/Leapfrog)
      z2dt = 2. * rdt
      IF( neuler == 0 .AND. kt == nit000 )   z2dt =rdt

      z1_rau0 = 0.5_wp / rau0

      IF ( .NOT. PRESENT(kdum) ) THEN             ! jump ssh computing
         !                                        !------------------------------!
         !                                        !   After Sea Surface Height   !
         !                                        !------------------------------!
         !                                                ! Sea Surface Height at u-,v- and f-points (vvl case only)
         IF( lk_vvl ) THEN                                ! (required only in key_vvl case)
            IF (lwp) WRITE(numout,*) 'lk_vvl not available yet'
            CALL abort
         ENDIF
         !                                                ! Sea surface elevation time stepping
         ssha_ad( :,:) = ssha_ad( :,:) * tmask(:,:,1)
         sshb_ad( :,:) = sshb_ad( :,:) + ssha_ad(:,:)
         ssha_ad( :,:) = ssha_ad( :,:) * z2dt
         zhdivad( :,:) = zhdivad( :,:) - ssha_ad(:,:)
         ssha_ad( :,:) = ssha_ad( :,:) * z1_rau0 
         emp_ad(  :,:) = emp_ad(  :,:) - ssha_ad(:,:)
         emp_b_ad(:,:) = emp_b_ad(:,:) - ssha_ad(:,:)
         ssha_ad(:,:) = 0.0_wp

         DO jk = 1, jpkm1                                 ! Horizontal divergence of barotropic transports
            hdivn_ad(:,:,jk) = hdivn_ad(:,:,jk) + e3t(:,:,jk) * zhdivad(:,:)
         END DO
         zhdivad(:,:) = 0._wp

      ENDIF

      CALL div_cur_adj( kt )            ! Horizontal divergence & Relative vorticity
      !
      !
      !
      IF( kt == nit000 ) wn_ad(:,:,jpk) = 0._wp
      !
      !
      !
      CALL wrk_dealloc( jpi, jpj, z2d, zhdivad )
      !
      IF( nn_timing == 1 )  CALL timing_stop('ssh_wzv_adj')
      !
   END SUBROUTINE ssh_wzv_adj

   SUBROUTINE ssh_nxt_adj( kt )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE ssh_nxt_adj  ***
      !!
      !! ** Purpose of the direct :
      !!              achieve the sea surface  height time stepping by
      !!              applying Asselin time filter and swapping the arrays
      !!              ssha  already computed in ssh_wzv
      !!
      !! ** Method  : - apply Asselin time fiter to now ssh and swap :
      !!             sshn = ssha + atfp * ( sshb -2 sshn + ssha )
      !!             sshn = ssha
      !!
      !! ** action  : - sshb, sshn   : before & now sea surface height
      !!                               ready for the next time step
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!
      INTEGER  ::   ji, jj               ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('ssh_nxt_adj')
      !

      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'ssh_nxt_adj : next sea surface height (Asselin time filter + swap)'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF
      ! Time filter and swap of the ssh
      ! -------------------------------
      IF( lk_vvl ) THEN      ! Variable volume levels :   ssh at t-, u-, v, f-points
         !                   ! ---------------------- !
         IF (lwp) WRITE(numout,*) 'lk_vvl not available yet'
         CALL abort
         !
      ELSE                   ! fixed levels :   ssh at t-point only

         !                   ! ------------ !
         IF( neuler == 0 .AND. kt == nit000 ) THEN      ! Euler time-stepping at first time-step : no filter
            ssha_ad(:,:) = ssha_ad(:,:) + sshn_ad(:,:)
            sshn_ad(:,:) = 0.0_wp
            !
         ELSE                                           ! Leap-Frog time-stepping: Asselin filter + swap
            DO jj = jpj, 1, -1
               DO ji = jpi, 1, -1                                ! before <-- now filtered
                  ssha_ad(ji,jj) = ssha_ad(ji,jj) + sshn_ad(ji,jj)
                  sshn_ad(ji,jj) = 0._wp
                  sshn_ad(ji,jj) = sshn_ad(ji,jj) + (1.0_wp - 2.0 * atfp) * sshb_ad(ji,jj)
                  ssha_ad(ji,jj) = ssha_ad(ji,jj) + atfp * sshb_ad(ji,jj)
                  sshb_ad(ji,jj) = atfp * sshb_ad(ji,jj)
               END DO
            END DO
         ENDIF
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('ssh_nxt_adj')
      !
   END SUBROUTINE ssh_nxt_adj

   SUBROUTINE ssh_wzv_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!          ***  ROUTINE ssh_wzv_adj_tst : TEST OF wzv_adj  ***
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
      !!            dx = ( un_tl, vn_tl, hdivn_tl, rotn_tl, emp_tl, sshb_tl ) and
      !!            dy = ( hdivn_tl, hdivb_tl, rotn_tl, rotb_tl, wn_tl, ssha_tl )
      !!
      !! History :
      !!        ! 2010-04 (F. Vigilant)
      !!-----------------------------------------------------------------------

      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zun_tlin,     & ! Tangent input: now u-velocity
         & zvn_tlin,     & ! Tangent input: now v-velocity
         & zhdivn_tlin,  & ! Tangent input: now horizontal divergence
         & zrotn_tlin,   & ! Tangent input: now horizontal divergence
         & zhdivn_tlout, & ! Tangent output: now horizontal divergence
         & zrotn_tlout,  & ! Tangent output: now horizontal divergence
         & zrotb_tlout,  & ! Tangent output: now horizontal divergence
         & zhdivb_tlout, & ! Tangent output: now horizontal divergence
         & zwn_tlout,    & ! Tangent output: now w-velocity
         & zwn_adin,     & ! Adjoint input: now w-velocity
         & zhdivn_adout, & ! Adjoint output: now horizontal divergence
         & zrotn_adin,   & ! Adjoint input: now horizontal divergence
         & zrotn_adout,  & ! Adjoint output: now horizontal divergence
         & zrotb_adin,   & ! Adjoint input: now horizontal divergence
         & zhdivn_adin,  & ! Adjoint input: now horizontal divergence
         & zhdivb_adin,  & ! Adjoint output: now horizontal divergence
         & zun_adout,    & ! Adjoint output: now horizontal divergence
         & zvn_adout,    & ! Adjoint output: now horizontal divergence
         & znu,          & ! 3D random field for u
         & znv             ! 3D random field for v
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zsshb_tlin,   & ! Tangent input: before SSH
         & zssha_tlout,  & ! Tangent input: before SSH
         & zsshb_adout,  & ! Adjoint output: before SSH
         & zssha_adin,   & ! Adjoint output: before SSH
         & zemp_tlin,    & ! Tangent input: EmP
         & zemp_adout,   & ! Adjoint output: EmP
         & znssh,        & ! 2D random field for SSH
         & znemp           ! 2D random field for EmP

      INTEGER :: &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d           ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
                              ! random field standard deviation for:
         & zstdssh,         & !   SSH
         & zstdemp,         & !   EMP
         & zsp1,            & ! scalar product involving the tangent routine
         & zsp2,            & ! scalar product involving the adjoint routine
         & zsp2_1,          & !   scalar product components
         & zsp2_2,          &
         & zsp2_3,          &
         & zsp2_4,          &
         & zsp2_5,          &
         & zsp2_6,          &
         & z2dt,            & ! temporary scalars
         & zraur
      CHARACTER (LEN=14) :: &
         & cl_name

      ! Allocate memory

      ALLOCATE( &
         & zhdivn_tlin(jpi,jpj,jpk),  &
         & zhdivb_tlout(jpi,jpj,jpk), &
         & zhdivn_tlout(jpi,jpj,jpk), &
         & zrotn_tlin(jpi,jpj,jpk),   &
         & zrotn_tlout(jpi,jpj,jpk),  &
         & zrotb_tlout(jpi,jpj,jpk),  &
         & zwn_tlout(jpi,jpj,jpk),    &
         & zwn_adin(jpi,jpj,jpk),     &
         & zhdivn_adout(jpi,jpj,jpk), &
         & zhdivb_adin(jpi,jpj,jpk),  &
         & zrotn_adin(jpi,jpj,jpk),   &
         & zrotn_adout(jpi,jpj,jpk),  &
         & zrotb_adin(jpi,jpj,jpk),   &
         & zhdivn_adin(jpi,jpj,jpk),  &
         & zun_tlin(jpi,jpj,jpk),     &
         & zvn_tlin(jpi,jpj,jpk),     &
         & zun_adout(jpi,jpj,jpk),    &
         & zvn_adout(jpi,jpj,jpk),    &
         & znu(jpi,jpj,jpk),          &
         & znv(jpi,jpj,jpk)           &
         & )
      ALLOCATE( &
         & zsshb_tlin(jpi,jpj),       &
         & zsshb_adout(jpi,jpj),      &
         & zssha_tlout(jpi,jpj),      &
         & zssha_adin(jpi,jpj),       &
         & zemp_tlin(jpi,jpj),        &
         & zemp_adout(jpi,jpj),       &
         & znssh(jpi,jpj),            &
         & znemp(jpi,jpj)             &
         & )


      ! Initialize constants

      z2dt  = 2.0_wp * rdt       ! time step: leap-frog
      zraur = 1.0_wp / rau0      ! inverse density of pure water (m3/kg)

      zhdivn_tlin(:,:,:) = 0.0_wp
      zrotn_tlin(:,:,:) = 0.0_wp
      zemp_tlin(:,:) = 0.0_wp
      zsshb_tlin(:,:) = 0.0_wp
      zun_tlin (:,:,:) = 0.0_wp
      zvn_tlin (:,:,:) = 0.0_wp

      zhdivn_tlout(:,:,:) = 0.0_wp
      zhdivb_tlout(:,:,:) = 0.0_wp
      zrotn_tlout(:,:,:)  = 0.0_wp
      zrotb_tlout(:,:,:)  = 0.0_wp
      zwn_tlout(:,:,:) = 0.0_wp
      zssha_tlout(:,:) = 0.0_wp

      zhdivn_adin(:,:,:) = 0.0_wp
      zhdivb_adin(:,:,:) = 0.0_wp
      zrotn_adin(:,:,:)  = 0.0_wp
      zrotb_adin(:,:,:)  = 0.0_wp
      zwn_adin(:,:,:) = 0.0_wp
      zssha_adin(:,:) = 0.0_wp

      zhdivn_adout(:,:,:) = 0.0_wp
      zrotn_adout(:,:,:) = 0.0_wp
      zemp_adout(:,:) = 0.0_wp
      zsshb_adout(:,:) = 0.0_wp
      zun_adout (:,:,:) = 0.0_wp
      zvn_adout (:,:,:) = 0.0_wp

      un_tl   (:,:,:) = 0.0_wp
      vn_tl   (:,:,:) = 0.0_wp
      hdivn_tl(:,:,:) = 0.0_wp
      hdivb_tl(:,:,:) = 0.0_wp
      rotn_tl (:,:,:) = 0.0_wp
      rotb_tl (:,:,:) = 0.0_wp
      wn_tl(:,:,:) = 0.0_wp
      ssha_tl(:,:) = 0.0_wp
      sshb_tl(:,:) = 0.0_wp
      emp_tl(:,:) = 0.0_wp

      un_ad   (:,:,:) = 0.0_wp
      vn_ad   (:,:,:) = 0.0_wp
      hdivn_ad(:,:,:) = 0.0_wp
      hdivb_ad(:,:,:) = 0.0_wp
      rotn_ad (:,:,:) = 0.0_wp
      rotb_ad (:,:,:) = 0.0_wp
      wn_ad(:,:,:) = 0.0_wp
      sshb_ad(:,:) = 0.0_wp
      ssha_ad(:,:) = 0.0_wp
      emp_ad(:,:) = 0.0_wp

      !=============================================================
      ! 1) dx = ( un_tl, vn_tl, emp_tl, sshb_tl ) and dy = ( wn_tl )
      !                   - or -
      ! 2) dx = ( hdivn_tl ) and dy = ( wn_tl )
      !=============================================================

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  znu, 'U', 0.0_wp, stdu )
      CALL grid_random(  znv, 'V', 0.0_wp, stdv )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zun_tlin(ji,jj,jk) = znu(ji,jj,jk)
               zvn_tlin(ji,jj,jk) = znv(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  znssh, 'T', 0.0_wp, stdssh )
      CALL grid_random(  znemp, 'T', 0.0_wp, stdssh )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsshb_tlin(ji,jj) = znssh(ji,jj)
            zemp_tlin (ji,jj) = znemp(ji,jj) / ( z2dt * zraur )
         END DO
      END DO

      un_tl(:,:,:) = zun_tlin(:,:,:)
      vn_tl(:,:,:) = zvn_tlin(:,:,:)
      CALL div_cur_tan( nit000 )    ! Generate noise hdiv/rot fields

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zhdivn_tlin(ji,jj,jk) = 0.5_wp * hdivn_tl(ji,jj,jk)
               zrotn_tlin (ji,jj,jk) = 0.5_wp * rotn_tl (ji,jj,jk)
            END DO
         END DO
      END DO

      ! re-initialization to zero
      un_tl   (:,:,:) = 0.0_wp
      vn_tl   (:,:,:) = 0.0_wp
      hdivb_tl(:,:,:) = 0.0_wp
      hdivn_tl(:,:,:) = 0.0_wp
      rotb_tl (:,:,:) = 0.0_wp
      rotn_tl (:,:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      hdivn_tl(:,:,:) = zhdivn_tlin(:,:,:)
      rotn_tl(:,:,:) = zrotn_tlin(:,:,:)
      sshb_tl(:,:) = zsshb_tlin(:,:)
      emp_tl (:,:) = zemp_tlin (:,:)
      un_tl(:,:,:) = zun_tlin(:,:,:)
      vn_tl(:,:,:) = zvn_tlin(:,:,:)

      CALL ssh_wzv_tan( nit000+1 )

      zwn_tlout(:,:,:) = wn_tl(:,:,:)
      zssha_tlout(:,: ) = ssha_tl(:,:)
      zhdivb_tlout(:,:,:) = hdivb_tl(:,:,:)
      zhdivn_tlout(:,:,:) = hdivn_tl(:,:,:)
      zrotb_tlout(:,:,:) = rotb_tl(:,:,:)
      zrotn_tlout(:,:,:) = rotn_tl(:,:,:)
      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zwn_adin(ji,jj,jk) = zwn_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zssha_adin(ji,jj) = zssha_tlout(ji,jj) &
               &                   * e1t(ji,jj) * e2t(ji,jj) * wesp_ssh &
               &                   * tmask(ji,jj,1)
         END DO
      END DO
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zhdivb_adin(ji,jj,jk) = zhdivb_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk)
              zhdivn_adin(ji,jj,jk) = zhdivn_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zrotb_adin(ji,jj,jk) = zrotb_tlout(ji,jj,jk) &
                 &               * e1f(ji,jj) * e2f(ji,jj) * e3f(ji,jj,jk)
              zrotn_adin(ji,jj,jk) = zrotn_tlout(ji,jj,jk) &
                 &               * e1f(ji,jj) * e2f(ji,jj) * e3f(ji,jj,jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------


      zsp1 = DOT_PRODUCT( zwn_tlout, zwn_adin ) + DOT_PRODUCT( zssha_tlout, zssha_adin ) &
           & +  DOT_PRODUCT( zhdivb_tlout, zhdivb_adin ) + DOT_PRODUCT( zhdivn_tlout, zhdivn_adin ) &
           & +  DOT_PRODUCT( zrotb_tlout, zrotb_adin ) + DOT_PRODUCT( zrotn_tlout, zrotn_adin )
      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      wn_ad(:,:,:) = zwn_adin(:,:,:)
      ssha_ad(:,:) = zssha_adin(:,:)
      hdivb_ad(:,:,:) = zhdivb_adin(:,:,:)
      hdivn_ad(:,:,:) = zhdivn_adin(:,:,:)
      rotb_ad(:,:,:) = zrotb_adin(:,:,:)
      rotn_ad(:,:,:) = zrotn_adin(:,:,:)

      CALL ssh_wzv_adj( nit000+1 )

      zrotn_adout(:,:,:) = rotn_ad(:,:,:)
      zhdivn_adout(:,:,:) = hdivn_ad(:,:,:)
      zsshb_adout(:,:) = sshb_ad(:,:)
      zemp_adout (:,:) = emp_ad (:,:)
      zun_adout(:,:,:) = un_ad(:,:,:)
      zvn_adout(:,:,:) = vn_ad(:,:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1 = DOT_PRODUCT( zun_tlin,    zun_adout    )
      zsp2_2 = DOT_PRODUCT( zvn_tlin,    zvn_adout    )
      zsp2_3 = DOT_PRODUCT( zhdivn_tlin, zhdivn_adout )
      zsp2_4 = DOT_PRODUCT( zemp_tlin,   zemp_adout   )
      zsp2_5 = DOT_PRODUCT( zsshb_tlin,  zsshb_adout  )
      zsp2_6 = DOT_PRODUCT( zrotn_tlin,  zrotn_adout  )

      zsp2 = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5 + zsp2_6

      ! Compare the scalar products
      ! 14 char:'12345678901234'
      cl_name = 'sshwzv_adj    '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

   END SUBROUTINE ssh_wzv_adj_tst

   SUBROUTINE ssh_nxt_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!          ***  ROUTINE ssh_nxt_adj_tst : TEST OF nxt_adj  ***
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
      !!            dx = ( sshb_tl, sshn_tl, ssha_tl ) and
      !!            dy = ( ssb_tl, sshn_tl )
      !!
      !! History :
      !!        ! 2010-05 (F. Vigilant)
      !!-----------------------------------------------------------------------

      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zsshb_tlin,   & ! Tangent input: before SSH
         & zsshn_tlin,   & ! Tangent input: before SSH
         & zssha_tlin,   & ! Tangent input: before SSH
         & zsshb_tlout,  & ! Tangent output: before SSH
         & zsshn_tlout,  & ! Tangent output: before SSH
         & zsshb_adin,   & ! Adjoint input: before SSH
         & zsshn_adin,   & ! Adjoint input: before SSH
         & zsshb_adout,  & ! Adjoint output: before SSH
         & zsshn_adout,  & ! Adjoint output: before SSH
         & zssha_adout,  & ! Adjoint output: before SSH
         & znssh           ! 2D random field for EmP

      INTEGER :: &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d           ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
                              ! random field standard deviation for:
         & zstdssh,         & !   SSH
         & zsp1,            & ! scalar product involving the tangent routine
         & zsp2,            & ! scalar product involving the adjoint routine
         & zsp1_1,          & !   scalar product components
         & zsp1_2,          &
         & zsp2_1,          & !   scalar product components
         & zsp2_2,          &
         & zsp2_3,          &
         & zsp2_4
      CHARACTER (LEN=14) :: &
         & cl_name

      ! Allocate memory

      ALLOCATE( &
         & zsshb_tlin(jpi,jpj),       &
         & zsshn_tlin(jpi,jpj),       &
         & zssha_tlin(jpi,jpj),       &
         & zsshb_tlout(jpi,jpj),      &
         & zsshn_tlout(jpi,jpj),      &
         & zsshb_adin(jpi,jpj),       &
         & zsshn_adin(jpi,jpj),       &
         & zsshb_adout(jpi,jpj),      &
         & zsshn_adout(jpi,jpj),      &
         & zssha_adout(jpi,jpj),      &
         & znssh(jpi,jpj)             &
         & )


      ! Initialize constants

      zsshb_tlin(:,:)  = 0.0_wp
      zsshn_tlin(:,:)  = 0.0_wp
      zssha_tlin(:,:)  = 0.0_wp

      zsshb_tlout(:,:) = 0.0_wp
      zsshn_tlout(:,:) = 0.0_wp

      zsshb_adout(:,:) = 0.0_wp
      zsshn_adout(:,:) = 0.0_wp
      zssha_adout(:,:) = 0.0_wp

      zsshb_adin(:,:)  = 0.0_wp
      zsshn_adin(:,:)  = 0.0_wp

      sshb_tl(:,:)     = 0.0_wp
      sshn_tl(:,:)     = 0.0_wp
      ssha_tl(:,:)     = 0.0_wp

      sshb_ad(:,:)     = 0.0_wp
      sshn_ad(:,:)     = 0.0_wp
      ssha_ad(:,:)     = 0.0_wp

      !=============================================================
      ! dx = ( sshb_tl, sshn_tl, ssha_tl ) and dy = ( ssb_tl, sshn_tl )
      !=============================================================

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random( znssh, 'T', 0.0_wp, stdssh )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsshb_tlin(ji,jj) = znssh(ji,jj)
         END DO
      END DO

      CALL grid_random( znssh, 'T', 0.0_wp, stdssh )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsshn_tlin(ji,jj) = znssh(ji,jj)
         END DO
      END DO

      CALL grid_random( znssh, 'T', 0.0_wp, stdssh )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zssha_tlin(ji,jj) = znssh(ji,jj)
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      sshb_tl(:,:) = zsshb_tlin(:,:)
      sshn_tl(:,:) = zsshn_tlin(:,:)
      ssha_tl(:,:) = zssha_tlin(:,:)

      CALL ssh_nxt_tan( nit000+1 )

      zsshb_tlout(:,: ) = sshb_tl(:,:)
      zsshn_tlout(:,: ) = sshn_tl(:,:)
      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsshb_adin(ji,jj) = zsshb_tlout(ji,jj) &
               &                   * e1t(ji,jj) * e2t(ji,jj) * wesp_ssh &
               &                   * tmask(ji,jj,1)
            zsshn_adin(ji,jj) = zsshn_tlout(ji,jj) &
               &                   * e1t(ji,jj) * e2t(ji,jj) * wesp_ssh &
               &                   * tmask(ji,jj,1)
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------
      zsp1_1 = DOT_PRODUCT( zsshb_tlout, zsshb_adin )
      zsp1_2 = DOT_PRODUCT( zsshn_tlout, zsshn_adin )

      zsp1 = zsp1_1 + zsp1_2
      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      sshb_ad(:,:) = zsshb_adin(:,:)
      sshn_ad(:,:) = zsshn_adin(:,:)

      CALL ssh_nxt_adj( nit000+1 )

      zsshb_adout(:,:) = sshb_ad(:,:)
      zsshn_adout(:,:) = sshn_ad(:,:)
      zssha_adout(:,:) = ssha_ad(:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1 = DOT_PRODUCT( zsshb_tlin,  zsshb_adout  )
      zsp2_2 = DOT_PRODUCT( zsshn_tlin,  zsshn_adout  )
      zsp2_3 = DOT_PRODUCT( zssha_tlin,  zssha_adout  )

      zsp2 = zsp2_1 + zsp2_2 + zsp2_3

      ! Compare the scalar products
      ! 14 char:'12345678901234'
      cl_name = 'sshnxt_adj    '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

   END SUBROUTINE ssh_nxt_adj_tst

   !!======================================================================

END MODULE sshwzv_tam
