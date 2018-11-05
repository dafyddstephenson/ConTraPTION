MODULE dynadv_tam
   !!==============================================================================
   !!                       ***  MODULE  dynadv_tam  ***
   !! Ocean active tracers:  advection scheme control
   !!                        Tangent and Adjoint module
   !!==============================================================================
   !! History of the direct module:
   !!           9.0  ! 2006-11  (G. Madec)  Original code
   !! History of the TAM module:
   !!           9.0  ! 2008-08  (A. Vidard) first version
   !!   NEMO    3.2  ! 2010-04 (F. Vigilant) 3.2 version
   !!   NEMO    3.4  ! 2012-07 (P.-A. Bouttier) 3.4 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_adv      : compute the momentum advection trend
   !!   dyn_adv_ctl  : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE oce
   USE dom_oce
   USE oce_tam
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE dynadv
   USE dynadv_cen2_tam
   USE dynadv_ubs_tam
   USE dynkeg_tam
   USE dynzad_tam
   USE sshwzv_tam
   USE sshwzv
   USE divcur
   USE divcur_tam
   USE in_out_manager
   USE lib_mpp
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_adv_tan     ! routine called by steptan module
   PUBLIC dyn_adv_adj     ! routine called by stepadj module
   PUBLIC dyn_adv_adj_tst ! routine called by the tst module
   PUBLIC dyn_adv_init_tam

   INTEGER    ::   nadv   ! choice of the formulation and scheme for the advection
   LOGICAL    ::   lfirst=.TRUE.
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

   SUBROUTINE dyn_adv_tan( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the ocean momentum advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following nadv
      !!      NB: in flux form advection (ln_dynadv_cen2 or ln_dynadv_ubs=T)
      !!      a metric term is add to the coriolis term while in vector form
      !!      it is the relative vorticity which is added to coriolis term
      !!      (see dynvor module).
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_adv_tan')
      !
      SELECT CASE ( nadv )                     ! compute advection trend and add it to general trend
      CASE ( 0 )
                      CALL dyn_keg_tan     ( kt )    ! vector form : horizontal gradient of kinetic energy
                      CALL dyn_zad_tan     ( kt )    ! vector form : vertical advection
      CASE ( 1 )
                      CALL dyn_adv_cen2_tan( kt )    ! 2nd order centered scheme
      CASE ( 2 )
                      CALL dyn_adv_ubs_tan ( kt )    ! 3rd order UBS      scheme
      !
      CASE (-1 )                                 ! esopa: test all possibility with control print
                      CALL dyn_keg_tan     ( kt )
                      CALL dyn_zad_tan     ( kt )
                      CALL dyn_adv_cen2_tan( kt )
                      CALL dyn_adv_ubs_tan ( kt )
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_adv_tan')
      !
   END SUBROUTINE dyn_adv_tan

   SUBROUTINE dyn_adv_adj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!            compute the ocean momentum advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following nadv
      !!      NB: in flux form advection (ln_dynadv_cen2 or ln_dynadv_ubs=T)
      !!      a metric term is add to the coriolis term while in vector form
      !!      it is the relative vorticity which is added to coriolis term
      !!      (see dynvor module).
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_adv_adj')
      !
      SELECT CASE ( nadv )                     ! compute advection trend and add it to general trend
      CASE ( 0 )
                      CALL dyn_zad_adj     ( kt )    ! vector form : vertical advection
                      CALL dyn_keg_adj     ( kt )    ! vector form : horizontal gradient of kinetic energy
      CASE ( 1 )
         IF (lwp) WRITE(numout,*) 'dyn_adv_cen2_adj not available yet'
         CALL abort
      CASE ( 2 )
         IF (lwp) WRITE(numout,*) 'dyn_adv_ubs_adj not available yet'
         CALL abort
      CASE (-1 )                                 ! esopa: test all possibility with control print
                      CALL dyn_zad_adj     ( kt )
                      CALL dyn_keg_adj     ( kt )
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_adv_adj')
      !
   END SUBROUTINE dyn_adv_adj

   SUBROUTINE dyn_adv_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dyn_adv_adj_tst ***
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
         & jk,    &
         & jt

      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zun_tlin,     & ! Tangent input:  now  u-velocity
         & zvn_tlin,     & ! Tangent input:  now  v-velocity
         & zwn_tlin,     & ! Tangent input:  now  w-velocity
         & zua_tlin,     & ! Tangent input: after u-velocity
         & zva_tlin,     & ! Tangent input: after u-velocity
         & zua_tlout,    & ! Tangent output:after u-velocity
         & zva_tlout,    & ! Tangent output:after v-velocity
         & zua_adin,     & ! adjoint input: after u-velocity
         & zva_adin,     & ! adjoint input: after v-velocity
         & zun_adout,    & ! adjoint output: now  u-velocity
         & zvn_adout,    & ! adjoint output: now  v-velocity
         & zwn_adout,    & ! adjoint output: now  u-velocity
         & zua_adout,    & ! adjoint output:after v-velocity
         & zva_adout,    & ! adjoint output:after u-velocity
         & zuvw            ! 3D random field for u, v and w

      REAL(KIND=wp) ::   &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4,       &
         & zsp2_5

      CHARACTER(LEN=14) :: cl_name


      ! Allocate memory

      ALLOCATE( &
         & zun_tlin(jpi,jpj,jpk),     &
         & zvn_tlin(jpi,jpj,jpk),     &
         & zwn_tlin(jpi,jpj,jpk),     &
         & zua_tlin(jpi,jpj,jpk),     &
         & zva_tlin(jpi,jpj,jpk),     &
         & zua_tlout(jpi,jpj,jpk),    &
         & zva_tlout(jpi,jpj,jpk),    &
         & zua_adin(jpi,jpj,jpk),     &
         & zva_adin(jpi,jpj,jpk),     &
         & zun_adout(jpi,jpj,jpk),    &
         & zvn_adout(jpi,jpj,jpk),    &
         & zwn_adout(jpi,jpj,jpk),    &
         & zua_adout(jpi,jpj,jpk),    &
         & zva_adout(jpi,jpj,jpk),    &
         & zuvw(jpi,jpj,jpk)           &
         & )


      !===================================================================================
      ! 1)  dx = ( un_tl, vn_tl, ua_tl, va_tl )        --> dynkeg
      ! and dy = ( ua_tl, va_tl )
      ! 2)  dx = ( un_tl, vn_tl, wn_tl, ua_tl, va_tl ) --> dynkeg
      ! and dy = ( ua_tl, va_tl )
      ! 3)  dx = ( un_tl, vn_tl, wn_tl, ua_tl, va_tl ) --> dynadv
      ! and dy = ( ua_tl, va_tl )
      !======================================================================

      DO jt = 1, 3
      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          zun_tlin(:,:,:)  = 0.0_wp
          zvn_tlin(:,:,:)  = 0.0_wp
          zwn_tlin(:,:,:)  = 0.0_wp
          zua_tlin(:,:,:)  = 0.0_wp
          zva_tlin(:,:,:)  = 0.0_wp
          zua_tlout(:,:,:) = 0.0_wp
          zva_tlout(:,:,:) = 0.0_wp
          zua_adin(:,:,:)  = 0.0_wp
          zva_adin(:,:,:)  = 0.0_wp
          zun_adout(:,:,:) = 0.0_wp
          zvn_adout(:,:,:) = 0.0_wp
          zwn_adout(:,:,:) = 0.0_wp
          zua_adout(:,:,:) = 0.0_wp
          zva_adout(:,:,:) = 0.0_wp
          zuvw(:,:,:)      = 0.0_wp

          un_tl(:,:,:) = 0.0_wp
          vn_tl(:,:,:) = 0.0_wp
          wn_tl(:,:,:) = 0.0_wp
          ua_tl(:,:,:) = 0.0_wp
          va_tl(:,:,:) = 0.0_wp
          un_ad(:,:,:) = 0.0_wp
          vn_ad(:,:,:) = 0.0_wp
          wn_ad(:,:,:) = 0.0_wp
          ua_ad(:,:,:) = 0.0_wp
          va_ad(:,:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

         CALL grid_random(  zuvw, 'U', 0.0_wp, stdu )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zun_tlin(ji,jj,jk) = zuvw(ji,jj,jk)
               END DO
            END DO
         END DO
         CALL grid_random(  zuvw, 'V', 0.0_wp, stdv )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zvn_tlin(ji,jj,jk) = zuvw(ji,jj,jk)
               END DO
            END DO
         END DO

         CALL grid_random(  zuvw, 'W', 0.0_wp, stdw )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zwn_tlin(ji,jj,jk) = zuvw(ji,jj,jk)
               END DO
            END DO
         END DO

         CALL grid_random(  zuvw, 'U', 0.0_wp, stdu )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zua_tlin(ji,jj,jk) = zuvw(ji,jj,jk)
               END DO
            END DO
         END DO

         CALL grid_random(  zuvw, 'V', 0.0_wp, stdv )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zva_tlin(ji,jj,jk) = zuvw(ji,jj,jk)
               END DO
            END DO
         END DO

         un_tl(:,:,:) = zun_tlin(:,:,:)
         vn_tl(:,:,:) = zvn_tlin(:,:,:)

         wn_tl(:,:,:) = zwn_tlin(:,:,:)
         ua_tl(:,:,:) = zua_tlin(:,:,:)
         va_tl(:,:,:) = zva_tlin(:,:,:)

         SELECT CASE ( jt )
         CASE ( 1 )
            CALL dyn_adv_init_tam
            CALL dyn_keg_tan( nit000 )
         CASE ( 2 )
            CALL dyn_adv_init_tam
            CALL dyn_zad_tan( nit000 )
         CASE ( 3 )
            CALL dyn_adv_tan ( nit000 )
         END SELECT

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

         SELECT CASE ( jt )
         CASE ( 1 )
            CALL dyn_keg_adj( nitend )
         CASE ( 2 )
            CALL dyn_zad_adj( nitend )
         CASE ( 3 )
            CALL dyn_adv_adj( nitend )
         END SELECT

         zun_adout(:,:,:) = un_ad(:,:,:)
         zvn_adout(:,:,:) = vn_ad(:,:,:)
         zwn_adout(:,:,:) = wn_ad(:,:,:)
         zua_adout(:,:,:) = ua_ad(:,:,:)
         zva_adout(:,:,:) = va_ad(:,:,:)

         zsp2_1 = DOT_PRODUCT( zun_tlin, zun_adout )
         zsp2_2 = DOT_PRODUCT( zvn_tlin, zvn_adout )
         IF (jt .EQ. 1) THEN
            zsp2_3 = 0.0_wp
         ELSE
            zsp2_3 = DOT_PRODUCT( zwn_tlin, zwn_adout )
         ENDIF
         zsp2_4 = DOT_PRODUCT( zua_tlin, zua_adout )
         zsp2_5 = DOT_PRODUCT( zva_tlin, zva_adout )
         zsp2   = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5

         ! Compare the scalar products

         SELECT CASE ( jt )
         CASE ( 1 )
            cl_name = 'dyn_keg_adj   '
         CASE ( 2 )
            cl_name = 'dyn_zad_adj   '
         CASE ( 3 )
            cl_name = 'dyn_adv_adj   '
         END SELECT

         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      ENDDO

      DEALLOCATE( &
         & zun_tlin,     &
         & zvn_tlin,     &
         & zwn_tlin,     &
         & zua_tlin,     &
         & zva_tlin,     &
         & zua_tlout,    &
         & zva_tlout,    &
         & zua_adin,     &
         & zva_adin,     &
         & zun_adout,    &
         & zvn_adout,    &
         & zwn_adout,    &
         & zua_adout,    &
         & zva_adout,    &
         & zuvw          &
         & )

   END SUBROUTINE dyn_adv_adj_tst
  !!======================================================================
   SUBROUTINE dyn_adv_init_tam
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_adv_ctl  ***
      !!
      !! ** Purpose :   Control the consistency between namelist options for
      !!              momentum advection formulation & scheme and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio

      NAMELIST/namdyn_adv/ ln_dynadv_vec, ln_dynadv_cen2 , ln_dynadv_ubs
      !!----------------------------------------------------------------------

      IF (lfirst) THEN

         REWIND ( numnam )               ! Read Namelist namdyn_adv : momentum advection scheme
         READ   ( numnam, namdyn_adv )

         IF(lwp) THEN                    ! Namelist print
            WRITE(numout,*)
            WRITE(numout,*) 'dyn_adv_init : choice/control of the momentum advection scheme'
            WRITE(numout,*) '~~~~~~~~~~~'
            WRITE(numout,*) '       Namelist namdyn_adv : chose a advection formulation & scheme for momentum'
            WRITE(numout,*) '          Vector/flux form (T/F)             ln_dynadv_vec  = ', ln_dynadv_vec
            WRITE(numout,*) '          2nd order centred advection scheme ln_dynadv_cen2 = ', ln_dynadv_cen2
            WRITE(numout,*) '          3rd order UBS advection scheme     ln_dynadv_ubs  = ', ln_dynadv_ubs
         ENDIF

         ioptio = 0                      ! Parameter control
         IF( ln_dynadv_vec  )   ioptio = ioptio + 1
         IF( ln_dynadv_cen2 )   ioptio = ioptio + 1
         IF( ln_dynadv_ubs  )   ioptio = ioptio + 1
         IF( lk_esopa       )   ioptio =          1

         IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE advection scheme in namelist namdyn_adv' )

      !                               ! Set nadv
         IF( ln_dynadv_vec  )   nadv =  0
         IF( ln_dynadv_cen2 )   nadv =  1
         IF( ln_dynadv_ubs  )   nadv =  2
         IF( lk_esopa       )   nadv = -1

         IF(lwp) THEN                    ! Print the choice
            WRITE(numout,*)
            IF( nadv ==  0 )   WRITE(numout,*) '         vector form : keg + zad + vor is used'
            IF( nadv ==  1 )   WRITE(numout,*) '         flux form   : 2nd order scheme is used'
            IF( nadv ==  2 )   WRITE(numout,*) '         flux form   : UBS       scheme is used'
            IF( nadv == -1 )   WRITE(numout,*) '         esopa test: use all advection formulation'
         ENDIF
         !
         lfirst = .FALSE.
      END IF
   END SUBROUTINE dyn_adv_init_tam
END MODULE dynadv_tam
