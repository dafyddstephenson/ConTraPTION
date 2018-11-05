MODULE traadv_tam
   !!==============================================================================
   !!                       ***  MODULE  traadv_tam  ***
   !! Ocean active tracers:  advection trend -
   !!                        Tangent and Adjoint Module
   !!==============================================================================
   !! History of the direct routine:
   !!            9.0  !  05-11  (G. Madec)  Original code
   !! History of the TAM:
   !!                 !  08-06  (A. Vidard) Skeleton
   !!                 !  09-03  (F. Vigilant) Add tra_adv_ctl_tam routine
   !!                                         Allow call to tra_eiv_tam/adj
   !!                 ! 12-07   (P.-A. Bouttier) Phase with 3.4 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_adv      : compute ocean tracer advection trend
   !!   tra_adv_init  : control the different options of advection scheme
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce_tam
   USE oce
   USE ldftra_oce
   USE dom_oce         ! ocean space and time domain
   USE traadv_cen2_tam
   USE traadv_tvd ! 2016-11-29 added TVD advection scheme as an option for passive tracer transport 
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE cla_tam
   USE iom
   USE lib_mpp
   USE wrk_nemo
   USE timing
!!!later   ! 2016-05-16 Addition of EIV to background trajectory (optional)
!!!later   USE traadv_eiv
   USE tamctl

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_tan     ! routine called by steptan module
   PUBLIC   tra_adv_adj     ! routine called by stepadj module
   PUBLIC   tra_adv_init_tam ! routine called by stepadj module
   PUBLIC   tra_adv_adj_tst ! routine called by tst module

   !!* Namelist nam_traadv
   LOGICAL, PUBLIC ::   ln_traadv_cen2   = .TRUE.       ! 2nd order centered scheme flag
   LOGICAL, PUBLIC ::   ln_traadv_tvd    = .FALSE.      ! TVD scheme flag
!!!later   LOGICAL, PUBLIC ::   ln_traadv_muscl  = .FALSE.      ! MUSCL scheme flag
!!!later   LOGICAL, PUBLIC ::   ln_traadv_muscl2 = .FALSE.      ! MUSCL2 scheme flag
!!!later   LOGICAL, PUBLIC ::   ln_traadv_ubs    = .FALSE.      ! UBS scheme flag
!!!later   LOGICAL, PUBLIC ::   ln_traadv_qck    = .FALSE.      ! QUICKEST scheme flag
!!!later   LOGICAL, PUBLIC ::   ln_traadv_tvdpt  = .FALSE. ! 2016-11-29 added TVD advection scheme as an option for passive tracer transport 

   REAL,    PUBLIC :: rn_traadv_weight_h = 0._wp
   REAL,    PUBLIC :: rn_traadv_weight_v = 0._wp

   INTEGER ::   nadv   ! choice of the type of advection scheme
   

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

   SUBROUTINE tra_adv_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!              compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following nadv
      !!----------------------------------------------------------------------
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zuntl, zvntl, zwntl   ! effective velocity
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zun, zvn, zwn         ! effective velocity
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      INTEGER               ::   jk   ! dummy loop index

      !!----------------------------------------------------------------------
      !

      IF( nn_timing == 1 )  CALL timing_start('tra_adv_tan')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zun, zvn, zwn )
      CALL wrk_alloc( jpi, jpj, jpk, zuntl, zvntl, zwntl )
      !
      IF( neuler == 0 .AND. kt == nit000 ) THEN     !at nit000
         r2dtra(:) =  rdttra(:)                     !     = rdtra (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1) THEN                !at nit000 or nit000+1
         r2dtra(:) = 2._wp * rdttra(:)              !     = 2 rdttra (leapfrog)
      ENDIF
      !
      IF( nn_cla == 1 )   CALL cla_traadv_tan( kt )
      !
      !                                               !==  effective transport  ==!
      DO jk = 1, jpkm1
         zun(:,:,jk) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)                  ! eulerian transport only
         zvn(:,:,jk) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)
         zwn(:,:,jk) = e1t(:,:) * e2t(:,:)      * wn(:,:,jk)
         zuntl(:,:,jk) = e2u(:,:) * e3u(:,:,jk) * un_tl(:,:,jk)                  ! eulerian transport only
         zvntl(:,:,jk) = e1v(:,:) * e3v(:,:,jk) * vn_tl(:,:,jk)
         zwntl(:,:,jk) = e1t(:,:) * e2t(:,:)      * wn_tl(:,:,jk)
      END DO
      zun(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zvn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zwn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zuntl(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zvntl(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zwntl(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      !
      IF ( kt == nit000 ) THEN
         ! 2016-11-29 added TVD advection scheme as an option for passive tracer transport 
         IF (ln_traadv_tvd) THEN
            if (lwp) write(numout,*) 'TVD advection scheme in use for passive tracer transport: CAUTION - NONLINEAR ADVECTION SCHEME'
            nadv = 2
         ELSE
            IF(lwp) WRITE(numout,*) ' tra_adv_tam: 2nd order scheme is forced in TAM'
            nadv = 1 ! force tra_adv_cen2 for tangent
         END IF
      END IF 

      SELECT CASE ( nadv )                           ! compute advection trend and add it to general trend
      CASE ( 1 ) ;
         CALL tra_adv_cen2_tan ( kt, nit000, zun, zvn, zwn, tsn, zuntl, zvntl, zwntl, tsn_tl, tsa_tl, jpts, rn_traadv_weight_h,rn_traadv_weight_v )    ! 2nd order centered scheme
      ! 2016-11-29 added call to TVD advection scheme of nonlinear model as an option for passive tracer transport 
      CASE ( 2 );
         CALL tra_adv_tvd      ( kt, nit000, 'TRA', r2dtra, zun, zvn, zwn, tsb_tl, tsn_tl, tsa_tl, jpts )
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_tan')
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zun, zvn, zwn, zuntl, zvntl, zwntl )
      !
   END SUBROUTINE tra_adv_tan

   SUBROUTINE tra_adv_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!              compute the ocean tracer advection trend.
      !!
      !! ** Method  : - Update (ua,va) with the advection term following nadv
      !!----------------------------------------------------------------------
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zunad, zvnad, zwnad   ! effective velocity
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   zun, zvn, zwn   ! effective velocity
      INTEGER :: jk
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('tra_adv_adj')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zun, zvn, zwn )
      CALL wrk_alloc( jpi, jpj, jpk, zunad, zvnad, zwnad )
      !
      zunad(:,:,:) = 0._wp                                                     ! no transport trough the bottom
      zvnad(:,:,:) = 0._wp                                                     ! no transport trough the bottom
      zwnad(:,:,:) = 0._wp
      !
      IF( neuler == 0 .AND. kt == nit000 ) THEN     ! at nit000
         r2dtra(:) =  rdttra(:)                          ! = rdtra (restarting with Euler time stepping)
      ELSEIF( kt <= nit000 + 1) THEN                ! at nit000 or nit000+1
         r2dtra(:) = 2._wp * rdttra(:)                   ! = 2 rdttra (leapfrog)
      ENDIF
      !
      DO jk = 1, jpkm1
         zun(:,:,jk) = e2u(:,:) * e3u(:,:,jk) * un(:,:,jk)                  ! eulerian transport only
         zvn(:,:,jk) = e1v(:,:) * e3v(:,:,jk) * vn(:,:,jk)
         zwn(:,:,jk) = e1t(:,:) * e2t(:,:)      * wn(:,:,jk)
      END DO
      zun(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zvn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      zwn(:,:,jpk) = 0._wp                                                     ! no transport trough the bottom
      !
      IF ( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*) ' tra_adv_tam: 2nd order scheme is forced in TAM'
         nadv = 1 ! force tra_adv_cen2 for adjoint
      END IF
      !
      SELECT CASE ( nadv )                           ! compute advection trend and add it to general trend
      CASE ( 1 ) ; CALL tra_adv_cen2_adj ( kt, nit000, zun, zvn, zwn, tsn, zunad, zvnad, zwnad, tsn_ad, tsa_ad, jpts, rn_traadv_weight_h,rn_traadv_weight_v )    ! 2nd order centered scheme
      CASE ( 2 ); !!!2017-08-21 Added TVD advection scheme to adjoint with negative velocity fields 
         CALL tra_adv_tvd      ( kt, nit000, 'TRA', r2dtra,-1*zun,-1*zvn,-1*zwn, tsb_tl, tsn_tl, tsa_tl, jpts )
      END SELECT
      !
      DO jk = jpkm1, 1, -1
         un_ad(:,:,jk) = un_ad(:,:,jk) + e2u(:,:) * e3u(:,:,jk) *  zunad(:,:,jk)
         vn_ad(:,:,jk) = vn_ad(:,:,jk) + e1v(:,:) * e3v(:,:,jk) *  zvnad(:,:,jk)
         wn_ad(:,:,jk) = wn_ad(:,:,jk) + e1t(:,:) * e2t(:,:) *  zwnad(:,:,jk)
      END DO
      !
      IF( nn_cla == 1 )   CALL cla_traadv_adj( kt )
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_adv_adj')
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zun, zvn, zwn, zunad, zvnad, zwnad )
      !
   END SUBROUTINE tra_adv_adj
   SUBROUTINE tra_adv_adj_tst( kumadt )
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
      !!
      !! History :
      !!        ! 08-08 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      CALL tra_adv_cen2_adj_tst( kumadt )

   END SUBROUTINE tra_adv_adj_tst

   SUBROUTINE tra_adv_init_tam
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_ctl_tam  ***
      !!
      !! ** Purpose :   Control the consistency between namelist options for
      !!              tracer advection schemes and set nadv
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio
      NAMELIST/namtra_adv_tam/ ln_traadv_cen2 , &
         & ln_traadv_tvd,    &
!!!later         &                 ln_traadv_muscl, ln_traadv_muscl2, &
!!!later         &                 ln_traadv_ubs  , ln_traadv_qck, &
         &                 rn_traadv_weight_h, rn_traadv_weight_v
!!!later         &                 ln_traadv_tvdpt ! 2016-11-29 added TVD advection scheme as an option for passive tracer transport 
      
      !!----------------------------------------------------------------------

      REWIND( numnam )               ! Read Namelist nam_traadv : tracer advection scheme
      READ(numnam, namtra_adv_tam)

      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_adv_init_tam : choice/control of the tracer advection scheme'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '       Namelist nam_traadv_tam : chose a advection scheme for tracers'
         WRITE(numout,*) '          2nd order advection scheme     ln_traadv_cen2   = ', ln_traadv_cen2
         WRITE(numout,*) '          TVD advection scheme           ln_traadv_tvd    = ', ln_traadv_tvd
!!!later         WRITE(numout,*) '          MUSCL  advection scheme        ln_traadv_muscl  = ', ln_traadv_muscl
!!!later         WRITE(numout,*) '          MUSCL2 advection scheme        ln_traadv_muscl2 = ', ln_traadv_muscl2
!!!later         WRITE(numout,*) '          UBS    advection scheme        ln_traadv_ubs    = ', ln_traadv_ubs
!!!later         WRITE(numout,*) '          QUICKEST advection scheme      ln_traadv_qck    = ', ln_traadv_qck
!!!later         WRITE(numout,*) ' TVD advection scheme for passive tracer transport = ', ln_traadv_tvdpt ! 2016-11-29 added TVD advection scheme as an option for passive tracer transport 
         IF( rn_traadv_weight_h < 0.0_wp) THEN 
            WRITE(numout,*) '         Weighted-mean horizontal advection scheme selected '
         ELSE
            WRITE(numout,*) 'weighted centred v. upstream advection scheme, horizontal upstream weighting = ', rn_traadv_weight_h*100, '%'
         ENDIF
         IF( rn_traadv_weight_v < 0.0_wp) THEN 
            WRITE(numout,*) '         Weighted-mean vertical advection scheme selected '
         ELSE
            WRITE(numout,*) 'weighted centred v. upstream advection scheme, vertical upstream weighting = ', rn_traadv_weight_v*100, '%'
         ENDIF

      ENDIF

      ioptio = 0                      ! Parameter control
      IF( ln_traadv_cen2   )   ioptio = ioptio + 1
      IF( ln_traadv_tvd    )   ioptio = ioptio + 1
!!!later      IF( ln_traadv_muscl  )   ioptio = ioptio + 1
!!!later      IF( ln_traadv_muscl2 )   ioptio = ioptio + 1
!!!later      IF( ln_traadv_ubs    )   ioptio = ioptio + 1
!!!later      IF( ln_traadv_qck    )   ioptio = ioptio + 1
!!!later      IF( ln_traadv_tvdpt  )   ioptio =   ioptio + 1
      IF( lk_esopa         )   ioptio =          1

      IF( ioptio /= 1 )   CALL ctl_stop( 'Choose ONE advection scheme in namelist nam_traadv' )

      IF( nn_cla == 1 .AND. .NOT. ln_traadv_cen2 )   &
         &                CALL ctl_stop( 'cross-land advection only with 2nd order advection scheme' )

      !                              ! Set nadv
      IF( ln_traadv_cen2   )   nadv =  1
      IF( ln_traadv_tvd    )   nadv =  2
!!!later      ! IF( ln_traadv_muscl  )   nadv =  3
!!!later      ! IF( ln_traadv_muscl2 )   nadv =  4
!!!later      ! IF( ln_traadv_ubs    )   nadv =  5
!!!later      ! IF( ln_traadv_qck    )   nadv =  6
      IF( lk_esopa         )   nadv = -1

      IF(lwp) THEN                   ! Print the choice
         WRITE(numout,*)
         IF( nadv ==  1 )   WRITE(numout,*) '         2nd order scheme is used'
         IF( nadv ==  2 )   WRITE(numout,*) '         TVD scheme: WARNING non-linear advection scheme selected'
!!!later         IF( nadv ==  3 )   WRITE(numout,*) '         MUSCL     Not Available in NEMO TAM'
!!!later         IF( nadv ==  4 )   WRITE(numout,*) '         MUSCL2    Not Available in NEMO TAM'
!!!later         IF( nadv ==  5 )   WRITE(numout,*) '         UBS       Not Available in NEMO TAM'
!!!later         IF( nadv ==  6 )   WRITE(numout,*) '         QUICKEST  Not Available in NEMO TAM'
         IF( nadv == -1 )   WRITE(numout,*) '         esopa test: Not Available in NEMO TAM'
      ENDIF
      !
   END SUBROUTINE tra_adv_init_tam


  !!======================================================================
END MODULE traadv_tam
