MODULE ldftra
   !!======================================================================
   !!                       ***  MODULE  ldftra  ***
   !! Ocean physics:  lateral diffusivity coefficient 
   !!=====================================================================
   !! History :        ! 1997-07  (G. Madec)  from inimix.F split in 2 routines
   !!   NEMO      1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!             2.0  ! 2005-11  (G. Madec)  
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_tra_init : initialization, namelist read, and parameters control
   !!   ldf_tra_c3d   : 3D eddy viscosity coefficient initialization
   !!   ldf_tra_c2d   : 2D eddy viscosity coefficient initialization
   !!   ldf_tra_c1d   : 1D eddy viscosity coefficient initialization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE ldfslp          ! ???
   USE in_out_manager  ! I/O manager
   USE ioipsl
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_tra_init   ! called by opa.F90

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra.F90 3294 2012-01-28 16:44:18Z rblod $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_tra_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_tra_init  ***
      !! 
      !! ** Purpose :   initializations of the tracer lateral mixing coeff.
      !!
      !! ** Method  :   the Eddy diffusivity and eddy induced velocity ceoff.
      !!      are defined as follows:
      !!         default option   : constant coef. aht0, aeiv0 (namelist)
      !!        'key_traldf_c1d': depth dependent coef. defined in 
      !!                            in ldf_tra_c1d routine
      !!        'key_traldf_c2d': latitude and longitude dependent coef.
      !!                            defined in ldf_tra_c2d routine
      !!        'key_traldf_c3d': latitude, longitude, depth dependent coef.
      !!                            defined in ldf_tra_c3d routine
      !!
      !!      N.B. User defined include files.  By default, 3d and 2d coef.
      !!      are set to a constant value given in the namelist and the 1d
      !!      coefficients are initialized to a hyperbolic tangent vertical
      !!      profile.
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio               ! temporary integer
      LOGICAL ::   ll_print = .FALSE.   ! =T print eddy coef. in numout
      !! 
      NAMELIST/namtra_ldf/ ln_traldf_lap  , ln_traldf_bilap,                  &
         &                 ln_traldf_level, ln_traldf_hor  , ln_traldf_iso,   &
         &                 ln_traldf_grif , ln_traldf_gdia,                   &
         &                 ln_triad_iso   , ln_botmix_grif,                   &
         &                 rn_aht_0       , rn_ahtb_0      , rn_aeiv_0,       &
         &                 rn_slpmax
      !!----------------------------------------------------------------------

      !  Define the lateral tracer physics parameters
      ! =============================================
    
      REWIND( numnam )                  ! Read Namelist namtra_ldf : Lateral physics on tracers
      READ  ( numnam, namtra_ldf )

      IF(lwp) THEN                      ! control print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_tra_init : lateral tracer physics'
         WRITE(numout,*) '~~~~~~~~~~~~ '
         WRITE(numout,*) '   Namelist namtra_ldf : lateral mixing parameters (type, direction, coefficients)'
         WRITE(numout,*) '      laplacian operator            ln_traldf_lap   = ', ln_traldf_lap
         WRITE(numout,*) '      bilaplacian operator          ln_traldf_bilap = ', ln_traldf_bilap
         WRITE(numout,*) '      iso-level                     ln_traldf_level = ', ln_traldf_level
         WRITE(numout,*) '      horizontal (geopotential)     ln_traldf_hor   = ', ln_traldf_hor
         WRITE(numout,*) '      iso-neutral                   ln_traldf_iso   = ', ln_traldf_iso
         WRITE(numout,*) '      iso-neutral (Griffies)        ln_traldf_grif  = ', ln_traldf_grif
         WRITE(numout,*) '      Griffies strmfn diagnostics   ln_traldf_gdia  = ', ln_traldf_gdia
         WRITE(numout,*) '      lateral eddy diffusivity      rn_aht_0        = ', rn_aht_0
         WRITE(numout,*) '      background hor. diffusivity   rn_ahtb_0       = ', rn_ahtb_0
         WRITE(numout,*) '      eddy induced velocity coef.   rn_aeiv_0       = ', rn_aeiv_0
         WRITE(numout,*) '      maximum isoppycnal slope      rn_slpmax       = ', rn_slpmax
         WRITE(numout,*) '      pure lateral mixing in ML     ln_triad_iso    = ', ln_triad_iso
         WRITE(numout,*) '      lateral mixing on bottom      ln_botmix_grif  = ', ln_botmix_grif
         WRITE(numout,*)
      ENDIF

      !                                ! convert DOCTOR namelist names into OLD names
      aht0  = rn_aht_0
      ahtb0 = rn_ahtb_0
      aeiv0 = rn_aeiv_0

      !                                ! Parameter control

      ! ... Check consistency for type and direction :
      !           ==> will be done in traldf module

      ! ... Space variation of eddy coefficients
      ioptio = 0
      IF(lwp) WRITE(numout,*) '          tracer mixing coef. = F( latitude, longitude)'
      ioptio = ioptio + 1
      IF( ioptio == 0 ) THEN
          IF(lwp) WRITE(numout,*) '          tracer mixing coef. = constant (default option)'
        ELSEIF( ioptio > 1 ) THEN
           CALL ctl_stop('          use only one of the following keys:',   &
             &           ' key_traldf_c3d, key_traldf_c2d, key_traldf_c1d' )
      ENDIF

      IF( ln_traldf_bilap ) THEN
         IF(lwp) WRITE(numout,*) '          biharmonic tracer diffusion'
         IF( aht0 > 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal diffusivity coef. aht0 must be negative' )
      ELSE
         IF(lwp) WRITE(numout,*) '          harmonic tracer diffusion (default)'
         IF( aht0 < 0 .AND. .NOT. lk_esopa )   CALL ctl_stop('The horizontal diffusivity coef. aht0 must be positive' )
      ENDIF


      !  Lateral eddy diffusivity and eddy induced velocity coefficients
      ! ================================================================
      CALL ldf_tra_c2d( ll_print )      ! aht = 2D coef. = F( longitude, latitude )
      !
   END SUBROUTINE ldf_tra_init

   !!----------------------------------------------------------------------
   !!                      ***  ldftra_c2d.h90  ***
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_c2d.h90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_tra_c2d( ld_print )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldftra_c2d  ***
      !!              
      !! ** Purpose :   initializations of horizontally non uniform eddy 
      !!      diffusivity coefficients
      !!
      !! ** Method :
      !!       biharmonic operator    : ahtt = defined at T-level
      !!                                ahtu,ahtv,ahtw never used
      !!       harmonic operator (ahtt never used)
      !!           iso-model level   : ahtu, ahtv defined at u-, v-points
      !!         isopycnal         : ahtu, ahtv, ahtw defined at u-, v-, w-pts
      !!         or geopotential   
      !!       eddy induced velocity
      !!           always harmonic   : aeiu, aeiv, aeiw defined at u-, v-, w-pts
      !!----------------------------------------------------------------------
      LOGICAL, INTENT (in) ::   ld_print   ! If true, print arrays in numout
      !
      INTEGER ::   ji, jj   ! dummy loop indices
      REAL(wp) ::   za00, zd_max, zeumax, zevmax, zetmax
      !!----------------------------------------------------------------------

      IF( lk_traldf_eiv ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' ldf_tra_c2d : 2D eddy diffusivity and eddy'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~   --  induced velocity coefficients'
      ELSE
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) ' ldf_tra2d : 2D eddy diffusivity coefficient'
         IF(lwp) WRITE(numout,*) ' ~~~~~~~~~~~   --'
      ENDIF

      zd_max = MAX( MAXVAL( e1t(:,:) ), MAXVAL( e2t(:,:) ) )
      IF( lk_mpp ) CALL mpp_max( zd_max )   ! max over the global domain

      ! harmonic operator : (U-, V-, W-points)
      ! ==================
      IF( ln_traldf_lap ) THEN
         !
         za00 = aht0 / zd_max
         !
         DO jj = 1, jpj 
            DO ji = 1, jpi 
               zeumax = MAX( e1u(ji,jj), e2u(ji,jj) ) 
               zevmax = MAX( e1v(ji,jj), e2v(ji,jj) ) 
               zetmax = MAX( e1t(ji,jj), e2t(ji,jj) )
               ahtu(ji,jj) = za00 * zeumax ! set ahtu = ahtv at u- and v-points, 
               ahtv(ji,jj) = za00 * zevmax ! and ahtw at w-point (idem T-point) 
               ahtw(ji,jj) = za00 * zetmax ! 
            END DO
         END DO

         CALL lbc_lnk( ahtu, 'U', 1. )   ! Lateral boundary conditions
         CALL lbc_lnk( ahtv, 'V', 1. )   ! (no change of sign)
         CALL lbc_lnk( ahtw, 'W', 1. )

         ! Special case for ORCA R2 and R4 configurations (overwrite the value of ahtu ahtv ahtw)
         ! ==============================================
         IF( cp_cfg == "orca" .AND. ( jp_cfg == 2 .OR. jp_cfg == 4 ) )   THEN
            ahtu(:,:) = aht0              ! set ahtu = ahtv at u- and v-points,
            ahtv(:,:) = aht0              ! and ahtw at w-point
            ahtw(:,:) = aht0              ! (here : no space variation)
            IF(lwp) WRITE(numout,*) '               ORCA R2 or R4 case'
            IF(lwp) WRITE(numout,*) '               Constant values used for eddy diffusivity coefficients'
            IF(lwp) WRITE(numout,*) '               Variation lat/lon only for eddy induced velocity coefficients'
         ENDIF

         ! Control print
         IF( lwp .AND. ld_print ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: ahtu array'
            CALL prihre( ahtu, jpi, jpj, 1, jpi, 1,   &
               &                         1, jpj, 1, 1.e-3, numout )
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: ahtv array'
            CALL prihre( ahtv, jpi, jpj, 1, jpi, 1,   &
               &                         1, jpj, 1, 1.e-3, numout )
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: ahtw array'
            CALL prihre( ahtw, jpi, jpj, 1, jpi, 1,   &
               &                         1, jpj, 1, 1.e-3, numout )
         ENDIF
      ENDIF
      
      ! biharmonic operator : (T-point)
      ! ====================
      IF( ln_traldf_bilap ) THEN
         ! (USER: modify ahtt following your desiderata)
         ! Here: ahm is proportional to the cube of the maximum of the gridspacing
         !       in the to horizontal direction

         zd_max = MAX( MAXVAL( e1t(:,:) ), MAXVAL( e2t(:,:) ) )
         IF( lk_mpp )   CALL mpp_max( zd_max )   ! max over the global domain

         za00 = aht0 / ( zd_max * zd_max * zd_max )
         DO jj = 1, jpj
            DO ji = 1, jpi
               zetmax = MAX( e1t(ji,jj), e2t(ji,jj) )
               ahtt(ji,jj) = za00 * zetmax * zetmax * zetmax      ! set ahtt at T-point
            END DO
         END DO

         CALL lbc_lnk( ahtt, 'T', 1. )   ! Lateral boundary conditions on ( ahtt )

         ! Control print
         IF( lwp .AND. ld_print ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: 2D ahtt array'
            CALL prihre( ahtt, jpi, jpj, 1, jpi, 1,   &
               &                         1, jpj, 1, 1.e-3, numout )
         ENDIF
      ENDIF

      ! set aeiu = aeiv at u- and v-points, and aeiw at w-point (idem T-point)
      ! (here no space variation)
      aeiu(:,:) = aeiv0
      aeiv(:,:) = aeiv0
      aeiw(:,:) = aeiv0
      
      IF( cp_cfg == "orca" .AND. jp_cfg == 4 ) THEN
         !                                 ! Cancel eiv in Gibraltar strait
         aeiu( mi0(68):mi1(71) , mj0(50):mj1(53) ) = 0.e0
         aeiv( mi0(68):mi1(71) , mj0(50):mj1(53) ) = 0.e0
         aeiw( mi0(68):mi1(71) , mj0(50):mj1(53) ) = 0.e0
         !                                 ! Cancel eiv in Mediterrannean sea
         aeiu( mi0(70):mi1(90) , mj0(49):mj1(56) ) = 0.e0
         aeiv( mi0(70):mi1(90) , mj0(49):mj1(56) ) = 0.e0
         aeiw( mi0(70):mi1(90) , mj0(49):mj1(56) ) = 0.e0
      ENDIF

      ! Lateral boundary conditions on ( aeiu, aeiv, aeiw )
      CALL lbc_lnk( aeiu, 'U', 1. )
      CALL lbc_lnk( aeiv, 'V', 1. )
      CALL lbc_lnk( aeiw, 'W', 1. )

      ! Control print
      IF( lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: aeiu array'
         CALL prihre(aeiu,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: aeiv array'
         CALL prihre(aeiv,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         WRITE(numout,*)
         WRITE(numout,*) 'inildf: aeiw array'
         CALL prihre(aeiw,jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
      ENDIF
      !
   END SUBROUTINE ldf_tra_c2d

   !!======================================================================
END MODULE ldftra
