MODULE ldfdyn
   !!======================================================================
   !!                       ***  MODULE  ldfdyn  ***
   !! Ocean physics:  lateral viscosity coefficient 
   !!=====================================================================
   !! History :  OPA  ! 1997-07  (G. Madec)  multi dimensional coefficients
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ldf_dyn_init : initialization, namelist read, and parameters control
   !!   ldf_dyn_c3d   : 3D eddy viscosity coefficient initialization
   !!   ldf_dyn_c2d   : 2D eddy viscosity coefficient initialization
   !!   ldf_dyn_c1d   : 1D eddy viscosity coefficient initialization
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers   
   USE dom_oce         ! ocean space and time domain 
   USE ldfdyn_oce      ! ocean dynamics lateral physics
   USE phycst          ! physical constants
   USE ldfslp          ! ???
   USE ioipsl
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation

   IMPLICIT NONE
   PRIVATE

   PUBLIC   ldf_dyn_init   ! called by opa.F90

  INTERFACE ldf_zpf
     MODULE PROCEDURE ldf_zpf_1d, ldf_zpf_1d_3d, ldf_zpf_3d
  END INTERFACE

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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfdyn.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE ldf_dyn_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_init  ***
      !!                   
      !! ** Purpose :   set the horizontal ocean dynamics physics
      !!
      !! ** Method  :  
      !!      -  default option : ahm = constant coef. = rn_ahm_0 (namelist)
      !!      - 'key_dynldf_c1d': ahm = F(depth)                     see ldf_dyn_c1d.h90
      !!      - 'key_dynldf_c2d': ahm = F(latitude,longitude)        see ldf_dyn_c2d.h90
      !!      - 'key_dynldf_c3d': ahm = F(latitude,longitude,depth)  see ldf_dyn_c3d.h90
      !!
      !!      N.B. User defined include files.  By default, 3d and 2d coef.
      !!      are set to a constant value given in the namelist and the 1d
      !!      coefficients are initialized to a hyperbolic tangent vertical
      !!      profile.
      !!
      !! Reference :   Madec, G. and M. Imbard, 1996: Climate Dynamics, 12, 381-388.
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio         ! ???
      LOGICAL ::   ll_print = .FALSE.    ! Logical flag for printing viscosity coef.
      !! 
      NAMELIST/namdyn_ldf/ ln_dynldf_lap  , ln_dynldf_bilap,                  &
         &                 ln_dynldf_level, ln_dynldf_hor  , ln_dynldf_iso,   &
         &                 rn_ahm_0_lap   , rn_ahmb_0      , rn_ahm_0_blp
      !!----------------------------------------------------------------------

      REWIND( numnam )                  ! Read Namelist namdyn_ldf : Lateral physics
      READ  ( numnam, namdyn_ldf )

      IF(lwp) THEN                      ! Parameter print
         WRITE(numout,*)
         WRITE(numout,*) 'ldf_dyn : lateral momentum physics'
         WRITE(numout,*) '~~~~~~~'
         WRITE(numout,*) '   Namelist nam_dynldf : set lateral mixing parameters'
         WRITE(numout,*) '      laplacian operator                      ln_dynldf_lap   = ', ln_dynldf_lap
         WRITE(numout,*) '      bilaplacian operator                    ln_dynldf_bilap = ', ln_dynldf_bilap
         WRITE(numout,*) '      iso-level                               ln_dynldf_level = ', ln_dynldf_level
         WRITE(numout,*) '      horizontal (geopotential)               ln_dynldf_hor   = ', ln_dynldf_hor
         WRITE(numout,*) '      iso-neutral                             ln_dynldf_iso   = ', ln_dynldf_iso
         WRITE(numout,*) '      horizontal laplacian eddy viscosity     rn_ahm_0_lap    = ', rn_ahm_0_lap
         WRITE(numout,*) '      background viscosity                    rn_ahmb_0       = ', rn_ahmb_0
         WRITE(numout,*) '      horizontal bilaplacian eddy viscosity   rn_ahm_0_blp    = ', rn_ahm_0_blp
      ENDIF

      ahm0     = rn_ahm_0_lap              ! OLD namelist variables defined from DOCTOR namelist variables
      ahmb0    = rn_ahmb_0
      ahm0_blp = rn_ahm_0_blp

      ! ... check of lateral diffusive operator on tracers
      !           ==> will be done in trazdf module

      ! ... Space variation of eddy coefficients
      ioptio = 0
      IF(lwp) WRITE(numout,*) '   momentum mixing coef. = F( latitude, longitude, depth)'
      ioptio = ioptio+1
      IF( ioptio == 0 ) THEN
          IF(lwp) WRITE(numout,*) '   momentum mixing coef. = constant  (default option)'
        ELSEIF( ioptio > 1 ) THEN
           CALL ctl_stop( 'use only one of the following keys: key_dynldf_c3d, key_dynldf_c2d, key_dynldf_c1d' )
      ENDIF


      IF( ln_dynldf_bilap ) THEN
         IF(lwp) WRITE(numout,*) '   biharmonic momentum diffusion'
         IF( .NOT. ln_dynldf_lap ) ahm0 = ahm0_blp   ! Allow spatially varying coefs, which use ahm0 as input
         IF( ahm0_blp > 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal viscosity coef. ahm0 must be negative' )
      ELSE
         IF(lwp) WRITE(numout,*) '   harmonic momentum diff. (default)'
         IF( ahm0 < 0 .AND. .NOT. lk_esopa )   CALL ctl_stop( 'The horizontal viscosity coef. ahm0 must be positive' )
      ENDIF


      ! Lateral eddy viscosity
      ! ======================
      CALL ldf_dyn_c3d( ll_print )   ! ahm = 3D coef. = F( longitude, latitude, depth )
      !
   END SUBROUTINE ldf_dyn_init

   !!----------------------------------------------------------------------
   !!                        ***  ldfdyn_c3d.h90  ***
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfdyn_c3d.h90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   'key_dynldf_c3d'             3D lateral eddy viscosity coefficients
   !!----------------------------------------------------------------------

   SUBROUTINE ldf_dyn_c3d( ld_print )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_c3d  ***
      !!                   
      !! ** Purpose :   initializations of the horizontal ocean physics
      !!
      !! ** Method  :   3D eddy viscosity coef. ( longitude, latitude, depth )
      !!       laplacian operator   : ahm1, ahm2 defined at T- and F-points
      !!                              ahm2, ahm4 never used
      !!       bilaplacian operator : ahm1, ahm2 never used
      !!                           :  ahm3, ahm4 defined at U- and V-points
      !!       ??? explanation of the default is missing
      !!----------------------------------------------------------------------
      USE ldftra_oce, ONLY :   aht0
      !!
      LOGICAL, INTENT (in) ::   ld_print   ! If true, output arrays on numout
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zr = 0.2     ! maximum of the reduction factor at the bottom ocean ( 0 < zr < 1 )
      REAL(wp) ::   zh = 500.    ! depth of at which start the reduction ( > dept(1) )
      REAL(wp) ::   zd_max       ! maximum grid spacing over the global domain
      REAL(wp) ::   za00, zc, zd, zetmax, zefmax, zeumax, zevmax   ! local scalars
      REAL(wp), POINTER, DIMENSION(:) :: zcoef   
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpk, zcoef )
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ldf_dyn_c3d : 3D lateral eddy viscosity coefficient'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'

      
      ! Set ahm1 and ahm2  ( T- and F- points) (used for laplacian operators
      ! =================                       whatever its orientation is)
      IF( ln_dynldf_lap ) THEN
         ! define ahm1 and ahm2 at the right grid point position
         ! (USER: modify ahm1 and ahm2 following your desiderata)

         zd_max = MAX( MAXVAL( e1t(:,:) ), MAXVAL( e2t(:,:) ) )
         IF( lk_mpp )   CALL mpp_max( zd_max )   ! max over the global domain

         IF(lwp) WRITE(numout,*) '              laplacian operator: ahm proportional to e1'
         IF(lwp) WRITE(numout,*) '              maximum grid-spacing = ', zd_max, ' maximum value for ahm = ', ahm0

         za00 = ahm0 / zd_max

         IF( ln_dynldf_iso ) THEN
            IF(lwp) WRITE(numout,*) '              Caution, as implemented now, the isopycnal part of momentum'
            IF(lwp) WRITE(numout,*) '                 mixing use aht0 as eddy viscosity coefficient. Thus, it is'
            IF(lwp) WRITE(numout,*) '                 uniform and you must be sure that your ahm is greater than'
            IF(lwp) WRITE(numout,*) '                 aht0 everywhere in the model domain.'
         ENDIF

         CALL ldf_zpf( .TRUE. , 1000., 500., 0.25, gdept(:,:,:), ahm1 )   ! vertical profile
         CALL ldf_zpf( .TRUE. , 1000., 500., 0.25, gdept(:,:,:), ahm2 )   ! vertical profile
         DO jk = 1,jpk
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zetmax = MAX( e1t(ji,jj), e2t(ji,jj) )
                  zefmax = MAX( e1f(ji,jj), e2f(ji,jj) )
                  ahm1(ji,jj,jk) = za00 * zetmax * ahm1(ji,jj,jk)
                  ahm2(ji,jj,jk) = za00 * zefmax * ahm2(ji,jj,jk)
               END DO
            END DO
         END DO


         ! Special case for ORCA R1, R2 and R4 configurations (overwrite the value of ahm1 ahm2)
         ! ==============================================
         IF( cp_cfg == "orca" .AND. ( jp_cfg == 1 .OR. jp_cfg == 2 .OR. jp_cfg == 4 ) ) THEN 
            IF(lwp) WRITE(numout,*)
            IF(lwp) WRITE(numout,*) '              ORCA R1, R2 or R4: overwrite the previous definition of ahm'
            IF(lwp) WRITE(numout,*) '              ================='
            CALL ldf_dyn_c3d_orca( ld_print )
         ENDIF

      ENDIF
      
      ! Control print
      IF(lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '         3D ahm1 array (k=1)'
         CALL prihre( ahm1(:,:,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1.e-3, numout )
         WRITE(numout,*)
         WRITE(numout,*) '         3D ahm2 array (k=1)'
         CALL prihre( ahm2(:,:,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1.e-3, numout )
      ENDIF


      ! ahm3 and ahm4 at U- and V-points (used for bilaplacian operator
      ! ================================  whatever its orientation is)
      ! (USER: modify ahm3 and ahm4 following your desiderata)
      ! Here: ahm is proportional to the cube of the maximum of the gridspacing
      !       in the to horizontal direction

      IF( ln_dynldf_bilap ) THEN

         zd_max = MAX( MAXVAL( e1u(:,:) ), MAXVAL( e2u(:,:) ) )
         IF( lk_mpp )   CALL mpp_max( zd_max )   ! max over the global domain

         IF(lwp) WRITE(numout,*) '              bi-laplacian operator: ahm proportional to e1**3 '
         IF(lwp) WRITE(numout,*) '              maximum grid-spacing = ', zd_max, ' maximum value for ahm = ', ahm0

         za00 = ahm0_blp / ( zd_max * zd_max * zd_max )
         DO jj = 1, jpj
            DO ji = 1, jpi
               zeumax = MAX( e1u(ji,jj), e2u(ji,jj) )
               zevmax = MAX( e1v(ji,jj), e2v(ji,jj) )
               ahm3(ji,jj,1) = za00 * zeumax * zeumax * zeumax
               ahm4(ji,jj,1) = za00 * zevmax * zevmax * zevmax
            END DO
         END DO

         zh = MAX( zh, gdept(1,1,1) )   ! at least the first reach ahm0
         IF( ln_zco ) THEN               ! z-coordinate, same profile everywhere
            IF(lwp) WRITE(numout,'(36x," ahm ", 7x)')
            DO jk = 1, jpk
               IF( gdept(1,1,jk) <= zh ) THEN
                  zcoef(jk) = 1.e0
               ELSE
                  zcoef(jk) = 1.e0 + ( zr - 1.e0 )   &
                     &               * (  1. - EXP( ( gdept(1,1,jk   ) - zh ) / zh )  )   &
                     &               / (  1. - EXP( ( gdept(1,1,jpkm1) - zh ) / zh )  )
               ENDIF
               ahm3(:,:,jk) = ahm3(:,:,1) * zcoef(jk)
               ahm4(:,:,jk) = ahm4(:,:,1) * zcoef(jk)
               IF(lwp) WRITE(numout,'(34x,E7.2,8x,i3)') zcoef(jk) * ahm0, jk
            END DO
         ELSE                            ! partial steps or s-ccordinate
            zc = MAXVAL( gdept(:,:,jpkm1) )
            IF( lk_mpp )   CALL mpp_max( zc )   ! max over the global domain

            zc = 1. / (  1. - EXP( ( zc - zh ) / zh )  )
            DO jk = 2, jpkm1
               DO jj = 1, jpj
                  DO ji = 1, jpi
                     IF( gdept(ji,jj,jk) <= zh ) THEN
                        ahm3(ji,jj,jk) = ahm3(ji,jj,1)
                        ahm4(ji,jj,jk) = ahm4(ji,jj,1)
                     ELSE
                        zd = 1.e0 + ( zr - 1.e0 ) * (  1. - EXP( ( gdept(ji,jj,jk) - zh ) / zh )  ) * zc
                        ahm3(ji,jj,jk) = ahm3(ji,jj,1) * zd
                        ahm4(ji,jj,jk) = ahm4(ji,jj,1) * zd
                     ENDIF
                  END DO
               END DO
            END DO
            ahm3(:,:,jpk) = ahm3(:,:,jpkm1)
            ahm4(:,:,jpk) = ahm4(:,:,jpkm1)
            IF(lwp) WRITE(numout,'(36x," ahm ", 7x)')
            DO jk = 1, jpk
               IF(lwp) WRITE(numout,'(30x,E10.2,8x,i3)') ahm3(1,1,jk), jk
            END DO
         ENDIF

         ! Control print
         IF( lwp .AND. ld_print ) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: ahm3 array at level 1'
            CALL prihre(ahm3(:,:,1  ),jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
            WRITE(numout,*)
            WRITE(numout,*) 'inildf: ahm4 array at level 1'
            CALL prihre(ahm4(:,:,jpk),jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         ENDIF
      ENDIF
      !
      CALL wrk_dealloc( jpk, zcoef )
      !
   END SUBROUTINE ldf_dyn_c3d


   SUBROUTINE ldf_dyn_c3d_orca( ld_print )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_dyn_c3d  ***
      !!                   
      !! ** Purpose :   ORCA R1, R2 and R4 only
      !!
      !! ** Method  :   blah blah blah ....
      !!----------------------------------------------------------------------
      USE ldftra_oce, ONLY:   aht0
      !!
      LOGICAL, INTENT(in) ::   ld_print   ! If true, output arrays on numout
      !!
      INTEGER ::   ji, jj, jk, jn      ! dummy loop indices
      INTEGER ::   ii0, ii1, ij0, ij1  ! local integers
      INTEGER ::   inum, iim, ijm      ! 
      INTEGER ::   ifreq, il1, il2, ij, ii
      REAL(wp) ::   zahmeq, zcoff, zcoft, zmsk   ! local scalars
      REAL(wp) ::   zemax , zemin, zeref, zahmm
      CHARACTER (len=15) ::   clexp
      INTEGER , POINTER, DIMENSION(:,:)  :: icof
      INTEGER , POINTER, DIMENSION(:,:)  :: idata
      REAL(wp), POINTER, DIMENSION(:  )  :: zcoef   
      REAL(wp), POINTER, DIMENSION(:,:)  :: zahm0
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi   , jpj   , icof  )
      CALL wrk_alloc( jpidta, jpjdta, idata )
      CALL wrk_alloc( jpk   ,         zcoef )
      CALL wrk_alloc( jpi   , jpj   , zahm0 )
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'ldfdyn_c3d_orca : 3D eddy viscosity coefficient'
      IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~'
      IF(lwp) WRITE(numout,*) '        orca R1, R2 or R4 configuration: reduced in the surface Eq. strip '

      ! Read 2d integer array to specify western boundary increase in the
      ! ===================== equatorial strip (20N-20S) defined at t-points

      CALL ctl_opn( inum, 'ahmcoef', 'OLD', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp )
      READ(inum,9101) clexp, iim, ijm
      READ(inum,'(/)')
      ifreq = 40
      il1 = 1
      DO jn = 1, jpidta/ifreq+1
         READ(inum,'(/)')
         il2 = MIN( jpidta, il1+ifreq-1 )
         READ(inum,9201) ( ii, ji = il1, il2, 5 )
         READ(inum,'(/)')
         DO jj = jpjdta, 1, -1
            READ(inum,9202) ij, ( idata(ji,jj), ji = il1, il2 )
         END DO
         il1 = il1 + ifreq
      END DO
      
      DO jj = 1, nlcj
         DO ji = 1, nlci
            icof(ji,jj) = idata( mig(ji), mjg(jj) )
         END DO
      END DO
      DO jj = nlcj+1, jpj
         DO ji = 1, nlci
            icof(ji,jj) = icof(ji,nlcj)
         END DO
      END DO
      DO jj = 1, jpj
         DO ji = nlci+1, jpi
            icof(ji,jj) = icof(nlci,jj)
         END DO
      END DO
      
9101  FORMAT(1x,a15,2i8)
9201  FORMAT(3x,13(i3,12x))
9202  FORMAT(i3,41i3)
      
      ! Set ahm1 and ahm2
      ! =================
      
      ! define ahm1 and ahm2 at the right grid point position
      ! (USER: modify ahm1 and ahm2 following your desiderata)
      ! biharmonic : ahm1 (ahm2) defined at u- (v-) point
      ! harmonic   : ahm1 (ahm2) defined at t- (f-) point
      
      ! first level : as for 2D coefficients
      
      ! Decrease ahm to zahmeq m2/s in the tropics
      ! (from 90 to 20 degre: ahm = constant
      ! from 20 to  2.5 degre: ahm = decrease in (1-cos)/2
      ! from  2.5 to  0 degre: ahm = constant
      ! symmetric in the south hemisphere)
      
      IF( jp_cfg == 4 )   THEN
         zahmeq = 5.0 * aht0
         zahmm  = min( 160000.0, ahm0)
         zemax = MAXVAL ( e1t(:,:) * e2t(:,:), tmask(:,:,1) .GE. 0.5 )
         zemin = MINVAL ( e1t(:,:) * e2t(:,:), tmask(:,:,1) .GE. 0.5 )
         zeref = MAXVAL ( e1t(:,:) * e2t(:,:),   &
             &   tmask(:,:,1) .GE. 0.5 .AND. ABS(gphit(:,:)) .GT. 50. )
 
         DO jj = 1, jpj
           DO ji = 1, jpi
              zmsk = e1t(ji,jj) * e2t(ji,jj)
              IF( abs(gphit(ji,jj)) .LE. 15 ) THEN 
                 zahm0(ji,jj) = ahm0
              ELSE 
                 IF( zmsk .GE. zeref ) THEN 
                    zahm0(ji,jj) = ahm0
                 ELSE 
                    zahm0(ji,jj) = zahmm + (ahm0-zahmm)*(1.0 -   &
                        &          cos((rpi*0.5*(zmsk-zemin)/(zeref-zemin))))
                 ENDIF 
              ENDIF 
           END DO 
         END DO 
      ENDIF

      IF( jp_cfg == 2 )   THEN
         zahmeq     = aht0
         zahmm      = ahm0
         zahm0(:,:) = ahm0
      ENDIF

      IF( jp_cfg == 1 )   THEN
         zahmeq     = aht0  ! reduced to aht0 on equator; set to ahm0 if no tropical reduction is required
         zahmm      = ahm0
         zahm0(:,:) = ahm0
      ENDIF

      DO jj = 1, jpj
         DO ji = 1, jpi
            IF( ABS(gphif(ji,jj)) >= 20.) THEN
               ahm2(ji,jj,1) =  zahm0(ji,jj)
            ELSEIF( ABS(gphif(ji,jj)) <= 2.5) THEN
               ahm2(ji,jj,1) =  zahmeq
            ELSE
               ahm2(ji,jj,1) = zahmeq + (zahm0(ji,jj)-zahmeq)/2.   &
                  &            *(1.-COS( rad*(ABS(gphif(ji,jj))-2.5)*180./17.5 ) )
            ENDIF
            IF( ABS(gphit(ji,jj)) >= 20.) THEN
               ahm1(ji,jj,1) =  zahm0(ji,jj)
            ELSEIF( ABS(gphit(ji,jj)) <= 2.5) THEN
               ahm1(ji,jj,1) =  zahmeq
            ELSE
               ahm1(ji,jj,1) = zahmeq + (zahm0(ji,jj)-zahmeq)/2.   &
                  &            *(1.-COS( rad*(ABS(gphit(ji,jj))-2.5)*180./17.5 ) )
            ENDIF
         END DO
      END DO
      
      ! increase along western boundaries of equatorial strip
      ! t-point
      DO jj = 1, jpjm1
         DO ji = 1, jpim1
            zcoft = float( icof(ji,jj) ) / 100.
            ahm1(ji,jj,1) = zcoft * zahm0(ji,jj) + (1.-zcoft) * ahm1(ji,jj,1)
         END DO
      END DO
      ! f-point
      icof(:,:) = icof(:,:) * tmask(:,:,1)
      DO jj = 1, jpjm1
         DO ji = 1, jpim1   ! NO vector opt.
            zmsk = tmask(ji,jj+1,1) + tmask(ji+1,jj+1,1) + tmask(ji,jj,1) + tmask(ji,jj+1,1)
            IF( zmsk == 0. ) THEN
               zcoff = 1.
            ELSE
               zcoff = FLOAT( icof(ji,jj+1) + icof(ji+1,jj+1) + icof(ji,jj) + icof(ji,jj+1) )   &
                     / (zmsk * 100.)
            ENDIF
            ahm2(ji,jj,1) = zcoff * zahm0(ji,jj) + (1.-zcoff) * ahm2(ji,jj,1)
         END DO
      END DO

      ! other level: re-increase the coef in the deep ocean
      !==================================================================
      ! Prior to v3.3, zcoeff was hardwired according to k-index jk.
      !
      ! From v3.3 onwards this has been generalised to a function of 
      ! depth so that it can be used with any number of levels.
      !
      ! The function has been chosen to match the original values (shown
      ! in the following comments) when using the standard 31 ORCA levels.  
      !     DO jk = 1, 21
      !        zcoef(jk) = 1._wp
      !     END DO
      !     zcoef(22) = 2._wp
      !     zcoef(23) = 3._wp
      !     zcoef(24) = 5._wp
      !     zcoef(25) = 7._wp
      !     zcoef(26) = 9._wp
      !     DO jk = 27, jpk
      !        zcoef(jk) = 10._wp
      !     END DO
      !==================================================================

       IF(lwp) THEN
          WRITE(numout,*)
          WRITE(numout,*) '         1D zcoef array '
          WRITE(numout,*) '         ~~~~~~~~~~~~~~ '
          WRITE(numout,*)
          WRITE(numout,*) '    jk        zcoef '
       ENDIF

      DO jk=1, jpk
         zcoef(jk) = 1.0_wp + NINT(9.0_wp*(gdept_0(jk)-800.0_wp)/(3000.0_wp-800.0_wp))
         zcoef(jk) = MIN(10.0_wp, MAX(1.0_wp, zcoef(jk)))
         IF(lwp) WRITE(numout,'(4x,i3,6x,f7.3)') jk,zcoef(jk)
      END DO

      DO jk = 2, jpk
         ahm1(:,:,jk) = MIN( zahm0(:,:), zcoef(jk) * ahm1(:,:,1) )
         ahm2(:,:,jk) = MIN( zahm0(:,:), zcoef(jk) * ahm2(:,:,1) )
      END DO

      IF( jp_cfg == 4 )   THEN         ! Limit AHM in Gibraltar strait
         ij0 = 50   ;   ij1 = 53
         ii0 = 69   ;   ii1 = 71
         DO jk = 1, jpk
            ahm1(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),jk) = MIN( zahmm, ahm1(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),jk) ) 
            ahm2(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),jk) = MIN( zahmm, ahm2(mi0(ii0):mi1(ii1),mj0(ij0):mj1(ij1),jk) )
         END DO
      ENDIF
      CALL lbc_lnk( ahm1, 'T', 1. )   ! Lateral boundary conditions (unchanged sign)
      CALL lbc_lnk( ahm2, 'F', 1. )


      IF(lwp) THEN                    ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '         3D ahm1 array (k=1)'
         CALL prihre( ahm1(:,:,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1.e-3, numout )
         WRITE(numout,*)
         WRITE(numout,*) '         3D ahm2 array (k=1)'
         CALL prihre( ahm2(:,:,1), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1.e-3, numout )
         WRITE(numout,*)
         WRITE(numout,*) '         3D ahm2 array (k=jpk)'
         CALL prihre( ahm2(:,:,jpk), jpi, jpj, 1, jpi, 20, 1, jpj, 20, 1.e-3, numout )
      ENDIF


      ! Set ahm3 and ahm4
      ! =================

      ! define ahm3 and ahm4 at the right grid point position
      ! initialization to a constant value
      !     (USER: modify ahm3 and ahm4 following your desiderata)
      !     harmonic isopycnal or geopotential:
      !                          ahm3 (ahm4) defined at u- (v-) point
      DO jk = 1, jpk
         DO jj = 2, jpj
            DO ji = 2, jpi
               ahm3(ji,jj,jk) = 0.5 * ( ahm2(ji,jj,jk) + ahm2(ji  ,jj-1,jk) )
               ahm4(ji,jj,jk) = 0.5 * ( ahm2(ji,jj,jk) + ahm2(ji-1,jj  ,jk) )
            END DO
         END DO
      END DO
      ahm3 ( :, 1, :) = ahm3 ( :, 2, :)
      ahm4 ( :, 1, :) = ahm4 ( :, 2, :)
      
      CALL lbc_lnk( ahm3, 'U', 1. )    ! Lateral boundary conditions (unchanged sign)
      CALL lbc_lnk( ahm4, 'V', 1. )

      ! Control print

      IF( lwp .AND. ld_print ) THEN
         WRITE(numout,*)
         WRITE(numout,*) '         ahm3 array level 1'
         CALL prihre(ahm3(:,:,1),jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
         WRITE(numout,*)
         WRITE(numout,*) '         ahm4 array level 1'
         CALL prihre(ahm4(:,:,1),jpi,jpj,1,jpi,1,1,jpj,1,1.e-3,numout)
      ENDIF
      !
      CALL wrk_dealloc( jpi   , jpj   , icof  )
      CALL wrk_dealloc( jpidta, jpjdta, idata )
      CALL wrk_dealloc( jpk   ,         zcoef )
      CALL wrk_dealloc( jpi   , jpj   , zahm0 )
      !
   END SUBROUTINE ldf_dyn_c3d_orca


   SUBROUTINE ldf_zpf_1d( ld_print, pdam, pwam, pbot, pdep, pah )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_zpf  ***
      !!                   
      !! ** Purpose :   vertical adimensional profile for eddy coefficient
      !!
      !! ** Method  :   1D eddy viscosity coefficients ( depth )
      !!----------------------------------------------------------------------
      LOGICAL , INTENT(in   )                 ::   ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT(in   )                 ::   pdam       ! depth of the inflection point
      REAL(wp), INTENT(in   )                 ::   pwam       ! width of inflection
      REAL(wp), INTENT(in   )                 ::   pbot       ! bottom value (0<pbot<= 1)
      REAL(wp), INTENT(in   ), DIMENSION(jpk) ::   pdep       ! depth of the gridpoint (T, U, V, F)
      REAL(wp), INTENT(inout), DIMENSION(jpk) ::   pah        ! adimensional vertical profile
      !!
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
      !!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept_0(1    ) ) / pwam )
      zm01 = TANH( ( pdam - gdept_0(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         pah(jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
      END DO

      IF(lwp .AND. ld_print ) THEN      ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(jk), pdep(jk)
         END DO
      ENDIF
      !
   END SUBROUTINE ldf_zpf_1d


   SUBROUTINE ldf_zpf_1d_3d( ld_print, pdam, pwam, pbot, pdep, pah )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_zpf  ***
      !!
      !! ** Purpose :   vertical adimensional profile for eddy coefficient
      !!
      !! ** Method  :   1D eddy viscosity coefficients ( depth )
      !!----------------------------------------------------------------------
      LOGICAL , INTENT(in   )                         ::   ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT(in   )                         ::   pdam       ! depth of the inflection point
      REAL(wp), INTENT(in   )                         ::   pwam       ! width of inflection
      REAL(wp), INTENT(in   )                         ::   pbot       ! bottom value (0<pbot<= 1)
      REAL(wp), INTENT(in   ), DIMENSION          (:) ::   pdep       ! depth of the gridpoint (T, U, V, F)
      REAL(wp), INTENT(inout), DIMENSION      (:,:,:) ::   pah        ! adimensional vertical profile
      !!
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs, zcf  ! temporary scalars
      !!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept_0(1    ) ) / pwam )
      zm01 = TANH( ( pdam - gdept_0(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         zcf = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(jk) ) / pwam )  )
         pah(:,:,jk) = zcf
      END DO

      IF(lwp .AND. ld_print ) THEN      ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(jk)
         END DO
      ENDIF
      !
   END SUBROUTINE ldf_zpf_1d_3d


   SUBROUTINE ldf_zpf_3d( ld_print, pdam, pwam, pbot, pdep, pah )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE ldf_zpf  ***
      !!
      !! ** Purpose :   vertical adimensional profile for eddy coefficient
      !!
      !! ** Method  :   3D for partial step or s-coordinate
      !!----------------------------------------------------------------------
      LOGICAL , INTENT(in   )                         ::   ld_print   ! If true, output arrays on numout
      REAL(wp), INTENT(in   )                         ::   pdam       ! depth of the inflection point
      REAL(wp), INTENT(in   )                         ::   pwam       ! width of inflection
      REAL(wp), INTENT(in   )                         ::   pbot       ! bottom value (0<pbot<= 1)
      REAL(wp), INTENT(in   ), DIMENSION      (:,:,:) ::   pdep       ! dep of the gridpoint (T, U, V, F)
      REAL(wp), INTENT(inout), DIMENSION      (:,:,:) ::   pah        ! adimensional vertical profile
      !!
      INTEGER  ::   jk           ! dummy loop indices
      REAL(wp) ::   zm00, zm01, zmhb, zmhs       ! temporary scalars
      !!----------------------------------------------------------------------

      zm00 = TANH( ( pdam - gdept_0(1    ) ) / pwam )   
      zm01 = TANH( ( pdam - gdept_0(jpkm1) ) / pwam )
      zmhs = zm00 / zm01
      zmhb = ( 1.e0 - pbot ) / ( 1.e0 - zmhs ) / zm01

      DO jk = 1, jpk
         pah(:,:,jk) = 1.e0 + zmhb * ( zm00 - TANH( ( pdam - pdep(:,:,jk) ) / pwam )  )
      END DO

      IF(lwp .AND. ld_print ) THEN      ! Control print
         WRITE(numout,*)
         WRITE(numout,*) '         ahm profile : '
         WRITE(numout,*)
         WRITE(numout,'("  jk      ahm       ","  depth t-level " )')
         DO jk = 1, jpk
            WRITE(numout,'(i6,2f12.4,3x,2f12.4)') jk, pah(1,1,jk), pdep(1,1,jk)
         END DO
      ENDIF
      !
   END SUBROUTINE ldf_zpf_3d

   !!======================================================================
END MODULE ldfdyn
