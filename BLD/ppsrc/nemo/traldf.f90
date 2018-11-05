MODULE traldf
   !!======================================================================
   !!                       ***  MODULE  traldf  ***
   !! Ocean Active tracers : lateral diffusive trends 
   !!=====================================================================
   !! History :  9.0  ! 2005-11 (G. Madec)  Original code
   !!       NEMO 3.0  ! 2008-01  (C. Ethe, G. Madec)  merge TRC-TRA 
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_ldf      : update the tracer trend with the lateral diffusion
   !!   tra_ldf_init : initialization, namelist read, and parameters control
   !!       ldf_ano  : compute lateral diffusion for constant T-S profiles
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE ldftra_oce      ! ocean tracer   lateral physics
   USE ldfslp          ! ???
   USE traldf_bilapg   ! lateral mixing            (tra_ldf_bilapg routine)
   USE traldf_bilap    ! lateral mixing             (tra_ldf_bilap routine)
   USE traldf_iso      ! lateral mixing               (tra_ldf_iso routine)
   USE traldf_iso_grif ! lateral mixing          (tra_ldf_iso_grif routine)
   USE traldf_lap      ! lateral mixing               (tra_ldf_lap routine)
   USE trdmod_oce      ! ocean space and time domain
   USE trdtra          ! ocean active tracers trends
   USE prtctl          ! Print control
   USE in_out_manager  ! I/O manager
   USE lib_mpp         ! distribued memory computing library
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_ldf         ! called by step.F90 
   PUBLIC   tra_ldf_init    ! called by opa.F90 
   !
   INTEGER ::   nldf = 0   ! type of lateral diffusion used defined from ln_traldf_... namlist logicals)

   REAL, SAVE, ALLOCATABLE, DIMENSION(:,:,:) ::   t0_ldf, s0_ldf   !: lateral diffusion trends of T & S for a cst profile
   !                                                               !  (key_traldf_ano only)

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
   !! $Id: traldf.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_ldf( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf  ***
      !! 
      !! ** Purpose :   compute the lateral ocean tracer physics.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      REAL(wp), POINTER, DIMENSION(:,:,:) ::  ztrdt, ztrds
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_ldf')
      !
      rldf = 1     ! For active tracers the 

      IF( l_trdtra )   THEN                    !* Save ta and sa trends
         CALL wrk_alloc( jpi, jpj, jpk, ztrdt, ztrds ) 
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) 
         ztrds(:,:,:) = tsa(:,:,:,jp_sal)
      ENDIF

      SELECT CASE ( nldf )                       ! compute lateral mixing trend and add it to the general trend
      CASE ( 0 )   ;   CALL tra_ldf_lap     ( kt, nit000, 'TRA', gtsu, gtsv, tsb, tsa, jpts        )  ! iso-level laplacian
      CASE ( 1 )                                                                              ! rotated laplacian
         IF( ln_traldf_grif ) THEN                                                          
                       CALL tra_ldf_iso_grif( kt, nit000,'TRA', gtsu, gtsv, tsb, tsa, jpts, ahtb0 )      ! Griffies operator
         ELSE                                                                                
                       CALL tra_ldf_iso     ( kt, nit000, 'TRA', gtsu, gtsv, tsb, tsa, jpts, ahtb0 )      ! Madec operator
         ENDIF
      CASE ( 2 )   ;   CALL tra_ldf_bilap   ( kt, nit000, 'TRA', gtsu, gtsv, tsb, tsa, jpts        )  ! iso-level bilaplacian
      CASE ( 3 )   ;   CALL tra_ldf_bilapg  ( kt, nit000, 'TRA',             tsb, tsa, jpts        )  ! s-coord. geopot. bilap.
         !
      CASE ( -1 )                                ! esopa: test all possibility with control print
         CALL tra_ldf_lap   ( kt, nit000, 'TRA', gtsu, gtsv, tsb, tsa, jpts        ) 
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' ldf0 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         IF( ln_traldf_grif ) THEN
            CALL tra_ldf_iso_grif( kt, nit000, 'TRA', gtsu, gtsv, tsb, tsa, jpts, ahtb0 )
         ELSE
            CALL tra_ldf_iso     ( kt, nit000, 'TRA', gtsu, gtsv, tsb, tsa, jpts, ahtb0 )  
         ENDIF
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' ldf1 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_ldf_bilap ( kt, nit000, 'TRA', gtsu, gtsv, tsb, tsa, jpts        ) 
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' ldf2 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
         CALL tra_ldf_bilapg( kt, nit000, 'TRA',             tsb, tsa, jpts        ) 
         CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' ldf3 - Ta: ', mask1=tmask,               &
         &             tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      END SELECT


      IF( l_trdtra )   THEN                      ! save the horizontal diffusive trends for further diagnostics
         ztrdt(:,:,:) = tsa(:,:,:,jp_tem) - ztrdt(:,:,:)
         ztrds(:,:,:) = tsa(:,:,:,jp_sal) - ztrds(:,:,:)
         CALL trd_tra( kt, 'TRA', jp_tem, jptra_trd_ldf, ztrdt )
         CALL trd_tra( kt, 'TRA', jp_sal, jptra_trd_ldf, ztrds )
         CALL wrk_dealloc( jpi, jpj, jpk, ztrdt, ztrds ) 
      ENDIF
      !                                          ! print mean trends (used for debugging)
      IF(ln_ctl)   CALL prt_ctl( tab3d_1=tsa(:,:,:,jp_tem), clinfo1=' ldf  - Ta: ', mask1=tmask,               &
         &                       tab3d_2=tsa(:,:,:,jp_sal), clinfo2=       ' Sa: ', mask2=tmask, clinfo3='tra' )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_ldf')
      !
   END SUBROUTINE tra_ldf


   SUBROUTINE tra_ldf_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_ldf_init  ***
      !! 
      !! ** Purpose :   Choice of the operator for the lateral tracer diffusion
      !!
      !! ** Method  :   set nldf from the namtra_ldf logicals
      !!      nldf == -1   ESOPA test: ALL operators are used
      !!      nldf ==  0   laplacian operator
      !!      nldf ==  1   Rotated laplacian operator
      !!      nldf ==  2   bilaplacian operator
      !!      nldf ==  3   Rotated bilaplacian
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio, ierr         ! temporary integers 
      !!----------------------------------------------------------------------

      !  Define the lateral mixing oparator for tracers
      ! ===============================================
    
      IF(lwp) THEN                    ! Namelist print
         WRITE(numout,*)
         WRITE(numout,*) 'tra_ldf_init : lateral tracer diffusive operator'
         WRITE(numout,*) '~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namtra_ldf already read in ldftra module'
         WRITE(numout,*) '   see ldf_tra_init report for lateral mixing parameters'
         WRITE(numout,*)
      ENDIF

      !                               ! control the input
      ioptio = 0
      IF( ln_traldf_lap   )   ioptio = ioptio + 1
      IF( ln_traldf_bilap )   ioptio = ioptio + 1
      IF( ioptio >  1 )   CALL ctl_stop( '          use ONE or NONE of the 2 lap/bilap operator type on tracer' )
      IF( ioptio == 0 )   nldf = -2   ! No lateral diffusion
      ioptio = 0
      IF( ln_traldf_level )   ioptio = ioptio + 1
      IF( ln_traldf_hor   )   ioptio = ioptio + 1
      IF( ln_traldf_iso   )   ioptio = ioptio + 1
      IF( ioptio >  1 )   CALL ctl_stop( '          use only ONE direction (level/hor/iso)' )

      ! defined the type of lateral diffusion from ln_traldf_... logicals
      ! CAUTION : nldf = 1 is used in trazdf_imp, change it carefully
      ierr = 0
      IF( ln_traldf_lap ) THEN       ! laplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_traldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_traldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_traldf_level )   ierr = 1      ! iso-level not allowed
            IF ( ln_traldf_hor   )   nldf = 0      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! z-coordinate
            IF ( ln_traldf_level )   nldf = 0      ! iso-level  (no rotation)
            IF ( ln_traldf_hor   )   nldf = 1      ! horizontal (   rotation)
            IF ( ln_traldf_iso   )   nldf = 1      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ln_traldf_bilap ) THEN      ! bilaplacian operator
         IF ( ln_zco ) THEN                ! z-coordinate
            IF ( ln_traldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_traldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_zps ) THEN             ! z-coordinate
            IF ( ln_traldf_level )   ierr = 1      ! iso-level not allowed 
            IF ( ln_traldf_hor   )   nldf = 2      ! horizontal (no rotation)
            IF ( ln_traldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
         IF ( ln_sco ) THEN             ! z-coordinate
            IF ( ln_traldf_level )   nldf = 2      ! iso-level  (no rotation)
            IF ( ln_traldf_hor   )   nldf = 3      ! horizontal (   rotation)
            IF ( ln_traldf_iso   )   ierr = 2      ! isoneutral (   rotation)
         ENDIF
      ENDIF

      IF( ierr == 1 )   CALL ctl_stop( ' iso-level in z-coordinate - partial step, not allowed' )
      IF( ierr == 2 )   CALL ctl_stop( ' isoneutral bilaplacian operator does not exist' )
      IF( lk_traldf_eiv .AND. .NOT.ln_traldf_iso )   &
           CALL ctl_stop( '          eddy induced velocity on tracers',   &
           &              ' the eddy induced velocity on tracers requires isopycnal laplacian diffusion' )
      IF( nldf == 1 .OR. nldf == 3 ) THEN      ! rotation
         IF( .NOT.lk_ldfslp )   CALL ctl_stop( '          the rotation of the diffusive tensor require key_ldfslp' )
         l_traldf_rot = .TRUE.                 ! needed for trazdf_imp
      ENDIF

      IF( lk_esopa ) THEN
         IF(lwp) WRITE(numout,*) '          esopa control: use all lateral physics options'
         nldf = -1
      ENDIF

      IF(lwp) THEN
         WRITE(numout,*)
         IF( nldf == -2 )   WRITE(numout,*) '          NO lateral diffusion'
         IF( nldf == -1 )   WRITE(numout,*) '          ESOPA test All scheme used'
         IF( nldf ==  0 )   WRITE(numout,*) '          laplacian operator'
         IF( nldf ==  1 )   WRITE(numout,*) '          Rotated laplacian operator'
         IF( nldf ==  2 )   WRITE(numout,*) '          bilaplacian operator'
         IF( nldf ==  3 )   WRITE(numout,*) '          Rotated bilaplacian'
      ENDIF

      ! Reference T & S diffusivity (if necessary)
      ! ===========================
      CALL ldf_ano
      !
   END SUBROUTINE tra_ldf_init

   !!----------------------------------------------------------------------
   !!   default option :   Dummy code   NO T & S background profiles
   !!----------------------------------------------------------------------
   SUBROUTINE ldf_ano
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'tra:ldf_ano : lateral diffusion acting on the full fields'
         WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
   END SUBROUTINE ldf_ano

   !!======================================================================
END MODULE traldf
