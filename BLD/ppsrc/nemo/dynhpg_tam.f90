MODULE dynhpg_tam
   !!======================================================================
   !!                       ***  MODULE  dynhpg_tam  ***
   !! Ocean dynamics:  hydrostatic pressure gradient trend
   !!                  Tangent and Adjoint module
   !!======================================================================
   !! History of the direct module:
   !!            1.0  !  87-09  (P. Andrich, M.-A. Foujols)  hpg_zco: Original code
   !!            5.0  !  91-11  (G. Madec)
   !!            7.0  !  96-01  (G. Madec)  hpg_sco: Original code for s-coordinates
   !!            8.0  !  97-05  (G. Madec)  split dynber into dynkeg and dynhpg
   !!            8.5  !  02-07  (G. Madec)  F90: Free form and module
   !!            8.5  !  02-08  (A. Bozec)  hpg_zps: Original code
   !!            9.0  !  05-10  (A. Beckmann, B.W. An)  various s-coordinate options
   !!                           Original code for hpg_ctl, hpg_hel hpg_wdj, hpg_djc, hpg_rot
   !!            9.0  !  05-11  (G. Madec) style & small optimisation
   !! History of the TAM module:
   !!            9.0  !  08-06  (A. Vidard) Skeleton
   !!                 !  08-11  (A. Vidard) Nemo v3
   !!                 !  12-07  (P.-A. Bouttier) 3.4 version
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_hpg      : update the momentum trend with the now horizontal
   !!                  gradient of the hydrostatic pressure
   !!       hpg_ctl  : initialisation and control of options
   !!       hpg_zco  : z-coordinate scheme
   !!       hpg_zps  : z-coordinate plus partial steps (interpolation)
   !!       hpg_sco  : s-coordinate (standard jacobian formulation)
   !!       hpg_hel  : s-coordinate (helsinki modification)
   !!       hpg_wdj  : s-coordinate (weighted density jacobian)
   !!       hpg_djc  : s-coordinate (Density Jacobian with Cubic polynomial)
   !!       hpg_rot  : s-coordinate (ROTated axes scheme)
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE dynhpg
   USE phycst
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE lib_mpp
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   dyn_hpg_tan    ! routine called by step_tam module
   PUBLIC   dyn_hpg_adj    ! routine called by step_tam module
   PUBLIC   dyn_hpg_init_tam
   PUBLIC   dyn_hpg_adj_tst! routine called by test module

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

   SUBROUTINE dyn_hpg_tan( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg_tan  ***
      !!
      !! ** Method of the direct routine:
      !!              Call the hydrostatic pressure gradient routine
      !!              using the scheme defined in the namelist
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_hpg_tan')
      !
      SELECT CASE ( nhpg )      ! Hydrastatic pressure gradient computation
      CASE (  0 )   ;   CALL hpg_zco_tan    ( kt )      ! z-coordinate
      CASE (  1 )   ;   CALL hpg_zps_tan    ( kt )      ! z-coordinate plus partial steps (interpolation)
      CASE (  2 )   ;   CALL hpg_sco_tan    ( kt )      ! s-coordinate (standard jacobian formulation)
      CASE (  3 )   ;   CALL hpg_djc_tan    ( kt )      ! s-coordinate (helsinki modification)
      CASE (  4 )   ;   CALL hpg_prj_tan    ( kt )      ! s-coordinate (weighted density jacobian)
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_hpg_tan')
      !
   END SUBROUTINE dyn_hpg_tan
   SUBROUTINE dyn_hpg_adj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE dyn_hpg_adj  ***
      !!
      !! ** Method of the direct routine:
      !!              call the hydrostatic pressure gradient routine
      !!              using the scheme defined in the namelist
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend (l_trddyn=T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_hpg_adj')
      !
      SELECT CASE ( nhpg )      ! Hydrastatic pressure gradient computation
      CASE (  0 )   ;   CALL hpg_zco_adj    ( kt )      ! z-coordinate
      CASE (  1 )   ;   CALL hpg_zps_adj    ( kt )      ! z-coordinate plus partial steps (interpolation)
      CASE (  2 )   ;   CALL hpg_sco_adj    ( kt )      ! s-coordinate (standard jacobian formulation)
      CASE (  3 )   ;   CALL hpg_djc_adj    ( kt )      ! s-coordinate (helsinki modification)
      CASE (  4 )   ;   CALL hpg_prj_adj    ( kt )      ! s-coordinate (weighted density jacobian)
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('dyn_hpg_adj')
      !
   END SUBROUTINE dyn_hpg_adj

   SUBROUTINE dyn_hpg_init_tam
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE hpg_ctl_tam  ***
      !!
      !! ** Purpose :   initializations for the hydrostatic pressure gradient
      !!              computation and consistency control
      !!
      !! ** Action  :   Read the namelist namdyn_hpg and check the consistency
      !!      with the type of vertical coordinate used (zco, zps, sco)
      !!----------------------------------------------------------------------
      INTEGER ::   ioptio = 0      ! temporary integer

      NAMELIST/namdyn_hpg/ ln_hpg_zco, ln_hpg_zps, ln_hpg_sco,     &
         &                 ln_hpg_djc, ln_hpg_prj, ln_dynhpg_imp
      !!----------------------------------------------------------------------

      REWIND ( numnam )               ! Read Namelist nam_dynhpg : pressure gradient calculation options
      READ   ( numnam, namdyn_hpg )

      IF(lwp) THEN                    ! Control print
         WRITE(numout,*)
         WRITE(numout,*) 'dyn_ctl_tam : hydrostatic pressure gradient'
         WRITE(numout,*) '~~~~~~~~~~~~~~~'
         WRITE(numout,*) '   Namelist namdyn_hpg : choice of hpg scheme'
         WRITE(numout,*) '      z-coord. - full steps                             ln_hpg_zco    = ', ln_hpg_zco
         WRITE(numout,*) '      z-coord. - partial steps (interpolation)          ln_hpg_zps    = ', ln_hpg_zps
         WRITE(numout,*) '      s-coord. (standard jacobian formulation)          ln_hpg_sco    = ', ln_hpg_sco
         WRITE(numout,*) '      s-coord. (Density Jacobian: Cubic polynomial)     ln_hpg_djc    = ', ln_hpg_djc
         WRITE(numout,*) '      s-coord. (Pressure Jacobian: Cubic polynomial)    ln_hpg_prj    = ', ln_hpg_prj
         WRITE(numout,*) '      time stepping: centered (F) or semi-implicit (T)  ln_dynhpg_imp = ', ln_dynhpg_imp
      ENDIF
      !
      IF( ln_hpg_djc )   &
         &   CALL ctl_stop('dyn_hpg_init_tam : Density Jacobian: Cubic polynominal method &
                           & currently disabled (bugs under investigation). Please select &
                           & either  ln_hpg_sco or  ln_hpg_prj instead')
      IF( lk_vvl .AND. .NOT. ln_hpg_sco )   THEN
         CALL ctl_stop( 'dyn_hpg_init_tam : variable volume key_vvl compatible only with the standard jacobian formulation hpg_sco')
      ENDIF
      !                               ! Set nhpg from ln_hpg_... flags
      IF( ln_hpg_zco )   nhpg = 0
      IF( ln_hpg_zps )   nhpg = 1
      IF( ln_hpg_sco )   nhpg = 2
      IF( ln_hpg_djc )   nhpg = 3
      IF( ln_hpg_prj )   nhpg = 4
      !                               ! Consitency check
      ioptio = 0
      IF( ln_hpg_zco )   ioptio = ioptio + 1
      IF( ln_hpg_zps )   ioptio = ioptio + 1
      IF( ln_hpg_sco )   ioptio = ioptio + 1
      IF( ln_hpg_djc )   ioptio = ioptio + 1
      IF( ln_hpg_prj )   ioptio = ioptio + 1
      IF ( ioptio /= 1 )   CALL ctl_stop( ' NO or several hydrostatic pressure gradient options used' )
      !
   END SUBROUTINE dyn_hpg_init_tam
   SUBROUTINE hpg_zco_tan( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_zco_tan  ***
      !!
      !! ** Method of the direct routine:
      !!      z-coordinate case, levels are horizontal surfaces.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      add it to the general momentum trend (ua,va).
      !!            ua = ua - 1/e1u * zhpi
      !!            va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua_tl,va_tl) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk       ! dummy loop indices
      REAL(wp) ::   zcoef0, zcoef1   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zhpitl, zhpjtl
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj,jpk, zhpitl, zhpjtl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zco_tan : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   z-coordinate case '
      ENDIF
      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp
      ! Surface value
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            zcoef1 = zcoef0 * e3w(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpitl(ji,jj,1) = zcoef1 * ( rhd_tl(ji+1,jj  ,1) - rhd_tl(ji,jj,1) ) / e1u(ji,jj)
            zhpjtl(ji,jj,1) = zcoef1 * ( rhd_tl(ji  ,jj+1,1) - rhd_tl(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua_tl(ji,jj,1) = ua_tl(ji,jj,1) + zhpitl(ji,jj,1)
            va_tl(ji,jj,1) = va_tl(ji,jj,1) + zhpjtl(ji,jj,1)
         END DO
      END DO
      !
      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zcoef1 = zcoef0 * e3w(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpitl(ji,jj,jk) = zhpitl(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd_tl(ji+1,jj,jk)+rhd_tl(ji+1,jj,jk-1) )   &
                  &                       - ( rhd_tl(ji  ,jj,jk)+rhd_tl(ji  ,jj,jk-1) )  ) / e1u(ji,jj)

               zhpjtl(ji,jj,jk) = zhpjtl(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd_tl(ji,jj+1,jk)+rhd_tl(ji,jj+1,jk-1) )   &
                  &                       - ( rhd_tl(ji,jj,  jk)+rhd_tl(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) + zhpitl(ji,jj,jk)
               va_tl(ji,jj,jk) = va_tl(ji,jj,jk) + zhpjtl(ji,jj,jk)
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zhpitl, zhpjtl )
      !
   END SUBROUTINE hpg_zco_tan
   SUBROUTINE hpg_zco_adj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_zco_tan  ***
      !!
      !! ** Method of the direct routine:
      !!      z-coordinate case, levels are horizontal surfaces.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level:    zhpi = grav .....
      !!                zhpj = grav .....
      !!      add it to the general momentum trend (ua,va).
      !!            ua = ua - 1/e1u * zhpi
      !!            va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua_tl,va_tl) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk       ! dummy loop indices
      REAL(wp) ::   zcoef0, zcoef1   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zhpiad, zhpjad
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj,jpk, zhpiad, zhpjad )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zco_adj : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   z-coordinate case '
      ENDIF
      ! adjoint variables initialization
      zhpiad = 0.0_wp
      zhpjad = 0.0_wp
      ! Local constant initialization
      zcoef0 = - grav * 0.5

      ! interior value (2=<jk=<jpkm1)
      DO jk = jpkm1, 2, -1
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               zcoef1 = zcoef0 * e3w(ji,jj,jk)
               ! add to the general momentum trend
               zhpiad(ji,jj,jk) = zhpiad(ji,jj,jk) + ua_ad(ji,jj,jk)
               zhpjad(ji,jj,jk) = zhpjad(ji,jj,jk) + va_ad(ji,jj,jk)
               ! hydrostatic pressure gradient
               rhd_ad(ji,jj+1,jk  ) = rhd_ad(ji,jj+1,jk  ) + zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               rhd_ad(ji,jj+1,jk-1) = rhd_ad(ji,jj+1,jk-1) + zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               rhd_ad(ji,jj  ,jk  ) = rhd_ad(ji,jj  ,jk  ) - zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               rhd_ad(ji,jj  ,jk-1) = rhd_ad(ji,jj  ,jk-1) - zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               zhpjad(ji,jj  ,jk-1) = zhpjad(ji,jj  ,jk-1) + zhpjad(ji,jj,jk)
               zhpjad(ji,jj  ,jk  ) = 0.0_wp
               !
               rhd_ad(ji+1,jj,jk  ) = rhd_ad(ji+1,jj,jk  ) + zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               rhd_ad(ji+1,jj,jk-1) = rhd_ad(ji+1,jj,jk-1) + zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               rhd_ad(ji  ,jj,jk  ) = rhd_ad(ji  ,jj,jk  ) - zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               rhd_ad(ji  ,jj,jk-1) = rhd_ad(ji  ,jj,jk-1) - zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               zhpiad(ji  ,jj,jk-1) = zhpiad(ji  ,jj,jk-1) + zhpiad(ji,jj,jk)
               zhpiad(ji  ,jj,jk  ) = 0.0_wp
               !
            END DO
         END DO
      END DO
      ! Surface value
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            zcoef1 = zcoef0 * e3w(ji,jj,1)
            ! add to the general momentum trend
            zhpiad(ji,jj,1) = zhpiad(ji,jj,1) + ua_ad(ji,jj,1)
            zhpjad(ji,jj,1) = zhpjad(ji,jj,1) + va_ad(ji,jj,1)
            ! hydrostatic pressure gradient
            rhd_ad(ji,jj+1,1) = rhd_ad(ji,jj+1,1) + zhpjad(ji,jj,1) * zcoef1 / e2v(ji,jj)
            rhd_ad(ji,jj  ,1) = rhd_ad(ji,jj  ,1) - zhpjad(ji,jj,1) * zcoef1 / e2v(ji,jj)
            zhpjad(ji,jj,1) = 0.0_wp
            !
            rhd_ad(ji+1,jj,1) = rhd_ad(ji+1,jj,1) + zhpiad(ji,jj,1) * zcoef1 / e1u(ji,jj)
            rhd_ad(ji  ,jj,1) = rhd_ad(ji  ,jj,1) - zhpiad(ji,jj,1) * zcoef1 / e1u(ji,jj)
            zhpiad(ji,jj,1) = 0.0_wp
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zhpiad, zhpjad )
      !
   END SUBROUTINE hpg_zco_adj
   SUBROUTINE hpg_zps_tan( kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE hpg_zps  ***
      !!
      !! ** Method of the direct routine:
      !!              z-coordinate plus partial steps case.  blahblah...
      !!
      !! ** Action  : - Update (ua_tl,va_tl) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk                       ! dummy loop indices
      INTEGER  ::   iku, ikv                         ! temporary integers
      REAL(wp) ::   zcoef0, zcoef1, zcoef2, zcoef3   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:) :: zhpitl, zhpjtl
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj,jpk, zhpitl, zhpjtl )
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zps_tan : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   z-coordinate with partial steps - vector optimization'
      ENDIF

      ! Local constant initialization
      zcoef0 = - grav * 0.5_wp

      !  Surface value
      DO jj = 2, jpjm1
         DO ji = 2, jpim1   ! vector opt.
            zcoef1 = zcoef0 * e3w(ji,jj,1)
            ! hydrostatic pressure gradient
            zhpitl(ji,jj,1) = zcoef1 * ( rhd_tl(ji+1,jj  ,1) - rhd_tl(ji,jj,1) ) / e1u(ji,jj)
            zhpjtl(ji,jj,1) = zcoef1 * ( rhd_tl(ji  ,jj+1,1) - rhd_tl(ji,jj,1) ) / e2v(ji,jj)
            ! add to the general momentum trend
            ua_tl(ji,jj,1) = ua_tl(ji,jj,1) + zhpitl(ji,jj,1)
            va_tl(ji,jj,1) = va_tl(ji,jj,1) + zhpjtl(ji,jj,1)
         END DO
      END DO

      ! interior value (2=<jk=<jpkm1)
      DO jk = 2, jpkm1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               zcoef1 = zcoef0 * e3w(ji,jj,jk)
               ! hydrostatic pressure gradient
               zhpitl(ji,jj,jk) = zhpitl(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd_tl(ji+1,jj,jk) + rhd_tl(ji+1,jj,jk-1) )   &
                  &                       - ( rhd_tl(ji  ,jj,jk) + rhd_tl(ji  ,jj,jk-1) )  ) / e1u(ji,jj)

               zhpjtl(ji,jj,jk) = zhpjtl(ji,jj,jk-1)   &
                  &           + zcoef1 * (  ( rhd_tl(ji,jj+1,jk) + rhd_tl(ji,jj+1,jk-1) )   &
                  &                       - ( rhd_tl(ji,jj,  jk) + rhd_tl(ji,jj  ,jk-1) )  ) / e2v(ji,jj)
               ! add to the general momentum trend
               ua_tl(ji,jj,jk) = ua_tl(ji,jj,jk) + zhpitl(ji,jj,jk)
               va_tl(ji,jj,jk) = va_tl(ji,jj,jk) + zhpjtl(ji,jj,jk)
            END DO
         END DO
      END DO

      ! partial steps correction at the last level  (new gradient with  intgrd.F)
      DO jj = 2, jpjm1
         DO ji = 2, jpim1
            iku = mbku(ji,jj)
            ikv = mbkv(ji,jj)
            zcoef2 = zcoef0 * MIN( e3w(ji,jj,iku), e3w(ji+1,jj  ,iku) )
            zcoef3 = zcoef0 * MIN( e3w(ji,jj,ikv), e3w(ji  ,jj+1,ikv) )
            ! on i-direction
            IF ( iku > 1 ) THEN             ! on i-direction (level 2 or more)
               ua_tl (ji,jj,iku) = ua_tl(ji,jj,iku) - zhpitl(ji,jj,iku)   ! subtract old value
               zhpitl(ji,jj,iku) = zhpitl(ji,jj,iku-1)               &    ! compute the new one
                  &              + zcoef2 * ( rhd_tl(ji+1,jj,iku-1) - rhd_tl(ji,jj,iku-1) + gru_tl(ji,jj) ) / e1u(ji,jj)
               ua_tl (ji,jj,iku) = ua_tl(ji,jj,iku) + zhpitl(ji,jj,iku)   ! add the new one to the general momentum trend
            ENDIF

            IF ( ikv > 1 ) THEN            ! on j-direction
               va_tl(ji,jj,ikv) = va_tl(ji,jj,ikv) - zhpjtl(ji,jj,ikv)    ! subtract old value
               zhpjtl (ji,jj,ikv) = zhpjtl(ji,jj,ikv-1)   &               ! compute the new one
                  &               + zcoef3 * ( rhd_tl(ji,jj+1,ikv-1) - rhd_tl(ji,jj,ikv-1) + grv_tl(ji,jj) ) / e2v(ji,jj)
               va_tl(ji,jj,ikv) = va_tl(ji,jj,ikv) + zhpjtl(ji,jj,ikv)    ! add the new one to the general momentum trend
            ENDIF
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zhpitl, zhpjtl )
      !
   END SUBROUTINE hpg_zps_tan
   SUBROUTINE hpg_zps_adj( kt )
      !!---------------------------------------------------------------------
      !!                 ***  ROUTINE hpg_zps  ***
      !!
      !! ** Method of the direct routine:
      !!              z-coordinate plus partial steps case.  blahblah...
      !!
      !! ** Action  : - Update (ua_tl,va_tl) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk                       ! dummy loop indices
      INTEGER  ::   iku, ikv                         ! temporary integers
      REAL(wp) ::   zcoef0, zcoef1, zcoef2, zcoef3   ! temporary scalars
      REAL(wp), POINTER, DIMENSION(:,:,:):: zhpiad, zhpjad
      !!----------------------------------------------------------------------
      !
      CALL wrk_alloc( jpi,jpj,jpk, zhpiad, zhpjad )
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn:hpg_zps_adj : hydrostatic pressure gradient trend'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~~~~~   z-coordinate with partial steps - vector optimization'
      ENDIF
      zhpiad(:,:,:) = 0.0_wp
      zhpjad(:,:,:) = 0.0_wp
      ! Local constant initialization
      zcoef0 = - grav * 0.5

      ! partial steps correction at the last level  (new gradient with  intgrd.F)
      DO jj = jpjm1, 2, -1
         DO ji = jpim1, 2, -1
            iku = mbku(ji,jj)
            ikv = mbkv(ji,jj)
            zcoef2 = zcoef0 * MIN( e3w(ji,jj,iku), e3w(ji+1,jj  ,iku) )
            zcoef3 = zcoef0 * MIN( e3w(ji,jj,ikv), e3w(ji  ,jj+1,ikv) )
            ! on i-direction
            IF ( iku > 2 ) THEN
               ! add the new one to the general momentum trend
               zhpiad(ji,jj,iku) = zhpiad(ji,jj,iku) + ua_ad(ji,jj,iku)
               ! compute the new one
               rhd_ad(ji+1,jj,iku-1) = rhd_ad(ji+1,jj,iku-1) + zhpiad (ji,jj,iku) * zcoef2 / e1u(ji,jj)
               rhd_ad(ji,jj,iku-1) = rhd_ad(ji,jj,iku-1) - zhpiad (ji,jj,iku) * zcoef2 / e1u(ji,jj)
               gru_ad(ji,jj) = gru_ad(ji,jj) + zhpiad (ji,jj,iku) * zcoef2 / e1u(ji,jj)
               zhpiad(ji,jj,iku-1) = zhpiad(ji,jj,iku-1) + zhpiad (ji,jj,iku)
               zhpiad (ji,jj,iku) = 0.0_wp
               ! subtract old value
               zhpiad(ji,jj,iku) = zhpiad(ji,jj,iku) - ua_ad(ji,jj,iku)
            ENDIF
            ! on j-direction
            IF ( ikv > 2 ) THEN
               ! add the new one to the general momentum trend
               zhpjad(ji,jj,ikv) = zhpjad(ji,jj,ikv) + va_ad(ji,jj,ikv)
               ! compute the new one
               rhd_ad(ji,jj+1,ikv-1) = rhd_ad(ji,jj+1,ikv-1) + zhpjad (ji,jj,ikv) * zcoef3 / e2v(ji,jj)
               rhd_ad(ji,jj,ikv-1) = rhd_ad(ji,jj,ikv-1) -zhpjad (ji,jj,ikv) * zcoef3 / e2v(ji,jj)
               grv_ad(ji,jj) = grv_ad(ji,jj) +zhpjad (ji,jj,ikv) * zcoef3 / e2v(ji,jj)
               zhpjad(ji,jj,ikv-1) = zhpjad(ji,jj,ikv-1) + zhpjad(ji,jj,ikv)
               zhpjad (ji,jj,ikv) = 0.0_wp
               ! subtract old value
               zhpjad(ji,jj,ikv) = zhpjad(ji,jj,ikv) - va_ad(ji,jj,ikv)
            ENDIF
         END DO
      END DO
      !
      ! interior value (2=<jk=<jpkm1)
      DO jk = jpkm1, 2, -1
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               zcoef1 = zcoef0 * e3w(ji,jj,jk)
               ! add to the general momentum trend
               zhpiad(ji,jj,jk) = zhpiad(ji,jj,jk) + ua_ad(ji,jj,jk)
               zhpjad(ji,jj,jk) = zhpjad(ji,jj,jk) + va_ad(ji,jj,jk)
               ! hydrostatic pressure gradient
               rhd_ad(ji,jj+1,jk  ) = rhd_ad(ji,jj+1,jk  ) + zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               rhd_ad(ji,jj+1,jk-1) = rhd_ad(ji,jj+1,jk-1) + zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               rhd_ad(ji,jj  ,jk  ) = rhd_ad(ji,jj  ,jk  ) - zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               rhd_ad(ji,jj  ,jk-1) = rhd_ad(ji,jj  ,jk-1) - zhpjad(ji,jj,jk) * zcoef1 / e2v(ji,jj)
               zhpjad(ji,jj  ,jk-1) = zhpjad(ji,jj  ,jk-1) + zhpjad(ji,jj,jk)
               zhpjad(ji,jj  ,jk  ) = 0.0_wp
               !
               rhd_ad(ji+1,jj,jk  ) = rhd_ad(ji+1,jj,jk  ) + zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               rhd_ad(ji+1,jj,jk-1) = rhd_ad(ji+1,jj,jk-1) + zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               rhd_ad(ji  ,jj,jk  ) = rhd_ad(ji  ,jj,jk  ) - zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               rhd_ad(ji  ,jj,jk-1) = rhd_ad(ji  ,jj,jk-1) - zhpiad(ji,jj,jk) * zcoef1 / e1u(ji,jj)
               zhpiad(ji  ,jj,jk-1) = zhpiad(ji  ,jj,jk-1) + zhpiad(ji,jj,jk)
               zhpiad(ji  ,jj,jk  ) = 0.0_wp
            END DO
         END DO
      END DO
      !  Surface value
      DO jj = jpjm1, 2, -1
         DO ji = jpim1, 2, -1   ! vector opt.
            zcoef1 = zcoef0 * e3w(ji,jj,1)
            ! add to the general momentum trend
            zhpiad(ji,jj,1) = zhpiad(ji,jj,1) + ua_ad(ji,jj,1)
            zhpjad(ji,jj,1) = zhpjad(ji,jj,1) + va_ad(ji,jj,1)
            ! hydrostatic pressure gradient
            rhd_ad(ji+1,jj  ,1) = rhd_ad(ji+1,jj  ,1) + zhpiad(ji,jj,1) * zcoef1 / e1u(ji,jj)
            rhd_ad(ji  ,jj  ,1) = rhd_ad(ji  ,jj  ,1) - zhpiad(ji,jj,1) * zcoef1 / e1u(ji,jj)
            rhd_ad(ji  ,jj+1,1) = rhd_ad(ji  ,jj+1,1) + zhpjad(ji,jj,1) * zcoef1 / e2v(ji,jj)
            rhd_ad(ji  ,jj  ,1) = rhd_ad(ji  ,jj  ,1) - zhpjad(ji,jj,1) * zcoef1 / e2v(ji,jj)
            zhpiad(ji  ,jj  ,1) = 0.0_wp
            zhpjad(ji  ,jj  ,1) = 0.0_wp
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi,jpj,jpk, zhpiad, zhpjad )
      !
   END SUBROUTINE hpg_zps_adj
   SUBROUTINE hpg_sco_tan( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_sco_tan  ***
      !!
      !! ** Method of the direct routine:   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      CALL ctl_stop( 'hpg_sco_tan not available yet')
   END SUBROUTINE hpg_sco_tan
   SUBROUTINE hpg_sco_adj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_sco_adj  ***
      !!
      !! ** Method of the direct routine:   s-coordinate case. Jacobian scheme.
      !!      The now hydrostatic pressure gradient at a given level, jk,
      !!      is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      CALL ctl_stop( 'hpg_sco_adj not available yet')
   END SUBROUTINE hpg_sco_adj
   SUBROUTINE hpg_djc_tan( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_hel_tan  ***
      !!
      !! ** Method of the direct routine:   s-coordinate case.
      !!      The now hydrostatic pressure gradient at a given level
      !!      jk is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend (l_trddyn=T)
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      CALL ctl_stop( 'hpg_djc_tan not available yet')
   END SUBROUTINE hpg_djc_tan
   SUBROUTINE hpg_djc_adj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_hel_adj  ***
      !!
      !! ** Method of the direct routine:   s-coordinate case.
      !!      The now hydrostatic pressure gradient at a given level
      !!      jk is computed by taking the vertical integral of the in-situ
      !!      density gradient along the model level from the suface to that
      !!      level. s-coordinates (ln_sco): a corrective term is added
      !!      to the horizontal pressure gradient :
      !!         zhpi = grav .....  + 1/e1u mi(rhd) di[ grav dep3w ]
      !!         zhpj = grav .....  + 1/e2v mj(rhd) dj[ grav dep3w ]
      !!      add it to the general momentum trend (ua,va).
      !!         ua = ua - 1/e1u * zhpi
      !!         va = va - 1/e2v * zhpj
      !!
      !! ** Action : - Update (ua,va) with the now hydrastatic pressure trend
      !!             - Save the trend (l_trddyn=T)
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      CALL ctl_stop( 'hpg_djc_adj not available yet')
   END SUBROUTINE hpg_djc_adj
   SUBROUTINE hpg_prj_tan( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_wdj_tan  ***
      !!
      !! ** Method of the direct roiutine:
      !!      Weighted Density Jacobian (wdj) scheme (song 1998)
      !!      The weighting coefficients from the namelist parameter gamm
      !!      (alpha=0.5-gamm ; beta=1-alpha=0.5+gamm)
      !!
      !! Reference : Song, Mon. Wea. Rev., 126, 3213-3230, 1998.
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      CALL ctl_stop( 'hpg_prj_tan not available yet')
   END SUBROUTINE hpg_prj_tan
   SUBROUTINE hpg_prj_adj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE hpg_wdj_adj  ***
      !!
      !! ** Method of the direct roiutine:
      !!      Weighted Density Jacobian (wdj) scheme (song 1998)
      !!      The weighting coefficients from the namelist parameter gamm
      !!      (alpha=0.5-gamm ; beta=1-alpha=0.5+gamm)
      !!
      !! Reference : Song, Mon. Wea. Rev., 126, 3213-3230, 1998.
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      CALL ctl_stop( 'hpg_prj_adj not available yet')
   END SUBROUTINE hpg_prj_adj

   SUBROUTINE dyn_hpg_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE dynhpg_adj_tst ***
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
      !!        ! 08-07 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   &
         & zrhd_tlin,              &  ! in situ density anomalie
         & zua_tlin,               &  ! after u- velocity
         & zva_tlin,               &  ! after v- velocity
         & zua_tlout,              &  ! after u- velocity
         & zva_tlout                  ! after v- velocity
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE ::   &
         & zrhd_adout,             &  ! in situ density anomalie
         & zua_adout,              &  ! after u- velocity
         & zva_adout,              &  ! after v- velocity
         & zua_adin,               &  ! after u- velocity
         & zva_adin                   ! after v- velocity
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   &
         & zgru_tlin,              &
         & zgrv_tlin,              &
         & zgru_adout,             &
         & zgrv_adout

      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zrh,          & ! 3D random field for rhd
         & zau,          & ! 3D random field for u
         & zav             ! 3D random field for v
      REAL(wp), DIMENSION(:,:), ALLOCATABLE ::   &
         & zgru,         & ! 2D random field for gru
         & zgrv            ! 2D random field for grv
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
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d           ! 2D seed for the random number generator
      INTEGER :: &
         & iseed, &
         & ji, &
         & jj, &
         & jk
      CHARACTER(LEN=14) :: cl_name

      ! Allocate memory
      ALLOCATE( &
         & zrhd_tlin(jpi,jpj,jpk),  &
         & zua_tlin(jpi,jpj,jpk),   &
         & zva_tlin(jpi,jpj,jpk),   &
         & zgru_tlin(jpi,jpj),      &
         & zgrv_tlin(jpi,jpj),      &
         & zua_tlout(jpi,jpj,jpk),  &
         & zva_tlout(jpi,jpj,jpk),  &
         & zrhd_adout(jpi,jpj,jpk), &
         & zua_adout(jpi,jpj,jpk),  &
         & zva_adout(jpi,jpj,jpk),  &
         & zgru_adout(jpi,jpj),     &
         & zgrv_adout(jpi,jpj),     &
         & zua_adin(jpi,jpj,jpk),   &
         & zva_adin(jpi,jpj,jpk),   &
         & zrh(jpi,jpj,jpk),        &
         & zau(jpi,jpj,jpk),        &
         & zav(jpi,jpj,jpk),        &
         & zgru(jpi,jpj),           &
         & zgrv(jpi,jpj)            &
         &     )


      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      zrhd_tlin(:,:,:)  = 0.0_wp
      zua_tlin(:,:,:)   = 0.0_wp
      zva_tlin(:,:,:)   = 0.0_wp
      zgru_tlin(:,:)    = 0.0_wp
      zgrv_tlin(:,:)    = 0.0_wp
      zua_tlout(:,:,:)  = 0.0_wp
      zva_tlout(:,:,:)  = 0.0_wp
      zgru_adout(:,:)   = 0.0_wp
      zgrv_adout(:,:)   = 0.0_wp
      zrhd_adout(:,:,:) = 0.0_wp
      zua_adout(:,:,:)  = 0.0_wp
      zva_adout(:,:,:)  = 0.0_wp
      zua_adin(:,:,:)   = 0.0_wp
      zva_adin(:,:,:)   = 0.0_wp
      zrh(:,:,:)        = 0.0_wp
      zau(:,:,:)        = 0.0_wp
      zav(:,:,:)        = 0.0_wp
      zgru(:,:)         = 0.0_wp
      zgrv(:,:)         = 0.0_wp


      gru_tl(:,:)   = 0.0_wp
      grv_tl(:,:)   = 0.0_wp
      gru_ad(:,:)   = 0.0_wp
      grv_ad(:,:)   = 0.0_wp
      ua_tl(:,:,:)  = 0.0_wp
      va_tl(:,:,:)  = 0.0_wp
      rhd_tl(:,:,:) = 0.0_wp
      ua_ad(:,:,:)  = 0.0_wp
      va_ad(:,:,:)  = 0.0_wp
      rhd_ad(:,:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  zau, 'U', 0.0_wp, stdu )
      CALL grid_random(  zav, 'V', 0.0_wp, stdv )
      CALL grid_random(  zrh, 'W', 0.0_wp, stdr )
      CALL grid_random(  zgru, 'U', 0.0_wp, stdu )
      CALL grid_random(  zgrv, 'V', 0.0_wp, stdv )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zrhd_tlin(ji,jj,jk) = zrh(ji,jj,jk)
               zua_tlin(ji,jj,jk)  = zau(ji,jj,jk)
               zva_tlin(ji,jj,jk)  = zav(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgru_tlin(ji,jj)   = zgru(ji,jj)
            zgrv_tlin(ji,jj)   = zgrv(ji,jj)
         END DO
      END DO
      ua_tl(:,:,:)  = zua_tlin(:,:,:)
      va_tl(:,:,:)  = zva_tlin(:,:,:)
      rhd_tl(:,:,:) = zrhd_tlin(:,:,:)
      gru_tl(:,:)   = zgru_tlin(:,:)
      grv_tl(:,:)   = zgrv_tlin(:,:)

      CALL dyn_hpg_tan ( nit000 )

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

      CALL dyn_hpg_adj ( nit000 )

      zgru_adout(:,:)   = gru_ad(:,:)
      zgrv_adout(:,:)   = grv_ad(:,:)
      zrhd_adout(:,:,:) = rhd_ad(:,:,:)
      zua_adout(:,:,:)  = ua_ad(:,:,:)
      zva_adout(:,:,:)  = va_ad(:,:,:)

      zsp2_1 = DOT_PRODUCT( zgru_tlin, zgru_adout )
      zsp2_2 = DOT_PRODUCT( zgrv_tlin, zgrv_adout )
      zsp2_3 = DOT_PRODUCT( zrhd_tlin, zrhd_adout )
      zsp2_4 = DOT_PRODUCT( zua_tlin, zua_adout )
      zsp2_5 = DOT_PRODUCT( zva_tlin, zva_adout )
      zsp2   = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5
      ! Compare the scalar products

      cl_name = 'dyn_hpg_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE( &
         & zrhd_tlin,  &
         & zua_tlin,   &
         & zva_tlin,   &
         & zgru_tlin,  &
         & zgrv_tlin,  &
         & zua_tlout,  &
         & zva_tlout,  &
         & zrhd_adout, &
         & zua_adout,  &
         & zva_adout,  &
         & zgru_adout, &
         & zgrv_adout, &
         & zua_adin,   &
         & zva_adin,   &
         & zrh,        &
         & zau,        &
         & zav,        &
         & zgru,       &
         & zgrv        &
         &              )
   END SUBROUTINE dyn_hpg_adj_tst
END MODULE dynhpg_tam
