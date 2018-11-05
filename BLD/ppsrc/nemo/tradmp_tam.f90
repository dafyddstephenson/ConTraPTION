MODULE tradmp_tam
   !!======================================================================
   !!                       ***  MODULE  tradmp_tam  ***
   !! Ocean physics: internal restoring trend on active tracers (T and S)
   !!                Tangent and Adjoint Module
   !!======================================================================
   !! History of the direct module:
   !!            OPA  ! 1991-03  (O. Marti, G. Madec)  Original code
   !!                 ! 1992-06  (M. Imbard)  doctor norme
   !!                 ! 1996-01  (G. Madec)  statement function for e3
   !!                 ! 1997-05  (G. Madec)  macro-tasked on jk-slab
   !!                 ! 1998-07  (M. Imbard, G. Madec) ORCA version
   !!            7.0  ! 2001-02  (M. Imbard)  cofdis, Original code
   !!            8.1  ! 2001-02  (G. Madec, E. Durand)  cleaning
   !!  NEMO      1.0  ! 2002-08  (G. Madec, E. Durand)  free form + modules
   !!                 ! 2003-08  (M. Balmaseda)  mods to allow eq. damping
   !!            3.2  ! 2009-08  (G. Madec, C. Talandier)  DOCTOR norm for namelist parameter
   !! History of the TAM:
   !!                 ! 2008-09  (A. Vidard) tangent and adjoint module
   !!                                       of the 03-08 version
   !!        NEMO 3.2 ! 2010-04 (F. Vigilant) 3.2 version
   !!        NEMO 3.4 ! 2012-07 (P.-A. Bouttier) 3.4 Version@
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   tra_dmp      : update the tracer trend with the internal damping
   !!   tra_dmp_init : initialization, namlist read, parameters control
   !!   dtacof_zoom  : restoring coefficient for zoom domain
   !!   dtacof       : restoring coefficient for global domain
   !!   cofdis       : compute the distance to the coastline
   !!----------------------------------------------------------------------
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE zdf_oce
   USE in_out_manager
   USE phycst
   USE lib_mpp
   USE tradmp
   USE prtctl
   USE gridrandom
   USE dotprodfld
   USE dtatsd
   USE zdfmxl
   USE tstool_tam
   USE paresp
   USE lib_mpp
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC tra_dmp_tan      ! routine called by step_tam.F90
   PUBLIC tra_dmp_adj      ! routine called by step_tam.F90
   PUBLIC tra_dmp_init_tam
   PUBLIC tra_dmp_adj_tst  ! routine called by tst.F90
   !LOGICAL, PUBLIC            ::   lk_tradmp = .TRUE.     !: internal damping flag
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) :: &
     & strdmp_tl,   & !: damping salinity trend (psu/s)
     & ttrdmp_tl,   & !: damping temperature trend (psu/s)
     & ttrdmp_ad,   & !: damping temperature trend (psu/s)
     & strdmp_ad      !: damping salinity trend (psu/s)


   LOGICAL :: lfirst = .TRUE.        !: flag for initialisation
   !                                 !!* Namelist namtra_dmp : T & S newtonian damping *
   INTEGER  ::   nn_hdmp =   -1      ! = 0/-1/'latitude' for damping over T and S
   INTEGER  ::   nn_zdmp =    0      ! = 0/1/2 flag for damping in the mixed layer
   REAL(wp) ::   rn_surf =   50._wp  ! surface time scale for internal damping        [days]
   REAL(wp) ::   rn_bot  =  360._wp  ! bottom time scale for internal damping         [days]
   REAL(wp) ::   rn_dep  =  800._wp  ! depth of transition between rn_surf and rn_bot [meters]
   INTEGER  ::   nn_file =    2      ! = 1 create a damping.coeff NetCDF file
   LOGICAL, PRIVATE, SAVE :: ll_alloctl = .FALSE.
   LOGICAL, PRIVATE, SAVE :: ll_allocad = .FALSE.
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

   INTEGER FUNCTION tra_dmp_alloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_dmp_alloc_tam  ***
      !!----------------------------------------------------------------------
      INTEGER, OPTIONAL :: kmode
      INTEGER :: ierr(2)
      integer :: jmode

      IF ( PRESENT(kmode) ) THEN 
         jmode = kmode 
      ELSE
         jmode = 0
      END IF

      IF ( ( jmode == 0 ) .OR. ( jmode == 1 ) .AND. ( .NOT. ll_alloctl) ) THEN
         ALLOCATE( strdmp_tl(jpi,jpj,jpk) , ttrdmp_tl(jpi,jpj,jpk), STAT= ierr(1) )
      END IF
      IF ( ( jmode == 0 ) .OR. ( jmode == 2 ) .AND. ( .NOT. ll_allocad) ) THEN
         ALLOCATE( strdmp_ad(jpi,jpj,jpk) , ttrdmp_ad(jpi,jpj,jpk), STAT= ierr(2) )
      END IF
      tra_dmp_alloc_tam = SUM( ierr )
      !
      IF( lk_mpp            )   CALL mpp_sum ( tra_dmp_alloc_tam )
      IF( tra_dmp_alloc_tam > 0 )   CALL ctl_warn('tra_dmp_alloc_tam: allocation of arrays failed')
      !
   END FUNCTION tra_dmp_alloc_tam

   INTEGER FUNCTION tra_dmp_dealloc_tam( kmode )
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION tra_dmp_alloc_tam  ***
      !!----------------------------------------------------------------------
      INTEGER, OPTIONAL :: kmode
      INTEGER :: ierr(2)
      integer :: jmode

      IF ( PRESENT( kmode) ) THEN 
         jmode = kmode 
      ELSE
         jmode = 0
      END IF

      IF ( ( jmode == 0 ) .OR. ( jmode == 1 ) .AND. ( ll_alloctl) ) THEN
         DEALLOCATE( strdmp_tl, ttrdmp_tl, STAT= ierr(1) )
      END IF
      IF ( ( jmode == 0 ) .OR. ( jmode == 2 ) .AND. ( ll_allocad) ) THEN
         DEALLOCATE( strdmp_ad, ttrdmp_ad, STAT= ierr(2) )
      END IF
      tra_dmp_dealloc_tam = SUM( ierr )
      !
      IF( lk_mpp            )   CALL mpp_sum ( tra_dmp_dealloc_tam )
      IF( tra_dmp_dealloc_tam > 0 )   CALL ctl_warn('tra_dmp_dealloc_tam: deallocation of arrays failed')
      !
   END FUNCTION tra_dmp_dealloc_tam

   SUBROUTINE tra_dmp_tan( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp_tan  ***
      !!
      !! ** Purpose of direct routine:
      !!      Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - (ta,sa)   tracer trends updated with the damping trend
      !!
      !! ASSUME key_zdfcst_tam
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jn            ! dummy loop indices
      REAL(wp) ::   ztest, ztatl, zsatl       ! temporary scalars
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_dmp_tan')
      !
      ! 1. Newtonian damping trends on tracer fields
      ! --------------------------------------------
      !    compute the newtonian damping trends depending on nmldmp

      SELECT CASE ( nn_zdmp )
      !
      CASE( 0 )                ! newtonian damping throughout the water column
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ztatl = - resto(ji,jj,jk) *  tsb_tl(ji,jj,jk,jp_tem)
                  zsatl = - resto(ji,jj,jk) *  tsb_tl(ji,jj,jk,jp_sal)
                  ! add the trends to the general tracer trends
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem) + ztatl
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal) + zsatl
                  ! save the salinity trend (used in flx to close the salt budget)
                  strdmp_tl(ji,jj,jk) = zsatl
                  ttrdmp_tl(ji,jj,jk) = ztatl
               END DO
            END DO
         END DO
         !
      CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ztest = avt(ji,jj,jk) - 5.e-4_wp
      !! ASSUME key_zdfcst_tam
                  IF( ztest < 0. ) THEN
                     ztatl = - resto(ji,jj,jk) * tsb_tl(ji,jj,jk,jp_tem)
                     zsatl = - resto(ji,jj,jk) * tsb_tl(ji,jj,jk,jp_sal)
                  ELSE
                     ztatl = 0.0_wp
                     zsatl = 0.0_wp
                  ENDIF
                  ! add the trends to the general tracer trends
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem) + ztatl
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal) + zsatl
                  ! save the salinity trend (used in flx to close the salt budget)
                  strdmp_tl(ji,jj,jk) = zsatl
                  ttrdmp_tl(ji,jj,jk) = ztatl
               END DO
            END DO
         END DO
         !
      CASE ( 2 )                ! no damping in the mixed layer
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  IF( gdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     ztatl = - resto(ji,jj,jk) * tsb_tl(ji,jj,jk,jp_tem)
                     zsatl = - resto(ji,jj,jk) * tsb_tl(ji,jj,jk,jp_sal)
                  ELSE
                     ztatl = 0.e0
                     zsatl = 0.e0
                  ENDIF
                  ! add the trends to the general tracer trends
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem) + ztatl
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal) + zsatl
                  ! save the salinity trend (used in flx to close the salt budget)
                  strdmp_tl(ji,jj,jk) = zsatl
                  ttrdmp_tl(ji,jj,jk) = ztatl
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_dmp_tan')
      !
   END SUBROUTINE tra_dmp_tan
   SUBROUTINE tra_dmp_adj( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_dmp_adj  ***
      !!
      !! ** Purpose of direct routine:
      !!      Compute the tracer trend due to a newtonian damping
      !!      of the tracer field towards given data field and add it to the
      !!      general tracer trends.
      !!
      !! ** Method  :   Newtonian damping towards t_dta and s_dta computed
      !!      and add to the general tracer trends:
      !!                     ta = ta + resto * (t_dta - tb)
      !!                     sa = sa + resto * (s_dta - sb)
      !!         The trend is computed either throughout the water column
      !!      (nlmdmp=0) or in area of weak vertical mixing (nlmdmp=1) or
      !!      below the well mixed layer (nlmdmp=2)
      !!
      !! ** Action  : - (ta,sa)   tracer trends updated with the damping trend
      !!              - save the trends in (ttrd,strd) ('key_trdtra')
      !! ASSUME key_zdfcst_tam
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jn            ! dummy loop indices
      REAL(wp) ::   ztest, ztaad, zsaad       ! temporary scalars
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_dmp_adj')
      !
      ! 1. Newtonian damping trends on tracer fields
      ! --------------------------------------------
      !    compute the newtonian damping trends depending on nmldmp
      ztaad = 0.0_wp
      zsaad = 0.0_wp

      SELECT CASE ( nn_zdmp )
      !
      CASE( 0 )                ! newtonian damping throughout the water column
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ! save the salinity trend (used in flx to close the salt budget)
                  ztaad = ztaad + ttrdmp_ad(ji,jj,jk)
                  ttrdmp_ad(ji,jj,jk) = 0.0_wp
                  zsaad = zsaad + strdmp_ad(ji,jj,jk)
                  strdmp_ad(ji,jj,jk) = 0.0_wp
                  ! add the trends to the general tracer trends
                  ztaad = tsa_ad(ji,jj,jk,jp_tem) + ztaad
                  zsaad = tsa_ad(ji,jj,jk,jp_sal) + zsaad

                  tsb_ad(ji,jj,jk,jp_tem)  = tsb_ad(ji,jj,jk,jp_tem) - ztaad * resto(ji,jj,jk)
                  tsb_ad(ji,jj,jk,jp_sal)  = tsb_ad(ji,jj,jk,jp_sal) - zsaad * resto(ji,jj,jk)
                  ztaad = 0.0_wp
                  zsaad = 0.0_wp
                END DO
            END DO
         END DO
         !
      CASE ( 1 )                ! no damping in the turbocline (avt > 5 cm2/s)
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ! save the salinity trend (used in flx to close the salt budget)
                  ztaad = ztaad + ttrdmp_ad(ji,jj,jk)
                  ttrdmp_ad(ji,jj,jk) = 0.0_wp
                  zsaad = zsaad + strdmp_ad(ji,jj,jk)
                  strdmp_ad(ji,jj,jk) = 0.0_wp
                  ! add the trends to the general tracer trends
                  ztaad = tsa_ad(ji,jj,jk,jp_tem) + ztaad
                  zsaad = tsa_ad(ji,jj,jk,jp_sal) + zsaad
                  ztest = avt(ji,jj,jk) - 5.e-4
      !! ASSUME key_zdfcst_tam
                  IF( ztest < 0. ) THEN
                     tsb_ad(ji,jj,jk,jp_tem)  = tsb_ad(ji,jj,jk,jp_tem) - ztaad * resto(ji,jj,jk)
                     tsb_ad(ji,jj,jk,jp_sal)  = tsb_ad(ji,jj,jk,jp_sal) - zsaad * resto(ji,jj,jk)
                     ztaad = 0.0_wp
                     zsaad = 0.0_wp
                  ELSE
                     ztaad = 0.0_wp
                     zsaad = 0.0_wp
                  ENDIF
               END DO
            END DO
         END DO
         !
      CASE ( 2 )                ! no damping in the mixed layer
         DO jk = 1, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  ! save the salinity trend (used in flx to close the salt budget)
                  ztaad = ztaad + ttrdmp_ad(ji,jj,jk)
                  ttrdmp_ad(ji,jj,jk) = 0.0_wp
                  zsaad = zsaad + strdmp_ad(ji,jj,jk)
                  strdmp_ad(ji,jj,jk) = 0.0_wp
                  ! add the trends to the general tracer trends
                  ztaad = tsa_ad(ji,jj,jk,jp_tem) + ztaad
                  zsaad = tsa_ad(ji,jj,jk,jp_sal) + zsaad
                  IF( gdept(ji,jj,jk) >= hmlp (ji,jj) ) THEN
                     tsb_ad(ji,jj,jk,jp_tem)  = tsb_ad(ji,jj,jk,jp_tem) - ztaad * resto(ji,jj,jk)
                     tsb_ad(ji,jj,jk,jp_sal)  = tsb_ad(ji,jj,jk,jp_sal) - zsaad * resto(ji,jj,jk)
                     ztaad = 0.0_wp
                     zsaad = 0.0_wp
                  ELSE
                     ztaad = 0.0_wp
                     zsaad = 0.0_wp
                  ENDIF
               END DO
            END DO
         END DO
         !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_dmp_adj')
      !
   END SUBROUTINE tra_dmp_adj
   SUBROUTINE tra_dmp_init_tam
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_dmp_init_tam  ***
      !!
      !! ** Purpose :   Initialization for the newtonian damping
      !!
      !! ** Method  :   read the nammbf namelist and check the parameters
      !!      called by tra_dmp at the first timestep (nit000)
      !!----------------------------------------------------------------------

      !NAMELIST/namtra_dmp/ nn_hdmp, nn_file, nn_zdmp, rn_surf, rn_bot, rn_dep
      !!----------------------------------------------------------------------

      IF (lfirst) THEN

         !REWIND ( numnam )                  ! Read Namelist namtdp : temperature and salinity damping term
         !READ   ( numnam, namtra_dmp )
         !IF( lzoom )   nn_zdmp = 0           ! restoring to climatology at closed north or south boundaries

         !IF(lwp) THEN                       ! Namelist print
            !WRITE(numout,*)
            !WRITE(numout,*) 'tra_dmp_init_tam : T and S newtonian damping'
            !WRITE(numout,*) '~~~~~~~~~~~'
            !WRITE(numout,*) '      Namelist namtra_dmp : set damping parameter'
            !WRITE(numout,*) '      T and S damping option         nn_hdmp = ', nn_hdmp
            !WRITE(numout,*) '      mixed layer damping option     nn_zdmp = ', nn_zdmp, '(zoom: forced to 0)'
            !WRITE(numout,*) '      surface time scale (days)      rn_surf = ', rn_surf
            !WRITE(numout,*) '      bottom time scale (days)       rn_bot  = ', rn_bot
            !WRITE(numout,*) '      depth of transition (meters)   rn_dep  = ', rn_dep
            !WRITE(numout,*) '      create a damping.coeff file    nn_file = ', nn_file
         !ENDIF

         !! ** modified to allow damping at the equator
         IF( ln_tradmp ) THEN               ! initialization for T-S damping
         !
            IF( tra_dmp_alloc_tam() /= 0 )   CALL ctl_stop( 'STOP', 'tra_dmp_init_tam: unable to allocate arrays' )
            SELECT CASE ( nn_hdmp )
            CASE (  -1  )   ;   IF(lwp) WRITE(numout,*) '   tracer damping in the Med & Red seas only'
            CASE ( 0:90 )   ;   IF(lwp) WRITE(numout,*) '   tracer damping poleward of', nn_hdmp, ' degrees'
            CASE DEFAULT
               IF(lwp) WRITE(numout,*) '          tracer damping disabled', nn_hdmp
               !         WRITE(ctmp1,*) '          bad flag value for ndmp = ', ndmp
               !         CALL ctl_stop(ctmp1)
            END SELECT

            SELECT CASE ( nn_zdmp )
            CASE ( 0 )   ;   IF(lwp) WRITE(numout,*) '   tracer damping throughout the water column'
            CASE ( 1 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the turbocline (avt > 5 cm2/s)'
            CASE ( 2 )   ;   IF(lwp) WRITE(numout,*) '   no tracer damping in the mixed layer'
            CASE DEFAULT
               WRITE(ctmp1,*) '   bad flag value for nn_zdmp = ', nn_zdmp
               CALL ctl_stop(ctmp1)
            END SELECT
            !
            !* We don't need data for TAM *!
            !IF( .NOT. ln_tsd_tradmp )
               !CALL ctl_warn( 'tra_dmp_init_tam: read T-S data not initialized, we force ln_tsd_tradmp=T' )
            !ENDIF
            !                          ! Damping coefficients initialization
            !IF( lzoom ) THEN   ;   CALL dtacof_zoom( resto )
            !ELSE               ;   CALL dtacof( nn_hdmp, rn_surf, rn_bot, rn_dep, nn_file, 'TRA', resto )
            !ENDIF
            !
            lfirst = .FALSE.
         END IF
         ttrdmp_tl(:,:,:) = 0.0_wp       ! internal damping salinity trend (used in ocesbc)
         ttrdmp_ad(:,:,:) = 0.0_wp       ! internal damping salinity trend (used in ocesbc)
         strdmp_tl(:,:,:) = 0.0_wp       ! internal damping salinity trend (used in ocesbc)
         strdmp_ad(:,:,:) = 0.0_wp       ! internal damping salinity trend (used in ocesbc)
      END IF
   END SUBROUTINE tra_dmp_init_tam

   SUBROUTINE tra_dmp_adj_tst ( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!          ***  ROUTINE tra_sbc_adj_tst : TEST OF tra_sbc_adj  ***
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
       !! History :
      !!        ! 08-08 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER ::  &
         & istp,  &
         & jstp,  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zsb_tlin,     &! Tangent input : before salinity
         & ztb_tlin,     &! Tangent input : before temperature
         & zsa_tlin,     &! Tangent input : after salinity
         & zta_tlin,     &! Tangent input : after temperature
         & zsa_tlout,    &! Tangent output: after salinity
         & zta_tlout,    &! Tangent output: after temperature
         & zstrdmp_tlout,&! Tangent output: salinity trend
         & zttrdmp_tlout,&! Tangent output: temperature trend
         & zsb_adout,    &! Adjoint output : before salinity
         & ztb_adout,    &! Adjoint output : before temperature
         & zsa_adout,    &! Adjoint output : after salinity
         & zta_adout,    &! Adjoint output : after temperature
         & zsa_adin,     &! Adjoint input : after salinity
         & zta_adin,     &! Adjoint input : after temperature
         & zstrdmp_adin, &! Adjoint input : salinity trend
         & zttrdmp_adin, &! Tangent output: temperature trend
         & z3r            ! 3D field

      REAL(KIND=wp) ::       &
         & zsp1,             & ! scalar product involving the tangent routine
         & zsp1_1,           & ! scalar product involving the tangent routine
         & zsp1_2,           & ! scalar product involving the tangent routine
         & zsp1_3,           & ! scalar product involving the tangent routine
         & zsp1_4,           & ! scalar product involving the tangent routine
         & zsp2,             & ! scalar product involving the adjoint routine
         & zsp2_1,           & ! scalar product involving the adjoint routine
         & zsp2_2,           & ! scalar product involving the adjoint routine
         & zsp2_3,           & ! scalar product involving the adjoint routine
         & zsp2_4              ! scalar product involving the adjoint routine
      CHARACTER(LEN=14) ::   &
         & cl_name

      ALLOCATE( &
         & ztb_tlin(jpi,jpj,jpk),     &
         & zsb_tlin(jpi,jpj,jpk),     &
         & zsa_tlin(jpi,jpj,jpk),     &
         & zta_tlin(jpi,jpj,jpk),     &
         & zsa_tlout(jpi,jpj,jpk),    &
         & zta_tlout(jpi,jpj,jpk),    &
         & zstrdmp_tlout(jpi,jpj,jpk),&
         & zttrdmp_tlout(jpi,jpj,jpk),&
         & ztb_adout(jpi,jpj,jpk),    &
         & zsb_adout(jpi,jpj,jpk),    &
         & zsa_adout(jpi,jpj,jpk),    &
         & zta_adout(jpi,jpj,jpk),    &
         & zsa_adin(jpi,jpj,jpk),     &
         & zta_adin(jpi,jpj,jpk),     &
         & zttrdmp_adin(jpi,jpj,jpk), &
         & zstrdmp_adin(jpi,jpj,jpk), &
         & z3r         (jpi,jpj,jpk)  &
         & )
      ! Test for time steps nit000 and nit000 + 1 (the matrix changes)

      DO jstp = nit000, nit000 + 2
         istp = jstp
         IF ( jstp == nit000 +2 ) istp = nitend

      ! Initialize the reference state
      avt(:,:,:) = 1.e-1
      !=============================================================
      ! 1) dx = ( T ) and dy = ( T )
      !=============================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      zsb_tlin(:,:,:)      = 0.0_wp
      ztb_tlin(:,:,:)      = 0.0_wp
      zsa_tlin(:,:,:)      = 0.0_wp
      zta_tlin(:,:,:)      = 0.0_wp
      zsa_tlout(:,:,:)     = 0.0_wp
      zta_tlout(:,:,:)     = 0.0_wp
      zstrdmp_tlout(:,:,:) = 0.0_wp
      zttrdmp_tlout(:,:,:) = 0.0_wp
      zsb_adout(:,:,:)     = 0.0_wp
      ztb_adout(:,:,:)     = 0.0_wp
      zsa_adout(:,:,:)     = 0.0_wp
      zta_adout(:,:,:)     = 0.0_wp
      zsa_adin(:,:,:)      = 0.0_wp
      zta_adin(:,:,:)      = 0.0_wp
      zstrdmp_adin(:,:,:)  = 0.0_wp
      zttrdmp_adin(:,:,:)  = 0.0_wp

      tsa_ad(:,:,:,:)      = 0.0_wp
      tsb_ad(:,:,:,:)      = 0.0_wp
      strdmp_ad(:,:,:)  = 0.0_wp
      ttrdmp_ad(:,:,:)  = 0.0_wp
      tsa_tl(:,:,:,:)      = 0.0_wp
      tsb_tl(:,:,:,:)      = 0.0_wp
      strdmp_tl(:,:,:)  = 0.0_wp
      ttrdmp_tl(:,:,:)  = 0.0_wp

      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zsb_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               ztb_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  z3r, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zsa_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      CALL grid_random(  z3r, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zta_tlin(ji,jj,jk) = z3r(ji,jj,jk)
            END DO
         END DO
      END DO

      tsb_tl(:,:,:,jp_sal) = zsb_tlin(:,:,:)
      tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
      tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)
      tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)

      CALL tra_dmp_tan( istp )

      zsa_tlout(:,:,:)     = tsa_tl(:,:,:,jp_sal)
      zta_tlout(:,:,:)     = tsa_tl(:,:,:,jp_tem)
      zstrdmp_tlout(:,:,:) = strdmp_tl(:,:,:)
      zttrdmp_tlout(:,:,:) = ttrdmp_tl(:,:,:)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsa_adin(ji,jj,jk)     = zsa_tlout(ji,jj,jk) &
                 &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &                   * tmask(ji,jj,jk) * wesp_s(jk)
              zta_adin(ji,jj,jk)     = zta_tlout(ji,jj,jk) &
                 &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &                   * tmask(ji,jj,jk) * wesp_t(jk)
              zstrdmp_adin(ji,jj,jk) = zstrdmp_tlout(ji,jj,jk) &
                 &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &                   * tmask(ji,jj,jk) * wesp_s(jk)
              zttrdmp_adin(ji,jj,jk) = zttrdmp_tlout(ji,jj,jk) &
                 &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &                   * tmask(ji,jj,jk) * wesp_t(jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1_1 = DOT_PRODUCT( zsa_tlout    , zsa_adin     )
      zsp1_2 = DOT_PRODUCT( zta_tlout    , zta_adin     )
      zsp1_3 = DOT_PRODUCT( zstrdmp_tlout, zstrdmp_adin )
      zsp1_4 = DOT_PRODUCT( zttrdmp_tlout, zttrdmp_adin )
      zsp1   = zsp1_1 + zsp1_2 + zsp1_3 + zsp1_4

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      tsa_ad(:,:,:,jp_sal)     = zsa_adin(:,:,:)
      tsa_ad(:,:,:,jp_tem)     = zta_adin(:,:,:)
      strdmp_ad(:,:,:) = zstrdmp_adin(:,:,:)
      ttrdmp_ad(:,:,:) = zttrdmp_adin(:,:,:)

      CALL tra_dmp_adj( istp )

      zsb_adout(:,:,:) = tsb_ad(:,:,:,jp_sal)
      ztb_adout(:,:,:) = tsb_ad(:,:,:,jp_tem)
      zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)
      zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1 = DOT_PRODUCT( zsb_tlin  , zsb_adout   )
      zsp2_2 = DOT_PRODUCT( ztb_tlin  , ztb_adout   )
      zsp2_3 = DOT_PRODUCT( zsa_tlin  , zsa_adout   )
      zsp2_4 = DOT_PRODUCT( zta_tlin  , zta_adout   )

      zsp2 = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4

      ! Compare the scalar products

      ! 14 char:'12345678901234'
      IF ( istp == nit000 ) THEN
         cl_name = 'tra_dmp_adj T1'
      ELSEIF ( istp == nit000 +1 ) THEN
         cl_name = 'tra_dmp_adj T2'
      ELSEIF ( istp == nitend ) THEN
         cl_name = 'tra_dmp_adj T3'
      END IF
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      END DO

      DEALLOCATE( &
         & zsb_tlin,     &
         & ztb_tlin,     &
         & zsa_tlin,     &
         & zta_tlin,     &
         & zsa_tlout,    &
         & zta_tlout,    &
         & zstrdmp_tlout,&
         & zttrdmp_tlout,&
         & zsb_adout,    &
         & ztb_adout,    &
         & zsa_adout,    &
         & zta_adout,    &
         & zsa_adin,     &
         & zta_adin,     &
         & zstrdmp_adin, &
         & zttrdmp_adin, &
         & z3r           &
         & )

   END SUBROUTINE tra_dmp_adj_tst

END MODULE tradmp_tam
