MODULE closea_tam
   !!======================================================================
   !!                       ***  MODULE  closea_tam  ***
   !! Closed Seas  : specific treatments associated with closed seas
   !!                Tangent and adjoint module
   !!======================================================================
   !! History of the direct module:
   !!             8.2  !  00-05  (O. Marti)  Original code
   !!             8.5  !  02-06  (E. Durand, G. Madec)  F90
   !!             9.0  !  06-07  (G. Madec)  add clo_rnf, clo_ups, clo_bat
   !! History of the tangent module:
   !!             9.0  !  08-11  (A. Vidard)   skeleton tam of sbc_clo
   !!             9.0  !  09-08  (F. Vigilant) tangent and adjoint version
   !!        NEMO 3.4  !  07-12  (P.-A. Bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_clo    : Special handling of closed seas
   !!----------------------------------------------------------------------

   USE par_oce
   USE par_kind
   USE sbc_oce_tam
   USE lbclnk_tam
   USE lbclnk
   USE dom_oce
   USE sbc_oce
   USE closea
   USE lib_mpp                 ! distribued memory computing library
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE tstool_tam

   IMPLICIT NONE
   PRIVATE

   PUBLIC sbc_clo_tan      ! routine called by sbcmod_tam module
   PUBLIC sbc_clo_adj      ! routine called by sbcmod_tam module
   PUBLIC sbc_clo_adj_tst  ! routine called by tst module

   REAL(wp), DIMENSION (jpncs+1)       ::   surf             ! closed sea surface

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
   !!  OPA 9.0 , LOCEAN-IPSL (2006)
   !! $Id: closea.F90 1146 2008-06-25 11:42:56Z rblod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS
   SUBROUTINE sbc_clo_tan( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_clo_tan  ***
      !!
      !! ** Purpose :   Special handling of closed seas (Tangent version)
      !!
      !! ** Method  :   Water flux is forced to zero over closed sea
      !!      Excess is shared between remaining ocean, or
      !!      put as run-off in open ocean.
      !!
      !! ** Action  :   emp, emps   updated surface freshwater fluxes at kt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean model time step
      !
      INTEGER                     ::   ji, jj, jc, jn   ! dummy loop indices
      REAL(wp)                    ::   zze2
      REAL(wp), DIMENSION (jpncs) ::   zfwftl
      !!----------------------------------------------------------------------
      !
      !                                                   !------------------!
      IF( kt == nit000 ) THEN                             !  Initialisation  !
         !                                                !------------------!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'sbc_clo_tan : closed seas '
         IF(lwp) WRITE(numout,*)'~~~~~~~'

         ! Total surface of ocean
         surf(jpncs+1) = SUM( e1t(:,:) * e2t(:,:) * tmask_i(:,:) )

         DO jc = 1, jpncs
            surf(jc) =0.0_wp
            DO jj = ncsj1(jc), ncsj2(jc)
               DO ji = ncsi1(jc), ncsi2(jc)
                  surf(jc) = surf(jc) + e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)      ! surface of closed seas
               END DO
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum ( surf, jpncs+1 )       ! mpp: sum over all the global domain

         IF(lwp) WRITE(numout,*)'     Closed sea surfaces'
         DO jc = 1, jpncs
            IF(lwp)WRITE(numout,FMT='(1I3,4I4,5X,F16.2)') jc, ncsi1(jc), ncsi2(jc), ncsj1(jc), ncsj2(jc), surf(jc)
         END DO

         ! jpncs+1 : surface of sea, closed seas excluded
         DO jc = 1, jpncs
            surf(jpncs+1) = surf(jpncs+1) - surf(jc)
         END DO
         !
      ENDIF
      !                                                   !--------------------!
      !                                                   !  update emp, emps  !
      zfwftl = 0.0_wp                                     !--------------------!
      DO jc = 1, jpncs
         DO jj = ncsj1(jc), ncsj2(jc)
            DO ji = ncsi1(jc), ncsi2(jc)
               zfwftl(jc) = zfwftl(jc) + e1t(ji,jj) * e2t(ji,jj) * ( emp_tl(ji,jj) - rnf(ji,jj) ) * tmask_i(ji,jj)
            END DO
         END DO
      END DO
      IF( lk_mpp )   CALL mpp_sum ( zfwftl(:) , jpncs )       ! mpp: sum over all the global domain

      IF( cp_cfg == "orca" .AND. jp_cfg == 2 ) THEN      ! Black Sea case for ORCA_R2 configuration
         zze2    = ( zfwftl(3) + zfwftl(4) ) / 2.
         zfwftl(3) = zze2
         zfwftl(4) = zze2
      ENDIF
      DO jc = 1, jpncs
         !
         IF( ncstt(jc) == 0 ) THEN
            ! water/evap excess is shared by all open ocean
            emp_tl (:,:) = emp_tl (:,:) + zfwftl(jc) / surf(jpncs+1)
            emps_tl(:,:) = emps_tl(:,:) + zfwftl(jc) / surf(jpncs+1)
         ELSEIF( ncstt(jc) == 1 ) THEN
            ! Excess water in open sea, at outflow location, excess evap shared
            IF ( zfwftl(jc) <= 0.e0 ) THEN
                DO jn = 1, ncsnr(jc)
                  ji = mi0(ncsir(jc,jn))
                  jj = mj0(ncsjr(jc,jn)) ! Location of outflow in open ocean
                  IF (      ji > 1 .AND. ji < jpi   &
                      .AND. jj > 1 .AND. jj < jpj ) THEN
                      emp_tl (ji,jj) = emp_tl (ji,jj) + zfwftl(jc) /   &
                         (FLOAT(ncsnr(jc)) * e1t(ji,jj) * e2t(ji,jj))
                      emps_tl(ji,jj) = emps_tl(ji,jj) + zfwftl(jc) /   &
                          (FLOAT(ncsnr(jc)) * e1t(ji,jj) * e2t(ji,jj))
                  END IF
                END DO
            ELSE
                emp_tl (:,:) = emp_tl (:,:) + zfwftl(jc) / surf(jpncs+1)
                emps_tl(:,:) = emps_tl(:,:) + zfwftl(jc) / surf(jpncs+1)
            ENDIF
         ELSEIF( ncstt(jc) == 2 ) THEN
            ! Excess e-p+r (either sign) goes to open ocean, at outflow location
            IF(      ji > 1 .AND. ji < jpi    &
               .AND. jj > 1 .AND. jj < jpj ) THEN
                DO jn = 1, ncsnr(jc)
                  ji = mi0(ncsir(jc,jn))
                  jj = mj0(ncsjr(jc,jn)) ! Location of outflow in open ocean
                  emp_tl (ji,jj) = emp_tl (ji,jj) + zfwftl(jc)   &
                      / (FLOAT(ncsnr(jc)) *  e1t(ji,jj) * e2t(ji,jj) )
                  emps_tl(ji,jj) = emps_tl(ji,jj) + zfwftl(jc)   &
                      / (FLOAT(ncsnr(jc)) *  e1t(ji,jj) * e2t(ji,jj) )
                END DO
            ENDIF
         ENDIF
         !
         DO jj = ncsj1(jc), ncsj2(jc)
            DO ji = ncsi1(jc), ncsi2(jc)
               emp_tl (ji,jj) = emp_tl (ji,jj) - zfwftl(jc) / surf(jc)
               emps_tl(ji,jj) = emps_tl(ji,jj) - zfwftl(jc) / surf(jc)
            END DO
         END DO
         !
      END DO
      !
      CALL lbc_lnk( emp_tl , 'T', 1. )
      CALL lbc_lnk( emps_tl, 'T', 1. )
      !
   END SUBROUTINE sbc_clo_tan
   SUBROUTINE sbc_clo_adj( kt )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_clo_adj  ***
      !!
      !! ** Purpose :   Special handling of closed seas (Adjoint version)
      !!
      !! ** Method  :   Water flux is forced to zero over closed sea
      !!      Excess is shared between remaining ocean, or
      !!      put as run-off in open ocean.
      !!
      !! ** Action  :   emp, emps   updated surface freshwater fluxes at kt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean model time step
      !
      !
      INTEGER                     ::   ji, jj, jc, jn   ! dummy loop indices
      REAL(wp)                    ::   zze2
      REAL(wp), DIMENSION (jpncs) ::   zfwfad
      !!----------------------------------------------------------------------
      !                                                   !------------------!
      IF( kt == nit000 ) THEN                             !  Initialisation  !
         !                                                !------------------!
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*)'sbc_clo_adj : closed seas '
         IF(lwp) WRITE(numout,*)'~~~~~~~'

         ! Total surface of ocean
         surf(jpncs+1) = SUM( e1t(:,:) * e2t(:,:) * tmask_i(:,:) )

         DO jc = 1, jpncs
            surf(jc) =0.0_wp
            DO jj = ncsj1(jc), ncsj2(jc)
               DO ji = ncsi1(jc), ncsi2(jc)
                  surf(jc) = surf(jc) + e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)      ! surface of closed seas
               END DO
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum ( surf, jpncs+1 )       ! mpp: sum over all the global domain

         IF(lwp) WRITE(numout,*)'     Closed sea surfaces'
         DO jc = 1, jpncs
            IF(lwp)WRITE(numout,FMT='(1I3,4I4,5X,F16.2)') jc, ncsi1(jc), ncsi2(jc), ncsj1(jc), ncsj2(jc), surf(jc)
         END DO

         ! jpncs+1 : surface of sea, closed seas excluded
         DO jc = 1, jpncs
            surf(jpncs+1) = surf(jpncs+1) - surf(jc)
         END DO
         !
      ENDIF
      zfwfad = 0.0_wp
      CALL lbc_lnk_adj( emp_ad , 'T', 1. )
      CALL lbc_lnk_adj( emps_ad, 'T', 1. )
      DO jc = jpncs, 1, -1
         !
         DO jj = ncsj1(jc), ncsj2(jc)
            DO ji = ncsi1(jc), ncsi2(jc)
               zfwfad(jc) = zfwfad(jc) - emps_ad(ji,jj) / surf(jc)
               zfwfad(jc) = zfwfad(jc) - emp_ad(ji,jj)  / surf(jc)
            END DO
         END DO
         !
         IF( ncstt(jc) == 0 ) THEN
            ! water/evap excess is shared by all open ocean
            DO jj = 1, jpj
               DO ji = 1, jpi
                  zfwfad(jc) = zfwfad(jc) + emps_ad(ji,jj) / surf(jpncs+1)
                  zfwfad(jc) = zfwfad(jc) + emp_ad(ji,jj)  / surf(jpncs+1)
               END DO
            END DO
         ELSEIF( ncstt(jc) == 1 ) THEN
            ! Excess water in open sea, at outflow location, excess evap shared
            IF ( zfwfad(jc) <= 0.e0 ) THEN
                DO jn = 1, ncsnr(jc)
                  ji = mi0(ncsir(jc,jn))
                  jj = mj0(ncsjr(jc,jn)) ! Location of outflow in open ocean
                  IF (      ji > 1 .AND. ji < jpi   &
                      .AND. jj > 1 .AND. jj < jpj ) THEN
                      zfwfad(jc) = zfwfad(jc) + emps_ad(ji,jj) /   &
                          (FLOAT(ncsnr(jc)) * e1t(ji,jj) * e2t(ji,jj))
                      zfwfad(jc) = zfwfad(jc) + emp_ad(ji,jj)  /   &
                          (FLOAT(ncsnr(jc)) * e1t(ji,jj) * e2t(ji,jj))
                  END IF
                END DO
            ELSE
                DO jj = 1, jpj
                   DO ji = 1, jpi
                      zfwfad(jc) = zfwfad(jc) + emps_ad(ji,jj) / surf(jpncs+1)
                      zfwfad(jc) = zfwfad(jc) + emp_ad(ji,jj)  / surf(jpncs+1)
                   END DO
                END DO
            ENDIF
         ELSEIF( ncstt(jc) == 2 ) THEN
            ! Excess e-p+r (either sign) goes to open ocean, at outflow location
            IF(      ji > 1 .AND. ji < jpi    &
               .AND. jj > 1 .AND. jj < jpj ) THEN
                DO jn = 1, ncsnr(jc)
                  ji = mi0(ncsir(jc,jn))
                  jj = mj0(ncsjr(jc,jn)) ! Location of outflow in open ocean
                  zfwfad(jc) = zfwfad(jc) + emps_ad(ji,jj)  &
                      / (FLOAT(ncsnr(jc)) *  e1t(ji,jj) * e2t(ji,jj) )
                  zfwfad(jc) = zfwfad(jc) + emp_ad(ji,jj)   &
                      / (FLOAT(ncsnr(jc)) *  e1t(ji,jj) * e2t(ji,jj) )
                END DO
            ENDIF
         ENDIF
         !
      END DO
      !                                                   !--------------------!
      !                                                   !  update emp, emps  !
                                                          !--------------------!
      IF( lk_mpp )   CALL mpp_sum ( zfwfad(:) , jpncs )       ! mpp: sum over all the global domain
      DO jc = 1, jpncs
         DO jj = ncsj1(jc), ncsj2(jc)
            DO ji = ncsi1(jc), ncsi2(jc)
               emp_ad(ji,jj) = emp_ad(ji,jj) + e1t(ji,jj) * e2t(ji,jj) * zfwfad(jc) * tmask_i(ji,jj)
            END DO
         END DO
      END DO
      zfwfad = 0.0_wp
   END SUBROUTINE sbc_clo_adj
   SUBROUTINE sbc_clo_adj_tst ( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE sbc_clo_adj_tst ***
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
      !!        ! 09-08 (F. Vigilant)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      !! * Local declarations
      INTEGER ::  &
         & istp,  &
         & jstp,  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zemp_tlin ,     & ! Tangent input
         & zemps_tlin ,    & ! Tangent input
         & zemp_tlout,     & ! Tangent output
         & zemps_tlout,    & ! Tangent output
         & zemp_adin ,     & ! Adjoint input
         & zemps_adin ,    & ! Adjoint input
         & zemp_adout,     & ! Adjoint output
         & zemps_adout,    & ! Adjoint output
         & zr                ! 2D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & zemp_tlin(  jpi,jpj),     &
         & zemps_tlin( jpi,jpj),     &
         & zemp_tlout( jpi,jpj),     &
         & zemps_tlout(jpi,jpj),     &
         & zemp_adin(  jpi,jpj),     &
         & zemps_adin( jpi,jpj),     &
         & zemp_adout( jpi,jpj),     &
         & zemps_adout(jpi,jpj),     &
         & zr(         jpi,jpj)      &
         & )
      !==================================================================
      ! 1) dx = ( emp_tl, emps_tl, ssh_tl ) and
      !    dy = ( emp_tl, emps_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          zemp_tlin  (:,:) = 0.0_wp
          zemps_tlin (:,:) = 0.0_wp
          zemp_tlout (:,:) = 0.0_wp
          zemps_tlout(:,:) = 0.0_wp
          zemp_adin  (:,:) = 0.0_wp
          zemps_adin (:,:) = 0.0_wp
          zemp_adout (:,:) = 0.0_wp
          zemps_adout(:,:) = 0.0_wp
          zr(:,:)          = 0.0_wp
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random( zr, 'T', 0.0_wp, stdemp )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zemp_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stdemp )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zemps_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      emps_tl(:,:) = zemps_tlin(:,:)
      emp_tl (:,:) = zemp_tlin (:,:)

      CALL sbc_clo_tan( nit000 )

      zemps_tlout(:,:) = emps_tl(:,:)
      zemp_tlout (:,:) = emp_tl (:,:)
      !-----------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !-----------------------------------------------------------------

      DO jj = 1, jpi
         DO ji = 1, jpj
            zemp_adin( ji,jj) = zemp_tlout( ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1)
            zemps_adin(ji,jj) = zemps_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1)
         END DO
      END DO
      !-----------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !-----------------------------------------------------------------
      zsp1 = DOT_PRODUCT( zemp_tlout,  zemp_adin  ) &
         & + DOT_PRODUCT( zemps_tlout, zemps_adin )

      !-----------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !-----------------------------------------------------------------

      emp_ad (:,:) = zemp_adin (:,:)
      emps_ad(:,:) = zemps_adin(:,:)

      CALL sbc_clo_adj ( nitend )

      zemps_adout(:,:) = emps_ad(:,:)
      zemp_adout (:,:) = emp_ad (:,:)

      zsp2 = DOT_PRODUCT( zemp_tlin,  zemp_adout  ) &
         & + DOT_PRODUCT( zemps_tlin, zemps_adout )

      ! 14 char:'12345678901234'
      cl_name = 'sbc_clo_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE(       &
         & zemp_tlin,   &
         & zemps_tlin,  &
         & zemp_tlout,  &
         & zemps_tlout, &
         & zemp_adin,   &
         & zemps_adin,  &
         & zemp_adout,  &
         & zemps_adout, &
         & zr           &
         & )
   END SUBROUTINE sbc_clo_adj_tst
   !!======================================================================
END MODULE closea_tam
