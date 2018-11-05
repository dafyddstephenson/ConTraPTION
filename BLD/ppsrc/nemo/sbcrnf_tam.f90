MODULE sbcrnf_tam
   !!======================================================================
   !!                       ***  MODULE  sbcrnf_tam  ***
   !! Ocean forcing:  river runoff
   !!=====================================================================
   !! History :  OPA  ! 2000-11  (R. Hordoir, E. Durand)  NetCDF FORMAT
   !!   NEMO     1.0  ! 2002-09  (G. Madec)  F90: Free form and module
   !!            3.0  ! 2006-07  (G. Madec)  Surface module
   !!            3.2  ! 2009-04  (B. Lemaire)  Introduce iom_put
   !!            3.3  ! 2010-10  (R. Furner, G. Madec) runoff distributed over ocean levels
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_rnf      : monthly runoffs read in a NetCDF file
   !!   sbc_rnf_init : runoffs initialisation
   !!   rnf_mouth    : set river mouth mask
   !!----------------------------------------------------------------------
   USE dom_oce         ! ocean space and time domain
   USE phycst          ! physical constants
   USE sbc_oce         ! surface boundary condition variables
   USE sbc_oce_tam         ! surface boundary condition variables
   USE closea          ! closed seas
   USE fldread         ! read input field at current time step
   USE restart         ! restart
   USE in_out_manager  ! I/O manager
   USE iom             ! I/O module
   USE lib_mpp         ! MPP library
   USE sbcrnf

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_rnf_tan       ! routine call in sbcmod module
   PUBLIC   sbc_rnf_div_tan   ! routine called in sshwzv module
   PUBLIC   sbc_rnf_adj       ! routine call in sbcmod module
   PUBLIC   sbc_rnf_div_adj   ! routine called in sshwzv module
   PUBLIC   sbc_rnf_alloc_tam ! routine call in sbcmod module
   PUBLIC   sbc_rnf_init_tam

   !                                                     !!* namsbc_rnf namelist *

   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rnf_tsc_b_tl, rnf_tsc_tl  !: before and now T & S runoff contents   [K.m/s & PSU.m/s]
   REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:,:) ::   rnf_tsc_b_ad, rnf_tsc_ad  !: before and now T & S runoff contents   [K.m/s & PSU.m/s

   LOGICAL, PRIVATE, SAVE :: ll_alloctl = .FALSE.
   LOGICAL, PRIVATE, SAVE :: ll_allocad = .FALSE.
   REAL(wp) ::   r1_rau0   ! = 1 / rau0

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
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   INTEGER FUNCTION sbc_rnf_alloc_tam( kmode )
      INTEGER, OPTIONAL :: kmode
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_rnf_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: jmode
      INTEGER, DIMENSION(2) :: ierr

      IF ( PRESENT( kmode ) ) THEN
         jmode = kmode
      ELSE
         jmode = 0
      END IF


      IF ( ( jmode == 0 ) .OR. ( jmode == 1 ) .AND. ( .NOT. ll_alloctl ) ) THEN
         ALLOCATE(  rnf_tsc_b_tl(jpi,jpj,jpts) , rnf_tsc_tl (jpi,jpj,jpts), STAT = ierr(1) )
         ll_alloctl = .TRUE.
      END IF
      IF ( ( jmode == 0 ) .OR. ( jmode == 2 ) .AND. ( .NOT. ll_allocad ) ) THEN
         ALLOCATE(  rnf_tsc_b_ad(jpi,jpj,jpts) , rnf_tsc_ad (jpi,jpj,jpts) , STAT=ierr(2) )
         ll_allocad = .TRUE.
      END IF
      sbc_rnf_alloc_tam = SUM( ierr )
         !
      IF( lk_mpp            )   CALL mpp_sum ( sbc_rnf_alloc_tam )
      IF( sbc_rnf_alloc_tam > 0 )   CALL ctl_warn('sbc_rnf_alloc: allocation of arrays failed')
   END FUNCTION sbc_rnf_alloc_tam

   INTEGER FUNCTION sbc_rnf_dealloc_tam( kmode )
      INTEGER, OPTIONAL :: kmode
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE sbc_rnf_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER :: jmode
      INTEGER, DIMENSION(2) :: ierr

      IF ( PRESENT( kmode ) ) THEN
         jmode = kmode
      ELSE
         jmode = 0
      END IF


      IF ( ( jmode == 0 ) .OR. ( jmode == 1 ) .AND. ( ll_alloctl ) ) THEN
         DEALLOCATE(  rnf_tsc_b_tl, rnf_tsc_tl, STAT = ierr(1) )
         ll_alloctl = .FALSE.
      END IF
      IF ( ( jmode == 0 ) .OR. ( jmode == 2 ) .AND. ( ll_allocad ) ) THEN
         DEALLOCATE(  rnf_tsc_b_ad, rnf_tsc_ad, STAT=ierr(2) )
         ll_allocad = .FALSE.
      END IF
      sbc_rnf_dealloc_tam = SUM( ierr )
         !
      IF( lk_mpp            )   CALL mpp_sum ( sbc_rnf_dealloc_tam )
      IF( sbc_rnf_dealloc_tam > 0 )   CALL ctl_warn('sbc_rnf_dealloc: deallocation of arrays failed')
   END FUNCTION sbc_rnf_dealloc_tam

   SUBROUTINE sbc_rnf_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf_tan  ***
      !!
      !! ** Purpose :   Introduce a climatological run off forcing
      !!
      !! ** Method  :   Set each river mouth with a monthly climatology
      !!                provided from different data.
      !!                CAUTION : upward water flux, runoff forced to be < 0
      !!
      !! ** Action  :   runoff updated runoff field at time-step kt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      !                                            ! ---------------------------------------- !
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         rnf_b_tl    (:,:  ) = rnf_tl    (:,:  )               ! Swap the ocean forcing fields except at nit000
         rnf_tsc_b_tl(:,:,:) = rnf_tsc_tl(:,:,:)               ! where before fields are set at the end of the routine
         !
      ENDIF

      !                                                   !-------------------!
      IF( .NOT. ln_rnf_emp ) THEN                         !   Update runoff   !
         !                                                !-------------------!
         !
                             CALL fld_read ( kt, nn_fsbc, sf_rnf   )    ! Read Runoffs data and provide it at kt
         IF( ln_rnf_tem  )   CALL fld_read ( kt, nn_fsbc, sf_t_rnf )    ! idem for runoffs temperature if required
         IF( ln_rnf_sal  )   CALL fld_read ( kt, nn_fsbc, sf_s_rnf )    ! idem for runoffs salinity    if required
         !
         ! Runoff reduction only associated to the ORCA2_LIM configuration
         ! when reading the NetCDF file runoff_1m_nomask.nc
         IF( cp_cfg == 'orca' .AND. jp_cfg == 2 )   THEN
            WHERE( 40._wp < gphit(:,:) .AND. gphit(:,:) < 65._wp )
               sf_rnf(1)%fnow(:,:,1) = 0.85 * sf_rnf(1)%fnow(:,:,1)
            END WHERE
         ENDIF
         !
         IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
            rnf_tl(:,:) = 0.0_wp
            !
            r1_rau0 = 1._wp / rau0
            !                                                     ! set temperature & salinity content of runoffs
            IF( ln_rnf_tem ) THEN                                       ! use runoffs temperature data
               rnf_tsc_tl(:,:,jp_tem) = ( sf_t_rnf(1)%fnow(:,:,1) ) * rnf_tl(:,:) * r1_rau0
               WHERE( sf_t_rnf(1)%fnow(:,:,1) == -999 )                 ! if missing data value use SST as runoffs temperature
               rnf_tsc_tl(:,:,jp_tem) = (sst_m_tl(:,:) * rnf(:,:) &
                  &                      + sst_m(:,:) * rnf_tl(:,:)) * r1_rau0
               END WHERE
            ELSE                                                        ! use SST as runoffs temperature
               rnf_tsc_tl(:,:,jp_tem) = (sst_m_tl(:,:) * rnf(:,:) &
                  &                      + sst_m(:,:) * rnf_tl(:,:)) * r1_rau0
            ENDIF
            !                                                           ! use runoffs salinity data
            IF( ln_rnf_sal )   rnf_tsc_tl(:,:,jp_sal) = ( sf_s_rnf(1)%fnow(:,:,1) ) * rnf_tl(:,:) * r1_rau0
            !                                                           ! else use S=0 for runoffs (done one for all in the init)
            !
            IF( ln_rnf_tem .OR. ln_rnf_sal ) THEN                 ! runoffs as outflow: use ocean SST and SSS
               WHERE( rnf(:,:) < 0._wp )                                 ! example baltic model when flow is out of domain
               rnf_tsc_tl(:,:,jp_tem) = (sst_m_tl(:,:) * rnf(:,:) &
                  &                      + sst_m(:,:) * rnf_tl(:,:)) * r1_rau0
               rnf_tsc_tl(:,:,jp_sal) = (sss_m_tl(:,:) * rnf(:,:) &
                  &                      + sss_m(:,:) * rnf_tl(:,:)) * r1_rau0
               END WHERE
            ENDIF
            !
         ENDIF
         !
      ENDIF
      !
      IF( kt == nit000 ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         IF( ln_rstart ) THEN
            rnf_b_tl(:,:) = 0.0_wp
            rnf_tsc_b_tl(:,:,:) = 0.0_wp
         ELSE                                                   !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields set to nit000'
             rnf_b_tl    (:,:  ) = rnf_tl    (:,:  )
             rnf_tsc_b_tl(:,:,:) = rnf_tsc_tl(:,:,:)
         ENDIF
      ENDIF
      !
   END SUBROUTINE sbc_rnf_tan


   SUBROUTINE sbc_rnf_div_tan( phdivn_tl )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :
      !!                CAUTION : rnf is positive (inflow) decreasing the
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phdivn_tl   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   r1_rau0   ! local scalar
      REAL(wp) ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact = 0.5_wp
      !
      r1_rau0 = 1._wp / rau0
      IF( ln_rnf_depth ) THEN      !==   runoff distributed over several levels   ==!
         IF( lk_vvl ) THEN             ! variable volume case
            !DO jj = 1, jpj                   ! update the depth over which runoffs are distributed
               !DO ji = 1, jpi
                  !h_rnf(ji,jj) = 0._wp
                  !DO jk = 1, nk_rnf(ji,jj)                           ! recalculates h_rnf to be the depth in metres
                     !h_rnf(ji,jj) = h_rnf(ji,jj) + e3t(ji,jj,jk)   ! to the bottom of the relevant grid box
                  !END DO
                  !!                          ! apply the runoff input flow
                  !DO jk = 1, nk_rnf(ji,jj)
                     !phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact * r1_rau0 / h_rnf(ji,jj)
                  !END DO
               !END DO
            !END DO
            CALL ctl_stop('key_vvl not implemented in TAM yet')
         ELSE                          ! constant volume case : just apply the runoff input flow
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     phdivn_tl(ji,jj,jk) = phdivn_tl(ji,jj,jk) - ( rnf_tl(ji,jj) + rnf_b_tl(ji,jj) )   &
                                         &                       * zfact * r1_rau0 / h_rnf(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
      ELSE                       !==   runoff put only at the surface   ==!
         IF( lk_vvl ) THEN              ! variable volume case
            CALL ctl_stop('key_vvl not implemented in TAM yet')
            !h_rnf(:,:) = e3t(:,:,1)   ! recalculate h_rnf to be depth of top box
         ENDIF
         phdivn_tl(:,:,1) = phdivn_tl(:,:,1) - ( rnf_tl(:,:) + rnf_b_tl(:,:) ) * zfact * r1_rau0 / e3t(:,:,1)
      ENDIF
      !
   END SUBROUTINE sbc_rnf_div_tan

   SUBROUTINE sbc_rnf_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf_adj  ***
      !!
      !! ** Purpose :   Introduce a climatological run off forcing
      !!
      !! ** Method  :   Set each river mouth with a monthly climatology
      !!                provided from different data.
      !!                CAUTION : upward water flux, runoff forced to be < 0
      !!
      !! ** Action  :   runoff updated runoff field at time-step kt
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt          ! ocean time step
      !!
      INTEGER  ::   ji, jj   ! dummy loop indices
      !!----------------------------------------------------------------------
      !
      IF( kt == nitend ) THEN                          !   set the forcing field at nit000 - 1    !
         !                                             ! ---------------------------------------- !
         IF( ln_rstart ) THEN
            rnf_b_ad(:,:) = 0.0_wp
            rnf_tsc_b_ad(:,:,:) = 0.0_wp
         ELSE                                                   !* no restart: set from nit000 values
            IF(lwp) WRITE(numout,*) '          nit000-1 runoff forcing fields set to nit000'
             rnf_ad    (:,:  ) = rnf_ad    (:,:  )  + rnf_b_ad(:,:)
             rnf_tsc_ad(:,:,:) = rnf_tsc_ad(:,:,:)  + rnf_tsc_b_ad(:,:,:)
         ENDIF
      ENDIF
      !                                            ! ---------------------------------------- !
      !                                                   !-------------------!
      IF( .NOT. ln_rnf_emp ) THEN                         !   Update runoff   !
         !                                                !-------------------!
         !
                             CALL fld_read ( kt, nn_fsbc, sf_rnf   )    ! Read Runoffs data and provide it at kt
         IF( ln_rnf_tem  )   CALL fld_read ( kt, nn_fsbc, sf_t_rnf )    ! idem for runoffs temperature if required
         IF( ln_rnf_sal  )   CALL fld_read ( kt, nn_fsbc, sf_s_rnf )    ! idem for runoffs salinity    if required
         !
         ! Runoff reduction only associated to the ORCA2_LIM configuration
         ! when reading the NetCDF file runoff_1m_nomask.nc
         IF( cp_cfg == 'orca' .AND. jp_cfg == 2 )   THEN
            WHERE( 40._wp < gphit(:,:) .AND. gphit(:,:) < 65._wp )
               sf_rnf(1)%fnow(:,:,1) = 0.85 * sf_rnf(1)%fnow(:,:,1)
            END WHERE
         ENDIF
         !
         IF( MOD( kt - 1, nn_fsbc ) == 0 ) THEN
            rnf_ad(:,:) = 0.0_wp
            !
            r1_rau0 = 1._wp / rau0
            !                                                     ! set temperature & salinity content of runoffs
            !
            IF( ln_rnf_tem .OR. ln_rnf_sal ) THEN                 ! runoffs as outflow: use ocean SST and SSS
               WHERE( rnf(:,:) < 0._wp )                                 ! example baltic model when flow is out of domain
                  sss_m_ad(:,:) = sst_m_ad(:,:) + r1_rau0 * (rnf(:,:) * rnf_tsc_ad(:,:,jp_sal))
                  rnf_ad(:,:)   = rnf_ad(:,:)   + r1_rau0 * (sss_m(:,:) * rnf_tsc_ad(:,:,jp_sal))
                  sst_m_ad(:,:) = sst_m_ad(:,:) + r1_rau0 * (rnf(:,:) * rnf_tsc_ad(:,:,jp_tem))
                  rnf_ad(:,:)   = rnf_ad(:,:)   + r1_rau0 * (sst_m(:,:) * rnf_tsc_ad(:,:,jp_tem))
               END WHERE
            ENDIF
            IF( ln_rnf_sal ) THEN
               rnf_ad(:,:) = rnf_ad(:,:) + ( sf_s_rnf(1)%fnow(:,:,1) ) * rnf_tsc_ad(:,:,jp_sal) * r1_rau0
            END IF
            !                                                           ! else use S=0 for runoffs (done one for all in the init)
            IF( ln_rnf_tem ) THEN                                       ! use runoffs temperature data
               rnf_ad(:,:) = rnf_ad(:,:) + ( sf_s_rnf(1)%fnow(:,:,1) ) * rnf_tsc_ad(:,:,jp_tem) * r1_rau0
               WHERE( sf_t_rnf(1)%fnow(:,:,1) == -999 )                 ! if missing data value use SST as runoffs temperature
                  sst_m_ad(:,:) = sst_m_ad(:,:) + r1_rau0 * (rnf(:,:) * rnf_tsc_ad(:,:,jp_tem))
                  rnf_ad(:,:)   = rnf_ad(:,:)   + r1_rau0 * (sst_m(:,:) * rnf_tsc_ad(:,:,jp_tem))
               END WHERE
            ELSE                                                        ! use SST as runoffs temperature
               sst_m_ad(:,:) = sst_m_ad(:,:) + r1_rau0 * (rnf(:,:) * rnf_tsc_ad(:,:,jp_tem))
               rnf_ad(:,:)   = rnf_ad(:,:)   + r1_rau0 * (sst_m(:,:) * rnf_tsc_ad(:,:,jp_tem))
            ENDIF
            !                                                           ! use runoffs salinity data
            !
         ENDIF
         !
      ENDIF
      IF( kt /= nit000 ) THEN                      !          Swap of forcing fields          !
         !                                         ! ---------------------------------------- !
         rnf_ad    (:,:  ) = rnf_ad    (:,:  ) + rnf_b_ad(:,:)               ! Swap the ocean forcing fields except at nit000
         rnf_tsc_ad(:,:,:) = rnf_tsc_ad(:,:,:) +rnf_tsc_b_ad(:,:,:)              ! where before fields are set at the end of the routine
         !
      ENDIF

      !
   END SUBROUTINE sbc_rnf_adj


   SUBROUTINE sbc_rnf_div_adj( phdivn_ad )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf  ***
      !!
      !! ** Purpose :   update the horizontal divergence with the runoff inflow
      !!
      !! ** Method  :
      !!                CAUTION : rnf is positive (inflow) decreasing the
      !!                          divergence and expressed in m/s
      !!
      !! ** Action  :   phdivn   decreased by the runoff inflow
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   phdivn_ad   ! horizontal divergence
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   r1_rau0   ! local scalar
      REAL(wp) ::   zfact     ! local scalar
      !!----------------------------------------------------------------------
      !
      zfact = 0.5_wp
      !
      r1_rau0 = 1._wp / rau0
      IF( ln_rnf_depth ) THEN      !==   runoff distributed over several levels   ==!
         IF( lk_vvl ) THEN             ! variable volume case
            !DO jj = 1, jpj                   ! update the depth over which runoffs are distributed
               !DO ji = 1, jpi
                  !h_rnf(ji,jj) = 0._wp
                  !DO jk = 1, nk_rnf(ji,jj)                           ! recalculates h_rnf to be the depth in metres
                     !h_rnf(ji,jj) = h_rnf(ji,jj) + e3t(ji,jj,jk)   ! to the bottom of the relevant grid box
                  !END DO
                  !!                          ! apply the runoff input flow
                  !DO jk = 1, nk_rnf(ji,jj)
                     !phdivn(ji,jj,jk) = phdivn(ji,jj,jk) - ( rnf(ji,jj) + rnf_b(ji,jj) ) * zfact * r1_rau0 / h_rnf(ji,jj)
                  !END DO
               !END DO
            !END DO
            CALL ctl_stop('key_vvl not implemented in TAM yet')
         ELSE                          ! consadjt volume case : just apply the runoff input flow
            DO jj = 1, jpj
               DO ji = 1, jpi
                  DO jk = 1, nk_rnf(ji,jj)
                     rnf_ad(:,:) = rnf_ad(:,:) - phdivn_ad(:,:,1) * zfact * r1_rau0 / h_rnf(ji,jj)
                     rnf_b_ad(:,:) = rnf_b_ad(:,:) - phdivn_ad(:,:,1) * zfact * r1_rau0 / h_rnf(ji,jj)
                  END DO
               END DO
            END DO
         ENDIF
      ELSE                       !==   runoff put only at the surface   ==!
         IF( lk_vvl ) THEN              ! variable volume case
            CALL ctl_stop('key_vvl not implemented in TAM yet')
            !h_rnf(:,:) = e3t(:,:,1)   ! recalculate h_rnf to be depth of top box
         ENDIF
         rnf_ad(:,:) = rnf_ad(:,:) - phdivn_ad(:,:,1) * zfact * r1_rau0 / e3t(:,:,1)
         rnf_b_ad(:,:) = rnf_b_ad(:,:) - phdivn_ad(:,:,1) * zfact * r1_rau0 / e3t(:,:,1)
      ENDIF
      !
   END SUBROUTINE sbc_rnf_div_adj

   SUBROUTINE sbc_rnf_init_tam
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_rnf_init  ***
      !!
      !! ** Purpose :   Initialisation of the runoffs if (ln_rnf=T)
      !!
      !! ** Method  : - read the runoff namsbc_rnf namelist
      !!
      !! ** Action  : - read parameters
      !!----------------------------------------------------------------------
      CHARACTER(len=32) ::   rn_dep_file   ! runoff file name
      INTEGER           ::   ji, jj, jk    ! dummy loop indices
      INTEGER           ::   ierror, inum  ! temporary integer
      !!
      IF(lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'sbc_rnf_tam : runoff '
         WRITE(numout,*) '~~~~~~~ '
         WRITE(numout,*) '   Namelist namsbc_rnf'
         WRITE(numout,*) '      runoff in a file to be read                ln_rnf_emp   = ', ln_rnf_emp
         WRITE(numout,*) '      specific river mouths treatment            ln_rnf_mouth = ', ln_rnf_mouth
         WRITE(numout,*) '      river mouth additional Kz                  rn_avt_rnf   = ', rn_avt_rnf
         WRITE(numout,*) '      depth of river mouth additional mixing     rn_hrnf      = ', rn_hrnf
         WRITE(numout,*) '      multiplicative factor for runoff           rn_rfact     = ', rn_rfact
      ENDIF

      !                                   ! ==================
      !                                   !   Type of runoff
      !                                   ! ==================
      !                                         !==  allocate runoff arrays
      IF( sbc_rnf_alloc_tam() /= 0 )   CALL ctl_stop( 'STOP', 'sbc_rnf_alloc_tam : unable to allocate arrays' )
      !
      !
      rnf_tl(:,:) =  0._wp                         ! runoff initialisation
      rnf_ad(:,:) =  0._wp                         ! runoff initialisation
      rnf_tsc_tl(:,:,:) = 0._wp                    ! runoffs temperature & salinty contents initilisation
      rnf_tsc_ad(:,:,:) = 0._wp                    ! runoffs temperature & salinty contents initilisation
      rnf_tsc_b_tl(:,:,:) = 0._wp                    ! runoffs temperature & salinty contents initilisation
      rnf_tsc_b_ad(:,:,:) = 0._wp                    ! runoffs temperature & salinty contents initilisation
      !
   END SUBROUTINE sbc_rnf_init_tam
   !!======================================================================
END MODULE sbcrnf_tam
