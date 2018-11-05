MODULE trasbc_tam
   !!==============================================================================
   !!                       ***  MODULE  trasbc  ***
   !! Ocean active tracers:  surface boundary condition
   !!==============================================================================
   !! History :  OPA  !  1998-10  (G. Madec, G. Roullet, M. Imbard)  Original code
   !!            8.2  !  2001-02  (D. Ludicone)  sea ice and free surface
   !!  NEMO      1.0  !  2002-06  (G. Madec)  F90: Free form and module
   !!            3.3  !  2010-04  (M. Leclair, G. Madec)  Forcing averaged over 2 time steps
   !!             -   !  2010-09  (C. Ethe, G. Madec) Merge TRA-TRC
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_sbc      : update the tracer trend at ocean surface
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and active tracers
   USE oce_tam
   USE sbc_oce         ! surface boundary condition: ocean
   USE sbc_oce_tam
   USE dom_oce         ! ocean space domain variables
   USE phycst          ! physical constant
   USE traqsr          ! solar radiation penetration
   USE traqsr_tam
   USE trdmod_oce      ! ocean trends
   USE trdtra          ! ocean trends
   USE in_out_manager  ! I/O manager
   USE prtctl          ! Print control
   USE restart         ! ocean restart
   USE sbcrnf          ! River runoff
   USE sbcrnf_tam          ! River runoff
   USE sbcmod          ! ln_rnf
   USE iom
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lbclnk_tam          ! ocean lateral boundary conditions (or mpp link)
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE tstool_tam
   USE paresp
   USE dotprodfld
   USE gridrandom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_sbc_tan    ! routine called by step.F90
   PUBLIC   tra_sbc_adj    ! routine called by step.F90
   PUBLIC   tra_sbc_adj_tst

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
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_sbc_tan ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sbc_tan  ***
      !!
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (flux through the interface, concentration/dilution effect)
      !!      and add it to the general trend of tracer equations.
      !!
      !! ** Method :
      !!      Following Roullet and Madec (2000), the air-sea flux can be divided
      !!      into three effects: (1) Fext, external forcing;
      !!      (2) Fwi, concentration/dilution effect due to water exchanged
      !!         at the surface by evaporation, precipitations and runoff (E-P-R);
      !!      (3) Fwe, tracer carried with the water that is exchanged.
      !!
      !!      Fext, flux through the air-sea interface for temperature and salt:
      !!            - temperature : heat flux q (w/m2). If penetrative solar
      !!         radiation q is only the non solar part of the heat flux, the
      !!         solar part is added in traqsr.F routine.
      !!            ta = ta + q /(rau0 rcp e3t)  for k=1
      !!            - salinity    : no salt flux
      !!
      !!      The formulation for Fwb and Fwi vary according to the free
      !!      surface formulation (linear or variable volume).
      !!      * Linear free surface
      !!            The surface freshwater flux modifies the ocean volume
      !!         and thus the concentration of a tracer and the temperature.
      !!         First order of the effect of surface freshwater exchange
      !!         for salinity, it can be neglected on temperature (especially
      !!         as the temperature of precipitations and runoffs is usually
      !!         unknown).
      !!            - temperature : we assume that the temperature of both
      !!         precipitations and runoffs is equal to the SST, thus there
      !!         is no additional flux since in this case, the concentration
      !!         dilution effect is balanced by the net heat flux associated
      !!         to the freshwater exchange (Fwe+Fwi=0):
      !!            (Tp P - Te E) + SST (P-E) = 0 when Tp=Te=SST
      !!            - salinity    : evaporation, precipitation and runoff
      !!         water has a zero salinity (Fwe=0), thus only Fwi remains:
      !!            sa = sa + emp * sn / e3t   for k=1
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation minus runoff) given in kg/m2/s is divided
      !!         by 1035 kg/m3 (density of ocena water) to obtain m/s.
      !!         Note: even though Fwe does not appear explicitly for
      !!         temperature in this routine, the heat carried by the water
      !!         exchanged through the surface is part of the total heat flux
      !!         forcing and must be taken into account in the global heat
      !!         balance).
      !!      * nonlinear free surface (variable volume, lk_vvl)
      !!         contrary to the linear free surface case, Fwi is properly
      !!         taken into account by using the true layer thicknesses to
      !!         calculate tracer content and advection. There is no need to
      !!         deal with it in this routine.
      !!           - temperature: Fwe=SST (P-E+R) is added to Fext.
      !!           - salinity:  Fwe = 0, there is no surface flux of salt.
      !!
      !! ** Action  : - Update the 1st level of (ta,sa) with the trend associated
      !!                with the tracer surface boundary condition
      !!              - save the trend it in ttrd ('key_trdtra')
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      REAL(wp) ::   zfact, z1_e3t, zsrau, zdep
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_sbc_tan')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_sbc_tan : TRAcer Surface Boundary Condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      zsrau = 1. / rau0             ! initialization

      IF( .NOT.ln_traqsr ) THEN     ! no solar radiation penetration
         qns_tl(:,:) = qns_tl(:,:) + qsr_tl(:,:)      ! total heat flux in qns
         qsr_tl(:,:) = 0.e0                     ! qsr set to zero
      ENDIF
      !
      !----------------------------------------
      !        EMP, EMPS and QNS effects
      !----------------------------------------
      !                                             Set before sbc tracer content fields
      !                                             ************************************
      IF( kt == nit000 ) THEN                      ! Set the forcing field at nit000 - 1
         !                                         ! -----------------------------------
         IF( ln_rstart )  THEN
            zfact = 0.5e0
            sbc_tsc_b_tl(:,:,:) = 0.0_wp
         ELSE                                         ! No restart or restart not found: Euler forward time stepping
            zfact = 1.e0
            sbc_tsc_b_tl(:,:,:) = 0.0_wp
         ENDIF
      ELSE                                         ! Swap of forcing fields
         !                                         ! ----------------------
         zfact = 0.5e0
         sbc_tsc_b_tl(:,:,:) = sbc_tsc_tl(:,:,:)
      ENDIF
      !                                             Compute now sbc tracer content fields
      !                                             *************************************

                                                   ! Concentration dilution effect on (t,s) due to
                                                   ! evaporation, precipitation and qns, but not river runoff

      IF( lk_vvl ) THEN                            ! Variable Volume case
         !DO jj = 1, jpj
            !DO ji = 1, jpi
               !! temperature : heat flux + cooling/heating effet of EMP flux
               !sbc_tsc(ji,jj,jp_tem) = ro0cpr * qns(ji,jj) - zsrau * emp(ji,jj) * tsn(ji,jj,1,jp_tem)
               !! concent./dilut. effect due to sea-ice melt/formation and (possibly) SSS restoration
               !sbc_tsc(ji,jj,jp_sal) = ( emps(ji,jj) - emp(ji,jj) ) * zsrau * tsn(ji,jj,1,jp_sal)
            !END DO
         !END DO
         CALL ctl_stop('key_vvl not implemented in TAM yet')
      ELSE                                         ! Constant Volume case
         DO jj = 2, jpj
            DO ji = 2, jpim1   ! vector opt.
               ! temperature : heat flux
               sbc_tsc_tl(ji,jj,jp_tem) = ro0cpr * qns_tl(ji,jj)
               ! salinity    : salt flux + concent./dilut. effect (both in emps)
               sbc_tsc_tl(ji,jj,jp_sal) = zsrau * ( emps_tl(ji,jj) * tsn(ji,jj,1,jp_sal) &
                  &                          + emps(ji,jj) * tsn_tl(ji,jj,1,jp_sal) )
            END DO
         END DO
      ENDIF
      ! Concentration dilution effect on (t,s) due to evapouration, precipitation and qns, but not river runoff
      DO jn = 1, jpts
         DO jj = 2, jpj
            DO ji = 2, jpim1   ! vector opt.
               z1_e3t = zfact / e3t(ji,jj,1)
               tsa_tl(ji,jj,1,jn) = tsa_tl(ji,jj,1,jn) + ( sbc_tsc_b_tl(ji,jj,jn) + sbc_tsc_tl(ji,jj,jn) ) * z1_e3t
            END DO
         END DO
      END DO
      !
      !----------------------------------------
      !        River Runoff effects
      !----------------------------------------
      !
      zfact = 0.5e0
      !
      ! Effect on (t,s) due to river runoff (dilution effect automatically applied via vertical tracer advection)
      IF( ln_rnf ) THEN
         DO jj = 2, jpj
            DO ji = 2, jpim1
               zdep = 1. / h_rnf(ji,jj)
               zdep = zfact * zdep
               IF ( rnf(ji,jj) /= 0._wp ) THEN
                  DO jk = 1, nk_rnf(ji,jj)
                                        tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)   &
                                          &               +  ( rnf_tsc_b_tl(ji,jj,jp_tem) + rnf_tsc_tl(ji,jj,jp_tem) ) * zdep
                     IF( ln_rnf_sal )   tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal)   &
                                          &               +  ( rnf_tsc_b_tl(ji,jj,jp_sal) + rnf_tsc_tl(ji,jj,jp_sal) ) * zdep
                  END DO
               ENDIF
            END DO
         END DO
      ENDIF
!!gm  It should be useless
      CALL lbc_lnk( tsa_tl(:,:,:,jp_tem), 'T', 1. )    ;    CALL lbc_lnk( tsa_tl(:,:,:,jp_sal), 'T', 1. )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_sbc_tan')
      !
   END SUBROUTINE tra_sbc_tan

   SUBROUTINE tra_sbc_adj ( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_sbc_adj  ***
      !!
      !! ** Purpose :   Compute the tracer surface boundary condition trend of
      !!      (flux through the interface, concentration/dilution effect)
      !!      and add it to the general trend of tracer equations.
      !!
      !! ** Method :
      !!      Following Roullet and Madec (2000), the air-sea flux can be divided
      !!      into three effects: (1) Fext, external forcing;
      !!      (2) Fwi, concentration/dilution effect due to water exchanged
      !!         at the surface by evaporation, precipitations and runoff (E-P-R);
      !!      (3) Fwe, tracer carried with the water that is exchanged.
      !!
      !!      Fext, flux through the air-sea interface for temperature and salt:
      !!            - temperature : heat flux q (w/m2). If penetrative solar
      !!         radiation q is only the non solar part of the heat flux, the
      !!         solar part is added in traqsr.F routine.
      !!            ta = ta + q /(rau0 rcp e3t)  for k=1
      !!            - salinity    : no salt flux
      !!
      !!      The formulation for Fwb and Fwi vary according to the free
      !!      surface formulation (linear or variable volume).
      !!      * Linear free surface
      !!            The surface freshwater flux modifies the ocean volume
      !!         and thus the concentration of a tracer and the temperature.
      !!         First order of the effect of surface freshwater exchange
      !!         for salinity, it can be neglected on temperature (especially
      !!         as the temperature of precipitations and runoffs is usually
      !!         unknown).
      !!            - temperature : we assume that the temperature of both
      !!         precipitations and runoffs is equal to the SST, thus there
      !!         is no additional flux since in this case, the concentration
      !!         dilution effect is balanced by the net heat flux associated
      !!         to the freshwater exchange (Fwe+Fwi=0):
      !!            (Tp P - Te E) + SST (P-E) = 0 when Tp=Te=SST
      !!            - salinity    : evaporation, precipitation and runoff
      !!         water has a zero salinity (Fwe=0), thus only Fwi remains:
      !!            sa = sa + emp * sn / e3t   for k=1
      !!         where emp, the surface freshwater budget (evaporation minus
      !!         precipitation minus runoff) given in kg/m2/s is divided
      !!         by 1035 kg/m3 (density of ocena water) to obtain m/s.
      !!         Note: even though Fwe does not appear explicitly for
      !!         temperature in this routine, the heat carried by the water
      !!         exchanged through the surface is part of the total heat flux
      !!         forcing and must be taken into account in the global heat
      !!         balance).
      !!      * nonlinear free surface (variable volume, lk_vvl)
      !!         contrary to the linear free surface case, Fwi is properly
      !!         taken into account by using the true layer thicknesses to
      !!         calculate tracer content and advection. There is no need to
      !!         deal with it in this routine.
      !!           - temperature: Fwe=SST (P-E+R) is added to Fext.
      !!           - salinity:  Fwe = 0, there is no surface flux of salt.
      !!
      !! ** Action  : - Update the 1st level of (ta,sa) with the trend associated
      !!                with the tracer surface boundary condition
      !!              - save the trend it in ttrd ('key_trdtra')
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      REAL(wp) ::   zfact, z1_e3t, zsrau, zdep
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_sbc_adj')
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_sbc_adj : TRAcer Surface Boundary Condition'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF

      zsrau = 1. / rau0             ! initialization
      !!gm  It should be useless
      CALL lbc_lnk_adj( tsa_ad(:,:,:,jp_tem), 'T', 1. )    ;    CALL lbc_lnk_adj( tsa_ad(:,:,:,jp_sal), 'T', 1. )
      !
      !----------------------------------------
      !        River Runoff effects
      !----------------------------------------
      !
      zfact = 0.5e0

      ! Effect on (t,s) due to river runoff (dilution effect automatically applied via vertical tracer advection)
      IF( ln_rnf ) THEN
         DO jj = 2, jpj
            DO ji = 2, jpim1
               zdep = 1. / h_rnf(ji,jj)
               zdep = zfact * zdep
               IF ( rnf(ji,jj) /= 0._wp ) THEN
                  DO jk = 1, nk_rnf(ji,jj)
                     rnf_tsc_b_ad(ji,jj,jp_tem) = rnf_tsc_b_ad(ji,jj,jp_tem) + tsa_ad(ji,jj,jk,jp_tem) * zdep
                     rnf_tsc_ad(ji,jj,jp_tem) = rnf_tsc_ad(ji,jj,jp_tem) + tsa_ad(ji,jj,jk,jp_tem) * zdep
                     IF( ln_rnf_sal ) THEN
                        rnf_tsc_b_ad(ji,jj,jp_sal) = rnf_tsc_b_ad(ji,jj,jp_sal) + tsa_ad(ji,jj,jk,jp_sal) * zdep
                        rnf_tsc_ad(ji,jj,jp_sal) = rnf_tsc_ad(ji,jj,jp_sal) + tsa_ad(ji,jj,jk,jp_sal) * zdep
                     ENDIF
                  END DO
               ENDIF
            END DO
         END DO
      ENDIF
      !
      IF( kt == nit000 ) THEN                      ! Set the forcing field at nit000 - 1
         !                                         ! -----------------------------------
         IF( ln_rstart )  THEN
            zfact = 0.5e0
         ELSE                                         ! No restart or restart not found: Euler forward time stepping
            zfact = 1.e0
         ENDIF
      ENDIF
      !                                          Set before sbc tracer content fields
      !                                          ************************************
      ! Concentration dilution effect on (t,s) due to evapouration, precipitation and qns, but not river runoff
      DO jn = 1, jpts
         DO jj = 2, jpj
            DO ji = 2, jpim1   ! vector opt.
               z1_e3t = zfact / e3t(ji,jj,1)
               sbc_tsc_b_ad(ji,jj,jn) = sbc_tsc_b_ad(ji,jj,jn) + tsa_ad(ji,jj,1,jn) * z1_e3t
               sbc_tsc_ad(ji,jj,jn) = sbc_tsc_ad(ji,jj,jn) + tsa_ad(ji,jj,1,jn) * z1_e3t
            END DO
         END DO
      END DO
      !
      !                                          Compute now sbc tracer content fields
      !                                          *************************************

                                                   ! Concentration dilution effect on (t,s) due to
                                                   ! evaporation, precipitation and qns, but not river runoff

      IF( lk_vvl ) THEN                            ! Variable Volume case
         !DO jj = 1, jpj
            !DO ji = 1, jpi
               !! temperature : heat flux + cooling/heating effet of EMP flux
               !sbc_tsc(ji,jj,jp_tem) = ro0cpr * qns(ji,jj) - zsrau * emp(ji,jj) * tsn(ji,jj,1,jp_tem)
               !! concent./dilut. effect due to sea-ice melt/formation and (possibly) SSS restoration
               !sbc_tsc(ji,jj,jp_sal) = ( emps(ji,jj) - emp(ji,jj) ) * zsrau * tsn(ji,jj,1,jp_sal)
            !END DO
         !END DO
         CALL ctl_stop('key_vvl not implemented in TAM yet')
      ELSE                                         ! Constant Volume case
         DO jj = 2, jpj
            DO ji = 2, jpim1   ! vector opt.
               ! salinity    : salt flux + concent./dilut. effect (both in emps)
               emps_ad(ji,jj) = emps_ad(ji,jj) + zsrau * (sbc_tsc_ad(ji,jj,jp_sal) * tsn(ji,jj,1,jp_sal))
               tsn_ad(ji,jj,1,jp_sal) = tsn_ad(ji,jj,1,jp_sal) + zsrau * (sbc_tsc_ad(ji,jj,jp_sal) * emps(ji,jj))
               ! temperature : heat flux
               qns_ad(ji,jj) = qns_ad(ji,jj) + ro0cpr * sbc_tsc_ad(ji,jj,jp_tem)
               sbc_tsc_ad(ji,jj,jp_sal) = 0._wp
               sbc_tsc_ad(ji,jj,jp_tem) = 0._wp
            END DO
         END DO
      ENDIF

      IF (kt /= nit000 ) THEN                       ! Swap of forcing fields
         !                                         ! ----------------------
         sbc_tsc_ad(:,:,:) = sbc_tsc_ad(:,:,:) + sbc_tsc_b_ad(:,:,:)
         sbc_tsc_b_ad(:,:,:) = 0._wp
      ELSE
         sbc_tsc_b_ad(:,:,:) = 0._wp
      ENDIF
      !
      !----------------------------------------
      !        EMP, EMPS and QNS effects
      !----------------------------------------
      !
      IF( .NOT.ln_traqsr ) THEN     ! no solar radiation penetration
         qsr_ad(:,:) = qsr_ad(:,:) + qns_ad(:,:)
         qns_ad(:,:) = 0._wp
      ENDIF
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_sbc_adj')
      !
   END SUBROUTINE tra_sbc_adj
   SUBROUTINE tra_sbc_adj_tst ( kumadt )
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
      USE trj_tam
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER :: &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zsn_tlin,       &! Tangent input : now salinity
         & zsa_tlin,       &! Tangent input : after salinity
         & zta_tlin,       &! Tangent input : after temperature
         & zqns_tlin,      &! Tangent input : solar radiation (w/m2)
         & zqsr_tlin,      &! Tangent input : evaporation - precipitation (free surface)
         & zemps_tlin,     &! Tangent input : evaporation - precipitation (free surface)
         & zsbc_tc_tlin,   &! Tangent input : evaporation - precipitation (free surface)
         & zsbc_sc_tlin,   &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_tc_tlin,   &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_sc_tlin,   &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_tc_b_tlin, &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_sc_b_tlin, &! Tangent input : evaporation - precipitation (free surface)
         & zsa_tlout,      &! Tangent output: after salinity
         & zta_tlout,      &! Tangent output: after temperature
         & zqns_tlout,     &! Tangent output: after temperature
         & zqsr_tlout,     &! Tangent output: after temperature
         & zsbc_tc_tlout,  &! Tangent output: after temperature
         & zsbc_sc_tlout,  &! Tangent output: after temperature
         & zsbc_tc_b_tlout,&! Tangent output: after temperature
         & zsbc_sc_b_tlout,&! Tangent output: after temperature
         & zsa_adin,       &! Adjoint input : after salinity
         & zta_adin,       &! Adjoint input : after temperature
         & zqns_adin,     &! Tangent output: after temperature
         & zqsr_adin,     &! Tangent output: after temperature
         & zsbc_tc_adin,  &! Tangent output: after temperature
         & zsbc_sc_adin,  &! Tangent output: after temperature
         & zsbc_tc_b_adin,&! Tangent output: after temperature
         & zsbc_sc_b_adin,&! Tangent output: after temperature
         & zsn_adout,      &! Adjoint output: now salinity
         & zsa_adout,      &! Adjoint output: after salinity
         & zta_adout,      &! Adjoint output: after temperature
         & zqns_adout,     &! Adjoint output: solar radiation (w/m2)
         & zqsr_adout,      &! Tangent input : evaporation - precipitation (free surface)
         & zemps_adout,     &! Tangent input : evaporation - precipitation (free surface)
         & zsbc_tc_adout,   &! Tangent input : evaporation - precipitation (free surface)
         & zsbc_sc_adout,   &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_tc_adout,   &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_sc_adout,   &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_tc_b_adout, &! Tangent input : evaporation - precipitation (free surface)
         & zrnf_sc_b_adout, &! Tangent input : evaporation - precipitation (free surface)
         & zsn,            &! temporary now salinity
         & zsa,            &! temporary after salinity
         & zta,            &! temporary after temperature
         & zqns,           &! temporary solar radiation (w/m2)
         & zemps          ! temporary evaporation - precipitation (free surface)
      REAL(KIND=wp) ::       &
         & zsp1,             & ! scalar product involving the tangent routine
         & zsp1_1,           & ! scalar product involving the tangent routine
         & zsp1_2,           & ! scalar product involving the tangent routine
         & zsp1_3,           & ! scalar product involving the tangent routine
         & zsp1_4,           & ! scalar product involving the tangent routine
         & zsp1_5,           & ! scalar product involving the tangent routine
         & zsp1_6,           & ! scalar product involving the tangent routine
         & zsp1_7,           & ! scalar product involving the tangent routine
         & zsp1_8,           & ! scalar product involving the tangent routine
         & zsp2,             & ! scalar product involving the adjoint routine
         & zsp2_1,           & ! scalar product involving the adjoint routine
         & zsp2_2,           & ! scalar product involving the adjoint routine
         & zsp2_3,           & ! scalar product involving the adjoint routine
         & zsp2_4,           & ! scalar product involving the adjoint routine
         & zsp2_5,           & ! scalar product involving the adjoint routine
         & zsp2_6,           & ! scalar product involving the adjoint routine
         & zsp2_7,           & ! scalar product involving the adjoint routine
         & zsp2_8,           & ! scalar product involving the adjoint routine
         & zsp2_9,           & ! scalar product involving the adjoint routine
         & zsp2_10,           & ! scalar product involving the adjoint routine
         & zsp2_11,           & ! scalar product involving the adjoint routine
         & zsp2_12,           & ! scalar product involving the adjoint routine
         & z2dt,             & ! temporary scalars
         & zraur
      CHARACTER (LEN=14) :: &
         & cl_name

      ALLOCATE( &
         & zsn_tlin(jpi,jpj),     &
         & zsa_tlin(jpi,jpj),     &
         & zta_tlin(jpi,jpj),     &
         & zqns_tlin(jpi,jpj),     &
         & zqsr_tlin(jpi,jpj),     &
         & zemps_tlin(jpi,jpj),     &
         & zsbc_tc_tlin(jpi,jpj),     &
         & zsbc_sc_tlin(jpi,jpj),     &
         & zrnf_tc_tlin(jpi,jpj),     &
         & zrnf_sc_tlin(jpi,jpj),     &
         & zrnf_tc_b_tlin(jpi,jpj),     &
         & zrnf_sc_b_tlin(jpi,jpj),     &
         & zsa_tlout(jpi,jpj),    &
         & zta_tlout(jpi,jpj),    &
         & zqns_tlout(jpi,jpj),    &
         & zqsr_tlout(jpi,jpj),    &
         & zsbc_tc_tlout(jpi,jpj),    &
         & zsbc_sc_tlout(jpi,jpj),    &
         & zsbc_tc_b_tlout(jpi,jpj),    &
         & zsbc_sc_b_tlout(jpi,jpj),    &
         & zsn_adout(jpi,jpj),     &
         & zsa_adout(jpi,jpj),     &
         & zta_adout(jpi,jpj),     &
         & zqns_adout(jpi,jpj),     &
         & zqsr_adout(jpi,jpj),     &
         & zemps_adout(jpi,jpj),     &
         & zsbc_tc_adout(jpi,jpj),     &
         & zsbc_sc_adout(jpi,jpj),     &
         & zrnf_tc_adout(jpi,jpj),     &
         & zrnf_sc_adout(jpi,jpj),     &
         & zrnf_tc_b_adout(jpi,jpj),     &
         & zrnf_sc_b_adout(jpi,jpj),     &
         & zsa_adin(jpi,jpj),    &
         & zta_adin(jpi,jpj),    &
         & zqns_adin(jpi,jpj),    &
         & zqsr_adin(jpi,jpj),    &
         & zsbc_tc_adin(jpi,jpj),    &
         & zsbc_sc_adin(jpi,jpj),    &
         & zsbc_tc_b_adin(jpi,jpj),    &
         & zsbc_sc_b_adin(jpi,jpj),    &
         & zsn(jpi,jpj),          &
         & zsa(jpi,jpj),          &
         & zta(jpi,jpj),          &
         & zqns(jpi,jpj),         &
         & zemps(jpi,jpj)         &
         & )


      ! Initialize constants
      z2dt  = 2.0_wp * rdt       ! time step: leap-frog
      zraur = 1.0_wp / rau0      ! inverse density of pure water (m3/kg)

      ! Initialize the reference state

      !===========================================================================
      ! 1) dx = ( qns_tl, sn_tl, emps_tl, ta_tl, sa_tl ) and dy = ( ta_tl, sa_tl )
      !===========================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      tsn_tl  (:,:,:,:) = 0.0_wp
      tsa_tl  (:,:,:,:) = 0.0_wp
      emps_tl(:,:)   = 0.0_wp
      qns_tl (:,:)   = 0.0_wp
      qsr_tl (:,:)   = 0.0_wp
      sbc_tsc_tl (:,:,:)   = 0.0_wp
      sbc_tsc_b_tl (:,:,:)   = 0.0_wp
      rnf_tsc_tl (:,:,:)   = 0.0_wp
      rnf_tsc_b_tl (:,:,:)   = 0.0_wp
      tsn_ad  (:,:,:,:) = 0.0_wp
      tsa_ad  (:,:,:,:) = 0.0_wp
      emps_ad(:,:)   = 0.0_wp
      qns_ad (:,:)   = 0.0_wp
      qsr_ad (:,:)   = 0.0_wp
      sbc_tsc_ad (:,:,:)   = 0.0_wp
      sbc_tsc_b_ad (:,:,:)   = 0.0_wp
      rnf_tsc_ad (:,:,:)   = 0.0_wp
      rnf_tsc_b_ad (:,:,:)   = 0.0_wp

      zsn_tlin   (:,:) = 0.0_wp
      zsa_tlin   (:,:) = 0.0_wp
      zta_tlin   (:,:) = 0.0_wp
      zqsr_tlin   (:,:) = 0.0_wp
      zqns_tlin   (:,:) = 0.0_wp
      zemps_tlin   (:,:) = 0.0_wp
      zsbc_tc_tlin   (:,:) = 0.0_wp
      zsbc_sc_tlin   (:,:) = 0.0_wp
      zrnf_tc_tlin   (:,:) = 0.0_wp
      zrnf_sc_tlin   (:,:) = 0.0_wp
      zrnf_tc_b_tlin   (:,:) = 0.0_wp
      zrnf_sc_b_tlin   (:,:) = 0.0_wp
      zsa_tlout  (:,:) = 0.0_wp
      zta_tlout  (:,:) = 0.0_wp
      zqns_tlout  (:,:) = 0.0_wp
      zqsr_tlout  (:,:) = 0.0_wp
      zsbc_tc_tlout  (:,:) = 0.0_wp
      zsbc_sc_tlout  (:,:) = 0.0_wp
      zsbc_tc_b_tlout  (:,:) = 0.0_wp
      zsbc_sc_b_tlout  (:,:) = 0.0_wp
      zsn_adout   (:,:) = 0.0_wp
      zsa_adout   (:,:) = 0.0_wp
      zta_adout   (:,:) = 0.0_wp
      zqsr_adout   (:,:) = 0.0_wp
      zqns_adout   (:,:) = 0.0_wp
      zemps_adout   (:,:) = 0.0_wp
      zsbc_tc_adout   (:,:) = 0.0_wp
      zsbc_sc_adout   (:,:) = 0.0_wp
      zrnf_tc_adout   (:,:) = 0.0_wp
      zrnf_sc_adout   (:,:) = 0.0_wp
      zrnf_tc_b_adout   (:,:) = 0.0_wp
      zrnf_sc_b_adout   (:,:) = 0.0_wp
      zsa_adin  (:,:) = 0.0_wp
      zta_adin  (:,:) = 0.0_wp
      zqns_adin  (:,:) = 0.0_wp
      zqsr_adin  (:,:) = 0.0_wp
      zsbc_tc_adin  (:,:) = 0.0_wp
      zsbc_sc_adin  (:,:) = 0.0_wp
      zsbc_tc_b_adin  (:,:) = 0.0_wp
      zsbc_sc_b_adin  (:,:) = 0.0_wp

      CALL grid_random(  zemps, 'T', 0.0_wp, stdemp )
      CALL grid_random(  zqns, 'T', 0.0_wp, stdqns )
      CALL grid_random(  zsn, 'T', 0.0_wp, stds )
      CALL grid_random(  zsa, 'T', 0.0_wp, stds )
      CALL grid_random(  zta, 'T', 0.0_wp, stdt )

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsn_tlin  (ji,jj) = zsn  (ji,jj)
            zsa_tlin  (ji,jj) = zsa  (ji,jj)
            zta_tlin  (ji,jj) = zta  (ji,jj)
            zemps_tlin(ji,jj) = zemps(ji,jj) / ( z2dt * zraur )
            zqns_tlin (ji,jj) = zqns (ji,jj)
            zqsr_tlin (ji,jj) = zqns (ji,jj)
            zsbc_tc_tlin (ji,jj) = zqns (ji,jj)
            zsbc_sc_tlin (ji,jj) = zqns (ji,jj)
            zrnf_tc_tlin (ji,jj) = zqns (ji,jj)
            zrnf_sc_tlin (ji,jj) = zqns (ji,jj)
            zrnf_tc_b_tlin (ji,jj) = zqns (ji,jj)
            zrnf_sc_b_tlin (ji,jj) = zqns (ji,jj)
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      tsn_tl  (:,:,1,jp_sal) = zsn_tlin  (:,:)
      tsa_tl  (:,:,1,jp_sal) = zsa_tlin  (:,:)
      tsa_tl  (:,:,1,jp_tem) = zta_tlin  (:,:)
      emps_tl(:,:)   = zemps_tlin(:,:)
      qns_tl (:,:)   = zqns_tlin (:,:)
      qsr_tl (:,:)   = zqsr_tlin (:,:)
      sbc_tsc_tl (:,:,jp_tem)   = zsbc_tc_tlin (:,:)
      sbc_tsc_tl (:,:,jp_sal)   = zsbc_sc_tlin (:,:)
      rnf_tsc_tl (:,:,jp_tem)   = zrnf_tc_tlin (:,:)
      rnf_tsc_tl (:,:,jp_sal)   = zrnf_sc_tlin (:,:)
      rnf_tsc_b_tl (:,:,jp_tem)   = zrnf_tc_b_tlin (:,:)
      rnf_tsc_b_tl (:,:,jp_sal)   = zrnf_sc_b_tlin (:,:)

      CALL tra_sbc_tan( nit000 )

      zsa_tlout(:,:) = tsa_tl(:,:,1,jp_sal)
      zta_tlout(:,:) = tsa_tl(:,:,1,jp_tem)
      zqns_tlout(:,:) = qns_tl(:,:)
      zqsr_tlout(:,:) = qsr_tl(:,:)
      zsbc_tc_tlout(:,:) = sbc_tsc_tl(:,:,jp_tem)
      zsbc_sc_tlout(:,:) = sbc_tsc_tl(:,:,jp_sal)
      zsbc_tc_b_tlout(:,:) = sbc_tsc_b_tl(:,:,jp_tem)
      zsbc_sc_b_tlout(:,:) = sbc_tsc_b_tl(:,:,jp_sal)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zsa_adin(ji,jj) = zsa_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_s(1)
            zta_adin(ji,jj) = zta_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
            zqns_adin(ji,jj) = zqns_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
            zqsr_adin(ji,jj) = zqsr_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
            zsbc_tc_adin(ji,jj) = zsbc_tc_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
            zsbc_sc_adin(ji,jj) = zsbc_sc_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
            zsbc_tc_b_adin(ji,jj) = zsbc_tc_b_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
            zsbc_sc_b_adin(ji,jj) = zsbc_sc_b_tlout(ji,jj) &
               &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
               &               * tmask(ji,jj,1) * wesp_t(1)
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1_1 = DOT_PRODUCT( zsa_tlout, zsa_adin )
      zsp1_2 = DOT_PRODUCT( zta_tlout, zta_adin )
      zsp1_3 = DOT_PRODUCT( zqns_tlout, zqns_adin )
      zsp1_4 = DOT_PRODUCT( zqsr_tlout, zqsr_adin )
      zsp1_5 = DOT_PRODUCT( zsbc_tc_tlout, zsbc_tc_adin )
      zsp1_6 = DOT_PRODUCT( zsbc_sc_tlout, zsbc_sc_adin )
      zsp1_7 = DOT_PRODUCT( zsbc_tc_b_tlout, zsbc_tc_b_adin )
      zsp1_8 = DOT_PRODUCT( zsbc_sc_b_tlout, zsbc_sc_b_adin )
      zsp1   = zsp1_1 + zsp1_2 + zsp1_3 + zsp1_4 + zsp1_5 + zsp1_6 + zsp1_7 + zsp1_8

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      tsa_ad(:,:,1,jp_sal) = zsa_adin(:,:)
      tsa_ad(:,:,1,jp_tem) = zta_adin(:,:)
      qns_ad(:,:) = zqns_adin(:,:)
      qsr_ad(:,:) = zqsr_adin(:,:)
      sbc_tsc_ad(:,:,jp_tem) = zsbc_tc_adin(:,:)
      sbc_tsc_ad(:,:,jp_sal) = zsbc_sc_adin(:,:)
      sbc_tsc_b_ad(:,:,jp_tem) = zsbc_tc_b_adin(:,:)
      sbc_tsc_b_ad(:,:,jp_sal) = zsbc_sc_b_adin(:,:)


      CALL tra_sbc_adj( nit000 )

      zsn_adout  (:,:) = tsn_ad(:,:,1,jp_sal)
      zsa_adout  (:,:) = tsa_ad(:,:,1,jp_sal)
      zta_adout  (:,:) = tsa_ad(:,:,1,jp_tem)
      zqns_adout (:,:) = qns_ad(:,: )
      zqsr_adout (:,:) = qsr_ad(:,: )
      zemps_adout(:,:) = emps_ad(:,:)
      zsbc_tc_adout(:,:) = sbc_tsc_ad(:,:,jp_tem)
      zsbc_sc_adout(:,:) = sbc_tsc_ad(:,:,jp_sal)
      zrnf_tc_adout(:,:) = rnf_tsc_ad(:,:,jp_tem)
      zrnf_sc_adout(:,:) = rnf_tsc_ad(:,:,jp_sal)
      zrnf_tc_b_adout(:,:) = rnf_tsc_b_ad(:,:,jp_tem)
      zrnf_sc_b_adout(:,:) = rnf_tsc_b_ad(:,:,jp_sal)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1 = DOT_PRODUCT( zsn_tlin  , zsn_adout   )
      zsp2_2 = DOT_PRODUCT( zsa_tlin  , zsa_adout   )
      zsp2_3 = DOT_PRODUCT( zta_tlin  , zta_adout   )
      zsp2_4 = DOT_PRODUCT( zqns_tlin , zqns_adout  )
      zsp2_5 = DOT_PRODUCT( zemps_tlin, zemps_adout )
      zsp2_6 = DOT_PRODUCT( zqsr_tlin, zqsr_adout )
      zsp2_7 = DOT_PRODUCT( zsbc_tc_tlin, zsbc_tc_adout )
      zsp2_8 = DOT_PRODUCT( zsbc_sc_tlin, zsbc_sc_adout )
      zsp2_9 = DOT_PRODUCT( zrnf_tc_tlin, zrnf_tc_adout )
      zsp2_10 = DOT_PRODUCT( zrnf_sc_tlin, zrnf_sc_adout )
      zsp2_11 = DOT_PRODUCT( zrnf_tc_b_tlin, zrnf_tc_b_adout )
      zsp2_12 = DOT_PRODUCT( zrnf_sc_b_tlin, zrnf_sc_b_adout )

      zsp2 = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5 + zsp2_6 + zsp2_7 + zsp2_8 + zsp2_9 + zsp2_10 + zsp2_11 + zsp2_12

      ! Compare the scalar products

      ! 14 char:'12345678901234'
      cl_name = 'tra_sbc_adj   '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      DEALLOCATE( &
         & zsn_tlin,     &
         & zsa_tlin,     &
         & zta_tlin,     &
         & zqns_tlin,     &
         & zqsr_tlin,     &
         & zemps_tlin,     &
         & zsbc_tc_tlin,     &
         & zsbc_sc_tlin,     &
         & zrnf_tc_tlin,     &
         & zrnf_sc_tlin,     &
         & zrnf_tc_b_tlin,     &
         & zrnf_sc_b_tlin,     &
         & zsa_tlout,    &
         & zta_tlout,    &
         & zqns_tlout,    &
         & zqsr_tlout,    &
         & zsbc_tc_tlout,    &
         & zsbc_sc_tlout,    &
         & zsbc_tc_b_tlout,    &
         & zsbc_sc_b_tlout,    &
         & zsn_adout,     &
         & zsa_adout,     &
         & zta_adout,     &
         & zqns_adout,     &
         & zqsr_adout,     &
         & zemps_adout,     &
         & zsbc_tc_adout,     &
         & zsbc_sc_adout,     &
         & zrnf_tc_adout,     &
         & zrnf_sc_adout,     &
         & zrnf_tc_b_adout,     &
         & zrnf_sc_b_adout,     &
         & zsa_adin,    &
         & zta_adin,    &
         & zqns_adin,    &
         & zqsr_adin,    &
         & zsbc_tc_adin,    &
         & zsbc_sc_adin,    &
         & zsbc_tc_b_adin,    &
         & zsbc_sc_b_adin,    &
         & zsn,          &
         & zsa,          &
         & zta,          &
         & zqns,         &
         & zemps         &
         & )
   END SUBROUTINE tra_sbc_adj_tst
   !!======================================================================
END MODULE trasbc_tam
