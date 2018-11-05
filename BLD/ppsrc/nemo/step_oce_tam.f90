MODULE step_oce_tam
   !!======================================================================
   !!                       ***  MODULE step_oce  ***
   !! Ocean time-stepping : module used in both initialisation phase and time stepping
   !!======================================================================
   !! History :   3.3  ! 2010-08  (C. Ethe)  Original code - reorganisation of the initial phase
   !!----------------------------------------------------------------------
   USE oce              ! ocean dynamics and tracers variables
   USE dom_oce          ! ocean space and time domain variables
   USE par_oce
   USE zdf_oce          ! ocean vertical physics variables
   USE ldftra_oce       ! ocean tracer   - trends
   USE ldfdyn_oce       ! ocean dynamics - trends
   USE in_out_manager   ! I/O manager
   USE iom              !
   USE lbclnk

   USE daymod           ! calendar                         (day     routine)

   USE sbcmod           ! surface boundary condition       (sbc     routine)
   USE sbcrnf           ! surface boundary condition: runoff variables
   USE sbccpl           ! surface boundary condition: coupled formulation (call send at end of step)
   USE cpl_oasis3, ONLY : lk_cpl
   USE sbctide          ! Tide initialisation

   USE traqsr           ! solar radiation penetration      (tra_qsr routine)
   USE trasbc           ! surface boundary condition       (tra_sbc routine)
   USE trabbc           ! bottom boundary condition        (tra_bbc routine)
   USE trabbl           ! bottom boundary layer            (tra_bbl routine)
   USE tradmp           ! internal damping                 (tra_dmp routine)
   USE traadv           ! advection scheme control     (tra_adv_ctl routine)
   USE traldf           ! lateral mixing                   (tra_ldf routine)
   !   zdfkpp           ! KPP non-local tracer fluxes      (tra_kpp routine)
   USE trazdf           ! vertical mixing                  (tra_zdf routine)
   USE tranxt           ! time-stepping                    (tra_nxt routine)
   USE tranpc           ! non-penetrative convection       (tra_npc routine)

   USE eosbn2           ! equation of state                (eos_bn2 routine)

   USE dynadv           ! advection                        (dyn_adv routine)
   USE dynbfr           ! Bottom friction terms            (dyn_bfr routine)
   USE dynvor           ! vorticity term                   (dyn_vor routine)
   USE dynhpg           ! hydrostatic pressure grad.       (dyn_hpg routine)
   USE dynldf           ! lateral momentum diffusion       (dyn_ldf routine)
   USE dynzdf           ! vertical diffusion               (dyn_zdf routine)
   USE dynspg_oce       ! surface pressure gradient        (dyn_spg routine)
   USE dynspg           ! surface pressure gradient        (dyn_spg routine)
   USE dynnept          ! simp. form of Neptune effect(dyn_nept_cor routine)

   USE dynnxt           ! time-stepping                    (dyn_nxt routine)

   USE obc_par          ! open boundary condition variables
   USE obcdta           ! open boundary condition data     (obc_dta routine)
   USE obcrst           ! open boundary cond. restart      (obc_rst routine)
   USE obcrad           ! open boundary cond. radiation    (obc_rad routine)

   USE bdy_par          ! for lk_bdy
   USE bdydta           ! open boundary condition data     (bdy_dta routine)

   USE sshwzv           ! vertical velocity and ssh        (ssh_wzv routine)

   USE ldfslp           ! iso-neutral slopes               (ldf_slp routine)
   USE ldfeiv           ! eddy induced velocity coef.      (ldf_eiv routine)

   USE zdftmx           ! tide-induced vertical mixing     (zdf_tmx routine)
   USE zdfbfr           ! bottom friction                  (zdf_bfr routine)
   USE zdftke           ! TKE vertical mixing              (zdf_tke routine)
   USE zdfgls           ! GLS vertical mixing              (zdf_gls routine)
   USE zdfkpp           ! KPP vertical mixing              (zdf_kpp routine)
   USE zdfddm           ! double diffusion mixing          (zdf_ddm routine)
   USE zdfevd           ! enhanced vertical diffusion      (zdf_evd routine)
   USE zdfric           ! Richardson vertical mixing       (zdf_ric routine)
   USE zdfmxl           ! Mixed-layer depth                (zdf_mxl routine)

   USE zpshde           ! partial step: hor. derivative     (zps_hde routine)

   USE diawri           ! Standard run outputs             (dia_wri routine)
   USE trdicp           ! Ocean momentum/tracers trends    (trd_wri routine)
   USE trdmld           ! mixed-layer trends               (trd_mld routine)
   USE trdmld_rst       ! restart for mixed-layer trends
   USE trdmod_oce       ! ocean momentum/tracers trends
   USE trdmod           ! momentum/tracers trends
   USE trdvor           ! vorticity budget                 (trd_vor routine)
   USE diaptr           ! poleward transports              (dia_ptr routine)
   USE diadct           ! sections transports              (dia_dct routine)
   USE diaar5           ! AR5 diagnosics                   (dia_ar5 routine)
   USE diahth           ! thermocline depth                (dia_hth routine)
   USE diafwb           ! freshwater budget                (dia_fwb routine)
   USE diahsb           ! heat, salt and volume budgets    (dia_hsb routine)
   USE diaharm
   USE flo_oce          ! floats variables
   USE floats           ! floats computation               (flo_stp routine)

   USE asminc           ! assimilation increments      (tra_asm_inc routine)
   !                                                   (dyn_asm_inc routine)

   USE stpctl           ! time stepping control            (stp_ctl routine)
   USE restart          ! ocean restart                    (rst_wri routine)
   USE prtctl           ! Print control                    (prt_ctl routine)

   USE diaobs           ! Observation operator

   USE timing           ! Timing

   USE oce_tam
   USE lbclnk_tam
   USE daymod_tam      ! calendar                         (adjoint of day     routine)
   USE sbc_oce_tam
   USE sbcmod_tam
   USE traqsr_tam      ! solar radiation penetration      (adjoint of tra_qsr routine)
   USE trasbc_tam      ! surface boundary condition       (adjoint of tra_sbc routine)
   USE trabbl_tam      ! bottom boundary layer            (adjoint of tra_bbl routine)
   USE trabbc_tam      ! bottom boundary condition        (adjoint of tra_bbc routine)
   USE tradmp_tam      ! internal damping                 (adjoint of tra_dmp routine)
   USE traadv_tam      ! advection scheme control     (adjoint of tra_adv_ctl routine)
   USE traldf_tam      ! lateral mixing                   (adjoint of tra_ldf routine)
   USE cla_tam         ! cross land advection             (adjoint of tra_cla routine)
   USE trazdf_tam      ! vertical mixing                  (adjoint of tra_zdf routine)
   USE tranxt_tam      ! time-stepping                    (adjoint of tra_nxt routine)
   USE eosbn2_tam      ! equation of state                (adjoint of eos_bn2 routine)
   USE dynadv_tam      ! advection                        (adjoint of dyn_adv routine)
   USE dynvor_tam      ! vorticity term                   (adjoint of dyn_vor routine)
   USE dynhpg_tam      ! hydrostatic pressure grad.       (adjoint of dyn_hpg routine)
   USE dynldf_tam      ! lateral momentum diffusion       (adjoint of dyn_ldf routine)
   USE dynzdf_tam      ! vertical diffusion               (adjoint of dyn_zdf routine)
   USE dynspg_tam      ! surface pressure gradient        (adjoint of dyn_spg routine)
   USE dynnxt_tam      ! time-stepping                    (adjoint of dyn_nxt routine)
   USE dynbfr_tam      ! time-stepping                    (adjoint of dyn_nxt routine)
   USE sshwzv_tam      ! vertical velocity and ssh        (ssh_wzv routine)
   USE divcur_tam      ! hor. divergence and curl      (adjoint of div & cur routines)
   USE cla_tam     ! cross land: hor. divergence      (adjoint of div_cla routine)
   USE zdfbfr_tam
   USE zpshde_tam      ! partial step: hor. derivative     (adjoint of zps_hde routine)
   USE trj_tam
   USE stpctl_tam      ! time stepping control            (adjoint of stp_ctl routine)
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE paresp
   USE istate_tam      !: Initial state setting          (istate_init routine)
   USE sol_oce
   USE sol_oce_tam
   USE trc_oce_tam
   USE sbcrnf_tam
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE step_oce_tam
