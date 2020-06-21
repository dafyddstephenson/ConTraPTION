MODULE pttam
#ifdef key_tam
!2017-02-10/16 separated TAM  passive tracer time step loops from nemogcm_tam.F90

  USE par_oce
  USE domain
  USE oce
  USE oce_tam
  USE sbc_oce_tam
  USE sbc_oce, ONLY:rnf_b !!! 2017-03-24 fix to avoid problems with uninitialised field rnf_b
  USE sol_oce_tam
  USE tamctl
  USE iom
  USE trj_tam
  USE wrk_nemo
  USE step_tam, ONLY: stp_tan, stp_adj
  USE step_oce_tam

  USE sbcssr_tam, ONLY: qrp_tl, erp_tl, qrp_ad, erp_ad, nn_sstr
  USE sbcfwb_tam, ONLY: a_fwb_tl, a_fwb_ad

  USE sbcmod_tam, ONLY: sbc_tmp_rm

#  include "domzgr_substitute.h90"

  IMPLICIT NONE

  PRIVATE

  ! Namelisted variables
  CHARACTER(len=132),PUBLIC :: cn_pttam_init
  INTEGER                   :: nn_pttam_out_freq

  !!!! 2016-06-20: added adjoint time stepping loop
  INTEGER :: jk
  
  ! Variable declaration
  REAL(KIND=wp), POINTER, DIMENSION(:,:,:) :: ztn_tlin
  INTEGER:: ncid


  REAL(wp), PUBLIC, ALLOCATABLE, SAVE, DIMENSION(:,:):: tmp_rm,sal_rm !: 2017-03-24  mean tracer volume removed in top layer
  PUBLIC pt_init 
  PUBLIC pt_finalise
  PUBLIC pt_run
  PUBLIC pt_tam_wri !2017-03-03 added subroutine to write output variables to file 

CONTAINS 

  SUBROUTINE pt_init

    NAMELIST/nampttam/cn_pttam_init, nn_pttam_out_freq

    cn_pttam_init = 'PTTAM_init.nc'
    nn_pttam_out_freq = 15

    REWIND(numnam)
    READ(numnam, nampttam)
    IF (lwp) THEN
       WRITE(numout,*) "pttam - filename for initialisation field: ", TRIM(cn_pttam_init)
       WRITE(numout,*) "pttam - 1/(output frequency):              ", nn_pttam_out_freq
    END IF

    IF (.NOT.ALLOCATED(tmp_rm))    ALLOCATE(tmp_rm(jpi,jpj), sal_rm(jpi,jpj))
    
    tmp_rm(:,:) = 0._wp
    sal_rm(:,:) = 0._wp !2017-03-24 storage surface removal of tracer by damping

    !!! 2017-03-24 fix to avoid problems with uninitialised field rnf_b
                  ! force flux boundary condition
    rnf_b(:,:) = 0.0_wp

  END SUBROUTINE pt_init

  SUBROUTINE pt_finalise
 
     IF(ALLOCATED(tmp_rm)) DEALLOCATE(tmp_rm,sal_rm)

  END SUBROUTINE pt_finalise

  SUBROUTINE pt_run(pn_swi)

    INTEGER::pn_swi

    IF (pn_swi == 200) THEN
       CALL pt_tan
    ELSEIF (pn_swi == 201) THEN
       CALL pt_adj
    END IF
    
  END SUBROUTINE pt_run

  SUBROUTINE pt_tan

    INTEGER::istep

    ! Initialisation as in test subroutine of OPATAM_SRC/step.F90
    CALL     oce_tam_init(1)
    CALL sbc_oce_tam_init(1)
    CALL sol_oce_tam_init(1)
#if defined key_tradmp
    CALL trc_oce_tam_init(1)
    strdmp_tl = 0.0_wp
    ttrdmp_tl = 0.0_wp
#endif
    qrp_tl = 0.0_wp
    erp_tl = 0.0_wp   
    emp_tl(:,:) = 0.0_wp
    a_fwb_tl = 0.0_wp

    ! Variable allocation and initialisation
    CALL wrk_alloc(jpi,jpj,jpk,ztn_tlin)
    ztn_tlin(:,:,:) = 0.0_wp

    ! Reading in of initial perturbation
    !!!! 2016-05-16 added variable for name of initial tracer distribution
    CALL iom_open(cn_pttam_init,ncid,kiolib = jpnf90)
    CALL iom_get(ncid,jpdom_autoglo,"Tinit",ztn_tlin,0)

    ztn_tlin(:,:,:) = ztn_tlin(:,:,:)*tmask(:,:,:)

#if defined key_mpp_mpi
    CALL lbc_lnk(ztn_tlin(:,:,:), 'T', 1.0_wp)
#endif

    ! Initialisation of TL model
    istep = nit000 - 1
    CALL trj_rea( istep, 1)
    istep = nit000

    CALL day_tam(nit000, 0)
    un_tl(:,:,:) = 0.0_wp
    vn_tl(:,:,:) = 0.0_wp
    sshn_tl(:,:) = 0.0_wp
    
    tsn_tl(:,:,:,jp_tem) = ztn_tlin(:,:,:)
    tsn_tl(:,:,:,jp_sal) = ztn_tlin(:,:,:)

    CALL iom_close(ncid)

    ub_tl(:,:,:) = 0.0_wp
    vb_tl(:,:,:) = 0.0_wp
    sshb_tl(:,:) = 0.0_wp
    tsb_tl(:,:,:,jp_tem) = ztn_tlin(:,:,:)
    tsb_tl(:,:,:,jp_sal) = ztn_tlin(:,:,:)

    CALL pt_tam_wri(nit000 - 1,0) !write to output file for initial step
    CALL istate_init_tan
    ! Time step loop
    DO istep = nit000, nitend, 1
       CALL stp_tan( istep )
       IF (nn_sstr == 1) THEN
          tmp_rm(:,:) = tmp_rm(:,:) - qrp_tl(:,:)*rdttra(1)*ro0cpr*e1t(:,:)*e2t(:,:)
       END IF
       
       un_tl(:,:,:) = 0.0_wp
       vn_tl(:,:,:) = 0.0_wp
       sshn_tl(:,:) = 0.0_wp 
       ub_tl(:,:,:) = 0.0_wp
       vb_tl(:,:,:) = 0.0_wp
       sshb_tl(:,:) = 0.0_wp

       ! write output ocasionally...

       IF (MOD(istep - nit000 + 1, nn_pttam_out_freq) ==0) THEN
          CALL pt_tam_wri( istep,0 ) !2017-03-03 added output writing to separate subroutine
       END IF

    END DO
    ! Variable deallocation
    CALL wrk_dealloc(jpi, jpj, jpk, ztn_tlin)

    IF (lwp) THEN
       WRITE(numout,*)
       WRITE(numout,*) ' TL_PASSIVE: Finished!'
       WRITE(numout,*) ' ---------------------'
       WRITE(numout,*)
    ENDIF
    CALL flush(numout)

    ! 2016-06-20: added adjoint time stepping loop
  END SUBROUTINE pt_tan

SUBROUTINE pt_adj

  INTEGER::istep
  
  CALL trj_rea(nit000-1, 1)
  DO istep = nit000, nitend
     CALL day_tam(istep, 0)
  END DO
  CALL trj_rea(istep-1, -1)

  call oce_tam_init(2)
  call sbc_oce_tam_init(2)
  call sol_oce_tam_init(2)
  call trc_oce_tam_init(2)
#if defined key_tradmp
  strdmp_ad       = 0.0_wp
  ttrdmp_ad       = 0.0_wp
#endif
  qrp_ad          = 0.0_wp
  erp_ad          = 0.0_wp
  emp_ad(:,:)     = 0.0_wp
  a_fwb_ad        = 0.0_wp

!####################################################################################################
!!! Read in passive tracer initial dye injection and initialise other variables

! Variable allocation and initialisation
  CALL wrk_alloc(jpi,jpj,jpk,ztn_tlin)
  ztn_tlin(:,:,:) = 0.0_wp
  ! Reading in of initial perturbation
  CALL iom_open(cn_pttam_init,ncid,kiolib = jpnf90)
  CALL iom_get(ncid,jpdom_autoglo,"Tinit",ztn_tlin,0)

#if defined key_mpp_mpi
  CALL lbc_lnk(ztn_tlin(:,:,:), 'T', 1.0_wp)
#endif

  tsn_ad(:,:,:,:) = 0.0_wp
     DO jk=1,jpk
        tsn_ad(:,:,jk,jp_tem) = 1.0_wp*tmsk_i(:,:,jk)*ztn_tlin(:,:,jk)*e1t(:,:)*e2t(:,:)*fse3t(:,:,jk)
        tsn_ad(:,:,jk,jp_sal) = 1.0_wp*tmsk_i(:,:,jk)*ztn_tlin(:,:,jk)*e1t(:,:)*e2t(:,:)*fse3t(:,:,jk)
     END DO
     tsn_ad(:,:,1,jp_tem) = tsn_ad(:,:,1,jp_tem)  +  1.0_wp*tmsk_i(:,:,1)*ztn_tlin(:,:,1)*e1t(:,:)*e2t(:,:)*sshn(:,:)
     tsn_ad(:,:,1,jp_sal) = tsn_ad(:,:,1,jp_sal)  +  1.0_wp*tmsk_i(:,:,1)*ztn_tlin(:,:,1)*e1t(:,:)*e2t(:,:)*sshn(:,:)

  un_ad(:,:,:)    = 0.0_wp
  vn_ad(:,:,:)    = 0.0_wp
  sshn_ad(:,:)    = 0.0_wp

  tsb_ad(:,:,:,:) = 0.0_wp
  ub_ad(:,:,:)    = 0.0_wp
  vb_ad(:,:,:)    = 0.0_wp
  sshb_ad(:,:)    = 0.0_wp
!####################################################################################################

  DO istep = nitend, nit000, -1

     IF (MOD(istep - nit000 + 1, nn_pttam_out_freq) ==0) THEN
        CALL pt_tam_wri( istep,1) !2017-03-03 added output writing to separate subroutine
     END IF

     un_ad(:,:,:) = 0.0_wp
     vn_ad(:,:,:) = 0.0_wp
     sshn_ad(:,:) = 0.0_wp 
     ub_ad(:,:,:) = 0.0_wp
     vb_ad(:,:,:) = 0.0_wp
     sshb_ad(:,:) = 0.0_wp

     CALL stp_adj(istep)

     tmp_rm(:,:) = sbc_tmp_rm(:,:)

  END DO

  
  CALL istate_init_adj

 CALL pt_tam_wri(nit000 - 1,1) !write to output file for initial step


END SUBROUTINE pt_adj

SUBROUTINE pt_tam_wri( kstp , wri_swi )


  REAL(wp), POINTER, DIMENSION(:,:,:) :: zicapprox3d 
  INTEGER, INTENT( in ) :: wri_swi
  INTEGER, INTENT( in ) :: kstp
  CHARACTER(LEN=132)::zfname

  CALL wrk_alloc(jpi, jpj, jpk, zicapprox3d)


!pt_conc,pt_vol and pt_vent_vol?
IF (wri_swi==0) THEN

   WRITE(zfname, FMT='(A,I0.8,A)') 'PTTAM_output_', kstp, '.nc'
   
   CALL iom_open(zfname, ncid, ldwrt=.TRUE., kiolib = jprstlib)
   !! link tsn_tl and tsb_tl into one variable of tangent-linear tracer concentration
   zicapprox3d(:,:,:) = tsn_tl(:,:,:,jp_tem) + tsb_tl(:,:,:,jp_tem)
   CALL lbc_lnk(zicapprox3d(:,:,:), 'T', 1.0_wp)
   CALL iom_rstput(kstp, kstp, ncid,    'pt_conc_tl', zicapprox3d)
   !! Output ventilation volume
   CALL iom_rstput(kstp, kstp, ncid, 'pt_vent_tl', tmp_rm(:,:))


ELSEIF (wri_swi==1) THEN
   WRITE(zfname, FMT='(A,I0.8,A)') 'PTTAM_output_', kstp, '.nc'

   CALL iom_open(zfname, ncid, ldwrt=.TRUE., kiolib = jprstlib)
   !! join tsn_ad and tsb_ad into one variable of adjoint tracer volume
   zicapprox3d(:,:,:) = tsn_ad(:,:,:,jp_tem) + tsb_ad(:,:,:,jp_tem)
   CALL lbc_lnk_adj(zicapprox3d(:,:,:), 'T', 1.0_wp)
   CALL iom_rstput(kstp, kstp, ncid,    'pt_vol_ad' , zicapprox3d)
   !! Output ventilation volume
   CALL iom_rstput(kstp, kstp, ncid, 'pt_vent_ad', tmp_rm(:,:))

END IF
   CALL iom_rstput(kstp, kstp, ncid, 'tn'   , tsn(:,:,:,jp_tem)   )
   CALL iom_rstput(kstp, kstp, ncid, 'sn'   , tsn(:,:,:,jp_sal)   )
   CALL iom_close(ncid)
   CALL wrk_dealloc(jpi, jpj, jpk, zicapprox3d) 

END SUBROUTINE pt_tam_wri

#endif
END MODULE pttam
