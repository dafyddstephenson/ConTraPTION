MODULE stpctl_tam
   !!==============================================================================
   !!                       ***  MODULE  stpctl_tam  ***
   !! Ocean run control :  gross check of the ocean time stepping
   !!                      Tangent and adjoint module
   !!==============================================================================
   !! History of the direct module:
   !!            OPA  ! 1991-03  (G. Madec) Original code
   !!            6.0  ! 1992-06  (M. Imbard)
   !!            8.0  ! 1997-06  (A.M. Treguier)
   !!   NEMO     1.0  ! 2002-06  (G. Madec)  F90: Free form and module
   !!            2.0  ! 2009-07  (G. Madec)  Add statistic for time-spliting
   !! History of the T&A module:
   !!            3.0  ! 2009-??  (A. Vidard) orignal version
   !!            3.2  ! 2010-04  (A. Vidard) nemo3.2 update
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   stp_ctl_tan      : Control the run
   !!   stp_ctl_adj      : Control the run
   !!----------------------------------------------------------------------
   !! * Modules used
   USE in_out_manager
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE sol_oce
   USE lib_mpp
   USE dynspg_oce

   IMPLICIT NONE
   PRIVATE

   !! * Accessibility
   PUBLIC stp_ctl_tan           ! routine called by steptan.F90
   PUBLIC stp_ctl_adj           ! routine called by stepadj.F90
   !!----------------------------------------------------------------------

CONTAINS

   SUBROUTINE stp_ctl_tan( kt, kindic, ksign )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl_tan  ***
      !!
      !! ** Purpose of the direct routine:   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Print solver statistics in numsol
      !!              - Stop the run IF problem for the solver ( indec < 0 )
      !!
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      INTEGER, INTENT( in ) ::   ksign
      INTEGER, INTENT( inout ) ::   kindic  ! indicator of solver convergence

      !! * local declarations
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      INTEGER  ::   ii, ij, ik                 ! temporary integers
      REAL(wp) ::   zumax, zvmax, ztmax, zsmax ! temporary scalars
      INTEGER, DIMENSION(3) ::   ilocu         !
      INTEGER, DIMENSION(3) ::   ilocv         !
      INTEGER, DIMENSION(3) ::   iloct         !
      INTEGER, DIMENSION(3) ::   ilocs         !
      CHARACTER(len=80) :: clname
      LOGICAL :: lfirst =.TRUE.
      IF( kt == nit000 .AND. lwp .AND. lfirst ) THEN
         WRITE(numout,*)
         WRITE(numout,*) 'stp_ctl_tan : time-stepping control'
         WRITE(numout,*) '~~~~~~~~~~~'
         ! open time.step file
         clname = 'time_tan.step'
         CALL ctl_opn( numstp, clname, 'REPLACE', 'FORMATTED', 'SEQUENTIAL', -1, numout, lwp, narea )
         lfirst = .FALSE.
      ENDIF

      IF(lwp) WRITE(numstp, '(1x, i8)' ) kt        !* save the current time step in numstp
      IF(lwp) REWIND(numstp)                       !  --------------------------

      ! elliptic solver statistics (if required)
      ! --------------------------
      IF( lk_dynspg_flt ) THEN
      ! Solver
      IF(lwp) WRITE(numsol,9200) kt, niter, res, SQRT(epsr)/eps

      ! Islands (if exist)
!!!!AV TO DO      IF( lk_isl )   CALL isl_stp_ctl_tan( kt, kindic )


      ! Output in numwso and numwvo IF kindic<0
      ! ---------------------------------------
      !    (i.e. problem for the solver)
      IF( kindic < 0 ) THEN
         IF(lwp) THEN
            WRITE(numout,*) ' stpctl_tan: the elliptic solver DO not converge or explode'
            WRITE(numout,*) ' ========= '
            WRITE(numout,9200) kt, niter, res, sqrt(epsr)/eps
            WRITE(numout,*)
           WRITE(numout,*) ' =========  *******************************'
         ENDIF
      ENDIF
      ENDIF

9200  FORMAT(' it :', i8, ' niter :', i4, ' res :',e30.18,' b :',e30.18)

      ! Test maximum of tangent fields
      ! ------------------------

      ! 1. zonal velocity
      ! -----------------

      !! zumax = MAXVAL( ABS( un_tl(:,:,:) ) )   ! slower than the following loop on NEC SX5
      zumax = 0.e0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zumax = MAX(zumax,ABS(un_tl(ji,jj,jk)))
          END DO
        END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( zumax )       ! max over the global domain

      IF( MOD( kt, nwrite ) == 1 .AND. lwp ) WRITE(numout,*) ' ==>> time-step= ',kt,' abs(U) max: ', zumax
      IF( zumax > 50._wp) THEN
         IF( lk_mpp ) THEN
            CALL mpp_maxloc(ABS(un_tl),umask,zumax,ii,ij,ik)
         ELSE
            ilocu = MAXLOC( ABS( un_tl(:,:,:) ) )
            ii = ilocu(1) + nimpp - 1
            ij = ilocu(2) + njmpp - 1
            ik = ilocu(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stpctl_tan: the zonal velocity is larger than 50 m/s'
            WRITE(numout,*) ' ========= '
            WRITE(numout,9400) kt, zumax, ii, ij, ik
            WRITE(numout,*)
         ENDIF
         kindic = -3

      ENDIF
9400  FORMAT (' kt=',i6,' max abs(U): ',1pg11.4,', i j k: ',3i5)
      ! 2. meridional velocity
      ! ----------------------
      zvmax = 0.e0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zvmax = MAX(zvmax,ABS(vn_tl(ji,jj,jk)))
          END DO
        END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( zvmax )   ! max over the global domain

      IF( MOD( kt, nwrite ) == 1 .AND. lwp ) WRITE(numout,*) ' ==>> time-step= ',kt,' abs(V) max: ', zvmax
      IF( zvmax > 50.) THEN
         IF( lk_mpp ) THEN
            CALL mpp_maxloc(ABS(vn_tl),vmask,zvmax,ii,ij,ik)
         ELSE
            ilocv = MAXLOC( ABS( vn_tl(:,:,:) ) )
            ii = ilocv(1) + nimpp - 1
            ij = ilocv(2) + njmpp - 1
            ik = ilocv(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stpctl_tan: the meridional  velocity is larger than 50 m/s'
            WRITE(numout,*) ' ========= '
            WRITE(numout,9410) kt, zvmax, ii, ij, ik
            WRITE(numout,*)
         ENDIF
         kindic = -3
      ENDIF
9410  FORMAT (' kt=',i6,' max abs(V): ',1pg11.4,', i j k: ',3i5)
      ! 3. Temperature
      ! ----------------------
      ztmax = 0.e0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               ztmax = MAX(ztmax,ABS(tsn_tl(ji,jj,jk,jp_tem)))
          END DO
        END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( ztmax )   ! max over the global domain

      IF( MOD( kt, nwrite ) == 1 .AND. lwp ) WRITE(numout,*) ' ==>> time-step= ',kt,' abs(T) max: ', zTmax
      IF( ztmax > 80.) THEN
         IF( lk_mpp ) THEN
            CALL mpp_maxloc(ABS(tsn_tl(:,:,:,jp_tem)),tmask,ztmax,ii,ij,ik)
         ELSE
            iloct = MAXLOC( ABS( tsn_tl(:,:,:,jp_tem) ) )
            ii = iloct(1) + nimpp - 1
            ij = iloct(2) + njmpp - 1
            ik = iloct(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stpctl_tan: the temperature is larger than 80 K'
            WRITE(numout,*) ' ========= '
            WRITE(numout,9420) kt, ztmax, ii, ij, ik
            WRITE(numout,*)
         ENDIF
         kindic = -3
      ENDIF
9420  FORMAT (' kt=',i6,' max abs(T): ',1pg11.4,', i j k: ',3i5)

      ! 3. Salinity
      ! ----------------------
      zsmax = 0.e0
      DO jk = 1, jpk
         DO jj = 1, jpj
            DO ji = 1, jpi
               zsmax = MAX(zsmax,ABS(tsn_tl(ji,jj,jk,jp_sal)))
          END DO
        END DO
      END DO
      IF( lk_mpp )   CALL mpp_max( zsmax )   ! max over the global domain

      IF( MOD( kt, nwrite ) == 1 .AND. lwp ) WRITE(numout,*) ' ==>> time-step= ',kt,' abs(S) max: ', zsmax
      IF( zsmax > 100.) THEN
         IF( lk_mpp ) THEN
            CALL mpp_maxloc(ABS(tsn_tl(:,:,:,jp_sal)),tmask,zsmax,ii,ij,ik)
         ELSE
            ilocs = MAXLOC( ABS( tsn_tl(:,:,:,jp_sal) ) )
            ii = ilocs(1) + nimpp - 1
            ij = ilocs(2) + njmpp - 1
            ik = ilocs(3)
         ENDIF
         IF(lwp) THEN
            WRITE(numout,cform_err)
            WRITE(numout,*) ' stpctl_tan: the Salinity is larger than 100 o/oo'
            WRITE(numout,*) ' ========== '
            WRITE(numout,9430) kt, zsmax, ii, ij, ik
            WRITE(numout,*)
         ENDIF
         kindic = -3
      ENDIF
9430  FORMAT (' kt=',i6,' max abs(S): ',1pg11.4,', i j k: ',3i5)

   END SUBROUTINE stp_ctl_tan
   SUBROUTINE stp_ctl_adj( kt, kindic, ksign )
      !!----------------------------------------------------------------------
      !!                    ***  ROUTINE stp_ctl_adj  ***
      !!
      !! ** Purpose of the direct routine:   Control the run
      !!
      !! ** Method  : - Save the time step in numstp
      !!              - Print it each 50 time steps
      !!              - Print solver statistics in numsol
      !!              - Stop the run IF problem for the solver ( indec < 0 )
      !!
      !! History :
      !!        !  91-03  ()
      !!        !  91-11  (G. Madec)
      !!        !  92-06  (M. Imbard)
      !!        !  97-06  (A.M. Treguier)
      !!   8.5  !  02-06  (G. Madec)  F90: Free form and module
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      INTEGER, INTENT( in ) ::   ksign
      INTEGER, INTENT( inout ) ::   kindic  ! indicator of solver convergence

      !! * local declarations

   END SUBROUTINE stp_ctl_adj

   !!======================================================================
END MODULE stpctl_tam
