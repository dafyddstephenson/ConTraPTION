MODULE sbcfwb_tam
   !!======================================================================
   !!                       ***  MODULE  sbcfwb  ***
   !! Ocean fluxes   : domain averaged freshwater budget
   !!======================================================================
   !! History of the direct module:
   !!            8.2  !  01-02  (E. Durand)  Original code
   !!            8.5  !  02-06  (G. Madec)  F90: Free form and module
   !!            9.0  !  06-08  (G. Madec)  Surface module
   !!            9.2  !  09-07  (C. Talandier) emp mean s spread over erp area
   !! History of the T&A module:
   !!            9.0  !  06-08  (A. Vidard)  Surface module
   !!            9.2  !  09-07  (A. Vidard)  NEMO3.2 update
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   sbc_fwb      : freshwater budget for global ocean configurations
   !!                  in free surface and forced mode
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE oce_tam
   USE sbc_oce_tam
   USE dom_oce
   USE tstool_tam
   USE in_out_manager
   USE lib_mpp                 ! distribued memory computing library
   USE gridrandom
   USE dotprodfld
   USE wrk_nemo
   USE lib_fortran
   USE timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   sbc_fwb_adj      ! routine called by sbcmod_tam
   PUBLIC   sbc_fwb_tan      ! routine called by sbcmod_tam
   PUBLIC   sbc_fwb_adj_tst  ! routine called by tst

   REAL(wp), PUBLIC ::   a_fwb_tl           ! for before year.
   REAL(wp), PUBLIC ::   a_fwb_ad           ! for before year.
   REAL(wp) ::   area               ! global mean ocean surface (interior domain)

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
   !! $Id: sbcfwb.F90 1168 2008-08-11 10:21:06Z rblod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sbc_fwb_tan( kt, kn_fwb, kn_fsbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_fwb_tan  ***
      !!
      !! ** Purpose :   Control the mean sea surface drift
      !!
      !! ** Method  :   several ways  depending on kn_fwb
      !!                =0 no control
      !!                =1 global mean of emp set to zero at each nn_fsbc time step
      !!                =2 annual global mean corrected from previous year
      !!                =3 variable relaxation time to zero or to a time varying field
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( in ) ::   kn_fsbc  !
      INTEGER, INTENT( in ) ::   kn_fwb   ! ocean time-step index
      !!
      INTEGER  ::   inum                  ! temporary logical unit
      INTEGER  ::   ikty, iyear           !
      CHARACTER (len=32) ::   clname
      REAL(wp) ::   z_emptl               ! temporary scalars
      !!----------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_fwb_tan')
      !
      IF( kt == nit000 ) THEN
         !
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_fwb_tan : FreshWater Budget correction.'
            WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
         !
         area = glob_sum( e1e2t(:,:) )           ! interior global domain surface
         !
      ENDIF

      SELECT CASE ( kn_fwb )
      !
      CASE ( 1 )                               ! global mean emp set to zero
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            z_emptl = glob_sum( e1e2t(:,:) * emp_tl(:,:) ) / area
            emp_tl (:,:) = emp_tl (:,:) - z_emptl
            emps_tl(:,:) = emps_tl(:,:) - z_emptl
         ENDIF
         !
      CASE ( 2 )                               ! emp budget adjusted from the previous year
         ! initialisation
         IF( kt == nit000 ) THEN
            a_fwb_tl = 0.0_wp
         ENDIF
         !
         ikty = 365 * 86400 / rdttra(1)    !!bug  use of 365 days leap year or 360d year !!!!!!!
         IF( MOD( kt, ikty ) == 0 ) THEN
            a_fwb_tl   = glob_sum( e1e2t(:,:) * sshn_tl(:,:) )
            a_fwb_tl   = a_fwb_tl * 1.e+3 / ( area * 86400. * 365. )     ! convert in Kg/m3/s = mm/s
!!gm        !                                                      !!bug 365d year
         ENDIF
         !
         ! correct the freshwater fluxes
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            emp_tl (:,:) = emp_tl (:,:) + a_fwb_tl
            emps_tl(:,:) = emps_tl(:,:) + a_fwb_tl
         ENDIF
         !
      CASE ( 3 )
         !
         WRITE(ctmp1,*)'sbc_fwb_tan: nn_fwb=', kn_fwb, ' is not available for TAM yet , choose either 0/1/2'
         CALL ctl_stop( ctmp1 )

         !
      CASE DEFAULT    ! you should never be there
         WRITE(ctmp1,*)'sbc_fwb_tan: nn_fwb=', kn_fwb, ' is not permitted for the FreshWater Budget correction, choose either 0/1/2'
         CALL ctl_stop( ctmp1 )
         !
      END SELECT
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_fwb_tan')
      !
   END SUBROUTINE sbc_fwb_tan


   SUBROUTINE sbc_fwb_adj( kt, kn_fwb, kn_fsbc )
      !!---------------------------------------------------------------------
      !!                  ***  ROUTINE sbc_fwb_adj  ***
      !!
      !! ** Purpose :   Control the mean sea surface drift
      !!
      !! ** Method  :   several ways  depending on kn_fwb
      !!                =0 no control
      !!                =1 global mean of emp set to zero at each nn_fsbc time step
      !!                =2 annual global mean corrected from previous year
      !!                =3 variable relaxation time to zero or to a time varying field
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt       ! ocean time-step index
      INTEGER, INTENT( in ) ::   kn_fsbc  !
      INTEGER, INTENT( in ) ::   kn_fwb   ! ocean time-step index
      !!
      INTEGER  ::   inum                  ! temporary logical unit
      INTEGER  ::   ikty, iyear           !
      INTEGER  ::   ji, jj           !
      CHARACTER (len=32) ::   clname
      REAL(wp) ::   z_empad               ! temporary scalars
      !!----------------------------------------------------------------------
      !
      !
      IF( nn_timing == 1 )  CALL timing_start('sbc_fwb_adj')
      !
      IF( kt == nitend ) THEN
         !
         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) 'sbc_fwb_adj : FreshWater Budget correction.'
            WRITE(numout,*) '~~~~~~~~~~~'
         ENDIF
         !
         area = glob_sum( e1e2t(:,:) )
         !
         ! initialisation
         a_fwb_ad = 0.0_wp
      ENDIF
      z_empad  = 0.0_wp

      SELECT CASE ( kn_fwb )
      !
      CASE ( 1 )                               ! global mean emp set to zero
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            DO jj =  1, jpj
               DO ji =  1, jpi
                  z_empad = z_empad - emps_ad(ji,jj) - emp_ad(ji,jj)
               END DO
            END DO
            IF( lk_mpp )   CALL  mpp_sum( z_empad    )   ! sum over the global domain
            z_empad = z_empad / area
            emp_ad(:,:) = emp_ad(:,:) +  e1e2t(:,:) * z_empad*tmask_i(:,:)
            z_empad = 0._wp
         ENDIF
         !
      CASE ( 2 )                               ! emp budget adjusted from the previous year
         !
         ! initialisation
         ikty = 365 * 86400 / rdttra(1)    !!bug  use of 365 days leap year or 360d year !!!!!!!
         ! correct the freshwater fluxes
         IF( MOD( kt-1, kn_fsbc ) == 0 ) THEN
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  a_fwb_ad = a_fwb_ad + emps_ad(ji,jj) + emp_ad (ji,jj)
               END DO
            END DO
         ENDIF
         !
         IF( MOD( kt, ikty ) == 0 ) THEN
            a_fwb_ad   = a_fwb_ad * 1.e+3 / ( area * 86400. * 365. )     ! convert in Kg/m3/s = mm/s
            IF( lk_mpp )   CALL  mpp_sum( a_fwb_ad )   ! sum over the global domain
            DO jj = 1, jpj
               DO ji = 1, jpi
                  sshn_ad(ji,jj) = sshn_ad(ji,jj) + e1e2t(ji,jj) * a_fwb_ad
               END DO
            END DO
            a_fwb_ad = 0.0_wp
         ENDIF
         !
      CASE ( 3 )
         !
         WRITE(ctmp1,*)'sbc_fwb_tan: nn_fwb=', kn_fwb, ' is not available for TAM yet , choose either 0/1/2'
         CALL ctl_stop( ctmp1 )
         !
      CASE DEFAULT    ! you should never be there
         WRITE(ctmp1,*)'sbc_fwb_adj: nn_fwb=', kn_fwb, ' is not permitted for the FreshWater Budget correction, choose either 0/1/2'
         CALL ctl_stop( ctmp1 )
         !
      END SELECT
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('sbc_fwb_adj')
      !
   END SUBROUTINE sbc_fwb_adj

   SUBROUTINE sbc_fwb_adj_tst ( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE example_adj_tst ***
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

      !! * Local declarations
      INTEGER ::  &
         & istp,  &
         & jstp,  &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & jn_fwb
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zemp_tlin ,     & ! Tangent input
         & zemps_tlin ,    & ! Tangent input
         & zssh_tlin ,     & ! Tangent input
         & zssh_tlout, zssh_adin, &
         & zemp_tlout,     & ! Tangent output
         & zemps_tlout,    & ! Tangent output
         & zemp_adin ,     & ! Adjoint input
         & zemps_adin ,    & ! Adjoint input
         & zemp_adout,     & ! Adjoint output
         & zemps_adout,    & ! Adjoint output
         & zssh_adout,     & ! Adjoint output
         & zr                ! 3D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & zssh_tlout( jpi,jpj),  zssh_adin ( jpi,jpj), &
         & zemp_tlin(  jpi,jpj),     &
         & zemps_tlin( jpi,jpj),     &
         & zssh_tlin(  jpi,jpj),     &
         & zemp_tlout( jpi,jpj),     &
         & zemps_tlout(jpi,jpj),     &
         & zemp_adin(  jpi,jpj),     &
         & zemps_adin( jpi,jpj),     &
         & zemp_adout( jpi,jpj),     &
         & zemps_adout(jpi,jpj),     &
         & zssh_adout( jpi,jpj),     &
         & zr(         jpi,jpj)      &
         & )
      !==================================================================
      ! 1) dx = ( emp_tl, emps_tl, ssh_tl ) and
      !    dy = ( emp_tl, emps_tl )
      !==================================================================
      ! Test for time steps nit000 and nit000 + 1 (the matrix changes)

      DO jstp = nit000, nit000 + 3
         istp = jstp
         IF ( jstp == nit000 + 3 ) istp = nitend

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          zemp_tlin  (:,:) = 0.0_wp
          zemps_tlin (:,:) = 0.0_wp
          zssh_tlin  (:,:) = 0.0_wp
          zemp_tlout (:,:) = 0.0_wp
          zemps_tlout(:,:) = 0.0_wp
          zemp_adin  (:,:) = 0.0_wp
          zemps_adin (:,:) = 0.0_wp
          zemp_adout (:,:) = 0.0_wp
          zemps_adout(:,:) = 0.0_wp
          zssh_adout (:,:) = 0.0_wp
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
      CALL grid_random( zr, 'T', 0.0_wp, stdssh )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zssh_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO

      DO jn_fwb = 1, 2

         a_fwb_tl         = 0.0_wp
         a_fwb_ad         = 0.0_wp
         sshn_ad    (:,:) = 0.0_wp

         sshn_tl(:,:) = zssh_tlin (:,:)
         emps_tl(:,:) = zemps_tlin(:,:)
         emp_tl (:,:) = zemp_tlin (:,:)

         CALL sbc_fwb_tan( istp, jn_fwb, 1 )

         zemps_tlout(:,:) = emps_tl(:,:)
         zemp_tlout (:,:) = emp_tl (:,:)
         zssh_tlout (:,:) = sshn_tl (:,:)
         !-----------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !-----------------------------------------------------------------

         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zemp_adin( ji,jj) = zemp_tlout( ji,jj) &
                  &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
                  &               * tmask(ji,jj,1)
               zemps_adin(ji,jj) = zemps_tlout(ji,jj) &
                  &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) &
                  &               * tmask(ji,jj,1)
               zssh_adin( ji,jj) = zssh_tlout( ji,jj) * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,1) * tmask(ji,jj,1)
            END DO
         END DO

         !-----------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !-----------------------------------------------------------------

         zsp1 = DOT_PRODUCT( zemp_tlout,  zemp_adin  ) &
            & + DOT_PRODUCT( zemps_tlout, zemps_adin )

         zsp1 = zsp1 + DOT_PRODUCT( zssh_tlout, zssh_adin )
         !-----------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !-----------------------------------------------------------------

         emp_ad (:,:) = zemp_adin (:,:)
         emps_ad(:,:) = zemps_adin(:,:)
         sshn_ad(:,:) = zssh_adin(:,:)
         CALL sbc_fwb_adj ( istp, jn_fwb, 1 )

         zemps_adout(:,:) = emps_ad(:,:)
         zemp_adout (:,:) = emp_ad (:,:)
         zssh_adout (:,:) = sshn_ad(:,:)

         zsp2 = DOT_PRODUCT( zemp_tlin,  zemp_adout  ) &
            & + DOT_PRODUCT( zssh_tlin,  zssh_adout  ) &
            & + DOT_PRODUCT( zemps_tlin, zemps_adout )

         ! 14 char:'12345678901234'
         IF ( istp == nit000 ) THEN
             WRITE (cl_name,"(A10,1x,i1,1x,A1)") 'sbcfwb_adj',jn_fwb,'1'
         ELSEIF ( istp == nit000 + 1 ) THEN
             WRITE (cl_name,"(A10,1x,i1,1x,A1)") 'sbcfwb_adj',jn_fwb,'2'
         ELSEIF ( istp == nit000 + 2 ) THEN
             WRITE (cl_name,"(A10,1x,i1,1x,A1)") 'sbcfwb_adj',jn_fwb,'3'
         ELSEIF ( istp == nitend ) THEN
             WRITE (cl_name,"(A10,1x,i1,1x,A1)") 'sbcfwb_adj',jn_fwb,'4'
         END IF
!         WRITE (cl_name,"(A11,2x,i1)") 'sbc_fwb_adj',jn_fwb
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      END DO

      END DO

      DEALLOCATE(       &
         & zemp_tlin,   &
         & zemps_tlin,  &
         & zssh_tlin,   &
         & zemp_tlout,  &
         & zemps_tlout, &
         & zemp_adin,   &
         & zemps_adin,  &
         & zemp_adout,  &
         & zemps_adout, &
         & zssh_adout,  &
         & zr           &
         & )

   END SUBROUTINE sbc_fwb_adj_tst
   !!======================================================================
END MODULE sbcfwb_tam
