MODULE divcur_tam
   !!==============================================================================
   !!       ***  MODULE  divcur_tam : TANGENT/ADJOINT OF MODULE divcur  ***
   !!
   !! Ocean diagnostic variable : horizontal divergence and relative vorticity
   !!
   !! History :  OPA  ! 1987-06  (P. Andrich, D. L Hostis)  Original code
   !!            4.0  ! 1991-11  (G. Madec)
   !!            6.0  ! 1993-03  (M. Guyon)  symetrical conditions
   !!            7.0  ! 1996-01  (G. Madec)  s-coordinates
   !!            8.0  ! 1997-06  (G. Madec)  lateral boundary cond., lbc
   !!            8.1  ! 1997-08  (J.M. Molines)  Open boundaries
   !!            8.2  ! 2000-03  (G. Madec)  no slip accurate
   !!  NEMO      1.0  ! 2002-09  (G. Madec, E. Durand)  Free form, F90
   !!             -   ! 2005-01  (J. Chanut) Unstructured open boundaries
   !!             -   ! 2003-08  (G. Madec)  merged of cur and div, free form, F90
   !!             -   ! 2005-01  (J. Chanut, A. Sellar) unstructured open boundaries
   !!            3.3  ! 2010-09  (D.Storkey and E.O'Dea) bug fixes for BDY module
   !!             -   ! 2010-10  (R. Furner, G. Madec) runoff and cla added directly here
   !!==============================================================================
   !! History of the TAM module:
   !!             7.0  !  95-01  (F. Van den Berghe)
   !!             8.0  !  96-04  (A. Weaver)
   !!             8.1  !  98-02  (A. Weaver)
   !!             8.2  !  00-08  (A. Weaver)
   !!             9.0  !  08-06  (A. Vidard) Skeleton
   !!             9.0  !  08-07  (A. Weaver)
   !!             9.0  !  08-11  (A. Vidard) Nemo v3
   !!             9.0  !  09-02  (A. Vidard) cleanup
   !!             9.0  !  07-12  (P.-A. Bouttier) Nemo v3.4
   !!----------------------------------------------------------------------
   !!   div_cur_tan     : Compute the horizontal divergence and relative
   !!                     vorticity fields (tangent routine)
   !!   div_cur_adj     : Compute the horizontal divergence and relative
   !!                     vorticity fields (adjoint routine)
   !!   div_cur_adj_tst : Test of the adjoint routine
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind
   USE par_oce
   USE in_out_manager
   USE dom_oce
   USE sbc_oce
   USE lbclnk
   USE lbclnk_tam
   USE oce_tam
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing
   USE cla_tam
   USE sbcrnf_tam

   PRIVATE

   !! * Accessibility
   PUBLIC div_cur_tan,    &   ! routine called by steptan.F90
      &   div_cur_adj,    &   ! routine called by stepadj.F90
      &   div_cur_adj_tst     ! adjoint test routine

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

   !!----------------------------------------------------------------------
   !!   Default option                           2nd order centered schemes
   !!----------------------------------------------------------------------
   SUBROUTINE div_cur_tan( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_cur_tan  ***
      !!
      !! ** Purpose of direct routine :
      !!      compute the horizontal divergence and the relative
      !!      vorticity at before and now time-step
      !!
      !! ** Method of direct routine :
      !!              - Divergence:
      !!      - save the divergence computed at the previous time-step
      !!      (note that the Asselin filter has not been applied on hdivb)
      !!      - compute the now divergence given by :
      !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!      Note: if lk_zco=T, e3u=e3v=e3t, they are simplified in the
      !!      above expression
      !!      - apply lateral boundary conditions on hdivn
      !!              - Relavtive Vorticity :
      !!      - save the curl computed at the previous time-step (rotb = rotn)
      !!      (note that the Asselin time filter has not been applied to rotb)
      !!      - compute the now curl in tensorial formalism:
      !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
      !!      - apply lateral boundary conditions on rotn through a call to
      !!      routine lbc_lnk routine.
      !!      Note: Coastal boundary condition: lateral friction set through
      !!      the value of fmask along the coast (see dommsk.F90) and shlat
      !!      (namelist parameter)
      !!
      !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence
      !!              - update rotb , rotn , the before & now rel. vorticity
      !!
      !! History of the direct routine:
      !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
      !!   4.0  !  91-11  (G. Madec)
      !!   6.0  !  93-03  (M. Guyon)  symetrical conditions
      !!   7.0  !  96-01  (G. Madec)  s-coordinates
      !!   8.0  !  97-06  (G. Madec)  lateral boundary cond., lbc
      !!   8.1  !  97-08  (J.M. Molines)  Open boundaries
      !!   9.0  !  02-09  (G. Madec, E. Durand)  Free form, F90
      !!        !  05-01  (J. Chanut) Unstructured open boundaries
      !! History of the TAM routine:
      !!   9.0  !  08-06  (A. Vidard) Skeleton
      !!        !  08-07  (A. Weaver) tangent of the 02-09 version
      !!        !  08-11  (A. Vidard) tangent of the 05-01 version
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) :: &
         & kt     ! ocean time-step index

      !! * Local declarations
      INTEGER :: &
         & ji,  & ! dummy loop indices
         & jj,  &
         & jk
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('div_cur_tan')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cur_tan : horizontal velocity divergence and'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   relative vorticity'
      ENDIF
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         hdivb_tl(:,:,jk) = hdivn_tl(:,:,jk)    ! time swap of div arrays
         rotb_tl (:,:,jk) = rotn_tl (:,:,jk)    ! time swap of rot arrays
         !                                             ! --------
         ! Horizontal divergence                       !   div
         !                                             ! --------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               hdivn_tl(ji,jj,jk) =   &
                  & (   e2u(ji  ,jj  ) * e3u(ji  ,jj  ,jk) * un_tl(ji  ,jj  ,jk) &
                  &   - e2u(ji-1,jj  ) * e3u(ji-1,jj  ,jk) * un_tl(ji-1,jj  ,jk) &
                  &   + e1v(ji  ,jj  ) * e3v(ji  ,jj  ,jk) * vn_tl(ji  ,jj  ,jk) &
                  &   - e1v(ji  ,jj-1) * e3v(ji  ,jj-1,jk) * vn_tl(ji  ,jj-1,jk) &
                  & ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO
         END DO
         !                                             ! --------
         ! relative vorticity                          !   rot
         !                                             ! --------
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               rotn_tl(ji,jj,jk) = (   e2v(ji+1,jj  ) * vn_tl(ji+1,jj  ,jk) &
                  &                  - e2v(ji  ,jj  ) * vn_tl(ji  ,jj  ,jk) &
                  &                  - e1u(ji  ,jj+1) * un_tl(ji  ,jj+1,jk) &
                  &                  + e1u(ji  ,jj  ) * un_tl(ji  ,jj  ,jk) &
                  &                ) * fmask(ji,jj,jk) / ( e1f(ji,jj) * e2f(ji,jj) )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      IF( ln_rnf      )   CALL sbc_rnf_div_tan( hdivn_tl )       ! runoffs (update hdivn field)
      IF( nn_cla == 1 )   CALL cla_div_tan    ( kt )             ! Cross Land Advection (update hdivn field)
      !!
      CALL lbc_lnk( hdivn_tl, 'T', 1. )
      CALL lbc_lnk( rotn_tl , 'F', 1. )     ! lateral boundary cond. (no sign change)
      !
      IF( nn_timing == 1 )  CALL timing_stop('div_cur_tan')
   END SUBROUTINE div_cur_tan

   SUBROUTINE div_cur_adj( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_cur_adj  ***
      !!
      !! ** Purpose of direct routine :
      !!      compute the horizontal divergence and the relative
      !!      vorticity at before and now time-step
      !!
      !! ** Method of direct routine :
      !!              - Divergence:
      !!      - save the divergence computed at the previous time-step
      !!      (note that the Asselin filter has not been applied on hdivb)
      !!      - compute the now divergence given by :
      !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!      Note: if lk_zco=T, e3u=e3v=e3t, they are simplified in the
      !!      above expression
      !!      - apply lateral boundary conditions on hdivn
      !!              - Relavtive Vorticity :
      !!      - save the curl computed at the previous time-step (rotb = rotn)
      !!      (note that the Asselin time filter has not been applied to rotb)
      !!      - compute the now curl in tensorial formalism:
      !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
      !!      - apply lateral boundary conditions on rotn through a call to
      !!      routine lbc_lnk routine.
      !!      Note: Coastal boundary condition: lateral friction set through
      !!      the value of fmask along the coast (see dommsk.F90) and shlat
      !!      (namelist parameter)
      !!
      !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence
      !!              - update rotb , rotn , the before & now rel. vorticity
      !!
      !! History of the direct routine:
      !!   1.0  !  87-06  (P. Andrich, D. L Hostis)  Original code
      !!   4.0  !  91-11  (G. Madec)
      !!   6.0  !  93-03  (M. Guyon)  symetrical conditions
      !!   7.0  !  96-01  (G. Madec)  s-coordinates
      !!   8.0  !  97-06  (G. Madec)  lateral boundary cond., lbc
      !!   8.1  !  97-08  (J.M. Molines)  Open boundaries
      !!   9.0  !  02-09  (G. Madec, E. Durand)  Free form, F90
      !! History of the TAM routine:
      !!   9.0  !  08-06  (A. Vidard) Skeleton
      !!   9.0  !  08-07  (A. Weaver)
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in ) :: &
         & kt     ! ocean time-step index

      !! * Local declarations
      INTEGER :: &
         & ji,  & ! dummy loop indices
         & jj,  &
         & jk
      !!----------------------------------------------------------------------
      !
      if( nn_timing == 1 )  call timing_start('div_cur_adj')
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cur_adj : horizontal velocity divergence and'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   relative vorticity'
      ENDIF
      ! 4. Lateral boundary conditions on hdivn and rotn
      ! ---------------------------------=======---======
      CALL lbc_lnk_adj( rotn_ad , 'F', 1.0_wp )     ! F-point, no sign change
      CALL lbc_lnk_adj( hdivn_ad, 'T', 1.0_wp )     ! T-point, no sign change
      !!
      IF( nn_cla == 1 )   CALL cla_div_adj    ( kt )             ! Cross Land Advection (update hdivn field)
      IF( ln_rnf      )   CALL sbc_rnf_div_adj( hdivn_ad )       ! runoffs (update hdivn field)
      !                                                ! ===============
      DO jk = jpkm1, 1, -1                             ! Horizontal slab
         !                                             ! ===============
         !                                             ! --------
         ! relative vorticity                          !   rot
         !                                             ! --------
         DO jj = jpjm1, 1, -1
            DO ji = jpim1, 1, -1   ! vector opt.
               rotn_ad(ji,jj,jk) =  rotn_ad(ji,jj,jk) * fmask(ji,jj,jk) &
                  &                       / ( e1f(ji,jj) * e2f(ji,jj) )
               un_ad(ji  ,jj  ,jk) = un_ad(ji  ,jj  ,jk) &
                  &                  + e1u(ji  ,jj  ) * rotn_ad(ji,jj,jk)
               un_ad(ji  ,jj+1,jk) = un_ad(ji  ,jj+1,jk) &
                  &                  - e1u(ji  ,jj+1) * rotn_ad(ji,jj,jk)
               vn_ad(ji  ,jj  ,jk) = vn_ad(ji  ,jj  ,jk) &
                  &                  - e2v(ji  ,jj  ) * rotn_ad(ji,jj,jk)
               vn_ad(ji+1,jj  ,jk) = vn_ad(ji+1,jj  ,jk) &
                  &                  + e2v(ji+1,jj  ) * rotn_ad(ji,jj,jk)
               rotn_ad(ji,jj,jk) = 0.0_wp
            END DO
         END DO
         !                                             ! --------
         ! Horizontal divergence                       !   div
         !                                             ! --------
         DO jj = jpjm1, 2, -1
            DO ji = jpim1, 2, -1   ! vector opt.
               hdivn_ad(ji,jj,jk) =  hdivn_ad(ji,jj,jk)  &
                  &                  / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               un_ad(ji  ,jj  ,jk) = un_ad(ji  ,jj  ,jk) &
                  &                  + e2u(ji  ,jj  ) * e3u(ji  ,jj  ,jk) &
                  &                                   * hdivn_ad(ji,jj,jk)
               un_ad(ji-1,jj  ,jk) = un_ad(ji-1,jj  ,jk) &
                  &                  - e2u(ji-1,jj  ) * e3u(ji-1,jj  ,jk) &
                  &                                   * hdivn_ad(ji,jj,jk)
               vn_ad(ji  ,jj  ,jk) = vn_ad(ji  ,jj  ,jk) &
                  &                  + e1v(ji  ,jj  ) * e3v(ji  ,jj  ,jk) &
                  &                                   * hdivn_ad(ji,jj,jk)
               vn_ad(ji  ,jj-1,jk) = vn_ad(ji  ,jj-1,jk) &
                  &                  - e1v(ji  ,jj-1) * e3v(ji  ,jj-1,jk) &
                  &                                   * hdivn_ad(ji,jj,jk)
               hdivn_ad(ji,jj,jk) = 0.0_wp
            END DO
         END DO
         !
         rotn_ad (:,:,jk) = rotn_ad (:,:,jk) + rotb_ad (:,:,jk)  ! time swap
         rotb_ad (:,:,jk) = 0.0_wp
         hdivn_ad(:,:,jk) = hdivn_ad(:,:,jk) + hdivb_ad(:,:,jk)  ! time swap
         hdivb_ad(:,:,jk) = 0.0_wp
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      if( nn_timing == 1 )  call timing_stop('div_cur_adj')
      !
   END SUBROUTINE div_cur_adj


   SUBROUTINE div_cur_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!          ***  ROUTINE div_cur_adj_tst : TEST OF div_cur_adj  ***
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
      !!      1) dx = ( un_tl, vn_tl ) and
      !!         dy = ( hdivn_tl )
      !!      2) dx = ( un_tl, vn_tl ) and
      !!         dy = ( rotntl )
      !!
      !! History :
      !!        ! 08-07 (A. Weaver)
      !!-----------------------------------------------------------------------

      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      INTEGER :: &
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zun_tlin,     & ! Tangent input: now u-velocity
         & zvn_tlin,     & ! Tangent input: now v-velocity
         & zhdivn_tlin,  & ! Tangent input: now horizontal divergence
         & zrotn_tlin,   & ! Tangent input: now relative vorticity
         & zhdivb_tlout, & ! Tangent output: before horizontal divergence
         & zhdivn_tlout, & ! Tangent output: now horizontal divergence
         & zrotb_tlout,  & ! Tangent output: before relative vorticity
         & zrotn_tlout,  & ! Tangent output: now relative vorticity
         & zhdivb_adin,  & ! Adjoint input: before horizontal divergence
         & zhdivn_adin,  & ! Adjoint input: now horizontal divergence
         & zrotb_adin,   & ! Adjoint input: before relative vorticity
         & zrotn_adin,   & ! Adjoint input: now relative vorticity
         & zun_adout,    & ! Adjoint output: now u-velocity
         & zvn_adout,    & ! Adjoint output: now v-velocity
         & zhdivn_adout, & ! Adjoint output: now horizontal divergence
         & zrotn_adout,  & ! Adjoint output: now relative vorticity
         & znu,          & ! 3D random field for u
         & znv             ! 3D random field for v

      REAL(KIND=wp) :: &
                           ! random field standard deviation for:
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp1_3,       & !
         & zsp1_4,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4

      CHARACTER(LEN=14) :: cl_name

      ! Allocate memory

      ALLOCATE( &
         & zun_tlin(jpi,jpj,jpk),     &
         & zvn_tlin(jpi,jpj,jpk),     &
         & zhdivn_tlin(jpi,jpj,jpk),  &
         & zrotn_tlin(jpi,jpj,jpk),   &
         & zhdivb_tlout(jpi,jpj,jpk), &
         & zhdivn_tlout(jpi,jpj,jpk), &
         & zrotb_tlout(jpi,jpj,jpk),  &
         & zrotn_tlout(jpi,jpj,jpk),  &
         & zhdivb_adin(jpi,jpj,jpk),  &
         & zhdivn_adin(jpi,jpj,jpk),  &
         & zrotb_adin(jpi,jpj,jpk),   &
         & zrotn_adin(jpi,jpj,jpk),   &
         & zun_adout(jpi,jpj,jpk),    &
         & zvn_adout(jpi,jpj,jpk),    &
         & zhdivn_adout(jpi,jpj,jpk), &
         & zrotn_adout(jpi,jpj,jpk),  &
         & znu(jpi,jpj,jpk),          &
         & znv(jpi,jpj,jpk)           &
         & )


      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      zun_tlin    (:,:,:) = 0.0_wp
      zvn_tlin    (:,:,:) = 0.0_wp
      zhdivn_tlin (:,:,:) = 0.0_wp
      zrotn_tlin  (:,:,:) = 0.0_wp
      zhdivb_tlout(:,:,:) = 0.0_wp
      zhdivn_tlout(:,:,:) = 0.0_wp
      zrotb_tlout (:,:,:) = 0.0_wp
      zrotn_tlout (:,:,:) = 0.0_wp
      zhdivb_adin (:,:,:) = 0.0_wp
      zhdivn_adin (:,:,:) = 0.0_wp
      zrotb_adin  (:,:,:) = 0.0_wp
      zrotn_adin  (:,:,:) = 0.0_wp
      zrotn_adout (:,:,:) = 0.0_wp
      zhdivn_adout(:,:,:) = 0.0_wp
      zun_adout   (:,:,:) = 0.0_wp
      zvn_adout   (:,:,:) = 0.0_wp

      un_tl   (:,:,:) = 0.0_wp
      vn_tl   (:,:,:) = 0.0_wp
      hdivb_tl(:,:,:) = 0.0_wp
      hdivn_tl(:,:,:) = 0.0_wp
      rotb_tl (:,:,:) = 0.0_wp
      rotn_tl (:,:,:) = 0.0_wp
      hdivb_ad(:,:,:) = 0.0_wp
      hdivn_ad(:,:,:) = 0.0_wp
      rotb_ad (:,:,:) = 0.0_wp
      rotn_ad (:,:,:) = 0.0_wp
      un_ad   (:,:,:) = 0.0_wp
      vn_ad   (:,:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random( znu, 'U', 0.0_wp, stdu )

      CALL grid_random( znv, 'V', 0.0_wp, stdv )

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zun_tlin(ji,jj,jk) = znu(ji,jj,jk)
               zvn_tlin(ji,jj,jk) = znv(ji,jj,jk)
            END DO
         END DO
      END DO

      un_tl(:,:,:) = zun_tlin(:,:,:)
      vn_tl(:,:,:) = zvn_tlin(:,:,:)

      CALL div_cur_tan( nit000 )  ! Generate noise for before hdiv/rot fields

      DO jk = 1, jpk
         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zhdivn_tlin(ji,jj,jk) = 0.5_wp * hdivn_tl(ji,jj,jk)
               zrotn_tlin (ji,jj,jk) = 0.5_wp * rotn_tl (ji,jj,jk)
            END DO
         END DO
      END DO

      un_tl   (:,:,:) = 0.0_wp
      vn_tl   (:,:,:) = 0.0_wp
      hdivb_tl(:,:,:) = 0.0_wp
      hdivn_tl(:,:,:) = 0.0_wp
      rotb_tl (:,:,:) = 0.0_wp
      rotn_tl (:,:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      un_tl   (:,:,:) = zun_tlin   (:,:,:)
      vn_tl   (:,:,:) = zvn_tlin   (:,:,:)
      hdivn_tl(:,:,:) = zhdivn_tlin(:,:,:)

      CALL div_cur_tan( nit000 )

      zhdivb_tlout(:,:,:) = hdivb_tl(:,:,:)
      zhdivn_tlout(:,:,:) = hdivn_tl(:,:,:)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zhdivb_adin(ji,jj,jk) = zhdivb_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk)
              zhdivn_adin(ji,jj,jk) = zhdivn_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1_1 = DOT_PRODUCT( zhdivb_tlout, zhdivb_adin )
      zsp1_2 = DOT_PRODUCT( zhdivn_tlout, zhdivn_adin )
      zsp1   = zsp1_1 + zsp1_2

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      hdivb_ad(:,:,:) = zhdivb_adin(:,:,:)
      hdivn_ad(:,:,:) = zhdivn_adin(:,:,:)
      rotb_ad (:,:,:) = 0.0_wp
      rotn_ad (:,:,:) = 0.0_wp

      CALL div_cur_adj( nit000 )

      zun_adout   (:,:,:) = un_ad   (:,:,:)
      zvn_adout   (:,:,:) = vn_ad   (:,:,:)
      zhdivn_adout(:,:,:) = hdivn_ad(:,:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1 = DOT_PRODUCT( zun_tlin,    zun_adout    )
      zsp2_2 = DOT_PRODUCT( zvn_tlin,    zvn_adout    )
      zsp2_3 = DOT_PRODUCT( zhdivn_tlin, zhdivn_adout )
      zsp2   = zsp2_1 + zsp2_2 + zsp2_3

      cl_name = 'div_cur_adj T1'
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      !=============================================================
      ! 2) dx = ( un_tl, vn_tl, rotn_tl ) and
      !    dy = ( rotb_tl, rotn_tl )
      !=============================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      un_tl   (:,:,:) = 0.0_wp
      vn_tl   (:,:,:) = 0.0_wp
      hdivb_tl(:,:,:) = 0.0_wp
      hdivn_tl(:,:,:) = 0.0_wp
      rotb_tl (:,:,:) = 0.0_wp
      rotn_tl (:,:,:) = 0.0_wp
      hdivb_ad(:,:,:) = 0.0_wp
      hdivn_ad(:,:,:) = 0.0_wp
      rotb_ad (:,:,:) = 0.0_wp
      rotn_ad (:,:,:) = 0.0_wp
      un_ad   (:,:,:) = 0.0_wp
      vn_ad   (:,:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      un_tl  (:,:,:) = zun_tlin  (:,:,:)
      vn_tl  (:,:,:) = zvn_tlin  (:,:,:)
      rotn_tl(:,:,:) = zrotn_tlin(:,:,:)

      CALL div_cur_tan( nit000 )

      zrotb_tlout(:,:,:) = rotb_tl(:,:,:)
      zrotn_tlout(:,:,:) = rotn_tl(:,:,:)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zrotb_adin(ji,jj,jk) = zrotb_tlout(ji,jj,jk) &
                 &               * e1f(ji,jj) * e2f(ji,jj) * e3f(ji,jj,jk)
              zrotn_adin(ji,jj,jk) = zrotn_tlout(ji,jj,jk) &
                 &               * e1f(ji,jj) * e2f(ji,jj) * e3f(ji,jj,jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1_1 = DOT_PRODUCT( zrotb_tlout, zrotb_adin )
      zsp1_2 = DOT_PRODUCT( zrotn_tlout, zrotn_adin )
      zsp1   = zsp1_1 + zsp1_2

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      rotb_ad (:,:,:) = zrotb_adin(:,:,:)
      rotn_ad (:,:,:) = zrotn_adin(:,:,:)
      hdivb_ad(:,:,:) = 0.0_wp
      hdivn_ad(:,:,:) = 0.0_wp

      CALL div_cur_adj( nit000 )

      zun_adout  (:,:,:) = un_ad  (:,:,:)
      zvn_adout  (:,:,:) = vn_ad  (:,:,:)
      zrotn_adout(:,:,:) = rotn_ad(:,:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1 = DOT_PRODUCT( zun_tlin,   zun_adout   )
      zsp2_2 = DOT_PRODUCT( zvn_tlin,   zvn_adout   )
      zsp2_3 = DOT_PRODUCT( zrotn_tlin, zrotn_adout )
      zsp2   = zsp2_1 + zsp2_2 + zsp2_3

      cl_name = 'div_cur_adj T2'
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      !=============================================================
      ! 3) dx = ( un_tl, vn_tl, rotn_tl, hdin_tl ) and
      !    dy = ( rotb_tl, rotn_tl, hdivn_tl, hdivb_tl )
      !=============================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------

      un_tl   (:,:,:) = 0.0_wp
      vn_tl   (:,:,:) = 0.0_wp
      hdivb_tl(:,:,:) = 0.0_wp
      hdivn_tl(:,:,:) = 0.0_wp
      rotb_tl (:,:,:) = 0.0_wp
      rotn_tl (:,:,:) = 0.0_wp
      hdivb_ad(:,:,:) = 0.0_wp
      hdivn_ad(:,:,:) = 0.0_wp
      rotb_ad (:,:,:) = 0.0_wp
      rotn_ad (:,:,:) = 0.0_wp
      un_ad   (:,:,:) = 0.0_wp
      vn_ad   (:,:,:) = 0.0_wp

      !--------------------------------------------------------------------
      ! Call the tangent routine: dy = L dx
      !--------------------------------------------------------------------

      un_tl  (:,:,:) = zun_tlin  (:,:,:)
      vn_tl  (:,:,:) = zvn_tlin  (:,:,:)
      rotn_tl(:,:,:) = zrotn_tlin(:,:,:)
      hdivn_tl(:,:,:) = zhdivn_tlin(:,:,:)

      CALL div_cur_tan( nit000 )

      zhdivb_tlout(:,:,:) = hdivb_tl(:,:,:)
      zhdivn_tlout(:,:,:) = hdivn_tl(:,:,:)
      zrotb_tlout(:,:,:) = rotb_tl(:,:,:)
      zrotn_tlout(:,:,:) = rotn_tl(:,:,:)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zhdivb_adin(ji,jj,jk) = zhdivb_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk)
              zhdivn_adin(ji,jj,jk) = zhdivn_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk)
            END DO
         END DO
      END DO
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zrotb_adin(ji,jj,jk) = zrotb_tlout(ji,jj,jk) &
                 &               * e1f(ji,jj) * e2f(ji,jj) * e3f(ji,jj,jk)
              zrotn_adin(ji,jj,jk) = zrotn_tlout(ji,jj,jk) &
                 &               * e1f(ji,jj) * e2f(ji,jj) * e3f(ji,jj,jk)
            END DO
         END DO
      END DO

      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1_1 = DOT_PRODUCT( zhdivb_tlout, zhdivb_adin )
      zsp1_2 = DOT_PRODUCT( zhdivn_tlout, zhdivn_adin )
      zsp1_3 = DOT_PRODUCT( zrotb_tlout, zrotb_adin )
      zsp1_4 = DOT_PRODUCT( zrotn_tlout, zrotn_adin )
      zsp1   = zsp1_1 + zsp1_2 + zsp1_3 + zsp1_4

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      hdivb_ad(:,:,:) = zhdivb_adin(:,:,:)
      hdivn_ad(:,:,:) = zhdivn_adin(:,:,:)
      rotb_ad (:,:,:) = zrotb_adin(:,:,:)
      rotn_ad (:,:,:) = zrotn_adin(:,:,:)

      CALL div_cur_adj( nit000 )

      zun_adout  (:,:,:) = un_ad  (:,:,:)
      zvn_adout  (:,:,:) = vn_ad  (:,:,:)
      zrotn_adout(:,:,:) = rotn_ad(:,:,:)
      zhdivn_adout(:,:,:) = hdivn_ad(:,:,:)

      !--------------------------------------------------------------------
      ! Compute the scalar product: dx^T L^T W dy
      !--------------------------------------------------------------------

      zsp2_1 = DOT_PRODUCT( zun_tlin,   zun_adout   )
      zsp2_2 = DOT_PRODUCT( zvn_tlin,   zvn_adout   )
      zsp2_3 = DOT_PRODUCT( zrotn_tlin, zrotn_adout )
      zsp2_4 = DOT_PRODUCT( zhdivn_tlin, zhdivn_adout )
      zsp2   = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4

      cl_name = 'div_cur_adj T3'
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )


      DEALLOCATE( &
         & zun_tlin,     &
         & zvn_tlin,     &
         & zhdivn_tlin,  &
         & zrotn_tlin,   &
         & zhdivb_tlout, &
         & zhdivn_tlout, &
         & zrotb_tlout,  &
         & zrotn_tlout,  &
         & zhdivb_adin,  &
         & zhdivn_adin,  &
         & zrotb_adin,   &
         & zrotn_adin,   &
         & zun_adout,    &
         & zvn_adout,    &
         & zhdivn_adout, &
         & zrotn_adout,  &
         & znu,          &
         & znv           &
         & )

   END SUBROUTINE div_cur_adj_tst

   !!======================================================================

END MODULE divcur_tam
