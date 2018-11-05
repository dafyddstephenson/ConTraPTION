MODULE divcur
   !!==============================================================================
   !!                       ***  MODULE  divcur  ***
   !! Ocean diagnostic variable : horizontal divergence and relative vorticity
   !!==============================================================================
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
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   div_cur    : Compute the horizontal divergence and relative
   !!                vorticity fields
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE sbc_oce, ONLY : ln_rnf   ! surface boundary condition: ocean
   USE sbcrnf          ! river runoff 
   USE cla             ! cross land advection             (cla_div routine)
   USE in_out_manager  ! I/O manager
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE lib_mpp         ! MPP library
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   div_cur    ! routine called by step.F90 and istate.F90

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
   !! $Id: divcur.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   !!----------------------------------------------------------------------
   !!   Default option                           2nd order centered schemes
   !!----------------------------------------------------------------------

   SUBROUTINE div_cur( kt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE div_cur  ***
      !!                    
      !! ** Purpose :   compute the horizontal divergence and the relative
      !!      vorticity at before and now time-step
      !!
      !! ** Method  : - Divergence:
      !!      - save the divergence computed at the previous time-step
      !!      (note that the Asselin filter has not been applied on hdivb)
      !!      - compute the now divergence given by :
      !!         hdivn = 1/(e1t*e2t*e3t) ( di[e2u*e3u un] + dj[e1v*e3v vn] )
      !!      correct hdiv with runoff inflow (div_rnf) and cross land flow (div_cla) 
      !!              - Relavtive Vorticity :
      !!      - save the curl computed at the previous time-step (rotb = rotn)
      !!      (note that the Asselin time filter has not been applied to rotb)
      !!      - compute the now curl in tensorial formalism:
      !!            rotn = 1/(e1f*e2f) ( di[e2v vn] - dj[e1u un] )
      !!      Note: Coastal boundary condition: lateral friction set through
      !!      the value of fmask along the coast (see dommsk.F90) and shlat
      !!      (namelist parameter)
      !!
      !! ** Action  : - update hdivb, hdivn, the before & now hor. divergence
      !!              - update rotb , rotn , the before & now rel. vorticity
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk    ! dummy loop indices
      REAL(wp) ::   zraur, zdep   ! local scalars
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('div_cur')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cur : horizontal velocity divergence and'
         IF(lwp) WRITE(numout,*) '~~~~~~~   relative vorticity'
      ENDIF

      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         !
         hdivb(:,:,jk) = hdivn(:,:,jk)    ! time swap of div arrays
         rotb (:,:,jk) = rotn (:,:,jk)    ! time swap of rot arrays
         !
         !                                             ! --------
         ! Horizontal divergence                       !   div 
         !                                             ! --------
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               hdivn(ji,jj,jk) =   &
                  (  e2u(ji,jj)*e3u(ji,jj,jk) * un(ji,jj,jk) - e2u(ji-1,jj)*e3u(ji-1,jj,jk) * un(ji-1,jj,jk)       &
                   + e1v(ji,jj)*e3v(ji,jj,jk) * vn(ji,jj,jk) - e1v(ji,jj-1)*e3v(ji,jj-1,jk) * vn(ji,jj-1,jk)  )    &
                  / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO  
         END DO  

         IF( .NOT. AGRIF_Root() ) THEN
            IF ((nbondi ==  1).OR.(nbondi == 2)) hdivn(nlci-1 , :     ,jk) = 0.e0      ! east
            IF ((nbondi == -1).OR.(nbondi == 2)) hdivn(2      , :     ,jk) = 0.e0      ! west
            IF ((nbondj ==  1).OR.(nbondj == 2)) hdivn(:      ,nlcj-1 ,jk) = 0.e0      ! north
            IF ((nbondj == -1).OR.(nbondj == 2)) hdivn(:      ,2      ,jk) = 0.e0      ! south
         ENDIF

         !                                             ! --------
         ! relative vorticity                          !   rot 
         !                                             ! --------
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               rotn(ji,jj,jk) = (  e2v(ji+1,jj  ) * vn(ji+1,jj  ,jk) - e2v(ji,jj) * vn(ji,jj,jk)    &
                  &              - e1u(ji  ,jj+1) * un(ji  ,jj+1,jk) + e1u(ji,jj) * un(ji,jj,jk)  ) &
                  &           * fmask(ji,jj,jk) / ( e1f(ji,jj) * e2f(ji,jj) )
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============

      IF( ln_rnf      )   CALL sbc_rnf_div( hdivn )          ! runoffs (update hdivn field)
      IF( nn_cla == 1 )   CALL cla_div    ( kt )             ! Cross Land Advection (update hdivn field)
      !
      CALL lbc_lnk( hdivn, 'T', 1. )   ;   CALL lbc_lnk( rotn , 'F', 1. )     ! lateral boundary cond. (no sign change)
      !
      IF( nn_timing == 1 )  CALL timing_stop('div_cur')
      !
   END SUBROUTINE div_cur
   
   !!======================================================================
END MODULE divcur
