MODULE dynldf_lap
   !!======================================================================
   !!                       ***  MODULE  dynldf_lap  ***
   !! Ocean dynamics:  lateral viscosity trend
   !!======================================================================
   !! History :  OPA  ! 1990-09 (G. Madec) Original code
   !!            4.0  ! 1991-11 (G. Madec)
   !!            6.0  ! 1996-01 (G. Madec) statement function for e3 and ahm
   !!   NEMO     1.0  ! 2002-06 (G. Madec)  F90: Free form and module
   !!             -   ! 2004-08 (C. Talandier) New trends organization
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   dyn_ldf_lap  : update the momentum trend with the lateral diffusion
   !!                  using an iso-level harmonic operator
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers
   USE dom_oce         ! ocean space and time domain
   USE ldfdyn_oce      ! ocean dynamics: lateral physics
   USE zdf_oce         ! ocean vertical physics
   USE in_out_manager  ! I/O manager
   USE trdmod          ! ocean dynamics trends 
   USE trdmod_oce      ! ocean variables trends
   USE ldfslp          ! iso-neutral slopes 
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC dyn_ldf_lap  ! called by step.F90

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
   !!                    ***  ldfdyn_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsahm., the lateral eddy viscosity coeff. 
   !!      with a constant, or 1D, or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfdyn_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   !!
   !! fsahmt, fsahmf - used for laplaian operator only
   !! fsahmu, fsahmv - used for bilaplacian operator only
   !!
!   ' key_dynldf_c3d' :                  3D coefficient
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
   !! $Id: dynldf_lap.F90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE dyn_ldf_lap( kt )
      !!----------------------------------------------------------------------
      !!                     ***  ROUTINE dyn_ldf_lap  ***
      !!                       
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The before horizontal momentum diffusion trend is an
      !!      harmonic operator (laplacian type) which separates the divergent
      !!      and rotational parts of the flow.
      !!      Its horizontal components are computed as follow:
      !!         difu = 1/e1u di[ahmt hdivb] - 1/(e2u*e3u) dj-1[e3f ahmf rotb]
      !!         difv = 1/e2v dj[ahmt hdivb] + 1/(e1v*e3v) di-1[e3f ahmf rotb]
      !!      in the rotational part of the diffusion.
      !!      Add this before trend to the general trend (ua,va):
      !!            (ua,va) = (ua,va) + (diffu,diffv)
      !!      'key_trddyn' activated: the two components of the horizontal
      !!                                 diffusion trend are saved.
      !!
      !! ** Action : - Update (ua,va) with the before iso-level harmonic 
      !!               mixing trend.
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt   ! ocean time-step index
      !
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      REAL(wp) ::   zua, zva, ze2u, ze1v   ! local scalars
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('dyn_ldf_lap')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'dyn_ldf : iso-level harmonic (laplacian) operator'
         IF(lwp) WRITE(numout,*) '~~~~~~~ '
      ENDIF
      !                                                ! ===============
      DO jk = 1, jpkm1                                 ! Horizontal slab
         !                                             ! ===============
         DO jj = 2, jpjm1
            DO ji = 2, jpim1   ! vector opt.
               ze2u = rotb (ji,jj,jk) * ahm2(ji,jj,jk) * e3f(ji,jj,jk)
               ze1v = hdivb(ji,jj,jk) * ahm1(ji,jj,jk)
               ! horizontal diffusive trends
               zua = - ( ze2u - rotb (ji,jj-1,jk)*ahm2(ji,jj-1,jk)*e3f(ji,jj-1,jk) ) / ( e2u(ji,jj) * e3u(ji,jj,jk) )   &
                     + ( hdivb(ji+1,jj,jk)*ahm1(ji+1,jj,jk) - ze1v                   ) / e1u(ji,jj)

               zva = + ( ze2u - rotb (ji-1,jj,jk)*ahm2(ji-1,jj,jk)*e3f(ji-1,jj,jk) ) / ( e1v(ji,jj) * e3v(ji,jj,jk) )   &
                     + ( hdivb(ji,jj+1,jk)*ahm1(ji,jj+1,jk) - ze1v                   ) / e2v(ji,jj)

               ! add it to the general momentum trends
               ua(ji,jj,jk) = ua(ji,jj,jk) + zua
               va(ji,jj,jk) = va(ji,jj,jk) + zva
            END DO
         END DO
         !                                             ! ===============
      END DO                                           !   End of slab
      !                                                ! ===============
      IF( nn_timing == 1 )  CALL timing_stop('dyn_ldf_lap')
      !
   END SUBROUTINE dyn_ldf_lap

   !!======================================================================
END MODULE dynldf_lap
