MODULE traadv_eiv
   !!======================================================================
   !!                    ***  MODULE  traadv_eiv  ***
   !! Ocean tracers:  advection trend - eddy induced velocity
   !!======================================================================
   !! History :  1.0  !  2005-11 (G. Madec)  Original code, from traldf and zdf _iso
   !!            3.3  !  2010-05 (C. Ethe, G. Madec)  merge TRC-TRA 
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_traldf_eiv'                  rotation of the lateral mixing tensor
   !!----------------------------------------------------------------------
   !!   tra_ldf_iso : update the tracer trend with the horizontal component
   !!                 of iso neutral laplacian operator or horizontal 
   !!                 laplacian operator in s-coordinate
   !!----------------------------------------------------------------------
   USE oce             ! ocean dynamics and tracers variables
   USE dom_oce         ! ocean space and time domain variables
   USE ldftra_oce      ! ocean active tracers: lateral physics
   USE ldfslp          ! iso-neutral slopes
   USE in_out_manager  ! I/O manager
   USE iom
   USE trc_oce         ! share passive tracers/Ocean variables
   USE phycst          ! physical constants
   USE lbclnk          ! ocean lateral boundary conditions (or mpp link)
   USE diaar5, ONLY:   lk_diaar5
   USE wrk_nemo        ! Memory Allocation
   USE timing          ! Timing

   IMPLICIT NONE
   PRIVATE

   PUBLIC   tra_adv_eiv   ! routine called by step.F90

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
   !!                    *** ldftra_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldftra_substitute.h90 3294 2012-01-28 16:44:18Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'key_traldf_c2d' :                 aht: 2D coefficient
   !!----------------------------------------------------------------------
   !!                   ***  ldfeiv_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaei. the eddy induced velocity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: ldfeiv_substitute.h90 2528 2010-12-27 17:33:53Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
!   'traldf_c2d' :                           eiv: 2D coefficient
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
   !! $Id: traadv_eiv.F90 3788 2013-02-10 12:14:59Z gm $
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_adv_eiv( kt, kit000, pun, pvn, pwn, cdtype )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_adv_eiv  ***
      !! 
      !! ** Purpose :   Compute the before horizontal tracer (t & s) diffusive 
      !!      trend and add it to the general trend of tracer equation.
      !!
      !! ** Method  :   The eddy induced advection is computed from the slope
      !!      of iso-neutral surfaces computed in routine ldf_slp as follows:
      !!         zu_eiv =  1/(e2u e3u)   dk[ aeiu e2u mi(wslpi) ]
      !!         zv_eiv =  1/(e1v e3v)   dk[ aeiv e1v mj(wslpj)
      !!         zw_eiv = -1/(e1t e2t) { di[ aeiu e2u mi(wslpi) ]
      !!                               + dj[ aeiv e1v mj(wslpj) ] }
      !!      add the eiv component to the model velocity:
      !!         p.n = p.n + z._eiv
      !!
      !! ** Action  : - add to p.n the eiv component
      !!----------------------------------------------------------------------
      INTEGER                         , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                         , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3)                , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pun      ! in : 3 ocean velocity components 
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pvn      ! out: 3 ocean velocity components
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pwn      ! increased by the eiv
      !!
      INTEGER  ::   ji, jj, jk                 ! dummy loop indices
      REAL(wp) ::   zuwk, zuwk1, zuwi, zuwi1   ! local scalars
      REAL(wp) ::   zvwk, zvwk1, zvwj, zvwj1   !   -      -
      REAL(wp) ::   zztmp                      ! local scalar
      REAL(wp), POINTER, DIMENSION(:,:) :: zu_eiv, zv_eiv, zw_eiv, z2d
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_adv_eiv')
      !
      CALL wrk_alloc( jpi, jpj, zu_eiv, zv_eiv, zw_eiv, z2d )

      IF( kt == kit000 )  THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_adv_eiv : eddy induced advection on ', cdtype,' :'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~   add to velocity fields the eiv component'
         IF( cdtype == 'TRA') THEN
            u_eiv(:,:,:) = 0.e0
            v_eiv(:,:,:) = 0.e0
            w_eiv(:,:,:) = 0.e0
         END IF
      ENDIF

      zu_eiv(:,:) = 0.e0   ;   zv_eiv(:,:) = 0.e0   ;    zw_eiv(:,:) = 0.e0  
      
                                                    ! =================
      DO jk = 1, jpkm1                              !  Horizontal slab
         !                                          ! =================
         DO jj = 1, jpjm1
            DO ji = 1, jpim1   ! vector opt.
               zuwk = ( wslpi(ji,jj,jk  ) + wslpi(ji+1,jj,jk  ) ) * aeiu(ji,jj) * umask(ji,jj,jk  )
               zuwk1= ( wslpi(ji,jj,jk+1) + wslpi(ji+1,jj,jk+1) ) * aeiu(ji,jj) * umask(ji,jj,jk+1)
               zvwk = ( wslpj(ji,jj,jk  ) + wslpj(ji,jj+1,jk  ) ) * aeiv(ji,jj) * vmask(ji,jj,jk  )
               zvwk1= ( wslpj(ji,jj,jk+1) + wslpj(ji,jj+1,jk+1) ) * aeiv(ji,jj) * vmask(ji,jj,jk+1)

               zu_eiv(ji,jj) = 0.5 * umask(ji,jj,jk) * ( zuwk - zuwk1 ) 
               zv_eiv(ji,jj) = 0.5 * vmask(ji,jj,jk) * ( zvwk - zvwk1 ) 
   
               pun(ji,jj,jk) = pun(ji,jj,jk) + e2u(ji,jj) * zu_eiv(ji,jj)
               pvn(ji,jj,jk) = pvn(ji,jj,jk) + e1v(ji,jj) * zv_eiv(ji,jj)
            END DO
         END DO
         IF( cdtype == 'TRA') THEN
            u_eiv(:,:,jk) = zu_eiv(:,:) / e3u(:,:,jk)
            v_eiv(:,:,jk) = zv_eiv(:,:) / e3v(:,:,jk)
         END IF
         IF( jk >=2 ) THEN                             ! jk=1 zw_eiv=0, not computed
            DO jj = 2, jpjm1
               DO ji = 2, jpim1   ! vector opt.
                  zuwi  = ( wslpi(ji,jj,jk)+wslpi(ji-1,jj,jk) ) * aeiu(ji-1,jj) * e2u(ji-1,jj) * umask(ji-1,jj,jk)
                  zuwi1 = ( wslpi(ji,jj,jk)+wslpi(ji+1,jj,jk) ) * aeiu(ji  ,jj) * e2u(ji  ,jj) * umask(ji  ,jj,jk)
                  zvwj  = ( wslpj(ji,jj,jk)+wslpj(ji,jj-1,jk) ) * aeiv(ji,jj-1) * e1v(ji,jj-1) * vmask(ji,jj-1,jk)
                  zvwj1 = ( wslpj(ji,jj,jk)+wslpj(ji,jj+1,jk) ) * aeiv(ji,jj  ) * e1v(ji  ,jj) * vmask(ji  ,jj,jk)
  
                  zw_eiv(ji,jj) = - 0.5 * tmask(ji,jj,jk) * ( zuwi1 - zuwi + zvwj1 - zvwj ) 
                  pwn(ji,jj,jk) = pwn(ji,jj,jk) + zw_eiv(ji,jj)
               END DO
            END DO
            IF( cdtype == 'TRA')  w_eiv(:,:,jk) = zw_eiv(:,:) / ( e1t(:,:) * e2t(:,:) )
         ENDIF
         !                                          ! =================
      END DO                                        !    End of slab  
      !                                             ! =================

      IF( cdtype == 'TRA') THEN
         CALL iom_put( "uoce_eiv", u_eiv )    ! i-eiv current
         CALL iom_put( "voce_eiv", v_eiv )    ! j-eiv current
         CALL iom_put( "woce_eiv", w_eiv )    ! vert. eiv current
         IF( lk_diaar5 ) THEN
            zztmp = 0.5 * rau0 * rcp 
            z2d(:,:) = 0.e0 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z2d(ji,jj) = z2d(ji,jj) + zztmp * u_eiv(ji,jj,jk) &
                       &         * (tsn(ji,jj,jk,jp_tem)+tsn(ji+1,jj,jk,jp_tem)) * e2u(ji,jj) * e3u(ji,jj,jk) 
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z2d, 'U', -1. )
            CALL iom_put( "ueiv_heattr", z2d )                  ! heat transport in i-direction
            z2d(:,:) = 0.e0 
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     z2d(ji,jj) = z2d(ji,jj) + zztmp * v_eiv(ji,jj,jk) &
                     &           * (tsn(ji,jj,jk,jp_tem)+tsn(ji,jj+1,jk,jp_tem)) * e1v(ji,jj) * e3v(ji,jj,jk) 
                  END DO
               END DO
            END DO
            CALL lbc_lnk( z2d, 'V', -1. )
            CALL iom_put( "veiv_heattr", z2d )                  !  heat transport in i-direction
         ENDIF
    END IF
      ! 
      CALL wrk_dealloc( jpi, jpj, zu_eiv, zv_eiv, zw_eiv, z2d )
      !
      IF( nn_timing == 1 )  CALL timing_stop( 'tra_adv_eiv')
      !
    END SUBROUTINE tra_adv_eiv


   !!==============================================================================
END MODULE traadv_eiv
