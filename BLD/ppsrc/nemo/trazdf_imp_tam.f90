MODULE trazdf_imp_tam
   !!==============================================================================
   !!                 ***  MODULE  trazdf_imp_tam  ***
   !! Ocean active tracers:  vertical component of the tracer mixing trend
   !!                        Tangent and Adjoint Module
   !!==============================================================================
   !! History of the direct module:
   !!            OPA  !  1990-10  (B. Blanke)  Original code
   !!            7.0  !  1991-11  (G. Madec)
   !!                 !  1992-06  (M. Imbard) correction on tracer trend loops
   !!                 !  1996-01  (G. Madec) statement function for e3
   !!                 !  1997-05  (G. Madec) vertical component of isopycnal
   !!                 !  1997-07  (G. Madec) geopotential diffusion in s-coord
   !!                 !  2000-08  (G. Madec) double diffusive mixing
   !!   NEMO     1.0  !  2002-08  (G. Madec) F90: Free form and module
   !!            2.0  !  2006-11  (G. Madec) New step reorganisation
   !!            3.2  !  2009-03  (G. Madec)  heat and salt content trends

   !! History of the T&A module:
   !!        !  09-01 (A. Vidard) tam of the 06-11 version
   !!----------------------------------------------------------------------
   !!   tra_zdf_imp_tan : Update the tracer trend with the diagonal vertical
   !!                 part of the mixing tensor (tangent).
   !!   tra_zdf_imp_adj : Update the tracer trend with the diagonal vertical
   !!                 part of the mixing tensor (adjoint).
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind
   USE par_oce
   USE oce_tam
   USE dom_oce
   USE oce
   USE zdf_oce
   USE ldftra_oce
   USE zdfddm
   USE traldf_tam
   USE in_out_manager
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE trc_oce
   USE trc_oce_tam
   USE ldftra
   USE lib_mpp
   USE wrk_nemo
   USE timing
   USE ldfslp
   USE paresp

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC tra_zdf_imp_tan       !  routine called by tra_zdf_tan.F90
   PUBLIC tra_zdf_imp_adj       !  routine called by tra_zdf_adj.F90
   PUBLIC tra_zdf_imp_adj_tst   !  routine called by tst.F90

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
   !!                    *** zdfddm_substitute.h90  ***
   !!----------------------------------------------------------------------
   !! ** purpose :   substitute fsaht. the eddy diffusivity coeff.
   !!      with a constant or 1D or 2D or 3D array, using CPP macro.
   !!----------------------------------------------------------------------
!   Defautl option :                     avs = avt
   !!----------------------------------------------------------------------
   !! NEMO/OPA 4.0 , NEMO Consortium (2011)
   !! $Id: zdfddm_substitute.h90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
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
   !!----------------------------------------------------------------------
   !!  OPA 9.0 , LOCEAN-IPSL (2005)
   !! $Id: trazdf_imp.F90 1156 2008-06-26 16:06:45Z rblod $
   !! Software governed by the CeCILL licence (modipsl/doc/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE tra_zdf_imp_tan( kt, kit000, cdtype, p2dt, ptb_tl, pta_tl, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_imp_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!     Compute the trend due to the vertical tracer diffusion
      !!     including the vertical component of lateral mixing (only for 2nd
      !!     order operator, for fourth order it is already computed and add
      !!     to the general trend in traldf.F) and add it to the general trend
      !!     of the tracer equations.
      !!
      !! ** Method of the direct routine :
      !!      The vertical component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neutral or geo-
      !!      potential surfaces to which an eddy induced advection can be
      !!      added. It is computed using before fields (forward in time) and
      !!      isopycnal or geopotential slopes computed in routine ldfslp.
      !!
      !!      Second part: vertical trend associated with the vertical physics
      !!      ===========  (including the vertical flux proportional to dk[t]
      !!                  associated with the lateral mixing, through the
      !!                  update of avt)
      !!      The vertical diffusion of tracers (t & s) is given by:
      !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme (t=ta).
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend ta,sa :
      !!         ta = ta + dz( avt dz(t) )
      !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
      !!
      !!      Third part: recover avt resulting from the vertical physics
      !!      ==========  alone, for further diagnostics (for example to
      !!                  compute the turbocline depth in zdfmxl.F90).
      !!         avt = zavt
      !!         (avs = zavs if lk_zdfddm=T )
      !!
      !! ** Remarks on the tangent routine : - key_vvl is not available in tangent yet.
      !!    Once it will be this routine wil need to be rewritten
      !!                                     - simplified version, slopes (wslp[ij])
      !!     assumed to be constant (read from the trajectory). same for av[ts]
      !!
      !!---------------------------------------------------------------------
      !!
      INTEGER                              , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000   ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt     ! number of tracers
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt     ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(in   ) ::   ptb_tl   ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_tl   ! tracer trend
      !! * Local declarations
      INTEGER  ::   ji, jj, jk, jn           ! dummy loop indices
      REAL(wp) ::   zavi, zrhstl, znvvl,   & ! temporary scalars
         ze3tb, ze3tn, ze3ta, zvsfvvl        ! variable vertical scale factors
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   &
         zwi, zwt, zwd, zws                     ! workspace arrays
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_imp_tan')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwi, zwt, zwd, zws )
      !
      IF( kt == kit000 ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'tra_zdf_imp_tan : implicit vertical mixing on ', cdtype
         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~ '
      ENDIF

      ! I.1 Variable volume : to take into account vertical variable vertical scale factors
      ! -------------------
      ! ... not available  in tangent yet
      ! II. Vertical trend associated with the vertical physics
      ! =======================================================
      !     (including the vertical flux proportional to dk[t] associated
      !      with the lateral mixing, through the avt update)
      !     dk[ avt dk[ (t,s) ] ] diffusive trends
      DO jn = 1, kjpt                                 !  tracer loop  !
         !                                            ! ============= !
         !
         !  Matrix construction
         ! --------------------
         ! Build matrix if temperature or salinity (only in double diffusion case) or first passive tracer
         !
         IF(  ( cdtype == 'TRA' .AND. ( jn == jp_tem .OR. ( jn == jp_sal .AND. lk_zdfddm ) ) ) .OR.   &
            & ( cdtype == 'TRC' .AND. jn == 1 )  )  THEN
            !
            ! vertical mixing coef.: avt for temperature, avs for salinity and passive tracers
            IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN   ;   zwt(:,:,2:jpk) = avt  (:,:,2:jpk)
            ELSE                                            ;   zwt(:,:,2:jpk) = avt(:,:,2:jpk)
            ENDIF
            zwt(:,:,1) = 0._wp
            !
            ! II.0 Matrix construction
            ! ------------------------
            ! update and save of avt (and avs if double diffusive mixing)
            IF ( ln_traldf_grif ) THEN
               IF ( lwp ) WRITE(numout, *) 'Griffies operator for lateral tracer diffusion not avaible in TAM yet'
               CALL abort
            ELSE IF( l_traldf_rot ) THEN
               DO jk = 2, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        zwt(ji,jj,jk) = zwt(ji,jj,jk) + rldf * ahtw(ji,jj)                       &   ! vertical mixing coef. due to lateral mixing
                           & * (  wslpi(ji,jj,jk) * wslpi(ji,jj,jk)   &
                           &    + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
                     END DO
                  END DO
               END DO
            ENDIF
            ! Diagonal, inferior, superior  (including the bottom boundary condition via avt masked)
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ze3ta = 1._wp                                ! after scale factor at T-point
                     ze3tn = e3t(ji,jj,jk)                      ! now   scale factor at T-point
                     zwi(ji,jj,jk) = - p2dt(jk) * zwt(ji,jj,jk  ) / ( ze3tn * e3w(ji,jj,jk  ) )
                     zws(ji,jj,jk) = - p2dt(jk) * zwt(ji,jj,jk+1) / ( ze3tn * e3w(ji,jj,jk+1) )
                     zwd(ji,jj,jk) = ze3ta - zwi(ji,jj,jk) - zws(ji,jj,jk)
                  END DO
               END DO
            END DO
            !
            ! II.1. Vertical diffusion on t
            ! ---------------------------
            !
            !! Matrix inversion from the first level
            !!----------------------------------------------------------------------
            !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
            !
            !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
            !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
            !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
            !        (        ...               )( ...  ) ( ...  )
            !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
            !
            !   m is decomposed in the product of an upper and lower triangular matrix
            !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
            !   The second member is in 2d array zwy
            !   The solution is in 2d array zwx
            !   The 3d arry zwt is a work space array
            !   zwy is used and then used as a work space array : its value is modified!

            ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zwt(ji,jj,1) = zwd(ji,jj,1)
               END DO
            END DO
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
                  END DO
               END DO
            END DO
         END IF
         ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ze3tb = 1._wp
               ze3tn = 1._wp
               pta_tl(ji,jj,1,jn) = ze3tb * ptb_tl(ji,jj,1,jn) + p2dt(1) * ze3tn * pta_tl(ji,jj,1,jn)
            END DO
         END DO

         DO jk = 2, jpkm1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ze3tb = 1._wp
                  ze3tn = 1._wp
                  zrhstl = ze3tb * ptb_tl(ji,jj,jk,jn) + p2dt(jk) * ze3tn * pta_tl(ji,jj,jk,jn)   ! zrhs=right hand side
                  pta_tl(ji,jj,jk,jn) = zrhstl - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) * pta_tl(ji,jj,jk-1,jn)
               END DO
            END DO
         END DO

         ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
         ! Save the masked temperature after in ta
         ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               pta_tl(ji,jj,jpkm1,jn) = pta_tl(ji,jj,jpkm1,jn) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
            END DO
         END DO
         DO jk = jpk-2, 1, -1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  pta_tl(ji,jj,jk,jn) = ( pta_tl(ji,jj,jk,jn) - zws(ji,jj,jk) * pta_tl(ji,jj,jk+1,jn) )   &
                                  &   / zwt(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
      !                                            ! ================= !
      END DO                                          !  end tracer loop  !
      !                                               ! ================= !
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwi, zwt, zwd, zws )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_imp_tan')
      !
   END SUBROUTINE tra_zdf_imp_tan
   SUBROUTINE tra_zdf_imp_adj( kt, kit000, cdtype, p2dt, ptb_ad, pta_ad, kjpt )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE tra_zdf_imp_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!     Compute the trend due to the vertical tracer diffusion
      !!     including the vertical component of lateral mixing (only for 2nd
      !!     order operator, for fourth order it is already computed and add
      !!     to the general trend in traldf.F) and add it to the general trend
      !!     of the tracer equations.
      !!
      !! ** Method of the direct routine :
      !!      The vertical component of the lateral diffusive trends
      !!      is provided by a 2nd order operator rotated along neutral or geo-
      !!      potential surfaces to which an eddy induced advection can be
      !!      added. It is computed using before fields (forward in time) and
      !!      isopycnal or geopotential slopes computed in routine ldfslp.
      !!
      !!      Second part: vertical trend associated with the vertical physics
      !!      ===========  (including the vertical flux proportional to dk[t]
      !!                  associated with the lateral mixing, through the
      !!                  update of avt)
      !!      The vertical diffusion of tracers (t & s) is given by:
      !!             difft = dz( avt dz(t) ) = 1/e3t dk+1( avt/e3w dk(t) )
      !!      It is computed using a backward time scheme (t=ta).
      !!      Surface and bottom boundary conditions: no diffusive flux on
      !!      both tracers (bottom, applied through the masked field avt).
      !!      Add this trend to the general trend ta,sa :
      !!         ta = ta + dz( avt dz(t) )
      !!         (sa = sa + dz( avs dz(t) ) if lk_zdfddm=T )
      !!
      !!      Third part: recover avt resulting from the vertical physics
      !!      ==========  alone, for further diagnostics (for example to
      !!                  compute the turbocline depth in zdfmxl.F90).
      !!         avt = zavt
      !!         (avs = zavs if lk_zdfddm=T )
      !!
      !! ** Remarks on the adjoint routine : - key_vvl is not available in adjoint yet.
      !!    Once it will be this routine wil need to be rewritten
      !!                                     - simplified version, slopes (wslp[ij])
      !!     assumed to be constant (read from the trajectory). same for av[ts]
      !!
      !!---------------------------------------------------------------------
      INTEGER                              , INTENT(in   ) ::   kt       ! ocean time-step index
      INTEGER                              , INTENT(in   ) ::   kit000          ! first time step index
      CHARACTER(len=3)                     , INTENT(in   ) ::   cdtype   ! =TRA or TRC (tracer indicator)
      INTEGER                              , INTENT(in   ) ::   kjpt     ! number of tracers
      REAL(wp), DIMENSION(        jpk     ), INTENT(in   ) ::   p2dt     ! vertical profile of tracer time-step
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   ptb_ad      ! before and now tracer fields
      REAL(wp), DIMENSION(jpi,jpj,jpk,kjpt), INTENT(inout) ::   pta_ad      ! tracer trend
      !! * Local declarations
      INTEGER  ::   ji, jj, jk, jn               ! dummy loop indices
      REAL(wp) ::   zavi, zrhsad, znvvl,   & ! temporary scalars
         ze3tb, ze3tn, ze3ta, zvsfvvl        ! variable vertical scale factors
      REAL(wp), POINTER, DIMENSION(:,:,:) ::   &
         zwi, zwt, zws, zwd                     ! workspace arrays
      !!---------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_zdf_imp_adj')
      !
      CALL wrk_alloc( jpi, jpj, jpk, zwi, zwt, zws, zwd )
      !
      IF( kt == nitend ) THEN
         IF(lwp)WRITE(numout,*)
         IF(lwp)WRITE(numout,*) 'tra_zdf_imp_adj : implicit vertical mixing on', cdtype
         IF(lwp)WRITE(numout,*) '~~~~~~~~~~~~~~~ '
      ENDIF

      ! I.1 Variable volume : to take into account vertical variable vertical scale factors
      ! -------------------
      ! ... not available  in tangent yet
      ! II. Vertical trend associated with the vertical physics
      ! =======================================================
      !     (including the vertical flux proportional to dk[t] associated
      !      with the lateral mixing, through the avt update)
      !     dk[ avt dk[ (t,s) ] ] diffusive trends
      DO jn = 1, kjpt                                 !  tracer loop  !
         !                                            ! ============= !
         !
         !  Matrix construction
         ! --------------------
         ! Build matrix if temperature or salinity (only in double diffusion case) or first passive tracer
         !
         IF(  ( cdtype == 'TRA' .AND. ( jn == jp_tem .OR. ( jn == jp_sal .AND. lk_zdfddm ) ) ) .OR.   &
            & ( cdtype == 'TRC' .AND. jn == 1 )  )  THEN
            !
            ! vertical mixing coef.: avt for temperature, avs for salinity and passive tracers
            IF( cdtype == 'TRA' .AND. jn == jp_tem ) THEN   ;   zwt(:,:,2:jpk) = avt  (:,:,2:jpk)
            ELSE                                            ;   zwt(:,:,2:jpk) = avt(:,:,2:jpk)
            ENDIF
            zwt(:,:,1) = 0._wp

      ! update and save of avt (and avs if double diffusive mixing)
            IF ( ln_traldf_grif ) THEN
               IF ( lwp ) WRITE(numout, *) 'Griffies operator for lateral tracer diffusion not avaible in TAM yet'
               CALL abort
            ELSE IF( l_traldf_rot ) THEN
               DO jk = 2, jpkm1
                  DO jj = 2, jpjm1
                     DO ji = 2, jpim1   ! vector opt.
                        zwt(ji,jj,jk) = zwt(ji,jj,jk) + rldf * ahtw(ji,jj)                       &   ! vertical mixing coef. due to lateral mixing
                           & * (  wslpi(ji,jj,jk) * wslpi(ji,jj,jk)   &
                           &    + wslpj(ji,jj,jk) * wslpj(ji,jj,jk)  )
                     END DO
                  END DO
               END DO
            ENDIF
            ! Diagonal, inferior, superior  (including the bottom boundary condition via avt masked)
            DO jk = 1, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1   ! vector opt.
                     ze3ta = 1._wp                                ! after scale factor at T-point
                     ze3tn = e3t(ji,jj,jk)                      ! now   scale factor at T-point
                     zwi(ji,jj,jk) = - p2dt(jk) * zwt(ji,jj,jk  ) / ( ze3tn * e3w(ji,jj,jk  ) )
                     zws(ji,jj,jk) = - p2dt(jk) * zwt(ji,jj,jk+1) / ( ze3tn * e3w(ji,jj,jk+1) )
                     zwd(ji,jj,jk) = ze3ta - zwi(ji,jj,jk) - zws(ji,jj,jk)
                  END DO
               END DO
            END DO

            !! Matrix inversion from the first level
            !!----------------------------------------------------------------------
            !   solve m.x = y  where m is a tri diagonal matrix ( jpk*jpk )
            !
            !        ( zwd1 zws1   0    0    0  )( zwx1 ) ( zwy1 )
            !        ( zwi2 zwd2 zws2   0    0  )( zwx2 ) ( zwy2 )
            !        (  0   zwi3 zwd3 zws3   0  )( zwx3 )=( zwy3 )
            !        (        ...               )( ...  ) ( ...  )
            !        (  0    0    0   zwik zwdk )( zwxk ) ( zwyk )
            !
            !   m is decomposed in the product of an upper and lower triangular matrix
            !   The 3 diagonal terms are in 2d arrays: zwd, zws, zwi
            !   The second member is in 2d array zwy
            !   The solution is in 2d array zwx
            !   The 3d arry zwt is a work space array
            !   zwy is used and then used as a work space array : its value is modified!
            ! first recurrence:   Tk = Dk - Ik Sk-1 / Tk-1   (increasing k)
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  zwt(ji,jj,1) = zwd(ji,jj,1)
               END DO
            END DO
            DO jk = 2, jpkm1
               DO jj = 2, jpjm1
                  DO ji = 2, jpim1
                     zwt(ji,jj,jk) = zwd(ji,jj,jk) - zwi(ji,jj,jk) * zws(ji,jj,jk-1)  /zwt(ji,jj,jk-1)
                  END DO
               END DO
            END DO
         END IF
         ! third recurrence: Xk = (Zk - Sk Xk+1 ) / Tk
         ! Save the masked temperature after in ta
         ! (c a u t i o n:  temperature not its trend, Leap-frog scheme done it will not be done in tranxt)
         DO jk = 1, jpk-2
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  pta_ad(ji,jj,jk+1,jn) = pta_ad(ji,jj,jk+1,jn) - zws(ji,jj,jk) * pta_ad(ji,jj,jk,jn)   &
                                    &   / zwt(ji,jj,jk) * tmask(ji,jj,jk)
                  pta_ad(ji,jj,jk,jn)   = pta_ad(ji,jj,jk,jn) / zwt(ji,jj,jk) * tmask(ji,jj,jk)
               END DO
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               pta_ad(ji,jj,jpkm1,jn) = pta_ad(ji,jj,jpkm1,jn) / zwt(ji,jj,jpkm1) * tmask(ji,jj,jpkm1)
            END DO
         END DO
         ! second recurrence:    Zk = Yk - Ik / Tk-1  Zk-1
         DO jk = jpkm1, 2, -1
            DO jj = 2, jpjm1
               DO ji = 2, jpim1
                  ze3tb = 1._wp
                  ze3tn = 1._wp
                  zrhsad = pta_ad(ji,jj,jk,jn)
                  pta_ad(ji,jj,jk-1,jn) = pta_ad(ji,jj,jk-1,jn) - zwi(ji,jj,jk) / zwt(ji,jj,jk-1) * pta_ad(ji,jj,jk,jn)
                  pta_ad(ji,jj,jk,jn) = 0.0_wp
                  ptb_ad(ji,jj,jk,jn) = ptb_ad(ji,jj,jk,jn) + ze3tb * zrhsad
                  pta_ad(ji,jj,jk,jn) = pta_ad(ji,jj,jk,jn) + p2dt(jk) * ze3tn * zrhsad
               END DO
            END DO
         END DO
         DO jj = 2, jpjm1
            DO ji = 2, jpim1
               ze3tb = 1._wp
               ze3tn = 1._wp
               ptb_ad(ji,jj,1,jn) = ptb_ad(ji,jj,1,jn) + ze3tb * pta_ad(ji,jj,1,jn)
               pta_ad(ji,jj,1,jn) = pta_ad(ji,jj,1,jn) * p2dt(1) * ze3tn
            END DO
         END DO
      END DO
      !
      CALL wrk_dealloc( jpi, jpj, jpk, zwi, zwt, zws, zwd )
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_zdf_imp_adj')
      !
   END SUBROUTINE tra_zdf_imp_adj
   SUBROUTINE tra_zdf_imp_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE tra_zdf_imp_adj_tst ***
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
         & jk
      REAL(KIND=wp) ::   &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zta_tlin ,     & ! Tangent input
         & ztb_tlin ,     & ! Tangent input
         & zsa_tlin ,     & ! Tangent input
         & zsb_tlin ,     & ! Tangent input
         & zta_tlout,     & ! Tangent output
         & zsa_tlout,     & ! Tangent output
         & zta_adin ,     & ! Adjoint input
         & zsa_adin ,     & ! Adjoint input
         & zta_adout,     & ! Adjoint output
         & ztb_adout,     & ! Adjoint output
         & zsa_adout,     & ! Adjoint output
         & zsb_adout,     & ! Adjoint output
         & zr             ! 3D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory

      ALLOCATE( &
         & zta_tlin( jpi,jpj,jpk),     &
         & zsa_tlin( jpi,jpj,jpk),     &
         & ztb_tlin( jpi,jpj,jpk),     &
         & zsb_tlin( jpi,jpj,jpk),     &
         & zta_tlout(jpi,jpj,jpk),     &
         & zsa_tlout(jpi,jpj,jpk),     &
         & zta_adin( jpi,jpj,jpk),     &
         & zsa_adin( jpi,jpj,jpk),     &
         & zta_adout(jpi,jpj,jpk),     &
         & zsa_adout(jpi,jpj,jpk),     &
         & ztb_adout(jpi,jpj,jpk),     &
         & zsb_adout(jpi,jpj,jpk),     &
         & zr(       jpi,jpj,jpk)      &
         & )
      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      ! initialization (normally done in traldf)
      l_traldf_rot = .TRUE.

      ! Test for time steps nit000 and nit000 + 1 (the matrix changes)

      DO jstp = nit000, nit000 + 2
         istp = jstp
         IF ( jstp == nit000+2 ) istp = nitend

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
          zta_tlin( :,:,:) = 0.0_wp
          ztb_tlin( :,:,:) = 0.0_wp
          zsa_tlin( :,:,:) = 0.0_wp
          zsb_tlin( :,:,:) = 0.0_wp
          zta_tlout(:,:,:) = 0.0_wp
          zsa_tlout(:,:,:) = 0.0_wp
          zta_adin( :,:,:) = 0.0_wp
          zsa_adin( :,:,:) = 0.0_wp
          zta_adout(:,:,:) = 0.0_wp
          zsa_adout(:,:,:) = 0.0_wp
          ztb_adout(:,:,:) = 0.0_wp
          zsb_adout(:,:,:) = 0.0_wp
          zr(       :,:,:) = 0.0_wp
          tsb_ad(:,:,:,:)     = 0.0_wp
          tsb_ad(:,:,:,:)     = 0.0_wp

          r2dtra(:) =  2.* rdttra(:)
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------

      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stdt )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              ztb_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsa_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO
      CALL grid_random(  zr, 'T', 0.0_wp, stds )
      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zsb_tlin(ji,jj,jk) = zr(ji,jj,jk)
            END DO
         END DO
      END DO


      tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
      tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)
      tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
      tsb_tl(:,:,:,jp_sal) = zsb_tlin(:,:,:)
      CALL tra_zdf_imp_tan ( istp, nit000, 'TRA', r2dtra, tsb_tl, tsa_tl, jpts )
      zta_tlout(:,:,:) = tsa_tl(:,:,:,jp_tem)
      zsa_tlout(:,:,:) = tsa_tl(:,:,:,jp_sal)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jk = 1, jpk
        DO jj = nldj, nlej
           DO ji = nldi, nlei
              zta_adin(ji,jj,jk) = zta_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_t(jk)
              zsa_adin(ji,jj,jk) = zsa_tlout(ji,jj,jk) &
                 &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                 &               * tmask(ji,jj,jk) * wesp_s(jk)
            END DO
         END DO
      END DO
      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------

      zsp1 = DOT_PRODUCT( zta_tlout, zta_adin ) &
         & + DOT_PRODUCT( zsa_tlout, zsa_adin )

      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------

      tsa_ad(:,:,:,jp_tem) = zta_adin(:,:,:)
      tsa_ad(:,:,:,jp_sal) = zsa_adin(:,:,:)

      CALL tra_zdf_imp_adj ( istp, nit000, 'TRA', r2dtra, tsb_ad, tsa_ad, jpts )

      zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)
      zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)
      ztb_adout(:,:,:) = tsb_ad(:,:,:,jp_tem)
      zsb_adout(:,:,:) = tsb_ad(:,:,:,jp_sal)
      zsp2 = DOT_PRODUCT( zta_tlin, zta_adout ) &
         & + DOT_PRODUCT( zsa_tlin, zsa_adout ) &
         & + DOT_PRODUCT( ztb_tlin, ztb_adout ) &
         & + DOT_PRODUCT( zsb_tlin, zsb_adout )

      ! 14 char:'12345678901234'
      IF ( istp == nit000 ) THEN
         cl_name = 'trazdfimpadjT1'
      ELSEIF ( istp == nit000 +1 ) THEN
         cl_name = 'trazdfimpadjT2'
      ELSEIF ( istp == nitend ) THEN
         cl_name = 'trazdfimpadjT3'
      END IF
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      END DO

      DEALLOCATE(   &
         & zta_tlin,  &
         & ztb_tlin,  &
         & zsa_tlin,  &
         & zsb_tlin,  &
         & zta_tlout, &
         & zsa_tlout, &
         & zta_adin,  &
         & zsa_adin,  &
         & zta_adout, &
         & ztb_adout, &
         & zsa_adout, &
         & zsb_adout, &
         & zr       &
         & )



   END SUBROUTINE tra_zdf_imp_adj_tst
END MODULE trazdf_imp_tam
