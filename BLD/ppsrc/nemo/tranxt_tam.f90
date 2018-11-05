MODULE tranxt_tam
   !!======================================================================
   !!                       ***  MODULE  tranxt_tam  ***
   !! Ocean active tracers:  time stepping on temperature and salinity
   !!               Tangent and Adjoint module
   !!======================================================================
   !! History of the direct module:
   !!            7.0  !  91-11  (G. Madec)  Original code
   !!                 !  93-03  (M. Guyon)  symetrical conditions
   !!                 !  96-02  (G. Madec & M. Imbard)  opa release 8.0
   !!            8.0  !  96-04  (A. Weaver)  Euler forward step
   !!            8.2  !  99-02  (G. Madec, N. Grima)  semi-implicit pressure grad.
   !!  NEMO      1.0  !  2002-08  (G. Madec)  F90: Free form and module
   !!             -   !  2002-11  (C. Talandier, A-M Treguier) Open boundaries
   !!             -   !  2005-04  (C. Deltel) Add Asselin trend in the ML budget
   !!            2.0  !  2006-02  (L. Debreu, C. Mazauric) Agrif implementation
   !!            3.0  !  2008-06  (G. Madec)  time stepping always done in trazd
   !!            3.1  !  2009-02  (G. Madec, R. Benshila)  re-introduce the vvl option
   !! History of the TAM module:
   !!            2.0  !  2008-09  (A. Vidard) tam of the 2006-02 version
   !!            3.0  !  2008-11  (A. Vidard) tam of the 2008-06 version
   !!             -   !  2009-01  (A. Weaver) corrections to test
   !!            3.2  !  2010-04  (F. Vigilant)  version 3.2
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   tra_nxt_tan    : time stepping on temperature and salinity (tangent)
   !!   tra_nxt_adj    : time stepping on temperature and salinity (adjoint)
   !!----------------------------------------------------------------------
   USE par_kind
   USE par_oce
   USE dynhpg
   USE oce_tam
   USE zdf_oce
   USE dom_oce
   USE tranxt
   USE in_out_manager
   USE lbclnk
   USE lbclnk_tam
   USE gridrandom
   USE dotprodfld
   USE paresp
   USE tstool_tam                   !           salinity
   USE traqsr
   USE wrk_nemo
   USE timing

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC   tra_nxt_tan     ! routine called by step_tam.F90
   PUBLIC   tra_nxt_adj     ! routine called by step_tam.F90
   PUBLIC   tra_nxt_adj_tst ! routine called by tst.F90

   REAL(wp) :: rbcp
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

CONTAINS

   SUBROUTINE tra_nxt_tan( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tranxt_tan  ***
      !!
      !! ** Purpose of the direct routine:
      !!             Apply the boundary condition on the after temperature
      !!             and salinity fields, achieved the time stepping by adding
      !!             the Asselin filter on now fields and swapping the fields.
      !!
      !! ** Method  :   At this stage of the computation, ta and sa are the
      !!             after temperature and salinity as the time stepping has
      !!             been performed in trazdf_imp or trazdf_exp module.
      !!
      !!              - Apply lateral boundary conditions on (ta,sa)
      !!             at the local domain   boundaries through lbc_lnk call,
      !!             at the radiative open boundaries (lk_obc=T),
      !!             at the relaxed   open boundaries (lk_bdy=T), and
      !!             at the AGRIF zoom     boundaries (lk_agrif=T)
      !!
      !!              - Update lateral boundary conditions on AGRIF children
      !!             domains (lk_agrif=T)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      INTEGER             :: jn, jk
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('tra_nxt_tan')
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_tan : achieve the time stepping by Asselin filter and array swap'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
      ENDIF
      rbcp = 0.25 * (1. + atfp) * (1. + atfp) * ( 1. - atfp)

      ! Update after tracer on domain lateral boundaries
      !
      CALL lbc_lnk( tsa_tl(:,:,:,jp_tem), 'T', 1.0_wp )      ! local domain boundaries  (T-point, unchanged sign)
      CALL lbc_lnk( tsa_tl(:,:,:,jp_sal), 'T', 1.0_wp )
      !
      ! set time step size (Euler/Leapfrog)
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dtra(:) =     rdttra(:)      ! at nit000             (Euler)
      ELSEIF( kt <= nit000 + 1 )           THEN   ;   r2dtra(:) = 2.* rdttra(:)      ! at nit000 or nit000+1 (Leapfrog)
      ENDIF

      IF ( neuler == 0 .AND. kt == nit000 ) THEN
         DO jn = 1, jpts
            DO jk = 1, jpkm1
               tsn_tl(:,:,jk,jn) = tsa_tl(:,:,jk,jn)
            END DO
         END DO
      ELSE
         ! Leap-Frog + Asselin filter time stepping
         IF( lk_vvl )   THEN   ;   CALL tra_nxt_vvl_tan( kt )      ! variable volume level (vvl)
         ELSE                  ;   CALL tra_nxt_fix_tan( kt, nit000, 'TRA', tsb_tl, tsn_tl, tsa_tl, jpts )      ! fixed    volume level
         ENDIF
      END IF
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_nxt_tan')
      !
   END SUBROUTINE tra_nxt_tan

   SUBROUTINE tra_nxt_fix_tan( kt, kit000, cdtype, ptb_tl, ptn_tl, pta_tl, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_fix_tan  ***
      !!
      !! ** Purpose :   fixed volume: apply the Asselin time filter and
      !!                swap the tracer fields.
      !!
      !! ** Method  : - Apply a Asselin time filter on now fields.
      !!              - save in (ta,sa) an average over the three time levels
      !!             which will be used to compute rdn and thus the semi-implicit
      !!             hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = (ta+2tn+tb)/4       ln_dynhpg_imp = T
      !!             ztm = 0                   otherwise
      !!             tb  = tn + atfp*[ tb - 2 tn + ta ]
      !!             tn  = ta
      !!             ta  = ztm       (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      INTEGER         , INTENT(in   )                               ::   kit000      ! first time step index
      CHARACTER(len=3), INTENT(in   )                               ::   cdtype      ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   )                               ::   kjpt        ! number of tracers
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptb_tl      ! before tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptn_tl      ! now tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   pta_tl      ! tracer trend
      !!
      INTEGER  ::   ji, jj, jk, jn            ! dummy loop indices
      REAL(wp) ::   ztntl, ztdtl, ztn     ! temporary scalars
      LOGICAL  ::   ll_tra_hpg
      !!----------------------------------------------------------------------

      IF( kt == kit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_fix_tan : time stepping ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      ztntl = 0._wp
      ztdtl = 0._wp
      !
      IF( cdtype == 'TRA' )  THEN   ;   ll_tra_hpg = ln_dynhpg_imp    ! active  tracers case  and  semi-implicit hpg
      ELSE                          ;   ll_tra_hpg = .FALSE.          ! passive tracers case or NO semi-implicit hpg
      ENDIF
      !
      DO jn = 1, kjpt
         !
         DO jk = 1, jpkm1
            DO jj = 1, jpj
               DO ji = 1, jpi
                  ztntl = ptn_tl(ji,jj,jk,jn)
                  ztdtl = pta_tl(ji,jj,jk,jn) - 2. * ztntl + ptb_tl(ji,jj,jk,jn)      !  time laplacian on tracers
                  !
                  ptb_tl(ji,jj,jk,jn) = ztntl + atfp * ztdtl                          ! ptb <-- filtered ptn
                  ptn_tl(ji,jj,jk,jn) = pta_tl(ji,jj,jk,jn)                           ! ptn <-- pta
                  !
                  IF( ll_tra_hpg )   pta_tl(ji,jj,jk,jn) = ztntl + rbcp * ztdtl       ! pta <-- Brown & Campana average
               END DO
           END DO
         END DO
         !
      END DO
      !
   END SUBROUTINE tra_nxt_fix_tan

   SUBROUTINE tra_nxt_vvl_tan( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_vvl_tan  ***
      !!
      !! ** Purpose :   Time varying volume: apply the Asselin time filter
      !!                and swap the tracer fields.
      !!
      !! ** Method  : - Apply a thickness weighted Asselin time filter on now fields.
      !!              - save in (ta,sa) a thickness weighted average over the three
      !!             time levels which will be used to compute rdn and thus the semi-
      !!             implicit hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = (e3t_a*ta+2*e3t_n*tn+e3t_b*tb)   ln_dynhpg_imp = T
      !!                  /(e3t_a   +2*e3t_n   +e3t_b   )
      !!             ztm = 0                                otherwise
      !!             tb  = ( e3t_n*tn + atfp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )
      !!                  /( e3t_n    + atfp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )
      !!             tn  = ta
      !!             ta  = zt        (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      REAL(wp) ::   ztm , ztc_f , ztf , ztca, ztcn, ztcb   ! temporary scalar
      REAL(wp) ::   zsm , zsc_f , zsf , zsca, zscn, zscb   !    -         -
      REAL(wp) ::   ze3mr, ze3fr                           !    -         -
      REAL(wp) ::   ze3t_b, ze3t_n, ze3t_a, ze3t_f         !    -         -
      !!----------------------------------------------------------------------

      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_vvl_tan : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      IF(lwp) WRITE(numout,*) "key_vvl net available in tangent yet"
      CALL abort
      !
   END SUBROUTINE tra_nxt_vvl_tan

   SUBROUTINE tra_nxt_adj( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tranxt_adj  ***
      !!
      !! ** Purpose of the direct routine:
      !!             Apply the boundary condition on the after temperature
      !!             and salinity fields, achieved the time stepping by adding
      !!             the Asselin filter on now fields and swapping the fields.
      !!
      !! ** Method  :   At this stage of the computation, ta and sa are the
      !!             after temperature and salinity as the time stepping has
      !!             been performed in trazdf_imp or trazdf_exp module.
      !!
      !!              - Apply lateral boundary conditions on (ta,sa)
      !!             at the local domain   boundaries through lbc_lnk call,
      !!             at the radiative open boundaries (lk_obc=T),
      !!             at the relaxed   open boundaries (lk_bdy=T), and
      !!             at the AGRIF zoom     boundaries (lk_agrif=T)
      !!
      !!              - Update lateral boundary conditions on AGRIF children
      !!             domains (lk_agrif=T)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      !!
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      INTEGER             ::   jn, jk
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start( 'tra_nxt_adj')
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_adj : achieve the time stepping by Asselin filter and array swap'
         IF(lwp) WRITE(numout,*) '~~~~~~~'
         rbcp = 0.25 * (1. + atfp) * (1. + atfp) * ( 1. - atfp)
      ENDIF

      ! set time step size (Euler/Leapfrog)
      r2dtra(:) =  2.* rdttra(:) ! initialization
      IF( neuler == 0 .AND. kt == nit000 ) THEN   ;   r2dtra(:) =     rdttra(:)      ! at nit000             (Euler)
      ELSEIF( kt <= nit000 + 1 )           THEN   ;   r2dtra(:) = 2.* rdttra(:)      ! at nit000 or nit000+1 (Leapfrog)
      ENDIF


      IF( neuler == 0 .AND. kt == nit000 ) THEN
         DO jn = 1, jpts
            DO jk = jpkm1, 1, -1
               tsa_ad(:,:,jk,jn) = tsa_ad(:,:,jk,jn) + tsn_ad(:,:,jk,jn)
               tsn_ad(:,:,jk,jn) = 0._wp
            END DO
         END DO
      ELSE
         !! Leap-Frog + Asselin filter time stepping
         IF( lk_vvl )   THEN   ;   CALL tra_nxt_vvl_adj( kt )      ! variable volume level (vvl)
         ELSE                  ;   CALL tra_nxt_fix_adj( kt, nit000, 'TRA', tsb_ad, tsn_ad, tsa_ad, jpts )      ! fixed    volume level
         ENDIF
      ENDIF

      ! Update after tracer on domain lateral boundaries
      !
      CALL lbc_lnk_adj( tsa_ad(:,:,:,jp_sal), 'T', 1.0_wp )
      CALL lbc_lnk_adj( tsa_ad(:,:,:,jp_tem), 'T', 1.0_wp )      ! local domain boundaries  (T-point, unchanged sign)
      !
      !
      IF( nn_timing == 1 )  CALL timing_stop('tra_nxt_adj')
      !
   END SUBROUTINE tra_nxt_adj

   SUBROUTINE tra_nxt_fix_adj( kt, kit000, cdtype, ptb_ad, ptn_ad, pta_ad, kjpt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_fix_adj  ***
      !!
      !! ** Purpose :   fixed volume: apply the Asselin time filter and
      !!                swap the tracer fields.
      !!
      !! ** Method  : - Apply a Asselin time filter on now fields.
      !!              - save in (ta,sa) an average over the three time levels
      !!             which will be used to compute rdn and thus the semi-implicit
      !!             hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = (ta+2tn+tb)/4       ln_dynhpg_imp = T
      !!             ztm = 0                   otherwise
      !!             tb  = tn + atfp*[ tb - 2 tn + ta ]
      !!             tn  = ta
      !!             ta  = ztm       (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      INTEGER         , INTENT(in   )                               ::   kit000      ! first time step index
      CHARACTER(len=3), INTENT(in   )                               ::   cdtype      ! =TRA or TRC (tracer indicator)
      INTEGER         , INTENT(in   )                               ::   kjpt        ! number of tracers
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptb_ad      ! before tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   ptn_ad      ! now tracer fields
      REAL(wp)        , INTENT(inout), DIMENSION(jpi,jpj,jpk,kjpt)  ::   pta_ad      ! tracer trend
      !!
      INTEGER  ::   ji, jj, jk, jn            ! dummy loop indices
      REAL(wp) ::   ztnad, ztdad, ztn     ! temporary scalars
      LOGICAL  ::   ll_tra_hpg
      !!----------------------------------------------------------------------

      IF( kt == kit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_fix_adj : time stepping ', cdtype
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF
      !
      IF( cdtype == 'TRA' )  THEN   ;   ll_tra_hpg = ln_dynhpg_imp    ! active  tracers case  and  semi-implicit hpg
      ELSE                          ;   ll_tra_hpg = .FALSE.          ! passive tracers case or NO semi-implicit hpg
      ENDIF
      ztnad = 0._wp
      ztdad = 0._wp
      !
      DO jn = 1, kjpt
         !
         DO jk = jpkm1, 1, -1
            DO jj = jpj, 1, -1
               DO ji = jpi, 1, -1
                  IF( ll_tra_hpg ) THEN
                     ztnad = ztnad + pta_ad(ji,jj,jk,jn)
                     ztdad = ztdad + rbcp * pta_ad(ji,jj,jk,jn)
                     pta_ad(ji,jj,jk,jn) = 0._wp
                  END IF
                  pta_ad(ji,jj,jk,jn) = pta_ad(ji,jj,jk,jn) + ptn_ad(ji,jj,jk,jn)
                  ptn_ad(ji,jj,jk,jn) = 0._wp
                  ztdad = ztdad + atfp * ptb_ad(ji,jj,jk,jn)
                  ztnad = ztnad + ptb_ad(ji,jj,jk,jn)
                  ptb_ad(ji,jj,jk,jn) = 0._wp
                  ptb_ad(ji,jj,jk,jn) = ztdad
                  pta_ad(ji,jj,jk,jn) = pta_ad(ji,jj,jk,jn) + ztdad
                  ztnad =  ztnad - 2._wp * ztdad
                  ptn_ad(ji,jj,jk,jn) = ptn_ad(ji,jj,jk,jn) + ztnad
                  ztdad = 0._wp
                  ztnad = 0._wp
               END DO
           END DO
         END DO
         !
      END DO
      !
   END SUBROUTINE tra_nxt_fix_adj

   SUBROUTINE tra_nxt_vvl_adj( kt )
      !!----------------------------------------------------------------------
      !!                   ***  ROUTINE tra_nxt_vvl_adj  ***
      !!
      !! ** Purpose :   Time varying volume: apply the Asselin time filter
      !!                and swap the tracer fields.
      !!
      !! ** Method  : - Apply a thickness weighted Asselin time filter on now fields.
      !!              - save in (ta,sa) a thickness weighted average over the three
      !!             time levels which will be used to compute rdn and thus the semi-
      !!             implicit hydrostatic pressure gradient (ln_dynhpg_imp = T)
      !!              - swap tracer fields to prepare the next time_step.
      !!                This can be summurized for tempearture as:
      !!             ztm = (e3t_a*ta+2*e3t_n*tn+e3t_b*tb)   ln_dynhpg_imp = T
      !!                  /(e3t_a   +2*e3t_n   +e3t_b   )
      !!             ztm = 0                                otherwise
      !!             tb  = ( e3t_n*tn + atfp*[ e3t_b*tb - 2 e3t_n*tn + e3t_a*ta ] )
      !!                  /( e3t_n    + atfp*[ e3t_b    - 2 e3t_n    + e3t_a    ] )
      !!             tn  = ta
      !!             ta  = zt        (NB: reset to 0 after eos_bn2 call)
      !!
      !! ** Action  : - (tb,sb) and (tn,sn) ready for the next time step
      !!              - (ta,sa) time averaged (t,s)   (ln_dynhpg_imp = T)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt    ! ocean time-step index
      !!
      INTEGER  ::   ji, jj, jk             ! dummy loop indices
      REAL(wp) ::   ztm , ztc_f , ztf , ztca, ztcn, ztcb   ! temporary scalar
      REAL(wp) ::   zsm , zsc_f , zsf , zsca, zscn, zscb   !    -         -
      REAL(wp) ::   ze3mr, ze3fr                           !    -         -
      REAL(wp) ::   ze3t_b, ze3t_n, ze3t_a, ze3t_f         !    -         -
      !!----------------------------------------------------------------------

      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_nxt_vvl_adj : time stepping'
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~~'
      ENDIF

      IF(lwp) WRITE(numout,*) "key_vvl net available in tangent yet"
      CALL abort
      !
   END SUBROUTINE tra_nxt_vvl_adj

   SUBROUTINE tra_nxt_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!          ***  ROUTINE tra_nxt_adj_tst : TEST OF tra_nxt_adj  ***
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
         & kumadt            ! Output unit

      INTEGER ::  &
         & ji,    &          ! dummy loop indices
         & jj,    &
         & jk,    &
         & jn

      LOGICAL ::  &          ! local variable for time scheme
         & ll_dynhpg_imp

      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d        ! 2D seed for the random number generator

      !! * Local declarations
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zsb_tlin,     &! Tangent input : before salinity
         & ztb_tlin,     &! Tangent input : before temperature
         & zsa_tlin,     &! Tangent input : after salinity
         & zta_tlin,     &! Tangent input : after temperature
         & zsn_tlin,     &! Tangent input : now salinity
         & ztn_tlin,     &! Tangent input : now temperature
         & zsb_tlout,    &! Tangent output: before salinity
         & ztb_tlout,    &! Tangent output: before temperature
         & zsa_tlout,    &! Tangent output: after salinity
         & zta_tlout,    &! Tangent output: after temperature
         & zsn_tlout,    &! Tangent output: now salinity
         & ztn_tlout,    &! Tangent output: now temperature
         & zsb_adin,     &! Adjoint input : before salinity
         & ztb_adin,     &! Adjoint input : before temperature
         & zsa_adin,     &! Adjoint input : after salinity
         & zta_adin,     &! Adjoint input : after temperature
         & zsn_adin,     &! Adjoint input : now salinity
         & ztn_adin,     &! Adjoint input : now temperature
         & zsb_adout,    &! Adjoint output: before salinity
         & ztb_adout,    &! Adjoint output: before temperature
         & zsa_adout,    &! Adjoint output: after salinity
         & zta_adout,    &! Adjoint output: after temperature
         & zsn_adout,    &! Adjoint output: now salinity
         & ztn_adout,    &! Adjoint output: now temperature
         & zr             ! 3D field

      REAL(KIND=wp) ::       &
         & zsp1,             & ! scalar product involving the tangent routine
         & zsp1_1,           & ! scalar product involving the tangent routine
         & zsp1_2,           & ! scalar product involving the tangent routine
         & zsp1_3,           & ! scalar product involving the tangent routine
         & zsp1_4,           & ! scalar product involving the tangent routine
         & zsp1_5,           & ! scalar product involving the tangent routine
         & zsp1_6,           & ! scalar product involving the tangent routine
         & zsp2,             & ! scalar product involving the adjoint routine
         & zsp2_1,           & ! scalar product involving the adjoint routine
         & zsp2_2,           & ! scalar product involving the adjoint routine
         & zsp2_3,           & ! scalar product involving the adjoint routine
         & zsp2_4,           & ! scalar product involving the adjoint routine
         & zsp2_5,           & ! scalar product involving the adjoint routine
         & zsp2_6              ! scalar product involving the adjoint routine
      CHARACTER(LEN=14) ::   &
         & cl_name

      ALLOCATE( &
         & zsb_tlin(jpi,jpj,jpk),    &! Tangent input : before salinity
         & ztb_tlin(jpi,jpj,jpk),    &! Tangent input : before temperature
         & zsa_tlin(jpi,jpj,jpk),    &! Tangent input : after salinity
         & zta_tlin(jpi,jpj,jpk),    &! Tangent input : after temperature
         & zsn_tlin(jpi,jpj,jpk),    &! Tangent input : now salinity
         & ztn_tlin(jpi,jpj,jpk),    &! Tangent input : now temperature
         & zsb_tlout(jpi,jpj,jpk),   &! Tangent output: before salinity
         & ztb_tlout(jpi,jpj,jpk),   &! Tangent output: before temperature
         & zsa_tlout(jpi,jpj,jpk),   &! Tangent output: after salinity
         & zta_tlout(jpi,jpj,jpk),   &! Tangent output: after temperature
         & zsn_tlout(jpi,jpj,jpk),   &! Tangent output: now salinity
         & ztn_tlout(jpi,jpj,jpk),   &! Tangent output: now temperature
         & zsb_adin(jpi,jpj,jpk),    &! Adjoint input : before salinity
         & ztb_adin(jpi,jpj,jpk),    &! Adjoint input : before temperature
         & zsa_adin(jpi,jpj,jpk),    &! Adjoint input : after salinity
         & zta_adin(jpi,jpj,jpk),    &! Adjoint input : after temperature
         & zsn_adin(jpi,jpj,jpk),    &! Adjoint input : now salinity
         & ztn_adin(jpi,jpj,jpk),    &! Adjoint input : now temperature
         & zsb_adout(jpi,jpj,jpk),   &! Adjoint output: before salinity
         & ztb_adout(jpi,jpj,jpk),   &! Adjoint output: before temperature
         & zsa_adout(jpi,jpj,jpk),   &! Adjoint output: after salinity
         & zta_adout(jpi,jpj,jpk),   &! Adjoint output: after temperature
         & zsn_adout(jpi,jpj,jpk),   &! Adjoint output: now salinity
         & ztn_adout(jpi,jpj,jpk),   &! Adjoint output: now temperature
         & zr       (jpi,jpj,jpk)    &! 3D field
         & )

      ll_dynhpg_imp = ln_dynhpg_imp     ! store namelist define time scheme

      DO jn = 1, 2

         IF ( jn .EQ. 1) ln_dynhpg_imp = .TRUE.
         IF ( jn .EQ. 2) ln_dynhpg_imp = .FALSE.

         !==================================================================
         ! 1) dx = ( tb_tl, tn_tl, ta_tl,       dy = ( tb_tl, tn_tl, ta_tl,
         !           sb_tl, sn_tl, sa_tl  ) and        sb_tl, sn_tl, sa_tl )
         !==================================================================

         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------
         tsb_tl(:,:,:,:) = 0.0_wp
         tsa_tl(:,:,:,:) = 0.0_wp
         tsn_tl(:,:,:,:) = 0.0_wp
         tsb_ad(:,:,:,:) = 0.0_wp
         tsa_ad(:,:,:,:) = 0.0_wp
         tsn_ad(:,:,:,:) = 0.0_wp
         zsb_tlin(:,:,:) = 0.0_wp
         ztb_tlin(:,:,:) = 0.0_wp
         zsa_tlin(:,:,:) = 0.0_wp
         zta_tlin(:,:,:) = 0.0_wp
         zsn_tlin(:,:,:) = 0.0_wp
         ztn_tlin(:,:,:) = 0.0_wp

         r2dtra(:) =  2.* rdttra(:) ! initialization

         CALL grid_random(  zr, 'T', 0.0_wp, stds )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zsb_tlin(ji,jj,jk) = zr(ji,jj,jk)
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

         CALL grid_random(  zr, 'T', 0.0_wp, stdt )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zta_tlin(ji,jj,jk) = zr(ji,jj,jk)
               END DO
            END DO
         END DO

         CALL grid_random(  zr, 'T', 0.0_wp, stds )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zsn_tlin(ji,jj,jk) = zr(ji,jj,jk)
               END DO
            END DO
         END DO

         CALL grid_random(  zr, 'T', 0.0_wp, stdt )
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  ztn_tlin(ji,jj,jk) = zr(ji,jj,jk)
               END DO
            END DO
         END DO

         tsb_tl(:,:,:,jp_sal) = zsb_tlin(:,:,:)
         tsb_tl(:,:,:,jp_tem) = ztb_tlin(:,:,:)
         tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)
         tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
         tsn_tl(:,:,:,jp_sal) = zsn_tlin(:,:,:)
         tsn_tl(:,:,:,jp_tem) = ztn_tlin(:,:,:)

         CALL tra_nxt_tan( nit000 + 1 )

         zsa_tlout(:,:,:) = tsa_tl(:,:,:,jp_sal)
         zta_tlout(:,:,:) = tsa_tl(:,:,:,jp_tem)
         zsb_tlout(:,:,:) = tsb_tl(:,:,:,jp_sal)
         ztb_tlout(:,:,:) = tsb_tl(:,:,:,jp_tem)
         zsn_tlout(:,:,:) = tsn_tl(:,:,:,jp_sal)
         ztn_tlout(:,:,:) = tsn_tl(:,:,:,jp_tem)

         !--------------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !--------------------------------------------------------------------

         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zsa_adin(ji,jj,jk)     = zsa_tlout(ji,jj,jk) &
                       &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                       &                   * tmask(ji,jj,jk) * wesp_s(jk)
                  zta_adin(ji,jj,jk)     = zta_tlout(ji,jj,jk) &
                       &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                       &                   * tmask(ji,jj,jk) * wesp_t(jk)
                  zsb_adin(ji,jj,jk)     = zsb_tlout(ji,jj,jk) &
                       &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                       &                   * tmask(ji,jj,jk) * wesp_s(jk)
                  ztb_adin(ji,jj,jk)     = ztb_tlout(ji,jj,jk) &
                       &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                       &                   * tmask(ji,jj,jk) * wesp_t(jk)
                  zsn_adin(ji,jj,jk)     = zsn_tlout(ji,jj,jk) &
                       &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                       &                   * tmask(ji,jj,jk) * wesp_s(jk)
                  ztn_adin(ji,jj,jk)     = ztn_tlout(ji,jj,jk) &
                       &                   * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                       &                   * tmask(ji,jj,jk) * wesp_t(jk)
               END DO
            END DO
         END DO

         !--------------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !--------------------------------------------------------------------

         zsp1_1 = DOT_PRODUCT( zsa_tlout    , zsa_adin     )
         zsp1_2 = DOT_PRODUCT( zta_tlout    , zta_adin     )
         zsp1_3 = DOT_PRODUCT( zsb_tlout    , zsb_adin     )
         zsp1_4 = DOT_PRODUCT( ztb_tlout    , ztb_adin     )
         zsp1_5 = DOT_PRODUCT( zsn_tlout    , zsn_adin     )
         zsp1_6 = DOT_PRODUCT( ztn_tlout    , ztn_adin     )
         zsp1   = zsp1_1 + zsp1_2 + zsp1_3 + zsp1_4 + zsp1_5 + zsp1_6

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------

         tsa_ad(:,:,:,jp_sal)     = zsa_adin(:,:,:)
         tsa_ad(:,:,:,jp_tem)     = zta_adin(:,:,:)
         tsb_ad(:,:,:,jp_sal)     = zsb_adin(:,:,:)
         tsb_ad(:,:,:,jp_tem)     = ztb_adin(:,:,:)
         tsn_ad(:,:,:,jp_sal)     = zsn_adin(:,:,:)
         tsn_ad(:,:,:,jp_tem)     = ztn_adin(:,:,:)

         CALL tra_nxt_adj ( nit000 + 1 )

         zsb_adout(:,:,:) = tsb_ad(:,:,:,jp_sal)
         ztb_adout(:,:,:) = tsb_ad(:,:,:,jp_tem)
         zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)
         zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)
         zsn_adout(:,:,:) = tsn_ad(:,:,:,jp_sal)
         ztn_adout(:,:,:) = tsn_ad(:,:,:,jp_tem)

         !--------------------------------------------------------------------
         ! Compute the scalar product: dx^T L^T W dy
         !--------------------------------------------------------------------

         zsp2_1 = DOT_PRODUCT( zsb_tlin  , zsb_adout   )
         zsp2_2 = DOT_PRODUCT( ztb_tlin  , ztb_adout   )
         zsp2_3 = DOT_PRODUCT( zsa_tlin  , zsa_adout   )
         zsp2_4 = DOT_PRODUCT( zta_tlin  , zta_adout   )
         zsp2_5 = DOT_PRODUCT( zsn_tlin  , zsn_adout   )
         zsp2_6 = DOT_PRODUCT( ztn_tlin  , ztn_adout   )

         zsp2 = zsp2_1 + zsp2_2 + zsp2_3 + zsp2_4 + zsp2_5 + zsp2_6

         ! Compare the scalar products

         ! 14 char:'12345678901234'
         IF ( jn .EQ. 1) cl_name = 'tra_nxt_adj T1'
         IF ( jn .EQ. 2) cl_name = 'tra_nxt_adj T2'
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

      ENDDO

      ln_dynhpg_imp = ll_dynhpg_imp   ! restore initial value of ln_dynhpg_imp

      DEALLOCATE( &
         & zsb_tlin,     &
         & ztb_tlin,     &
         & zsa_tlin,     &
         & zta_tlin,     &
         & zsn_tlin,     &
         & ztn_tlin,     &
         & zsb_tlout,    &
         & ztb_tlout,    &
         & zsa_tlout,    &
         & zta_tlout,    &
         & zsn_tlout,    &
         & ztn_tlout,    &
         & zsb_adin,     &
         & ztb_adin,     &
         & zsa_adin,     &
         & zta_adin,     &
         & zsn_adin,     &
         & ztn_adin,     &
         & zsb_adout,    &
         & ztb_adout,    &
         & zsa_adout,    &
         & zta_adout,    &
         & zsn_adout,    &
         & ztn_adout,    &
         & zr            &
         & )

     END SUBROUTINE tra_nxt_adj_tst
   !!======================================================================
END MODULE tranxt_tam
