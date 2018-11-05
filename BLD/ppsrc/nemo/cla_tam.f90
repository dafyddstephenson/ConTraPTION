MODULE cla_tam
   !!======================================================================
   !!                    ***  MODULE  cla  ***
   !! Cross Land Advection : specific update of the horizontal divergence,
   !!                        tracer trends and after velocity
   !!
   !!                 ---   Specific to ORCA_R2   ---
   !!
   !!======================================================================
   !! History :  1.0  ! 2002-11 (A. Bozec)  Original code
   !!            3.2  ! 2009-07 (G. Madec)  merge cla, cla_div, tra_cla, cla_dynspg
   !!                 !                     and correct a mpp bug reported by A.R. Porter
   !! History of TAM :
   !!            3.4  ! 2012-07 (P.-A. Bouttier)
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_orca_r2'                                 global ocean model R2
   !!----------------------------------------------------------------------
   !!   cla_div           : update of horizontal divergence at cla straits
   !!   tra_cla           : update of tracers at cla straits
   !!   cla_dynspg        : update of after horizontal velocities at cla straits
   !!   cla_init          : initialisation - control check
   !!   cla_bab_el_mandeb : cross land advection for Bab-el-mandeb strait
   !!   cla_gibraltar     : cross land advection for Gibraltar strait
   !!   cla_hormuz        : cross land advection for Hormuz strait
   !!----------------------------------------------------------------------
   USE oce            ! ocean dynamics and tracers
   USE dom_oce        ! ocean space and time domain
   USE sbc_oce        ! surface boundary condition: ocean
   USE dynspg_oce     ! ocean dynamics: surface pressure gradient variables
   USE oce_tam        ! ocean dynamics and tracers
   USE sbc_oce_tam    ! surface boundary condition: ocean
   !USE dynspg_oce_tam ! ocean dynamics: surface pressure gradient variables
   USE in_out_manager ! I/O manager
   USE lib_mpp        ! distributed memory computing library
   USE lbclnk         ! ocean lateral boundary conditions (or mpp link)
   USE dotprodfld
   USE tstool_tam
   USE paresp
   USE gridrandom

   IMPLICIT NONE
   PRIVATE

   PUBLIC   cla_init_tam         ! routine called by opatam.F90
   PUBLIC   cla_div_tan          ! routine called by divcur_tan.F90
   PUBLIC   cla_traadv_tan       ! routine called by traadv_tan.F90
   PUBLIC   cla_dynspg_tan       ! routine called by dynspg_flt_tan.F90
   PUBLIC   cla_div_adj          ! routine called by divcur_adj.F90
   PUBLIC   cla_traadv_adj       ! routine called by traadv_adj.F90
   PUBLIC   cla_dynspg_adj       ! routine called by dynspg_flt_adj.F90
   PUBLIC   cla_div_adj_tst      ! routine called by tamtst.F90
   PUBLIC   cla_traadv_adj_tst   ! routine called by traadv_adj.F90
   PUBLIC   cla_dynspg_adj_tst   ! routine called by dynspg_flt_adj.F90

   INTEGER  ::   nbab, ngib, nhor   ! presence or not of required grid-point on local domain
   !                                ! for Bab-el-Mandeb, Gibraltar, and Hormuz straits

   !                                           !   fixed part  !  time evolving    !!! profile of hdiv for some straits
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_139_101, hdiv_139_101_kt    ! Gibraltar     (i,j)=(172,101)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_139_102                     ! Gibraltar     (i,j)=(139,102)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_141_102, hdiv_141_102_kt    ! Gibraltar     (i,j)=(141,102)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_161_88 , hdiv_161_88_kt     ! Bab-el-Mandeb (i,j)=(161,88)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_161_87                      ! Bab-el-Mandeb (i,j)=(161,87)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_160_89 , hdiv_160_89_kt     ! Bab-el-Mandeb (i,j)=(160,89)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_172_94                      ! Hormuz        (i,j)=(172, 94)

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   t_171_94_hor, s_171_94_hor       ! Temperature, salinity in Hormuz strait

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_139_101_tl, hdiv_139_101_kt_tl    ! Gibraltar     (i,j)=(172,101)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_139_102_tl                        ! Gibraltar     (i,j)=(139,102)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_141_102_tl, hdiv_141_102_kt_tl    ! Gibraltar     (i,j)=(141,102)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_161_88_tl , hdiv_161_88_kt_tl     ! Bab-el-Mandeb (i,j)=(161,88)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_161_87_tl                         ! Bab-el-Mandeb (i,j)=(161,87)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_160_89_tl , hdiv_160_89_kt_tl     ! Bab-el-Mandeb (i,j)=(160,89)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_172_94_tl                         ! Hormuz        (i,j)=(172, 94)

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   t_171_94_hor_tl, s_171_94_hor_tl       ! Temperature, salinity in Hormuz strait

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_139_101_ad, hdiv_139_101_kt_ad    ! Gibraltar     (i,j)=(172,101)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_139_102_ad                         ! Gibraltar     (i,j)=(139,102)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_141_102_ad, hdiv_141_102_kt_ad    ! Gibraltar     (i,j)=(141,102)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_161_88_ad , hdiv_161_88_kt_ad     ! Bab-el-Mandeb (i,j)=(161,88)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_161_87_ad                          ! Bab-el-Mandeb (i,j)=(161,87)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_160_89_ad , hdiv_160_89_kt_ad     ! Bab-el-Mandeb (i,j)=(160,89)
   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   hdiv_172_94_ad                          ! Hormuz        (i,j)=(172, 94)

   REAL(wp), ALLOCATABLE, SAVE, DIMENSION (:) ::   t_171_94_hor_ad, s_171_94_hor_ad       ! Temperature, salinity in Hormuz strait
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
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE cla_div_tan( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE div_cla_tan  ***
      !!
      !! ** Purpose :   update the horizontal divergence of the velocity field
      !!              at some straits ( Gibraltar, Bab el Mandeb and Hormuz ).
      !!
      !! ** Method  : - first time-step: initialisation of cla
      !!              - all   time-step: using imposed transport at each strait,
      !!              the now horizontal divergence is updated
      !!
      !! ** Action  :   phdivn   updted now horizontal divergence at cla straits
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cla_tan : cross land advection on hdiv '
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         !
         IF( nbab == 1 )   CALL cla_bab_el_mandeb_tan('ini')    ! Bab el Mandeb ( Red Sea - Indian ocean )
         IF( ngib == 1 )   CALL cla_gibraltar_tan    ('ini')    ! Gibraltar strait (Med Sea - Atlantic ocean)
         IF( nhor == 1 )   CALL cla_hormuz_tan       ('ini')    ! Hormuz Strait ( Persian Gulf - Indian ocean )
         !
      ENDIF
      !
      IF( nbab == 1    )   CALL cla_bab_el_mandeb_tan('div')    ! Bab el Mandeb ( Red Sea - Indian ocean )
      IF( ngib == 1    )   CALL cla_gibraltar_tan    ('div')    ! Gibraltar strait (Med Sea - Atlantic ocean)
      IF( nhor == 1    )   CALL cla_hormuz_tan       ('div')    ! Hormuz Strait ( Persian Gulf - Indian ocean )
      !
!!gm  lbc useless here, no?
!!gm      CALL lbc_lnk( hdivn, 'T', 1. )                    ! Lateral boundary conditions on hdivn
      !
   END SUBROUTINE cla_div_tan

   SUBROUTINE cla_div_adj( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE div_cla_adj  ***
      !!
      !! ** Purpose :   update the horizontal divergence of the velocity field
      !!              at some straits ( Gibraltar, Bab el Mandeb and Hormuz ).
      !!
      !! ** Method  : - first time-step: initialisation of cla
      !!              - all   time-step: using imposed transport at each strait,
      !!              the now horizontal divergence is updated
      !!
      !! ** Action  :   phdivn   updted now horizontal divergence at cla straits
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( kt == nitend ) THEN
         !
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'div_cla_adj : cross land advection on hdiv '
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
         !
         IF( nbab == 1 )   CALL cla_bab_el_mandeb_adj('ini')    ! Bab el Mandeb ( Red Sea - Indian ocean )
         IF( ngib == 1 )   CALL cla_gibraltar_adj    ('ini')    ! Gibraltar strait (Med Sea - Atlantic ocean)
         IF( nhor == 1 )   CALL cla_hormuz_adj       ('ini')    ! Hormuz Strait ( Persian Gulf - Indian ocean )
         !
      ENDIF
      !
      IF( nbab == 1    )   CALL cla_bab_el_mandeb_adj('div')    ! Bab el Mandeb ( Red Sea - Indian ocean )
      IF( ngib == 1    )   CALL cla_gibraltar_adj    ('div')    ! Gibraltar strait (Med Sea - Atlantic ocean)
      IF( nhor == 1    )   CALL cla_hormuz_adj       ('div')    ! Hormuz Strait ( Persian Gulf - Indian ocean )
      !
!!gm  lbc useless here, no?
!!gm      CALL lbc_lnk( hdivn, 'T', 1. )                    ! Lateral boundary conditions on hdivn
      !
   END SUBROUTINE cla_div_adj

   SUBROUTINE cla_traadv_tan( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tra_cla_tan  ***
      !!
      !! ** Purpose :   Update the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations
      !!      at some straits ( Bab el Mandeb, Gibraltar, Hormuz ).
      !!
      !! ** Method  :   using both imposed transport at each strait and T & S
      !!              budget, the now tracer trends is updated
      !!
      !! ** Action  :   (ta,sa)   updated now tracer trends at cla straits
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_cla_tan : cross land advection on tracers '
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
      ENDIF
      !
      IF( nbab == 1    )   CALL cla_bab_el_mandeb_tan('tra')      ! Bab el Mandeb strait
      IF( ngib == 1    )   CALL cla_gibraltar_tan    ('tra')      ! Gibraltar strait
      IF( nhor == 1    )   CALL cla_hormuz_tan       ('tra')      ! Hormuz Strait ( Persian Gulf)
      !
   END SUBROUTINE cla_traadv_tan

   SUBROUTINE cla_traadv_adj( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE tra_cla_adj  ***
      !!
      !! ** Purpose :   Update the now trend due to the advection of tracers
      !!      and add it to the general trend of passive tracer equations
      !!      at some straits ( Bab el Mandeb, Gibraltar, Hormuz ).
      !!
      !! ** Method  :   using both imposed transport at each strait and T & S
      !!              budget, the now tracer trends is updated
      !!
      !! ** Action  :   (ta,sa)   updated now tracer trends at cla straits
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'tra_cla_adj : cross land advection on tracers '
         IF(lwp) WRITE(numout,*) '~~~~~~~~'
      ENDIF
      !
      IF( nbab == 1    )   CALL cla_bab_el_mandeb_adj('tra')      ! Bab el Mandeb strait
      IF( ngib == 1    )   CALL cla_gibraltar_adj    ('tra')      ! Gibraltar strait
      IF( nhor == 1    )   CALL cla_hormuz_adj       ('tra')      ! Hormuz Strait ( Persian Gulf)
      !
   END SUBROUTINE cla_traadv_adj

   SUBROUTINE cla_dynspg_tan( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE cla_dynspg_tan  ***
      !!
      !! ** Purpose :   Update the after velocity at some straits
      !!              (Bab el Mandeb, Gibraltar, Hormuz).
      !!
      !! ** Method  :   required to compute the filtered surface pressure gradient
      !!
      !! ** Action  :   (ua,va)   after velocity at the cla straits
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( kt == nit000 ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'cla_dynspg_tan : cross land advection on (ua,va) '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !
      IF( nbab == 1    )   CALL cla_bab_el_mandeb_tan('spg')      ! Bab el Mandeb strait
      IF( ngib == 1    )   CALL cla_gibraltar_tan    ('spg')      ! Gibraltar strait
      IF( nhor == 1    )   CALL cla_hormuz_tan       ('spg')      ! Hormuz Strait ( Persian Gulf)
      !
!!gm lbc is needed here, not?
!!gm      CALL lbc_lnk( hdivn, 'U', -1. )   ;   CALL lbc_lnk( hdivn, 'V', -1. )      ! Lateral boundary conditions
      !
   END SUBROUTINE cla_dynspg_tan

   SUBROUTINE cla_dynspg_adj( kt )
      !!----------------------------------------------------------------------
      !!                 ***  ROUTINE cla_dynspg_adj  ***
      !!
      !! ** Purpose :   Update the after velocity at some straits
      !!              (Bab el Mandeb, Gibraltar, Hormuz).
      !!
      !! ** Method  :   required to compute the filtered surface pressure gradient
      !!
      !! ** Action  :   (ua,va)   after velocity at the cla straits
      !!----------------------------------------------------------------------
      INTEGER, INTENT( in ) ::   kt         ! ocean time-step index
      !!----------------------------------------------------------------------
      !
      IF( kt == nitend ) THEN
         IF(lwp) WRITE(numout,*)
         IF(lwp) WRITE(numout,*) 'cla_dynspg_adj : cross land advection on (ua,va) '
         IF(lwp) WRITE(numout,*) '~~~~~~~~~~'
      ENDIF
      !
      IF( nbab == 1    )   CALL cla_bab_el_mandeb_adj('spg')      ! Bab el Mandeb strait
      IF( ngib == 1    )   CALL cla_gibraltar_adj    ('spg')      ! Gibraltar strait
      IF( nhor == 1    )   CALL cla_hormuz_adj       ('spg')      ! Hormuz Strait ( Persian Gulf)
      !
!!gm lbc is needed here, not?
!!gm      CALL lbc_lnk( hdivn, 'U', -1. )   ;   CALL lbc_lnk( hdivn, 'V', -1. )      ! Lateral boundary conditions
      !
   END SUBROUTINE cla_dynspg_adj

   SUBROUTINE cla_init_tam
      !! -------------------------------------------------------------------
      !!                   ***  ROUTINE cla_init_tam  ***
      !!
      !! ** Purpose :   control check for mpp computation
      !!
      !! ** Method  : - All the strait grid-points must be inside one of the
      !!              local domain interior for the cla advection to work
      !!              properly in mpp (i.e. inside (2:jpim1,2:jpjm1) ).
      !!              Define the corresponding indicators (nbab, ngib, nhor)
      !!              - The profiles of cross-land fluxes are currently hard
      !!              coded for L31 levels. Stop if jpk/=31
      !!
      !! ** Action  :   nbab, ngib, nhor   strait inside the local domain or not
      !!---------------------------------------------------------------------
      REAL(wp) ::   ztemp   ! local scalar
      INTEGER  ::   ierr    ! local integer
      !!---------------------------------------------------------------------
      !
      IF(lwp) WRITE(numout,*)
      IF(lwp) WRITE(numout,*) 'cla_init_tam : cross land advection initialisation '
      IF(lwp) WRITE(numout,*) '~~~~~~~~~'
      !
      !                           ! Allocate arrays for this module
      ALLOCATE( hdiv_139_101(jpk) , hdiv_139_101_kt(jpk) ,     &    ! Gibraltar
         &      hdiv_139_102(jpk) ,                            &
         &      hdiv_141_102(jpk) , hdiv_141_102_kt(jpk) ,     &
         &      hdiv_161_88 (jpk) , hdiv_161_88_kt (jpk) ,     &    ! Bab-el-Mandeb
         &      hdiv_161_87 (jpk) ,                            &
         &      hdiv_160_89 (jpk) , hdiv_160_89_kt (jpk) ,     &     ! Hormuz
         &      hdiv_172_94 (jpk) ,                            &
         &      t_171_94_hor(jpk) , s_171_94_hor   (jpk) , STAT=ierr )
      ALLOCATE( hdiv_139_101_tl(jpk) , hdiv_139_101_kt_tl(jpk) ,     &    ! Gibraltar
         &      hdiv_139_102_tl(jpk) ,                               &
         &      hdiv_141_102_tl(jpk) , hdiv_141_102_kt_tl(jpk) ,     &
         &      hdiv_161_88_tl (jpk) , hdiv_161_88_kt_tl (jpk) ,     &    ! Bab-el-Mandeb
         &      hdiv_161_87_tl (jpk) ,                               &
         &      hdiv_160_89_tl (jpk) , hdiv_160_89_kt_tl (jpk) ,     &     ! Hormuz
         &      hdiv_172_94_tl (jpk) ,                               &
         &      t_171_94_hor_tl(jpk) , s_171_94_hor_tl   (jpk) , STAT=ierr )
      ALLOCATE( hdiv_139_101_ad(jpk) , hdiv_139_101_kt_ad(jpk) ,     &    ! Gibraltar
         &      hdiv_139_102_ad(jpk) ,                               &
         &      hdiv_141_102_ad(jpk) , hdiv_141_102_kt_ad(jpk) ,     &
         &      hdiv_161_88_ad (jpk) , hdiv_161_88_kt_ad (jpk) ,     &    ! Bab-el-Mandeb
         &      hdiv_161_87_ad (jpk) ,                               &
         &      hdiv_160_89_ad (jpk) , hdiv_160_89_kt_ad (jpk) ,     &     ! Hormuz
         &      hdiv_172_94_ad (jpk) ,                               &
         &      t_171_94_hor_ad(jpk) , s_171_94_hor_ad   (jpk) , STAT=ierr )
      IF( lk_mpp    )   CALL mpp_sum( ierr )
      IF( ierr /= 0 )   CALL ctl_stop( 'STOP', 'cla_init_tam: unable to allocate arrays' )
      !
      IF( .NOT.lk_dynspg_flt )   CALL ctl_stop( 'cla_init_tam: Cross Land Advection works only with lk_dynspg_flt=T ' )
      !
      IF( lk_vvl             )   CALL ctl_stop( 'cla_init_tam: Cross Land Advection does not work with lk_vvl=T option' )
      !
      IF( jpk /= 31          )   CALL ctl_stop( 'cla_init_tam: Cross Land Advection hard coded for ORCA_R2_L31' )
      !
      !                                        _|_______|_______|_
      !                                     89  |       |///////|
      !                                        _|_______|_______|_
      ! ------------------------ !          88  |///////|       |
      !   Bab el Mandeb strait   !             _|_______|_______|_
      ! ------------------------ !          87  |///////|       |
      !                                        _|_______|_______|_
      !                                         |  160  |  161  |
      !
      ! The 6 Bab el Mandeb grid-points must be inside one of the interior of the
      ! local domain for the cla advection to work properly (i.e. (2:jpim1,2:jpjm1)
      nbab = 0
      IF(  ( 1 <= mj0( 88) .AND. mj1( 89) <= jpj ) .AND.    &  !* (161,89), (161,88) and (161,88) on the local pocessor
         & ( 1 <= mi0(160) .AND. mi1(161) <= jpi )       )    nbab = 1
      !
      ! test if there is no local domain that includes all required grid-points
      ztemp = REAL( nbab )
      IF( lk_mpp )   CALL mpp_sum( ztemp )      ! sum with other processors value
      IF( ztemp == 0 ) THEN                     ! Only 2 points in each direction, this should never be a problem
         CALL ctl_stop( ' cross land advection at Bab-el_Mandeb does not work with your processor cutting: change it' )
      ENDIF
      !                                        ___________________________
      ! ------------------------ !         102  |       |///////|       |
      !     Gibraltar strait     !             _|_______|_______|_______|_
      ! ------------------------ !         101  |       |///////|       |
      !                                        _|_______|_______|_______|_
      !                                         |  139  |  140  |  141  |
      !
      ! The 6 Gibraltar grid-points must be inside one of the interior of the
      ! local domain for the cla advection to work properly (i.e. (2:jpim1,2:jpjm1)
      ngib = 0
      IF(  ( 2 <= mj0(101) .AND. mj1(102) <= jpjm1 ) .AND.    &  !* (139:141,101:102) on the local pocessor
         & ( 2 <= mi0(139) .AND. mi1(141) <= jpim1 )       )    ngib = 1
      !
      ! test if there is no local domain that includes all required grid-points
      ztemp = REAL( ngib )
      IF( lk_mpp )   CALL mpp_sum( ztemp )      ! sum with other processors value
      IF( ztemp == 0 ) THEN                     ! 3 points in i-direction, this may be a problem with some cutting
           CALL ctl_stop( ' cross land advection at Gibraltar does not work with your processor cutting: change it' )
      ENDIF
      !                                        _______________
      ! ------------------------ !          94  |/////|     |
      !       Hormuz strait      !             _|_____|_____|_
      ! ------------------------ !                171   172
      !
      ! The 2 Hormuz grid-points must be inside one of the interior of the
      ! local domain for the cla advection to work properly (i.e. (2:jpim1,2:jpjm1)
      nhor = 0
      IF(    2 <= mj0( 94) .AND. mj1( 94) <= jpjm1  .AND.  &
         &   2 <= mi0(171) .AND. mi1(172) <= jpim1         )   nhor = 1
      !
      ! test if there is no local domain that includes all required grid-points
      ztemp = REAL( nhor )
      IF( lk_mpp )   CALL mpp_sum( ztemp )      ! sum with other processors value
      IF( ztemp == 0 ) THEN                     ! 3 points in i-direction, this may be a problem with some cutting
           CALL ctl_stop( ' cross land advection at Hormuz does not work with your processor cutting: change it' )
      ENDIF
      !
   END SUBROUTINE cla_init_tam

   SUBROUTINE cla_bab_el_mandeb_tan( cd_td )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE cla_bab_el_mandeb_tan  ***
      !!
      !! ** Purpose :   update the now horizontal divergence, the tracer tendancy
      !!              and the after velocity in vicinity of Bab el Mandeb ( Red Sea - Indian ocean).
      !!
      !! ** Method  :   compute the exchanges at each side of the strait :
      !!
      !!       surf. zio_flow
      !! (+ balance of emp) /\  |\\\\\\\\\\\|
      !!                    ||  |\\\\\\\\\\\|
      !!    deep zio_flow   ||  |\\\\\\\\\\\|
      !!            |  ||   ||  |\\\\\\\\\\\|
      !!        89  |  ||   ||  |\\\\\\\\\\\|
      !!            |__\/_v_||__|____________
      !!            !\\\\\\\\\\\|          surf. zio_flow
      !!            |\\\\\\\\\\\|<===    (+ balance of emp)
      !!            |\\\\\\\\\\\u
      !!        88  |\\\\\\\\\\\|<---      deep  zrecirc (upper+deep at 2 different levels)
      !!            |___________|__________
      !!            !\\\\\\\\\\\|
      !!            |\\\\\\\\\\\| ---\     deep  zrecirc (upper+deep)
      !!        87  !\\\\\\\\\\\u ===/   + deep  zio_flow   (all at the same level)
      !!            !\\\\\\\\\\\|
      !!            !___________|__________
      !!                160         161
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in) ::   cd_td   ! ='div' update the divergence
      !                                         ! ='tra' update the tracers
      !                                         ! ='spg' update after velocity
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zemp_red_tl     ! temporary scalar
      REAL(wp) ::   zio_flow, zrecirc_upp, zrecirc_mid, zrecirc_bot
      !!---------------------------------------------------------------------
      !
      SELECT CASE( cd_td )
      !                     ! ---------------- !
      CASE( 'ini' )         !  initialisation  !
         !                  ! ---------------- !
         !
         zio_flow    = 0.4e6                       ! imposed in/out flow
         zrecirc_upp = 0.2e6                       ! imposed upper recirculation water
         zrecirc_bot = 0.5e6                       ! imposed bottom  recirculation water

         hdiv_161_88(:) = 0.e0                     ! (161,88) Gulf of Aden side, north point
         hdiv_161_87(:) = 0.e0                     ! (161,87) Gulf of Aden side, south point
         hdiv_160_89(:) = 0.e0                     ! (160,89) Red sea side
         hdiv_161_88_tl(:) = 0.e0                     ! (161,88) Gulf of Aden side, north point
         hdiv_161_87_tl(:) = 0.e0                     ! (161,87) Gulf of Aden side, south point
         hdiv_160_89_tl(:) = 0.e0                     ! (160,89) Red sea side

         DO jj = mj0(88), mj1(88)              !** profile of hdiv at (161,88)   (Gulf of Aden side, north point)
            DO ji = mi0(161), mi1(161)         !------------------------------
               DO jk = 1, 8                        ! surface in/out flow   (Ind -> Red)   (div >0)
                  hdiv_161_88(jk) = + zio_flow / ( 8. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                   ! recirculation water   (Ind -> Red)   (div >0)
               hdiv_161_88(20) =                 + zrecirc_upp   / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,20) )
               hdiv_161_88(21) = + ( zrecirc_bot - zrecirc_upp ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,21) )
            END DO
         END DO
         !
         DO jj = mj0(87), mj1(87)              !** profile of hdiv at (161,88)   (Gulf of Aden side, south point)
            DO ji = mi0(161), mi1(161)         !------------------------------
               !                                   ! deep out flow + recirculation   (Red -> Ind)   (div <0)
               hdiv_161_87(21) = - ( zio_flow + zrecirc_bot ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,21) )
            END DO
         END DO
         !
         DO jj = mj0(89), mj1(89)              !** profile of hdiv at (161,88)   (Red sea side)
            DO ji = mi0(160), mi1(160)         !------------------------------
               DO jk = 1, 8                        ! surface inflow    (Ind -> Red)   (div <0)
                  hdiv_160_89(jk) = - zio_flow / ( 8. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                   ! deep    outflow   (Red -> Ind)   (div >0)
               hdiv_160_89(16)    = + zio_flow / (      e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,16) )
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'div' )         !   update hdivn   ! (call by divcur module)
         !                  ! ---------=====-- !
         !                                     !** emp on the Red Sea   (div >0)
         zemp_red_tl = 0.e0                       !---------------------
         DO jj = mj0(87), mj1(96)                  ! sum over the Red sea
            DO ji = mi0(148), mi1(160)
               zemp_red_tl = zemp_red_tl + emp_tl(ji,jj) * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum( zemp_red_tl )   ! sum with other processors value
         zemp_red_tl = zemp_red_tl * 1.e-3               ! convert in m3
         !
         !                                     !** Correct hdivn (including emp adjustment)
         !                                     !-------------------------------------------
         DO jj = mj0(88), mj1(88)                  !* profile of hdiv at (161,88)   (Gulf of Aden side, north point)
            DO ji = mi0(161), mi1(161)
               hdiv_161_88_kt_tl(:) = hdiv_161_88_tl(:)
               DO jk = 1, 8                              ! increase the inflow from the Indian   (div >0)
                  hdiv_161_88_kt_tl(jk) = hdiv_161_88_tl(jk) + zemp_red_tl / (8. * e2u(ji,jj) * e3u(ji,jj,jk) )
               END DO
               hdivn_tl(ji,jj,:) = hdivn_tl(ji,jj,:) + hdiv_161_88_kt_tl(:)
            END DO
         END DO
         DO jj = mj0(87), mj1(87)                  !* profile of divergence at (161,87)   (Gulf of Aden side, south point)
            DO ji = mi0(161), mi1(161)
               hdivn_tl(ji,jj,:) = hdivn_tl(ji,jj,:) + hdiv_161_87_tl(:)
            END DO
         END DO
         DO jj = mj0(89), mj1(89)                  !* profile of divergence at (160,89)   (Red sea side)
            DO ji = mi0(160), mi1(160)
               hdiv_160_89_kt_tl(:) = hdiv_160_89_tl(:)
               DO jk = 1, 18                              ! increase the inflow from the Indian   (div <0)
                  hdiv_160_89_kt_tl(jk) = hdiv_160_89_tl(jk) - zemp_red_tl / (10. * e1v(ji,jj) * e3v(ji,jj,jk) )
               END DO
               hdivn_tl(ji, jj,:) = hdivn_tl(ji, jj,:) + hdiv_160_89_kt_tl(:)
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'tra' )         !  update (ta,sa)  ! (call by traadv module)
         !                  ! --------=======- !
         !
         DO jj = mj0(88), mj1(88)              !** (161,88)   (Gulf of Aden side, north point)
            DO ji = mi0(161), mi1(161)
               DO jk = 1, jpkm1                         ! surf inflow + reciculation (from Gulf of Aden)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)                         &
                     &                       - hdiv_161_88_kt_tl(jk) * tsn(ji,jj,jk,jp_tem) &
                     &                       - hdiv_161_88_kt(jk) * tsn_tl(ji,jj,jk,jp_tem)
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal)                         &
                     &                       - hdiv_161_88_kt_tl(jk) * tsn(ji,jj,jk,jp_sal) &
                     &                       - hdiv_161_88_kt(jk) * tsn_tl(ji,jj,jk,jp_sal)
               END DO
            END DO
         END DO
         DO jj = mj0(87), mj1(87)              !** (161,87)   (Gulf of Aden side, south point)
            DO ji = mi0(161), mi1(161)
               jk =  21                                 ! deep outflow + recirulation (combined flux)
               tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)                             &
                  &                        + hdiv_161_88_tl(20) * tsn(ji  ,jj+1,20,jp_tem)   &  ! upper recirculation from Gulf of Aden
                  &                        + hdiv_161_88(20) * tsn_tl(ji  ,jj+1,20,jp_tem)   &  ! upper recirculation from Gulf of Aden
                  &                        + hdiv_161_88_tl(21) * tsn(ji  ,jj+1,21,jp_tem)   &  ! deep  recirculation from Gulf of Aden
                  &                        + hdiv_161_88(21) * tsn_tl(ji  ,jj+1,21,jp_tem)   &  ! deep  recirculation from Gulf of Aden
                  &                        + hdiv_160_89_tl(16) * tsn(ji-1,jj+2,16,jp_tem)   &  ! deep inflow from Red sea
                  &                        + hdiv_160_89(16) * tsn_tl(ji-1,jj+2,16,jp_tem)      ! deep inflow from Red sea
               tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal) &
                  &                        + hdiv_161_88_tl(20) * tsn(ji  ,jj+1,20,jp_sal)   &
                  &                        + hdiv_161_88(20) * tsn_tl(ji  ,jj+1,20,jp_sal)   &
                  &                        + hdiv_161_88_tl(21) * tsn(ji  ,jj+1,21,jp_sal)   &
                  &                        + hdiv_161_88(21) * tsn_tl(ji  ,jj+1,21,jp_sal)   &
                  &                        + hdiv_160_89_tl(16) * tsn(ji-1,jj+2,16,jp_sal)   &
                  &                        + hdiv_160_89(16) * tsn_tl(ji-1,jj+2,16,jp_sal)
            END DO
         END DO
         DO jj = mj0(89), mj1(89)              !** (161,88)   (Red sea side)
            DO ji = mi0(160), mi1(160)
               DO jk = 1, 14                            ! surface inflow (from Gulf of Aden)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)                                &
                     &                    - hdiv_160_89_kt_tl(jk) * tsn(ji+1,jj-1,jk,jp_tem)       &
                     &                    - hdiv_160_89_kt(jk) * tsn_tl(ji+1,jj-1,jk,jp_tem)
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal)                                &
                     &                       - hdiv_160_89_kt_tl(jk) * tsn(ji+1,jj-1,jk,jp_sal)    &
                     &                       - hdiv_160_89_kt(jk) * tsn_tl(ji+1,jj-1,jk,jp_sal)
               END DO
               !                                  ! deep    outflow (from Red sea)
               tsa_tl(ji,jj,16,jp_tem) = tsa_tl(ji,jj,16,jp_tem)                                   &
                  &                       - hdiv_160_89_tl(16) * tsn(ji,jj,16,jp_tem)              &
                  &                       - hdiv_160_89(16) * tsn_tl(ji,jj,16,jp_tem)
               tsa_tl(ji,jj,16,jp_sal) = tsa_tl(ji,jj,16,jp_sal)                                   &
                  &                       - hdiv_160_89_tl(16) * tsn(ji,jj,16,jp_sal)              &
                  &                       - hdiv_160_89(16) * tsn_tl(ji,jj,16,jp_sal)
            END DO
         END DO
         !
         !                  ! ---------------- !
      CASE( 'spg' )         !  update (ua,va)  ! (call by dynspg module)
         !                  ! --------=======- !
         ! at this stage, (ua,va) are the after velocity, not the tendancy
         ! compute the velocity from the divergence at T-point
         !
         DO jj = mj0(88), mj1(88)              !** (160,88)   (Gulf of Aden side, north point)
            DO ji = mi0(160), mi1(160)                   ! 160, not 161 as it is a U-point)
               ua_tl(ji,jj,:) = - hdiv_161_88_kt_tl(:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                              * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         DO jj = mj0(87), mj1(87)              !** (160,87)   (Gulf of Aden side, south point)
            DO ji = mi0(160), mi1(160)                   ! 160, not 161 as it is a U-point)
               ua_tl(ji,jj,:) = - hdiv_161_87_tl(:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                           * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         DO jj = mj0(88), mj1(88)              !** profile of divergence at (160,89)   (Red sea side)
            DO ji = mi0(160), mi1(160)                   ! 88, not 89 as it is a V-point)
               va_tl(ji,jj,:) = - hdiv_160_89_kt_tl(:) / ( e1t(ji,jj+1) * e2t(ji,jj+1) * e3t(ji,jj+1,:) )   &
                  &                              * e1v(ji,jj) * e3v(ji,jj,:)
            END DO
         END DO
      END SELECT
      !
   END SUBROUTINE cla_bab_el_mandeb_tan

   SUBROUTINE cla_bab_el_mandeb_adj( cd_td )
      !!----------------------------------------------------------------------
      !!                ***  ROUTINE cla_bab_el_mandeb_adj  ***
      !!
      !! ** Purpose :   update the now horizontal divergence, the tracer tendancy
      !!              and the after velocity in vicinity of Bab el Mandeb ( Red Sea - Indian ocean).
      !!
      !! ** Method  :   compute the exchanges at each side of the strait :
      !!
      !!       surf. zio_flow
      !! (+ balance of emp) /\  |\\\\\\\\\\\|
      !!                    ||  |\\\\\\\\\\\|
      !!    deep zio_flow   ||  |\\\\\\\\\\\|
      !!            |  ||   ||  |\\\\\\\\\\\|
      !!        89  |  ||   ||  |\\\\\\\\\\\|
      !!            |__\/_v_||__|____________
      !!            !\\\\\\\\\\\|          surf. zio_flow
      !!            |\\\\\\\\\\\|<===    (+ balance of emp)
      !!            |\\\\\\\\\\\u
      !!        88  |\\\\\\\\\\\|<---      deep  zrecirc (upper+deep at 2 different levels)
      !!            |___________|__________
      !!            !\\\\\\\\\\\|
      !!            |\\\\\\\\\\\| ---\     deep  zrecirc (upper+deep)
      !!        87  !\\\\\\\\\\\u ===/   + deep  zio_flow   (all at the same level)
      !!            !\\\\\\\\\\\|
      !!            !___________|__________
      !!                160         161
      !!
      !!----------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in) ::   cd_td   ! ='div' update the divergence
      !                                         ! ='tra' update the tracers
      !                                         ! ='spg' update after velocity
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zemp_red_ad     ! temporary scalar
      REAL(wp) ::   zio_flow, zrecirc_upp, zrecirc_mid, zrecirc_bot
      !!---------------------------------------------------------------------
      !
      SELECT CASE( cd_td )
      !                     ! ---------------- !
      CASE( 'ini' )         !  initialisation  !
         !                  ! ---------------- !
         !
         zio_flow    = 0.4e6                       ! imposed in/out flow
         zrecirc_upp = 0.2e6                       ! imposed upper recirculation water
         zrecirc_bot = 0.5e6                       ! imposed bottom  recirculation water

         hdiv_161_88(:) = 0.e0                     ! (161,88) Gulf of Aden side, north point
         hdiv_161_87(:) = 0.e0                     ! (161,87) Gulf of Aden side, south point
         hdiv_160_89(:) = 0.e0                     ! (160,89) Red sea side
         hdiv_161_88_ad(:) = 0.e0                     ! (161,88) Gulf of Aden side, north point
         hdiv_161_87_ad(:) = 0.e0                     ! (161,87) Gulf of Aden side, south point
         hdiv_160_89_ad(:) = 0.e0                     ! (160,89) Red sea side

         DO jj = mj0(88), mj1(88)              !** profile of hdiv at (161,88)   (Gulf of Aden side, north point)
            DO ji = mi0(161), mi1(161)         !------------------------------
               DO jk = 1, 8                        ! surface in/out flow   (Ind -> Red)   (div >0)
                  hdiv_161_88(jk) = + zio_flow / ( 8. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                   ! recirculation water   (Ind -> Red)   (div >0)
               hdiv_161_88(20) =                 + zrecirc_upp   / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,20) )
               hdiv_161_88(21) = + ( zrecirc_bot - zrecirc_upp ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,21) )
            END DO
         END DO
         !
         DO jj = mj0(87), mj1(87)              !** profile of hdiv at (161,88)   (Gulf of Aden side, south point)
            DO ji = mi0(161), mi1(161)         !------------------------------
               !                                   ! deep out flow + recirculation   (Red -> Ind)   (div <0)
               hdiv_161_87(21) = - ( zio_flow + zrecirc_bot ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,21) )
            END DO
         END DO
         !
         DO jj = mj0(89), mj1(89)              !** profile of hdiv at (161,88)   (Red sea side)
            DO ji = mi0(160), mi1(160)         !------------------------------
               DO jk = 1, 8                        ! surface inflow    (Ind -> Red)   (div <0)
                  hdiv_160_89(jk) = - zio_flow / ( 8. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                   ! deep    outflow   (Red -> Ind)   (div >0)
               hdiv_160_89(16)    = + zio_flow / (      e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,16) )
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'div' )         !   update hdivn   ! (call by divcur module)
         !                  ! ---------=====-- !
         !                                     !** emp on the Red Sea   (div >0)
         zemp_red_ad = 0.e0                       !---------------------
         DO jj = mj1(89), mj0(89), -1                 !* profile of divergence at (160,89)   (Red sea side)
            DO ji = mi1(160), mi0(160), -1
               hdiv_160_89_kt_ad(:) = hdiv_160_89_kt_ad(:) + hdivn_ad(ji, jj,:)
               DO jk = 18, 1, -1                              ! increase the inflow from the Indian   (div <0)
                  zemp_red_ad =  zemp_red_ad - hdiv_160_89_ad(jk) / (10. * e1v(ji,jj) * e3v(ji,jj,jk) )
               END DO
               hdiv_160_89_ad(:) = hdiv_160_89_ad(:) + hdiv_160_89_kt_ad(:)
            END DO
         END DO
         DO jj = mj0(87), mj1(87)                  !* profile of divergence at (161,87)   (Gulf of Aden side, south point)
            DO ji = mi0(161), mi1(161)
               hdiv_161_87_ad(:) = hdiv_161_87_ad(:) + hdivn_ad(ji,jj,:)
            END DO
         END DO
         !                                     !** Correct hdivn (including emp adjustment)
         !                                     !-------------------------------------------
         DO jj = mj1(88), mj0(88), -1                  !* profile of hdiv at (161,88)   (Gulf of Aden side, north point)
            DO ji = mi1(161), mi0(161), -1
               hdiv_161_88_kt_ad(:) = hdiv_161_88_kt_ad(:) + hdivn_ad(ji,jj,:)
               DO jk = 8, 1, -1                              ! increase the inflow from the Indian   (div >0)
                  zemp_red_ad = zemp_red_ad + hdiv_161_88_ad(jk) / (8. * e2u(ji,jj) * e3u(ji,jj,jk) )
               END DO
               hdiv_161_88_ad(:) = hdiv_161_88_ad(:) + hdiv_161_88_kt_ad(:)
            END DO
         END DO
         zemp_red_ad = zemp_red_ad * 1.e-3               ! convert in m3
         IF( lk_mpp )   CALL mpp_sum( zemp_red_ad )   ! sum with other processors value

         DO jj = mj1(96), mj0(87), -1                  ! sum over the Red sea
            DO ji = mi1(160), mi0(148), -1
               emp_ad(ji,jj) = emp_ad(ji,jj) + zemp_red_ad * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
         !
         !                  ! ---------------- !
      CASE( 'tra' )         !  update (ta,sa)  ! (call by traadv module)
         !                  ! --------=======- !
         !===========================
         ! Direct model recomputation
         !===========================
         !
         DO jj = mj0(88), mj1(88)              !** (161,88)   (Gulf of Aden side, north point)
            DO ji = mi0(161), mi1(161)
               DO jk = 1, jpkm1                         ! surf inflow + reciculation (from Gulf of Aden)
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) - hdiv_161_88_kt(jk) * tsn(ji,jj,jk,jp_tem)
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) - hdiv_161_88_kt(jk) * tsn(ji,jj,jk,jp_sal)
               END DO
            END DO
         END DO
         DO jj = mj0(87), mj1(87)              !** (161,87)   (Gulf of Aden side, south point)
            DO ji = mi0(161), mi1(161)
               jk =  21                                 ! deep outflow + recirulation (combined flux)
               tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) + hdiv_161_88(20) * tsn(ji  ,jj+1,20,jp_tem)   &  ! upper recirculation from Gulf of Aden
                  &                        + hdiv_161_88(21) * tsn(ji  ,jj+1,21,jp_tem)   &  ! deep  recirculation from Gulf of Aden
                  &                        + hdiv_160_89(16) * tsn(ji-1,jj+2,16,jp_tem)      ! deep inflow from Red sea
               tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) + hdiv_161_88(20) * tsn(ji  ,jj+1,20,jp_sal)   &
                  &                        + hdiv_161_88(21) * tsn(ji  ,jj+1,21,jp_sal)   &
                  &                        + hdiv_160_89(16) * tsn(ji-1,jj+2,16,jp_sal)
            END DO
         END DO
         DO jj = mj0(89), mj1(89)              !** (161,88)   (Red sea side)
            DO ji = mi0(160), mi1(160)
               DO jk = 1, 14                            ! surface inflow (from Gulf of Aden)
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) - hdiv_160_89_kt(jk) * tsn(ji+1,jj-1,jk,jp_tem)
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) - hdiv_160_89_kt(jk) * tsn(ji+1,jj-1,jk,jp_sal)
               END DO
               !                                  ! deep    outflow (from Red sea)
               tsa(ji,jj,16,jp_tem) = tsa(ji,jj,16,jp_tem) - hdiv_160_89(16) * tsn(ji,jj,16,jp_tem)
               tsa(ji,jj,16,jp_sal) = tsa(ji,jj,16,jp_sal) - hdiv_160_89(16) * tsn(ji,jj,16,jp_sal)
            END DO
         END DO
         !=============
         ! Adjoint part
         !=============
         !
         DO jj = mj1(89), mj0(89), -1              !** (161,88)   (Red sea side)
            DO ji = mi1(160), mi0(160), -1
               hdiv_160_89_ad(16) = hdiv_160_89_ad(16) - tsa_ad(ji,jj,16,jp_sal) * tsn(ji,jj,16,jp_sal)
               tsn_ad(ji,jj,16,jp_sal) = tsn_ad(ji,jj,16,jp_sal) - tsa_ad(ji,jj,16,jp_sal) * tsn(ji,jj,16,jp_sal)
               hdiv_160_89_ad(16) = hdiv_160_89_ad(16) - tsa_ad(ji,jj,16,jp_tem) * tsn(ji,jj,16,jp_tem)
               tsn_ad(ji,jj,16,jp_tem) = tsn_ad(ji,jj,16,jp_tem) - tsa_ad(ji,jj,16,jp_tem) * tsn(ji,jj,16,jp_tem)
               DO jk = 14, 1, -1                            ! surface inflow (from Gulf of Aden)
                  hdiv_160_89_kt_ad(jk) = hdiv_160_89_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_sal) * tsn(ji+1,jj-1,jk,jp_sal)
                  tsn_ad(ji+1,jj-1,jk,jp_sal) = tsn_ad(ji+1,jj-1,jk,jp_sal) - tsa_ad(ji,jj,jk,jp_sal) * hdiv_160_89_kt(jk)
                  hdiv_160_89_kt_ad(jk) = hdiv_160_89_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_tem) * tsn(ji+1,jj-1,jk,jp_tem)
                  tsn_ad(ji+1,jj-1,jk,jp_tem) = tsn_ad(ji+1,jj-1,jk,jp_tem) - tsa_ad(ji,jj,jk,jp_tem) * hdiv_160_89_kt(jk)
               END DO
            END DO
         END DO
         DO jj = mj1(87), mj0(87), -1              !** (161,87)   (Gulf of Aden side, south point)
            DO ji = mi1(161), mi0(161), -1
               jk =  21                                 ! deep outflow + recirulation (combined flux)
               hdiv_161_88_ad(20) = hdiv_161_88_ad(20) + tsa_ad(ji,jj,jk,jp_sal) * tsn(ji  ,jj+1,20,jp_sal)
               hdiv_161_88_ad(21) = hdiv_161_88_ad(21) + tsa_ad(ji,jj,jk,jp_sal) * tsn(ji  ,jj+1,21,jp_sal)
               hdiv_160_89_ad(16) = hdiv_160_89_ad(16) + tsa_ad(ji,jj,jk,jp_sal) * tsn(ji-1,jj+2,16,jp_sal)
               tsn_ad(ji  ,jj+1,20,jp_sal) = tsn_ad(ji  ,jj+1,20,jp_sal) + tsa_ad(ji,jj,jk,jp_sal) * hdiv_161_88_ad(20)
               tsn_ad(ji  ,jj+1,21,jp_sal) = tsn_ad(ji  ,jj+1,21,jp_sal) + tsa_ad(ji,jj,jk,jp_sal) * hdiv_161_88_ad(21)
               tsn_ad(ji-1,jj+2,16,jp_sal) = tsn_ad(ji-1,jj+2,16,jp_sal) + tsa_ad(ji,jj,jk,jp_sal) * hdiv_160_89_ad(16)
               hdiv_161_88_ad(20) = hdiv_161_88_ad(20) + tsa_ad(ji,jj,jk,jp_tem) * tsn(ji  ,jj+1,20,jp_tem)
               hdiv_161_88_ad(21) = hdiv_161_88_ad(21) + tsa_ad(ji,jj,jk,jp_tem) * tsn(ji  ,jj+1,21,jp_tem)
               hdiv_160_89_ad(16) = hdiv_160_89_ad(16) + tsa_ad(ji,jj,jk,jp_tem) * tsn(ji-1,jj+2,16,jp_tem)
               tsn_ad(ji  ,jj+1,20,jp_tem) = tsn_ad(ji  ,jj+1,20,jp_tem) + tsa_ad(ji,jj,jk,jp_tem) * hdiv_161_88_ad(20)
               tsn_ad(ji  ,jj+1,21,jp_tem) = tsn_ad(ji  ,jj+1,21,jp_tem) + tsa_ad(ji,jj,jk,jp_tem) * hdiv_161_88_ad(21)
               tsn_ad(ji-1,jj+2,16,jp_tem) = tsn_ad(ji-1,jj+2,16,jp_tem) + tsa_ad(ji,jj,jk,jp_tem) * hdiv_160_89_ad(16)
            END DO
         END DO
         DO jj = mj1(88), mj0(88), -1              !** (161,88)   (Gulf of Aden side, north point)
            DO ji = mi1(161), mi0(161), -1
               DO jk = jpkm1, 1, -1                         ! surf inflow + reciculation (from Gulf of Aden)
                  hdiv_161_88_kt_ad(jk) = hdiv_161_88_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_sal) * tsn(ji,jj,jk,jp_sal)
                  tsn_ad(ji,jj,jk,jp_sal) = tsn_ad(ji,jj,jk,jp_sal) - hdiv_161_88_kt(jk) * tsa_ad(ji,jj,jk,jp_sal)
                  hdiv_161_88_kt_ad(jk) = hdiv_161_88_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_tem) * tsn(ji,jj,jk,jp_tem)
                  tsn_ad(ji,jj,jk,jp_tem) = tsn_ad(ji,jj,jk,jp_tem) - hdiv_161_88_kt(jk) * tsa_ad(ji,jj,jk,jp_tem)
               END DO
            END DO
         END DO
         !
         !                  ! ---------------- !
      CASE( 'spg' )         !  update (ua,va)  ! (call by dynspg module)
         !                  ! --------=======- !
         ! at this stage, (ua,va) are the after velocity, not the tendancy
         ! compute the velocity from the divergence at T-point
         !
         DO jj = mj1(88), mj0(88), -1              !** profile of divergence at (160,89)   (Red sea side)
            DO ji = mi1(160), mi0(160), -1                   ! 88, not 89 as it is a V-point)
               hdiv_160_89_kt_ad(:) =  hdiv_160_89_kt_ad(:) - va_ad(ji,jj,:) / ( e1t(ji,jj+1) * e2t(ji,jj+1) * e3t(ji,jj+1,:) )   &
                  &                              * e1v(ji,jj) * e3v(ji,jj,:)
               va_ad(ji,jj,:) = 0.0_wp
            END DO
         END DO
         DO jj = mj1(87), mj0(87), -1              !** (160,87)   (Gulf of Aden side, south point)
            DO ji = mi1(160), mi0(160), -1                   ! 160, not 161 as it is a U-point)
               hdiv_161_87_ad(:) =  hdiv_161_87_ad(:) - ua_ad(ji,jj,:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                           * e2u(ji,jj) * e3u(ji,jj,:)
               ua_ad(ji,jj,:) = 0.0_wp
            END DO
         END DO
         DO jj = mj1(88), mj0(88), -1              !** (160,88)   (Gulf of Aden side, north point)
            DO ji = mi1(160), mi0(160), -1                   ! 160, not 161 as it is a U-point)
               hdiv_161_88_kt_ad(:) = hdiv_161_88_kt_ad(:) - ua_ad(ji,jj,:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                              * e2u(ji,jj) * e3u(ji,jj,:)
               ua_ad(ji,jj,:) = 0.0_wp
            END DO
         END DO
      END SELECT
      !
   END SUBROUTINE cla_bab_el_mandeb_adj

   SUBROUTINE cla_gibraltar_tan( cd_td )
      !! -------------------------------------------------------------------
      !!                 ***  ROUTINE cla_gibraltar_tan  ***
      !!
      !! ** Purpose :   update the now horizontal divergence, the tracer
      !!              tendancyand the after velocity in vicinity of Gibraltar
      !!              strait ( Persian Gulf - Indian ocean ).
      !!
      !! ** Method :
      !!                     _______________________
      !!      deep  zio_flow /====|///////|====> surf. zio_flow
      !!    + deep  zrecirc  \----|///////|     (+balance of emp)
      !! 102                      u///////u
      !!      mid.  recicul    <--|///////|<==== deep  zio_flow
      !!                     _____|_______|_____
      !!      surf. zio_flow ====>|///////|
      !!    (+balance of emp)     |///////|
      !! 101                      u///////|
      !!      mid.  recicul    -->|///////|               Caution: zrecirc split into
      !!      deep  zrecirc  ---->|///////|                  upper & bottom recirculation
      !!                   _______|_______|_______
      !!                     139     140     141
      !!
      !!---------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in) ::   cd_td   ! ='div' update the divergence
      !                                         ! ='tra' update the tracers
      !                                         ! ='spg' update after velocity
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zemp_med     ! temporary scalar
      REAL(wp) ::   zio_flow, zrecirc_upp, zrecirc_mid, zrecirc_bot
      !!---------------------------------------------------------------------
      !
      SELECT CASE( cd_td )
      !                     ! ---------------- !
      CASE( 'ini' )         !  initialisation  !
         !                  ! ---------------- !
         !                                     !** initialization of the velocity
         hdiv_139_101(:) = 0.e0                     !  139,101 (Atlantic side, south point)
         hdiv_139_102(:) = 0.e0                     !  139,102 (Atlantic side, north point)
         hdiv_141_102(:) = 0.e0                     !  141,102 (Med sea  side)
         hdiv_139_101_tl(:) = 0.e0                     !  139,101 (Atlantic side, south point)
         hdiv_139_102_tl(:) = 0.e0                     !  139,102 (Atlantic side, north point)
         hdiv_141_102_tl(:) = 0.e0                     !  141,102 (Med sea  side)

         !                                     !** imposed transport
         zio_flow    = 0.8e6                        ! inflow surface  water
         zrecirc_mid = 0.7e6                        ! middle recirculation water
         zrecirc_upp = 2.5e6                        ! upper  recirculation water
         zrecirc_bot = 3.5e6                        ! bottom recirculation water
         !
         DO jj = mj0(101), mj1(101)            !** profile of hdiv at 139,101 (Atlantic side, south point)
            DO ji = mi0(139), mi1(139)         !-----------------------------
               DO jk = 1, 14                        ! surface in/out flow (Atl -> Med)   (div >0)
                  hdiv_139_101(jk) = + zio_flow / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               DO jk = 15, 20                       ! middle  reciculation (Atl 101 -> Atl 102)   (div >0)
                  hdiv_139_101(jk) = + zrecirc_mid / ( 6. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                    ! upper reciculation (Atl 101 -> Atl 101)   (div >0)
               hdiv_139_101(21) =               + zrecirc_upp / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               !
               !                                    ! upper & bottom reciculation (Atl 101 -> Atl 101 & 102)   (div >0)
               hdiv_139_101(22) = ( zrecirc_bot - zrecirc_upp ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO
         END DO
         DO jj = mj0(102), mj1(102)            !** profile of hdiv at 139,102 (Atlantic side, north point)
            DO ji = mi0(139), mi1(139)         !-----------------------------
               DO jk = 15, 20                       ! middle reciculation (Atl 101 -> Atl 102)   (div <0)
                  hdiv_139_102(jk) = - zrecirc_mid / ( 6. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                    ! outflow of Mediterranean sea + deep recirculation   (div <0)
               hdiv_139_102(22) = - ( zio_flow + zrecirc_bot ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO
         END DO
         DO jj = mj0(102), mj1(102)            !** velocity profile at 141,102  (Med sea side)
            DO ji = mi0(141), mi1(141)         !------------------------------
               DO  jk = 1, 14                       ! surface inflow in the Med     (div <0)
                  hdiv_141_102(jk) = - zio_flow / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                    ! deep    outflow toward the Atlantic    (div >0)
               hdiv_141_102(21)    = + zio_flow / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'div' )         !   update hdivn   ! (call by divcur module)
         !                  ! ---------=====-- !
         !                                     !** emp on the Mediterranean Sea  (div >0)
         zemp_med = 0.e0                       !-------------------------------
         DO jj = mj0(96), mj1(110)                  ! sum over the Med sea
            DO ji = mi0(141),mi1(181)
               zemp_med = zemp_med + emp_tl(ji,jj) * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
         DO jj = mj0(96), mj1(96)                   ! minus 2 points in Red Sea
            DO ji = mi0(148),mi1(148)
               zemp_med = zemp_med - emp_tl(ji,jj) * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
            DO ji = mi0(149),mi1(149)
               zemp_med = zemp_med - emp_tl(ji,jj) * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
         IF( lk_mpp )   CALL mpp_sum( zemp_med )    ! sum with other processors value
         zemp_med = zemp_med * 1.e-3                ! convert in m3
         !
         !                                     !** Correct hdivn (including emp adjustment)
         !                                     !-------------------------------------------
         DO jj = mj0(101), mj1(101)                 !* 139,101 (Atlantic side, south point)
            DO ji = mi0(139), mi1(139)
               hdiv_139_101_kt_tl(:) = hdiv_139_101_tl(:)
               DO jk = 1, 14                              ! increase the inflow from the Atlantic   (div >0)
                  hdiv_139_101_kt_tl(jk) = hdiv_139_101_tl(jk) + zemp_med / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               hdivn_tl(ji, jj,:) = hdivn_tl(ji, jj,:) + hdiv_139_101_kt_tl(:)
            END DO
         END DO
         DO jj = mj0(102), mj1(102)                 !* 139,102 (Atlantic side, north point)
            DO ji = mi0(139), mi1(139)
               hdivn_tl(ji,jj,:) = hdivn_tl(ji,jj,:) + hdiv_139_102_tl(:)
            END DO
         END DO
         DO jj = mj0(102), mj1(102)                 !* 141,102 (Med side)
            DO ji = mi0(141), mi1(141)
               hdiv_141_102_tl(:) = hdiv_141_102_tl(:)
               DO jk = 1, 14                              ! increase the inflow from the Atlantic   (div <0)
                  hdiv_141_102_kt_tl(jk) = hdiv_141_102_tl(jk) - zemp_med / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               hdivn_tl(ji, jj,:) = hdivn_tl(ji, jj,:) + hdiv_141_102_kt_tl(:)
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'tra' )         !  update (ta,sa)  ! (call by traadv module)
         !                  ! --------=======- !
         !
         DO jj = mj0(101), mj1(101)            !** 139,101 (Atlantic side, south point)   (div >0)
            DO ji = mi0(139), mi1(139)
               DO jk = 1, jpkm1                         ! surf inflow + mid. & bottom reciculation (from Atlantic)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)       &
                     &                       - hdiv_139_101_kt_tl(jk) * tsn(ji,jj,jk,jp_tem)  &
                     &                       - hdiv_139_101_kt(jk) * tsn_tl(ji,jj,jk,jp_tem)
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal)       &
                     &                       - hdiv_139_101_kt_tl(jk) * tsn(ji,jj,jk,jp_sal)  &
                     &                       - hdiv_139_101_kt(jk) * tsn_tl(ji,jj,jk,jp_sal)
               END DO
            END DO
         END DO
         !
         DO jj = mj0(102), mj1(102)            !** 139,102 (Atlantic side, north point)   (div <0)
            DO ji = mi0(139), mi1(139)
               DO jk = 15, 20                            ! middle  reciculation (Atl 101 -> Atl 102)   (div <0)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)      &
                     &                       - hdiv_139_102_tl(jk) * tsn(ji,jj-1,jk,jp_tem)  & ! middle Atlantic recirculation
                     &                       - hdiv_139_102(jk) * tsn_tl(ji,jj-1,jk,jp_tem)    ! middle Atlantic recirculation
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal)      &
                     &                       - hdiv_139_102_tl(jk) * tsn(ji,jj-1,jk,jp_sal)  &
                     &                       - hdiv_139_102(jk) * tsn_tl(ji,jj-1,jk,jp_sal)
               END DO
               !                                         ! upper & bottom Atl. reciculation (Atl 101 -> Atl 102) - (div <0)
               !                                         ! deep Med flow                    (Med 102 -> Atl 102) - (div <0)
               tsa_tl(ji,jj,22,jp_tem) = tsa_tl(ji,jj,22,jp_tem)                           &
                  &                       + hdiv_141_102_tl(21) * tsn(ji+2,jj,21,jp_tem)   &  ! deep Med flow
                  &                       + hdiv_141_102(21) * tsn_tl(ji+2,jj,21,jp_tem)   &  ! deep Med flow
                  &                       + hdiv_139_101_tl(21) * tsn(ji,jj-1,21,jp_tem)   &  ! upper  Atlantic recirculation
                  &                       + hdiv_139_101(21) * tsn_tl(ji,jj-1,21,jp_tem)   &  ! upper  Atlantic recirculation
                  &                       + hdiv_139_101_tl(22) * tsn(ji,jj-1,22,jp_tem)   &  ! bottom Atlantic recirculation
                  &                       + hdiv_139_101(22) * tsn_tl(ji,jj-1,22,jp_tem)      ! bottom Atlantic recirculation
               tsa_tl(ji,jj,22,jp_sal) = tsa_tl(ji,jj,22,jp_sal)                                 &
                  &                       + hdiv_141_102_tl(21) * tsn(ji+2,jj,21,jp_sal)   &
                  &                       + hdiv_141_102(21) * tsn_tl(ji+2,jj,21,jp_sal)   &
                  &                       + hdiv_139_101_tl(21) * tsn(ji,jj-1,21,jp_sal)   &
                  &                       + hdiv_139_101(21) * tsn_tl(ji,jj-1,21,jp_sal)   &
                  &                       + hdiv_139_101_tl(22) * tsn(ji,jj-1,22,jp_sal)   &
                  &                       + hdiv_139_101(22) * tsn_tl(ji,jj-1,22,jp_sal)
            END DO
         END DO
         DO jj = mj0(102), mj1(102)                 !* 141,102 (Med side)   (div <0)
            DO ji = mi0(141), mi1(141)
               DO jk = 1, 14                             ! surface flow from Atlantic to Med sea
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)                                  &
                     &                       - hdiv_141_102_kt_tl(jk) * tsn(ji-2,jj-1,jk,jp_tem)     &
                     &                       - hdiv_141_102_kt(jk) * tsn_tl(ji-2,jj-1,jk,jp_tem)
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal)                                        &
                     &                       - hdiv_141_102_kt_tl(jk) * tsn(ji-2,jj-1,jk,jp_sal)     &
                     &                       - hdiv_141_102_kt(jk) * tsn_tl(ji-2,jj-1,jk,jp_sal)
               END DO
               !                                         ! deeper flow from Med sea to Atlantic
               tsa_tl(ji,jj,21,jp_tem) = tsa_tl(ji,jj,21,jp_tem)                                           &
                  &                    - hdiv_141_102_tl(21) * tsn(ji,jj,21,jp_tem)                  &
                  &                    - hdiv_141_102(21) * tsn_tl(ji,jj,21,jp_tem)
               tsa_tl(ji,jj,21,jp_sal) = tsa_tl(ji,jj,21,jp_sal)                                           &
                  &                    - hdiv_141_102_tl(21) * tsn(ji,jj,21,jp_sal)                  &
                  &                    - hdiv_141_102(21) * tsn_tl(ji,jj,21,jp_sal)
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'spg' )         !  update (ua,va)  ! (call by dynspg module)
         !                  ! --------=======- !
         ! at this stage, (ua,va) are the after velocity, not the tendancy
         ! compute the velocity from the divergence at T-point
         !
         DO jj = mj0(101), mj1(101)            !** 139,101 (Atlantic side, south point)
            DO ji = mi0(139), mi1(139)                    ! div >0 => ua >0, same sign
               ua_tl(ji,jj,:) = hdiv_139_101_kt_tl(:) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,:) )   &
                  &                             * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         DO jj = mj0(102), mj1(102)            !** 139,102 (Atlantic side, north point)
            DO ji = mi0(139), mi1(139)                    ! div <0 => ua <0, same sign
               ua_tl(ji,jj,:) = hdiv_139_102_tl(:) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,:) )   &
                  &                          * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         DO jj = mj0(102), mj1(102)            !** 140,102 (Med side) (140 not 141 as it is a U-point)
            DO ji = mi0(140), mi1(140)                    ! div >0 => ua <0, opposite sign
               ua_tl(ji,jj,:) = - hdiv_141_102_tl(:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                            * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE cla_gibraltar_tan

   SUBROUTINE cla_gibraltar_adj( cd_td )
      !! -------------------------------------------------------------------
      !!                 ***  ROUTINE cla_gibraltar_adj  ***
      !!
      !! ** Purpose :   update the now horizontal divergence, the tracer
      !!              tendancyand the after velocity in vicinity of Gibraltar
      !!              strait ( Persian Gulf - Indian ocean ).
      !!
      !! ** Method :
      !!                     _______________________
      !!      deep  zio_flow /====|///////|====> surf. zio_flow
      !!    + deep  zrecirc  \----|///////|     (+balance of emp)
      !! 102                      u///////u
      !!      mid.  recicul    <--|///////|<==== deep  zio_flow
      !!                     _____|_______|_____
      !!      surf. zio_flow ====>|///////|
      !!    (+balance of emp)     |///////|
      !! 101                      u///////|
      !!      mid.  recicul    -->|///////|               Caution: zrecirc split into
      !!      deep  zrecirc  ---->|///////|                  upper & bottom recirculation
      !!                   _______|_______|_______
      !!                     139     140     141
      !!
      !!---------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in) ::   cd_td   ! ='div' update the divergence
      !                                         ! ='tra' update the tracers
      !                                         ! ='spg' update after velocity
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zemp_med     ! temporary scalar
      REAL(wp) ::   zio_flow, zrecirc_upp, zrecirc_mid, zrecirc_bot
      !!---------------------------------------------------------------------
      !
      SELECT CASE( cd_td )
      !                     ! ---------------- !
      CASE( 'ini' )         !  initialisation  !
         !                  ! ---------------- !
         !                                     !** initialization of the velocity
         hdiv_139_101(:) = 0.e0                     !  139,101 (Atlantic side, south point)
         hdiv_139_102(:) = 0.e0                     !  139,102 (Atlantic side, north point)
         hdiv_141_102(:) = 0.e0                     !  141,102 (Med sea  side)
         hdiv_139_101_ad(:) = 0.e0                     !  139,101 (Atlantic side, south point)
         hdiv_139_102_ad(:) = 0.e0                     !  139,102 (Atlantic side, north point)
         hdiv_141_102_ad(:) = 0.e0                     !  141,102 (Med sea  side)

         !                                     !** imposed transport
         zio_flow    = 0.8e6                        ! inflow surface  water
         zrecirc_mid = 0.7e6                        ! middle recirculation water
         zrecirc_upp = 2.5e6                        ! upper  recirculation water
         zrecirc_bot = 3.5e6                        ! bottom recirculation water
         !
         DO jj = mj0(101), mj1(101)            !** profile of hdiv at 139,101 (Atlantic side, south point)
            DO ji = mi0(139), mi1(139)         !-----------------------------
               DO jk = 1, 14                        ! surface in/out flow (Atl -> Med)   (div >0)
                  hdiv_139_101(jk) = + zio_flow / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               DO jk = 15, 20                       ! middle  reciculation (Atl 101 -> Atl 102)   (div >0)
                  hdiv_139_101(jk) = + zrecirc_mid / ( 6. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                    ! upper reciculation (Atl 101 -> Atl 101)   (div >0)
               hdiv_139_101(21) =               + zrecirc_upp / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               !
               !                                    ! upper & bottom reciculation (Atl 101 -> Atl 101 & 102)   (div >0)
               hdiv_139_101(22) = ( zrecirc_bot - zrecirc_upp ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO
         END DO
         DO jj = mj0(102), mj1(102)            !** profile of hdiv at 139,102 (Atlantic side, north point)
            DO ji = mi0(139), mi1(139)         !-----------------------------
               DO jk = 15, 20                       ! middle reciculation (Atl 101 -> Atl 102)   (div <0)
                  hdiv_139_102(jk) = - zrecirc_mid / ( 6. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                    ! outflow of Mediterranean sea + deep recirculation   (div <0)
               hdiv_139_102(22) = - ( zio_flow + zrecirc_bot ) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO
         END DO
         DO jj = mj0(102), mj1(102)            !** velocity profile at 141,102  (Med sea side)
            DO ji = mi0(141), mi1(141)         !------------------------------
               DO  jk = 1, 14                       ! surface inflow in the Med     (div <0)
                  hdiv_141_102(jk) = - zio_flow / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               !                                    ! deep    outflow toward the Atlantic    (div >0)
               hdiv_141_102(21)    = + zio_flow / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'div' )         !   update hdivn   ! (call by divcur module)
         !                  ! ---------=====-- !
         !                                     !** Correct hdivn (including emp adjustment)
         !                                     !-------------------------------------------
         DO jj = mj1(102), mj0(102), -1                 !* 141,102 (Med side)
            DO ji = mi1(141), mi0(141), -1
               hdiv_141_102_kt_ad(:) =  hdiv_141_102_kt_ad(:) + hdivn_ad(ji, jj,:)
               DO jk = 14, 1, -1                              ! increase the inflow from the Atlantic   (div <0)
                  zemp_med = zemp_med - hdiv_141_102_kt_ad(jk) / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                  hdiv_141_102_ad(jk) = hdiv_141_102_ad(jk) + hdiv_141_102_kt_ad(jk)
               END DO
               hdiv_141_102_ad(:) = hdiv_141_102_ad(:) + hdiv_141_102_kt_ad(:)
            END DO
         END DO
         DO jj = mj1(102), mj0(102), -1                 !* 139,102 (Atlantic side, north point)
            DO ji = mi1(139), mi0(139), -1
               hdiv_139_102_ad(:) = hdiv_139_102_ad(:) +  hdivn_ad(ji,jj,:)
            END DO
         END DO
         DO jj = mj1(101), mj0(101), -1                 !* 139,101 (Atlantic side, south point)
            DO ji = mi1(139), mi0(139), -1
               hdiv_139_101_kt_ad(:) = hdiv_139_101_kt_ad(:) + hdivn_ad(ji, jj,:)
               DO jk = 14, 1, -1                              ! increase the inflow from the Atlantic   (div >0)
                  zemp_med = zemp_med + hdiv_139_101_kt_ad(jk) / ( 14. * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
                  hdiv_139_101_ad(jk) = hdiv_139_101_ad(jk) + hdiv_139_101_kt_ad(jk)
               END DO
               hdiv_139_101_kt_ad(:) = hdiv_139_101_ad(:) + hdiv_139_101_kt_ad(:)
            END DO
         END DO
         !                                     !** emp on the Mediterranean Sea  (div >0)
         zemp_med = zemp_med * 1.e-3                ! convert in m3
         IF( lk_mpp )   CALL mpp_sum( zemp_med )    ! sum with other processors value
         DO jj = mj1(96), mj0(96), -1               ! minus 2 points in Red Sea
            DO ji = mi1(149),mi0(149), -1
               emp_ad(ji,jj) = emp_ad(ji,jj) - zemp_med * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
            DO ji = mi1(148),mi0(148), -1
               emp_ad(ji,jj) = emp_ad(ji,jj) - zemp_med * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
         DO jj = mj1(110), mj0(96), -1                  ! sum over the Med sea
            DO ji = mi1(181), mi0(141), -1
               emp_ad(ji,jj) = emp_ad(ji,jj) - zemp_med * e1t(ji,jj) * e2t(ji,jj) * tmask_i(ji,jj)
            END DO
         END DO
         zemp_med = 0.e0                       !-------------------------------
         !                  ! ---------------- !
      CASE( 'tra' )         !  update (ta,sa)  ! (call by traadv module)
         !                  ! --------=======- !
         !===========================
         ! Direct model recomputation
         !===========================
         !
         DO jj = mj0(101), mj1(101)            !** 139,101 (Atlantic side, south point)   (div >0)
            DO ji = mi0(139), mi1(139)
               DO jk = 1, jpkm1                         ! surf inflow + mid. & bottom reciculation (from Atlantic)
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) - hdiv_139_101_kt(jk) * tsn(ji,jj,jk,jp_tem)
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) - hdiv_139_101_kt(jk) * tsn(ji,jj,jk,jp_sal)
               END DO
            END DO
         END DO
         !
         DO jj = mj0(102), mj1(102)            !** 139,102 (Atlantic side, north point)   (div <0)
            DO ji = mi0(139), mi1(139)
               DO jk = 15, 20                            ! middle  reciculation (Atl 101 -> Atl 102)   (div <0)
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) - hdiv_139_102(jk) * tsn(ji,jj-1,jk,jp_tem)  ! middle Atlantic recirculation
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) - hdiv_139_102(jk) * tsn(ji,jj-1,jk,jp_sal)
               END DO
               !                                         ! upper & bottom Atl. reciculation (Atl 101 -> Atl 102) - (div <0)
               !                                         ! deep Med flow                    (Med 102 -> Atl 102) - (div <0)
               tsa(ji,jj,22,jp_tem) = tsa(ji,jj,22,jp_tem) + hdiv_141_102(21) * tsn(ji+2,jj,21,jp_tem)   &  ! deep Med flow
                  &                        + hdiv_139_101(21) * tsn(ji,jj-1,21,jp_tem)   &  ! upper  Atlantic recirculation
                  &                        + hdiv_139_101(22) * tsn(ji,jj-1,22,jp_tem)      ! bottom Atlantic recirculation
               tsa(ji,jj,22,jp_sal) = tsa(ji,jj,22,jp_sal) + hdiv_141_102(21) * tsn(ji+2,jj,21,jp_sal)   &
                  &                        + hdiv_139_101(21) * tsn(ji,jj-1,21,jp_sal)   &
                  &                        + hdiv_139_101(22) * tsn(ji,jj-1,22,jp_sal)
            END DO
         END DO
         DO jj = mj0(102), mj1(102)                 !* 141,102 (Med side)   (div <0)
            DO ji = mi0(141), mi1(141)
               DO jk = 1, 14                             ! surface flow from Atlantic to Med sea
                  tsa(ji,jj,jk,jp_tem) = tsa(ji,jj,jk,jp_tem) - hdiv_141_102_kt(jk) * tsn(ji-2,jj-1,jk,jp_tem)
                  tsa(ji,jj,jk,jp_sal) = tsa(ji,jj,jk,jp_sal) - hdiv_141_102_kt(jk) * tsn(ji-2,jj-1,jk,jp_sal)
               END DO
               !                                         ! deeper flow from Med sea to Atlantic
               tsa(ji,jj,21,jp_tem) = tsa(ji,jj,21,jp_tem) - hdiv_141_102(21) * tsn(ji,jj,21,jp_tem)
               tsa(ji,jj,21,jp_sal) = tsa(ji,jj,21,jp_sal) - hdiv_141_102(21) * tsn(ji,jj,21,jp_sal)
            END DO
         END DO
         !=============
         ! Adjoint part
         !=============
         !
         DO jj = mj1(102), mj0(102), -1                 !* 141,102 (Med side)   (div <0)
            DO ji = mi1(141), mi0(141), -1
               !                                         ! deeper flow from Med sea to Atlantic
               hdiv_141_102_ad(21) = hdiv_141_102_ad(21) - tsa_ad(ji,jj,21,jp_sal) * tsn(ji,jj,21,jp_sal)
               tsn_ad(ji,jj,21,jp_sal) = tsn_ad(ji,jj,21,jp_sal) - tsa_ad(ji,jj,21,jp_sal) * hdiv_141_102(21)
               hdiv_141_102_ad(21) = hdiv_141_102_ad(21) - tsa_ad(ji,jj,21,jp_tem) * tsn(ji,jj,21,jp_tem)
               tsn_ad(ji,jj,21,jp_tem) = tsn_ad(ji,jj,21,jp_tem) - tsa_ad(ji,jj,21,jp_tem) * hdiv_141_102(21)
               DO jk = 14, 1, -1                             ! surface flow from Atlantic to Med sea
                  hdiv_141_102_kt_ad(jk) = hdiv_141_102_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_sal) * tsn(ji-2,jj-1,jk,jp_sal)
                  tsn_ad(ji-2,jj-1,jk,jp_sal) = tsn_ad(ji-2,jj-1,jk,jp_sal) - tsa_ad(ji,jj,jk,jp_sal) *  hdiv_141_102_kt(jk)
                  hdiv_141_102_kt_ad(jk) = hdiv_141_102_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_tem) * tsn(ji-2,jj-1,jk,jp_tem)
                  tsn_ad(ji-2,jj-1,jk,jp_tem) = tsn_ad(ji-2,jj-1,jk,jp_tem) - tsa_ad(ji,jj,jk,jp_tem) *  hdiv_141_102_kt(jk)
               END DO
            END DO
         END DO
         !
         DO jj = mj1(102), mj0(102), -1            !** 139,102 (Atlantic side, north point)   (div <0)
            DO ji = mi1(139), mi0(139), -1
               !                                         ! upper & bottom Atl. reciculation (Atl 101 -> Atl 102) - (div <0)
               !                                         ! deep Med flow                    (Med 102 -> Atl 102) - (div <0)
               hdiv_141_102_ad(21) = hdiv_141_102_ad(21) + tsa_ad(ji,jj,22,jp_sal) * tsn(ji+2,jj,21,jp_sal)
               hdiv_139_101_ad(21) = hdiv_139_101_ad(21) + tsa_ad(ji,jj,22,jp_sal) * tsn(ji,jj-1,21,jp_sal)
               hdiv_139_101_ad(22) = hdiv_139_101_ad(22) + tsa_ad(ji,jj,22,jp_sal) * tsn(ji,jj-1,22,jp_sal)
               tsn_ad(ji+2,jj,21,jp_sal) = tsn_ad(ji+2,jj,21,jp_sal) + tsa_ad(ji,jj,22,jp_sal) * hdiv_141_102(21)
               tsn_ad(ji,jj-1,21,jp_sal) = tsn_ad(ji,jj-1,21,jp_sal) + tsa_ad(ji,jj,22,jp_sal) * hdiv_139_101(21)
               tsn_ad(ji,jj-1,22,jp_sal) = tsn_ad(ji,jj-1,22,jp_sal) + tsa_ad(ji,jj,22,jp_sal) * hdiv_139_101(22)
               hdiv_141_102_ad(21) = hdiv_141_102_ad(21) + tsa_ad(ji,jj,22,jp_tem) * tsn(ji+2,jj,21,jp_tem)
               hdiv_139_101_ad(21) = hdiv_139_101_ad(21) + tsa_ad(ji,jj,22,jp_tem) * tsn(ji,jj-1,21,jp_tem)
               hdiv_139_101_ad(22) = hdiv_139_101_ad(22) + tsa_ad(ji,jj,22,jp_tem) * tsn(ji,jj-1,22,jp_tem)
               tsn_ad(ji+2,jj,21,jp_tem) = tsn_ad(ji+2,jj,21,jp_tem) + tsa_ad(ji,jj,22,jp_tem) * hdiv_141_102(21)
               tsn_ad(ji,jj-1,21,jp_tem) = tsn_ad(ji,jj-1,21,jp_tem) + tsa_ad(ji,jj,22,jp_tem) * hdiv_139_101(21)
               tsn_ad(ji,jj-1,22,jp_tem) = tsn_ad(ji,jj-1,22,jp_tem) + tsa_ad(ji,jj,22,jp_tem) * hdiv_139_101(22)
               DO jk = 20, 15, -1                            ! middle  reciculation (Atl 101 -> Atl 102)   (div <0)
                  hdiv_139_102_ad(jk) = hdiv_139_102_ad(jk) - tsa_ad(ji,jj,jk,jp_sal) * tsn(ji,jj-1,jk,jp_sal)
                  tsn_ad(ji,jj-1,jk,jp_sal) = tsn_ad(ji,jj-1,jk,jp_sal) - tsa_ad(ji,jj,jk,jp_sal) * hdiv_139_102(jk)
                  hdiv_139_102_ad(jk) = hdiv_139_102_ad(jk) - tsa_ad(ji,jj,jk,jp_tem) * tsn(ji,jj-1,jk,jp_tem)
                  tsn_ad(ji,jj-1,jk,jp_tem) = tsn_ad(ji,jj-1,jk,jp_tem) - tsa_ad(ji,jj,jk,jp_tem) * hdiv_139_102(jk)
               END DO
            END DO
         END DO
         DO jj = mj1(101), mj0(101), -1            !** 139,101 (Atlantic side, south point)   (div >0)
            DO ji = mi1(139), mi0(139), -1
               DO jk = jpkm1, 1, -1                         ! surf inflow + mid. & bottom reciculation (from Atlantic)
                  hdiv_139_101_kt_ad(jk) = hdiv_139_101_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_sal) * tsn(ji,jj,jk,jp_sal)
                  tsn_ad(ji,jj,jk,jp_sal) = tsn_ad(ji,jj,jk,jp_sal) -  tsa_ad(ji,jj,jk,jp_sal) * hdiv_139_101_kt(jk)
                  hdiv_139_101_kt_ad(jk) = hdiv_139_101_kt_ad(jk) - tsa_ad(ji,jj,jk,jp_tem) * tsn(ji,jj,jk,jp_tem)
                  tsn_ad(ji,jj,jk,jp_tem) = tsn_ad(ji,jj,jk,jp_tem) -  tsa_ad(ji,jj,jk,jp_tem) * hdiv_139_101_kt(jk)
               END DO
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'spg' )         !  update (ua,va)  ! (call by dynspg module)
         !                  ! --------=======- !
         ! at this stage, (ua,va) are the after velocity, not the tendancy
         ! compute the velocity from the divergence at T-point
         !
         DO jj = mj0(102), mj1(102)            !** 140,102 (Med side) (140 not 141 as it is a U-point)
            DO ji = mi0(140), mi1(140)                    ! div >0 => ua <0, opposite sign
               hdiv_141_102_ad(:) = hdiv_141_102_ad(:) - ua_ad(ji,jj,:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                            * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         DO jj = mj0(102), mj1(102)            !** 139,102 (Atlantic side, north point)
            DO ji = mi0(139), mi1(139)                    ! div <0 => ua <0, same sign
               hdiv_139_102_ad(:) = hdiv_139_102_ad(:) +  ua_ad(ji,jj,:) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,:) )   &
                  &                          * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         DO jj = mj0(101), mj1(101)            !** 139,101 (Atlantic side, south point)
            DO ji = mi0(139), mi1(139)                    ! div >0 => ua >0, same sign
               hdiv_139_101_kt_ad(:) = hdiv_139_101_kt_ad(:) +  ua_ad(ji,jj,:) / ( e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,:) )   &
                  &                             * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE cla_gibraltar_adj

   SUBROUTINE cla_hormuz_tan( cd_td )
      !! -------------------------------------------------------------------
      !!                   ***  ROUTINE div_hormuz_tan  ***
      !!
      !! ** Purpose :   update the now horizontal divergence, the tracer
      !!              tendancyand the after velocity in vicinity of Hormuz
      !!              strait ( Persian Gulf - Indian ocean ).
      !!
      !! ** Method  :   Hormuz strait
      !!            ______________
      !!            |/////|<==      surface inflow
      !!        94  |/////|
      !!            |/////|==>      deep    outflow
      !!            |_____|_______
      !!              171    172
      !!---------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in) ::   cd_td   ! ='ini' initialisation
      !!                                        ! ='div' update the divergence
      !!                                        ! ='tra' update the tracers
      !!                                        ! ='spg' update after velocity
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zio_flow     ! temporary scalar
      !!---------------------------------------------------------------------
      !
      SELECT CASE( cd_td )
      !                     ! ---------------- !
      CASE( 'ini' )         !  initialisation  !
         !                  ! ---------------- !
         !                                     !** profile of horizontal divergence due to cross-land advection
         zio_flow  = 1.e6                          ! imposed in/out flow
         !
         hdiv_172_94(:) = 0.e0
         hdiv_172_94_tl(:) = 0.e0
         !
         DO jj = mj0(94), mj1(94)                  ! in/out flow at (i,j) = (172,94)
            DO ji = mi0(172), mi1(172)
               DO jk = 1, 8                            ! surface inflow  (Indian ocean to Persian Gulf) (div<0)
                  hdiv_172_94(jk) = - ( zio_flow / 8.e0 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               DO jk = 16, 18                          ! deep    outflow (Persian Gulf to Indian ocean) (div>0)
                  hdiv_172_94(jk) = + ( zio_flow / 3.e0 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
            END DO
         END DO
         !                                     !** T & S profile in the Hormuz strait (use in deep outflow)
         !      Temperature       and         Salinity
         t_171_94_hor(:)  = 0.e0   ;   s_171_94_hor(:)  = 0.e0
         t_171_94_hor(16) = 18.4   ;   s_171_94_hor(16) = 36.27
         t_171_94_hor(17) = 17.8   ;   s_171_94_hor(17) = 36.4
         t_171_94_hor(18) = 16.    ;   s_171_94_hor(18) = 36.27
         !
         !                  ! ---------------- !
      CASE( 'div' )         !   update hdivn   ! (call by divcur module)
         !                  ! ---------=====-- !
         !
         DO jj = mj0(94), mj1(94)              !** 172,94 (Indian ocean side)
            DO ji = mi0(172), mi1(172)
               hdivn_tl(ji,jj,:) = hdivn_tl(ji,jj,:) + hdiv_172_94_tl(:)
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'tra' )         !  update (ta,sa)  ! (call by traadv module)
         !                  ! --------=======- !
         !
         DO jj = mj0(94), mj1(94)              !** 172,94 (Indian ocean side)
            DO ji = mi0(172), mi1(172)
               DO jk = 1, 8                          ! surface inflow   (Indian ocean to Persian Gulf) (div<0)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem)                        &
                     &                       - hdiv_172_94_tl(jk) * tsn(ji,jj,jk,jp_tem)   &
                     &                       - hdiv_172_94(jk) * tsn_tl(ji,jj,jk,jp_tem)
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal)                        &
                     &                       - hdiv_172_94_tl(jk) * tsn(ji,jj,jk,jp_sal)   &
                     &                       - hdiv_172_94(jk) * tsn_tl(ji,jj,jk,jp_sal)
               END DO
               DO jk = 16, 18                        ! deep outflow     (Persian Gulf to Indian ocean) (div>0)
                  tsa_tl(ji,jj,jk,jp_tem) = tsa_tl(ji,jj,jk,jp_tem) - hdiv_172_94_tl(jk) * t_171_94_hor(jk)
                  tsa_tl(ji,jj,jk,jp_sal) = tsa_tl(ji,jj,jk,jp_sal) - hdiv_172_94_tl(jk) * s_171_94_hor(jk)
               END DO
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'spg' )         !  update (ua,va)  ! (call by dynspg module)
         !                  ! --------=======- !
         ! No barotropic flow through Hormuz strait
         ! at this stage, (ua,va) are the after velocity, not the tendancy
         ! compute the velocity from the divergence at T-point
         DO jj = mj0(94), mj1(94)              !** 171,94 (Indian ocean side) (171 not 172 as it is the western U-point)
            DO ji = mi0(171), mi1(171)                ! div >0 => ua >0, opposite sign
               ua_tl(ji,jj,:) = - hdiv_172_94_tl(:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                           * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE cla_hormuz_tan


   SUBROUTINE cla_hormuz_adj( cd_td )
      !! -------------------------------------------------------------------
      !!                   ***  ROUTINE div_hormuz_adj  ***
      !!
      !! ** Purpose :   update the now horizontal divergence, the tracer
      !!              tendancyand the after velocity in vicinity of Hormuz
      !!              strait ( Persian Gulf - Indian ocean ).
      !!
      !! ** Method  :   Hormuz strait
      !!            ______________
      !!            |/////|<==      surface inflow
      !!        94  |/////|
      !!            |/////|==>      deep    outflow
      !!            |_____|_______
      !!              171    172
      !!---------------------------------------------------------------------
      CHARACTER(len=1), INTENT(in) ::   cd_td   ! ='ini' initialisation
      !!                                        ! ='div' update the divergence
      !!                                        ! ='tra' update the tracers
      !!                                        ! ='spg' update after velocity
      !!
      INTEGER  ::   ji, jj, jk   ! dummy loop indices
      REAL(wp) ::   zio_flow     ! temporary scalar
      !!---------------------------------------------------------------------
      !
      SELECT CASE( cd_td )
      !                     ! ---------------- !
      CASE( 'ini' )         !  initialisation  !
         !                  ! ---------------- !
         !                                     !** profile of horizontal divergence due to cross-land advection
         zio_flow  = 1.e6                          ! imposed in/out flow
         !
         hdiv_172_94(:) = 0.e0
         hdiv_172_94_ad(:) = 0.e0
         !
         DO jj = mj0(94), mj1(94)                  ! in/out flow at (i,j) = (172,94)
            DO ji = mi0(172), mi1(172)
               DO jk = 1, 8                            ! surface inflow  (Indian ocean to Persian Gulf) (div<0)
                  hdiv_172_94(jk) = - ( zio_flow / 8.e0 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
               DO jk = 16, 18                          ! deep    outflow (Persian Gulf to Indian ocean) (div>0)
                  hdiv_172_94(jk) = + ( zio_flow / 3.e0 * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) )
               END DO
            END DO
         END DO
         !                                     !** T & S profile in the Hormuz strait (use in deep outflow)
         !      Temperature       and         Salinity
         t_171_94_hor(:)  = 0.e0   ;   s_171_94_hor(:)  = 0.e0
         t_171_94_hor(16) = 18.4   ;   s_171_94_hor(16) = 36.27
         t_171_94_hor(17) = 17.8   ;   s_171_94_hor(17) = 36.4
         t_171_94_hor(18) = 16.    ;   s_171_94_hor(18) = 36.27
         !
         !                  ! ---------------- !
      CASE( 'div' )         !   update hdivn   ! (call by divcur module)
         !                  ! ---------=====-- !
         !
         DO jj = mj1(94), mj0(94), -1              !** 172,94 (Indian ocean side)
            DO ji = mi1(172), mi0(172), -1
               hdiv_172_94_ad(:) = hdiv_172_94_ad(:) + hdivn_ad(ji,jj,:)
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'tra' )         !  update (ta,sa)  ! (call by traadv module)
         !                  ! --------=======- !
         !
         DO jj = mj1(94), mj0(94), -1              !** 172,94 (Indian ocean side)
            DO ji = mi1(172), mi0(172), -1
               DO jk = 18, 16, -1                        ! deep outflow     (Persian Gulf to Indian ocean) (div>0)
                  hdiv_172_94_ad(jk) = hdiv_172_94_ad(jk) - tsa_ad(ji,jj,jk,jp_sal) * s_171_94_hor(jk)
                  hdiv_172_94_ad(jk) = hdiv_172_94_ad(jk) - tsa_ad(ji,jj,jk,jp_tem) * s_171_94_hor(jk)
               END DO
               DO jk = 8, 1, -1                          ! surface inflow   (Indian ocean to Persian Gulf) (div<0)
                  hdiv_172_94_ad(jk) = hdiv_172_94_ad(jk) - tsa_ad(ji,jj,jk,jp_sal) * tsn(ji,jj,jk,jp_sal)
                  tsn_ad(ji,jj,jk,jp_sal) = tsn_ad(ji,jj,jk,jp_sal) - tsa_ad(ji,jj,jk,jp_sal) * hdiv_172_94(jk)
                  hdiv_172_94_ad(jk) = hdiv_172_94_ad(jk) - tsa_ad(ji,jj,jk,jp_tem) * tsn(ji,jj,jk,jp_tem)
                  tsn_ad(ji,jj,jk,jp_tem) = tsn_ad(ji,jj,jk,jp_tem) - tsa_ad(ji,jj,jk,jp_tem) * hdiv_172_94(jk)
               END DO
            END DO
         END DO
         !                  ! ---------------- !
      CASE( 'spg' )         !  update (ua,va)  ! (call by dynspg module)
         !                  ! --------=======- !
         ! No barotropic flow through Hormuz strait
         ! at this stage, (ua,va) are the after velocity, not the tendancy
         ! compute the velocity from the divergence at T-point
         DO jj = mj0(94), mj1(94)              !** 171,94 (Indian ocean side) (171 not 172 as it is the western U-point)
            DO ji = mi0(171), mi1(171)                ! div >0 => ua >0, opposite sign
               hdiv_172_94_ad(:) = hdiv_172_94_ad(:) - ua_ad(ji,jj,:) / ( e1t(ji+1,jj) * e2t(ji+1,jj) * e3t(ji+1,jj,:) )   &
                  &                           * e2u(ji,jj) * e3u(ji,jj,:)
            END DO
         END DO
         !
      END SELECT
      !
   END SUBROUTINE cla_hormuz_adj

   SUBROUTINE cla_div_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE cla_divadj_tst ***
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
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & jt

      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zemp_tlin,       & ! Tangent input
         & zemp_tlout,      & ! Tangent output
         & zemp_adin,       & ! adjoint input
         & zemp_adout,      & ! adjoint output
         & zemp

      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zhdivn_tlin,     & ! Tangent input
         & zhdivn_tlout,    & ! Tangent output
         & zhdivn_adin,     & ! adjoint input
         & zhdivn_adout,    & ! adjoint output
         & zhdivn

      REAL(KIND=wp) ::   &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp1_1,       & !   scalar product components
         & zsp1_2,       &
         & zsp2,         & ! scalar product involving the adjoint routine
         & zsp2_1,       & !   scalar product components
         & zsp2_2,       &
         & zsp2_3,       &
         & zsp2_4,       &
         & zsp2_5

      CHARACTER(LEN=14) :: cl_name

      ALLOCATE( &
         & zhdivn_tlin(jpi,jpj,jpk),   &
         & zhdivn_tlout(jpi,jpj,jpk),  &
         & zhdivn_adin(jpi,jpj,jpk),   &
         & zhdivn_adout(jpi,jpj,jpk),  &
         & zemp_tlin(jpi,jpj),     &
         & zemp_tlout(jpi,jpj),    &
         & zemp_adin(jpi,jpj),     &
         & zemp_adout(jpi,jpj),    &
         & zhdivn(jpi,jpj,jpk),        &
         & zemp(jpi,jpj) )

      DO jt = 1, 3
         !==================================================================
         ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
         !    dy = ( hdivb_tl, hdivn_tl )
         !==================================================================

         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------

          zhdivn_tlin(:,:,:)  = 0._wp
          zhdivn_tlout(:,:,:) = 0._wp
          zhdivn_adin(:,:,:)  = 0._wp
          zhdivn_adout(:,:,:) = 0._wp
          zemp_tlin(:,:)      = 0._wp
          zemp_tlout(:,:)     = 0._wp
          zemp_adin(:,:)      = 0._wp
          zemp_adout(:,:)     = 0._wp
          zhdivn(:,:,:)       = 0._wp
          zemp(:,:)           = 0._wp

         hdivn_tl(:,:,:) = 0._wp
         emp_tl(:,:)     = 0._wp
         hdivn_ad(:,:,:) = 0._wp
         emp_ad(:,:)     = 0._wp

         CALL grid_random( zemp,   'T', 0.0_wp, stdssh )
         CALL grid_random( zhdivn, 'T', 0.0_wp, stdu )

         DO jj = nldj, nlej
            DO ji = nldi, nlei
               zemp_tlin(ji,jj) = zemp(ji,jj)
            END DO
         END DO
         DO jk = 1, jpk
            DO jj = nldj, nlej
               DO ji = nldi, nlei
                  zhdivn_tlin(ji,jj,jk) = zhdivn(ji,jj,jk)
               END DO
            END DO
         END DO
         !--------------------------------------------------------------------
         ! Call the tangent routine: dy = L dx
         !--------------------------------------------------------------------
         emp_tl(:,:)     = zemp_tlin(:,:)
         hdivn_tl(:,:,:) = zhdivn_tlin(:,:,:)

         SELECT CASE (jt)
         CASE(1)
            nbab = 1
            ngib = 0
            nhor = 0
            CALL cla_div_tan( nit000 )
         CASE(2)
            nbab = 0
            ngib = 1
            nhor = 0
            CALL cla_div_tan( nit000 )
         CASE(3)
            nbab = 0
            ngib = 0
            nhor = 1
            CALL cla_div_tan( nit000 )
         END SELECT

         zhdivn_tlout(:,:,:) = hdivn_tl(:,:,:)


         !--------------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !--------------------------------------------------------------------
         DO jk = 1, jpk
           DO jj = nldj, nlej
              DO ji = nldi, nlei
                 zhdivn_adin(ji,jj,jk) = zhdivn_tlout(ji,jj,jk) &
                    &               * e1t(ji,jj) * e2t(ji,jj) * e3t(ji,jj,jk) &
                    &               * tmask(ji,jj,jk)
               END DO
            END DO
         END DO

         !--------------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !--------------------------------------------------------------------

         zsp1   = DOT_PRODUCT( zhdivn_tlout, zhdivn_adin )

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------

         hdivn_ad(:,:,:) = zhdivn_adin(:,:,:)

         SELECT CASE (jt)
         CASE(1)
            nbab = 1
            ngib = 0
            nhor = 0
            CALL cla_div_adj( nit000 )
         CASE(2)
            nbab = 0
            ngib = 1
            nhor = 0
            CALL cla_div_adj( nit000 )
         CASE(3)
            nbab = 0
            ngib = 0
            nhor = 1
            CALL cla_div_adj( nit000 )
         END SELECT

         zemp_adout   (:,:) = emp_ad   (:,:)
         zhdivn_adout(:,:,:) = hdivn_ad(:,:,:)

         !--------------------------------------------------------------------
         ! Compute the scalar product: dx^T L^T W dy
         !--------------------------------------------------------------------

         zsp2_1 = DOT_PRODUCT( zhdivn_tlin,    zhdivn_adout    )
         zsp2_2 = DOT_PRODUCT( zemp_tlin,    zemp_adout    )
         zsp2   = zsp2_1 + zsp2_2

         SELECT CASE (jt)
         CASE(1)
            cl_name = 'cladivadj babm'
         CASE(2)
            cl_name = 'cladivadj gibr'
         CASE(3)
            cl_name = 'cladivadj horm'
         END SELECT
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )
      END DO

      DEALLOCATE( &
         & zhdivn_tlin,   &
         & zhdivn_tlout,  &
         & zhdivn_adin,   &
         & zhdivn_adout,  &
         & zemp_tlin,     &
         & zemp_tlout,    &
         & zemp_adin,     &
         & zemp_adout,    &
         & zhdivn,        &
         & zemp )

   END SUBROUTINE cla_div_adj_tst

   SUBROUTINE cla_traadv_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE cla_divadj_tst ***
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
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & jt

      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & ztn_tlin ,     & ! Tangent input
         & zsn_tlin ,     zta_tlin ,     zsa_tlin ,     & ! Tangent input
         & ztn_adout,     & ! Adjoint output
         & zsn_adout,     zta_adout,     zsa_adout,     & ! Adjoint output
         & zta_tlout,     zsa_tlout,     & ! Tangent output
         & zta_adin ,     zsa_adin ,     & ! Adjoint input
         & zr             ! 3D random field
      CHARACTER(LEN=14) ::&
         & cl_name
      ! Allocate memory

      ALLOCATE( &
         & ztn_tlin( jpi,jpj,jpk),     zsn_tlin( jpi,jpj,jpk),     zta_tlin( jpi,jpj,jpk),     &
         & zsa_tlin( jpi,jpj,jpk),     zta_tlout(jpi,jpj,jpk),     zsa_tlout(jpi,jpj,jpk),     &
         & zta_adin( jpi,jpj,jpk),     zsa_adin( jpi,jpj,jpk),                                 &
         & ztn_adout(jpi,jpj,jpk),                                                             &
         & zsn_adout(jpi,jpj,jpk),     zta_adout(jpi,jpj,jpk),     zsa_adout(jpi,jpj,jpk),     &
         & zr(       jpi,jpj,jpk)      &
         & )

      DO jt = 1, 3
         !==================================================================
         ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
         !    dy = ( hdivb_tl, hdivn_tl )
         !==================================================================

         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------
         ztn_tlin( :,:,:) = 0.0_wp
         zsn_tlin( :,:,:) = 0.0_wp
         zta_tlin( :,:,:) = 0.0_wp
         zsa_tlin( :,:,:) = 0.0_wp
         zta_tlout(:,:,:) = 0.0_wp
         zsa_tlout(:,:,:) = 0.0_wp
         zta_adin( :,:,:) = 0.0_wp
         zsa_adin( :,:,:) = 0.0_wp
         ztn_adout(:,:,:) = 0.0_wp
         zsn_adout(:,:,:) = 0.0_wp
         zta_adout(:,:,:) = 0.0_wp
         zsa_adout(:,:,:) = 0.0_wp
         zr(       :,:,:) = 0.0_wp

         tsn_ad(:,:,:,jp_tem) = 0.0_wp
         tsn_ad(:,:,:,jp_sal) = 0.0_wp

         CALL grid_random(  zr, 'T', 0.0_wp, stdt )
         DO jk = 1, jpk
           DO jj = nldj, nlej
              DO ji = nldi, nlei
                 ztn_tlin(ji,jj,jk) = zr(ji,jj,jk)
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
                 zta_tlin(ji,jj,jk) = zr(ji,jj,jk)
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

         tsn_tl(:,:,:,jp_tem) = ztn_tlin(:,:,:)
         tsn_tl(:,:,:,jp_sal) = zsn_tlin(:,:,:)
         tsa_tl(:,:,:,jp_tem) = zta_tlin(:,:,:)
         tsa_tl(:,:,:,jp_sal) = zsa_tlin(:,:,:)

         !--------------------------------------------------------------------
         ! Call the tangent routine: dy = L dx
         !--------------------------------------------------------------------

         SELECT CASE (jt)
         CASE(1)
            nbab = 1
            ngib = 0
            nhor = 0
            CALL cla_traadv_tan( nit000 )
         CASE(2)
            nbab = 0
            ngib = 1
            nhor = 0
            CALL cla_traadv_tan( nit000 )
         CASE(3)
            nbab = 0
            ngib = 0
            nhor = 1
            CALL cla_traadv_tan( nit000 )
         END SELECT

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

         SELECT CASE (jt)
         CASE(1)
            nbab = 1
            ngib = 0
            nhor = 0
            CALL cla_traadv_adj( nit000 )
         CASE(2)
            nbab = 0
            ngib = 1
            nhor = 0
            CALL cla_traadv_adj( nit000 )
         CASE(3)
            nbab = 0
            ngib = 0
            nhor = 1
            CALL cla_traadv_adj( nit000 )
         END SELECT

         ztn_adout(:,:,:) = tsn_ad(:,:,:,jp_tem)
         zsn_adout(:,:,:) = tsn_ad(:,:,:,jp_sal)
         zta_adout(:,:,:) = tsa_ad(:,:,:,jp_tem)
         zsa_adout(:,:,:) = tsa_ad(:,:,:,jp_sal)

         zsp2 = DOT_PRODUCT( ztn_tlin, ztn_adout ) &
            & + DOT_PRODUCT( zsn_tlin, zsn_adout ) &
            & + DOT_PRODUCT( zta_tlin, zta_adout ) &
            & + DOT_PRODUCT( zsa_tlin, zsa_adout )

         SELECT CASE (jt)
         CASE(1)
            cl_name = 'clatraadv babm'
         CASE(2)
            cl_name = 'clatraadv gibr'
         CASE(3)
            cl_name = 'clatraadv horm'
         END SELECT
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )
      END DO

      DEALLOCATE(         &
         & ztn_tlin ,     zsn_tlin ,     &
         & zta_tlin ,     zsa_tlin ,     zta_tlout,     zsa_tlout,     zta_adin ,     &
         & zsa_adin ,     ztn_adout,     &
         & zsn_adout,     zta_adout,     zsa_adout,     zr                            &
         & )

   END SUBROUTINE cla_traadv_adj_tst

   SUBROUTINE cla_dynspg_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE cla_divadj_tst ***
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
         & ji,    &        ! dummy loop indices
         & jj,    &
         & jk,    &
         & jt

      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zua_tlin ,     & ! Tangent input
         & zva_tlin ,     & ! Tangent input
         & zua_tlout,     & ! Tangent output
         & zva_tlout,     & ! Tangent output
         & zua_adin ,     & ! Adjoint input
         & zva_adin ,     & ! Adjoint input
         & zua_adout,     & ! Adjoint output
         & zva_adout,     & ! Adjoint output
         & zr3d           ! 3D random field
       CHARACTER(LEN=14) :: &
         & cl_name
      ! Allocate memory

      ALLOCATE( &
         & zua_tlin(   jpi,jpj,jpk),     &
         & zva_tlin(   jpi,jpj,jpk),     &
         & zua_tlout(  jpi,jpj,jpk),     &
         & zva_tlout(  jpi,jpj,jpk),     &
         & zua_adin(   jpi,jpj,jpk),     &
         & zva_adin(   jpi,jpj,jpk),     &
         & zua_adout(  jpi,jpj,jpk),     &
         & zva_adout(  jpi,jpj,jpk),     &
         & zr3d(       jpi,jpj,jpk)      )

      DO jt = 1, 3
         !==================================================================
         ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
         !    dy = ( hdivb_tl, hdivn_tl )
         !==================================================================

         !--------------------------------------------------------------------
         ! Reset the tangent and adjoint variables
         !--------------------------------------------------------------------
         ua_ad( :,:,:) = 0.0_wp
         va_ad( :,:,:) = 0.0_wp
         !--------------------------------------------------------------------
         ! Initialize the tangent input with random noise: dx
         !--------------------------------------------------------------------

         CALL grid_random( zr3d, 'U', 0.0_wp, stdu )
         zua_tlin(:,:,:) = zr3d(:,:,:)
         CALL grid_random( zr3d, 'V', 0.0_wp, stdv )
         zva_tlin(:,:,:) = zr3d(:,:,:)

         ua_tl = zua_tlin
         va_tl = zva_tlin
         !--------------------------------------------------------------------
         ! Call the tangent routine: dy = L dx
         !--------------------------------------------------------------------

         SELECT CASE (jt)
         CASE(1)
            nbab = 1
            ngib = 0
            nhor = 0
            CALL cla_traadv_tan( nit000 )
         CASE(2)
            nbab = 0
            ngib = 1
            nhor = 0
            CALL cla_traadv_tan( nit000 )
         CASE(3)
            nbab = 0
            ngib = 0
            nhor = 1
            CALL cla_traadv_tan( nit000 )
         END SELECT

         zua_tlout = ua_tl
         zva_tlout = va_tl

         !--------------------------------------------------------------------
         ! Initialize the adjoint variables: dy^* = W dy
         !--------------------------------------------------------------------
         DO jk = 1, jpk
           DO jj = nldj, nlej
              DO ji = nldi, nlei
                 zua_adin(ji,jj,jk) = zua_tlout(ji,jj,jk) &
                    &               * e1u(ji,jj) * e2u(ji,jj) * e3u(ji,jj,jk) &
                    &               * umask(ji,jj,jk)
                 zva_adin(ji,jj,jk) = zva_tlout(ji,jj,jk) &
                    &               * e1v(ji,jj) * e2v(ji,jj) * e3v(ji,jj,jk) &
                    &               * vmask(ji,jj,jk)
               END DO
            END DO
         END DO
         !--------------------------------------------------------------------
         ! Compute the scalar product: ( L dx )^T W dy
         !--------------------------------------------------------------------

         zsp1 = DOT_PRODUCT( zua_tlout, zua_adin ) &
            & + DOT_PRODUCT( zva_tlout, zva_adin )

         !--------------------------------------------------------------------
         ! Call the adjoint routine: dx^* = L^T dy^*
         !--------------------------------------------------------------------
         ua_ad = zua_adin
         va_ad = zva_adin

         SELECT CASE (jt)
         CASE(1)
            nbab = 1
            ngib = 0
            nhor = 0
            CALL cla_div_adj( nit000 )
         CASE(2)
            nbab = 0
            ngib = 1
            nhor = 0
            CALL cla_div_adj( nit000 )
         CASE(3)
            nbab = 0
            ngib = 0
            nhor = 1
            CALL cla_div_adj( nit000 )
         END SELECT

         zua_adout   = ua_ad
         zva_adout   = va_ad

         zsp2 = DOT_PRODUCT( zua_tlin  , zua_adout   ) &
            & + DOT_PRODUCT( zva_tlin  , zva_adout   )

         SELECT CASE (jt)
         CASE(1)
            cl_name = 'cladynspg babm'
         CASE(2)
            cl_name = 'cladynspg gibr'
         CASE(3)
            cl_name = 'cladynspg horm'
         END SELECT
         CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )
      END DO

      DEALLOCATE(         &
         & zua_tlin,      &
         & zva_tlin,      &
         & zua_tlout,     &
         & zva_tlout,     &
         & zua_adin,      &
         & zva_adin,      &
         & zua_adout,     &
         & zva_adout,     &
         & zr3d           &
         & )

   END SUBROUTINE cla_dynspg_adj_tst

   SUBROUTINE cla_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE cla_adj_tst ***
      !!
      !! ** Purpose : Test the adjoint routine.
      !!
      !! ** Method  : call the test routines of all the cla_ routines
      !!
      !!
      !! History :
      !!        ! 08-08 (A. Vidard)
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt             ! Output unit

      CALL cla_div_adj_tst( kumadt )

      CALL cla_traadv_adj_tst( kumadt )
      
      CALL cla_dynspg_adj_tst( kumadt )

   END SUBROUTINE cla_adj_tst

END MODULE cla_tam
