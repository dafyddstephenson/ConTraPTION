MODULE lib_mpp_tam
   !!======================================================================
   !!                       ***  MODULE  lib_mpp_tam  ***
   !! Ocean numerics:  massively parallel processing library
   !!=====================================================================
   !! History :  OPA  !  1994  (M. Guyon, J. Escobar, M. Imbard)  Original code
   !!            7.0  !  1997  (A.M. Treguier)  SHMEM additions
   !!            8.0  !  1998  (M. Imbard, J. Escobar, L. Colombet ) SHMEM and MPI
   !!                 !  1998  (J.M. Molines) Open boundary conditions
   !!   NEMO     1.0  !  2003  (J.-M. Molines, G. Madec)  F90, free form
   !!                 !  2003  (J.M. Molines) add mpp_ini_north(_3d,_2d)
   !!             -   !  2004  (R. Bourdalle Badie)  isend option in mpi
   !!                 !  2004  (J.M. Molines) minloc, maxloc
   !!             -   !  2005  (G. Madec, S. Masson)  npolj=5,6 F-point & ice cases
   !!             -   !  2005  (R. Redler) Replacement of MPI_COMM_WORLD except for MPI_Abort
   !!             -   !  2005  (R. Benshila, G. Madec)  add extra halo case
   !!             -   !  2008  (R. Benshila) add mpp_ini_ice
   !!            3.2  !  2009  (R. Benshila) SHMEM suppression, north fold in lbc_nfd
   !!            3.2  !  2009  (O. Marti)    add mpp_ini_znl
   !!            4.0  !  2011  (G. Madec)  move ctl_ routines from in_out_manager
   !!   NEMOTAM  2.?  !  2007  (K. Mogensen) Original code (lib_mppadj)
   !!            3.0  !  2009  (A. Vidard) nemo v3 update
   !!            3.2  !  2010  (A. Vidard) 3.2 version, complete rewrite
   !!            3.4  !  2012  (P.-A. Bouttier) v3.4 update
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   ctl_stop   : update momentum and tracer Kz from a tke scheme
   !!   ctl_warn   : initialization, namelist read, and parameters control
   !!   ctl_opn    : Open file and check if required file is available.
   !!   get_unit    : give the index of an unused logical unit
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   'key_mpp_mpi'             MPI massively parallel processing library
   !!----------------------------------------------------------------------
   !!   mpp_lnk_adj       : interface (defined in lbclnk) for message passing of 2d or 3d arrays (mpp_lnk_2d, mpp_lnk_3d)
   !!   mpp_lnk_3d_gather_adj :  Message passing manadgement for two 3D arrays
   !!   mpp_lnk_e_adj     : interface (defined in lbclnk) for message passing of 2d array with extra halo (mpp_lnk_2d_e)
   !!   mpp_lbc_north_adj : north fold processors gathering
   !!   mpp_lbc_north_e_adj : variant of mpp_lbc_north for extra outer halo
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE lbcnfd_tam     ! north fold treatment
   USE lib_mpp
   USE in_out_manager ! I/O manager

   IMPLICIT NONE
   PRIVATE

   PUBLIC   mpp_lbc_north_adj, mpp_lbc_north_e_adj
   PUBLIC   mpp_lnk_3d_adj, mpp_lnk_3d_gather_adj, mpp_lnk_2d_adj, mpp_lnk_2d_e_adj
   PUBLIC   lib_mpp_alloc_adj
   PUBLIC   mpp_sum_nfd
   !! * Interfaces
   !! define generic interface for these routine as they are called sometimes
   !! with scalar arguments instead of array arguments, which causes problems
   !! for the compilation on AIX system as well as NEC and SGI. Ok on COMPACQ
   INTERFACE mpp_lbc_north_adj
      MODULE PROCEDURE mpp_lbc_north_3d_adj, mpp_lbc_north_2d_adj
   END INTERFACE

   INTERFACE mpp_sum_nfd
      MODULE PROCEDURE mpp_sum_nfd_3d,mpp_sum_nfd_4d
   END INTERFACE
   !! ========================= !!
   !!  MPI  variable definition !!
   !! ========================= !!
!$AGRIF_DO_NOT_TREAT
INCLUDE "mpif.h"
!$AGRIF_END_DO_NOT_TREAT

   ! message passing arrays
   REAL(wp), DIMENSION(:,:,:,:,:), ALLOCATABLE, SAVE ::   t4ns_ad, t4sn_ad   ! 2 x 3d for north-south & south-north
   REAL(wp), DIMENSION(:,:,:,:,:), ALLOCATABLE, SAVE ::   t4ew_ad, t4we_ad   ! 2 x 3d for east-west & west-east
   REAL(wp), DIMENSION(:,:,:,:,:), ALLOCATABLE, SAVE ::   t4p1_ad, t4p2_ad   ! 2 x 3d for north fold
   REAL(wp), DIMENSION(:,:,:,:)  , ALLOCATABLE, SAVE ::   t3ns_ad, t3sn_ad   ! 3d for north-south & south-north
   REAL(wp), DIMENSION(:,:,:,:)  , ALLOCATABLE, SAVE ::   t3ew_ad, t3we_ad   ! 3d for east-west & west-east
   REAL(wp), DIMENSION(:,:,:,:)  , ALLOCATABLE, SAVE ::   t3p1_ad, t3p2_ad   ! 3d for north fold
   REAL(wp), DIMENSION(:,:,:)    , ALLOCATABLE, SAVE ::   t2ns_ad, t2sn_ad   ! 2d for north-south & south-north
   REAL(wp), DIMENSION(:,:,:)    , ALLOCATABLE, SAVE ::   t2ew_ad, t2we_ad   ! 2d for east-west & west-east
   REAL(wp), DIMENSION(:,:,:)    , ALLOCATABLE, SAVE ::   t2p1_ad, t2p2_ad   ! 2d for north fold
   REAL(wp), DIMENSION(:,:,:)    , ALLOCATABLE, SAVE ::   tr2ns_ad, tr2sn_ad ! 2d for north-south & south-north + extra outer halo
   REAL(wp), DIMENSION(:,:,:)    , ALLOCATABLE, SAVE ::   tr2ew_ad, tr2we_ad ! 2d for east-west   & west-east   + extra outer halo

   ! Arrays used in mpp_lbc_north_3d()
   REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE, SAVE   ::   ztabad, znorthlocad
   REAL(wp), DIMENSION(:,:,:,:), ALLOCATABLE, SAVE   ::   znorthgloioad
   REAL(wp), DIMENSION(:,:,:)  , ALLOCATABLE, SAVE   ::   zfoldwkad      ! Workspace for message transfers avoiding mpi_allgather

   ! Arrays used in mpp_lbc_north_2d()
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, SAVE    ::   ztabad_2d, znorthlocad_2d
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, SAVE    ::   znorthgloioad_2d
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, SAVE    ::   zfoldwkad_2d    ! Workspace for message transfers avoiding mpi_allgather

   ! Arrays used in mpp_lbc_north_e()
   REAL(wp), DIMENSION(:,:)  , ALLOCATABLE, SAVE    ::   ztabad_e, znorthlocad_e
   REAL(wp), DIMENSION(:,:,:), ALLOCATABLE, SAVE    ::   znorthgloioad_e

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

INTEGER FUNCTION lib_mpp_alloc_adj( kumout )
      !!----------------------------------------------------------------------
      !!              ***  routine lib_mpp_alloc  ***
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kumout   ! ocean.output logical unit
      !!----------------------------------------------------------------------
      !
      ALLOCATE( t4ns_ad(jpi,jprecj,jpk,2,2) , t4sn_ad(jpi,jprecj,jpk,2,2) ,                                            &
         &      t4ew_ad(jpj,jpreci,jpk,2,2) , t4we_ad(jpj,jpreci,jpk,2,2) ,                                            &
         &      t4p1_ad(jpi,jprecj,jpk,2,2) , t4p2_ad(jpi,jprecj,jpk,2,2) ,                                            &
         &      t3ns_ad(jpi,jprecj,jpk,2)   , t3sn_ad(jpi,jprecj,jpk,2)   ,                                            &
         &      t3ew_ad(jpj,jpreci,jpk,2)   , t3we_ad(jpj,jpreci,jpk,2)   ,                                            &
         &      t3p1_Ad(jpi,jprecj,jpk,2)   , t3p2_ad(jpi,jprecj,jpk,2)   ,                                            &
         &      t2ns_ad(jpi,jprecj    ,2)   , t2sn_ad(jpi,jprecj    ,2)   ,                                            &
         &      t2ew_ad(jpj,jpreci    ,2)   , t2we_ad(jpj,jpreci    ,2)   ,                                            &
         &      t2p1_ad(jpi,jprecj    ,2)   , t2p2_ad(jpi,jprecj    ,2)   ,                                            &
         !
         &      tr2ns_ad(1-jpr2di:jpi+jpr2di,jprecj+jpr2dj,2) ,                                                     &
         &      tr2sn_ad(1-jpr2di:jpi+jpr2di,jprecj+jpr2dj,2) ,                                                     &
         &      tr2ew_ad(1-jpr2dj:jpj+jpr2dj,jpreci+jpr2di,2) ,                                                     &
         &      tr2we_ad(1-jpr2dj:jpj+jpr2dj,jpreci+jpr2di,2) ,                                                     &
         !
         &      ztabad(jpiglo,4,jpk) , znorthlocad(jpi,4,jpk) , znorthgloioad(jpi,4,jpk,jpni) ,                        &
         &      zfoldwkad(jpi,4,jpk) ,                                                                             &
         !
         &      ztabad_2d(jpiglo,4)  , znorthlocad_2d(jpi,4)  , znorthgloioad_2d(jpi,4,jpni)  ,                        &
         &      zfoldwkad_2d(jpi,4)  ,                                                                             &
         !
         &      ztabad_e(jpiglo,4+2*jpr2dj) , znorthlocad_e(jpi,4+2*jpr2dj) , znorthgloioad_e(jpi,4+2*jpr2dj,jpni) ,   &
         !
         &      STAT=lib_mpp_alloc_adj )
         !
      IF( lib_mpp_alloc_adj /= 0 ) THEN
         WRITE(kumout,cform_war)
         WRITE(kumout,*) 'lib_mpp_alloc_adj : failed to allocate arrays'
      ENDIF
      !
   END FUNCTION lib_mpp_alloc_adj

   SUBROUTINE mpp_lnk_3d_adj( ptab_ad, cd_type, psgn, cd_mpp, pval )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_3d_adj  ***
      !!
      !! ** Purpose :   Adjoint of Message passing manadgement
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !! ** Action  :   ptab_ad with update value at its periphery
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptab_ad  ! 3D array on which the boundary condition is applied
      CHARACTER(len=1)                , INTENT(in   ) ::   cd_type  ! define the nature of ptab_ad array grid-points
      !                                                             ! = T , U , V , F , W points
      REAL(wp)                        , INTENT(in   ) ::   psgn     ! =-1 the sign change across the north fold boundary
      !                                                             ! =  1. , the sign is kept
      CHARACTER(len=3), OPTIONAL      , INTENT(in   ) ::   cd_mpp   ! fill the overlap area only
      REAL(wp)        , OPTIONAL      , INTENT(in   ) ::   pval     ! background value (used at closed boundaries)
      !!
      INTEGER  ::   ji, jj, jk, jl             ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! temporary integers
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for key_mpi_isend
      REAL(wp) ::   zland
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for key_mpi_isend
      !!----------------------------------------------------------------------

      t3ns_ad = 0.0_wp ; t3sn_ad = 0.0_wp
      t3we_ad = 0.0_wp ; t3ew_ad = 0.0_wp
      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value
      ELSE                         ;   zland = 0.e0      ! zero by default
      ENDIF
      !
      ! 4. north fold treatment
      ! -----------------------
      !
      IF( npolj /= 0 .AND. .NOT. PRESENT(cd_mpp) ) THEN
         !
         SELECT CASE ( jpni )
         CASE ( 1 )     ;   CALL lbc_nfd_adj      ( ptab_ad, cd_type, psgn )   ! only 1 northern proc, no mpp
         CASE DEFAULT   ;   CALL mpp_lbc_north_adj( ptab_ad, cd_type, psgn )   ! for all northern procs.
         END SELECT
         !
      ENDIF
      !
      !
      ! 3. North and south directions
      ! -----------------------------
      ! always closed : we play only with the neigbours
      !
      !                           ! Write Dirichlet lateral conditions
      ijhom = nlcj-jprecj
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, jprecj
            t3ns_ad(:,jl,:,2) = t3ns_ad(:,jl,:,2) + ptab_ad(:,ijhom+jl,:)
            ptab_ad(:,ijhom+jl,:) = 0.0_wp
         END DO
      CASE ( 0 )
         DO jl = 1, jprecj
            t3ns_ad(:,jl,:,2) = t3ns_ad(:,jl,:,2) + ptab_ad(:,ijhom+jl,:)
            ptab_ad(:,ijhom+jl,:) = 0.0_wp
            t3sn_ad(:,jl,:,2) = t3sn_ad(:,jl,:,2) + ptab_ad(:,jl      ,:)
            ptab_ad(:,jl      ,:) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, jprecj
            t3sn_ad(:,jl,:,2) = t3sn_ad(:,jl,:,2) + ptab_ad(:,jl,:)
            ptab_ad(:,jl,:) = 0.0_wp
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = jprecj * jpi * jpk
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, t3ns_ad(1,1,1,2), imigr, nono, ml_req1 )
         CALL mpprecv( 3, t3sn_ad(1,1,1,1), imigr, nono )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 3, t3sn_ad(1,1,1,2), imigr, noso, ml_req1 )
         CALL mppsend( 4, t3ns_ad(1,1,1,2), imigr, nono, ml_req2 )
         CALL mpprecv( 3, t3sn_ad(1,1,1,1), imigr, nono )
         CALL mpprecv( 4, t3ns_ad(1,1,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 )
         CALL mppsend( 3, t3sn_ad(1,1,1,2), imigr, noso, ml_req1 )
         CALL mpprecv( 4, t3ns_ad(1,1,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      !
      !
      IF( nbondj /= 2 ) THEN      ! Read Dirichlet lateral conditions
         ijhom = nlcj-nrecj
         DO jl = 1, jprecj
            ptab_ad(:,jprecj+jl,:) = ptab_ad(:,jprecj+jl,:) + t3ns_ad(:,jl,:,1)
            t3ns_ad(:,jl,:,1) = 0.0_wp
            ptab_ad(:,ijhom +jl,:) = ptab_ad(:,ijhom +jl,:) + t3sn_ad(:,jl,:,1)
            t3sn_ad(:,jl,:,1) = 0.0_wp
         END DO
      ENDIF
      !
      ! 2. East and west directions exchange
      ! ------------------------------------
      !                           ! Write Dirichlet lateral conditions
      iihom = nlci-jpreci
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, jpreci
            t3ew_ad(:,jl,:,2)     = t3ew_ad(:,jl,:,2) + ptab_ad(iihom+jl,:,:)
            ptab_ad(iihom+jl,:,:) = 0.0_wp
         END DO
      CASE ( 0 )
         DO jl = 1, jpreci
            t3ew_ad(:,jl,:,2)     = t3ew_ad(:,jl,:,2) + ptab_ad(iihom+jl,:,:)
            ptab_ad(iihom+jl,:,:) = 0.0_wp
            t3we_ad(:,jl,:,2)     = t3we_ad(:,jl,:,2) + ptab_ad(jl      ,:,:)
            ptab_ad(jl      ,:,:) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, jpreci
            t3we_ad(:,jl,:,2)     = t3we_ad(:,jl,:,2) + ptab_ad(jl      ,:,:)
            ptab_ad(jl      ,:,:) = 0.0_wp
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = jpreci * jpj * jpk
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, t3ew_ad(1,1,1,2), imigr, noea, ml_req1 )
         CALL mpprecv( 1, t3we_ad(1,1,1,1), imigr, noea )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 1, t3we_ad(1,1,1,2), imigr, nowe, ml_req1 )
         CALL mppsend( 2, t3ew_ad(1,1,1,2), imigr, noea, ml_req2 )
         CALL mpprecv( 1, t3we_ad(1,1,1,1), imigr, noea )
         CALL mpprecv( 2, t3ew_ad(1,1,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 )
         CALL mppsend( 1, t3we_ad(1,1,1,2), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, t3ew_ad(1,1,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      !
      ! we play with the neigbours AND the row number because of the periodicity
      !
      SELECT CASE ( nbondi )      ! Read Dirichlet lateral conditions
      CASE ( -1, 0, 1 )                ! all exept 2 (i.e. close case)
         iihom = nlci-nreci
         DO jl = 1, jpreci
            ptab_ad(iihom +jl,:,:) = ptab_ad(iihom +jl,:,:) + t3we_ad(:,jl,:,1)
            t3we_ad(:,jl,:,1) = 0.0_wp
            ptab_ad(jpreci+jl,:,:) = ptab_ad(jpreci+jl,:,:) + t3ew_ad(:,jl,:,1)
            t3ew_ad(:,jl,:,1) = 0.0_wp
         END DO
      END SELECT
      !
      ! 1. standard boundary treatment
      ! ------------------------------
      IF( PRESENT( cd_mpp ) ) THEN      ! only fill added line/raw with existing values
         !
         ! WARNING ptab_ad is defined only between nld and nle
         DO jk = 1, jpk
            DO ji = nlci+1, jpi                 ! added column(s) (full)
               ptab_ad(nlei,      nlej   ,jk) = ptab_ad(nlei,nlej     ,jk) + SUM(ptab_ad(ji,nlej+1:jpj   ,jk))
               ptab_ad(ji  ,nlej+1:jpj   ,jk) = 0.0_wp
               ptab_ad(nlei,nldj         ,jk) = ptab_ad(nlei,nldj     ,jk) + SUM(ptab_ad(ji,1     :nldj-1,jk))
               ptab_ad(ji  ,1     :nldj-1,jk) = 0.0_wp
               ptab_ad(nlei,nldj  :nlej  ,jk) = ptab_ad(nlei,nldj:nlej,jk) + ptab_ad(ji,nldj  :nlej  ,jk)
               ptab_ad(ji  ,nldj  :nlej  ,jk) = 0.0_wp
            END DO
            DO jj = nlcj+1, jpj                 ! added line(s)   (inner only)
               ptab_ad(nlei         ,nlej,jk) = ptab_ad(nlei     ,nlej,jk) + SUM(ptab_ad(nlei+1:nlci  ,jj,jk))
               ptab_ad(nlei+1:nlci  ,jj  ,jk) = 0.0_wp
               ptab_ad(nldi         ,nlej,jk) = ptab_ad(nldi     ,nlej,jk) + SUM(ptab_ad(1     :nldi-1,jj,jk))
               ptab_ad(1     :nldi-1,jj  ,jk) = 0.0_wp
               ptab_ad(nldi  :nlei  ,nlej,jk) = ptab_ad(nldi:nlei,nlej,jk) + ptab_ad(nldi  :nlei  ,jj,jk)
               ptab_ad(nldi  :nlei  ,jj  ,jk) = 0.0_wp
            END DO
         END DO
         !
      ELSE                              ! standard close or cyclic treatment
      !
      !                                   ! North-South boundaries (always closed)
                                   ptab_ad(:,nlcj-jprecj+1:jpj   ,:) = 0.0_wp       ! north
      IF( .NOT. cd_type == 'F' )   ptab_ad(:,     1       :jprecj,:) = 0.0_wp       ! south except F-point
      !                                   ! East-West boundaries
      !                                        !* Cyclic east-west
      IF( nbondi == 2 .AND. (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
         ptab_ad( 2   ,:,:) = ptab_ad(  2  ,:,:) + ptab_ad(jpi,:,:)
         ptab_ad(jpi  ,:,:) = 0.0_wp
         ptab_ad(jpim1,:,:) = ptab_ad(jpim1,:,:) + ptab_ad( 1 ,:,:)
         ptab_ad( 1   ,:,:) = 0.0_wp
      ELSE                                     !* closed
         IF( .NOT. cd_type == 'F' )   ptab_ad(     1       :jpreci,:,:) = 0.0_wp    ! south except F-point
                                      ptab_ad(nlci-jpreci+1:jpi   ,:,:) = 0.0_wp    ! north
      ENDIF
      !
   ENDIF
      !
   END SUBROUTINE mpp_lnk_3d_adj


   SUBROUTINE mpp_lnk_2d_adj( pt2d_ad, cd_type, psgn, cd_mpp, pval )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_adj  ***
      !!
      !! ** Purpose :   Adjoint of Message passing manadgement for 2d array
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pt2d_ad     ! 2D array on which the boundary condition is applied
      CHARACTER(len=1)            , INTENT(in   ) ::   cd_type  ! define the nature of ptab_ad array grid-points
      !                                                         ! = T , U , V , F , W and I points
      REAL(wp)                    , INTENT(in   ) ::   psgn     ! =-1 the sign change across the north fold boundary
      !                                                         ! =  1. , the sign is kept
      CHARACTER(len=3), OPTIONAL  , INTENT(in   ) ::   cd_mpp   ! fill the overlap area only
      REAL(wp)        , OPTIONAL  , INTENT(in   ) ::   pval     ! background value (used at closed boundaries)
      !!
      INTEGER  ::   ji, jj, jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! temporary integers
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for key_mpi_isend
      REAL(wp) ::   zland
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for key_mpi_isend
      !!----------------------------------------------------------------------
      t2ns_ad = 0.0_wp ; t2sn_ad = 0.0_wp
      t2we_ad = 0.0_wp ; t2ew_ad = 0.0_wp
      IF( PRESENT( pval ) ) THEN   ;   zland = pval      ! set land value
      ELSE                         ;   zland = 0.e0      ! zero by default
      ENDIF
      !4. north fold treatment
      !-----------------------

      IF( npolj /= 0 .AND. .NOT. PRESENT(cd_mpp) ) THEN
         !
         SELECT CASE ( jpni )
         CASE ( 1 )
            CALL lbc_nfd_adj      ( pt2d_ad, cd_type, psgn )   ! only 1 northern proc, no mpp
         CASE DEFAULT
            CALL mpp_lbc_north_adj( pt2d_ad, cd_type, psgn )   ! for all northern procs.
         END SELECT
         !
      ENDIF

      !3. North and south directions
      !-----------------------------
                                !! Write Dirichlet lateral conditions
      ijhom = nlcj - jprecj
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, jprecj
            t2ns_ad(:,jl,2)     = t2ns_ad(:,jl,2) + pt2d_ad(:,ijhom+jl)
            pt2d_ad(:,ijhom+jl) = t2ns_ad(:,jl,2)
         END DO
      CASE ( 0 )
         DO jl = 1, jprecj
            t2ns_ad(:,jl,2)     = t2ns_ad(:,jl,2) + pt2d_ad(:,ijhom+jl)
            pt2d_ad(:,ijhom+jl) = 0.0_wp
            t2sn_ad(:,jl,2)     = t2sn_ad(:,jl,2) + pt2d_ad(:,jl      )
            pt2d_ad(:,jl      ) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, jprecj
            t2sn_ad(:,jl,2) = t2sn_ad(:,jl,2) + pt2d_ad(:,jl      )
            pt2d_ad(:,jl      ) = 0.0_wp
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = jprecj * jpi
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, t2ns_ad(1,1,2), imigr, nono, ml_req1 )
         CALL mpprecv( 3, t2sn_ad(1,1,1), imigr, nono )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 3, t2sn_ad(1,1,2), imigr, noso, ml_req1 )
         CALL mppsend( 4, t2ns_ad(1,1,2), imigr, nono, ml_req2 )
         CALL mpprecv( 3, t2sn_ad(1,1,1), imigr, nono )
         CALL mpprecv( 4, t2ns_ad(1,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 3, t2sn_ad(1,1,2), imigr, noso, ml_req1 )
         CALL mpprecv( 4, t2ns_ad(1,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT

      !always closed : we play only with the neigbours

      IF( nbondj /= 2 ) THEN      ! Read Dirichlet lateral conditions
         ijhom = nlcj-nrecj
         DO jl = 1, jprecj
            pt2d_ad(:,jprecj+jl) = pt2d_ad(:,jprecj+jl) + t2ns_ad(:,jl,1)
            t2ns_ad(:,jl,1) = 0.0_wp
            pt2d_ad(:,ijhom +jl) = pt2d_ad(:,ijhom +jl) + t2sn_ad(:,jl,1)
            t2sn_ad(:,jl,1) = 0.0_wp
         END DO
      ENDIF
      !2. East and west directions exchange
      !------------------------------------
                                !! Write Dirichlet lateral conditions
      iihom = nlci - jpreci
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, jpreci
            t2ew_ad(:,jl,2)     = t2ew_ad(:,jl,2) + pt2d_ad(iihom+jl,:)
            pt2d_ad(iihom+jl,:) = 0.0_wp
         END DO
      CASE ( 0 )
         DO jl = 1, jpreci
            t2ew_ad(:,jl,2)     = t2ew_ad(:,jl,2) + pt2d_ad(iihom+jl,:)
            pt2d_ad(iihom+jl,:) = 0.0_wp
            t2we_ad(:,jl,2)     = t2we_ad(:,jl,2) + pt2d_ad(jl      ,:)
            pt2d_ad(jl      ,:) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, jpreci
            t2we_ad(:,jl,2) = t2we_ad(:,jl,2) + pt2d_ad(jl,:)
            pt2d_ad(jl,:)   = 0.0_wp
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = jpreci * jpj
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, t2ew_ad(1,1,2), imigr, noea, ml_req1 )
         CALL mpprecv( 1, t2we_ad(1,1,1), imigr, noea )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 1, t2we_ad(1,1,2), imigr, nowe, ml_req1 )
         CALL mppsend( 2, t2ew_ad(1,1,2), imigr, noea, ml_req2 )
         CALL mpprecv( 1, t2we_ad(1,1,1), imigr, noea )
         CALL mpprecv( 2, t2ew_ad(1,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 1, t2we_ad(1,1,2), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, t2ew_ad(1,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      ! we play with the neigbours AND the row number because of the periodicity
      !
      SELECT CASE ( nbondi )      ! Read Dirichlet lateral conditions
      CASE ( -1, 0, 1 )                ! all exept 2 (i.e. close case)
         iihom = nlci-nreci
         DO jl = 1, jpreci
            pt2d_ad(iihom +jl,:) = pt2d_ad(iihom +jl,:) + t2we_ad(:,jl,1)
            t2we_ad(:,jl,1)      = 0.0_wp
            pt2d_ad(jpreci+jl,:) = pt2d_ad(jpreci+jl,:) + t2ew_ad(:,jl,1)
            t2ew_ad(:,jl,1)      = 0.0_wp
         END DO
      END SELECT
      !
      ! 1. standard boundary treatment
      ! ------------------------------
      !
      IF( PRESENT( cd_mpp ) ) THEN      ! only fill added line/raw with existing values
         !
         ! WARNING pt2d is defined only between nld and nle
         DO ji = nlci+1, jpi                 ! added column(s) (full)
            pt2d_ad(nlei,       nlej  ) = pt2d_ad(nlei,     nlej) + SUM(pt2d_ad(ji,nlej+1:jpj   ))
            pt2d_ad(ji  ,nlej+1:jpj   ) = 0.0_wp
            pt2d_ad(nlei,nldj         ) = pt2d_ad(nlei,nldj     ) + SUM(pt2d_ad(ji,1     :nldj-1))
            pt2d_ad(ji  ,1     :nldj-1) = 0.0_wp
            pt2d_ad(nlei,nldj  :nlej  ) = pt2d_ad(nlei,nldj:nlej) + pt2d_ad(ji,nldj  :nlej  )
            pt2d_ad(ji  ,nldj  :nlej  ) = 0.0_wp
         END DO
         DO jj = nlcj+1, jpj                 ! added line(s)   (inner only)
            pt2d_ad(nlei         ,nlej) = pt2d_ad(     nlei,     nlej) + SUM(pt2d_ad(nlei+1:nlci  , jj))
            pt2d_ad(nlei+1:nlci  ,jj  ) = 0.0_wp
            pt2d_ad(nldi         ,nlej) = pt2d_ad(nldi     ,     nlej) + SUM(pt2d_ad(1     :nldi-1, jj))
            pt2d_ad(1     :nldi-1,jj  ) = 0.0_wp
            pt2d_ad(nldi  :nlei  ,nlej) = pt2d_ad(nldi:nlei,     nlej) + pt2d_ad(nldi  :nlei  , jj)
            pt2d_ad(nldi  :nlei  , jj ) = 0.0_wp
         END DO
         !
      ELSE                              ! standard close or cyclic treatment
         !
         !                                   ! North-South boundaries (always closed)
            IF( .NOT. cd_type == 'F' )   pt2d_ad(:,     1       :jprecj) = 0.0_wp !south except F-point
                                         pt2d_ad(:,nlcj-jprecj+1:jpj   ) = 0.0_wp ! north
         !                                   ! East-West boundaries
         IF( nbondi == 2 .AND.   &                ! Cyclic east-west
            &    (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
            pt2d_ad(  2  ,:) = pt2d_ad(  2  ,:) + pt2d_ad(jpi,:) ! east
            pt2d_ad(jpi  ,:) = 0.0_wp
            pt2d_ad(jpim1,:) = pt2d_ad(jpim1,:) + pt2d_ad( 1 ,:) ! west
            pt2d_ad( 1   ,:) = 0.0_wp
         ELSE                                     ! closed
            IF( .NOT. cd_type == 'F' )   pt2d_ad(     1       :jpreci,:) = 0.0_wp   ! south except F-point
                                         pt2d_ad(nlci-jpreci+1:jpi   ,:) = 0.0_wp   ! north
         ENDIF
         !
      ENDIF

   END SUBROUTINE mpp_lnk_2d_adj


   SUBROUTINE mpp_lnk_3d_gather_adj( ptab1_ad, cd_type1, ptab2_ad, cd_type2, psgn )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_3d_gather_adj  ***
      !!
      !! ** Purpose :   Adjoint of Message passing manadgement for two 3D arrays
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !! ** Action  :   ptab1_ad and ptab2_ad  with update value at its periphery
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptab1_ad     ! first and second 3D array on which
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   ptab2_ad     ! the boundary condition is applied
      CHARACTER(len=1)                , INTENT(in   ) ::   cd_type1  ! nature of ptab1_ad and ptab2_ad arrays
      CHARACTER(len=1)                , INTENT(in   ) ::   cd_type2  ! i.e. grid-points = T , U , V , F or W points
      REAL(wp)                        , INTENT(in   ) ::   psgn      ! =-1 the sign change across the north fold boundary
      !!                                                             ! =  1. , the sign is kept
      INTEGER  ::   jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom        ! temporary integers
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for key_mpi_isend
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for key_mpi_isend
      !!----------------------------------------------------------------------
      t4ns_ad = 0.0_wp ; t4sn_ad = 0.0_wp
      t4we_ad = 0.0_wp ; t4ew_ad = 0.0_wp
      ! 4. north fold treatment
      ! -----------------------
      IF( npolj /= 0 ) THEN
         !
         SELECT CASE ( jpni )
         CASE ( 1 )
            CALL lbc_nfd_adj      ( ptab2_ad, cd_type2, psgn )
            CALL lbc_nfd_adj      ( ptab1_ad, cd_type1, psgn )   ! only for northern procs.
         CASE DEFAULT
            CALL mpp_lbc_north_adj( ptab2_ad, cd_type2, psgn)
            CALL mpp_lbc_north_adj( ptab1_ad, cd_type1, psgn )   ! for all northern procs.
         END SELECT
         !
      ENDIF
      !
      ! 3. North and south directions
      ! -----------------------------
      !                           ! Write Dirichlet lateral conditions
      ijhom = nlcj - jprecj
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, jprecj
            t4ns_ad(:,jl,:,2,2)    = t4ns_ad(:,jl,:,2,2) + ptab2_ad(:,ijhom+jl,:)
            ptab2_ad(:,ijhom+jl,:) = 0.0_wp
            t4ns_ad(:,jl,:,1,2)    = t4ns_ad(:,jl,:,1,2) + ptab1_ad(:,ijhom+jl,:)
            ptab1_ad(:,ijhom+jl,:) = 0.0_wp
         END DO
      CASE ( 0 )
         DO jl = 1, jprecj
            t4ns_ad(:,jl,:,2,2)    = t4ns_ad(:,jl,:,2,2) + ptab2_ad(:,ijhom+jl,:)
            ptab2_ad(:,ijhom+jl,:) = 0.0_wp
            t4sn_ad(:,jl,:,2,2)    = t4sn_ad(:,jl,:,2,2) + ptab2_ad(:,jl      ,:)
            ptab2_ad(:,jl      ,:) = 0.0_wp
            t4ns_ad(:,jl,:,1,2)    = t4ns_ad(:,jl,:,1,2) + ptab1_ad(:,ijhom+jl,:)
            ptab1_ad(:,ijhom+jl,:) = 0.0_wp
            t4sn_ad(:,jl,:,1,2)    = t4sn_ad(:,jl,:,1,2) + ptab1_ad(:,jl      ,:)
            ptab1_ad(:,jl      ,:) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, jprecj
            t4sn_ad(:,jl,:,2,2) = t4sn_ad(:,jl,:,2,2) + ptab2_ad(:,jl,:)
            ptab2_ad(:,jl,:) = 0.0_wp
            t4sn_ad(:,jl,:,1,2) = t4sn_ad(:,jl,:,1,2) + ptab1_ad(:,jl,:)
            ptab1_ad(:,jl,:) = 0.0_wp
         END DO
      END SELECT
      !                           ! Migrations
      imigr = jprecj * jpi * jpk * 2
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, t4ns_ad(1,1,1,1,2), imigr, nono, ml_req1 )
         CALL mpprecv( 3, t4sn_ad(1,1,1,1,1), imigr, nono )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 3, t4sn_ad(1,1,1,1,2), imigr, noso, ml_req1 )
         CALL mppsend( 4, t4ns_ad(1,1,1,1,2), imigr, nono, ml_req2 )
         CALL mpprecv( 3, t4sn_ad(1,1,1,1,1), imigr, nono )
         CALL mpprecv( 4, t4ns_ad(1,1,1,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 )
         CALL mppsend( 3, t4sn_ad(1,1,1,1,2), imigr, noso, ml_req1 )
         CALL mpprecv( 4, t4ns_ad(1,1,1,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      !
      ! always closed : we play only with the neigbours
      !
      IF( nbondj /= 2 ) THEN      ! Read Dirichlet lateral conditions
         ijhom = nlcj - nrecj
         DO jl = 1, jprecj
            ptab2_ad(:,jprecj+jl,:) = ptab2_ad(:,jprecj+jl,:) + t4ns_ad(:,jl,:,2,1)
            t4ns_ad(:,jl,:,2,1)     = 0.0_wp
            ptab2_ad(:,ijhom +jl,:) = ptab2_ad(:,ijhom +jl,:) + t4sn_ad(:,jl,:,2,1)
            t4sn_ad(:,jl,:,2,1)     = 0.0_wp
            ptab1_ad(:,jprecj+jl,:) = ptab1_ad(:,jprecj+jl,:) + t4ns_ad(:,jl,:,1,1)
            t4ns_ad(:,jl,:,1,1)     = 0.0_wp
            ptab1_ad(:,ijhom +jl,:) = ptab1_ad(:,ijhom +jl,:) + t4sn_ad(:,jl,:,1,1)
            t4sn_ad(:,jl,:,1,1)     = 0.0_wp
         END DO
      ENDIF
      !
      ! 2. East and west directions exchange
      ! ------------------------------------
      !                           ! Write Dirichlet lateral conditions
      iihom = nlci - jpreci
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, jpreci
            t4ew_ad(:,jl,:,2,2)    = t4ew_ad(:,jl,:,2,2) + ptab2_ad(iihom+jl,:,:)
            ptab2_ad(iihom+jl,:,:) = 0.0_wp
            t4ew_ad(:,jl,:,1,2)    = t4ew_ad(:,jl,:,1,2) + ptab1_ad(iihom+jl,:,:)
            ptab1_ad(iihom+jl,:,:) = 0.0_wp
         END DO
      CASE ( 0 )
         DO jl = 1, jpreci
            t4ew_ad(:,jl,:,2,2)    = t4ew_ad(:,jl,:,2,2) + ptab2_ad(iihom+jl,:,:)
            ptab2_ad(iihom+jl,:,:) = 0.0_wp
            t4we_ad(:,jl,:,2,2)    = t4we_ad(:,jl,:,2,2) + ptab2_ad(jl      ,:,:)
            ptab2_ad(jl      ,:,:) = 0.0_wp
            t4ew_ad(:,jl,:,1,2)    = t4ew_ad(:,jl,:,1,2) + ptab1_ad(iihom+jl,:,:)
            ptab1_ad(iihom+jl,:,:) = 0.0_wp
            t4we_ad(:,jl,:,1,2)    = t4we_ad(:,jl,:,1,2) + ptab1_ad(jl      ,:,:)
            ptab1_ad(jl      ,:,:) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, jpreci
            t4we_ad(:,jl,:,2,2) = t4we_ad(:,jl,:,2,2) + ptab2_ad(jl      ,:,:)
            t4we_ad(:,jl,:,1,2) = t4we_ad(:,jl,:,1,2) + ptab1_ad(jl      ,:,:)
            ptab1_ad(jl      ,:,:) = 0.0_wp
            ptab2_ad(jl      ,:,:) = 0.0_wp
         END DO
      END SELECT
      !
      !                           ! Migrations
      imigr = jpreci * jpj * jpk *2
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, t4ew_ad(1,1,1,1,2), imigr, noea, ml_req1 )
         CALL mpprecv( 1, t4we_ad(1,1,1,1,1), imigr, noea )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      CASE ( 0 )
         CALL mppsend( 1, t4we_ad(1,1,1,1,2), imigr, nowe, ml_req1 )
         CALL mppsend( 2, t4ew_ad(1,1,1,1,2), imigr, noea, ml_req2 )
         CALL mpprecv( 1, t4we_ad(1,1,1,1,1), imigr, noea )
         CALL mpprecv( 2, t4ew_ad(1,1,1,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2, ml_stat, ml_err)
      CASE ( 1 )
         CALL mppsend( 1, t4we_ad(1,1,1,1,2), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, t4ew_ad(1,1,1,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1, ml_stat, ml_err)
      END SELECT
      !
      ! we play with the neigbours AND the row number because of the periodicity
      !
      SELECT CASE ( nbondi )      ! Read Dirichlet lateral conditions
      CASE ( -1, 0, 1 )                ! all exept 2 (i.e. close case)
         iihom = nlci-nreci
         DO jl = 1, jpreci
            ptab2_ad(iihom +jl,:,:) = ptab2_ad(iihom +jl,:,:) + t4we_ad(:,jl,:,2,1)
            t4we_ad(:,jl,:,2,1)     = 0.0_wp
            ptab2_ad(jpreci+jl,:,:) = ptab2_ad(jpreci+jl,:,:) + t4ew_ad(:,jl,:,2,1)
            t4ew_ad(:,jl,:,2,1)     = 0.0_wp
            ptab1_ad(iihom +jl,:,:) = ptab1_ad(iihom +jl,:,:) + t4we_ad(:,jl,:,1,1)
            t4we_ad(:,jl,:,1,1)     = 0.0_wp
            ptab1_ad(jpreci+jl,:,:) = ptab1_ad(jpreci+jl,:,:) + t4ew_ad(:,jl,:,1,1)
            t4ew_ad(:,jl,:,1,1)     = 0.0_wp
         END DO
      END SELECT
      ! 1. standard boundary treatment
      ! ------------------------------
      !                                      ! East-West boundaries
      !                                           !* Cyclic east-west
      !                                      ! North-South boundaries
                                    ptab2_ad(:,nlcj-jprecj+1:jpj   ,:) = 0.e0
                                    ptab1_ad(:,nlcj-jprecj+1:jpj   ,:) = 0.e0    ! north
      IF( .NOT. cd_type2 == 'F' )   ptab2_ad(:,     1       :jprecj,:) = 0.e0
      IF( .NOT. cd_type1 == 'F' )   ptab1_ad(:,     1       :jprecj,:) = 0.e0    ! south except at F-point
      IF( nbondi == 2 .AND. (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
         ptab2_ad(  2  ,:,:) = ptab2_ad(  2  ,:,:) + ptab2_ad(jpi,:,:)
         ptab2_ad(jpi,:,:)   = 0.0_wp
         ptab2_ad(jpim1,:,:) = ptab2_ad(jpim1,:,:) + ptab2_ad( 1 ,:,:)
         ptab2_ad( 1 ,:,:)   = 0.0_wp
         ptab1_ad(  2  ,:,:) = ptab1_ad(  2  ,:,:) + ptab1_ad(jpi,:,:)
         ptab1_ad(jpi,:,:)   = 0.0_wp
         ptab1_ad(jpim1,:,:) = ptab1_ad(jpim1,:,:) + ptab1_ad( 1 ,:,:)
         ptab1_ad( 1 ,:,:)   = 0.0_wp
      ELSE                                        !* closed
         IF( .NOT. cd_type1 == 'F' )   ptab1_ad(     1       :jpreci,:,:) = 0.e0    ! south except at F-point
         IF( .NOT. cd_type2 == 'F' )   ptab2_ad(     1       :jpreci,:,:) = 0.e0
                                       ptab1_ad(nlci-jpreci+1:jpi   ,:,:) = 0.e0    ! north
                                       ptab2_ad(nlci-jpreci+1:jpi   ,:,:) = 0.e0
      ENDIF
   END SUBROUTINE mpp_lnk_3d_gather_adj


   SUBROUTINE mpp_lnk_2d_e_adj( pt2d_ad, cd_type, psgn )
      !!----------------------------------------------------------------------
      !!                  ***  routine mpp_lnk_2d_e_adj ***
      !!
      !! ** Purpose :   Adjoint of Message passing manadgement for 2d array (with halo)
      !!
      !! ** Method  :   Use mppsend and mpprecv function for passing mask
      !!      between processors following neighboring subdomains.
      !!            domain parameters
      !!                    nlci   : first dimension of the local subdomain
      !!                    nlcj   : second dimension of the local subdomain
      !!                    jpr2di : number of rows for extra outer halo
      !!                    jpr2dj : number of columns for extra outer halo
      !!                    nbondi : mark for "east-west local boundary"
      !!                    nbondj : mark for "north-south local boundary"
      !!                    noea   : number for local neighboring processors
      !!                    nowe   : number for local neighboring processors
      !!                    noso   : number for local neighboring processors
      !!                    nono   : number for local neighboring processors
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT(inout) ::   pt2d_ad  ! 2D array with extra halo
      CHARACTER(len=1)                                            , INTENT(in   ) ::   cd_type  ! nature of ptab_ad array grid-points
      !                                                                                         ! = T , U , V , F , W and I points
      REAL(wp)                                                    , INTENT(in   ) ::   psgn     ! =-1 the sign change across the
      !!                                                                                        ! north boundary, =  1. otherwise
      INTEGER  ::   jl   ! dummy loop indices
      INTEGER  ::   imigr, iihom, ijhom, rank, ierr  ! temporary integers
      INTEGER  ::   ipreci, iprecj             ! temporary integers
      INTEGER  ::   ml_req1, ml_req2, ml_err   ! for key_mpi_isend
      INTEGER, DIMENSION(MPI_STATUS_SIZE) ::   ml_stat   ! for key_mpi_isend
      !!----------------------------------------------------------------------
      ipreci = jpreci + jpr2di      ! take into account outer extra 2D overlap area
      iprecj = jprecj + jpr2dj
      tr2ns_ad(:,:,:) = 0.0_wp ; tr2sn_ad(:,:,:) = 0.0_wp
      tr2we_ad(:,:,:) = 0.0_wp ; tr2ew_ad(:,:,:) = 0.0_wp
      ! 3. North and south directions
      ! -----------------------------
      !                           ! Write Dirichlet lateral conditions
      ijhom = nlcj - jprecj
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         DO jl = 1, iprecj
            tr2ns_ad(:,jl,2) = tr2ns_ad(:,jl,2) + pt2d_ad(:,ijhom+jl)
            pt2d_ad(:,ijhom+jl) = 0.0_wp
         END DO
      CASE ( 0 )
         DO jl = 1, iprecj
            tr2ns_ad(:,jl,2)     = tr2ns_ad(:,jl,2) + pt2d_ad(:,ijhom+jl )
            pt2d_ad(:,ijhom+jl ) = 0.0_wp
            tr2sn_ad(:,jl,2)     = tr2sn_ad(:,jl,2) + pt2d_ad(:,jl-jpr2dj)
            pt2d_ad(:,jl-jpr2dj) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, iprecj
            tr2sn_ad(:,jl,2) = tr2sn_ad(:,jl,2) + pt2d_ad(:,jl-jpr2dj)
            pt2d_ad(:,jl-jpr2dj) = 0.0_wp
         END DO
      END SELECT
      !                           ! Migrations
      imigr = iprecj * ( jpi + 2*jpr2di )
      !
      SELECT CASE ( nbondj )
      CASE ( -1 )
         CALL mppsend( 4, tr2ns_ad(1-jpr2di,1,2), imigr, nono, ml_req1 )
         CALL mpprecv( 3, tr2sn_ad(1-jpr2di,1,1), imigr, nono )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 3, tr2sn_ad(1-jpr2di,1,2), imigr, noso, ml_req1 )
         CALL mppsend( 4, tr2ns_ad(1-jpr2di,1,2), imigr, nono, ml_req2 )
         CALL mpprecv( 3, tr2sn_ad(1-jpr2di,1,1), imigr, nono )
         CALL mpprecv( 4, tr2ns_ad(1-jpr2di,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 3, tr2sn_ad(1-jpr2di,1,2), imigr, noso, ml_req1 )
         CALL mpprecv( 4, tr2ns_ad(1-jpr2di,1,1), imigr, noso )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      ! always closed : we play only with the neigbours
      !
      IF( nbondj /= 2 ) THEN      ! Read Dirichlet lateral conditions
         ijhom = nlcj-nrecj-jpr2dj
         DO jl = 1, iprecj
            pt2d_ad(:,jprecj+jl) = pt2d_ad(:,jprecj+jl) + tr2ns_ad(:,jl,1)
            tr2ns_ad(:,jl,1) = 0.0_wp
            pt2d_ad(:,ijhom +jl) = pt2d_ad(:,ijhom +jl) + tr2sn_ad(:,jl,1)
            tr2sn_ad(:,jl,1) = 0.0_wp
         END DO
      ENDIF
      !
      ! 2. East and west directions exchange
      ! ------------------------------------
      !                           ! Write Dirichlet lateral conditions
      iihom = nlci - jpreci
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         DO jl = 1, ipreci
            tr2ew_ad(:,jl,2)    = tr2ew_ad(:,jl,2) + pt2d_ad(iihom+jl,:)
            pt2d_ad(iihom+jl,:) = 0.0_wp
         END DO
      CASE ( 0 )
         DO jl = 1, ipreci
            tr2ew_ad(:,jl,2)  = tr2ew_ad(:,jl,2) + pt2d_ad( iihom+jl,:)
            pt2d_ad( iihom+jl,:) = 0.0_wp
            tr2we_ad(:,jl,2)  = tr2we_ad(:,jl,2) + pt2d_ad(jl-jpr2di,:)
            pt2d_ad(jl-jpr2di,:) = 0.0_wp
         END DO
      CASE ( 1 )
         DO jl = 1, ipreci
            tr2we_ad(:,jl,2)     = tr2we_ad(:,jl,2) + pt2d_ad(jl-jpr2di,:)
            pt2d_ad(jl-jpr2di,:) = 0.0_wp
         END DO
      END SELECT
      !                           ! Migrations
      imigr = ipreci * ( jpj + 2*jpr2dj)
      !
      SELECT CASE ( nbondi )
      CASE ( -1 )
         CALL mppsend( 2, tr2ew_ad(1-jpr2dj,1,2), imigr, noea, ml_req1 )
         CALL mpprecv( 1, tr2we_ad(1-jpr2dj,1,1), imigr, noea )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      CASE ( 0 )
         CALL mppsend( 1, tr2we_ad(1-jpr2dj,1,2), imigr, nowe, ml_req1 )
         CALL mppsend( 2, tr2ew_ad(1-jpr2dj,1,2), imigr, noea, ml_req2 )
         CALL mpprecv( 1, tr2we_ad(1-jpr2dj,1,1), imigr, noea )
         CALL mpprecv( 2, tr2ew_ad(1-jpr2dj,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
         IF(l_isend) CALL mpi_wait(ml_req2,ml_stat,ml_err)
      CASE ( 1 )
         CALL mppsend( 1, tr2we_ad(1-jpr2dj,1,2), imigr, nowe, ml_req1 )
         CALL mpprecv( 2, tr2ew_ad(1-jpr2dj,1,1), imigr, nowe )
         IF(l_isend) CALL mpi_wait(ml_req1,ml_stat,ml_err)
      END SELECT
      !
      ! we play with the neigbours AND the row number because of the periodicity
      !
      SELECT CASE ( nbondi )      ! Read Dirichlet lateral conditions
      CASE ( -1, 0, 1 )                ! all exept 2 (i.e. close case)
         iihom = nlci-nreci-jpr2di
         DO jl = 1, ipreci
            pt2d_ad(iihom +jl,:) = pt2d_ad(iihom +jl,:) + tr2we_ad(:,jl,1)
            tr2we_ad(:,jl,1)     = 0.0_wp
            pt2d_ad(jpreci+jl,:) = pt2d_ad(jpreci+jl,:) + tr2ew_ad(:,jl,1)
            tr2ew_ad(:,jl,1)     = 0.0_wp
         END DO
      END SELECT
      !
      ! 1. standard boundary treatment
      ! ------------------------------
      ! Order matters Here !!!!
      !
      ! north fold treatment
      ! -----------------------
      IF( npolj /= 0 ) THEN
         !
         SELECT CASE ( jpni )
         CASE ( 1 )     ;   CALL lbc_nfd_adj        ( pt2d_ad(1:jpi,1:jpj+jpr2dj), cd_type, psgn, pr2dj=jpr2dj )
         CASE DEFAULT   ;   CALL mpp_lbc_north_e_adj( pt2d_ad                    , cd_type, psgn               )
         END SELECT
         !
      ENDIF
      !                                      ! East-West boundaries
      !                                           !* Cyclic east-west
      IF( nbondi == 2 .AND. (nperio == 1 .OR. nperio == 4 .OR. nperio == 6) ) THEN
         pt2d_ad(     2      :2+jpr2di  ,:) = pt2d_ad(     2      :2+jpr2di  ,:) &
            &                               + pt2d_ad(    jpi     :jpi+jpr2di,:)! west
         pt2d_ad(   jpi      :jpi+jpr2di,:) = 0.0_wp
         pt2d_ad(jpim1-jpr2di:  jpim1   ,:) = pt2d_ad(jpim1-jpr2di:  jpim1   ,:) &
            &                               + pt2d_ad(1-jpr2di    :     1    ,:)! east
         pt2d_ad(1-jpr2di    :     1    ,:) = 0.0_wp
         !
      ELSE                                        !* closed
                                      pt2d_ad(nlci-jpreci+1:jpi+jpr2di,:) = 0.e0    ! north
         IF( .NOT. cd_type == 'F' )   pt2d_ad(  1-jpr2di   :jpreci    ,:) = 0.e0    ! south except at F-point
      ENDIF
      !
      !                                      !* North-South boundaries (always colsed)
                                   pt2d_ad(:,nlcj-jprecj+1:jpj+jpr2dj) = 0.e0    ! north
      IF( .NOT. cd_type == 'F' )   pt2d_ad(:,  1-jpr2dj   :  jprecj  ) = 0.e0    ! south except at F-point


   END SUBROUTINE mpp_lnk_2d_e_adj


   SUBROUTINE mpp_lbc_north_3d_adj( pt3d_ad, cd_type, psgn )
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_3d_adj  ***
      !!
      !! ** Purpose :  Adjoint of Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4 northern lines of the global domain on 1 processor
      !!              and apply lbc north-fold on this sub array. Then we
      !!              scatter the north fold array back to the processors.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj,jpk), INTENT(inout) ::   pt3d_ad      ! 3D array on which the b.c. is applied
      CHARACTER(len=1)                , INTENT(in   ) ::   cd_type   ! nature of pt3d_ad grid-points
      !                                                              !   = T ,  U , V , F or W  gridpoints
      REAL(wp)                        , INTENT(in   ) ::   psgn      ! = -1. the sign change across the north fold
      !!                                                             ! =  1. , the sign is kept
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille, ildi, ilei, iilb
      INTEGER ::   ijpj, ijpjm1, ij, iproc
      INTEGER, DIMENSION (jpmaxngh)          ::   ml_req_nf          ! for mpi_isend when avoiding mpi_allgather
      INTEGER                                ::   ml_err             ! for mpi_isend when avoiding mpi_allgather
      INTEGER, DIMENSION(MPI_STATUS_SIZE)    ::   ml_stat            ! for mpi_isend when avoiding mpi_allgather
      !!----------------------------------------------------------------------
      !
      ijpj   = 4
      ityp = -1
      ijpjm1 = 3
      ztabad(:,:,:) = 0.0_wp  ;  znorthlocad (:,:,:)= 0.0_wp  ;  znorthgloioad(:,:,:,:) = 0.0_wp
      !
      !
      DO jj = nlcj, nlcj-ijpj+1, -1          ! Scatter back to pt3d_ad
         ij = jj - nlcj + ijpj
         DO ji= nlci, 1, -1
            ztabad(ji+nimpp-1,ij,:) = ztabad(ji+nimpp-1,ij,:) + pt3d_ad(ji,jj,:)
            pt3d_ad(ji,jj,:) = 0.0_wp
         END DO
      END DO
      !
      ! The ztabad array has been either:
      !  a. Fully populated by the mpi_allgather operation or
      !  b. Had the active points for this domain and northern neighbours populated
      !     by peer to peer exchanges
      ! Either way the array may be folded by lbc_nfd and the result for the span of
      ! this domain will be identical.
      !
      CALL lbc_nfd_adj( ztabad, cd_type, psgn )   ! North fold boundary condition
      !
      !                                     ! Build in procs of ncomm_north the znorthgloioad
      itaille = jpi * jpk * ijpj
      IF ( l_north_nogather ) THEN
         !
         ! Set the exchange type in order to access the correct list of active neighbours
         !
         SELECT CASE ( cd_type )
            CASE ( 'T' , 'W' )
               ityp = 1
            CASE ( 'U' )
               ityp = 2
            CASE ( 'V' )
               ityp = 3
            CASE ( 'F' )
               ityp = 4
            CASE ( 'I' )
               ityp = 5
            CASE DEFAULT
               ityp = -1                    ! Set a default value for unsupported types which
                                            ! will cause a fallback to the mpi_allgather method
         END SELECT
         IF ( ityp .gt. 0 ) THEN
         !
            DO jr = 1, nsndto(ityp)
               iproc = isendto(jr,ityp) + 1
               ildi = nldit (iproc)
               ilei = nleit (iproc)
               iilb = nimppt(iproc)
               DO jj = ijpj, 1, -1
                  DO ji = ilei, ildi, -1
                     zfoldwkad(ji,jj,:) = ztabad(ji+iilb-1,jj,:)
                     ztabad(ji+iilb-1,jj,:) = 0.0_wp
                  END DO
               END DO
               CALL mppsend(5, zfoldwkad, itaille, isendto(jr,ityp), ml_req_nf(jr) )
            END DO
            DO jr = 1, nsndto(ityp)
               CALL mpprecv(5, znorthlocad, itaille, isendto(jr,ityp))
            END DO
            IF (l_isend) THEN
               DO jr = 1, nsndto(ityp)
                  CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
               END DO
            ENDIF
            !
         ENDIF
         !
         ! Avoid the use of mpi_allgather by exchanging only with the processes already identified
         ! (in nemo_northcomms) as being  involved in this process' northern boundary exchange
         !
         DO jj = nlcj, nlcj-ijpj+1, -1         ! First put local values into the global array
            ij = jj - nlcj + ijpj
            DO ji = nlci, 1, -1
               pt3d_ad(ji,jj,:) = pt3d_ad(ji,jj,:) + ztabad(ji+nimpp-1,ij,:)
               ztabad(ji+nimpp-1,ij,:) = 0.0_wp
            END DO
         END DO
         !
      ENDIF
      IF ( ityp .lt. 0 ) THEN
         DO jr = 1, ndim_rank_north         ! recover the global north array
            iproc = nrank_north(jr) + 1
            ildi  = nldit (iproc)
            ilei  = nleit (iproc)
            iilb  = nimppt(iproc)
            DO jj = 1, ijpj
               DO ji = ildi, ilei
                  znorthgloioad(ji,jj,:,jr) = znorthgloioad(ji,jj,:,jr) + ztabad(ji+iilb-1,jj,:)
                  ztabad(ji+iilb-1,jj,:) = 0.0_wp
               END DO
            END DO
         END DO
         !                                     ! Build in procs of ncomm_north the znorthgloio
         itaille = jpi * jpk * ijpj
         ! Specific treatment of adjoint of mpi_allgather
         znorthgloioad = mpp_sum_nfd(znorthgloioad,jpi,4,jpk,jpni,ncomm_north)
         jr=  ndim_rank_north-jpnij+nproc+1
         znorthlocad(:,:,:) = znorthgloioad(:,:,:,jr)
         !
      ENDIF
      !
      DO jj = nlcj, nlcj - ijpj +1, -1          ! put in znorthlocad the last 4 jlines of pt3d_ad
         ij = jj - nlcj + ijpj
         pt3d_ad(:,jj,:) = pt3d_ad(:,jj,:) + znorthlocad(:,ij,:)
         znorthlocad(:,ij,:) = 0.0_wp
      END DO
      !
   END SUBROUTINE mpp_lbc_north_3d_adj


   SUBROUTINE mpp_lbc_north_2d_adj( pt2d_ad, cd_type, psgn)
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_2d_adj  ***
      !!
      !! ** Purpose :   Adjoint of Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 (for 2d array )
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4 northern lines of the global domain on 1 processor
      !!              and apply lbc north-fold on this sub array. Then we
      !!              scatter the north fold array back to the processors.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(jpi,jpj), INTENT(inout) ::   pt2d_ad      ! 3D array on which the b.c. is applied
      CHARACTER(len=1)            , INTENT(in   ) ::   cd_type   ! nature of pt3d_ad grid-points
      !                                                          !   = T ,  U , V , F or W  gridpoints
      REAL(wp)                    , INTENT(in   ) ::   psgn      ! = -1. the sign change across the north fold
      !!                                                             ! =  1. , the sign is kept
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille, ildi, ilei, iilb
      INTEGER ::   ijpj, ijpjm1, ij, iproc
      INTEGER, DIMENSION (jpmaxngh)      ::   ml_req_nf          ! for mpi_isend when avoiding mpi_allgather
      INTEGER                            ::   ml_err             ! for mpi_isend when avoiding mpi_allgather
      INTEGER, DIMENSION(MPI_STATUS_SIZE)::   ml_stat            ! for mpi_isend when avoiding mpi_allgather
      !!----------------------------------------------------------------------
      !
      ijpj   = 4
      ityp = -1
      ijpjm1 = 3
      ztabad_2d(:,:) = 0.0_wp ; znorthlocad_2d = 0.0_wp ; znorthgloioad_2d = 0.0_wp
      !
      DO jj = nlcj, nlcj-ijpj+1, -1             ! Scatter back to pt2d_ad
         ij = jj - nlcj + ijpj
         DO ji = nlci, 1, -1
            ztabad_2d(ji+nimpp-1,ij) = ztabad_2d(ji+nimpp-1,ij) + pt2d_ad(ji,jj)
            pt2d_ad(ji,jj) = 0.0_wp
         END DO
      END DO
      !
      ! The ztabad array has been either:
      !  a. Fully populated by the mpi_allgather operation or
      !  b. Had the active points for this domain and northern neighbours populated
      !     by peer to peer exchanges
      ! Either way the array may be folded by lbc_nfd and the result for the span of
      ! this domain will be identical.
      !
      CALL lbc_nfd_adj( ztabad_2d, cd_type, psgn )   ! North fold boundary condition
      !
      !
      !                                     ! Build in procs of ncomm_north the znorthgloioad_2d
      itaille = jpi * ijpj
      IF ( l_north_nogather ) THEN
      !
      ! Set the exchange type in order to access the correct list of active neighbours
      !
         SELECT CASE ( cd_type )
            CASE ( 'T' , 'W' )
               ityp = 1
            CASE ( 'U' )
               ityp = 2
            CASE ( 'V' )
               ityp = 3
            CASE ( 'F' )
               ityp = 4
            CASE ( 'I' )
               ityp = 5
            CASE DEFAULT
               ityp = -1                    ! Set a default value for unsupported types which
         !                                  ! will cause a fallback to the mpi_allgather method
         END SELECT
         !
         IF ( ityp .gt. 0 ) THEN
            !
            DO jr = nsndto(ityp), 1, -1
               iproc = isendto(jr,ityp) + 1
               ildi = nldit (iproc)
               ilei = nleit (iproc)
               iilb = nimppt(iproc)
               DO jj = ijpj, 1, -1
                  DO ji = ilei, ildi, -1
                     zfoldwkad_2d(ji,jj) = ztabad_2d(ji+iilb-1,jj)
                     ztabad_2d(ji+iilb-1,jj) = 0.0_wp
                  END DO
               END DO
               CALL mppsend(5, zfoldwkad_2d, itaille, isendto(jr,ityp), ml_req_nf(jr))
            END DO
            DO jr = nsndto(ityp), 1, -1
               CALL mpprecv(5, znorthlocad_2d, itaille, isendto(jr,ityp) )
            END DO
            !
            IF (l_isend) THEN
               DO jr = nsndto(ityp), 1, -1
                  CALL mpi_wait(ml_req_nf(jr), ml_stat, ml_err)
               END DO
            ENDIF
            !
         ENDIF
         !
         ! Avoid the use of mpi_allgather by exchanging only with the processes already identified
         ! (in nemo_northcomms) as being  involved in this process' northern boundary exchange
         !
         DO jj = nlcj-ijpj+1, nlcj          ! First put local values into the global array
            ij = jj - nlcj + ijpj
            DO ji = 1, nlci
               pt2d_ad(ji,jj) = pt2d_ad(ji,jj) + ztabad_2d(ji+nimpp-1,ij)
               ztabad_2d(ji+nimpp-1,ij) = 0.0_wp
            END DO
         END DO
         !
      ENDIF
      !
      IF ( ityp .lt. 0 ) THEN
         DO jr = 1, ndim_rank_north            ! recover the global north array
            iproc = nrank_north(jr) + 1
            ildi = nldit (iproc)
            ilei = nleit (iproc)
            iilb = nimppt(iproc)
            DO jj = 1, ijpj
               DO ji = ildi, ilei
                  znorthgloioad_2d(ji,jj,jr) = znorthgloioad_2d(ji,jj,jr) + ztabad_2d(ji+iilb-1,jj)
                  ztabad_2d(ji+iilb-1,jj) = 0.0_wp
               END DO
            END DO
         END DO
         !
         znorthgloioad_2d = mpp_sum_nfd(znorthgloioad_2d,jpi,4,jpni,ncomm_north)
         jr=  ndim_rank_north-jpnij+nproc+1
         znorthlocad_2d(:,:) = znorthgloioad_2d(:,:,jr)
         !
      ENDIF
      DO jj = nlcj, nlcj - ijpj +1, -1          ! put in znorthlocad the last 4 jlines of pt2d_ad
         ij = jj - nlcj + ijpj
         pt2d_ad(:,jj) = pt2d_ad(:,jj) + znorthlocad_2d(:,ij)
         znorthlocad_2d(:,ij) = 0.0_wp
      END DO
      !
   END SUBROUTINE mpp_lbc_north_2d_adj


   SUBROUTINE mpp_lbc_north_e_adj( pt2d_ad, cd_type, psgn)
      !!---------------------------------------------------------------------
      !!                   ***  routine mpp_lbc_north_e_adj  ***
      !!
      !! ** Purpose : Adjoint of Ensure proper north fold horizontal bondary condition
      !!              in mpp configuration in case of jpn1 > 1 and for 2d
      !!              array with outer extra halo
      !!
      !! ** Method  :   North fold condition and mpp with more than one proc
      !!              in i-direction require a specific treatment. We gather
      !!              the 4+2*jpr2dj northern lines of the global domain on 1
      !!              processor and apply lbc north-fold on this sub array.
      !!              Then we scatter the north fold array back to the processors.
      !!
      !!----------------------------------------------------------------------
      REAL(wp), DIMENSION(1-jpr2di:jpi+jpr2di,1-jpr2dj:jpj+jpr2dj), INTENT(inout) ::   pt2d_ad     ! 2D array with extra halo
      CHARACTER(len=1)                                            , INTENT(in   ) ::   cd_type  ! nature of pt3d_ad grid-points
      !                                                                                         !   = T ,  U , V , F or W -points
      REAL(wp)                                                    , INTENT(in   ) ::   psgn     ! = -1. the sign change across the
      !!                                                                                        ! north fold, =  1. otherwise
      INTEGER ::   ji, jj, jr
      INTEGER ::   ierr, itaille, ildi, ilei, iilb
      INTEGER ::   ijpj, ij, iproc
      !!----------------------------------------------------------------------
      !
      ijpj=4
      ztabad_e = 0.e0  ;  znorthlocad_e = 0.0_wp  ;  znorthgloioad_e = 0.0_wp
      !
      ij = jpr2dj
      !! Scatter back to pt2d
      DO jj = nlcj - ijpj + 1 , nlcj +jpr2dj
      ij  = ij +1
         DO ji= 1, nlci
            ztabad_e(ji+nimpp-1,ij) = ztabad_e(ji+nimpp-1,ij) + pt2d_ad(ji,jj)
            pt2d_ad(ji,jj)        = 0.0_wp
         END DO
      END DO
      !
      ! 2. North-Fold boundary conditions
      ! ----------------------------------
      CALL lbc_nfd_adj( ztabad_e(:,:), cd_type, psgn, jpr2dj)!, pr2dj = jpr2dj )
      !
      DO jr = 1, ndim_rank_north            ! recover the global north array
         iproc = nrank_north(jr) + 1
         ildi = nldit (iproc)
         ilei = nleit (iproc)
         iilb = nimppt(iproc)
         DO jj = 1, ijpj+2*jpr2dj
            DO ji = ildi, ilei
               znorthgloioad_e(ji,jj,jr) = znorthgloioad_e(ji,jj,jr) + ztabad_e(ji+iilb-1,jj)
               ztabad_e(ji+iilb-1,jj) = 0.0_wp
            END DO
         END DO
      END DO
      !
      itaille = jpi * ( ijpj + 2 * jpr2dj )
      ! Specific treatment of adjoint of mpi_allgather
      znorthgloioad_e = mpp_sum_nfd(znorthgloioad_e,jpi,4,jpni,ncomm_north)
      jr=  ndim_rank_north-jpnij+nproc+1
      znorthlocad_e(:,:) = znorthgloioad_e(:,:,jr)
      !
      ij=0
      ! put in znorthloc the last 4 jlines of pt2d
      DO jj = nlcj - ijpj + 1 - jpr2dj, nlcj +jpr2dj
         ij = ij + 1
         DO ji = 1, jpi
            pt2d_ad(ji,jj)     = pt2d_ad(ji,jj) + znorthlocad_e(ji,ij)
            znorthlocad_e(ji,ij) = 0.0_wp
         END DO
      END DO
      !
   END SUBROUTINE mpp_lbc_north_e_adj
   FUNCTION mpp_sum_nfd_4d( pval, kn1, kn2, kn3, kn4, kcom )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_sum_nfd_4d ***
      !!
      !! ** Purpose : Summation of arrays across processors (of kcom group)
      !!
      !! ** Method  : Pack 4d array to vector, do vector sum and then unpack
      !!
      !! ** Action  : This does only work for MPI.
      !!              It does not work for SHMEM.
      !!
      !! History :
      !!        !  10-11  (F. Vigilant)  Original code
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp), DIMENSION(kn1, kn2, kn3, kn4) :: mpp_sum_nfd_4d
      !! * Arguments
      INTEGER, INTENT(IN) :: kn1, kn2, kn3, kn4
      REAL(wp), DIMENSION(kn1, kn2, kn3, kn4), INTENT(IN) :: pval
      INTEGER , INTENT( in ), OPTIONAL        :: kcom
      !! * Local declarations
      REAL(wp), DIMENSION(:), ALLOCATABLE     :: zvec
      LOGICAL, DIMENSION(kn1, kn2, kn3, kn4)  :: zmask
      REAL(wp), DIMENSION(kn1, kn2, kn3, kn4) :: zfd
      INTEGER :: zdim

      zdim = kn1 * kn2 * kn3 * kn4
      ALLOCATE( zvec(zdim) )

      zvec = PACK( pval(1:kn1,1:kn2,1:kn3,1:kn4),.TRUE.)

      CALL mpp_sum( zvec, zdim, kcom )

      zmask(:,:,:,:) = .TRUE.
      mpp_sum_nfd_4d = UNPACK( zvec, zmask, zfd )

      DEALLOCATE( zvec )

   END FUNCTION mpp_sum_nfd_4d

   FUNCTION mpp_sum_nfd_3d( pval, kn1, kn2, kn3, kcom )
      !!----------------------------------------------------------------------
      !!               ***  ROUTINE mpp_sum_nfd_3d ***
      !!
      !! ** Purpose : Summation of arrays across processors (of kcom group)
      !!
      !! ** Method  : Pack 3d array to vector, do vector sum and then unpack
      !!
      !! ** Action  : This does only work for MPI.
      !!              It does not work for SHMEM.
      !!
      !! History :
      !!        !  10-11  (F. Vigilant)  Original code
      !!----------------------------------------------------------------------
      !! * Function return
      REAL(wp), DIMENSION(kn1, kn2, kn3) :: mpp_sum_nfd_3d
      !! * Arguments
      INTEGER, INTENT(IN) ::kn1, kn2, kn3
      REAL(wp), DIMENSION(kn1, kn2, kn3), INTENT(IN) :: pval
      INTEGER , INTENT( in ), OPTIONAL        :: kcom
      !! * Local declarations
      REAL(wp), DIMENSION(:), ALLOCATABLE     :: zvec
      LOGICAL, DIMENSION(kn1, kn2, kn3)  :: zmask
      REAL(wp), DIMENSION(kn1, kn2, kn3) :: zfd
      INTEGER :: zdim

      zdim = kn1 * kn2 * kn3
      ALLOCATE( zvec(zdim) )

      zvec = PACK( pval(1:kn1,1:kn2,1:kn3),.TRUE.)

      CALL mpp_sum( zvec, zdim, kcom )

      zmask(:,:,:) = .TRUE.
      mpp_sum_nfd_3d = UNPACK( zvec, zmask, zfd )

      DEALLOCATE( zvec )

   END FUNCTION mpp_sum_nfd_3d
   !!----------------------------------------------------------------------
END MODULE lib_mpp_tam
