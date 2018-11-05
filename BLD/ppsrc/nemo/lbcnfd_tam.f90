MODULE lbcnfd_tam
   !!======================================================================
   !!                       ***  MODULE  lbcnfd_tam  ***
   !! Ocean        : TAM of north fold  boundary conditions
   !!======================================================================
   !! History :  3.2  ! 2009-03  (R. Benshila)  Original code
   !! History of TAM : 3.2 ! 2010-04 (F. Vigilant) Original Code
   !!                  3.4 ! 2012-03 (P.-A. Bouttier) Update
   !!----------------------------------------------------------------------

   !!----------------------------------------------------------------------
   !!   lbc_nfd_adj       : generic interface for lbc_nfd_3d_adj and lbc_nfd_2d_adj routines
   !!   lbc_nfd_3d_adj    : lateral boundary condition: North fold treatment for a 3D arrays   (lbc_nfd_adj)
   !!   lbc_nfd_2d_adj    : lateral boundary condition: North fold treatment for a 2D arrays   (lbc_nfd_adj)
   !!----------------------------------------------------------------------
   USE dom_oce        ! ocean space and time domain
   USE in_out_manager ! I/O manager
   USE lib_mpp

   IMPLICIT NONE
   PRIVATE

   INTERFACE lbc_nfd_adj
      MODULE PROCEDURE   lbc_nfd_3d_adj, lbc_nfd_2d_adj
   END INTERFACE

   PUBLIC   lbc_nfd_adj      ! north fold conditions
   PUBLIC   lbc_nfd_adj_tst  ! Adjoint test routine
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE lbc_nfd_3d_adj( pt3d, cd_type, psgn )
      !!----------------------------------------------------------------------
      !!                  ***  routine lbc_nfd_3d_adj  ***
      !!
      !! ** Purpose :   Adjoint of 3D lateral boundary condition : North fold treatment
      !!              without processor exchanges.
      !!
      !! ** Method  :
      !!
      !! ** Action  :   pt3d with updated values along the north fold
      !!----------------------------------------------------------------------
      CHARACTER(len=1)          , INTENT(in   ) ::   cd_type   ! define the nature of ptab array grid-points
      !                                                        !   = T , U , V , F , W points
      REAL(wp)                  , INTENT(in   ) ::   psgn      ! control of the sign change
      !                                                        !   = -1. , the sign is changed if north fold boundary
      !                                                        !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(:,:,:), INTENT(inout) ::   pt3d      ! 3D array on which the boundary condition is applied
      !
      INTEGER  ::   ji, jk
      INTEGER  ::   zijt, ziju, zijpj, zijpjm1
      ! 
      REAL(wp) :: ztmp
      !!----------------------------------------------------------------------

      SELECT CASE ( jpni )
      CASE ( 1 )     ;   zijpj = nlcj      ! 1 proc only  along the i-direction
      CASE DEFAULT   ;   zijpj = 4         ! several proc along the i-direction
      END SELECT
      zijpjm1 = zijpj-1

      DO jk = 1, jpk
         !
         SELECT CASE ( npolj )
         !
         CASE ( 3 , 4 )                        ! *  North fold  T-point pivot
            !
            SELECT CASE ( cd_type )
         CASE ( 'T' , 'W' )                         ! T-, W-point
            DO ji = jpiglo, jpiglo/2+1, -1
               zijt = jpiglo-ji+2
               ztmp = psgn * pt3d(ji,zijpjm1,jk)
               pt3d(ji ,zijpjm1,jk) = 0.0_wp
               pt3d(zijt,zijpjm1,jk) = pt3d(zijt,zijpjm1,jk) + ztmp
            END DO

            DO ji = jpiglo, 2, -1
               zijt = jpiglo-ji+2
               pt3d(zijt,zijpj-2,jk) = pt3d(zijt,zijpj-2,jk) + psgn * pt3d(ji,zijpj,jk)
               pt3d(ji ,zijpj  ,jk) = 0.0_wp
            END DO
         CASE ( 'U' )                               ! U-point
            DO ji = jpiglo-1, jpiglo/2, -1
               ziju = jpiglo-ji+1
               ztmp = psgn * pt3d(ji,zijpjm1,jk)
               pt3d(ji ,zijpjm1,jk) = 0.0_wp
               pt3d(ziju,zijpjm1,jk) = pt3d(ziju,zijpjm1,jk) + ztmp
            END DO

            DO ji = jpiglo-1, 1, -1
               ziju = jpiglo-ji+1
               pt3d(ziju,zijpj-2,jk) = pt3d(ziju,zijpj-2,jk) + psgn * pt3d(ji,zijpj,jk)
               pt3d(ji ,zijpj  ,jk) = 0.0_wp
            END DO
         CASE ( 'V' )                               ! V-point
            DO ji = jpiglo, 2, -1
               zijt = jpiglo-ji+2
               pt3d(zijt,zijpj-3,jk) = pt3d(zijt,zijpj-3,jk) + psgn * pt3d(ji,zijpj  ,jk)
               pt3d(ji,zijpj  ,jk) = 0.0_wp
               pt3d(zijt,zijpj-2,jk) = pt3d(zijt,zijpj-2,jk) + psgn * pt3d(ji,zijpj-1,jk)
               pt3d(ji,zijpj-1,jk) = 0.0_wp
            END DO
         CASE ( 'F' )                               ! F-point
            DO ji = jpiglo-1, 1, -1
               ziju = jpiglo-ji+1
               pt3d(ziju,zijpj-3,jk) = pt3d(ziju,zijpj-3,jk) + psgn * pt3d(ji,zijpj  ,jk)
               pt3d(ji ,zijpj  ,jk) = 0.0_wp
               pt3d(ziju,zijpj-2,jk) = pt3d(ziju,zijpj-2,jk) + psgn * pt3d(ji,zijpj-1,jk)
               pt3d(ji ,zijpj-1,jk) = 0.0_wp
            END DO
         END SELECT
            !
         CASE ( 5 , 6 )                        ! *  North fold  F-point pivot
            !
            SELECT CASE ( cd_type )
            CASE ( 'T' , 'W' )                         ! T-, W-point
            DO ji = jpiglo, 1, -1
               zijt = jpiglo-ji+1
               pt3d(zijt,zijpj-1,jk) = pt3d(zijt,zijpj-1,jk) + psgn * pt3d(ji,zijpj,jk)
               pt3d(ji ,zijpj  ,jk) = 0.0_wp
            END DO
         CASE ( 'U' )                               ! U-point
            DO ji = jpiglo-1, 1, -1
               ziju = jpiglo-ji
               pt3d(ziju,zijpj-1,jk) = pt3d(ziju,zijpj-1,jk) + psgn * pt3d(ji,zijpj,jk)
               pt3d(ji ,zijpj  ,jk) = 0.0_wp
            END DO
         CASE ( 'V' )                               ! V-point
            DO ji = jpiglo, jpiglo/2+1, -1
               zijt = jpiglo-ji+1
               ztmp = psgn * pt3d(ji,zijpjm1,jk)
               pt3d(ji ,zijpjm1,jk) = 0.0_wp
               pt3d(zijt,zijpjm1,jk) = pt3d(zijt,zijpjm1,jk) + ztmp
            END DO
            DO ji = jpiglo, 1, -1
               zijt = jpiglo-ji+1
               pt3d(zijt,zijpj-2,jk) = pt3d(zijt,zijpj-2,jk) + psgn * pt3d(ji,zijpj,jk)
               pt3d(ji ,zijpj  ,jk) = 0.0_wp
            END DO
         CASE ( 'F' )                               ! F-point
            DO ji = jpiglo-1, jpiglo/2+1, -1
               ziju = jpiglo-ji
               ztmp = psgn * pt3d(ji,zijpjm1,jk)
               pt3d(ji ,zijpjm1,jk) = 0.0_wp
               pt3d(ziju,zijpjm1,jk) = pt3d(ziju,zijpjm1,jk) + ztmp
            END DO
            DO ji = jpiglo-1, 1, -1
               ziju = jpiglo-ji
               pt3d(ziju,zijpj-2,jk) = pt3d(ziju,zijpj-2,jk) + psgn * pt3d(ji,zijpj  ,jk)
               pt3d(ji ,zijpj  ,jk) = 0.0_wp
            END DO
         END SELECT
            !
         CASE DEFAULT                           ! *  closed : the code probably never go through
            !
            SELECT CASE ( cd_type)
            CASE ( 'T' , 'U' , 'V' , 'W' )             ! T-, U-, V-, W-points
               pt3d(:, 1  ,jk) = 0.e0
               pt3d(:,zijpj,jk) = 0.e0
            CASE ( 'F' )                               ! F-point
               pt3d(:,zijpj,jk) = 0.e0
            END SELECT
            !
         END SELECT     !  npolj
         !
      END DO
      !
   END SUBROUTINE lbc_nfd_3d_adj


   SUBROUTINE lbc_nfd_2d_adj( pt2d, cd_type, psgn, pr2dj )
      !!----------------------------------------------------------------------
      !!                  ***  routine lbc_nfd_2d_adj  ***
      !!
      !! ** Purpose :   Adjoint of 2D lateral boundary condition : North fold treatment
      !!       without processor exchanges.
      !!
      !! ** Method  :
      !!
      !! ** Action  :   pt2d with updated values along the north fold
      !!----------------------------------------------------------------------
      CHARACTER(len=1)        , INTENT(in   ) ::   cd_type   ! define the nature of ptab array grid-points
      !                                                      ! = T , U , V , F , W points
      REAL(wp)                , INTENT(in   ) ::   psgn      ! control of the sign change
      !                                                      !   = -1. , the sign is changed if north fold boundary
      !                                                      !   =  1. , the sign is kept  if north fold boundary
      REAL(wp), DIMENSION(:,:), INTENT(inout) ::   pt2d      ! 2D array on which the boundary condition is applied
      INTEGER , OPTIONAL      , INTENT(in   ) ::   pr2dj     ! number of additional halos
      !
      INTEGER  ::   ji, jl, ipr2dj
      INTEGER  ::   zijt, ziju, zijpj, zijpjm1
      !
      REAL(wp) :: ztmp
      !!----------------------------------------------------------------------

      SELECT CASE ( jpni )
      CASE ( 1 )     ;   zijpj = nlcj      ! 1 proc only  along the i-direction
      CASE DEFAULT   ;   zijpj = 4         ! several proc along the i-direction
      END SELECT
      !
      IF( PRESENT(pr2dj) ) THEN           ! use of additional halos
         ipr2dj = pr2dj
         IF( jpni > 1 )   zijpj = zijpj + ipr2dj
      ELSE
         ipr2dj = 0
      ENDIF
      !
      zijpjm1 = zijpj-1
      SELECT CASE ( npolj )
      !
      CASE ( 3, 4 )                       ! *  North fold  T-point pivot
         !
         SELECT CASE ( cd_type )
         !
         CASE ( 'T', 'W' )
            DO ji = jpiglo, jpiglo/2+1, -1
               zijt=jpiglo-ji+2
               ztmp = psgn * pt2d(ji,zijpj-1)
               pt2d(ji,zijpj-1) = 0.0_wp
               pt2d(zijt,zijpj-1) = pt2d(zijt,zijpj-1) + ztmp
            END DO
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo, 2, -1
                  zijt=jpiglo-ji+2
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji ,zijpj+jl  ) = 0.0_wp
                  pt2d(zijt,zijpj-2-jl) = pt2d(zijt,zijpj-2-jl) + ztmp
               END DO
            END DO
         CASE ( 'U' )                                     ! U-point
            DO ji = jpiglo-1, jpiglo/2, -1
               ziju = jpiglo-ji+1
               ztmp = psgn * pt2d(ji,zijpjm1)
               pt2d(ji,zijpjm1) = 0.0_wp
               pt2d(ziju,zijpjm1) = pt2d(ziju,zijpjm1) + ztmp
            END DO
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo-1, 1, -1
                  ziju = jpiglo-ji+1
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji ,zijpj+jl  ) = 0.0_wp
                  pt2d(ziju,zijpj-2-jl) = pt2d(ziju,zijpj-2-jl) + ztmp
               END DO
            END DO
         CASE ( 'V' )                                     ! V-point
            DO jl = ipr2dj, -1, -1
               DO ji = jpiglo, 2, -1
                  zijt = jpiglo-ji+2
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji ,zijpj+jl  ) = 0.0_wp
                  pt2d(zijt,zijpj-3-jl) = pt2d(zijt,zijpj-3-jl) + ztmp
               END DO
            END DO
         CASE ( 'F' )                               ! F-point
            DO jl = ipr2dj, -1, -1
               DO ji = jpiglo-1, 1, -1
                  ziju = jpiglo-ji+1
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
                  pt2d(ziju,zijpj-3-jl) = pt2d(ziju,zijpj-3-jl) + ztmp
               END DO
            END DO
         CASE ( 'I' )                                     ! ice U-V point
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo, 3, -1
                  ziju = jpiglo - ji + 3
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
                  pt2d(ziju,zijpj-1-jl) = pt2d(ziju,zijpj-1-jl) + ztmp
               END DO
               ztmp = psgn * pt2d(2,zijpj+jl)
               pt2d(2,zijpj+jl) = 0.0_wp
               pt2d(3,zijpj-1+jl) = pt2d(3,zijpj-1+jl) + ztmp
            END DO
         CASE ( 'J' )                                     ! first ice U-V point
            DO jl =ipr2dj,0,-1
               DO ji = 3, jpiglo
                  ziju = jpiglo - ji + 3
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
                  pt2d(ziju,zijpj-1-jl) = pt2d(ziju,zijpj-1-jl) + ztmp
               END DO
               pt2d(3,zijpj-1+jl) = pt2d(3,zijpj-1+jl) + psgn * pt2d(2,zijpj+jl)
               pt2d(2,zijpj+jl) = 0.0_wp
            END DO
         CASE ( 'K' )                                     ! second ice U-V point
            DO jl =0, ipr2dj
               DO ji = 3, jpiglo
                  ziju = jpiglo - ji + 3
                  pt2d(3,zijpj-1+jl) = pt2d(3,zijpj-1+jl) + psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
               END DO
               ztmp = psgn * pt2d(2,zijpj+jl)
               pt2d(2,zijpj+jl) = 0.0_wp
               pt2d(3,zijpj-1+jl) = pt2d(3,zijpj-1+jl) + ztmp
            END DO
         END SELECT
         !
      CASE ( 5, 6 )                        ! *  North fold  F-point pivot
         !
         SELECT CASE ( cd_type )
         CASE ( 'T' , 'W' )                          ! T-, W-point
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo, 1, -1
                  zijt = jpiglo-ji+1
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji ,zijpj+jl  ) = 0.0_wp
                  pt2d(zijt,zijpj-1-jl) = pt2d(zijt,zijpj-1-jl) + ztmp
               END DO
            END DO
         CASE ( 'U' )                                     ! U-point
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo-1, 1, -1
                  ziju = jpiglo-ji
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
                  pt2d(ziju,zijpj-1-jl) = pt2d(ziju,zijpj-1-jl) + ztmp
               END DO
            END DO
         CASE ( 'V' )                                     ! V-point
            DO ji = jpiglo, jpiglo/2+1, -1
               zijt = jpiglo-ji+1
               ztmp = psgn * pt2d(ji,zijpjm1)
               pt2d(ji ,zijpjm1) = 0.0_wp
               pt2d(zijt,zijpjm1) = pt2d(zijt,zijpjm1) + ztmp
            END DO
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo, 1, -1
                  zijt = jpiglo-ji+1
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
                  pt2d(zijt,zijpj-2-jl) = pt2d(zijt,zijpj-2-jl) + ztmp
               END DO
            END DO
         CASE ( 'F' )                               ! F-point
            DO ji = jpiglo-1, jpiglo/2+1, -1
               ziju = jpiglo-ji
               ztmp = psgn * pt2d(ji,zijpjm1)
               pt2d(ji ,zijpjm1) = 0.0_wp
               pt2d(ziju,zijpjm1) = pt2d(ziju,zijpjm1) + ztmp
            END DO
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo-1, 1, -1
                  ziju = jpiglo-ji
                  ztmp = psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji ,zijpj+jl  ) = 0.0_wp
                  pt2d(ziju,zijpj-2-jl) = pt2d(ziju,zijpj-2-jl) + ztmp
               END DO
            END DO
         CASE ( 'I' )                                  ! ice U-V point
            DO jl = ipr2dj, 0, -1
               DO ji = jpiglo-1, 2, -1
                  zijt = jpiglo - ji + 2
                  pt2d(ji ,zijpj-1-jl) = pt2d(ji ,zijpj-1-jl) + 0.5 * pt2d(ji,zijpj+jl)
                  pt2d(zijt,zijpj-1-jl) = pt2d(zijt,zijpj-1-jl) + 0.5 * psgn * pt2d(ji,zijpj+jl)
                  pt2d(ji ,zijpj+jl  ) = 0.0_wp
               END DO
            END DO
            pt2d( 2 ,zijpj:zijpj+ipr2dj) = 0.0_wp
         CASE ( 'J' )                                  ! first ice U-V point
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  zijt = jpiglo - ji + 2
                  pt2d(ji,zijpj-1-jl) = pt2d(ji,zijpj-1-jl) + pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
               END DO
            END DO
            pt2d( 2 ,zijpj:zijpj+ipr2dj) = 0.e0
         CASE ( 'K' )                                  ! second ice U-V point
            DO jl = 0, ipr2dj
               DO ji = 2 , jpiglo-1
                  zijt = jpiglo - ji + 2
                  pt2d(zijt,zijpj-1-jl) = pt2d(zijt,zijpj-1-jl) + pt2d(ji,zijpj+jl)
                  pt2d(ji,zijpj+jl) = 0.0_wp
               END DO
            END DO
            pt2d( 2 ,zijpj:zijpj+ipr2dj) = 0.e0
         END SELECT
         !
      CASE DEFAULT                           ! *  closed : the code probably never go through
         !
         SELECT CASE ( cd_type)
         CASE ( 'T' , 'U' , 'V' , 'W' )                 ! T-, U-, V-, W-points
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,zijpj:zijpj+ipr2dj) = 0.e0
         CASE ( 'F' )                                   ! F-point
            pt2d(:,zijpj:zijpj+ipr2dj) = 0.e0
         CASE ( 'I' )                                   ! ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,zijpj:zijpj+ipr2dj) = 0.e0
         CASE ( 'J' )                                   ! first ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,zijpj:zijpj+ipr2dj) = 0.e0
         CASE ( 'K' )                                   ! second ice U-V point
            pt2d(:, 1:1-ipr2dj     ) = 0.e0
            pt2d(:,zijpj:zijpj+ipr2dj) = 0.e0
         END SELECT
         !
      END SELECT
      !
   END SUBROUTINE lbc_nfd_2d_adj

   SUBROUTINE lbc_nfd_3d_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE lbc_nfd_3d_adj_tst ***
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
      !! ** Action  :
      !!
      !! History :
      !!        ! 2010-09 (F. Vigilant)
      !!-----------------------------------------------------------------------
      !! * Modules used
   USE gridrandom    , ONLY: & ! Random Gaussian noise on grids
      & grid_random
   USE dotprodfld    , ONLY: & ! Computes dot product for 3D and 2D fields
      & dot_product
   USE lbcnfd,         ONLY: &
      & lbc_nfd
   USE tstool_tam    , ONLY: &
      & stdu,                &
      & stdv,                &
      & stdt,                &
      & prntst_adj
   USE dom_oce       , ONLY: & ! Ocean space and time domain
      & tmask,               &
      & umask,               &
      & vmask,               &
      & mig,                 &
      & mjg,                 &
      & nldi,                &
      & nldj,                &
      & nlei,                &
      & nlej
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt        ! Output unit

      !! * Local declarations
      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zu_tlin,    & ! Tangent input: ua_tl
         & zv_tlin,    & ! Tangent input: va_tl
         & zt_tlin,    & ! Tangent input: ub_tl
         & zu_tlout,   & ! Tangent output: ua_tl
         & zv_tlout,   & ! Tangent output: va_tl
         & zt_tlout,   & ! Tangent output: ua_tl
         & zu_adin,    & ! Tangent output: ua_ad
         & zv_adin,    & ! Tangent output: va_ad
         & zt_adin,    & ! Adjoint input: ua_ad
         & z_adin,    & ! Adjoint input: va_ad
         & zu_adout,   & ! Adjoint output: ua_ad
         & zv_adout,   & ! Adjoint output: va_ad
         & zt_adout,   & ! Adjoint oputput: ub_ad
         & znu            ! 3D random field for u

      REAL(wp), DIMENSION(:,:,:), ALLOCATABLE :: &
         & zu_tl,    & ! Tangent input: ua_tl
         & zv_tl,    & ! Tangent input: va_tl
         & zt_tl,    & ! Tangent input: ub_tl
         & zu_ad,    & ! Tangent output: ua_ad
         & zv_ad,    & ! Tangent output: va_ad
         & zt_ad
      REAL(wp) :: &
         & zsp1,    &   ! scalar product involving the tangent routine
         & zsp2         ! scalar product involving the adjoint routine
      CHARACTER(len=1) :: &
         & cd_type
      INTEGER :: &
         & indic, &
         & istp
      INTEGER :: &
         & ji, &
         & jj, &
         & jk, &
         & kmod, &
         & jstp
      CHARACTER (LEN=14) :: &
         & cl_name
      INTEGER ::             &
         & zijpj, ziglo, zip, zdiff

      ! Allocate memory

      zijpj = 4

      SELECT CASE ( jpni )
      CASE ( 1 )     ;   zijpj = nlcj      ! 1 proc only  along the i-direction
      CASE DEFAULT   ;   zijpj = 4         ! several proc along the i-direction
      END SELECT

      ALLOCATE( &
         & zu_tlin(jpiglo,zijpj,jpk),  &
         & zv_tlin(jpiglo,zijpj,jpk),  &
         & zt_tlin(jpiglo,zijpj,jpk),  &
         & zu_tlout(jpiglo,zijpj,jpk), &
         & zv_tlout(jpiglo,zijpj,jpk), &
         & zt_tlout(jpiglo,zijpj,jpk), &
         & zu_adin(jpiglo,zijpj,jpk),  &
         & zv_adin(jpiglo,zijpj,jpk),  &
         & zt_adin(jpiglo,zijpj,jpk),  &
         & zu_adout(jpiglo,zijpj,jpk), &
         & zv_adout(jpiglo,zijpj,jpk), &
         & zt_adout(jpiglo,zijpj,jpk), &
         & znu(jpi,jpj,jpk)         &
         & )

      ALLOCATE( &
         & zu_tl(jpiglo,zijpj,jpk),  &
         & zv_tl(jpiglo,zijpj,jpk),  &
         & zt_tl(jpiglo,zijpj,jpk),  &
         & zu_ad(jpiglo,zijpj,jpk),  &
         & zv_ad(jpiglo,zijpj,jpk),  &
         & zt_ad(jpiglo,zijpj,jpk)   &
         & )

        zu_tlin (:,:,:) = 0.0_wp
        zv_tlin (:,:,:) = 0.0_wp
        zt_tlin (:,:,:) = 0.0_wp
        zu_tlout(:,:,:) = 0.0_wp
        zv_tlout(:,:,:) = 0.0_wp
        zt_tlout(:,:,:) = 0.0_wp
        zu_adin (:,:,:) = 0.0_wp
        zv_adin (:,:,:) = 0.0_wp
        zt_adin (:,:,:) = 0.0_wp
        zu_adout(:,:,:) = 0.0_wp
        zv_adout(:,:,:) = 0.0_wp
        zt_adout(:,:,:) = 0.0_wp

        zu_tl (:,:,:) = 0.0_wp
        zv_tl (:,:,:) = 0.0_wp
        zt_tl (:,:,:) = 0.0_wp
        zu_ad (:,:,:) = 0.0_wp
        zv_ad (:,:,:) = 0.0_wp
        zt_ad (:,:,:) = 0.0_wp

        ziglo = INT(jpiglo / (nlei - nldi) )
        zdiff = nlei - nldi
        !--------------------------------------------------------------------
        ! Initialize the tangent input with random noise: dx
        !--------------------------------------------------------------------
        CALL grid_random( znu, 'U', 0.0_wp, stdu )
        DO jk = 1, jpk
           DO jj = 1, 4
              DO ji = nldi, nlei
                 DO zip = 1, ziglo
                    zu_tlin(ji+(zip-1)*zdiff,jj,jk) = znu(ji,jj,jk)
                 ENDDO
              END DO
           END DO
        END DO
        CALL grid_random( znu, 'V', 0.0_wp, stdu )
        DO jk = 1, jpk
           DO jj = 1, 4
              DO ji = nldi, nlei
                 DO zip = 1, ziglo
                    zv_tlin(ji+(zip-1)*zdiff,jj,jk) = znu(ji,jj,jk)
                 ENDDO
              END DO
           END DO
        END DO
        CALL grid_random( znu, 'T', 0.0_wp, stdt )
        DO jk = 1, jpk
           DO jj = 1, 4
              DO ji = nldi, nlei
                 DO zip = 1, ziglo
                    zt_tlin(ji+(zip-1)*zdiff,jj,jk) = znu(ji,jj,jk)
                 ENDDO
              END DO
           END DO
        END DO

        !--------------------------------------------------------------------
        ! Call the tangent routine: dy = L dx
        !--------------------------------------------------------------------
        zu_tl(:,:,:) = zu_tlin(:,:,:)
        zv_tl(:,:,:) = zv_tlin(:,:,:)
        zt_tl(:,:,:) = zt_tlin(:,:,:)
        CALL lbc_nfd( zu_tl, 'U',  -1.0_wp )
        CALL lbc_nfd( zv_tl, 'V',  -1.0_wp )
        CALL lbc_nfd( zt_tl, 'T',   1.0_wp )
        zu_tlout(:,:,:) = zu_tl(:,:,:)
        zv_tlout(:,:,:) = zv_tl(:,:,:)
        zt_tlout(:,:,:) = zt_tl(:,:,:)

        DO jk = 1, jpk
           DO jj = 1,4
              DO ji = nldi, nlei
                 DO zip = 1, ziglo
                    zu_adin(ji+(zip-1)*zdiff,jj,jk) = zu_tlout(ji+(zip-1)*ziglo,jj,jk) &
                       &               * e1u(ji,jj) * e2u(ji,jj) &!* fse3u(ji,jj,jk) &
                       &               * umask(ji,jj,jk)
                    zv_adin(ji+(zip-1)*zdiff,jj,jk) = zv_tlout(ji+(zip-1)*ziglo,jj,jk) &
                       &               * e1v(ji,jj) * e2v(ji,jj) &!* fse3v(ji,jj,jk) &
                       &               * vmask(ji,jj,jk)
                    zt_adin(ji+(zip-1)*zdiff,jj,jk) = zt_tlout(ji+(zip-1)*ziglo,jj,jk) &
                       &               * e1t(ji,jj) * e2t(ji,jj) &!* fse3t(ji,jj,jk) &
                       &               * tmask(ji,jj,jk)
                 ENDDO
              END DO
           END DO
        END DO

        ! DOT_PRODUCT
        zsp1 = sum( PACK(zu_tlout(:,:,:),.TRUE.) * &
         &                       PACK( zu_adin(:,:,:),.TRUE.) )

        zu_ad(:,:,:) = zu_adin(:,:,:)
        zv_ad(:,:,:) = zv_adin(:,:,:)
        zt_ad(:,:,:) = zt_adin(:,:,:)

        CALL lbc_nfd_adj( zu_ad, 'U',  -1.0_wp )
        CALL lbc_nfd_adj( zv_ad, 'V',  -1.0_wp )
        CALL lbc_nfd_adj( zt_ad, 'T',   1.0_wp )

        zu_adout(:,:,:) = zu_ad(:,:,:)
        zv_adout(:,:,:) = zv_ad(:,:,:)
        zt_adout(:,:,:) = zt_ad(:,:,:)

        zsp2 = sum( PACK(zu_tlin(:,:,:),.TRUE.) * &
         &                       PACK( zu_adout(:,:,:),.TRUE.) )

        CALL mpp_sum( zsp1 )
        CALL mpp_sum( zsp2 )

        cl_name = 'lbc_nfd  U  3d'
        CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

        zsp1 = sum( PACK(zv_tlout(:,:,:),.TRUE.) * &
         &                       PACK( zv_adin(:,:,:),.TRUE.) )

        zsp2 = sum( PACK(zv_tlin(:,:,:),.TRUE.) * &
         &                       PACK( zv_adout(:,:,:),.TRUE.) )

        CALL mpp_sum( zsp1 )
        CALL mpp_sum( zsp2 )
        cl_name = 'lbc_nfd  V  3d'
        CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

        zsp1 = sum( PACK(zt_tlout(:,:,:),.TRUE.) * &
         &                       PACK( zt_adin(:,:,:),.TRUE.) )

        zsp2 = sum( PACK(zt_tlin(:,:,:),.TRUE.) * &
         &                       PACK( zt_adout(:,:,:),.TRUE.) )

        CALL mpp_sum( zsp1 )
        CALL mpp_sum( zsp2 )

      cl_name = 'lbc_nfd  T  3d'
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

   END SUBROUTINE lbc_nfd_3d_adj_tst

   SUBROUTINE lbc_nfd_2d_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE lbc_nfd_2d_adj_tst ***
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
      !! ** Action  :
      !!
      !! History :
      !!        ! 2010-09 (F. Vigilant)
      !!-----------------------------------------------------------------------
      !! * Modules used
   USE gridrandom    , ONLY: & ! Random Gaussian noise on grids
      & grid_random
   USE dotprodfld    , ONLY: & ! Computes dot product for 3D and 2D fields
      & dot_product
   USE lbcnfd,         ONLY: &
      & lbc_nfd
   USE tstool_tam    , ONLY: &
      & stdu,                &
      & stdv,                &
      & stdt,                &
      & prntst_adj
   USE dom_oce       , ONLY: & ! Ocean space and time domain
      & tmask,               &
      & umask,               &
      & vmask,               &
      & mig,                 &
      & mjg,                 &
      & nldi,                &
      & nldj,                &
      & nlei,                &
      & nlej
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt        ! Output unit

      !! * Local declarations
      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zu_tlin,    & ! Tangent input: ua_tl
         & zv_tlin,    & ! Tangent input: va_tl
         & zt_tlin,    & ! Tangent input: ub_tl
         & zu_tlout,   & ! Tangent output: ua_tl
         & zv_tlout,   & ! Tangent output: va_tl
         & zt_tlout,   & ! Tangent output: ua_tl
         & zu_adin,    & ! Tangent output: ua_ad
         & zv_adin,    & ! Tangent output: va_ad
         & zt_adin,    & ! Adjoint input: ua_ad
         & z_adin,    & ! Adjoint input: va_ad
         & zu_adout,   & ! Adjoint output: ua_ad
         & zv_adout,   & ! Adjoint output: va_ad
         & zt_adout,   & ! Adjoint oputput: ub_ad
         & znu            ! 3D random field for u

      REAL(wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zu_tl,    & ! Tangent input: ua_tl
         & zv_tl,    & ! Tangent input: va_tl
         & zt_tl,    & ! Tangent input: ub_tl
         & zu_ad,    & ! Tangent output: ua_ad
         & zv_ad,    & ! Tangent output: va_ad
         & zt_ad
      REAL(wp) :: &
         & zsp1,    &   ! scalar product involving the tangent routine
         & zsp2         ! scalar product involving the adjoint routine
      INTEGER, DIMENSION(jpi,jpj) :: &
         & iseed_2d    ! 2D seed for the random number generator
      CHARACTER(len=1) :: &
         & cd_type
      INTEGER :: &
         & indic, &
         & istp
      INTEGER :: &
         & ji, &
         & jj, &
         & jk, &
         & kmod, &
         & jstp
      CHARACTER (LEN=14) :: &
         & cl_name
      INTEGER ::             &
         & zijpj, ziglo, zip, zdiff

      ! Allocate memory

      SELECT CASE ( jpni )
      CASE ( 1 )     ;   zijpj = nlcj      ! 1 proc only  along the i-direction
      CASE DEFAULT   ;   zijpj = 4         ! several proc along the i-direction
      END SELECT

      ALLOCATE( &
         & zu_tlin(jpiglo,zijpj),  &
         & zv_tlin(jpiglo,zijpj),  &
         & zt_tlin(jpiglo,zijpj),  &
         & zu_tlout(jpiglo,zijpj), &
         & zv_tlout(jpiglo,zijpj), &
         & zt_tlout(jpiglo,zijpj), &
         & zu_adin(jpiglo,zijpj),  &
         & zv_adin(jpiglo,zijpj),  &
         & zt_adin(jpiglo,zijpj),  &
         & zu_adout(jpiglo,zijpj), &
         & zv_adout(jpiglo,zijpj), &
         & zt_adout(jpiglo,zijpj), &
         & znu(jpi,jpj)         &
         & )

      ALLOCATE( &
         & zu_tl(jpiglo,zijpj),  &
         & zv_tl(jpiglo,zijpj),  &
         & zt_tl(jpiglo,zijpj),  &
         & zu_ad(jpiglo,zijpj),  &
         & zv_ad(jpiglo,zijpj),  &
         & zt_ad(jpiglo,zijpj)   &
         & )

        zu_tlin (:,:) = 0.0_wp
        zv_tlin (:,:) = 0.0_wp
        zt_tlin (:,:) = 0.0_wp
        zu_tlout(:,:) = 0.0_wp
        zv_tlout(:,:) = 0.0_wp
        zt_tlout(:,:) = 0.0_wp
        zu_adin (:,:) = 0.0_wp
        zv_adin (:,:) = 0.0_wp
        zt_adin (:,:) = 0.0_wp
        zu_adout(:,:) = 0.0_wp
        zv_adout(:,:) = 0.0_wp
        zt_adout(:,:) = 0.0_wp

        zu_tl (:,:) = 0.0_wp
        zv_tl (:,:) = 0.0_wp
        zt_tl (:,:) = 0.0_wp
        zu_ad (:,:) = 0.0_wp
        zv_ad (:,:) = 0.0_wp
        zt_ad (:,:) = 0.0_wp

        ziglo = INT(jpiglo / (nlei - nldi) )
        zdiff = nlei - nldi
        !--------------------------------------------------------------------
        ! Initialize the tangent input with random noise: dx
        !--------------------------------------------------------------------
        CALL grid_random( znu, 'U', 0.0_wp, stdu )
        DO jj = 1, 4
           DO ji = nldi, nlei
              DO zip = 1, ziglo
                 zu_tlin(ji+(zip-1)*zdiff,jj) = znu(ji,jj)
              END DO
           END DO
        END DO
        CALL grid_random( znu, 'V', 0.0_wp, stdu )
        DO jj = 1, 4
           DO ji = nldi, nlei
              DO zip = 1, ziglo
                 zv_tlin(ji+(zip-1)*zdiff,jj) = znu(ji,jj)
              END DO
           END DO
        END DO
        CALL grid_random(  znu, 'T', 0.0_wp, stdt )
        DO jj = 1, 4
           DO ji = nldi, nlei
              DO zip = 1, ziglo
                 zt_tlin(ji+(zip-1)*zdiff,jj) = znu(ji,jj)
              END DO
           END DO
        END DO

        !--------------------------------------------------------------------
        ! Call the tangent routine: dy = L dx
        !--------------------------------------------------------------------
        zu_tl(:,:) = zu_tlin(:,:)
        zv_tl(:,:) = zv_tlin(:,:)
        zt_tl(:,:) = zt_tlin(:,:)

        CALL lbc_nfd( zu_tl, 'U',  -1.0_wp )
        CALL lbc_nfd( zv_tl, 'V',  -1.0_wp )
        CALL lbc_nfd( zt_tl, 'T',   1.0_wp )
        zu_tlout(:,:) = zu_tl(:,:)
        zv_tlout(:,:) = zv_tl(:,:)
        zt_tlout(:,:) = zt_tl(:,:)

        DO jj = 1,4
           DO ji = nldi, nlei
              DO zip = 1, ziglo
              zu_adin(ji+(zip-1)*zdiff,jj) = zu_tlout(ji+(zip-1)*ziglo,jj) &
                 &               * e1u(ji,jj) * e2u(ji,jj) &!* fse3u(ji,jj) &
                 &               * umask(ji,jj,1)
              zv_adin(ji+(zip-1)*zdiff,jj) = zv_tlout(ji+(zip-1)*ziglo,jj) &
                 &               * e1v(ji,jj) * e2v(ji,jj) &!* fse3v(ji,jj) &
                 &               * vmask(ji,jj,1)
              zt_adin(ji+(zip-1)*zdiff,jj) = zt_tlout(ji+(zip-1)*ziglo,jj) &
                 &               * e1t(ji,jj) * e2t(ji,jj) &!* fse3t(ji,jj) &
                 &               * tmask(ji,jj,1)
              END DO
           END DO
        END DO

        zsp1 = sum( PACK(zu_tlout(:,:),.TRUE.) * &
         &                       PACK( zu_adin(:,:),.TRUE.) )
        zu_ad(:,:) = zu_adin(:,:)
        zv_ad(:,:) = zv_adin(:,:)
        zt_ad(:,:) = zt_adin(:,:)

        CALL lbc_nfd_adj( zu_ad, 'U',  -1.0_wp )
        CALL lbc_nfd_adj( zv_ad, 'V',  -1.0_wp )
        CALL lbc_nfd_adj( zt_ad, 'T',   1.0_wp )

        zu_adout(:,:) = zu_ad(:,:)
        zv_adout(:,:) = zv_ad(:,:)
        zt_adout(:,:) = zt_ad(:,:)

        zsp2 = sum( PACK(zu_tlin(:,:),.TRUE.) * &
         &                       PACK( zu_adout(:,:),.TRUE.) )

        CALL mpp_sum( zsp1 )
        CALL mpp_sum( zsp2 )

        cl_name = 'lbc_nfd  U  2d'
        CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

        zsp1 = sum( PACK(zv_tlout(:,:),.TRUE.) * &
         &                       PACK( zv_adin(:,:),.TRUE.) )

        zsp2 = sum( PACK(zv_tlin(:,:),.TRUE.) * &
         &                       PACK( zv_adout(:,:),.TRUE.) )
        
        CALL mpp_sum( zsp1 )
        CALL mpp_sum( zsp2 )

        cl_name = 'lbc_nfd  V  2d'
        CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

        zsp1 = sum( PACK(zt_tlout(:,:),.TRUE.) * &
         &                       PACK( zt_adin(:,:),.TRUE.) )

        zsp2 = sum( PACK(zt_tlin(:,:),.TRUE.) * &
         &                       PACK( zt_adout(:,:),.TRUE.) )

        CALL mpp_sum( zsp1 )
        CALL mpp_sum( zsp2 )

        cl_name = 'lbc_nfd  T  2d'
        CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )

   END SUBROUTINE lbc_nfd_2d_adj_tst

   SUBROUTINE lbc_nfd_adj_tst( kumadt )
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE lbc_nfd_adj_tst ***
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
      !! ** Action  :
      !!
      !! History :
      !!        ! 2010-09 (F. Vigilant)
      !!-----------------------------------------------------------------------
      !! * Modules used
      !! * Arguments
      INTEGER, INTENT(IN) :: &
         & kumadt        ! Output unit

      CALL lbc_nfd_3d_adj_tst( kumadt )
      CALL lbc_nfd_2d_adj_tst( kumadt )

   END SUBROUTINE lbc_nfd_adj_tst
   !!======================================================================
END MODULE lbcnfd_tam
