MODULE solsor_tam
   !!======================================================================
   !!                     ***  MODULE  solsor  ***
   !! Ocean solver :  Tangent linear and Adjoint of Successive Over-Relaxation solver
   !!=====================================================================
   !! History :  OPA  ! 1990-10  (G. Madec)  Original code
   !!            7.1  ! 1993-04  (G. Madec)  time filter
   !!                 ! 1996-05  (G. Madec)  merge sor and pcg formulations
   !!                 ! 1996-11  (A. Weaver)  correction to preconditioning
   !!   NEMO     1.0  ! 2003-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
   !!            2.0  ! 2005-09  (R. Benshila, G. Madec)  MPI optimization
   !!----------------------------------------------------------------------
   !!
   !!----------------------------------------------------------------------
   !!   sol_sor     : Red-Black Successive Over-Relaxation solver
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_oce
   USE in_out_manager
   USE sol_oce
   USE solmat
   USE lib_mpp
   USE lbclnk
   USE lbclnk_tam
   USE sol_oce_tam
   USE dom_oce
   USE gridrandom
   USE dotprodfld
   USE tstool_tam
   USE lib_fortran
   USE wrk_nemo
   USE timing
   !
   IMPLICIT NONE
   PRIVATE
   !
   !! * Routine accessibility
   PUBLIC sol_sor_adj          !
   PUBLIC sol_sor_tan          !
   PUBLIC sol_sor_adj_tst      ! called by tamtst.F90
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id$
   !! Software governed by the CeCILL licence     (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
CONTAINS

   SUBROUTINE sol_sor_tan( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_sor_tan : TL of sol_sor  ***
      !!
      !! ** Purpose :   Solve the ellipic equation for the barotropic stream
      !!      function system (lk_dynspg_rl=T) or the transport divergence
      !!      system (lk_dynspg_flt=T) using a red-black successive-over-
      !!      relaxation method.
      !!       In the former case, the barotropic stream function trend has a
      !!     zero boundary condition along all coastlines (i.e. continent
      !!     as well as islands) while in the latter the boundary condition
      !!     specification is not required.
      !!       This routine provides a MPI optimization to the existing solsor
      !!     by reducing the number of call to lbc.
      !!
      !! ** Method  :   Successive-over-relaxation method using the red-black
      !!      technique. The former technique used was not compatible with
      !!      the north-fold boundary condition used in orca configurations.
      !!      Compared to the classical sol_sor, this routine provides a
      !!      mpp optimization by reducing the number of calls to lnc_lnk
      !!      The solution is computed on a larger area and the boudary
      !!      conditions only when the inside domain is reached.
      !!
      !! References :
      !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
      !!      Beare and Stevens 1997 Ann. Geophysicae 15, 1369-1377
      !!
      !! History of the direct routine:
      !!        !  90-10  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!   7.1  !  93-04  (G. Madec)  time filter
      !!        !  96-05  (G. Madec)  merge sor and pcg formulations
      !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
      !!   9.0  !  05-09  (R. Benshila, G. Madec)  MPI optimization
      !! History of the  T&A routine:
      !!        !  96-11  (A. Weaver)  correction to preconditioning
      !!   8.2  !  03-02  (C. Deltel) OPAVAR tangent-linear version
      !!   9.0  !  07-09  (K. Mogensen) tangent of the 03-04 version
      !!   9.0  !  09-02  (A. Vidard)   tangent of the 05-09 version
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in    ) ::   kt       ! Current timestep.
      INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                      ! gence is not reached: the model is
      !                                      ! stopped in step
      !                                      ! set to zero before the call of solsor
      !! * Local declarations
      INTEGER  ::   ji, jj, jn               ! dummy loop indices
      INTEGER  ::   ishift, icount, istp
      REAL(wp) ::   ztmp, zres, zres2

      INTEGER  ::   ijmppodd, ijmppeven
      INTEGER  ::   ijpr2d
      REAL(wp), POINTER, DIMENSION(:,:) ::   ztab                 ! 2D workspace
      !!----------------------------------------------------------------------
      !
      IF( nn_timing == 1 )  CALL timing_start('sol_sor_tan')
      !
      CALL wrk_alloc( jpi, jpj, ztab )
      !
      ijmppeven = MOD(nimpp+njmpp+jpr2di+jpr2dj  ,2)
      ijmppodd  = MOD(nimpp+njmpp+jpr2di+jpr2dj+1,2)
      ijpr2d = MAX(jpr2di,jpr2dj)
      icount = 0
      !                                                       ! ==============
      DO jn = 1,  nn_nmax                                         ! Iterative loop
         !                                                    ! ==============
         ! applied the lateral boundary conditions
         IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx_tl, c_solver_pt, 1.0_wp )
         ! Residuals
         ! ---------
         ! Guess black update
         DO jj = 2-jpr2dj, nlcj - 1 + jpr2dj
            ishift = MOD( jj-ijmppodd-jpr2dj, 2 )
            DO ji = 2-jpr2di+ishift, nlci - 1 + jpr2di, 2
               ztmp =                  gcb_tl(ji  ,jj  )   &
                  &   - gcp(ji,jj,1) * gcx_tl(ji  ,jj-1)   &
                  &   - gcp(ji,jj,2) * gcx_tl(ji-1,jj  )   &
                  &   - gcp(ji,jj,3) * gcx_tl(ji+1,jj  )   &
                  &   - gcp(ji,jj,4) * gcx_tl(ji  ,jj+1)
               ! Estimate of the residual
               zres = ztmp - gcx_tl(ji,jj)
               gcr_tl(ji,jj) = zres * gcdmat(ji,jj) * zres
               ! Guess update
               gcx_tl(ji,jj) = rn_sor * ztmp + ( 1.0_wp - rn_sor ) * gcx_tl(ji,jj)
            END DO
         END DO
         icount = icount + 1
         ! applied the lateral boundary conditions
         IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e( gcx_tl, c_solver_pt, 1.0_wp )
         ! Guess red update
         DO jj = 2-jpr2dj, nlcj-1+jpr2dj
            ishift = MOD( jj-ijmppeven-jpr2dj, 2 )
            DO ji = 2-jpr2di+ishift, nlci-1+jpr2di, 2
               ztmp =                  gcb_tl(ji  ,jj  )   &
                  &   - gcp(ji,jj,1) * gcx_tl(ji  ,jj-1)   &
                  &   - gcp(ji,jj,2) * gcx_tl(ji-1,jj  )   &
                  &   - gcp(ji,jj,3) * gcx_tl(ji+1,jj  )   &
                  &   - gcp(ji,jj,4) * gcx_tl(ji  ,jj+1)
               ! Estimate of the residual
               zres = ztmp - gcx_tl(ji,jj)
               gcr_tl(ji,jj) = zres * gcdmat(ji,jj) * zres
               ! Guess update
               gcx_tl(ji,jj) = rn_sor * ztmp + ( 1.0_wp - rn_sor ) * gcx_tl(ji,jj)
            END DO
         END DO
         icount = icount + 1
         ! test of convergence
         IF ( (jn > nn_nmin .AND. MOD( jn-nn_nmin, nn_nmod ) == 0) .OR. jn==nn_nmax) THEN
            SELECT CASE ( nn_sol_arp )
            CASE ( 0 )
               ! absolute precision (maximum value of the residual)
               zres2 = MAXVAL( gcr_tl(2:nlci - 1,2:nlcj - 1) )
               IF( lk_mpp )  CALL mpp_max( zres2 ) ! max over the global domain
               ! test of convergence
               res = SQRT( zres2 )
               IF( zres2 < rn_resmax .OR. jn == nn_nmax ) THEN
                  niter = jn
                  ncut = 999
                  ! Store number of iterations for adjoint computation
                  istp = kt - nit000 + 1
                  nitsor(istp) = niter
               ENDIF
            CASE ( 1 )                 ! relative precision
               ztab(:,:) = 0.0_wp
               ztab(2:nlci-1,2:nlcj-1) = gcr_tl(2:nlci-1,2:nlcj-1)
               rnorme = glob_sum(ztab)
               ! test of convergence
               res = SQRT( rnorme )
               IF( rnorme < epsr .OR. jn == nn_nmax ) THEN
                  niter = jn
                  ncut = 999
                  ! Store number of iterations for adjoint computation
                  istp = kt - nit000 + 1
                  nitsor(istp) = niter
               ENDIF
            END SELECT
            !****
            IF(lwp)WRITE(numsol,9300) jn, res, sqrt( epsr ) / eps
9300        FORMAT('          niter :',i6,' res :',e20.10,' b :',e20.10)
            !****
         ENDIF
         ! indicator of non-convergence or explosion
	      IF( jn == nn_nmax ) nitsor(istp) = jn
         IF( jn == nn_nmax .OR. SQRT(epsr)/eps > 1.e+20 ) kindic = -2
         IF( ncut == 999 ) GOTO 999
         !                                              ! =====================
      END DO                                            ! END of iterative loop
      !                                                 ! =====================
999   CONTINUE
      !  Output in gcx_tl
      !  ----------------
      CALL lbc_lnk_e( gcx_tl, c_solver_pt, 1.0_wp )    ! Lateral BCs
      !
      CALL wrk_dealloc( jpi, jpj, ztab )
      !
      IF( nn_timing == 1 )  CALL timing_stop('sol_sor_tan')
      !
   END SUBROUTINE sol_sor_tan
   SUBROUTINE sol_sor_adj( kt, kindic )
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE sol_sor_adj : adjoint of sol_sor  ***
      !!
      !! ** Purpose :   Solve the ellipic equation for the barotropic stream
      !!      function system (lk_dynspg_rl=T) or the transport divergence
      !!      system (lk_dynspg_flt=T) using a red-black successive-over-
      !!      relaxation method.
      !!       In the former case, the barotropic stream function trend has a
      !!     zero boundary condition along all coastlines (i.e. continent
      !!     as well as islands) while in the latter the boundary condition
      !!     specification is not required.
      !!       This routine provides a MPI optimization to the existing solsor
      !!     by reducing the number of call to lbc.
      !!
      !! ** Method  :   Successive-over-relaxation method using the red-black
      !!      technique. The former technique used was not compatible with
      !!      the north-fold boundary condition used in orca configurations.
      !!      Compared to the classical sol_sor, this routine provides a
      !!      mpp optimization by reducing the number of calls to lbc_lnk
      !!      The solution is computed on a larger area and the boudary
      !!      conditions only when the inside domain is reached.
      !! ** Comments on adjoint routine :
      !!      When the step in a tangent-linear DO loop is an arbitrary
      !!      integer then care must be taken in computing the lower bound
      !!      of the adjoint DO loop; i.e.,
      !!
      !!      If the tangent-linear DO loop is:  low_tl, up_tl, step
      !!
      !!      then the adjoint DO loop is:  low_ad, up_ad, -step
      !!
      !!      where  low_ad = low_tl + step * INT( ( up_tl - low_tl ) / step )
      !!             up_ad  = low_tl
      !!
      !!      NB. If step = 1 then low_ad = up_tl
      !!
      !! References :
      !!      Madec et al. 1988, Ocean Modelling, issue 78, 1-6.
      !!      Beare and Stevens 1997 Ann. Geophysicae 15, 1369-1377
      !!
      !! History of the direct routine:
      !!        !  90-10  (G. Madec)  Original code
      !!        !  91-11  (G. Madec)
      !!   7.1  !  93-04  (G. Madec)  time filter
      !!        !  96-05  (G. Madec)  merge sor and pcg formulations
      !!   9.0  !  03-04  (C. Deltel, G. Madec)  Red-Black SOR in free form
      !!   9.0  !  05-09  (R. Benshila, G. Madec)  MPI optimization
      !! History of the  T&A routine:
      !!        !  96-11  (A. Weaver)  correction to preconditioning
      !!   8.2  !  03-02  (C. Deltel) OPAVAR tangent-linear version
      !!   9.0  !  07-09  (K. Mogensen, A. Weaver) adjoint of the 03-04 version
      !!   9.0  !  09-02  (A. Vidard)  adjoint of the 05-09 version
      !!----------------------------------------------------------------------
      !! * Arguments
      INTEGER, INTENT( in    ) ::   kt       ! Current timestep.
      INTEGER, INTENT( inout ) ::   kindic   ! solver indicator, < 0 if the conver-
      !                                      ! gence is not reached: the model is
      !                                      ! stopped in step
      !                                      ! set to zero before the call of solsor
      !! * Local declarations
      INTEGER  ::   ji, jj, jn               ! dummy loop indices
      INTEGER  ::   ishift, icount, istp, iter, ilower
      REAL(wp) ::   ztmpad

      INTEGER  ::   ijmppodd, ijmppeven
      INTEGER  ::   ijpr2d
      !!----------------------------------------------------------------------
      IF( nn_timing == 1 )  CALL timing_start('sol_sor_adj')
      !
      ijmppeven = MOD(nimpp+njmpp+jpr2di+jpr2dj,2)
      ijmppodd  = MOD(nimpp+njmpp+jpr2di+jpr2dj+1,2)
      ijpr2d = MAX(jpr2di,jpr2dj)
      !
      ! Fixed number of iterations
      istp = kt - nit000 + 1
      iter = nitsor(istp)
      icount = iter * 2
      !  Output in gcx_ad
      !  ----------------
      CALL lbc_lnk_e_adj( gcx_ad, c_solver_pt, 1.0_wp )    ! Lateral BCs
      !                                                    ! ==============
      DO jn = iter, 1, -1                                  ! Iterative loop
         !                                                 ! ==============
         ! Guess red update
         DO jj =  nlcj-1+jpr2dj, 2-jpr2dj, -1
            ishift = MOD( jj-ijmppeven-jpr2dj, 2 )
            ! this weird computation is to cope with odd end of loop in the tangent
            ilower = 2-jpr2dj+ishift + 2 * INT( ( ( nlci-1+jpr2dj )-( 2-jpr2dj+ishift ) ) / 2 )
            DO ji = ilower, 2-jpr2dj+ishift, -2
               ! Guess update
               ztmpad = rn_sor * gcx_ad(ji,jj)
               gcx_ad(ji  ,jj  ) = gcx_ad(ji  ,jj  ) * ( 1.0_wp - rn_sor )

               gcb_ad(ji  ,jj  ) = gcb_ad(ji  ,jj  ) + ztmpad
               gcx_ad(ji  ,jj-1) = gcx_ad(ji  ,jj-1) - ztmpad * gcp(ji,jj,1)
               gcx_ad(ji-1,jj  ) = gcx_ad(ji-1,jj  ) - ztmpad * gcp(ji,jj,2)
               gcx_ad(ji+1,jj  ) = gcx_ad(ji+1,jj  ) - ztmpad * gcp(ji,jj,3)
               gcx_ad(ji  ,jj+1) = gcx_ad(ji  ,jj+1) - ztmpad * gcp(ji,jj,4)
            END DO
         END DO
	      icount = icount - 1
         ! applied the lateral boundary conditions
         IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e_adj( gcx_ad, c_solver_pt, 1.0_wp )   ! Lateral BCs
         ! Residus
         ! -------
         ! Guess black update
         DO jj = nlcj-1+jpr2dj, 2-jpr2dj, -1
            ishift = MOD( jj-ijmppodd-jpr2dj, 2 )
            ilower = 2-jpr2dj+ishift + 2 * INT( ( ( nlci-1+jpr2dj )-( 2-jpr2dj+ishift ) ) / 2 )
            DO ji = ilower, 2-jpr2dj+ishift, -2
               ! Guess update
               ztmpad = rn_sor * gcx_ad(ji,jj)
               gcx_ad(ji  ,jj  ) = gcx_ad(ji  ,jj  ) * ( 1.0_wp - rn_sor )

               gcb_ad(ji  ,jj  ) = gcb_ad(ji  ,jj  ) + ztmpad
               gcx_ad(ji  ,jj-1) = gcx_ad(ji  ,jj-1) - ztmpad * gcp(ji,jj,1)
               gcx_ad(ji-1,jj  ) = gcx_ad(ji-1,jj  ) - ztmpad * gcp(ji,jj,2)
               gcx_ad(ji+1,jj  ) = gcx_ad(ji+1,jj  ) - ztmpad * gcp(ji,jj,3)
               gcx_ad(ji  ,jj+1) = gcx_ad(ji  ,jj+1) - ztmpad * gcp(ji,jj,4)
            END DO
         END DO
	     icount = icount - 1
         ! applied the lateral boundary conditions
         IF( MOD(icount,ijpr2d+1) == 0 ) CALL lbc_lnk_e_adj( gcx_ad, c_solver_pt, 1.0_wp )   ! Lateral BCs
         !                                              ! =====================
      END DO                                            ! END of iterative loop
      !                                                 ! =====================
      IF( nn_timing == 1 )  CALL timing_stop('sol_sor_adj')

   END SUBROUTINE sol_sor_adj
   SUBROUTINE sol_sor_adj_tst( kumadt )
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
         & jn,    &
         & jk,    &
         & kindic,&        ! flags fo solver convergence
         & kmod,  &        ! frequency of test for the SOR solver
         & kt              ! number of iteration
      REAL(KIND=wp) :: &
         & zsp1,         & ! scalar product involving the tangent routine
         & zsp2            ! scalar product involving the adjoint routine
      REAL(KIND=wp), DIMENSION(:,:), ALLOCATABLE :: &
         & zgcb_tlin ,     & ! Tangent input
         & zgcx_tlin ,     & ! Tangent input
         & zgcx_tlout,     & ! Tangent output
         & zgcx_adin ,     & ! Adjoint input
         & zgcb_adout,     & ! Adjoint output
         & zgcx_adout,     & ! Adjoint output
         & zr             ! 3D random field
      CHARACTER(LEN=14) :: cl_name
      ! Allocate memory


      ALLOCATE( &
         & zgcb_tlin( jpi,jpj),     &
         & zgcx_tlin( jpi,jpj),     &
         & zgcx_tlout(jpi,jpj),     &
         & zgcx_adin( jpi,jpj),     &
         & zgcx_adout(jpi,jpj),     &
         & zgcb_adout(jpi,jpj),     &
         & zr(        jpi,jpj)      &
         & )

      ! Initialize the matrix of the elliptic equation

      CALL sol_mat( nit000 + 1 )

      !==================================================================
      ! 1) dx = ( un_tl, vn_tl, hdivn_tl ) and
      !    dy = ( hdivb_tl, hdivn_tl )
      !==================================================================

      !--------------------------------------------------------------------
      ! Reset the tangent and adjoint variables
      !--------------------------------------------------------------------
      zgcb_tlin( :,:) = 0.0_wp
      zgcx_tlin( :,:) = 0.0_wp
      zgcx_tlout(:,:) = 0.0_wp
      zgcx_adin( :,:) = 0.0_wp
      zgcx_adout(:,:) = 0.0_wp
      zgcb_adout(:,:) = 0.0_wp
      zr(        :,:) = 0.0_wp
      !--------------------------------------------------------------------
      ! Initialize the tangent input with random noise: dx
      !--------------------------------------------------------------------
      kt=nit000
      kindic=0
!      kmod = nn_nmod  ! store frequency of test for the SOR solver
!      nn_nmod = 1     ! force frequency to one (remove adj_tst dependancy to nn_nmin)

      CALL grid_random( zr, c_solver_pt, 0.0_wp, stdgc )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgcb_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO
      CALL grid_random( zr, c_solver_pt, 0.0_wp, stdgc )
      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgcx_tlin(ji,jj) = zr(ji,jj)
         END DO
      END DO
      ncut = 1 ! reinitialize the solver convergence flag
      gcr_tl(:,:) = 0.0_wp
      gcb_tl(:,:) = zgcb_tlin(:,:)
      gcx_tl(:,:) = zgcx_tlin(:,:)
      CALL sol_sor_tan(kt, kindic)
      zgcx_tlout(:,:) = gcx_tl(:,:)

      !--------------------------------------------------------------------
      ! Initialize the adjoint variables: dy^* = W dy
      !--------------------------------------------------------------------

      DO jj = nldj, nlej
         DO ji = nldi, nlei
            zgcx_adin(ji,jj) = zgcx_tlout(ji,jj)         &
               &               * e1t(ji,jj) * e2t(ji,jj) &
               &               * tmask(ji,jj,1)
         END DO
      END DO
      !--------------------------------------------------------------------
      ! Compute the scalar product: ( L dx )^T W dy
      !--------------------------------------------------------------------
      zsp1 = DOT_PRODUCT( zgcx_tlout, zgcx_adin )
      !--------------------------------------------------------------------
      ! Call the adjoint routine: dx^* = L^T dy^*
      !--------------------------------------------------------------------
      gcb_ad(:,:) = 0.0_wp
      gcx_ad(:,:) = zgcx_adin(:,:)
      CALL sol_sor_adj(kt, kindic)
      zgcx_adout(:,:) = gcx_ad(:,:)
      zgcb_adout(:,:) = gcb_ad(:,:)

      zsp2 =  DOT_PRODUCT( zgcx_tlin, zgcx_adout ) &
         & + DOT_PRODUCT( zgcb_tlin, zgcb_adout )

      cl_name = 'sol_sor_adj  '
      CALL prntst_adj( cl_name, kumadt, zsp1, zsp2 )
!      nn_nmod = kmod  ! restore initial frequency of test for the SOR solver

      nitsor(:) = jp_it0adj

      DEALLOCATE(      &
         & zgcb_tlin,  &
         & zgcx_tlin,  &
         & zgcx_tlout, &
         & zgcx_adin,  &
         & zgcx_adout, &
         & zgcb_adout, &
         & zr          &
         & )
   END SUBROUTINE sol_sor_adj_tst
END MODULE solsor_tam

