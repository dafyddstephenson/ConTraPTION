MODULE paresp
   !!======================================================================
   !!                       ***  MODULE paresp ***
   !! NEMOVAR : Weights for an energy-type scalar product
   !!======================================================================

   !!----------------------------------------------------------------------
   !!  par_esp : Compute and store weights for an energy-type scalar product
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind
   USE par_oce
   USE eosbn2
   USE phycst
   USE zdf_oce
   USE dom_oce
   USE lib_mpp          ! MPP library
   USE in_out_manager   ! I/O stuff
   USE wrk_nemo         ! Memory Allocation

   IMPLICIT NONE

   !! * Routine accessibility
   PRIVATE

   PUBLIC &
      & par_esp       !: Compute and store energy weights

   !! * Share module variables
   REAL(wp), DIMENSION(:), ALLOCATABLE, PUBLIC ::   &
      & wesp_t,   &   !: Normalized energy weights for temperature
      & wesp_s        !: Normalized energy weights for salinity
   REAL(wp), PUBLIC ::   &
      & wesp_u,   &   !: Normalized energy weights for velocity
      & wesp_ssh, &   !: Normalized energy weights for SSH
      & wesp_tau, &   !: Normalized energy weights for windstress
      & wesp_q,   &   !: Normalized energy weights for heat flux
      & wesp_emp      !: Normalized energy weights for EmP

CONTAINS
   INTEGER FUNCTION par_esp_alloc()
      !!----------------------------------------------------------------------
      !!                ***  FUNCTION exa_mpl_alloc  ***
      !!----------------------------------------------------------------------
      ALLOCATE( wesp_t(jpk) , wesp_s(jpk) , STAT= par_esp_alloc )
      !
      IF( lk_mpp             )   CALL mpp_sum ( par_esp_alloc )
      IF( par_esp_alloc /= 0 )   CALL ctl_warn('par_esp_alloc: failed to allocate arrays')
      !
   END FUNCTION par_esp_alloc

   SUBROUTINE par_esp
      !!-----------------------------------------------------------------------
      !!
      !!                  ***  ROUTINE par_esp  ***
      !!
      !! ** Purpose : Compute and store weights for an energy-type scalar
      !!              product.
      !!
      !! ** Method  :
      !!      The weighting coefficients are computed from an analytical
      !!      linear density profile which is similar to the one used in
      !!      default option in OPA direct.
      !!
      !!      The weighting coefficients for the forcing fields are estimated
      !!      from the surface boundary conditions.
      !!
      !!      The weighting coefficients for the velocity variables are
      !!      normalized to one.
      !!
      !!       wesp_t(k) = ( g / N(k) )^2 x alpha^2
      !!       wesp_s(k) = ( g / N(k) )^2 x beta^2
      !!       wesp_u    = 1
      !!       wesp_ssh  = g
      !!       wesp_tau  = ( ( e3w(1) / ( avu x rho0 ) )^2
      !!       wesp_q    = ( ( g / N(1) x alpha
      !!                  x e3w(1) / ( avT x rho0 x cp ) )^2
      !!       wesp_emp  = ( ( g / N(1) x beta
      !!                  x e3w(1) x S(1) / avT )^2
      !!
      !!      where
      !!
      !!        N(k)   is the Brunt-Vaisala frequency (depth-dependent)
      !!        alpha  is the thermal expansion coefficient
      !!        beta   is the salinity expansion coefficient
      !!        g      is gravity
      !!        avu    is vertical eddy viscosity coefficent for dynamics
      !!        avT    is vertical eddy diffusivity coefficent for tracers
      !!        e3w(1) is the vertical scale factor at top w-point
      !!        rho0   is a reference density
      !!        S(1)   is a reference salinity at the surface
      !!
      !!      The complete scalar product < ; > for the state vector
      !!
      !!       dx = ( du , dv , dT , dS , deta , dq , demp , dtaux , dtauy )
      !!
      !!      involves the volume/area elements:
      !!
      !!      < dx ; dx > = sum_i,j (  deta^2  * e1T * e2T           * wesp_ssh
      !!                             + dq^2    * e1T * e2T * e3w(1)  * wesp_q
      !!                             + demp^2  * e1T * e2T * e3w(1)  * wesp_emp
      !!                             + dtaux^2 * e1u * e2u * e3uw(1) * wesp_tau
      !!                             + dtauy^2 * e1v * e2v * e3vw(1) * wesp_tau
      !!                    + sum_k (  du^2 * e1u * e2u * e3u * wesp_u
      !!                             + dv^2 * e1v * e2v * e3v * wesp_u
      !!                             + dT^2 * e1T * e2T * e3T * wesp_T
      !!                             + dS^2 * e1T * e2T * e3T * wesp_S )    )
      !!
      !! ** Action  :
      !!
      !! History :
      !!        ! 97-06 (A. Weaver, J. Vialard) OPAVAR version
      !!        ! 07-11 (A. Weaver) NEMOVAR version based on OPAVAR paresp.F
      !!        ! 2011-07 (J. While) Adapted to deal with 2D grids
      !!-----------------------------------------------------------------------
      !! * Modules used

      !! * Arguments

      !! * Local declarations
      REAL(wp) :: zgrau0,  zk, zt0, zs0, zrau0
      REAL(wp), POINTER, DIMENSION(:) :: ztan,  zsan
      REAL(wp), POINTER, DIMENSION(:) :: zrhan, zbn2an
      INTEGER :: jk, ierr
      !!----------------------------------------------------------------------
      !
      ierr = par_esp_alloc()
      CALL wrk_alloc( jpk, ztan, zsan, zrhan, zbn2an )
      !
      !--------------------------------------------------------------------
      ! Local constant initialization
      !--------------------------------------------------------------------

      zt0   = 19.0_wp    ! Reference temperature
      zs0   = 35.0_wp    ! Reference salinity
      zrau0 = 1025.0_wp  ! Reference density

      zgrau0 = -1.0_wp * grav / zrau0

      !--------------------------------------------------------------------
      ! Initialize rho with an analytical profile
      !--------------------------------------------------------------------

      !  Define analytical temperature profile

      DO jk = 1, jpk
         ztan(jk) =  7.5_wp * ( 1.0_wp - TANH( ( gdept_0(jk) - 80.0_wp ) &
            &                 / 30.0_wp ) ) &
            &      + 10.0_wp * ( 6000.0_wp - gdept_0(jk) ) / 6000.0_wp
         ztan(jk) = ABS( ztan(jk) )
      END DO

      !  Define analytical salinity profile

      DO jk = 1, jpk
         zsan(jk) = 35.50_wp
      END DO

      !  Define analytical density profile

      DO jk = 1, jpk
         zrhan(jk) =  zrau0 * ( 1.0_wp - rn_alpha * ( ztan(jk) - zt0 ) &
            &                          + rn_beta  * ( zsan(jk) - zs0 ) )
      END DO

      !--------------------------------------------------------------------
      ! Calculate N^2 at T and S points
      !--------------------------------------------------------------------

      IF ( jpk > 2 ) THEN
         DO jk = 2, jpkm1
            zbn2an(jk) = ( zgrau0 / 2.0_wp ) &
               & * (   ( ( zrhan(jk-1) - zrhan(jk  ) ) / e3w_0(jk  ) ) &
               &     + ( ( zrhan(jk  ) - zrhan(jk+1) ) / e3w_0(jk+1) ) )
         END DO
      ELSE
         zbn2an(2) = 1._wp
      ENDIF

      zbn2an(1)   = zbn2an(2)
      zbn2an(jpk) = zbn2an(jpkm1)

      !--------------------------------------------------------------------
      ! Calculate energy-type weights
      !--------------------------------------------------------------------

      DO jk = 1, jpk
         zk = ( grav * grav ) / zbn2an(jk)
         wesp_t(jk) = zk * rn_alpha * rn_alpha
         wesp_s(jk) = zk * rn_beta  * rn_beta
      END DO

      wesp_u   = 1.0_wp
      wesp_ssh = grav
      wesp_tau = ( e3w_0(1) / ( rn_avm0 * zrau0 ) )**2
      wesp_q   = wesp_t(1) * ( e3w_0(1) / ( rn_avt0 * zrau0 * rcp ) )**2
      wesp_emp = wesp_s(1) * ( e3w_0(1) * zsan(1) / rn_avt0  )**2

      !--------------------------------------------------------------------
      ! Print
      !--------------------------------------------------------------------

      IF (lwp) THEN
         WRITE(numout,*)
         WRITE(numout,*) ' par_esp: Analytical T, S, rho, N^2 profiles'
         WRITE(numout,*) ' -------'
         WRITE(numout,*)
         WRITE(numout,995)
         WRITE(numout,996)
         DO jk = 1, jpk
            IF(lwp)WRITE(numout,1007) jk, ztan(jk), zsan(jk), &
               &                      zrhan(jk), zbn2an(jk)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' par_esp: Energy scalar product weights', &
            &            ' for T and S'
         WRITE(numout,*) ' -------'
         WRITE(numout,*)
         WRITE(numout,998)
         WRITE(numout,999)
         DO jk = 1, jpk
            WRITE(numout,1006) jk, wesp_t(jk), wesp_s(jk)
         END DO
         WRITE(numout,*)
         WRITE(numout,*) ' par_esp: Energy scalar product weights', &
            &            ' for u, SSH, tau, Q and EmP'
         WRITE(numout,*) ' -------'
         WRITE(numout,*)
         WRITE(numout,*) '          wesp_u   = ', wesp_u
         WRITE(numout,*) '          wesp_ssh = ', wesp_ssh
         WRITE(numout,*) '          wesp_tau = ', wesp_tau
         WRITE(numout,*) '          wesp_q   = ', wesp_q
         WRITE(numout,*) '          wesp_emp = ', wesp_emp

 995     FORMAT('    level        T              S    ', &
            &   '        rho             N^2 ')
 996     FORMAT('    -----     --------       --------', &
            &   '     ----------     ----------')
 998     FORMAT('    level   wesp_T        wesp_S   ')
 999     FORMAT('    ----- ----------   ------------')
 1005    FORMAT(4X,I2,3X,E12.5,1X,E12.5)
 1006    FORMAT(4X,I2,3X,E12.5,1X,E12.5)
 1007    FORMAT(4X,I2,6X,F10.5,5X,F10.5,5X,F10.5,5X,E12.5)
         CALL flush( numout )

      ENDIF

      DO jk = 1, jpk
         IF ( zbn2an(jk) < 0.0_wp ) THEN
            IF(lwp)WRITE(numout,990) zbn2an(jk), jk
990         FORMAT(' paresp : Error unstable density profile;', &
               &   ' zbn2an = ',E12.5,' at level ',I3)
!            CALL ctl_stop( 'paresp; unstable density profile' )
             IF(lwp)WRITE(numout,*)'WARNING'
	     IF(lwp)WRITE(numout,*)'paresp; unstable density profile'
         ENDIF
      ENDDO
      CALL wrk_dealloc( jpk, ztan, zsan, zrhan, zbn2an )

   END SUBROUTINE par_esp

END MODULE paresp
