MODULE daymod_tam
   !!======================================================================
   !!                       ***  MODULE  daymod_tam  ***
   !! Ocean        :  calendar, tangent and adjoint model version
   !!=====================================================================
   !! History :  OPA  ! 1994-09  (M. Pontaud M. Imbard)  Original code
   !!                 ! 1997-03  (O. Marti)
   !!                 ! 1997-05  (G. Madec)
   !!                 ! 1997-08  (M. Imbard)
   !!   NEMO     1.0  ! 2003-09  (G. Madec)  F90 + nyear, nmonth, nday
   !!                 ! 2004-01  (A.M. Treguier) new calculation based on adatrj
   !!                 ! 2006-08  (G. Madec)  surface module major update
   !! History :
   !!            OPA  ! 1998-2004 (A. Weaver, N. Daget) daytam
   !!            NEMO ! 2005-08  (A. Vidard) skeleton
   !!                 ! 2008-08  (A. Vidard) 04-01 version, based on daytam
   !!                 ! 2009-02  (A. Vidard) 06-08 version
   !!            3.2  ! 2010-03  (F. Vigilant)
   !!            3.2  ! 2012-07  (P.-A. Bouttier) Phasing with 3.4
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   day_tam        : calendar
   !!
   !!           -------------------------------
   !!           ----------- WARNING -----------
   !!
   !!   we assume that the time step is a divisor of the number of second of in a day
   !!             ---> MOD( rday, rdttra(1) ) == 0
   !!
   !!           ----------- WARNING -----------
   !!           -------------------------------
   !!
   !!----------------------------------------------------------------------
   !! * Modules used
   USE par_kind
   USE phycst
   USE dom_oce
   USE in_out_manager
   USE daymod
   USE prtctl
   USE ioipsl, ONLY : ymds2ju

   IMPLICIT NONE
   PRIVATE

   !! * Routine accessibility
   PUBLIC day_tam        ! called by steptan.F90 and stepadj.F90

CONTAINS


   SUBROUTINE day_tam( kt, kindic )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE day_tam  ***
      !!
      !! ** Purpose :   Compute the date with a day iteration IF necessary.
      !!    	    	forward for the tangent linear, backward for the adjoint
      !!
      !! * Arguments
      INTEGER, INTENT( in ) ::   kt      ! ocean time-step indices
      INTEGER, INTENT( in ) ::   kindic  ! forward (0) or backward (1)
      !!----------------------------------------------------------------------




      SELECT CASE ( kindic )
      CASE ( 0 )
         ! ----------------------------------------------
         ! Forward running Calendar for the tangent model
         ! ----------------------------------------------
         CALL day_tan( kt )
      CASE ( 1 )
         ! ----------------------------------------------
         ! Backward running Calendar for the adjoint model
         ! ----------------------------------------------
         CALL day_adj( kt )
      CASE default
         IF (lwp) WRITE(numout,*) 'day_tam called with a wrong kindic: ',kindic
         CALL abort
      END SELECT
   END SUBROUTINE day_tam

   SUBROUTINE day_tan ( kt )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE day_tan  ***
      !!
      !! ** Purpose :   Compute the date with a day iteration IF necessary.
      !!
      !! ** Method  : - ???
      !!
      !! ** Action  : - nyear     : current year
      !!              - nmonth    : current month of the year nyear
      !!              - nday      : current day of the month nmonth
      !!              - nday_year : current day of the year nyear
      !!              - ndastp    : = nyear*10000 + nmonth*100 + nday
      !!              - adatrj    : date in days since the beginning of the run
      !!              - nsec_year : current time of the year (in second since 00h, jan 1st)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt        ! ocean time-step indices
      !
      CHARACTER (len=25) ::   charout
      REAL(wp)           ::   zprec      ! fraction of day corresponding to 0.1 second
      !!----------------------------------------------------------------------
      !
      zprec = 0.1 / rday
      !                                                 ! New time-step
      nsec_year  = nsec_year  + ndt
      nsec_month = nsec_month + ndt
      nsec_week  = nsec_week  + ndt
      nsec_day   = nsec_day   + ndt
      adatrj  = adatrj  + rdttra(1) / rday
      fjulday = fjulday + rdttra(1) / rday
      IF( ABS(fjulday - REAL(NINT(fjulday),wp)) < zprec )   fjulday = REAL(NINT(fjulday),wp)   ! avoid truncation error
      IF( ABS(adatrj  - REAL(NINT(adatrj ),wp)) < zprec )   adatrj  = REAL(NINT(adatrj ),wp)   ! avoid truncation error

      IF( nsec_day > nsecd ) THEN                       ! New day
         !
         nday      = nday + 1
         nday_year = nday_year + 1
         nsec_day  = ndt05
         !
         IF( nday == nmonth_len(nmonth) + 1 ) THEN      ! New month
            nday   = 1
            nmonth = nmonth + 1
            nsec_month = ndt05
            IF( nmonth == 13 ) THEN                     ! New year
               nyear     = nyear + 1
               nmonth    = 1
               nday_year = 1
               nsec_year = ndt05
               nsec1jan000 = nsec1jan000 + nsecd * nyear_len(1)
               IF( nleapy == 1 )   CALL day_mth
            ENDIF
         ENDIF
         !
         ndastp = nyear * 10000 + nmonth * 100 + nday   ! New date
         !
         !compute first day of the year in julian days
         CALL ymds2ju( nyear, 01, 01, 0.0, fjulstartyear )
         !
         IF(lwp) WRITE(numout,'(a,i8,a,i4.4,a,i2.2,a,i2.2,a,i3.3)') '======>> time-step =', kt,   &
              &   '      New day, DATE Y/M/D = ', nyear, '/', nmonth, '/', nday, '      nday_year = ', nday_year
         IF(lwp) WRITE(numout,'(a,i8,a,i7,a,i5)') '         nsec_year = ', nsec_year,   &
              &   '   nsec_month = ', nsec_month, '   nsec_day = ', nsec_day, '   nsec_week = ', nsec_week
      ENDIF

      IF( nsec_week > 7*nsecd )   nsec_week = ndt05     ! New week

      IF(ln_ctl) THEN
         WRITE(charout,FMT="('kt =', I4,'  d/m/y =',I2,I2,I4)") kt, nday, nmonth, nyear
         CALL prt_ctl_info(charout)
      ENDIF
      !
   END SUBROUTINE day_tan

   SUBROUTINE day_adj( kt )
      !!----------------------------------------------------------------------
      !!                      ***  ROUTINE day_adj  ***
      !!
      !! ** Purpose :   Compute the date with a day iteration backward
      !!                if necessary.
      !!
      !! ** Method  : - ???
      !!
      !! ** Action  : - nyear     : current year
      !!              - nmonth    : current month of the year nyear
      !!              - nday      : current day of the month nmonth
      !!              - nday_year : current day of the year nyear
      !!              - ndastp    : = nyear*10000 + nmonth*100 + nday
      !!              - adatrj    : date in days since the beginning of the run
      !!              - nsec_year : current time of the year (in second since 00h, jan 1st)
      !!----------------------------------------------------------------------
      INTEGER, INTENT(in) ::   kt        ! ocean time-step indices
      !
      CHARACTER (len=25) ::   charout
      REAL(wp)           ::   zprec      ! fraction of day corresponding to 0.1 second
      !!----------------------------------------------------------------------
      zprec = 0.1 / rday
      !                                                 ! New time-step
      nsec_year  = nsec_year  - ndt
      nsec_month = nsec_month - ndt
      nsec_day   = nsec_day   - ndt
      adatrj = adatrj + rdttra(1) / rday
      fjulday = fjulday + rdttra(1) / rday
      IF( ABS(fjulday - REAL(NINT(fjulday),wp)) < zprec )   fjulday = REAL(NINT(fjulday),wp)   ! avoid truncation error
      IF( ABS(adatrj  - REAL(NINT(adatrj ),wp)) < zprec )   adatrj  = REAL(NINT(adatrj ),wp)   ! avoid truncation error

      IF( nsec_day < 0 ) THEN                        ! NEW day
         !
         nday      = nday - 1
         nday_year = nday_year - 1
         nsec_day  = rday - ndt05
         !
         IF( nday == 0 ) THEN      ! NEW month
            nmonth = nmonth - 1
            IF( nmonth == 0 ) THEN                     ! NEW year
               nyear     = nyear - 1
               nmonth    = 12
               nday_year = nyear_len(0)
               nsec_year = nday_year * rday - ndt05
               nsec1jan000 = nsec1jan000 - nsecd * nyear_len(0)
               IF( nleapy == 1 )   CALL day_mth
            ENDIF
            nday   = nmonth_len(nmonth)
            nsec_month = nmonth_len(nmonth) * rday - ndt05
         ENDIF
         !
         ndastp = nyear * 10000 + nmonth * 100 + nday   ! NEW date
         !
         IF(lwp) WRITE(numout,'(a,i8,a,i4.4,a,i2.2,a,i2.2,a,i3.3)') '======>> time-step =', kt,   &
              &   '      New day, DATE Y/M/D = ', nyear, '/', nmonth, '/', nday, '      nday_year = ', nday_year
         IF(lwp) WRITE(numout,'(a,i8,a,i7,a,i5)') '         nsec_year = ', nsec_year,   &
              &   '   nsec_month = ', nsec_month, '   nsec_day = ', nsec_day
      ENDIF

      !
   END SUBROUTINE day_adj

   !!======================================================================
END MODULE daymod_tam
