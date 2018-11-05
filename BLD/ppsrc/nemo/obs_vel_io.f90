MODULE obs_vel_io
   !!======================================================================
   !!                       ***  MODULE obs_vel_io  ***
   !! Observation operators : I/O for TAO velocity data
   !!======================================================================
   !! History : 
   !!             !  09-01  (K. Mogensen) Initial version based on obs_read_taovel
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   read_taofile    :  Read a obfbdata structure from an TAO velocity file.
   !!----------------------------------------------------------------------
   USE par_kind
   USE obs_utils
   USE obs_fbm
   USE obs_conv
   USE in_out_manager
   USE netcdf
   IMPLICIT NONE

   INTEGER, PARAMETER :: imaxlev = 10000

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obs_vel_io.F90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

CONTAINS

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: obsvel_io.h90 2287 2010-10-18 07:53:52Z smasson $
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------

   SUBROUTINE read_taondbc( cdfilename, inpfile, kunit, ldwp, ldgrid )
      !!---------------------------------------------------------------------
      !!
      !!                     ** ROUTINE read_enactfile **
      !!
      !! ** Purpose : Read from file the TAO data fro NDBC.
      !!
      !! ** Method  : The data file is a NetCDF file. 
      !!
      !! ** Action  :
      !!
      !! ** Reference : http://tao.noaa.gov/tao/data_deliv/deliv_ndbc.shtml
      !! History : 
      !!          ! 09-01 (K. Mogensen) Original version.
      !!----------------------------------------------------------------------
      !! * Arguments
      CHARACTER(LEN=*) :: cdfilename ! Input filename
      TYPE(obfbdata)   :: inpfile    ! Output obfbdata structure
      INTEGER          :: kunit      ! Unit for output
      LOGICAL          :: ldwp       ! Print info
      LOGICAL          :: ldgrid     ! Save grid info in data structure
      !! * Local declarations
      INTEGER :: iobs                ! Number of observations
      INTEGER :: ilev                ! Number of levels
      INTEGER :: ilat                ! Number of latitudes
      INTEGER :: ilon                ! Number of longtudes
      INTEGER :: itim                ! Number of obs. times
      INTEGER :: i_file_id
      INTEGER :: i_dimid_id
      INTEGER :: i_phi_id
      INTEGER :: i_lam_id
      INTEGER :: i_depth_id
      INTEGER :: i_var_id
      INTEGER :: i_time_id
      INTEGER :: i_time2_id
      INTEGER :: i_qc_var_id
      CHARACTER(LEN=40) :: cl_fld_lam
      CHARACTER(LEN=40) :: cl_fld_phi
      CHARACTER(LEN=40) :: cl_fld_depth
      CHARACTER(LEN=40) :: cl_fld_var_u
      CHARACTER(LEN=40) :: cl_fld_var_v
      CHARACTER(LEN=40) :: cl_fld_var_qc_uv1
      CHARACTER(LEN=40) :: cl_fld_var_qc_uv2
      CHARACTER(LEN=40) :: cl_fld_time
      CHARACTER(LEN=40) :: cl_fld_time2
      INTEGER :: ja
      INTEGER :: jo
      INTEGER :: jk
      INTEGER :: jt
      REAL(wp), ALLOCATABLE, DIMENSION(:,:,:,:) :: &
         & zv,     &
         & zu,     &
         & zuv1qc, &
         & zuv2qc
      REAL(wp), ALLOCATABLE, DIMENSION(:) :: &
         & zdep, &
         & zlat, &
         & zlon, &
         & zjuld
      REAL(wp) :: zl
      INTEGER, ALLOCATABLE, DIMENSION(:) :: &
         & itime, &
         & itime2
      CHARACTER(LEN=50) :: cdjulref
      CHARACTER(LEN=12), PARAMETER :: cl_name = 'read_taondbc'
      CHARACTER(len=1) :: cns, cew

      !-----------------------------------------------------------------------
      ! Initialization
      !-----------------------------------------------------------------------
      cl_fld_lam                 = 'lon'
      cl_fld_phi                 = 'lat'
      cl_fld_depth               = 'depth'
      cl_fld_time                = 'time'
      cl_fld_time2               = 'time2'

      !-----------------------------------------------------------------------
      ! Open file
      !-----------------------------------------------------------------------

      CALL chkerr( nf90_open( TRIM( cdfilename ), nf90_nowrite, &
            &      i_file_id ),           cl_name, 88 )

      !-----------------------------------------------------------------------
      ! Read the heading of the file
      !-----------------------------------------------------------------------
      IF(ldwp) WRITE(kunit,*)
      IF(ldwp) WRITE(kunit,*) ' read_taondbc :' 
      IF(ldwp) WRITE(kunit,*) ' ~~~~~~~~~~~~'
      
      !---------------------------------------------------------------------
      ! Read the number of observations and of levels to allocate array
      !---------------------------------------------------------------------
      CALL chkerr( nf90_inq_dimid        ( i_file_id, 'time', i_dimid_id ),        &
         &         cl_name, 101 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_dimid_id, len = itim ),    &
         &         cl_name, 103 )
      CALL chkerr( nf90_inq_dimid        ( i_file_id, 'depth', i_dimid_id ),       &
         &         cl_name, 105 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_dimid_id, len = ilev ),    &
         &         cl_name, 107 )
      CALL chkerr( nf90_inq_dimid        ( i_file_id, 'lat', i_dimid_id ),           &
         &         cl_name, 109 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_dimid_id, len = ilat ),    &
         &         cl_name, 111 )
      CALL chkerr( nf90_inq_dimid        ( i_file_id, 'lon', i_dimid_id ),         &
         &         cl_name, 113 )
      CALL chkerr( nf90_inquire_dimension( i_file_id, i_dimid_id, len = ilon ),    &
         &         cl_name, 115 ) 

      iobs = itim * ilat * ilon
      IF(ldwp)WRITE(kunit,*) '         No. of data records = ', iobs
      IF(ldwp)WRITE(kunit,*) '         No. of levels       = ', ilev
      IF(ldwp)WRITE(kunit,*) 

      !---------------------------------------------------------------------
      ! Allocate arrays
      !---------------------------------------------------------------------

      CALL init_obfbdata( inpfile )
      CALL alloc_obfbdata( inpfile, 2, iobs, ilev, 0, 0, ldgrid )
      inpfile%cname(1) = 'UVEL'
      inpfile%cname(2) = 'VVEL'
      inpfile%coblong(1) = 'Zonal current'
      inpfile%coblong(2) = 'Meridional current'
      inpfile%cobunit(1) = 'Meters per second'
      inpfile%cobunit(2) = 'Meters per second'

      ALLOCATE( &
         & zu(ilon,ilat,ilev,itim),     &
         & zv(ilon,ilat,ilev,itim),     &
         & zdep(ilev),                  &
         & zuv1qc(ilon,ilat,ilev,itim), &
         & zuv2qc(ilon,ilat,ilev,itim), &
         & itime(itim),                 &
         & itime2(itim),                &
         & zlat(ilat),                  &
         & zlon(ilon),                  &
         & zjuld(itim)                  &
         & )

      !---------------------------------------------------------------------
      ! Read the time/position variables 
      !---------------------------------------------------------------------
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_time, i_time_id ),                               &
         &         cl_name, 153 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_time_id, itime ),                                     &
         &         cl_name, 155 )

      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_time2, i_time2_id ),                             &
         &         cl_name, 158 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_time2_id, itime2 ),                                   &
         &         cl_name, 160 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_depth, i_depth_id ),                             &
            &         cl_name, 163 )         
      CALL chkerr( nf90_get_var  ( i_file_id, i_depth_id, zdep ),                                     &
         &         cl_name, 165 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_phi, i_phi_id ),                                 &
         &         cl_name, 168 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_phi_id, zlat ),                                       &
         &         cl_name, 170 )
      
      CALL chkerr( nf90_inq_varid( i_file_id, cl_fld_lam, i_lam_id ),                                 &
         &         cl_name, 173 )
      CALL chkerr( nf90_get_var  ( i_file_id, i_lam_id, zlon ),                                       &
         &         cl_name, 175 )
      
      !---------------------------------------------------------------------
      ! Read the variables
      !---------------------------------------------------------------------

      ! ADCP format assumed
      cl_fld_var_u = 'u_1205'
      IF ( nf90_inq_varid( i_file_id, cl_fld_var_u, i_var_id ) /= nf90_noerr ) THEN
         ! Try again with current meter format
         cl_fld_var_u = 'U_320'
         IF ( nf90_inq_varid( i_file_id, cl_fld_var_u, i_var_id ) /= nf90_noerr ) THEN
            CALL fatal_error( 'Unknown format in read_taondbc', 187 )
         ENDIF
      ENDIF
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, zu ),                                         &
         &         cl_name, 191 )
      
      ! ADCP format assumed
      cl_fld_var_v = 'v_1206'
      IF ( nf90_inq_varid( i_file_id, cl_fld_var_v, i_var_id ) /= nf90_noerr ) THEN
         ! Try again with current meter format
         cl_fld_var_v = 'V_321'
         IF ( nf90_inq_varid( i_file_id, cl_fld_var_v, i_var_id ) /= nf90_noerr ) THEN
            CALL fatal_error( 'Unknown format in read_taondbc', 199 )
         ENDIF
      ENDIF
      CALL chkerr( nf90_get_var  ( i_file_id, i_var_id, zv ),                                         &
         &         cl_name, 203 )

      !---------------------------------------------------------------------
      ! Read the QC attributes
      !---------------------------------------------------------------------
      
      ! ADCP format assumed
      cl_fld_var_qc_uv1 = 'QU_5205'
      IF ( nf90_inq_varid( i_file_id, cl_fld_var_qc_uv1, i_qc_var_id ) /= nf90_noerr ) THEN
         ! Try again with current meter format
         cl_fld_var_qc_uv1 = 'QCS_5300'
         IF ( nf90_inq_varid( i_file_id, cl_fld_var_qc_uv1, i_qc_var_id ) /= nf90_noerr ) THEN
            ! Try again with high freq. current meter format
            cl_fld_var_qc_uv1 = 'QCU_5320'
            IF ( nf90_inq_varid( i_file_id, cl_fld_var_qc_uv1, i_qc_var_id ) /= nf90_noerr ) THEN
               CALL fatal_error( 'Unknown format in read_taondbc', 218 )
            ENDIF
         ENDIF
      ENDIF
      CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, zuv1qc),                                   &
         &         cl_name, 223 )

      ! ADCP format assumed
      cl_fld_var_qc_uv2 = 'QV_5206'
      IF ( nf90_inq_varid( i_file_id, cl_fld_var_qc_uv2, i_qc_var_id ) /= nf90_noerr ) THEN
         ! Try again with current meter format
         cl_fld_var_qc_uv2 = 'QCD_5310'
         IF ( nf90_inq_varid( i_file_id, cl_fld_var_qc_uv2, i_qc_var_id ) /= nf90_noerr ) THEN
            ! Try again with high freq. current meter format
            cl_fld_var_qc_uv2 = 'QCV_5321'
            IF ( nf90_inq_varid( i_file_id, cl_fld_var_qc_uv2, i_qc_var_id ) /= nf90_noerr ) THEN
               CALL fatal_error( 'Unknown format in read_taondbc', 234 )
            ENDIF
         ENDIF
      ENDIF
      CALL chkerr( nf90_get_var  ( i_file_id, i_qc_var_id, zuv2qc),                                   &
         &         cl_name, 239 )

      !---------------------------------------------------------------------
      ! Close file
      !---------------------------------------------------------------------

      CALL chkerr( nf90_close( i_file_id ),           cl_name, 245 )

      !---------------------------------------------------------------------
      ! Convert to to 19500101 based Julian date
      !---------------------------------------------------------------------
      DO jt = 1, itim
         zjuld(jt) = REAL(itime(jt),wp) + REAL(itime2(jt),wp)/86400000.0_wp &
            &           - 2433283.0_wp
      END DO
      inpfile%cdjuldref = '19500101000000'

      !---------------------------------------------------------------------
      ! Copy info to obfbdata structure
      !---------------------------------------------------------------------

      iobs = 0
      DO jt = 1, itim
         DO ja = 1, ilat
            DO jo = 1, ilon
               iobs = iobs + 1
               zl = zlon(jo)
               IF ( zl > 180.0_wp ) zl = zl - 360.0_wp
               IF ( zl < 0 ) THEN
                  cew = 'w'
               ELSE
                  cew = 'e'
               ENDIF
               IF ( zlat(jo) < 0 ) THEN
                  cns = 's'
               ELSE
                  cns = 'n'
               ENDIF
               WRITE(inpfile%cdwmo(iobs),'(A1,I2.2,A1,I3.3)') &
                  & cns, ABS(NINT(zlat(ja))), cew, ABS(NINT(zl))
               DO jk = 1, ilev
                  inpfile%pob(jk,iobs,1)     = zu(jo,ja,jk,jt)
                  inpfile%pob(jk,iobs,2)     = zv(jo,ja,jk,jt)
                  inpfile%pdep(jk,iobs)      = zdep(jk)
                  inpfile%ivlqc(jk,iobs,1:2) = INT( MAX( zuv1qc(jo,ja,jk,jt), zuv2qc(jo,ja,jk,jt) ) )
               END DO
               inpfile%plam(iobs) = zlon(jo)
               inpfile%pphi(iobs) = zlat(ja)
               inpfile%ptim(iobs) = zjuld(jt)
            END DO
         END DO
      END DO

      ! No position, time, depth and variable QC in input files
      DO jo = 1, iobs
         inpfile%ipqc(jo) = 1
         inpfile%itqc(jo) = 1
         inpfile%ivqc(jo,1:2) = 1
         DO jk = 1, ilev
            inpfile%idqc(jk,jo) = 1
         END DO
      END DO

      !---------------------------------------------------------------------
      ! Set the platform information
      !---------------------------------------------------------------------
      inpfile%cdtyp(:)=' 820'

      !---------------------------------------------------------------------
      ! Set QC flags for missing data and rescale to m/s
      !---------------------------------------------------------------------

      DO jo = 1, iobs
         DO jk = 1, ilev
            IF ( ( ABS(inpfile%pob(jk,jo,1)) > 10000.0_wp ) .OR. &
               & ( ABS(inpfile%pob(jk,jo,2)) > 10000.0_wp ) ) THEN
               inpfile%ivlqc(jk,jo,:) = 4
               inpfile%pob(jk,jo,1) = fbrmdi
               inpfile%pob(jk,jo,2) = fbrmdi
            ELSE
               inpfile%pob(jk,jo,1) = 0.01 * inpfile%pob(jk,jo,1)
               inpfile%pob(jk,jo,2) = 0.01 * inpfile%pob(jk,jo,2)
            ENDIF
         END DO
      END DO

      !---------------------------------------------------------------------
      ! Set file indexes
      !---------------------------------------------------------------------

      DO jo = 1, inpfile%nobs
         inpfile%kindex(jo) = jo
      END DO

      !---------------------------------------------------------------------
      ! Initialize flags since they are not in the TAO input files
      !---------------------------------------------------------------------

      inpfile%ioqcf(:,:)      = 0
      inpfile%ipqcf(:,:)      = 0
      inpfile%itqcf(:,:)      = 0
      inpfile%idqcf(:,:,:)    = 0
      inpfile%ivqcf(:,:,:)    = 0
      inpfile%ivlqcf(:,:,:,:) = 0

      !---------------------------------------------------------------------
      ! Deallocate data
      !---------------------------------------------------------------------
      DEALLOCATE( &
         & zu,     &
         & zv,     &
         & zdep,   &
         & zuv1qc, &
         & zuv2qc, &
         & itime,  &
         & itime2, &
         & zlat,   &
         & zlon,   &
         & zjuld   &
         & )

   END SUBROUTINE read_taondbc

END MODULE obs_vel_io
