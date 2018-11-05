MODULE par_oce
   !!======================================================================
   !!                        ***  par_oce  ***
   !! Ocean :   set the ocean parameters
   !!======================================================================
   !! History :  OPA  !  1991     (Imbard, Levy, Madec)  Original code
   !!   NEMO     1.0  !  2004-01  (G. Madec, J.-M. Molines)  Free form and module
   !!            3.3  !  2010-09  (C. Ethe) TRA-TRC merge: add jpts, jp_tem & jp_sal
   !!----------------------------------------------------------------------
   USE par_kind          ! kind parameters

   IMPLICIT NONE
   PUBLIC

   !!----------------------------------------------------------------------
   !!   Domain decomposition
   !!----------------------------------------------------------------------
   !! if we dont use massively parallel computer (parameters jpni=jpnj=1) so jpiglo=jpi and jpjglo=jpj
   INTEGER, PUBLIC            ::   jpni         !: number of processors following i 
   INTEGER, PUBLIC            ::   jpnj         !: number of processors following j
   INTEGER, PUBLIC            ::   jpnij        !: nb of local domain = nb of processors ( <= jpni x jpnj )
   INTEGER, PUBLIC, PARAMETER ::   jpr2di = 0   !: number of columns for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   jpr2dj = 0   !: number of rows    for extra outer halo 
   INTEGER, PUBLIC, PARAMETER ::   jpreci = 1   !: number of columns for overlap 
   INTEGER, PUBLIC, PARAMETER ::   jprecj = 1   !: number of rows    for overlap 

   !! Ocean Domain sizes
   !! ------------------
   !!   data           domain   (jpidta,jpjdta)
   !!   global or zoom domain   (jpiglo,jpjglo)
   !!   local          domain   ( jpi  , jpj  )
   
   !!---------------------------------------------------------------------
   !!   'key_orca_r2'   :                           global ocean : ORCA R4
   !!---------------------------------------------------------------------
   !!---------------------------------------------------------------------
   !!                     ***  par_ORCA_R2.h90  ***  
   !!   Ocean Domain : 2 degrees resolution global ocean
   !!                  (0RCA_R2 configuration)
   !!---------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: par_ORCA_R2.h90 2715 2011-03-30 15:58:35Z rblod $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!----------------------------------------------------------------------
   CHARACTER (len=16)      &
      , PARAMETER  &
      ::    &  
      cp_cfg = "orca"           !: name of the configuration 
   INTEGER     &
      , PARAMETER  &
      :: &
      jp_cfg = 2,            &  !: resolution of the configuration (degrees)

      ! data size              !!! * size of all input files *
      jpidta  = 182,         &  !: 1st lateral dimension ( >= jpiglo )
      jpjdta  = 149,         &  !: 2nd    "       "      ( >= jpjglo )
      jpkdta  = 31              !: number of levels      ( >= jpk    ) 

      ! global domain size     !!! *  global domain  *
   INTEGER    &
      , PARAMETER  &
      :: &
      jpiglo  = jpidta,      &  !: 1st dimension of global domain --> i
      jpjglo  = jpjdta,      &  !: 2nd    "                  "    --> j
      ! starting position of the zoom 
      jpizoom =   1   ,      &  !: left bottom (i,j) indices of the zoom
      jpjzoom =   1   ,      &  !: in data domain indices
      ! Domain characteristics
      jperio  =   4             !: lateral cond. type (between 0 and 6)


   !!  Values set to pp_not_used indicates that this parameter is not used in THIS config.
   !!  Values set to pp_to_be_computed  indicates that variables will be computed in domzgr
   REAL(wp), PARAMETER ::   &
      pp_not_used       = 999999_wp , &  !:
      pp_to_be_computed = 0._wp          !:

   !! Coefficients associated with the horizontal coordinate system (jphgr_msh /= 0 )

   INTEGER,PARAMETER   ::    & !
      jphgr_msh = 0            !: type of horizontal mesh
      !                        !  = 0 curvilinear coordinate on the sphere
      !                        !      read in coordinate.nc file
      !                        !  = 1 geographical mesh on the sphere
      !                        !      with regular grid-spacing
      !                        !  = 2 f-plane with regular grid-spacing
      !                        !  = 3 beta-plane with regular grid-spacing
      !                        !  = 4 Mercator grid with T/U point at the equator  with
      !                        !      isotropic resolution (e1_deg)

      !   ppglam0 , ppgphi0: coordinates of the lower leftmost T point of the grid.
      !   The mercator grid starts only approximately at gphi0 because
      !   of the constraint that the equator be a T point.
   REAL(wp) ,PARAMETER ::       &  !
      ppglam0  = pp_not_used,   &  !: longitude of first raw and column T-point (jphgr_msh = 1)
      ppgphi0  = pp_not_used,   &  !: latitude  of first raw and column T-point (jphgr_msh = 1)
      !                            !  latitude for the Coriolis or Beta parameter (jphgr_msh = 2 or 3)
      ppe1_deg = pp_not_used,   &  !: zonal      grid-spacing (degrees)
      ppe2_deg = pp_not_used,   &  !: meridional grid-spacing (degrees)
      !
      ppe1_m   = pp_not_used,   &  !: zonal      grid-spacing (meters )
      ppe2_m   = pp_not_used       !: meridional grid-spacing (meters )

   !!
   !! Vertical grid parameter for domzgr
   !! ==================================
   !!
   REAL(wp), PARAMETER  ::       &
      &     ppsur = -4762.96143546300_wp    ,  &  !: ORCA r4, r2 and r05 coefficients
      &     ppa0  =   255.58049070440_wp    ,  &  !: (default coefficients)
      &     ppa1  =   245.58132232490_wp    ,  &  !:
      &     ppkth =    21.43336197938_wp    ,  &  !: (non dimensional): gives the approximate
      !                                           !: layer number above which  stretching will
      !                                           !: be maximum. Usually of order jpk/2.
      &     ppacr =     3.00000000000_wp          !: (non dimensional): stretching factor
      !                                           !: for the grid. The highest zacr, the smallest
      !                                           !: the stretching.

   !!
   !!  If both ppa0 ppa1 and ppsur are specified to 0, then
   !!  they are computed from ppdzmin, pphmax , ppkth, ppacr in dom_zgr
   !!
   REAL(wp), PARAMETER ::        &
      &     ppdzmin = pp_not_used           ,  &  !: (meters) vertical thickness of the top layer
      &     pphmax  = pp_not_used                 !: (meters) Maximum depth of the ocean gdepw(jpk)
   LOGICAL,  PARAMETER ::        &
      &     ldbletanh = .FALSE.                   !: Use/do not use double tanf function for vertical coordinates
   REAL(wp), PARAMETER ::        &
      &     ppa2    = pp_not_used           ,  &  !: Double tanh function parameters
      &     ppkth2  = pp_not_used           ,  &  !:
      &     ppacr2  = pp_not_used                 !:
   !!---------------------------------------------------------------------


   !!---------------------------------------------------------------------
   !! Active tracer parameters
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC, PARAMETER ::   jpts   = 2    !: Number of active tracers (=2, i.e. T & S )
   INTEGER, PUBLIC, PARAMETER ::   jp_tem = 1    !: indice for temperature
   INTEGER, PUBLIC, PARAMETER ::   jp_sal = 2    !: indice for salinity

   !!---------------------------------------------------------------------
   !! Domain Matrix size  (if AGRIF, they are not all parameters)
   !!---------------------------------------------------------------------
   INTEGER, PUBLIC  ::   jpi   ! = ( jpiglo-2*jpreci + (jpni-1) ) / jpni + 2*jpreci   !: first  dimension
   INTEGER, PUBLIC  ::   jpj   ! = ( jpjglo-2*jprecj + (jpnj-1) ) / jpnj + 2*jprecj   !: second dimension
   INTEGER, PUBLIC  ::   jpk   ! = jpkdta
   INTEGER, PUBLIC  ::   jpim1 ! = jpi-1                                            !: inner domain indices
   INTEGER, PUBLIC  ::   jpjm1 ! = jpj-1                                            !:   -     -      -
   INTEGER, PUBLIC  ::   jpkm1 ! = jpk-1                                            !:   -     -      -
   INTEGER, PUBLIC  ::   jpij  ! = jpi*jpj                                          !:  jpi x jpj

   !!---------------------------------------------------------------------
   !! Optimization/control flags
   !!---------------------------------------------------------------------
   LOGICAL, PUBLIC, PARAMETER ::   lk_esopa     = .FALSE.  !: flag to activate the all options

   LOGICAL, PUBLIC, PARAMETER ::   lk_vopt_loop = .FALSE.  !: vector optimization flag

   !!----------------------------------------------------------------------
   !! NEMO/OPA 3.3 , NEMO Consortium (2010)
   !! $Id: par_oce.F90 5196 2015-04-07 08:28:07Z pabouttier $ 
   !! Software governed by the CeCILL licence (NEMOGCM/NEMO_CeCILL.txt)
   !!======================================================================
END MODULE par_oce
