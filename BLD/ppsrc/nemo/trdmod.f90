MODULE trdmod
   !!======================================================================
   !!                       ***  MODULE  trdmod  ***
   !! Ocean diagnostics:  ocean tracers and dynamic trends
   !!=====================================================================
   !! History :  1.0  !  2004-08  (C. Talandier) Original code
   !!             -   !  2005-04  (C. Deltel)    Add Asselin trend in the ML budget
   !!            3.3  ! 2010-10  (C. Ethe, G. Madec) reorganisation of initialisation phase
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Default case :                                         Empty module
   !!----------------------------------------------------------------------
   USE trdmod_oce      ! ocean variables trends
   USE trdvor          ! ocean vorticity trends 
   USE trdicp          ! ocean bassin integral constraints properties
   USE trdmld          ! ocean active mixed layer tracers trends 
   !!----------------------------------------------------------------------
CONTAINS
   SUBROUTINE trd_mod(ptrd3dx, ptrd3dy, ktrd , ctype, kt )   ! Empty routine
      REAL(wp) ::   ptrd3dx(:,:,:), ptrd3dy(:,:,:)
      INTEGER  ::   ktrd, kt                            
      CHARACTER(len=3) ::  ctype                  
      WRITE(*,*) 'trd_3d: You should not have seen this print! error ?', ptrd3dx(1,1,1), ptrd3dy(1,1,1)
      WRITE(*,*) ' "   ": You should not have seen this print! error ?', ktrd, ctype, kt
   END SUBROUTINE trd_mod

   SUBROUTINE trd_mod_init
      !!----------------------------------------------------------------------
      !!                  ***  ROUTINE trd_mod_init  ***
      !! 
      !! ** Purpose :   Initialization of activated trends
      !!----------------------------------------------------------------------
      USE in_out_manager          ! I/O manager
      !!    
      NAMELIST/namtrd/ nn_trd, nn_ctls, cn_trdrst_in, cn_trdrst_out, ln_trdmld_restart, rn_ucf, ln_trdmld_instant
      !!----------------------------------------------------------------------

      IF( l_trdtra .OR. l_trddyn )   THEN
         REWIND( numnam )
         READ  ( numnam, namtrd )      ! namelist namtrd : trends diagnostic

         IF(lwp) THEN
            WRITE(numout,*)
            WRITE(numout,*) ' trd_mod_init : Momentum/Tracers trends'
            WRITE(numout,*) ' ~~~~~~~~~~~~~'
            WRITE(numout,*) '   Namelist namtrd : set trends parameters'
            WRITE(numout,*) '      frequency of trends diagnostics   nn_trd             = ', nn_trd
            WRITE(numout,*) '      control surface type              nn_ctls            = ', nn_ctls
            WRITE(numout,*) '      restart for ML diagnostics        ln_trdmld_restart  = ', ln_trdmld_restart
            WRITE(numout,*) '      instantaneous or mean ML T/S      ln_trdmld_instant  = ', ln_trdmld_instant
            WRITE(numout,*) '      unit conversion factor            rn_ucf             = ', rn_ucf
        ENDIF
      ENDIF
      !
      IF( lk_trddyn .OR. lk_trdtra )    CALL trd_icp_init       ! integral constraints trends
      IF( lk_trdmld                )    CALL trd_mld_init       ! mixed-layer trends (active  tracers)  
      IF( lk_trdvor                )    CALL trd_vor_init       ! vorticity trends        
      !
   END SUBROUTINE trd_mod_init

   !!======================================================================
END MODULE trdmod
