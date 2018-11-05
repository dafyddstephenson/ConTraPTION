MODULE sbccpl
   !!======================================================================
   !!                       ***  MODULE  sbccpl  ***
   !! Surface Boundary Condition :  momentum, heat and freshwater fluxes in coupled mode
   !!======================================================================
   !! History :  2.0  ! 2007-06  (R. Redler, N. Keenlyside, W. Park) Original code split into flxmod & taumod
   !!            3.0  ! 2008-02  (G. Madec, C Talandier)  surface module
   !!            3.1  ! 2009_02  (G. Madec, S. Masson, E. Maisonave, A. Caubel) generic coupled interface
   !!            3.4  ! 2011_11  (C. Harris) more flexibility + multi-category fields
   !!----------------------------------------------------------------------
   !!----------------------------------------------------------------------
   !!   Dummy module                                            NO coupling
   !!----------------------------------------------------------------------
   USE par_kind        ! kind definition
CONTAINS
   SUBROUTINE sbc_cpl_snd( kt )
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?', kt
   END SUBROUTINE sbc_cpl_snd
   !
   SUBROUTINE sbc_cpl_rcv( kt, k_fsbc, k_ice )     
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?', kt, k_fsbc, k_ice
   END SUBROUTINE sbc_cpl_rcv
   !
   SUBROUTINE sbc_cpl_ice_tau( p_taui, p_tauj )     
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_taui   ! i- & j-components of atmos-ice stress [N/m2]
      REAL(wp), INTENT(out), DIMENSION(:,:) ::   p_tauj   ! at I-point (B-grid) or U & V-point (C-grid)
      p_taui(:,:) = 0.   ;   p_tauj(:,:) = 0. ! stupid definition to avoid warning message when compiling...
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?'
   END SUBROUTINE sbc_cpl_ice_tau
   !
   SUBROUTINE sbc_cpl_ice_flx( p_frld , palbi   , psst    , pist  )
      REAL(wp), INTENT(in   ), DIMENSION(:,:  ) ::   p_frld     ! lead fraction                [0 to 1]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   palbi   ! ice albedo
      REAL(wp), INTENT(in   ), DIMENSION(:,:  ), OPTIONAL ::   psst    ! sea surface temperature      [Celcius]
      REAL(wp), INTENT(in   ), DIMENSION(:,:,:), OPTIONAL ::   pist    ! ice surface temperature      [Kelvin]
      WRITE(*,*) 'sbc_cpl_snd: You should not have seen this print! error?', p_frld(1,1), palbi(1,1,1), psst(1,1), pist(1,1,1) 
   END SUBROUTINE sbc_cpl_ice_flx
   

   !!======================================================================
END MODULE sbccpl
