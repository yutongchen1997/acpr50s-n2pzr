! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: pvessel (pressurizer vessel).
! 
! ------------------------------------------------------------------
module pvessel_mod
  !======= Inclusions =========
  use ctlsgnl_mod

  !======= Declarations =========
  implicit none

  type :: pvessel
    !=== Member Variables === 
    ! required input parameters ------------------------------------
    character(len=20) :: id
    ! id string of liq and gas section, material of liq and gas
    character(len=20) :: id_l, id_g, mat_l, mat_g
    ! inner diameter, total vessel height
    real(c_double) :: dia_i, len_t
    ! other parameters ---------------------------------------------
    ! cross-sectional area, total vessel volume
    real(c_double) :: csa_z, vol_t 
    ! boundary conditions ------------------------------------------
    ! mass flow rates: gas inlet, gas exit, surge line
    real(c_double) :: wg_in, wg_ex, w_su
    ! enthalpy flow rates: gas inlet, gas exit, surge line
    real(c_double) :: wghg_in, wghg_ex, wh_su
    ! wall heat transfer rates
    real(c_double) :: qw_g, qw_l
    ! dynamic variables --------------------------------------------
    ! mass, enthalpy, temperature, volume, pressure, level
    real(c_double) :: m_g, m_l  ! gas, liq mass         - [kg]
    real(c_double) :: h_g, h_l  ! gas, liq enthalpy     - [J]
    real(c_double) :: t_g, t_l  ! gas, liq temperature  - [K]
    real(c_double) :: v_g, v_l  ! gas, liq volume       - [m^3]
    real(c_double) :: p_g       ! gas pressure          - [Pa]
    real(c_double) :: z_l       ! internal liquid level - [m]
    ! derivatives of dynamic variables
    real(c_double) :: dmgdt, dmldt   ! time derivative of m  - [kg/s]
    real(c_double) :: dhgdt, dhldt   ! time derivative of h  - [J/s]
    real(c_double) :: dtgdt, dtldt   ! time derivative of t  - [K/s]
    real(c_double) :: dvgdt, dvldt   ! time derivative of v  - [m^3/s]
    real(c_double) :: dpdt           ! time derivative of pg - [Pa/s]
    ! internal variables -------------------------------------------
    real(c_double) :: htc_i   ! liq - gas interface htc - [W/m^2/K]  
    real(c_double) :: qlg_i   ! liq - gas heating rate  - [W]
  contains  
    !=== Member Procedures === 
    procedure :: init_pvessel
    procedure :: calc_res_pvessel

  end type pvessel

contains

  subroutine init_pvessel(self, id_l, id_g, mat_l, mat_g, inp)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(pvessel), intent(inout) :: self
    character(len=20), intent(in) :: id_l, id_g, mat_l, mat_g
    ! inp array: 1-dia_i, 2-len_t, 3-len_i, 4-tg0, 5-tl0, 6-pg0, 7-htc_0
    real(c_double), intent(in) :: inp(:) 
    !=== local variables ===
    real(c_double) :: rg0, rl0, s_hg0, s_hl0
    !======= Internals ============
    ! assign the input strings
    self%id_l = id_l
    self%id_g = id_g
    self%mat_l = mat_l
    self%mat_g = mat_g
    ! assign the input parameters
    self%dia_i = inp(1)
    self%len_t = inp(2)
    self%z_l   = inp(3)
    self%t_g   = inp(4)
    self%t_l   = inp(5)
    self%p_g   = inp(6)
    self%htc_i = inp(7)
    ! calculate cross-sectional area
    self%csa_z = pi*inp(1)**2.d0/4.d0
    ! calculate total volume
    self%vol_t = self%csa_z*inp(2)
    ! calculate initial liquid and gas volume
    self%v_l = self%csa_z*self%z_l
    self%v_g = self%vol_t-self%v_l
    ! calcualte initial liquid and gas specific enthalpy [J/kg]
    call hgeos_n2g(eosopt, self%p_g, self%t_g, s_hg0)
    call hleos_h2o(self%p_g, self%t_l, s_hl0)        
    ! calculate initial liquid and gas density [kg/m^3]
    call rgeos_n2g(eosopt, self%p_g, s_hg0, rg0)
    call rleos_h2o(self%p_g, s_hl0, rl0)
    ! calculate initial liquid and gas mass [kg]
    self%m_l = rl0*self%v_l
    self%m_g = rg0*self%v_g
    ! calculate initial liquid and gas enthalpy [J]
    self%h_l = self%m_l*s_hl0
    self%h_g = self%m_g*s_hg0
    ! calculate initial liquid to gas heating rate [W]
    self%qlg_i = self%htc_i*self%csa_z*(self%t_l-self%t_g)
  end subroutine init_pvessel

  subroutine calc_res_pvessel(self, res)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(pvessel), intent(inout) :: self
    real(c_double), intent(out) :: res(:)
    !=== local variables ===
    real(c_double) :: rg, s_hg, s_hl
    real(c_double) :: pgeos, rleos, tgeos, tleos
    !======= Internals ============
    ! initialization
    res = 0.d0
    ! mass equations
    res(1) = self%dmldt-self%w_su              ! liq mass
    res(2) = self%dmgdt-self%wg_in+self%wg_ex  ! gas mass
    ! energy equations
    res(3) = self%dhldt-self%wh_su+self%qlg_i                 ! liq energy
    res(4) = self%dhgdt-self%wghg_in+self%wghg_ex-self%qlg_i  ! gas energy
    ! eos equations
    rg   = self%m_g/self%v_g
    s_hg = self%h_g/self%m_g
    s_hl = self%h_l/self%m_l
    call pgeos_n2g(eosopt, rg, s_hg, pgeos)
    call tleos_h2o(self%p_g, s_hl, tleos)
    call tgeos_n2g(eosopt, self%p_g, s_hg, tgeos)
    res(5) = self%p_g-pgeos ! gas pressure
    res(6) = self%t_l-tleos ! liq temperature
    res(7) = self%t_g-tgeos ! gas temperature
    ! volume equations
    call rleos_h2o(self%p_g, s_hl, rleos)
    res(8) = self%v_l-self%m_l/rleos ! liq volume
    res(9) = self%v_l+self%v_g-self%vol_t! total volume
    ! 
  end subroutine calc_res_pvessel


end module pvessel_mod