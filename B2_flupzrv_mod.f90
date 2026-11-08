! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: pvessel (pressurizer vessel).
! 
! ------------------------------------------------------------------
module flupzrv_mod
  !======= Inclusions =========
  use ctlsgnl_mod

  !======= Declarations =========
  implicit none

  type :: flupzrv
    !=== Member Variables === 
    ! required input parameters ------------------------------------
    character(len=20) :: id
    ! id string of liq and gas section, material of liq and gas
    character(len=20) :: id_l, id_g, mat_l, mat_g, shape
    ! inner diameter, total vessel height
    real(c_double) :: dia_i, len_t
    ! other parameters ---------------------------------------------
    ! cross-sectional area, total vessel volume
    real(c_double) :: csa_z, vol_t
    ! volume when freelevel fills lower head, reaches upper head 
    real(c_double) :: vol_l, vol_u
    ! boundary conditions ------------------------------------------
    ! mass flow rates: gas inlet, gas exit, liq inlet, liq exit
    real(c_double) :: wg_in, wg_ex, wl_in, wl_ex
    ! enthalpy flow rates: gas inlet, gas exit, liq inlet, liq exit
    real(c_double) :: wghg_in, wghg_ex, wlhl_in, wlhl_ex
    ! wall heat transfer rates
    real(c_double) :: qw_g, qw_l
    ! dynamic variables --------------------------------------------
    ! mass, enthalpy, temperature, volume, pressure, level
    real(c_double) :: m_g, m_l  ! gas, liq mass         - [kg]
    real(c_double) :: h_g, h_l  ! gas, liq enthalpy     - [J]
    real(c_double) :: t_g, t_l  ! gas, liq temperature  - [K]
    real(c_double) :: v_g, v_l  ! gas, liq volume       - [m^3]
    real(c_double) :: p_g       ! gas pressure          - [Pa]
    ! derivatives of dynamic variables
    real(c_double) :: dmldt, dmgdt   ! time derivative of m  - [kg/s]
    real(c_double) :: dhldt, dhgdt   ! time derivative of h  - [J/s]
    real(c_double) :: dtldt, dtgdt   ! time derivative of t  - [K/s]
    real(c_double) :: dvldt, dvgdt   ! time derivative of v  - [m^3/s]
    real(c_double) :: dpdt           ! time derivative of pg - [Pa/s]
    ! internal variables -------------------------------------------
    real(c_double) :: htc_i      ! liq - gas interface htc    - [W/m^2/K]  
    real(c_double) :: qlg_i      ! liq - gas heating rate     - [W]
    real(c_double) :: z_l        ! internal liquid level      - [m]
    real(c_double) :: r_g, r_l   ! gas, liq density           - [kg/m^3]
    real(c_double) :: s_hg, s_hl ! gas, liq specific enthalpy - [J/kg]
    ! temperature probe variables ----------------------------------
    integer(c_int) :: n_tcs     ! number of TCs
    real(c_double), allocatable :: z_tcs(:)  ! axial location of TCs
    
  contains  
    !=== Member Procedures === 
    procedure :: init_flupzrv
    procedure :: calc_liq_vol
    procedure :: get_liq_lvl
    procedure :: calc_res_flupzrv

  end type flupzrv

contains

  subroutine init_flupzrv(self, id, id_l, id_g, mat_l, mat_g, shape, inp)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flupzrv), intent(inout) :: self
    character(len=20), intent(in) :: id, id_l, id_g, mat_l, mat_g, shape
    ! inp array: 1-dia_i, 2-len_t, 3-len_i, 4-tg0, 5-tl0, 6-pg0, 7-htc_0
    real(c_double), intent(in) :: inp(:) 
    !=== local variables ===
    real(c_double) :: rg0, rl0, s_hg0, s_hl0
    !======= Internals ============
    ! assign the input strings
    self%id = id
    self%id_l = id_l
    self%id_g = id_g
    self%mat_l = mat_l
    self%mat_g = mat_g
    self%shape = shape
    ! assign the input parameters
    self%dia_i = inp(1)
    self%len_t = inp(2)
    self%z_l   = inp(3)
    self%t_g   = inp(4)
    self%t_l   = inp(5)
    self%p_g   = inp(6)
    self%htc_i = inp(7)
    ! calculate cross-sectional area
    self%csa_z = pi*self%dia_i**2.d0/4.d0
    ! calculate total volume, initial liquid and gas volume
    call calc_liq_vol(self)
    self%v_g = self%vol_t-self%v_l
    ! calcualte initial liquid and gas specific enthalpy [J/kg]
    call hgeos_n2g(eosopt, self%p_g, self%t_g, s_hg0)
    call hfeos_h2o(self%p_g, self%t_l, s_hl0)        
    self%s_hg= s_hg0
    self%s_hl= s_hl0 
    ! calculate initial liquid and gas density [kg/m^3]
    call rgeos_n2g(eosopt, self%p_g, s_hg0, rg0)
    call rfeos_h2o(self%p_g, s_hl0, rl0)
    self%r_g = rg0
    self%r_l = rl0
    ! calculate initial liquid and gas mass [kg]
    self%m_g = self%r_g*self%v_g
    self%m_l = self%r_l*self%v_l
    ! calculate initial liquid and gas enthalpy [J]
    self%h_l = self%s_hl*self%m_l
    self%h_g = self%s_hg*self%m_g
    ! calculate initial liquid to gas freelevel heating rate [W]
    self%qlg_i = self%htc_i*self%csa_z*(self%t_l-self%t_g)
  end subroutine init_flupzrv

  subroutine calc_liq_vol(self)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flupzrv), intent(inout) :: self
    !=== local variables ===
    !======= Internals ============
    select case (self%shape)
      case('cylindrical')
        self%vol_t = self%csa_z*self%len_t
        self%vol_l = self%vol_t
        self%vol_u = self%vol_t
        self%v_l = self%csa_z*self%z_l
      case('hemispherical')
        self%vol_l = (pi*self%dia_i**3.d0)/12.d0  ! volume of the hemisphere
        self%vol_u = self%vol_l + self%csa_z*(self%len_t-1.0*self%dia_i)
        self%vol_t = self%vol_u + self%vol_l  
        if(self%z_l<0.5*self%dia_i) then
          self%v_l = (pi*self%z_l**2.d0)*(1.5d0*self%dia_i-self%z_l)/3.d0
        else if(self%z_l<self%len_t-0.5*self%dia_i) then
          self%v_l = (pi*self%dia_i**3.d0)/12.d0
          self%v_l = self%v_l + self%csa_z*(self%z_l-0.5*self%dia_i)
        else
          self%v_l = (pi*self%dia_i**3.d0)/12.d0 + self%csa_z*(self%len_t-self%dia_i)
          self%v_l = self%v_l + pi/3.d0*(0.25d0*self%dia_i**3.d0 + &
          (self%len_t-self%z_l)**3.d0-1.5d0*(self%len_t-self%z_l)**2.d0*self%dia_i)
        end if
      case default
        print*, 'no match for PV shape = ', self%shape
        stop 'program stopped at init_pvessel'
    end select
  end subroutine calc_liq_vol

  subroutine get_liq_lvl(self)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flupzrv), intent(inout) :: self
    !=== local variables ===
    real(c_double) :: zg, vg
    !======= Internals ============
    select case (self%shape)
      case('cylindrical')
        self%z_l = self%v_l/self%csa_z
      case('hemispherical')
        if(self%v_l<self%vol_l) then
          self%z_l = (self%dia_i/2.d0)*(1.d0-2.d0*SIN(ASIN(1.d0-12.d0*self%v_l/pi/(self%dia_i**3.d0))/3.d0))
        else if(self%v_l<self%vol_u) then
          self%z_l = (self%v_l-self%vol_l)/self%csa_z + self%dia_i/2.d0
        else
          vg = self%vol_t - self%v_l
          zg = (self%dia_i/2.d0)*(1.d0-2.d0*SIN(ASIN(1.d0-12.d0*vg/pi/(self%dia_i**3.d0))/3.d0))
          self%z_l = self%len_t - zg
        end if
      case default
        print*, 'no match for PV shape = ', self%shape
        stop 'program stopped at init_pvessel'
    end select
  end subroutine get_liq_lvl

  subroutine calc_res_flupzrv(self, res)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flupzrv), intent(inout) :: self
    real(c_double), intent(out) :: res(:)
    !=== local variables ===
    real(c_double) :: pgeos, rleos, tgeos, tleos
    !======= Internals ============
    ! initialization
    res = 0.d0
    ! mass equations
    res(1) = self%dmldt-self%wl_in+self%wl_ex  ! liq mass
    res(2) = self%dmgdt-self%wg_in+self%wg_ex  ! gas mass
    print*, 'self%wg_in-pzrv', self%wg_in, self%wg_ex
    ! energy equations
    self%qlg_i = self%htc_i*self%csa_z*(self%t_l-self%t_g)    ! liq-gas freelevel heating rate
    res(3) = self%dhldt-self%wlhl_in+self%wlhl_ex+self%qlg_i  ! liq energy
    res(4) = self%dhgdt-self%wghg_in+self%wghg_ex-self%qlg_i  ! gas energy
    ! eos equations
    call pgeos_n2g(eosopt, self%r_g, self%s_hg, pgeos)
    call tfeos_h2o(self%p_g, self%s_hl, tleos)
    call tgeos_n2g(eosopt, self%p_g, self%s_hg, tgeos)
    res(5) = self%p_g-pgeos ! gas pressure
    res(6) = self%t_l-tleos ! liq temperature
    res(7) = self%t_g-tgeos ! gas temperature
    ! volume equations
    call rfeos_h2o(self%p_g, self%s_hl, rleos)
    res(8) = self%v_l-self%m_l/rleos       ! liq volume
    res(9) = self%v_l+self%v_g-self%vol_t  ! total volume
    ! 
  end subroutine calc_res_flupzrv




end module flupzrv_mod