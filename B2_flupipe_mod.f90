! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: flupipe (pipe assembly).
! 
! ------------------------------------------------------------------
module flupipe_mod
  !======= Inclusions =========
  use corrsub_mod

  !======= Declarations =========
  implicit none

  type :: flupipe
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
    procedure :: init_flupipe
    procedure :: calc_res_flupipe

  end type flupipe

contains

  subroutine init_flupipe(self)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flupipe), intent(inout) :: self
    !=== local variables ===

    !======= Internals ============



  end subroutine init_flupipe

 

  subroutine calc_res_flupipe(self, res)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flupipe), intent(inout) :: self
    real(c_double), intent(out) :: res(:)
    !=== local variables ===
      
    !======= Internals ============
    ! initialization
      
    ! 
  end subroutine calc_res_flupipe




end module flupipe_mod