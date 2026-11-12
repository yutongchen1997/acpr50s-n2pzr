! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: flubndc (fluid boundary cell).
! 
! ------------------------------------------------------------------
module flubndc_mod
  !======= Inclusions =========
  use corrsub_mod

  !======= Declarations =========
  implicit none

  type :: flubndc
    !=== Member Variables === 
    ! required input parameters ------------------------------------
    ! id string, fbc type, fluid material
    character(len=20) :: id, typ, mat
    ! signal id string of fbc pressure and temperature
    character(len=20) :: sig_p, sig_t
    ! fluid pressure [Pa], fluid temperature [K]
    real(c_double) :: p_f, t_f
    ! fluid density [kg/m^3]
    real(c_double) :: r_f
    ! derived variables --------------------------------------------
    ! fluid enthalpy [J/kg]
    real(c_double) :: s_hf

  contains
    !=== Member Procedures ===
    procedure :: init_flubndc
    procedure :: updt_flubndc

  end type flubndc

contains

  subroutine init_flubndc(self, id, typ, mat, sig_p, sig_t)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flubndc), intent(inout) :: self
    ! fluid bc id, fluid bc type, fluid material type
    character(len=20), intent(in) :: id, typ, mat
    ! pressure signal id, temperature signal id
    character(len=20), intent(in) :: sig_p, sig_t

    !=== local variables ===
    
    !======= Internals ============
    ! required parameters
    self%id = id
    ! four fbc types: mass_flow_inlet, pressure_inlet, mass_flow_outlet, pressure_outlet
    self%typ = typ ! 4 fbc types: 'w_in', 'p_in', 'w_out', p_out'
    self%mat = mat ! 2 flu types: 'n2g', 'h2o'
    self%sig_p = sig_p
    self%sig_t = sig_t
    !
    self%p_f = 0.d0
    self%t_f = 0.d0
    self%s_hf= 0.d0
    self%r_f = 0.d0
    
  end subroutine init_flubndc

  subroutine updt_flubndc(self)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flubndc), intent(inout) :: self
    ! update fluid enthalpy with updated self%pf and self%tf
    select case (self%mat)
      case('n2g')
        call hgeos_n2g(eosopt, self%p_f, self%t_f, self%s_hf)
        call rgeos_n2g(eosopt, self%p_f, self%s_hf, self%r_f)
      case('h2o')
        call hfeos_h2o(self%p_f, self%t_f, self%s_hf)
        call rfeos_h2o(self%p_f, self%s_hf, self%r_f)
      case default
        print*, 'no match for fluid material = ', self%mat
        stop 'program stopped at updt_flubndc'
    end select
  end subroutine updt_flubndc

end module flubndc_mod