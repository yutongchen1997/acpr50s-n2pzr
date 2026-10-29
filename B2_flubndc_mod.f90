! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: flubndc (fluid boundary cell).
! 
! ------------------------------------------------------------------
module flubndc_mod
  !======= Inclusions =========
  use ctlsgnl_mod

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
    real(c_double) :: pf, tf
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
    self%pf  = 0.d0
    self%tf  = 0.d0
    self%s_hf= 0.d0
    
  end subroutine init_flubndc

  subroutine updt_flubndc(self, pf, tf)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flubndc), intent(inout) :: self
    real(c_double), intent(in) :: pf, tf
    ! update fluid enthalpy with updated pf and tf
    self%pf = pf
    self%tf = tf
    select case (self%mat)
      case('n2g')
        call hgeos_n2g(eosopt, self%pf, self%tf, self%s_hf)
      case('h2o')
        call hleos_h2o(self%pf, self%tf, self%s_hf)
      case default
        print*, 'no match for fluid material = ', self%mat
        stop 'program stopped at updt_flubndc'
    end select
  end subroutine updt_flubndc

end module flubndc_mod