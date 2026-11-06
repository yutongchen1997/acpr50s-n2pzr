! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: flubrch (fluid bracn volume that 
! could be linked by multible junctions).
! 
! ------------------------------------------------------------------
module flubrch_mod
  !======= Inclusions =========
  use ctlsgnl_mod

  !======= Declarations =========
  implicit none

  type :: flubrch
    !=== Member Variables ===
    ! required input parameters ------------------------------------
    ! junction id, fluid material type
    character(len=20) :: id, mat
    ! branch fluid volume
    real(c_double) :: vol
    ! boundary i/o effects -----------------------------------------
    ! mass flow rates: inlet, outlet
    real(c_double) :: wf_in, wf_ex
    ! enthalpy flow rates: inlet, outlet
    real(c_double) :: wfhf_in, wfhf_ex
    ! dynamic variables --------------------------------------------
    ! mass, enthalpy, temperature, pressure
    real(c_double) :: m_f    ! fluid mass         - [kg]
    real(c_double) :: h_f    ! fluid enthalpy     - [J]
    real(c_double) :: t_f    ! fluid temperature  - [K]
    real(c_double) :: p_f    ! fluid pressure     - [Pa]
    ! derivatives of dynamic variables -----------------------------
    real(c_double) :: dmfdt  ! time derivative of m  - [kg/s]
    real(c_double) :: dhfdt  ! time derivative of h  - [J/s]
    real(c_double) :: dtfdt  ! time derivative of t  - [K/s]
    real(c_double) :: dpfdt  ! time derivative of pg - [Pa/s]
    ! internal variables -------------------------------------------
    real(c_double) :: q_htr  ! electric heater power - [W]

  contains
    !=== Member Procedures ===
    procedure :: init_flubrch ! initialize junction parameters
  end type flubrch

contains

  subroutine init_flubrch(self, id, mat, inp)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flubrch), intent(inout) :: self
    character(len=20), intent(in) :: id, mat
    ! inp array: 1-vol0, 2-tf0, 3-pf0
    real(c_double) :: inp(:)
    !=== local variables ===
    real(c_double) :: rf0, s_hf0
    !======= Internals ============
    ! assign the input strings
    self%id = id
    self%mat= mat
    ! assign the input parameters
    self%vol = inp(1)
    self%t_f = inp(2)
    self%p_f = inp(3)
    ! calculate initial fluid specific enthalpy [J/kg]
    select case (self%mat)
      case('h2o')
        call hfeos_h2o(self%p_f,self%t_f,s_hf0)
      case('n2g')
        call hgeos_n2g(eosopt,self%p_f,self%t_f,s_hf0)
      case default
        print*, 'no match for fluid material = ', mat
        stop 'program stopped at init_flubrch'
    end select
    ! calculate initial fluid density [kg/m^3]
    select case (self%mat)
      case('h2o')
        call rfeos_h2o(self%p_f, s_hf0, rf0)
      case('n2g')
        call rgeos_n2g(eosopt, self%p_f, s_hf0, rf0)
      case default
        print*, 'no match for fluid material = ', mat
        stop 'program stopped at init_flubrch'
    end select
    ! initialize fluid mass and enthalpy

  end subroutine init_flubrch

end module flubrch_mod