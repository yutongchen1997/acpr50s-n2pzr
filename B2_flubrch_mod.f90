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
    ! signal id string of electric heater power
    character(len=20) :: sig_ht
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
    real(c_double) :: q_htr  ! electric heater power   - [W]
    real(c_double) :: r_f    ! fluid density           - [kg/m^3]
    real(c_double) :: s_hf   ! fluid specific enthaply - [J/kg]

  contains
    !=== Member Procedures ===
    procedure :: init_flubrch ! initialize junction parameters
    procedure :: calc_res_flubrch
  end type flubrch

contains

  subroutine init_flubrch(self, id, mat, sig_ht, inp)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flubrch), intent(inout) :: self
    character(len=20), intent(in) :: id, mat, sig_ht
    ! inp array: 1-vol0, 2-tf0, 3-pf0
    real(c_double) :: inp(:)
    !=== local variables ===
    real(c_double) :: rf0, s_hf0
    !======= Internals ============
    ! assign the input strings
    self%id = id
    self%mat= mat ! 2 flu types: 'n2g', 'h2o'
    self%sig_ht = sig_ht
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
    ! initialize fluid density, mass and enthalpy
    self%r_f = rf0
    self%s_hf= s_hf0
    self%m_f = self%r_f*self%vol
    self%h_f = self%m_f*self%s_hf
    !
    self%q_htr = 0.d0

  end subroutine init_flubrch

  subroutine calc_res_flubrch(self, res)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flubrch), intent(inout) :: self
    real(c_double), intent(out) :: res(:)
    !=== local variables ===
    real(c_double) :: rf, s_hf, tfeos, pfeos
    !======= Internals ============
    ! initialization 
    res = 0.d0
    ! mass equation
    res(1) = self%dmfdt-self%wf_in+self%wfhf_ex
    ! energy equation
    res(2) = self%dhfdt-self%wfhf_in+self%wfhf_ex-self%q_htr
    ! eos equations
    rf   = self%m_f/self%vol
    s_hf = self%h_f/self%m_f
    call tfeos_h2o(self%p_f, s_hf, tfeos)
    res(3) = self%t_f-tfeos
    call pfeos_h2o(rf, s_hf, pfeos)
    res(4) = self%p_f-pfeos
  end subroutine calc_res_flubrch

end module flubrch_mod