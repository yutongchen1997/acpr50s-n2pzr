! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: flujunc (junctions that link fluid
! components).
! 
! ------------------------------------------------------------------
module flujunc_mod
  !======= Inclusions =========
  use ctlsgnl_mod

  !======= Declarations =========
  implicit none

  type :: flujunc
    !=== Member Variables ===
    ! required input parameters ------------------------------------
    ! junction id and type
    character(len=20) :: id, typ
    ! 'from' and 'to' fluid component id string
    character(len=20) :: fid_f, fid_t
    ! signal id (optional)
    character(len=20) :: sig
    ! junction values: pump head [Pa], local loss coefficient [-]
    real(c_double) :: phed, kfac
    ! junction variable: mass flow rate [kg/s]
    real(c_double) :: mdot 

  contains
    !=== Member Procedures ===
    procedure :: init_flujunc ! initialize junction parameters
  end type flujunc

contains

  subroutine init_flujunc(self, juntp, fid_f, fid_t, sig)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flujunc), intent(inout) :: self
    character(len=20), intent(in) :: juntp, fid_f, fid_t
    character(len=20), optional, intent(in) :: sig
    !======= Internals ============
    self%typ   = juntp  ! junction type: 'junc_flow', 'junc_valv', 'junc_pump', 'junc_norm'
    self%fid_f = fid_f
    self%fid_t = fid_t
    if(present(sig)) then
      self%sig = sig
    end if
    !
    self%phed = 0.d0
    self%kfac  = 0.d0
    self%mdot  = 0.d0
  end subroutine init_flujunc

end module flujunc_mod