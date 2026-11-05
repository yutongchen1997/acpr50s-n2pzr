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
    ! junction variable: time derivative of mass flow rate [kg/s^2]
    real(c_double) :: dmdotdt

  contains
    !=== Member Procedures ===
    procedure :: init_flubrch ! initialize junction parameters
  end type flubrch

contains

  subroutine init_flubrch(self, id)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(flubrch), intent(inout) :: self
    character(len=20), intent(in) :: id
    !======= Internals ============

  end subroutine init_flubrch

end module flubrch_mod