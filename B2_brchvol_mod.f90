! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: brchvol (bracn volume that could be
! linked by several junctions).
! 
! ------------------------------------------------------------------
module brchvol_mod
  !======= Inclusions =========
  use ctlsgnl_mod

  !======= Declarations =========
  implicit none

  type :: brchvol
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
    procedure :: init_brchvol ! initialize junction parameters
  end type brchvol

contains

  subroutine init_brchvol(self, id)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(brchvol), intent(inout) :: self
    character(len=20), intent(in) :: id
    !======= Internals ============

  end subroutine init_brchvol

end module brchvol_mod