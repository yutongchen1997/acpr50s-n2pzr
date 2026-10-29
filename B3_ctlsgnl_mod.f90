! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: ctlsgnl (control signal).
! 
! ------------------------------------------------------------------
module ctlsgnl_mod
  !======= Inclusions ===========
  use corrsub_mod

  !======= Declarations =========
  implicit none

  type :: ctlsgnl
    !=== Member Variables === 
    ! signal id, signal type
    character(len=20) :: id, typ
    ! latest signal value
    real(c_double) :: val
    ! signal look-up table (only for 'time-dependent' type signal)
    real(c_double), allocatable :: tab(:,:)
  contains
    !=== Member Procedures === 
    procedure :: init_sgnl
  end type ctlsgnl

contains

  subroutine init_sgnl(self, sigid, sigtp, sgntab)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(ctlsgnl), intent(inout) :: self
    character(len=20), intent(in) :: sigid, sigtp
    real(c_double),     optional :: sgntab(:,:)
    !== Local variables ==
    integer(c_int), allocatable :: lentab
    !======= Internals =========
    self%id = sigid
    self%typ= sigtp
    select case(self%typ)
      case('const')  ! constant
        continue
      case('tdsig')  ! time-dependent signal
        lentab = size(sgntab(1,:))
        allocate(self%tab(2,lentab))
        self%tab = sgntab
      case('lutab')  ! look-up table signal
        print*, 'now under development...'
        stop 'program stopped at ctlsgnl%init_sgnl'
      case default
        print*, 'no match for signal type = ', self%typ
        stop 'program stopped at ctlsgnl%init_sgnl'
    end select
  end subroutine init_sgnl

end module ctlsgnl_mod