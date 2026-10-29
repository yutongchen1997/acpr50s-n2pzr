! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: sysdomn (storing, reading and writing
! system data ).
! 
! ------------------------------------------------------------------
module sysdomn_mod
  !======= Inclusions =========
  use fludomn_mod
  use slddomn_mod

  !======= Declarations =========
  implicit none
  
  type :: sysdomn
    !=== Member Variables === 
    ! time variables
    real(c_double) :: tstart, tend, dtout
    ! total number of daes
    integer(c_int) :: neq_tot = 0
    ! total number of signals
    integer(c_int) :: numsgnl = 0
    ! fluid domain object
    type(fludomn) :: fluid
    ! solid domain object
    type(slddomn) :: solid
    ! core object...
    ! ------------------------------
    ! control signal objects array
    type(ctlsgnl), allocatable :: signals(:)
  contains
    !=== Member Procedures === 
    procedure :: write_to_y
    procedure :: read_from_y

  end type sysdomn
  ! build the system object
  ! system_glo is the sysdomn-type object that stores all system data
  type(sysdomn) :: system_glo

contains
! ----------------------------------------------------------------------
! write_to_y: write unknowns to yval-vector for initialization
! ----------------------------------------------------------------------
  subroutine write_to_y(self, yval)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(sysdomn), intent(inout) :: self
    real(c_double), intent(out)   :: yval(:)
    !== local variables ==
    integer(c_int) :: i, j, k

    !======= Internals ============

  end subroutine write_to_y

! ----------------------------------------------------------------------
! read_from_y: read unknowns from y-vector
! ----------------------------------------------------------------------
  subroutine read_from_y(self, yval)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(sysdomn), intent(inout) :: self
    real(c_double), intent(in)    :: yval(:)
    !== local variables ==
    integer(c_int) :: i, j, k

    !======= Internals ============

  end subroutine read_from_y


end module sysdomn_mod