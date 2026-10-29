! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! The following is the module for composing residual vector of the 
! Differential Algebraic Equation set.
! 
! ------------------------------------------------------------------

module rescomp_mod

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod

  !======= Declarations =========
  implicit none

  integer, parameter :: myindextype = selected_int_kind(16)

  integer(kind=myindextype), parameter :: neq = 3
  integer(kind=myindextype), parameter :: nout = 12

contains

  ! ----------------------------------------------------------------
  ! resrob: The DAE residual function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function resrob(tres, sunvec_y, sunvec_yp, sunvec_r, user_data) &
    result(ierr) bind(C, name='resrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: tres      ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_yp ! derivative N_Vector
    type(N_Vector)        :: sunvec_r  ! residual N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yval(:)
    real(c_double), pointer :: ypval(:)
    real(c_double), pointer :: rval(:)

    !======= Internals ============

    ! get data arrays from SUNDIALS vectors
    yval => FN_VGetArrayPointer(sunvec_y)
    ypval => FN_VGetArrayPointer(sunvec_yp)
    rval => FN_VGetArrayPointer(sunvec_r)

    ! fill residual vector
    rval(1) = -0.04d0*yval(1) + 1.0d4*yval(2)*yval(3)
    rval(2) = -rval(1) - 3.0d7*yval(2)**2 - ypval(2)
    rval(1) = rval(1) - ypval(1)
    rval(3) = yval(1) + yval(2) + yval(3) - 1.d0

    ! return success
    ierr = 0
    return

  end function resrob

end module rescomp_mod
! ------------------------------------------------------------------