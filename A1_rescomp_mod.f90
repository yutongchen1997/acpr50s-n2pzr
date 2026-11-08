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
  use control_mod
  use mfldcpl_mod
  use fsundials_core_mod

  !======= Declarations =========
  implicit none
  ! integers

contains

  ! ----------------------------------------------------------------
  ! resrob: The DAE residual function
  !
  ! Return values:
  !    0 = success,
  !    1 = recoverable error,
  !   -1 = non-recoverable error
  ! ----------------------------------------------------------------
  integer(c_int) function resrob(time, sunvec_y, sunvec_yp, sunvec_r, user_data) &
    result(ierr) bind(C, name='resrob')

    !======= Inclusions ===========
    use, intrinsic :: iso_c_binding

    use fnvector_serial_mod

    !======= Declarations =========
    implicit none

    ! calling variables
    real(c_double), value :: time      ! current time
    type(N_Vector)        :: sunvec_y  ! solution N_Vector
    type(N_Vector)        :: sunvec_yp ! derivative N_Vector
    type(N_Vector)        :: sunvec_r  ! residual N_Vector
    type(c_ptr), value :: user_data ! user-defined data

    ! pointers to data in SUNDIALS vectors
    real(c_double), pointer :: yval(:)
    real(c_double), pointer :: ypval(:)
    real(c_double), pointer :: rval(:)

    ! locally defined arrays/variables
    real(c_double), allocatable :: res(:)
    real(c_double), allocatable :: res_fluid(:)
    real(c_double), allocatable :: res_solid(:)
    real(c_double) :: t1, t2  ! CPU time

    !======= Internals ============
    call CPU_TIME(t1)
    ! update the signal values
    call evaluate_signals(time) 
    
    ! get data arrays from SUNDIALS vectors
    yval  => FN_VGetArrayPointer(sunvec_y)
    ypval => FN_VGetArrayPointer(sunvec_yp)
    rval  => FN_VGetArrayPointer(sunvec_r)
    
    ! read the system_glo data from SUNDIALS vectors
    call system_glo%read_from_y_and_yp(yval, ypval)
    
    ! allocate residual array
    allocate(res(neq_tot))
    res = 0.d0
    ! fill residual vector
    if(solve_fluid) then
      allocate(res_fluid(neq_mdot+neq_pzrv+neq_brch))
      call system_glo%fluid%calc_res_fludomn(res_fluid)
    end if
    res = (/res_fluid/)
    rval= res
    deallocate(res)
    deallocate(res_fluid)
    ! return success
    ierr = 0
    !
    call CPU_TIME(t2)
    !
    return

  end function resrob

! ----------------------------------------------------------------------
! set_id_vector: a subroutine that builds the IDA variable ID vector:  
! = 1.0 means differential variable - (ODE)
! = 0.0 means algebraic variable    - (AE)
! output - the processed id vector
! ----------------------------------------------------------------------
  subroutine set_id_vector(id)
  
  !======= Declarations =========
  implicit none
  !== output variables ==
  real(c_double), intent(out) :: id(:)
  !== local variables ==
  integer(c_int) :: i, j, k

  !======= Internals ============
  if(solve_fluid) then
    ! mdots - ODE
    do i = 1, neq_mdot
      id(i) = 0.0    
    end do
    ! equations relative to pzrvs
    do i = neq_mdot+1, neq_mdot+neq_pzrv
      select case(mod(i-neq_mdot, 9))
        case(1,2,3,4)
          id(i) = 1.d0  ! mass, energy - ODE
        case(5,6,7,8,0)
          id(i) = 0.d0  ! pressure, temperature, volume - AE
        case default
          print*, 'under aktiv development 2026/11/08-2030'
      end select
    end do
    ! equations relative to brchs
    do i = neq_mdot+neq_pzrv+1, neq_mdot+neq_pzrv+neq_brch
      select case(mod(i-neq_mdot-neq_pzrv, 4))
        case(1,2)
          id(i) = 1.d0
        case(3,4)
          id(i) = 0.d0
        case default
          print*, 'under aktiv development 2026/11/08-2039'
      end select
    end do
  end if

  end subroutine set_id_vector





end module rescomp_mod
! ------------------------------------------------------------------