! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! The following is the 'main' program for analyzing nitrogen gas 
! pressurizer dynamics.
! 
! ------------------------------------------------------------------

program acprpzr

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fsundials_core_mod            ! Access SUNDIALS core types, data structures, etc.
  use fida_mod                      ! Fortran interface to IDA
  use fnvector_serial_mod           ! Fortran interface to serial N_Vector
  use fsunmatrix_dense_mod          ! Fortran interface to dense SUNMatrix

  use fsunlinsol_dense_mod          ! 0. Fortran interface to dense SUNLinearSolver
  use fsunlinsol_pcg_mod            ! 1. Fortran interface to pcg SUNLinearSolver
  use fsunlinsol_spbcgs_mod         ! 2. Fortran interface to spbcgs SUNLinearSolver
  use fsunlinsol_spfgmr_mod         ! 3. Fortran interface to spfgmr SUNLinearSolver
  use fsunlinsol_spgmr_mod          ! 4. Fortran interface to spgmr SUNLinearSolver
  use fsunlinsol_sptfqmr_mod        ! 5. Fortran interface to sptfqmr SUNLinearSolver

  use fsunnonlinsol_newton_mod      ! Fortran interface to Newton SUNNonlinearSolver

  use rescomp_mod                   ! Compose residuals of DAE functions

  !======= Declarations =========
  implicit none

  ! local variables
  real(c_double) :: rtol, t0, tout1, tout, tret(1)
  ! CPU time variables
  real(c_double) :: tic, tac
  integer(c_int) :: iout, retval, retvalr, nrtfn, rootsfound(2)

  type(N_Vector), pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_yp     ! sundials derivative vector
  type(N_Vector), pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: ida_mem       ! IDA memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context

  ! solution and tolerance vectors, neq is set in the dae_mod module
  real(c_double) :: yval(neq), ypval(neq), avtol(neq)

  !======= Internals ============
  call cpu_time(tic)

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  ! initialize solution vectors and tolerances
  yval(1) = 1.d0
  yval(2) = 0.d0
  yval(3) = 0.d0

  ypval(1) = -0.040
  ypval(2) = 0.04d0
  ypval(3) = 0.d0

  rtol = 1.d-4

  avtol(1) = 1.d-8
  avtol(2) = 1.d-6
  avtol(3) = 1.d-6

  ! create serial vectors
  sunvec_y => FN_VMake_Serial(neq, yval, sunctx)
  if (.not. associated(sunvec_y)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  sunvec_yp => FN_VMake_Serial(neq, ypval, sunctx)
  if (.not. associated(sunvec_yp)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  sunvec_av => FN_VMake_Serial(neq, avtol, sunctx)
  if (.not. associated(sunvec_av)) then
    print *, 'ERROR: sunvec = NULL'
    stop 1
  end if

  ! set integration limits
  t0 = 0.d0
  tout1 = 0.4d0

  call PrintHeader(rtol, avtol, yval)

  ! Call FIDACreate and FIDAInit to initialize IDA memory
  ida_mem = FIDACreate(sunctx)
  if (.not. c_associated(ida_mem)) then
    print *, 'ERROR: ida_mem = NULL'
    stop 1
  end if

  retval = FIDAInit(ida_mem, c_funloc(resrob), t0, sunvec_y, sunvec_yp)
  if (retval /= 0) then
    print *, 'Error in FIDAInit, retval = ', retval, '; halting'
    stop 1
  end if

  ! Call FIDASVtolerances to set tolerances
  retval = FIDASVtolerances(ida_mem, rtol, sunvec_av)
  if (retval /= 0) then
    print *, 'Error in FIDASVtolerances, retval = ', retval, '; halting'
    stop 1
  end if

  ! Create dense SUNMatrix for use in linear solves
  sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
  if (.not. associated(sunmat_A)) then
    print *, 'ERROR: sunmat = NULL'
    stop 1
  end if

  ! Create dense SUNLinearSolver object
  sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
  if (.not. associated(sunlinsol_LS)) then
    print *, 'ERROR: sunlinsol = NULL'
    stop 1
  end if

  ! Attach the matrix and linear solver
  retval = FIDASetLinearSolver(ida_mem, sunlinsol_LS, sunmat_A); 
  if (retval /= 0) then
    print *, 'Error in FIDASetLinearSolver, retval = ', retval, '; halting'
    stop 1
  end if

  ! Create Newton SUNNonlinearSolver object. IDA uses a
  ! Newton SUNNonlinearSolver by default, so it is not necessary
  ! to create it and attach it. It is done in this example code
  ! solely for demonstration purposes.
  sunnonlin_NLS => FSUNNonlinSol_Newton(sunvec_y, sunctx)
  if (.not. associated(sunnonlin_NLS)) then
    print *, 'ERROR: sunnonlinsol = NULL'
    stop 1
  end if

  ! Attach the nonlinear solver
  retval = FIDASetNonlinearSolver(ida_mem, sunnonlin_NLS)
  if (retval /= 0) then
    print *, 'Error in FIDASetNonlinearSolver, retval = ', retval, '; halting'
    stop 1
  end if

  ! In loop, call IDASolve, print results, and test for error.
  ! Break out of loop when NOUT preset output times have been reached.

  iout = 0
  tout = tout1
  do

    retval = FIDASolve(ida_mem, tout, tret, sunvec_y, sunvec_yp, IDA_NORMAL)
    if (retval < 0) then
      print *, 'Error in FIDASolve, retval = ', retval, '; halting'
      stop 1
    end if

    call PrintOutput(ida_mem, tret(1), yval)

    if (retval == IDA_ROOT_RETURN) then
      retvalr = FIDAGetRootInfo(ida_mem, rootsfound)
      if (retvalr < 0) then
        print *, 'Error in FIDAGetRootInfo, retval = ', retval, '; halting'
        stop 1
      end if
      print '(a,2(i2,2x))', "    rootsfound[] = ", rootsfound(1), rootsfound(2)
    end if

    if (retval == IDA_SUCCESS) then
      iout = iout + 1
      tout = tout*10.d0
    end if

    if (iout == NOUT) exit

  end do

  call PrintFinalStats(ida_mem)

  ! free memory
  call FIDAFree(ida_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_av)
  call FN_VDestroy(sunvec_yp)
  retval = FSUNContext_Free(sunctx)

  call cpu_time(tac)
  print*, 'total CPU time = ', tac - tic
end program acprpzr

! ----------------------------------------------------------------
! PrintHeader: prints first lines of output (problem description)
! ----------------------------------------------------------------
subroutine PrintHeader(rtol, avtol, y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use rescomp_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  real(c_double) :: rtol
  real(c_double) :: avtol(neq)
  real(c_double) :: y(neq)

  !======= Internals ============

  print *, " "
  print *, "idaRoberts_dns_f2003: Robertson kinetics DAE serial example problem for IDA"
  print *, "         Three equation chemical kinetics problem."
  print *, " "
  print *, "Linear solver: DENSE, with user-supplied Jacobian."
  print '(a,f6.4,a,3(es7.0,1x))', "Tolerance parameters:  rtol = ", rtol, "   atol = ", avtol
  print '(a,3(f5.2,1x),a)', "Initial conditions y0 = (", y, ")"
  print *, "Constraints and id not used."
  print *, " "
  print *, "-----------------------------------------------------------------------"
  print *, "  t             y1           y2           y3     | nst  k      h"
  print *, "-----------------------------------------------------------------------"

  return
end subroutine PrintHeader

! ----------------------------------------------------------------
! PrintOutput
! ----------------------------------------------------------------
subroutine PrintOutput(ida_mem, t, y)

  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding
  use fida_mod
  use rescomp_mod

  !======= Declarations =========
  implicit none

  ! calling variable
  type(c_ptr)    :: ida_mem
  real(c_double) :: t, y(neq)

  ! internal variables
  integer(c_int)  :: retval, kused(1)
  integer(c_long) :: nst(1)
  real(c_double)  :: hused(1)

  !======= Internals ============

  retval = FIDAGetLastOrder(ida_mem, kused)
  if (retval /= 0) then
    print *, 'Error in FIDAGetLastOrder, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumSteps(ida_mem, nst)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumSteps, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetLastStep(ida_mem, hused)
  if (retval /= 0) then
    print *, 'Error in FIDAGetLastStep, retval = ', retval, '; halting'
    stop 1
  end if

  print '(es12.4,1x,3(es12.4,1x),a,i3,2x,i1,1x,es12.4)', &
    t, y(1), y(2), y(3), "| ", nst, kused(1), hused(1)

end subroutine PrintOutput

! ----------------------------------------------------------------
! PrintFinalStats
!
! Print KINSOL statstics to standard out
! ----------------------------------------------------------------
subroutine PrintFinalStats(ida_mem)

  !======= Inclusions ===========
  use iso_c_binding
  use fida_mod

  !======= Declarations =========
  implicit none

  type(c_ptr), intent(in) :: ida_mem
  integer(c_int)  :: retval
  integer(c_long) :: nst(1), nni(1), nje(1), nre(1), nreLS(1), nge(1), ncfn(1), netf(1)

  !======= Internals ============

  retval = FIDAGetNumSteps(ida_mem, nst)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumSteps, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumResEvals(ida_mem, nre)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumResEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumJacEvals(ida_mem, nje)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumJacEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumNonlinSolvIters(ida_mem, nni)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumNonlinSolvIters, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumErrTestFails(ida_mem, netf)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumErrTestFails, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumNonlinSolvConvFails(ida_mem, ncfn)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumNonlinSolvConvFails, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumLinResEvals(ida_mem, nreLS)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumLinResEvals, retval = ', retval, '; halting'
    stop 1
  end if

  retval = FIDAGetNumGEvals(ida_mem, nge)
  if (retval /= 0) then
    print *, 'Error in FIDAGetNumGEvals, retval = ', retval, '; halting'
    stop 1
  end if

  print *, " "
  print *, "Final Run Statistics: "
  print *, "Number of steps                    = ", nst
  print *, "Number of residual evaluations     = ", nre + nreLS
  print *, "Number of Jacobian evaluations     = ", nje
  print *, "Number of nonlinear iterations     = ", nni
  print *, "Number of error test failures      = ", netf
  print *, "Number of nonlinear conv. failures = ", ncfn
  print *, "Number of root fn. evaluations     = ", nge

  return

end subroutine PrintFinalStats
