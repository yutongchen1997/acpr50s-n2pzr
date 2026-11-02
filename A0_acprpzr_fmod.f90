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

  ! time variables
  real(c_double) :: tout, tret(1)
  ! CPU time variables
  real(c_double) :: tic, tac
  integer(c_int) :: iout, retval
  ! integers
  integer, parameter :: myindextype = selected_int_kind(16)
  integer(kind=myindextype), parameter :: NOUT = 100
  integer(kind=myindextype) :: neq = 0
  ! SUNDIALS types
  type(N_Vector), pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_yp     ! sundials derivative vector
  type(N_Vector), pointer :: sunvec_av     ! sundials tolerance vector
  type(SUNMatrix), pointer :: sunmat_A      ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: ida_mem       ! IDA memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context

  ! solution and tolerance vectors, neq is set in the dae_mod module
  real(c_double) :: rtol
  real(c_double), allocatable :: avtol(:)
  real(c_double), allocatable :: yval(:), ypval(:)

  !======= Internals ============
  call CPU_TIME(tic)

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  ! open the input file and read input data
  call construct_input()
  ! open output files
  call open_output_files()
  ! convert neq type
  neq = int(neq_tot, kind = myindextype)
  ! set tolerance vectors
  rtol = rt
  ! 
  allocate(avtol(neq))
  if(neq_mdot>0) then
    ! under active development ... 2026/11/02-1453
  end if
  if(neq_pzrv>0) then
  ! fluid mass
    avtol(1:2) = at(1)
  ! fluid energy
    avtol(3:4) = at(2)
  ! fluid pressure
    avtol(5)   = at(3)
  ! fluid temperature
    avtol(6:7) = at(4)
  ! fluid volume
    avtol(8:9) = at(5)
  end if
  ! ---
  call system_glo%write_to_y_and_yp(yval, ypval)

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



  ! Call FIDACreate and FIDAInit to initialize IDA memory
  ida_mem = FIDACreate(sunctx)
  if (.not. c_associated(ida_mem)) then
    print *, 'ERROR: ida_mem = NULL'
    stop 1
  end if

  retval = FIDAInit(ida_mem, c_funloc(resrob), tstart, sunvec_y, sunvec_yp)
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

  !
  select case (solopt)
    case(1) ! Dense direct linear solver (internal) (Dense)
      ! Create dense SUNMatrix for use in linear solves
      sunmat_A => FSUNDenseMatrix(neq, neq, sunctx)
      if (.not. associated(sunmat_A)) then
         print *, 'ERROR: sunmat = NULL'
         stop 1
      end if
      ! Create dense SUNLinearSolver object
      sunlinsol_LS => FSUNLinSol_Dense(sunvec_y, sunmat_A, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunlinsol = NULL, Dense'
         stop 1
      end if
    case(2) ! Preconditioned Conjugate Gradient iterative solver (PCG)
      ! SUNMatrix
      sunmat_A => null()
      ! Tell SUNDIALS to use a PCG linear solver.
      sunlinsol_LS => FSUNLinSol_PCG(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, PCG'
         stop 1
      end if
      !
    case(3) ! Scaled-Preconditioned BiCGtab iterative Solver (SPBCGS)
      ! SUNMatrix
      sunmat_A => null()
      ! Tell SUNDIALS to use a SPBCGS linear solver
      sunlinsol_LS => FSUNLinSol_SPBCGS(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPBCGS'
         stop 1
      end if
      !
    case(4) ! Scaled-Preconditioned FGMRES iterative solver (SPFGMR)
      ! SUNMatrix
      sunmat_A => null()
      ! Tell SUNDIALS to use a SPFGMR linear solver.
      sunlinsol_LS => FSUNLinSol_SPFGMR(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPFGMR'
         stop 1
      end if
      ! Set the type of Gram-Schmidt orthogonalization to use
      retval = FSUNLinSol_SPFGMRSetGSType(sunlinsol_LS,SUN_MODIFIED_GS)
      if (retval /= 0) then
         print *, 'Error in FSUNLinSol_SPGMRSetGSType'
         stop 1
      end if
    case(5) ! Scaled-Preconditioned GMRES iterative solver
      ! SUNMatrix
      sunmat_A => null()
      ! Tell SUNDIALS to use a SPGMR linear solver.
      sunlinsol_LS => FSUNLinSol_SPGMR(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPGMR'
         stop 1
      end if
      ! Set the type of Gram-Schmidt orthogonalization to use
      retval = FSUNLinSol_SPGMRSetGSType(sunlinsol_LS,SUN_MODIFIED_GS)
      if (retval /= 0) then
         print *, 'Error in FSUNLinSol_SPGMRSetGSType'
         stop 1
      end if
    case(6) ! Scaled-Preconditioned TFQMR iterative solver
      ! SUNMatrix
      sunmat_A => null()
      ! Tell CVODE to use a SPTFQMR linear solver.
      sunlinsol_LS => FSUNLinSol_SPTFQMR(sunvec_y,SUN_PREC_NONE, 5, sunctx)
      if (.not. associated(sunlinsol_LS)) then
         print *, 'ERROR: sunls = NULL, SPGMR'
         stop 1
      end if
    case default
      print*, 'no match for solopt = ', solopt
      stop 'program stopped at program acprpzr'
  end select


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
  tout = tstart
  call print_output_files(tout)
  ! push on to the next step
  tout = tout + dtotp

  do while(tout<=tend)
    retval = FIDASolve(ida_mem, tout, tret, sunvec_y, sunvec_yp, IDA_NORMAL)
    if (retval < 0) then
      print *, 'Error in FIDASolve, retval = ', retval, '; halting'
      stop 1
    end if
    ! output current data
    call print_output_files(tout)
    print*, iout, tout
    if (retval == IDA_SUCCESS) then
      iout = iout + 1
      tout = tout + dtotp
    end if

    if (iout == NOUT) exit

  end do

  ! free memory
  call FIDAFree(ida_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_av)
  call FN_VDestroy(sunvec_yp)
  retval = FSUNContext_Free(sunctx)

  call CPU_TIME(tac)
  print*, 'total CPU time = ', tac - tic
end program acprpzr

