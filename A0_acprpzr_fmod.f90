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
  real(c_double) :: tout, tout1, tret(1)
  ! CPU time variables
  real(c_double) :: tic, tac
  integer(c_int) :: iout, retval, ic_opt
  ! integers
  integer, parameter :: myindextype = selected_int_kind(16)
  integer(kind=myindextype) :: neq = 0
  integer(c_int) :: i, j, k
  ! SUNDIALS types
  type(N_Vector), pointer :: sunvec_y      ! sundials solution vector
  type(N_Vector), pointer :: sunvec_yp     ! sundials derivative vector
  type(N_Vector), pointer :: sunvec_av     ! sundials tolerance vector
  type(N_Vector), pointer :: sunvec_id     ! sundials variable-ID vector (for FIDACalcIC)
  type(SUNMatrix), pointer :: sunmat_A     ! sundials matrix
  type(SUNLinearSolver), pointer :: sunlinsol_LS  ! sundials linear solver
  type(SUNNonLinearSolver), pointer :: sunnonlin_NLS ! sundials nonlinear solver
  type(c_ptr)                       :: ida_mem       ! IDA memory
  type(c_ptr)                       :: sunctx        ! SUNDIALS simulation context

  ! solution and derivative vectors, neq is set in the dae_mod module
  real(c_double), allocatable :: yval(:), ypval(:)
  real(c_double), allocatable :: id_vals(:)  ! variable-ID vector for FIDACalcIC

  !======= Internals ============
  call CPU_TIME(tic)

  retval = FSUNContext_Create(SUN_COMM_NULL, sunctx)

  ! open the input file and read input data
  call construct_input()
  ! open output files
  call open_output_files()
  ! echo the input file into the output directory
  call echo_input()
  ! sleep for 1 seconds
  call SLEEP(1)
  ! convert neq type
  neq = int(neq_tot, kind = myindextype)

  ! set tolerance vectors
  rtol = rt
  !
  allocate(yval(neq))
  allocate(ypval(neq))
  allocate(avtol(neq))
  allocate(id_vals(neq))
  ! set avtols
  i = 1
  if(neq_mdot>0) then
  ! junction mass flowrate [kg/s]
    avtol(i:i+neq_mdot-1) = at(1)
  end if
  i = 1 + neq_mdot
  if(neq_pzrv>0) then
    do j = 1, system_glo%fluid%npzrv
      ! fluid mass        [kg]
      avtol(i+9*j-9:i+9*j-8)= at(1)
      ! fluid energy      [J]
      avtol(i+9*j-7:i+9*j-6)= at(2)
      ! fluid pressure    [Pa]
      avtol(i+9*j-5)        = at(3)
      ! fluid temperature [K]
      avtol(i+9*j-4:i+9*j-3)= at(4)
      ! fluid volume      [m^3]
      avtol(i+9*j-2:i+9*j-1)= at(5)
    end do
  end if
  i = 1 + neq_mdot + neq_pzrv
  if(neq_brch>0) then
    do j = 1, system_glo%fluid%nbrch
      ! fluid mass        [kg]
      avtol(i+4*j-4) = at(1)
      ! fluid energy      [J]
      avtol(i+4*j-3) = at(2)
      ! fluid temperature [K]
      avtol(i+4*j-2) = at(4)
      ! fluid pressure    [Pa]
      avtol(i+4*j-1) = at(3)
    end do
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
  !
  ! User-configurable, default (500)
  retval = FIDASetMaxNumSteps(ida_mem, mxstps)
  if (retval /= 0) then
    print *, 'Error in FIDASetMaxNumSteps, retval = ', retval, '; halting'
    stop 1
  end if
  !
  ! User-configurable, default (10)
  retval = FIDASetMaxErrTestFails(ida_mem, maxnef)
  if (retval /= 0) then
    print *, 'Error in FIDASetMaxErrTestFails, retval = ', retval, '; halting'
    stop 1
  end if
  !
  ! User-configurable, default (10)
  retval = FIDASetMaxNonlinIters(ida_mem, maxcor)
  if (retval /= 0) then
    print *, 'Error in FIDASetMaxNonlinIters, retval = ', retval, '; halting'
    stop 1
  end if
  !
  ! set initial timestep
  retval = FIDASetInitStep(ida_mem, dtini)
  if (retval /= 0) then
    print *, 'Error in FIDASetInitStep, retval = ', retval, '; halting'
    stop 1
  end if
  ! set maximum timestep
  retval = FIDASetMaxStep(ida_mem, dtmax)
  if (retval /= 0) then
    print *, 'Error in FIDASetMaxStep, retval = ', retval, '; halting'
    stop 1
  end if
  ! set minimum timestep
  retval = FIDASetMinStep(ida_mem, dtmin)
  if (retval /= 0) then
    print *, 'Error in FIDASetMinStep, retval = ', retval, '; halting'
    stop 1
  end if

  !!!  YTCHEN - 2026/11/08 - 0025
  !  id_vals( 1) = 0
  !  id_vals( 2) = 1
  !  id_vals( 3) = 1
  !  id_vals( 4) = 1
  !  id_vals( 5) = 1
  !  id_vals( 6) = 0
  !  id_vals( 7) = 0
  !  id_vals( 8) = 0
  !  id_vals( 9) = 0
  !  id_vals(10) = 0
  !  id_vals(11) = 1
  !  id_vals(12) = 1
  !  id_vals(13) = 0
  !  id_vals(14) = 0
  !  
  !  sunvec_id => FN_VMake_Serial(neq, id_vals, sunctx)
  !  if (.not. associated(sunvec_id)) then
  !     print *, 'ERROR: sunvec_id = NULL'
  !     stop
  !  end if
  !  retval = FIDASetId(ida_mem, sunvec_id)
  !  if (retval /= 0) then
  !    print *, 'Error in FIDASetId, retval = ', retval
  !    stop
  !  end if
  !  retval = FIDASetSuppressAlg(ida_mem, 1_c_int)
  !   if (retval /= 0) then
  !      print *, 'Warning: FIDASetSuppressAlg returned retval = ', retval
  !   end if
  !  
  !  ic_opt = IDA_YA_YDP_INIT              ! y 已知，ydot 待求
  !  tout1  = tstart + 1.0d-3              ! 一致化目标时刻
  !  retval = FIDACalcIC(ida_mem, ic_opt, tout1)
  !  if (retval /= 0) then
  !    print *, 'FIDACalcIC failed, retval = ', retval; stop
  !  else
  !    print *, 'FIDACalcIC succeeded: initial derivative vector calculated.'
  !  end if
  !  print*, 'ypval1', ypval
  !  retval = FIDAGetConsistentIC(ida_mem, sunvec_y, sunvec_yp)
  !  if (retval /= 0) then
  !    print *, 'Warning: IDAGetConsistentIC returned retval = ', retval
  !  else
  !    print *, 'FIDAGetConsistentIC succeeded.'
  !  end if
  !  print*, 'ypval2', ypval
  !  call FN_VDestroy(sunvec_id)

  !
  iout = 0
  tout = tstart
  !
  ! ### KONSTANTIN MIKITYUK - START
  ! Always set the id vector and supress algebraic vars from LTE test.
  !!!!!!!!call set_id_vector(id_vals)
  !!!!!!!!sunvec_id => FN_VMake_Serial(neq, id_vals, sunctx)
  !!!!!!!!if (.not. associated(sunvec_id)) then
  !!!!!!!!   print *, 'ERROR: sunvec_id = NULL'
  !!!!!!!!   stop
  !!!!!!!!end if
  !!!!!!!!retval = FIDASetId(ida_mem, sunvec_id)
  !!!!!!!!if (retval /= 0) then
  !!!!!!!!   print *, 'Error in FIDASetId, retval = ', retval
  !!!!!!!!   stop
  !!!!!!!!end if
  !!!!!!!!retval = FIDASetSuppressAlg(ida_mem, 1_c_int)
  !!!!!!!!if (retval /= 0) then
  !!!!!!!!   print *, 'Warning: IDASetSuppressAlg returned retval = ', retval
  !!!!!!!!end if
  ! IDA_YA_YDP_INIT: 
  ! Get the consistent initial values of y(t0) and y'(t0) vectors
  !!!!!!!!print*, 'Re-calculating initial derivative vectors via FIDACalcIC...'
  !!!!!!!!retval = FIDACalcIC(ida_mem, IDA_YA_YDP_INIT, tout + dtini)
  ! ### KONSTANTIN MIKITYUK - END
  !
  call print_output_files(tout)
  write(*,'(A,I0,A,A,G0)') 'iout = ',iout,'    ','tout = ',tout
  ! push on to the next step
  iout = 1
  tout = tout + dtotp
  !
  ! In loop, call IDASolve, print results, and test for error.
  ! Break out of loop when nout preset output times have been reached.
  do while(tout<=tend)
    retval = FIDASolve(ida_mem, tout, tret, sunvec_y, sunvec_yp, IDA_NORMAL)
    if (retval < 0) then
      print *, 'Error in FIDASolve, retval = ', retval, '; halting'
      call print_output_files(tret(1))
      stop 1
    end if
    ! output current data
    call print_output_files(tout)
    write(*,'(A,I0,A,A,G0)') 'iout = ',iout,'    ','tout = ',tout
    if (retval == IDA_SUCCESS) then
      iout = iout + 1
      tout = tout + dtotp
    end if
    !
    if (iout == nout) exit
    !
  end do
  ! close output files
  call close_output_files()
  ! free memory
  call FIDAFree(ida_mem)
  retval = FSUNNonlinSolFree(sunnonlin_NLS)
  retval = FSUNLinSolFree(sunlinsol_LS)
  call FSUNMatDestroy(sunmat_A)
  call FN_VDestroy(sunvec_y)
  call FN_VDestroy(sunvec_av)
  call FN_VDestroy(sunvec_yp)
  retval = FSUNContext_Free(sunctx)
  !
  call CPU_TIME(tac)
  print*, 'total CPU time = ', tac - tic

end program acprpzr

