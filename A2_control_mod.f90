! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines control module, controlling the workflow of 
! the execution of the program.
! 
! ------------------------------------------------------------------
module control_mod
  !======= Inclusions ===========
  use sysdomn_mod

  !======= Declarations =========
  implicit none
  !== SHARED variables ==
  character(len=20) :: id, mat, typ, fcpid, sigid
  !character(len=20) :: mats(maxnlyr) = 'end'
  real(c_double) :: temp0, pres0
  !== UNIQUE variables ==
  ! XTIME variables-------------------------------------------------01
  real(c_double) :: tstart= 0.0d00, tend = 1.d+02, dtotp = 1.d+01
  real(c_double) :: dtini = 1.d-08, dtmin= 1.d-08, dtmax = 1.d-01
  ! XSOLV variables-------------------------------------------------02
  character(len=20) :: solver(nsolvr) = 'end'
  real(c_double) :: rt    = 1.0d-4
  real(c_double) :: at(4) = 1.0d-6
  integer(c_long):: mxstps= 1000000
  integer(c_int) :: solopt= 1
  integer(c_int) :: maxnef= 7, maxcor = 3
  ! XPZRV variables-------------------------------------------------03
  character(len=20) :: id_l, id_g, mat_l, mat_g
  real(c_double) :: dia_i, len_t, level0, tgas0, tliq0
  real(c_double) :: htc_i0
  ! XPIPE variables-------------------------------------------------04
  ! temporarily left blank
  ! XBNDC variables-------------------------------------------------05
  character(len=20) :: sig_p, sig_t
  ! XJUNC variables-------------------------------------------------06
  character(len=20) :: fid_f, fid_t
  ! XCORE variables-------------------------------------------------07
  ! XFUEL variables-------------------------------------------------08
  ! XCLAD variables-------------------------------------------------09
  ! XFROD variables-------------------------------------------------10
  ! XTHBC variables-------------------------------------------------11
  ! temporarily left blank
  ! XHSTR variables-------------------------------------------------12
  ! temporarily left blank
  ! XSGNL variables-------------------------------------------------13
  real(c_double) :: val = 0.d0
  real(c_double) :: tab(2, maxtab) = endval
  ! ----------------------------------------------------------------
  ! == NAMELIST declarations ==
  ! XTIME namelist 01.----------------------------------------------
  NAMELIST /XTIME/ tstart, tend  , dtotp , dtini , dtmin , dtmax
  ! XSOLV namelist 02.----------------------------------------------
  NAMELIST /XSOLV/ solver, rt    , at    , maxnef, maxcor, mxstps, &
                   solopt, eosopt
  ! XPZRV namelist 03.----------------------------------------------
  NAMELIST /XPZRV/ id_l  , id_g  , mat_l , mat_g , dia_i , len_t , &
                   level0, htc_i0, tgas0 , tliq0 , pres0
  ! XPIPE namelist 04.----------------------------------------------
  ! XBNDC namelist 05.----------------------------------------------
  NAMELIST /XBNDC/ id    , mat   , typ   , temp0 , pres0 , sig_p , &
                   sig_t
  ! XJUNC namelist 06.----------------------------------------------
  NAMELIST /XJUNC/ typ   , fid_f , fid_t , sigid
  ! XCORE namelist 07.----------------------------------------------
  ! XFUEL namelist 08.----------------------------------------------
  ! XCLAD namelist 09.----------------------------------------------
  ! XFROD namelist 10.----------------------------------------------
  ! XTHBC namelist 11.----------------------------------------------
  ! XHSTR namelist 12.----------------------------------------------
  ! XSGNL namelist 13.----------------------------------------------
  NAMELIST /XSGNL/ id    , typ   , val   , tab
  !== NAMELIST information ==
  ! array storing namelist titles
  character(len=6):: nmlstr(nnmlst)
  ! array storing namelist counts
  integer(c_int)  :: numnml(nnmlst)
  !== Fluid information ==
  ! array storing fluid show-up counts
  integer(c_int)   :: numfld(nfluids)

  contains
! ----------------------------------------------------------------------
! construct_input: a subroutine that read the input file and store all 
! input data into object type(systdat) :: system_obj
! output - solid material properties
! ----------------------------------------------------------------------
  subroutine construct_input()
    !======= Declarations =========
    implicit none
    !== local variables ==
    character(len=6) :: nml_str
    real(c_double),allocatable :: inp_pzrv(:)
    integer(c_int) :: i, j
    ! the objects
    type(pvessel), allocatable :: pvessel_arr(:)
    type(flubndc), allocatable :: flubndc_arr(:)
    type(flujunc), allocatable :: flujunc_arr(:)
    type(ctlsgnl), allocatable :: ctlsgnl_arr(:)
    !======= Internals ============
    ! assign the namelist titles
    nmlstr(idxtime) = '&XTIME'
    nmlstr(idxsolv) = '&XSOLV'
    nmlstr(idxpzrv) = '&XPZRV'
    nmlstr(idxbndc) = '&XBNDC'
    nmlstr(idxjunc) = '&XJUNC'
    nmlstr(idxsgnl) = '&XSGNL'
    ! open input file
    open(unit = inpfu, status = 'old', file = 'input', form = 'formatted')
    ! initialize the counters
    numnml = 0
    numfld = 0
    ! cycle through the input file to count the hits of each namelist
    do while (.TRUE.)
      read(inpfu, '(1X,A6)', end=100) nml_str
      do i = 1, nnmlst
        if(nml_str == nmlstr(i)) then
          numnml(i) = numnml(i) + 1
        end if
      end do
    end do
100 continue

    ! /XTIME/-----------------------------------------------------------
    if(numnml(idxtime)>=1) then
      rewind(inpfu)
      write(0,*) '===/XTIME/==='
      read(inpfu, nml=XTIME) ! only the first /XTIME/ is read
    else
      print*, 'no namelist /XTIME/ detected'
      stop 'program stop at construct_input'
    end if

    ! /XSOLV/-----------------------------------------------------------
    if(numnml(idxsolv)>=1) then
      rewind(inpfu)
      write(0,*) '===/XSOLV/==='
      read(inpfu, nml=XSOLV) ! only the first /XSOLV/ is read
      ! cycle the solver array
      do i = 1, nsolvr
        select case (solver(i))
          case('fluid')
            solve_fluid  = .TRUE.
            print*, 'solves fluid'
          case('htstr')
            solve_htstr  = .TRUE.
            print*, 'solves htstr'
          case default
            if(solver(i)=='end') then
              print*, 'all involved solvers are presented above...'
              exit
            else
              print*, 'no match for solver = ', solver(i)
              stop 'program stopped at construct_input'
            end if
        end select
      end do
      ! print SUNLinSol information (linear solver information)
      select case (solopt)
        case(1)
          print*, 'Dense direct linear solver (internal) selected'
        case(2)
          print*, 'PCG iterative solver selected'
        case(3)
          print*, 'SPBCGS solver selected'
        case(4)
          print*, 'SPFGMR solver selected'
        case(5)
          print*, 'SPGMR solver selected'
        case(6)
          print*, 'SPTFQMR solver selected'
        case default
          print*, 'no match for solopt = ', solopt
          stop 'program stopped at construct_input'
      end select
    else 
      print*, 'no namelist /XSOLV/ detected'
      stop 'program stop at construct_input'
    end if
    ! start of solve_fluid section =====================================
    if(solve_fluid) then
    ! /XPZRV/-----------------------------------------------------------
      if(numnml(idxpzrv)>=1) then
        allocate(pvessel_arr(numnml(idxpzrv)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XPZRV/=== *', numnml(idxpzrv)
        do i = 1, numnml(idxpzrv)
          read(inpfu, nml=xpzrv)
          allocate(inp_pzrv(7))
          inp_pzrv(1) = dia_i
          inp_pzrv(2) = len_t
          inp_pzrv(3) = level0
          inp_pzrv(4) = tgas0
          inp_pzrv(5) = tliq0
          inp_pzrv(6) = pres0
          inp_pzrv(7) = htc_i0
          call pvessel_arr(i)%init_pvessel(id_l, id_g, mat_l, mat_g, inp_pzrv)
          deallocate(inp_pzrv)
          ! check and count the fluid show-up times
          ! liquid phase
          select case (mat_l)
            case ('h2o')
              numfld(idxh2o) = numfld(idxh2o) + 1
            case ('n2g')
              print*, mat_l, ' is not suitable for liquid in pzrv'
              stop 'program stopped at construct_input'
            case default
              print*, 'no match for liquid = ', mat_l
              stop 'program stopped at construct_input'
          end select
          ! gas phase
          select case (mat_g)
            case ('h2o')
              numfld(idxh2o) = numfld(idxh2o) + 1
            case ('n2g')
              numfld(idxn2g) = numfld(idxn2g) + 1
            case default
              print*, 'no match for gas = ', mat_g
              stop 'program stopped at construct_input'
          end select
        end do
        ! assign the pvessel_arr(:) to system_glo
        call system_glo%fluid%init_pzrv(pvessel_arr)
      else
        print*, 'no namelist /XPRZV/ detected for fluid solver'
        continue
      end if
      deallocate(pvessel_arr)
    ! /XBNDC/-----------------------------------------------------------
      if(numnml(idxbndc)>=1) then
        allocate(flubndc_arr(numnml(idxbndc)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XBNDC/=== *', numnml(idxbndc)
        do i = 1, numnml(idxbndc)
          read(inpfu, nml=xbndc)
          call flubndc_arr(i)%init_flubndc(id, typ, mat, sig_p, sig_t)
          ! check and count the fluid show-up times
          select case (mat)
            case ('h2o')
              numfld(idxh2o) = numfld(idxh2o) + 1
            case ('n2g')
              numfld(idxn2g) = numfld(idxn2g) + 1
            case default
              print*, 'no match for fluid material = ', mat
              stop 'program stopped at construct_input'
          end select
        end do
        ! assign the flubndc_arr(:) to system_glo
        call system_glo%fluid%init_bndcs(flubndc_arr)
      else
        print*, 'no namelist /XBNDC/ detected for fluid solver'
        continue
      end if 
      deallocate(flubndc_arr)
    ! /XJUNC/-----------------------------------------------------------
      if(numnml(idxjunc)>=1) then
        allocate(flujunc_arr(numnml(idxjunc)))
        rewind(inpfu)
        write(0,'(1X,A15,I3)') '===/XJUNC/=== *', numnml(idxjunc)
        do i = 1, numnml(idxjunc)
          read(inpfu, nml=xjunc)
          ! read different types of junctions seperately
          select case (typ)
            case('junc_flow', 'junc_valv', 'junc_pump')
              call flujunc_arr(i)%init_flujunc(typ, fid_f, fid_t, sigid)
            case('junc_norm')
              call flujunc_arr(i)%init_flujunc(typ, fid_f, fid_t)
            case default
              print*, 'no match for junction type = ', typ
              stop 'program stopped at construct_input'
          end select
        end do
        ! assign the flujunc_arr(:) to system_glo
        call system_glo%fluid%init_juncs(flujunc_arr)
      else
        print*, 'no namelist /XJUNC/ detected for fluid solver'
        stop 'program stopped at construct_input'
      end if 
      deallocate(flujunc_arr)
      ! then get the fluid thermophysical data
      if(numfld(idxh2o)>0) call read_h2o_dat()
      if(numfld(idxn2g)>0) call read_n2g_dat()
    end if
    ! end of solve_fluid section =========================================
    ! /XSGNL/-------------------------------------------------------------
    if(numnml(idxsgnl)>=1) then
      allocate(ctlsgnl_arr(numnml(idxsgnl)))
      allocate(system_glo%signals(numnml(idxsgnl)))
      rewind(inpfu)
      write(0,'(1X,A15,I3)') '===/XSGNL/=== *', numnml(idxsgnl)
      do i = 1, numnml(idxsgnl)
        tab = endval
        read(inpfu, nml=xsgnl)
        ! echo the signals with fixed name and purpose
        select case (id)
          case('dtprint')
            print*, '*** system signal "dtprint" specified ***'
            set_dtprint = .TRUE.
          case default
            print*, ADJUSTL(TRIM(id)), ' signal specified'
        end select
        ! store the signal data
        select case (typ)
          case('const') ! constant
            call ctlsgnl_arr(i)%init_sgnl(id, typ)
            ctlsgnl_arr(i)%val = val
          case('tdsig') ! time-dependent signal
            ! 'j' is the counter of actual input number of sgntab
            j = 1
            do while (.TRUE.)
              if(tab(1,j)/=endval .and. tab(2,j)/=endval) then
                j = j+1
              else
                if(tab(1,j)/=endval .or. tab(2,j)/=endval) then
                  print*, 'Warning in signal id = ', id, ', time-dependent table'
                  print*, 'the x-array and y-array sizes do not match'
                end if
                exit
              end if 
            end do
            ! check the input table length
            if(j<3) then
              print*, 'data table needs at least 2 columns for signal = ', id
              stop 'program stopped at construct_inputs'
            else
              call ctlsgnl_arr(i)%init_sgnl(id, typ, tab(1:2,1:j-1))
              ctlsgnl_arr(i)%val = tab(2,1)
            end if
          case('lutab') ! look-up table signal
            print*, 'now under development, please wait with patience ...'
            stop 'program stopped at construct_input'
          case default
            print*, 'no match for signal type = ', typ
            stop 'program stopped at construct_input'
        end select
      end do
      system_glo%signals = ctlsgnl_arr
      system_glo%numsgnl = size(ctlsgnl_arr)
      deallocate(ctlsgnl_arr)
    else
      print*, 'no namelist /XSGNL/ detected'
      continue
    end if
    ! close input file
    close(inpfu)

  end subroutine construct_input

! ----------------------------------------------------------------------
! evaluate_signals: a subroutine that evaluate all signals, which means
! updating and patching the signal values during transient calculation
! input - time
! ---------------------------------------------------------------------
  subroutine evaluate_signals(time)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    real(c_double), intent(in) :: time
    !== local variables ==
    real(c_double) :: value       ! signal value
    real(c_double) :: pbnd, tbnd  ! boundary cell p, t
    integer(c_int) :: i, j
    logical :: pfound, tfound, jfound


    !======= Internals ============
    do i = 1, system_glo%numsgnl
      ! update the signal values array
      select case (system_glo%signals(i)%typ)
        case('const')
          continue  ! no further action for constant signal
        case('tdsig')
          call interp1d(system_glo%signals(i)%tab(1,:), &
                        system_glo%signals(i)%tab(2,:), time, value)
          system_glo%signals(i)%val = value
        case('lutab')
          print*, 'now under development...'
          stop 'program stopped at evaluate_signals'
        case default
          print*, 'no match for signal type = ', system_glo%signals(i)%typ
          stop 'program stopped at evaluate_signals'
      end select
    end do
    ! patch the signal values relative to fluid domain
    if(solve_fluid) then
      ! update the fluid boundary cells
      do i = 1, system_glo%fluid%nbndc
        pfound = .FALSE. ! signal for p is found?
        tfound = .FALSE. ! signal for t is found?
        do j = 1, system_glo%numsgnl
          if(system_glo%fluid%bndcs(i)%sig_p==system_glo%signals(j)%id) then
            pbnd = system_glo%signals(j)%val
            pfound = .TRUE.
          end if
          if(system_glo%fluid%bndcs(i)%sig_t==system_glo%signals(j)%id) then
            tbnd = system_glo%signals(j)%val
            tfound = .TRUE.
          end if
          if(pfound .and. tfound) exit
        end do
        !
        if(pfound .and. tfound) then
          call system_glo%fluid%bndcs(i)%updt_flubndc(pbnd, tbnd)
        else
          print*, 'lacking p/t signal for bnd cell ', system_glo%fluid%bndcs(i)%id
          stop 'program stopped at evaluate_signals'
        end if
      end do
      ! update the fluid junctions
      do i = 1, system_glo%fluid%njunc
        jfound = .FALSE. ! signal for junction is found?
        do j = 1, system_glo%numsgnl
          if(system_glo%fluid%juncs(i)%sig==system_glo%signals(j)%id) then
            select case (system_glo%fluid%juncs(i)%typ)
              case ('jun-flow')  ! mass-flow-rate specified
                system_glo%fluid%juncs(i)%mdot = system_glo%signals(j)%val
              case ('jun-valv')  ! k-factor specified 
                system_glo%fluid%juncs(i)%kfac = system_glo%signals(j)%val
              case ('jun-pump')  ! pump-head specified
                system_glo%fluid%juncs(i)%phed = system_glo%signals(j)%val
              case default
                if(system_glo%fluid%juncs(i)%typ=='jun-norm') then
                  print*, 'warning, signal patched to normal junction would be neglected ...'
                  print*, system_glo%fluid%juncs(i)%fid_f,' -> ',system_glo%fluid%juncs(i)%fid_t
                  continue
                else
                  print*, 'no match for junction type = ', system_glo%fluid%juncs(i)%typ
                  stop 'program stopped at evaluate_signals'
                end if
            end select
            jfound = .TRUE.
            exit
          end if
        end do
      end do
    end if
    ! patch the signal values relative to solid domain
    ! now under development...
  end subroutine evaluate_signals

! ----------------------------------------------------------------------
! open_output_files: a subroutine that open the output files, assigning 
! the file units and name strings
! ----------------------------------------------------------------------
  subroutine open_output_files()
    !======= Declarations =========
    implicit none
    !== local variables ==
    ! output directory name
    character(len=20) :: outdir
    ! path for outputs
    character(len=50) :: path4outdir
    ! file id string
    character(len=250):: fileid
    ! column id string
    character(len=18) :: colid
    ! system command string
    character(len=100):: syscmd
    ! loop counters
    integer(c_int) :: i, j, k
    ! status variable
    integer(c_int) :: istatus
    ! integers recording system time
    integer(c_int) :: time_arr(8)
    integer(c_int) :: year, month, day, hour, minute, second

    !======= Internals ============




  end subroutine open_output_files

! ----------------------------------------------------------------------
! print_output_files: a subroutine that write calculation results to 
! output files every output time step
! ----------------------------------------------------------------------
  subroutine print_output_files(time)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    real(c_double), intent(in) :: time
    !== local variables ==
    character(len=18) :: colstr
    integer(c_int) :: i, j, k

    !======= Internals ============


  end subroutine print_output_files







end module control_mod
