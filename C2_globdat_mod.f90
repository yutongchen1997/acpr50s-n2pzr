! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module contains all global variables for the program.
! 
! ------------------------------------------------------------------

module globdat_mod
  !======= Inclusions ===========
  use, intrinsic :: iso_c_binding

  !======= Declarations =========
  implicit none
  !== Constant Parameters ==
  ! scientific constants
  real(c_double), parameter :: pi     = 3.141592653589793d0
  real(c_double), parameter :: grav   = 9.8d0              ! gravity acceleration
  real(c_double), parameter :: r_u    = 8.31446261815324d0 ! universal gas constant
  real(c_double), parameter :: gam_n2 = 1.4d0              ! gama coefficient of n2
  real(c_double), parameter :: mol_n2 = 28.0134d-3         ! molecular weight of n2
  real(c_double), parameter :: small  = 1.0d-10
  real(c_double), parameter :: minor  = 1.0d-20
  real(c_double), parameter :: large  = 1.0d+10
  real(c_double), parameter :: huge   = 1.0d+20
  real(c_double), parameter :: endval = 1145141919810.d0   ! end marker
  ! array sizes
  integer(c_int), parameter :: maxtab = 100    ! maximum look-up table length
  integer(c_int), parameter :: maxnlyr= 10     ! maximum number of htstr layers
  integer(c_int), parameter :: nnmlst = 13     ! total number of input file namelists
  integer(c_int), parameter :: nsolvr = 4      ! total number of available solvers
  integer(c_int), parameter :: nfluids= 2      ! total number of available fluids
  ! integers relative to file units
  integer(c_int), parameter :: inpfu  = 5      ! input file unit
  integer(c_int), parameter :: otpfu  = 10     ! output file unit baseline
  integer(c_int), parameter :: n2gfu  = 10001  ! N2 gas property data file unit
  integer(c_int), parameter :: h2ofu  = 10005  ! H2O property data file unit
  integer(c_int), parameter :: h2osu  = 10006  ! H2O saturated data file unit
  ! indexes of namelist titles at array nmlstr(:)
  integer(c_int), parameter :: idxtime= 1      ! index of &XTIME - Time
  integer(c_int), parameter :: idxsolv= 2      ! index of &XSOLV - Solver
  integer(c_int), parameter :: idxpzrv= 3      ! index of &XPZRV - Pressurizer vessel
  integer(c_int), parameter :: idxpipe= 4      ! index of &XPIPE - Pipe
  integer(c_int), parameter :: idxbndc= 5      ! index of &XBNDC - Fluid boundary cell
  integer(c_int), parameter :: idxjunc= 6      ! index of &XJUNC - Junction
  integer(c_int), parameter :: idxcore= 7      ! index of &XCORE - Reactor core
  integer(c_int), parameter :: idxfuel= 8      ! index of &XFUEL - Fuel pellet
  integer(c_int), parameter :: idxclad= 9      ! index of &XCLAD - Cladding
  integer(c_int), parameter :: idxfrod= 10     ! index of &XFROD - Fuel rod
  integer(c_int), parameter :: idxthbc= 11     ! index of &XTHBC - Thermal B.C.
  integer(c_int), parameter :: idxhstr= 12     ! index of &XHSTR - Heat Structure
  integer(c_int), parameter :: idxsgnl= 13     ! index of &XSGNL - Signal
  ! indexes of fluid name at fluid counter array numflu(:)
  integer(c_int), parameter :: idxn2g = 1      ! index of N2 fluid (gas phase)
  integer(c_int), parameter :: idxh2o = 2      ! index of H2O fluid
  !
  !== Global Variables ==
  ! logical variables
  logical :: solve_htstr  = .FALSE.
  logical :: solve_fluid  = .FALSE.
  logical :: set_dtprint  = .FALSE.
  ! N2 gas thermal property data
  real(c_double), allocatable :: p_n2g(:)
  real(c_double), allocatable :: h_n2g(:)
  real(c_double), allocatable :: tpf_n2g(:,:,:)
  ! H2O thermal property data
  real(c_double), allocatable :: p_h2o(:)
  real(c_double), allocatable :: h_h2o(:)
  real(c_double), allocatable :: tpf_h2o(:,:,:)
  real(c_double), allocatable :: sat_h2o(:,:)
  ! fluid index data
  integer(c_int), allocatable :: findindex(:,:) ! findindex(j,k) means the j-th pipe, k-th node's 'i' in 1-d array
  integer(c_int) :: pindex   ! index of the first pressure equation in Bx=b
  integer(c_int) :: neq_mdot = 0 !
  integer(c_int) :: neq_fvol = 0 ! fluid volume
  integer(c_int) :: neq_enth = 0 ! fluid enthalpy
  ! solid index data
  integer(c_int) :: neq_hstr = 0 ! number of ode equations of heat structure node temperature
  integer(c_int) :: neq_frod = 0 ! number of ode equations of fuelrod node temperature
  ! total number of DAEs 
  integer(c_int) :: neq_tot = 0
  ! selection of gas-phase equation-of-state model
  character(len=20) :: eosopt = 'igas'
  ! parallel calculation data
  integer(c_int) :: num_of_threads = 1 ! number of omp threads
  ! tolerances
  real(c_double) :: rtol
  real(c_double), allocatable :: avtol(:)
  
end module globdat_mod