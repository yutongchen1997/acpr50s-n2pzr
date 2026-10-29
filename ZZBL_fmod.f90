! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! The following is the 'main' program for analyzing nitrogen gas 
! pressurizer dynamics.
! 
! ------------------------------------------------------------------
program zzbl

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

  use corrsub_mod
  !======= Declarations =========
  implicit none

  real(c_double) :: pg, hg, hl
  real(c_double) :: tg, tl, rl
  character(len=20) :: eosmod  ! eos model type 

  call read_n2g_dat()
  call read_h2o_dat()

  eosmod = 'nist'
  pg = 15.500D6
  tg = 373.00D0
  tl = 373.00D0
  

  call hgeos_n2g(eosmod ,pg, tg, hg)
  call hleos_h2o(pg,tl,hl)
  !call tleos_h2o(pg, hl, tl)
  !call rleos_h2o(pg, hl, rl)
  
  print*, 'hg = ', hg
  print*, 'hl = ', hl
  !print*, 'tl = ', tl
  !print*, 'rl = ', rl

  


end program zzbl