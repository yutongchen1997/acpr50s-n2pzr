! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: sldfrod (fuel rod)
! 
! ------------------------------------------------------------------
module sldfrod_mod
  !======= Inclusions =========
  use corrsub_mod

  !======= Declarations =========
  implicit none

  type :: sldfrod
    !=== Member Variables ===
    ! required input parameters ------------------------------------
    ! id string of fuel rod, fuel pellet, gas gap, cladding 
    character(len=20) :: id, id_f, id_g, id_c
    ! id string of coupling fluid pipe component
    character(len=20) :: id_pipe
    ! coupling fluid pipe node, multiplier (fuel rod number)
    integer(c_int) :: nz, mltpl
    ! fuel pellet
    type(fuel) :: fuel
    ! gasgap
    type(ggap) :: ggap
    ! cladding
    type(clad) :: clad
  contains
    !=== Member Procedures ===
    procedure :: init_sldfrod
    procedure :: calc_res_sldfrod

  end type sldfrod

  
contains


end module sldfrod_mod