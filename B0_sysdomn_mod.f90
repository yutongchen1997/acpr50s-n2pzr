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
    procedure :: write_to_y_and_yp
    procedure :: read_from_y_and_yp

  end type sysdomn
  ! build the system object
  ! system_glo is the sysdomn-type object that stores all system data
  type(sysdomn) :: system_glo

contains
! ----------------------------------------------------------------------
! write_to_y: write unknowns to yval-vector for initialization
! ----------------------------------------------------------------------
  subroutine write_to_y_and_yp(self, yval, ypval)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(sysdomn), intent(inout) :: self
    real(c_double), intent(out)   :: yval(:)
    real(c_double), intent(out)   :: ypval(:)
    !== local variables ==
    integer(c_int) :: i, j, k

    !======= Internals ============
    k = 1
    ! fluid solver
    if(solve_fluid) then
      ! pvessel data----------------------------------------------------
      do i = 1, self%fluid%npzrv
        ! fluid masses to be solved by sundials
        yval(k) = self%fluid%pzrvs(i)%m_l  ! 01
        k = k + 1
        yval(k) = self%fluid%pzrvs(i)%m_g  ! 02
        k = k + 1
        ! fluid total ENTHALPIES to be solved by sundials
        yval(k) = self%fluid%pzrvs(i)%h_l  ! 03
        k = k + 1
        yval(k) = self%fluid%pzrvs(i)%h_g  ! 04
        k = k + 1
        ! fluid pressure to be solved by sundials
        yval(k) = self%fluid%pzrvs(i)%p_g  ! 05
        k = k + 1
        ! fluid  temperatures to be solved by sundials
        yval(k) = self%fluid%pzrvs(i)%t_l  ! 06
        k = k + 1
        yval(k) = self%fluid%pzrvs(i)%t_g  ! 07
        k = k + 1
        ! fluid volumes to be solved by sundials
        yval(k) = self%fluid%pzrvs(i)%v_l  ! 08
        k = k + 1
        yval(k) = self%fluid%pzrvs(i)%v_g  ! 09
        k = k + 1
      end do
      ! external junctions----------------------------------------------
      do i = 1, self%fluid%njunc
        ! mdots to be solved by sundials
        if(self%fluid%juncs(i)%typ=='junc_flow') then
          continue
        else
          yval(k) = self%fluid%juncs(i)%mdot
          k = k + 1
        end if
      end do
      ! ypvals----------------------------------------------------------
      ! initialize vector of time derivatives to zero
      ypval = 0.d0
    end if

  end subroutine write_to_y_and_yp

! ----------------------------------------------------------------------
! read_from_y: read unknowns from y-vector
! ----------------------------------------------------------------------
  subroutine read_from_y_and_yp(self, yval, ypval)
    !======= Declarations =========
    implicit none
    !== I/O variables ==
    class(sysdomn), intent(inout) :: self
    real(c_double), intent(in)    :: yval(:)
    real(c_double), intent(in)    :: ypval(:)
    !== local variables ==
    integer(c_int) :: i, j, k

    !======= Internals ============
    k = 1
    ! fluid solver
    if(solve_fluid) then
      ! pvessel data----------------------------------------------------
      do i = 1, self%fluid%npzrv
        ! fluid masses to be solved by sundials
        self%fluid%pzrvs(i)%m_l  = yval(k)
        self%fluid%pzrvs(i)%dmldt= ypval(k)
        k = k + 1
        self%fluid%pzrvs(i)%m_g  = yval(k)
        self%fluid%pzrvs(i)%dmgdt= ypval(k)
        k = k + 1
        ! fluid total ENTHALPIES to be solved by sundials
        self%fluid%pzrvs(i)%h_l  = yval(k)
        self%fluid%pzrvs(i)%dhldt= ypval(k)
        k = k + 1
        self%fluid%pzrvs(i)%h_g  = yval(k)
        self%fluid%pzrvs(i)%dhgdt= ypval(k)
        k = k + 1
        ! fluid pressure to be solved by sundials
        self%fluid%pzrvs(i)%p_g = yval(k)
        self%fluid%pzrvs(i)%dpdt= ypval(k)
        k = k + 1
        ! fluid  temperatures to be solved by sundials
        self%fluid%pzrvs(i)%t_l  = yval(k)
        self%fluid%pzrvs(i)%dtldt= ypval(k)
        k = k + 1
        self%fluid%pzrvs(i)%t_g  = yval(k)
        self%fluid%pzrvs(i)%dtgdt= ypval(k)
        k = k + 1
        ! fluid volumes to be solved by sundials
        self%fluid%pzrvs(i)%v_l  = yval(k)
        self%fluid%pzrvs(i)%dvldt= ypval(k)
        k = k + 1
        self%fluid%pzrvs(i)%v_g  = yval(k)
        self%fluid%pzrvs(i)%dvgdt= ypval(k)
        k = k + 1
        ! update liquid level
        call self%fluid%pzrvs(i)%get_liq_lvl()
      end do
      ! external fluid junctions----------------------------------------
      do i = 1, self%fluid%njunc
        ! mdots to be solved by sundials
        if(self%fluid%juncs(i)%typ=='junc_flow') then
          continue
        else
          self%fluid%juncs(i)%mdot   = yval(k)
          self%fluid%juncs(i)%dmdotdt= ypval(k)
          k = k + 1
        end if
      end do
    end if

  end subroutine read_from_y_and_yp


end module sysdomn_mod