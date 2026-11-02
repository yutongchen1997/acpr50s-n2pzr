! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: fludomn (fluid domain).
! 
! ------------------------------------------------------------------
module fludomn_mod
  !======= Inclusions =========
  use pvessel_mod
  use flujunc_mod
  use flubndc_mod

  !======= Declarations =========
  implicit none

  type :: fludomn
    !=== Member Variables ===
    ! number of boundary cells, external junctions
    integer(c_int) :: nbndc, njunc, npzrv
    ! pressurizer vessel object
    type(pvessel), allocatable :: pzrvs(:)
    ! array storing boundary cell objects
    type(flubndc), allocatable :: bndcs(:)
    ! array storing junction objects
    type(flujunc), allocatable :: juncs(:)
    ! array storing junction connection information
    character(len=20), allocatable:: f_ids(:), t_ids(:)
    character(len=20), allocatable:: f_typ(:), t_typ(:)
    integer(c_int), allocatable   :: f_idx(:), t_idx(:)

    
  contains
    !=== Member Procedures ===
    procedure :: init_bndcs
    procedure :: init_juncs
    procedure :: init_pzrv
    procedure :: init_fludomn
    procedure :: calc_res_fludomn

  end type fludomn

contains
  ! initialize fluid boundary cells
  subroutine init_bndcs(self, bndc_arr)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(fludomn), intent(inout) :: self
    class(flubndc), intent(in)    :: bndc_arr(:)
    !=== local variables ===

    !======= Internals ============
    ! allocate and assign boundary cell array
    self%nbndc = size(bndc_arr)
    allocate(self%bndcs(self%nbndc))
    self%bndcs = bndc_arr
  end subroutine init_bndcs

  ! initialize junctions
  subroutine init_juncs(self, junc_arr)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(fludomn), intent(inout) :: self
    class(flujunc), intent(in)    :: junc_arr(:)
    !=== local variables ===
    integer(c_int) :: i
    !======= Internals ============
    ! allocate and assign mdot and junction array
    self%njunc = size(junc_arr)
    allocate(self%juncs(self%njunc))
    self%juncs = junc_arr
    do i = 1, self%njunc
      if(self%juncs(i)%typ=='junc_flow') then
        continue
      else
        neq_mdot = neq_mdot + 1
      end if
    end do
  end subroutine init_juncs

  ! initialize pressurizer vessel
  subroutine init_pzrv(self, pzrv_arr)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(fludomn), intent(inout) :: self
    class(pvessel), intent(in)    :: pzrv_arr(:)
    !=== local variables ===

    !======= Internals ============
    ! allocate pressurizer vessel object
    self%npzrv = size(pzrv_arr)
    allocate(self%pzrvs(self%npzrv))
    self%pzrvs = pzrv_arr
    ! number of pzrv equations
    neq_pzrv = 9*self%npzrv ! 9 equations for each pvessel
  end subroutine init_pzrv

  ! initialize fluid domain
  subroutine init_fludomn(self)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(fludomn), intent(inout) :: self
    !=== local variables ===
    integer(c_int) :: i, j
    logical :: junc_f, junc_t
    !======= Internals ============
    ! allocate the arrays storing junction connection information
    allocate(self%f_ids(self%njunc))
    allocate(self%t_ids(self%njunc))
    allocate(self%f_typ(self%njunc))
    allocate(self%t_typ(self%njunc))
    allocate(self%f_idx(self%njunc))
    allocate(self%t_idx(self%njunc))
    ! check if there is a unconnected junction
    do j = 1, self%njunc
      junc_f = .FALSE.
      junc_t = .FALSE.
      ! if connected with fbcs
      do i = 1, self%nbndc  
        if(self%juncs(j)%fid_f==self%bndcs(i)%id) then
          junc_f = .TRUE.
          self%f_ids(j) = self%bndcs(i)%id
          self%f_typ(j) = 'bndc'
          self%f_idx(j) = i
        else if(self%juncs(j)%fid_t==self%bndcs(i)%id) then
          junc_t = .TRUE.
          self%t_ids(j) = self%bndcs(i)%id
          self%t_typ(j) = 'bndc'
          self%t_idx(j) = i
        end if
      end do
      ! if connected with pvessel
      if(.NOT.(junc_f)) then
        do i = 1, self%npzrv
          if(self%juncs(j)%fid_f==self%pzrvs(i)%id_g) then
            junc_f = .TRUE.
            self%f_ids(j) = self%pzrvs(i)%id_g
            self%f_typ(j) = 'pzrv-g'
            self%f_idx(j) = i
          else if(self%juncs(j)%fid_f==self%pzrvs(i)%id_l) then
            junc_f = .TRUE.
            self%f_ids(j) = self%pzrvs(i)%id_l
            self%f_typ(j) = 'pzrv-l'
            self%f_idx(j) = i
          end if
        end do
      end if 
      if(.not.(junc_t)) then
        do i = 1, self%npzrv
          if(self%juncs(j)%fid_t==self%pzrvs(i)%id_g) then
            junc_t = .TRUE.
            self%t_ids(j) = self%pzrvs(i)%id_g
            self%t_typ(j) = 'pzrv-g'
            self%t_idx(j) = i
          else if(self%juncs(j)%fid_t==self%pzrvs(i)%id_l) then
            junc_t = .TRUE.
            self%t_ids(j) = self%pzrvs(i)%id_l
            self%t_typ(j) = 'pzrv-l'
            self%t_idx(j) = i
          end if
        end do
      end if 
      ! if found unconnected junctions
      if((.not.junc_f).or.(.not.junc_t)) then
        print*, 'unconnected j-th junction detected, j = ', j
        print*, 'junciton id = ', self%juncs(j)%id
        stop 'program stopped at init_fludomn'
      end if
    end do

    ! check if the type of connected fbcs and junctions match
    do i = 1, self%nbndc
      select case (self%bndcs(i)%typ)
        case('w_in')  ! mass-flow-inlet fbc
          do j = 1, self%njunc
            if(self%juncs(j)%fid_t==self%bndcs(i)%id) then
              print*, 'error in fbc ', self%bndcs(i)%id
              print*, 'the "w_in" bc cannot serve as "to" of juncs'
              stop 'program stopped at init_fludomn'
            else if(self%juncs(j)%fid_f==self%bndcs(i)%id) then
              select case(self%juncs(j)%typ)
                case('junc_flow')
                  continue
                case default
                  print*, 'the "w_in" bc must serve as "from" of junc_flow'
                  stop 'program stopped at init_fludomn'
              end select
            end if 
          end do
        case('p_in')  ! pressure-inlet fbc
          do j = 1, self%njunc
            if(self%juncs(j)%fid_t==self%bndcs(i)%id) then
              print*, 'error in fbc ', self%bndcs(i)%id
              print*, 'the "p_in" bc cannot serve as "to" of juncs'
              stop 'program stopped at init_fludomn'
            else if(self%juncs(j)%fid_f==self%bndcs(i)%id) then
              select case(self%juncs(j)%typ)
                case('junc_norm','junc_valv')
                  continue
                case default
                  print*, 'the "p_in" bc must serve as "from" of junc_norm or junc_valv'
                  stop 'program stopped at init_fludomn'
              end select
            end if
          end do
        case('w_out') ! mass-flow-outlet fbc
          do j = 1, self%njunc
            if(self%juncs(j)%fid_f==self%bndcs(i)%id) then
              print*, 'error in fbc ', self%bndcs(i)%id
              print*, 'the "w_out" bc cannot serve as "from" of juncs'
              stop 'program stopped at init_fludomn'
            else if(self%juncs(j)%fid_t==self%bndcs(i)%id) then
              select case(self%juncs(j)%typ)
                case('junc_flow')
                  continue
                case default
                  print*, 'the "w_out" bc must serve as "from" of junc_flow'
                  stop 'program stopped at init_fludomn'
              end select
            end if
          end do
        case('p_out') ! pressure-outlet fbc
          do j = 1, self%njunc
            if(self%juncs(j)%fid_f==self%bndcs(i)%id) then
              print*, 'error in fbc ', self%bndcs(i)%id
              print*, 'the "p_out" bc cannot serve as "from" of juncs'
              stop 'program stopped at init_fludomn'
            else if(self%juncs(j)%fid_t==self%bndcs(i)%id) then
              select case(self%juncs(j)%typ)
                case('junc_norm','junc_valv')
                  continue
                case default
                  print*, 'the "p_out" bc must serve as "to" of junc_norm or junc_valv'
                  stop 'program stopped at init_fludomn'
              end select
            end if
          end do
        case default
          continue
      end select
    end do

  end subroutine init_fludomn

  subroutine calc_res_fludomn(self,res)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(fludomn), intent(inout) :: self
    ! residual array
    real(c_double), intent(out) :: res(:)
    !=== local variables ===
    real(c_double) :: mdot
    real(c_double) :: hg, hl, mg, ml
    integer(c_int) :: i, j, k
    !======= Internals ============
    ! calculate the mass and energy transport between fluid components
    do i = 1, self%npzrv
      ! initialization of pvessels inlets/outlets
      self%pzrvs(i)%wl_su  = 0.d0
      self%pzrvs(i)%wg_in  = 0.d0
      self%pzrvs(i)%wg_ex  = 0.d0
      self%pzrvs(i)%wlhl_su= 0.d0
      self%pzrvs(i)%wghg_in= 0.d0
      self%pzrvs(i)%wghg_ex= 0.d0
    end do
    do j = 1, self%njunc
      ! effects on the 'from' component
      select case(self%f_typ(j))
        case('bndc')    ! from 'bndc'
          continue
        case('pzrv-g')  ! from 'pzrv-g'
          i = self%f_idx(j)
          mdot = self%juncs(j)%mdot
          self%pzrvs(i)%wg_ex = self%pzrvs(i)%wg_ex + mdot
          hg = self%pzrvs(i)%h_g  ! [J]
          mg = self%pzrvs(i)%m_g  ! [kg]
          self%pzrvs(i)%wghg_ex=self%pzrvs(i)%wghg_ex + mdot*hg/mg
        case('pzrv-l')  ! from 'pzrv-l'
          i = self%f_idx(j)
          mdot = self%juncs(j)%mdot
          self%pzrvs(i)%wl_su = self%pzrvs(i)%wl_su - mdot
          hl = self%pzrvs(i)%h_l  ! [J]
          ml = self%pzrvs(i)%m_l  ! [kg]
          self%pzrvs(i)%wlhl_su=self%pzrvs(i)%wlhl_su - mdot*hl/ml
        case default
          stop 'under development - 2026/11/01-0010'
      end select
      ! effects on the 'to' component
      select case(self%t_typ(j))
        case('bndc')    ! to 'bndc'
          continue
        case('pzrv-g')  ! to 'pzrv-g'
          i = self%t_idx(j)
          k = self%f_idx(j)
          mdot = self%juncs(j)%mdot
          self%pzrvs(i)%wg_in = self%pzrvs(i)%wg_in + mdot
          hg = self%bndcs(k)%s_hf  ! bndc enthalpy [J/kg]
          self%pzrvs(i)%wghg_in=self%pzrvs(i)%wghg_in + mdot*hg
        case('pzrv-l')  ! to 'pzrv-l'
          i = self%t_idx(j)
          k = self%f_idx(j)
          mdot = self%juncs(j)%mdot
          self%pzrvs(i)%wl_su = self%pzrvs(i)%wl_su + mdot
          hl = self%bndcs(k)%s_hf  ! bndc enthalpy [J/kg]
          self%pzrvs(i)%wlhl_su=self%pzrvs(i)%wlhl_su + mdot*hl
        case default
          stop 'under development - 2026/11/01-0125'
      end select
    end do
    ! calculate residuals of equations related to pvessels
    do i = 1, self%npzrv
      call self%pzrvs(i)%calc_res_pvessel(res(9*i-8:9*i))
    end do
    ! calculate residuals of equations related to junctions
    ! now under active development...
  end subroutine calc_res_fludomn


end module fludomn_mod