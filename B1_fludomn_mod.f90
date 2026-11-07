! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines type :: fludomn (fluid domain).
! 
! ------------------------------------------------------------------
module fludomn_mod
  !======= Inclusions =========
  use flupzrv_mod
  use flubrch_mod
  use flujunc_mod
  use flubndc_mod

  !======= Declarations =========
  implicit none

  type :: fludomn
    !=== Member Variables ===
    ! number of fluid components
    integer(c_int) :: nbndc, nbrch, njunc, npzrv
    ! array storing boundary cell objects
    type(flubndc), allocatable :: bndcs(:)
    ! array storing branch cell objects
    type(flubrch), allocatable :: brchs(:)
    ! array storing junction objects
    type(flujunc), allocatable :: juncs(:)
    ! array storing pressurizer vessel object
    type(flupzrv), allocatable :: pzrvs(:)
    ! array storing junction connection information
    character(len=20), allocatable:: f_ids(:), t_ids(:)
    character(len=20), allocatable:: f_typ(:), t_typ(:)
    integer(c_int), allocatable   :: f_idx(:), t_idx(:)

  contains
    !=== Member Procedures ===
    procedure :: init_bndcs
    procedure :: init_brchs
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

  ! initialize fluid branch cells
  subroutine init_brchs(self, brch_arr)
    !======= Declarations =========
    implicit none
    !=== I/0 variables ===
    class(fludomn), intent(inout) :: self
    class(flubrch), intent(in)    :: brch_arr(:)
    !=== local variables ===

    !======= Internals ============
    ! allocate and assign branch cell array
    self%nbrch = size(brch_arr)
    allocate(self%brchs(self%nbrch))
    self%brchs = brch_arr
  end subroutine init_brchs

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
    class(flupzrv), intent(in)    :: pzrv_arr(:)
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
      ! if connected with fluid boundary cell
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
      ! if connected with branch
      do i = 1, self%nbrch

      end do

      ! under aktiv development ...

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
                case('junc_norm','junc_orfz')
                  continue
                case default
                  print*, 'the "p_in" bc must serve as "from" of junc_norm or junc_orfz'
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
                case('junc_norm','junc_orfz')
                  continue
                case default
                  print*, 'the "p_out" bc must serve as "to" of junc_norm or junc_orfz'
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
    real(c_double) :: p_1, p_2, r_f
    real(c_double) :: s_hf
    integer(c_int) :: i, j, k
    !======= Internals ============
    ! calculate the mass and energy transport between fluid components
    do i = 1, self%npzrv
      ! initialization of pvessels inlets/outlets
      self%pzrvs(i)%wl_in  = 0.d0
      self%pzrvs(i)%wl_ex  = 0.d0
      self%pzrvs(i)%wg_in  = 0.d0
      self%pzrvs(i)%wg_ex  = 0.d0
      self%pzrvs(i)%wlhl_in= 0.d0
      self%pzrvs(i)%wlhl_ex= 0.d0
      self%pzrvs(i)%wghg_in= 0.d0
      self%pzrvs(i)%wghg_ex= 0.d0
    end do
    k = 1
    do j = 1, self%njunc
      mdot = self%juncs(j)%mdot
      if(mdot>0.d0) then
        ! effects on the 'from' component
        i = self%f_idx(j)
        select case(self%f_typ(j))
          case('bndc')
            r_f  = self%bndcs(i)%r_f
            s_hf = self%bndcs(i)%s_hf
          case('brch')
            r_f = self%brchs(i)%r_f
            s_hf= self%brchs(i)%s_hf
            self%brchs(i)%wf_ex = self%brchs(i)%wf_ex + mdot
            self%brchs(i)%wfhf_ex=self%brchs(i)%wfhf_ex + mdot*s_hf
          case('pzrv-g')
            r_f = self%pzrvs(i)%r_g
            s_hf= self%pzrvs(i)%s_hg
            self%pzrvs(i)%wg_ex = self%pzrvs(i)%wg_ex + mdot
            self%pzrvs(i)%wghg_ex=self%pzrvs(i)%wghg_ex + mdot*s_hf 
          case('pzrv-l')
            r_f = self%pzrvs(i)%r_l
            s_hf= self%pzrvs(i)%s_hl
            self%pzrvs(i)%wl_ex = self%pzrvs(i)%wl_ex + mdot
            self%pzrvs(i)%wlhl_ex=self%pzrvs(i)%wlhl_ex + mdot*s_hf
          case default
            stop 'under development - 2026/11/07-0947'
        end select
        ! effects on the 'to' component
        i = self%t_idx(j)
        select case(self%t_typ(j))
          case('bndc')
            continue
          case('brch')
            self%brchs(i)%wf_in = self%brchs(i)%wf_in + mdot
            self%brchs(i)%wfhf_in=self%brchs(i)%wfhf_in + mdot*s_hf
          case('pzrv-g')
            self%pzrvs(i)%wg_in = self%pzrvs(i)%wg_in + mdot
            self%pzrvs(i)%wghg_in=self%pzrvs(i)%wghg_in + mdot*s_hf
          case('pzrv-l')
            self%pzrvs(i)%wl_in = self%pzrvs(i)%wl_in + mdot
            self%pzrvs(i)%wlhl_in=self%pzrvs(i)%wlhl_in + mdot*s_hf
          case default
            stop 'under development - 2026/11/07-1205'
        end select
      else
        ! here mdot<=0.d0
        ! effects on the 'to' component
        i = self%t_idx(j)
        select case(self%t_typ(j))
          case('bndc')
            r_f = self%bndcs(i)%r_f
            s_hf= self%bndcs(i)%s_hf
          case('brch')
            r_f = self%brchs(i)%r_f
            s_hf= self%brchs(i)%s_hf
            self%brchs(i)%wf_in = self%brchs(i)%wf_in + mdot
            self%brchs(i)%wfhf_in=self%brchs(i)%wfhf_in + mdot*s_hf
          case('pzrv-g')
            r_f = self%pzrvs(i)%r_g
            s_hf= self%pzrvs(i)%s_hg
            self%pzrvs(i)%wg_in = self%pzrvs(i)%wg_in + mdot
            self%pzrvs(i)%wghg_in=self%pzrvs(i)%wghg_in + mdot*s_hf
          case('pzrv-l')
            r_f = self%pzrvs(i)%r_l
            s_hf= self%pzrvs(i)%s_hl
            self%pzrvs(i)%wl_in = self%pzrvs(i)%wl_in - mdot
            self%pzrvs(i)%wlhl_in=self%pzrvs(i)%wlhl_in + mdot*s_hf
          case default
            stop 'under development - 2026/11/07-1249'
        end select
        ! effects on the 'from' component
        i = self%f_idx(j)
        select case(self%f_typ(j))
          case('bndc')
            continue
          case('brch')
            self%brchs(i)%wf_ex = self%brchs(i)%wf_ex + mdot
            self%brchs(i)%wfhf_ex=self%brchs(i)%wfhf_ex + mdot*s_hf
          case('pzrv-g')
            self%pzrvs(i)%wg_ex = self%pzrvs(i)%wg_ex + mdot
            self%pzrvs(i)%wghg_ex=self%pzrvs(i)%wghg_ex + mdot*s_hf
          case('pzrv-l')
            self%pzrvs(i)%wl_ex = self%pzrvs(i)%wl_ex + mdot
            self%pzrvs(i)%wlhl_ex=self%pzrvs(i)%wlhl_ex + mdot*s_hf
          case default
            stop 'under development - 2026/11/07-1620'
        end select
      end if
      ! calculate residuals of equations related to junctions
      select case(self%juncs(j)%typ)
        case('junc_flow')
          continue
        case('junc_pump')
          stop 'under development - 2026/11/07-1632'
        case('junc_norm')
          stop 'under development - 2026/11/07-1632'
        case('junc_orfz')
          ! get the pressure of 'from' component
          i = self%f_idx(j)
          select case(self%f_typ(j))
            case('bndc')
              p_1 = self%bndcs(i)%p_f
            case('brch')
              p_1 = self%brchs(i)%p_f
            case('pzrv-g', 'pzrv-l')
              p_1 = self%pzrvs(i)%p_g
            case default
              stop 'under development - 2026/11/07-1653'
          end select
          ! get the pressure of 'to' component
          i = self%t_idx(j)
          select case(self%t_typ(j))
            case('bndc')
              p_2 = self%bndcs(i)%p_f
            case('brch')
              p_2 = self%brchs(i)%p_f
            case('pzrv-g', 'pzrv-l')
              p_2 = self%pzrvs(i)%p_g
            case default
              stop 'under development - 2026/11/07-1707'
          end select
          ! calculate the residuals





          k = k + 1
        case default
          stop 'under development - 2026/11/07-1628'
      end select

    end do
    !
    ! calculate residuals of equations related to pvessels
    do i = 1, self%npzrv
      call self%pzrvs(i)%calc_res_flupzrv(res(neq_mdot+9*i-8:neq_mdot+9*i))
    end do
    ! caluclate residuals of equations related to branches
    do i = 1, self%nbrch
      ! ...
    end do

  end subroutine calc_res_fludomn


end module fludomn_mod