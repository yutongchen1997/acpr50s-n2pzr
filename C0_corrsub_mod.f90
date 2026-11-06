! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module contains all subroutines that calculate
! auxiliary correlations
! 
! ------------------------------------------------------------------
module corrsub_mod
  !======= Inclusions =========
  use interpo_mod
  use readdat_mod
  !======= Declarations =========
  implicit none

contains
! ----------------------------------------------------------------------
! matpro_s: a subroutine for solid material properties calculation
! input  - solid material id and temperature [K]
! output - solid material properties
! ----------------------------------------------------------------------
  subroutine matpro_s(matid, temp, rho, cp, kpa)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! material id
    character(len=20) :: matid
    ! temperature
    real(c_double), intent(in) :: temp(:)
    !== output variables ==
    ! density, heat_capacity, heat_conductivity
    real(c_double), intent(out) :: rho(:), cp(:), kpa(:)
    !== local variables ==
    real(c_double), allocatable :: tc(:), tf(:)
    integer(c_int) :: i, nlen
    ! coefficients for correlations
    real(c_double) :: a0, a1, a2, a3, a4

    !======= Internals ============
    ! check if input is correct
    nlen = size(temp)
    if(size(rho)/=nlen .or. size(cp)/=nlen .or. size(kpa)/=nlen) then
      stop 'Array size does not match in matpro_s'
    end if
    ! alllcate temporary arrays
    allocate(tc(nlen))
    allocate(tf(nlen))
    tc = temp - 273.15d0   ! Celsius temperature array
    tf = tc*1.8d0 + 32.d0  ! Fahrenheit temperature array
    select case (matid)
      case('ss316') ! Stainless Steel, Type 316, developed based on TRACE manual
        ! density
        rho = 7954.d0
        ! heat capacity
        do i = 1, nlen
          ! TRACE model limitations
          tf(i) = MIN(tf(i), 2400.d0)
        end do
        a0 = 426.17d0
        a1 = 4.3816d-1
        a2 =-6.3759d-4
        a3 = 4.4803d-7
        a4 =-1.0729d-10
        cp  = a0 + a1*tf + a2*tf**(2.d0) + a3*tf**(3.d0) + a4*tf**(4.d0)
        ! heat conductivity
        kpa = 9.248d0 + 1.571d-2*temp
      case default
        print*, 'no match for material name = ', matid
        stop 'program stopped at matpro_s'
    end select  
    ! deallocate arrays
    deallocate(tc)
    deallocate(tf)
  end subroutine matpro_s

! ----------------------------------------------------------------------
! pgeos_n2g: a subroutine to calculate partial pressure of n2 gas with 
! equation of state model 
! input  - eos type, gas density [kg/m^3], gas enthalpy [J/kg]
! output - gas pressure [Pa]
! ----------------------------------------------------------------------
  subroutine pgeos_n2g(eosmod, rg, hg, pg)
    !======= Declarations =========
    implicit none
    !== input variables ==
    character(len=20), intent(in):: eosmod  ! eos model type 
    ! gas density, gas enthalpy
    real(c_double), intent(in) :: rg, hg
    !== output variables ==
    ! gas partial pressure
    real(c_double), intent(out) :: pg
    !== local variables ==
    ! variables for Newton iterative algorithm
    real(c_double) :: d_rg, rg_1, rg_2, pg_1, pg_2
    real(c_double) :: k_rp    ! the slope
    integer(c_int) :: icount

    !======= Internals ============
    select case (eosmod)
      case('igas')    ! use ideal gas eos model
        pg = rg*hg*(gam_n2-1.d0)/gam_n2
      case('nist')    ! interpolate with NIST table
        ! initial guess of pg by ideal-gas model
        pg_1 = rg*hg*(gam_n2-1.d0)/gam_n2 
        ! then calcualte n2g density with NIST table
        call interp2d(p_n2g,h_n2g,tpf_n2g(1,:,:),pg_1,hg,rg_1)    !column 1 - density
        d_rg = rg-rg_1    ! difference between input and estimated density
        icount = 0
        do while(abs(d_rg)/rg>rtol)
          if(d_rg>0.d0) then
            ! the estimated pg is lower than true value
            pg_2 = pg_1+1.d0
            call interp2d(p_n2g,h_n2g,tpf_n2g(1,:,:),pg_2,hg,rg_2)
            k_rp = rg_2-rg_1           ! the slope of pg vs. rg curve
          else
            ! the estimated pg is higher than true value
            pg_2 = pg_1-1.d0
            call interp2d(p_n2g,h_n2g,tpf_n2g(1,:,:),pg_2,hg,rg_2)
            k_rp = rg_1-rg_2           ! the slope of pg vs. rg curve
          end if 
          pg_1 = pg_1 + d_rg/k_rp    ! correction of pg_1
          call interp2d(p_n2g,h_n2g,tpf_n2g(1,:,:),pg_1,hg,rg_1)
          d_rg = rg-rg_1
          ! max iterations limiter
          icount = icount+1
          if(icount>50) then
            print*, 'pgeos_n2g was not converged under nist mode!'
            exit
          end if 
        end do
        pg = pg_1
      case default
        print*, 'no match for eos model = ', eosmod
        stop 'program stopped at pgeos_n2g'
    end select
  end subroutine pgeos_n2g

! ----------------------------------------------------------------------
! tgeos_n2g: a subroutine to calculate n2 gas temperature with 
! equation of state model 
! input  - eos type, gas pressure [Pa], gas enthalpy [J/kg]
! output - gas temperature [K]
! ----------------------------------------------------------------------
  subroutine tgeos_n2g(eosmod, pg, hg, tg)
    !======= Declarations =========
    implicit none
    !== input variables ==
    character(len=20), intent(in):: eosmod  ! eos model type 
    ! gas pressure, gas enthalpy
    real(c_double), intent(in) :: pg, hg
    !== output variables ==
    ! gas temperature
    real(c_double), intent(out) :: tg
 
    !======= Internals ============
    select case (eosmod)
      case('igas')    ! use ideal gas eos model
        tg = hg*(gam_n2-1.d0)*mol_n2/gam_n2/r_u
      case('nist')    ! interpolate with NIST table
        call interp2d(p_n2g,h_n2g,tpf_n2g(6,:,:),pg,hg,tg)    !column 6 - temperature
      case default
        print*, 'no match for eos model = ', eosmod
        stop 'program stopped at tgeos_n2g'
    end select
  end subroutine tgeos_n2g

! ----------------------------------------------------------------------
! tleos_h2o: a subroutine to calculate temperature of h2o with 
! equation of state model 
! input  - gas pressure [Pa], water enthalpy [J/kg]
! output - water temperature [K]
! ----------------------------------------------------------------------
  subroutine tleos_h2o(pg, hl, tl)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! gas pressure, water enthalpy
    real(c_double), intent(in) :: pg, hl
    !== output variables ==
    ! water temperature 
    real(c_double), intent(out) :: tl
    !======= Internals ============
    call interp2d(p_h2o,h_h2o,tpf_h2o(6,:,:),pg,hl,tl)   !column 6 - temperature
  end subroutine tleos_h2o

! ----------------------------------------------------------------------
! rgeos_n2g: a subroutine to calculate density of n2g with 
! equation of state model 
! input  - eos type, gas pressure [Pa], gas enthalpy [J/kg]
! output - gas density [kg/m^3]
! ----------------------------------------------------------------------
  subroutine rgeos_n2g(eosmod, pg, hg, rg)
    !======= Declarations =========
    implicit none
    !== input variables ==
    character(len=20), intent(in):: eosmod  ! eos model type 
    ! gas pressure, gas enthalpy
    real(c_double), intent(in) :: pg, hg
    !== output variables ==
    ! gas density
    real(c_double), intent(out) :: rg
    !======= Internals ============
    select case (eosmod)
      case('igas')    ! use ideal gas eos model
        rg = gam_n2*pg/hg/(gam_n2-1.d0)
      case('nist')    ! interpolate with NIST table
        call interp2d(p_n2g,h_n2g,tpf_n2g(1,:,:),pg,hg,rg)
      case default
        print*, 'no match for eos model = ', eosmod
        stop 'program stopped at rgeos_n2g'
    end select
  end subroutine rgeos_n2g

! ----------------------------------------------------------------------
! rfeos_h2o: a subroutine to calculate density of h2o with 
! equation of state model 
! input  - water pressure [Pa], water enthalpy [J/kg]
! output - water density [kg/m^3]
! ----------------------------------------------------------------------
  subroutine rfeos_h2o(pf, hf, rf)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! gas pressure, water enthalpy
    real(c_double), intent(in) :: pf, hf
    !== output variables ==
    ! water density 
    real(c_double), intent(out) :: rf
    !======= Internals ============
    call interp2d(p_h2o,h_h2o,tpf_h2o(1,:,:),pf,hf,rf)   !column 1 - density
  end subroutine rfeos_h2o

! ----------------------------------------------------------------------
! hgeos_n2g: a subroutine to calculate enthalpy of n2 gas based on 
! equation of state model 
! input  - gas pressure [Pa], temperature [K]
! output - gas enthalpy [J/kg]
! ----------------------------------------------------------------------
  subroutine hgeos_n2g(eosmod, pg, tg, hg)
    !======= Declarations =========
    implicit none
    !== input variables ==
    character(len=20), intent(in):: eosmod  ! eos model type 
    ! gas pressure, gas temperature
    real(c_double), intent(in) :: pg, tg
    !== output variables ==
    ! gas enthalpy
    real(c_double), intent(out) :: hg
    !== local variables ==
    ! variables for Newton iterative algorithm
    real(c_double) :: d_tg, tg_1, tg_2, hg_1, hg_2
    real(c_double) :: k_th    ! the slope
    integer(c_int) :: icount
    !======= Internals ============
    select case (eosmod)
      case('igas')    ! use ideal gas eos model
        hg = gam_n2*r_u*tg/mol_n2/(gam_n2-1.d0)
      case('nist')    ! interpolate with NIST table
        ! initial guess of hg by ideal-gas model
        hg_1 = gam_n2*r_u*tg/mol_n2/(gam_n2-1.d0)
        ! then calcualte n2g temperature with NIST table
        call interp2d(p_n2g,h_n2g,tpf_n2g(6,:,:),pg,hg_1,tg_1)
        d_tg = tg-tg_1    ! difference between input and estimated tg
        icount = 0
        do while(abs(d_tg)/tg>1.0d-6)
          if(d_tg>0.d0) then
            ! the estimated hg is lower than true value
            hg_2 = hg_1+1.d0
            call interp2d(p_n2g,h_n2g,tpf_n2g(6,:,:),pg,hg_2,tg_2)
            k_th = tg_2-tg_1           ! the slope of tg vs. hg curve
          else
            ! the estimated hg is higher than true value
            hg_2 = hg_1-1.d0
            call interp2d(p_n2g,h_n2g,tpf_n2g(6,:,:),pg,hg_2,tg_2)
            k_th = tg_1-tg_2           ! the slope of tg vs. hg curve
          end if
          hg_1 = hg_1 + d_tg/k_th    ! correction of hg_1
          call interp2d(p_n2g,h_n2g,tpf_n2g(6,:,:),pg,hg_1,tg_1)
          d_tg = tg-tg_1
          ! max iterations limiter
          icount = icount+1
          if(icount>50) then
            print*, 'hgeos_n2g was not converged under nist mode!'
            exit
          end if 
        end do
        hg = hg_1
      case default
        print*, 'no match for eos model = ', eosmod
        stop 'program stopped at hgeos_n2g'
    end select
  end subroutine hgeos_n2g

! ----------------------------------------------------------------------
! hfeos_h2o: a subroutine to calculate enthalpy of h2o liquid based on 
! equation of state model 
! input  - h2o pressure [Pa], h2o temperature [K]
! output - h2o enthalpy [J/kg]
! ----------------------------------------------------------------------
  subroutine hfeos_h2o(pf, tf, hf)
    !======= Declarations =========
    implicit none
    !== input variables ==
    ! h2o pressure, h2o temperature
    real(c_double), intent(in) :: pf, tf
    !== output variables ==
    ! h2o enthalpy
    real(c_double), intent(out) :: hf
    !== local variables ==
    ! variables for Newton iterative algorithm
    real(c_double) :: d_tf, tf_1, tf_2, hf_1, hf_2
    real(c_double) :: k_th    ! the slope
    integer(c_int) :: icount
    !======= Internals ============
    ! initial guess of hl assuming saturation state - hlsat(tl)
    call interp1d(sat_h2o(9,:),sat_h2o(10,:),tf,hf_1)
    ! then calculate h2o temperature with NIST table
    call interp2d(p_h2o,h_h2o,tpf_h2o(6,:,:),pf,hf_1,tf_1)
    d_tf = tf - tf_1    ! difference between input and estimated tl
    icount = 0
    do while(abs(d_tf)/tf>rtol)
      if(d_tf>0.d0) then
        ! the estimated hl is lower than true value
        hf_2 = hf_1+1.d0
        call interp2d(p_h2o,h_h2o,tpf_h2o(6,:,:),pf,hf_2,tf_2)
        k_th = tf_2-tf_1           ! the slope
      else
        ! the estimated hf is higher than true value
        hf_2 = hf_1-1.d0
        call interp2d(p_h2o,h_h2o,tpf_h2o(6,:,:),pf,hf_2,tf_2)
        k_th = tf_1-tf_2           ! the slope
      end if
      hf_1 = hf_1 + d_tf/k_th    ! correction of hl_1
      call interp2d(p_h2o,h_h2o,tpf_h2o(6,:,:),pf,hf_1,tf_1)
      d_tf = tf - tf_1
      ! max interations limiter
      icount = icount+1
      if(icount>50) then
        print*, 'hfeos_h2o was not converged!'
        exit
      end if 
    end do 
    hf = hf_1
  end subroutine hfeos_h2o


  
end module corrsub_mod