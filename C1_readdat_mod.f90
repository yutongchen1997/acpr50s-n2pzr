! ------------------------------------------------------------------
! Programmer(s): Yutong. Chen @ CGN::CNPRI::SRD::HX
! ------------------------------------------------------------------
! 
! This module defines subroutines reading data of material 
! properties from '.dat' files.
! 
! ------------------------------------------------------------------
module readdat_mod
  !======= Inclusions ===========
  use globdat_mod
  
contains

! ----------------------------------------------------------------------
! read_n2g_dat: a subroutine reading thermal property datapoints for n2
! gas, filling global tpf_n2g(6,:,:) array containing list of rho, niu, 
! miu, kpa, cp, t
! ----------------------------------------------------------------------
  subroutine read_n2g_dat()
    !======= Declarations =========
    implicit none
    !== input variables ==
    !== output variables ==
    !== local variables ==
    real(c_double) :: c_real
    real(c_double) :: p1, p2
    integer(c_int) :: i, j, filelen
    integer(c_int) :: lenp, lenh     ! p, h array length
    !======= Internals ============
    print*, 'Reading N2 property data from tpf_n2g.dat file...'
    ! open n2 gas data file
    open(unit=n2gfu, status='old', file='tpf_n2g.dat', form='formatted')
    filelen = 0
    ! loop over the .dat file to get array length
    do while(.TRUE.)
      read(n2gfu, '(F12.5)', end=100) c_real
      filelen = filelen + 1
    end do
    100 continue
    ! get the length of h array
    lenh = 0
    rewind(n2gfu)
    do while(.TRUE.)
      read(n2gfu,'(F12.5)') p1
      if(lenh == 0) then
        p2 = p1
        lenh = lenh + 1
      else
        if(p1 == p2) then
          lenh = lenh + 1
        else
          exit
        end if
      end if
    end do
    ! the length of p array
    lenp = filelen/lenh
    ! allocate data arrays
    allocate(p_n2g(lenp))
    allocate(h_n2g(lenh))
    allocate(tpf_n2g(6,lenp,lenh))
    ! read data into global tpf_n2g(6,lenp,lenh) array
    ! tpf_n2g(6,lenp,lenh):
    ! column 1-density, 2-kinetic_viscosity, 3-dynamic_viscosity, 4-heat_conductivity,
    ! 5-heat_capacity, 6-temperature
    rewind(n2gfu)
    do i = 1, lenp
      do j = 1, lenh
        read(n2gfu,'(9F12.5)')p_n2g(i),h_n2g(j),tpf_n2g(1,i,j),&
        tpf_n2g(2,i,j), tpf_n2g(3,i,j), tpf_n2g(4,i,j), tpf_n2g(5,i,j),  &
        tpf_n2g(6,i,j)
      end do
    end do
    close(n2gfu)
    print*, 'Done.'
  end subroutine read_n2g_dat

! ----------------------------------------------------------------------
! read_h2o_dat: a subroutine reading thermal property datapoints for h2o 
! fluid, filling global tpf_h2o(7,:,:) array containing list of rho,
! niu, miu, kpa, cp, t, xe and global sat_h2o(11,:) 
! array containing list of saturated rhol, rhog, miul, miug, cpl, cpg, 
! kpal, kpag, tsat, hls, hgs, hgl, sigma
! ----------------------------------------------------------------------
  subroutine read_h2o_dat()
    !======= Declarations =========
    implicit none
    !== input variables ==
    !== output variables ==
    !== local variables ==
    real(c_double) :: c_real
    real(c_double) :: p1, p2
    integer(c_int) :: i, j, filelen
    integer(c_int) :: lenp, lenh     ! p, h array length
    !======= Internals ============
    print*, 'Reading H2O property data from tpf_h2o.dat & sat_h2o.dat file...'
    ! open h2o fluid data file
    open(unit=h2ofu, status='old', file='tpf_h2o.dat', form='formatted')
    open(unit=h2osu, status='old', file='sat_h2o.dat', form='formatted')
    filelen = 0
    ! loop over the .dat file to get array length
    do while(.TRUE.)
      read(h2ofu, '(F12.5)', end=200) c_real
      filelen = filelen + 1
    end do
    200 continue
    ! get the length of h array
    lenh = 0 
    rewind(h2ofu)
    do while(.TRUE.)
      read(h2ofu,'(F12.5)') p1
      if(lenh == 0) then
        p2 = p1
        lenh = lenh + 1
      else
        if(p1 == p2) then
          lenh = lenh + 1
        else
          exit
        end if
      end if
    end do
    ! the length of p array
    lenp = filelen/lenh
    ! allocate data arrays
    allocate(p_h2o(lenp))
    allocate(h_h2o(lenh))
    allocate(tpf_h2o( 7,lenp,lenh))
    allocate(sat_h2o(13,lenp))
    ! read data into global tpf_h2o(7,lenp,lenh) and sat_h2o(11,lenp) array
    ! tpf_h2o(7,lenp,lenh):
    ! column 1-density, 2-kinetic_viscosity, 3-dynamic_viscosity, 4-heat_conductivity,
    ! 5-heat_capacity, 6-temperature, 7-gas_mass_quality
    ! sat_h2o(11,lenp):
    ! column 1-rhol, 2-rhog, 3-miul, 4-miug, 5-cpl, 6-cpg, 7-kpal, 8-kpag, 9-tsat, 
    ! 10-hls, 11-hgs, 12-hgl, 13-sigma
    rewind(h2ofu)
    rewind(h2osu)
    do i = 1, lenp
      read(h2osu,'(13F12.5)') sat_h2o(1,i), sat_h2o(2,i), sat_h2o(3,i),     &
      sat_h2o(4,i), sat_h2o(5,i), sat_h2o(6,i), sat_h2o(7,i), sat_h2o(8,i), &
      sat_h2o(9,i), sat_h2o(10,i),sat_h2o(11,i),sat_h2o(12,i),sat_h2o(13,i)
      do j = 1, lenh
        read(h2ofu,'(9F12.5)')p_h2o(i),h_h2o(j),tpf_h2o(1,i,j),&
        tpf_h2o(2,i,j), tpf_h2o(3,i,j), tpf_h2o(4,i,j), tpf_h2o(5,i,j),  &
        tpf_h2o(6,i,j), tpf_h2o(7,i,j)
      end do
    end do
    close(h2ofu)
    close(h2osu)
    print*, 'Done.'
  end subroutine read_h2o_dat

end module readdat_mod