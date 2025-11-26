!! Taken from: https://github.com/moritzgubler/vc-sqnm

!! The variable cell shape optimization method is based on the following 
!! paper: https://arxiv.org/abs/2206.07339
!! More details about the SQNM optimization method are available here:
!! https://comphys.unibas.ch/publications/Schaefer2015.pdf
!! Author of this document: Moritz Gubler 

! Copyright (C) 2022 Moritz Gubler
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module historylist
  use iso_c_binding
  implicit none
  type hist_list
  !! Historylist that is used by the sqnm class.
  !! More informations about the SQNM algorithm can be found here: https://aip.scitation.org/doi/10.1063/1.4905665
  integer(c_int) :: nhistx
  integer(c_int) :: ndim
  integer(c_int), private :: icount
  logical, private :: is_initialized = .false.
  real(c_double), allocatable, dimension(:, :) :: list
  real(c_double), allocatable, dimension(:, :) :: diff_list
  real(c_double), allocatable, dimension(:, :) :: norm_diff_list
  real(c_double), allocatable, dimension(:) :: old_x
  contains
  procedure :: init
  procedure :: add
  procedure :: get_length
  procedure :: close_history_list

  end type hist_list
contains

  subroutine init(t, ndim, nhistx)
    !! initializes the historylist object.
    class(hist_list) :: t
    integer(c_int), intent(in) :: ndim
    !! dimension of the optimization problem
    integer(c_int), intent(in) :: nhistx
    !! maximal length of history list

    t%icount = 1
    t%ndim = ndim
    t%nhistx = nhistx
    allocate(t%list(ndim, nhistx), t%diff_list(ndim, nhistx))
    allocate(t%old_x(ndim), t%norm_diff_list(ndim, nhistx))
    t%is_initialized = .true.
  end subroutine init

  subroutine add(t, x)
    !! Add a vector to the history list.
    class(hist_list) :: t
    real(c_double), intent(in) :: x(t%ndim)
    !! Vector to add.
    integer :: i
   
    if ( t%icount <= t%nhistx) then ! list not yet full
      t%list(:, t%icount) = x
      do i = 2, t%icount
        t%diff_list(:, i - 1) = t%list(:, i) - t%list(:, i - 1)
        t%norm_diff_list(:, i-1) = t%diff_list(:, i-1) / norm2(t%diff_list(:, i-1))
      end do
      t%icount = t%icount + 1
    else ! list is full
      t%icount = t%nhistx + 2
      t%old_x = t%list(:, 1)
      do i = 1, t%nhistx - 1
        t%list(:, i) = t%list(:, i + 1)                
      end do
      t%list(:, t%nhistx) = x
      t%diff_list(:, 1) =  t%list(:, 1) - t%old_x
      t%norm_diff_list(:, 1) = t%diff_list(:, 1) / norm2(t%diff_list(:, 1))
      do i = 2, t%nhistx
        t%diff_list(:, i) = t%list(:, i) - t%list(:, i - 1)
        t%norm_diff_list(:, i) = t%diff_list(:, i) / norm2(t%diff_list(:, i))
      end do
    end if
  end subroutine add

  integer function get_length(t)
  !! returns the length of the historylist.
    class(hist_list) :: t
    get_length = t%icount - 2
  end function get_length

  subroutine close_history_list(t)
    !! closes the history list and deallocates the memory.
    class(hist_list) :: t
    if (t%is_initialized) then
      deallocate(t%list, t%diff_list, t%norm_diff_list, t%old_x)
      t%is_initialized = .false.
    end if
  end subroutine close_history_list

end module historylist

module sqnm
  use iso_c_binding
  use historylist
  implicit none

  type sqnm_optimizer
    !! This class is an implementation of the stabilized quasi newton optimization method.
    !! More informations about the algorithm can be found here: https://aip.scitation.org/doi/10.1063/1.4905665
    integer(c_int) :: nhistx
    integer(c_int) :: ndim
    real(c_double) :: eps_subsp
    real(c_double) :: alpha0
    type(hist_list) :: x_list
    type(hist_list) :: flist
    real(c_double) :: alpha
    real(c_double), allocatable, dimension(:) :: dir_of_descent
    real(c_double) :: prev_f
    !! previous value of target function
    real(c_double), allocatable, dimension(:) :: prev_df_dx
    !! previous derivative of target function
    real(c_double), allocatable, dimension(:, :) :: s_evec
    real(c_double), allocatable, dimension(:) :: s_eval
    real(c_double), allocatable, dimension(:, :) :: dr_subsp
    real(c_double), allocatable, dimension(:, :) :: df_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec_subsp
    real(c_double), allocatable, dimension(:, :) :: h_evec
    real(c_double), allocatable, dimension(:) :: h_eval
    real(c_double), allocatable, dimension(:) :: res
    real(c_double), allocatable, dimension(:) :: res_temp
    real(c_double) :: gainratio
    integer :: nhist
    real(c_double), allocatable, dimension(:) :: expected_positions
    logical :: estimate_step_size

    INTEGER :: lwork
    REAL(8), DIMENSION(:), ALLOCATABLE:: work

    contains
    procedure :: initialize_sqnm
    procedure :: sqnm_step
    procedure :: get_lower_bound
    procedure :: close_sqnm
  end type sqnm_optimizer
contains

subroutine initialize_sqnm(t, ndim, nhistx, alpha, alpha0, eps_subsp)
  class(sqnm_optimizer) :: t
  integer(c_int) :: ndim
  !! dimension of the optimization problem
  integer(c_int) :: nhistx
  !! maximal length of history list
  real(c_double) :: alpha
  !! Initial step size. Should be approximately the inverse of the largest eigenvalue of the Hessian matrix.
  real(c_double) :: alpha0
  !! Lowest step size that is allowed.
  real(c_double) :: eps_subsp
  !! Lower limit on linear dependencies in history list.
  
  t%ndim = ndim
  t%nhistx = nhistx
  t%estimate_step_size = .false.
  if ( alpha <= 0.d0 ) then
    t%estimate_step_size = .true.
    t%alpha = -alpha
  else
    t%alpha = alpha
  end if
  
  t%alpha0 = alpha0
  t%eps_subsp = eps_subsp
  call t%x_list%init(ndim, nhistx)
  call t%flist%init(ndim, nhistx)

  allocate(t%s_evec(t%nhistx, t%nhistx), t%s_eval(t%nhistx))
  allocate(t%prev_df_dx(ndim))
  allocate(t%dr_subsp(t%ndim, nhistx))
  allocate(t%df_subsp(t%ndim, nhistx))
  allocate(t%h_evec_subsp(t%nhistx, t%nhistx))
  allocate(t%h_eval(nhistx))
  allocate(t%h_evec(t%ndim, t%nhistx))
  allocate(t%res(nhistx))
  allocate(t%res_temp(ndim))
  allocate(t%dir_of_descent(ndim))
  allocate(t%expected_positions(ndim))

  t%lwork = 100 * t%nhistx
  allocate(t%work(t%lwork))
  
end subroutine initialize_sqnm

subroutine close_sqnm(t)
    !! closes the sqnm optimizer and frees all memory that was allocated.
    !! This function should be called at the end of the optimization.
    class(sqnm_optimizer) :: t
    
    call t%x_list%close_history_list()
    call t%flist%close_history_list()

    if (allocated(t%s_evec)) deallocate(t%s_evec)
    if (allocated(t%s_eval)) deallocate(t%s_eval)
    if (allocated(t%prev_df_dx)) deallocate(t%prev_df_dx)
    if (allocated(t%dr_subsp)) deallocate(t%dr_subsp)
    if (allocated(t%df_subsp)) deallocate(t%df_subsp)
    if (allocated(t%h_evec_subsp)) deallocate(t%h_evec_subsp)
    if (allocated(t%h_eval)) deallocate(t%h_eval)
    if (allocated(t%h_evec)) deallocate(t%h_evec)
    if (allocated(t%res)) deallocate(t%res)
    if (allocated(t%res_temp)) deallocate(t%res_temp)
    if (allocated(t%dir_of_descent)) deallocate(t%dir_of_descent)
    if (allocated(t%expected_positions)) deallocate(t%expected_positions)
    if (allocated(t%work)) deallocate(t%work)

end subroutine

subroutine sqnm_step(t, x, f_of_x, df_dx, dir_of_descent)
  !! Calculates a set of new coordinates based on the function value and derivatives provide on input.
  !! The idea behind this function is, that the user evaluates the function at the new point this method suggested and
  !! then calls this method again with the function value at the new point until convergence was reached.

  class(sqnm_optimizer) :: t
  real(c_double), intent(in) :: x(t%ndim)
  !! Array containg position vector x.
  real(c_double), intent(in) :: f_of_x
  !! Value of target function at point x.
  real(c_double), intent(in) :: df_dx(t%ndim)
  !! derivative of function f with respect to x.
  real(c_double), intent(out) :: dir_of_descent(t%ndim)
  !! Direction of descent x+dir_of_descent is the new point
  !! of the function that should be evaluated by the user.

  real(c_double) :: l1, l2


  integer :: dim_subsp
  integer :: i, ihist, k, j

  ! lapack variables
  INTEGER :: info

  !! check if gradient is already sufficiently small and return if this is the case.
  if ( maxval(abs( df_dx )) < 1.d-11 ) then
    dir_of_descent = 0.d0
    t%dir_of_descent = 0.d0
    return
  end if

  call t%x_list%add(x)
  call t%flist%add(df_dx)
  t%nhist = t%x_list%get_length()

  if ( t%nhist == 0 ) then !! first step
    t%dir_of_descent = - t%alpha * df_dx
  else

    ! check if positions have been changed and print a warning if they were.
    if ( maxval(abs(x - t%expected_positions)) > 1.d-9 ) then
      print*, "SQNM was not called with positions that were expected. If this was not done on purpose, it is probably a bug."
      print*, "Were atoms that left the simulation box put back into the cell? This is not allowed."
    end if

    if ( t%estimate_step_size ) then
      l1 = (f_of_x - t%prev_f + t%alpha * norm2(t%prev_df_dx)**2) / (.5d0 * (t%alpha**2) * (norm2(t%prev_df_dx)**2))
      l2 = norm2(df_dx - t%prev_df_dx) / (t%alpha * norm2(t%prev_df_dx))
      t%alpha = 1 / max(l1, l2)
      print'(a, g0.4)', 'Automatic initial step size guess: ', t%alpha
      t%estimate_step_size = .false.
    else
      ! calculate gainratio
      t%gainratio = (f_of_x - t%prev_f) / (.5d0 * dot_product(t%dir_of_descent, t%prev_df_dx))
      if (t%gainratio < 0.5d0 ) t%alpha = max(t%alpha0, t%alpha * 0.65d0)
      if (t%gainratio > 1.05d0) t%alpha = t%alpha * 1.05d0
    end if

    ! calculate overlab matrix of basis
    t%s_evec(:t%nhist, :t%nhist) =  matmul(transpose(t%x_list%norm_diff_list(:, :t%nhist)) &
    , t%x_list%norm_diff_list(:, :t%nhist))
    call dsyev('v', 'u', t%nhist, t%s_evec(:t%nhist, :t%nhist), t%nhist, t%s_eval(:t%nhist) &
      , t%work, t%lwork, info)
    if (info /= 0) stop 's dsyev'

    dim_subsp = 0
    do i = 1, t%nhist
      if (t%s_eval(i) / t%s_eval(t%nhist) > t%eps_subsp) then
        dim_subsp = dim_subsp + 1
      !else
      !  print*, 'remove dimension'
      end if
    end do
    t%s_eval(1:dim_subsp) = t%s_eval((t%nhist - dim_subsp + 1):t%nhist)
    t%s_evec(:, 1:dim_subsp) = t%s_evec(:, (t%nhist - dim_subsp + 1):t%nhist)

    ! compute eq. 11
    t%dr_subsp(:,:dim_subsp) = 0.d0
    t%df_subsp(:,:dim_subsp) = 0.d0
    do i = 1, dim_subsp
      do ihist = 1, t%nhist
        t%dr_subsp(:, i) = t%dr_subsp(:, i) + t%s_evec(ihist, i) * t%x_list%norm_diff_list(:, ihist)
        t%df_subsp(:, i) = t%df_subsp(:, i) + t%s_evec(ihist, i) * t%flist%diff_list(:, ihist) &
          / norm2(t%x_list%diff_list(:, ihist)) 
      end do
      t%dr_subsp(:, i) = t%dr_subsp(:, i) / sqrt(t%s_eval(i))
      t%df_subsp(:, i) = t%df_subsp(:, i) / sqrt(t%s_eval(i))
    end do

    !! compute eq. 13
    t%h_evec_subsp(:dim_subsp, :dim_subsp) = .5d0 * (matmul(transpose(t%df_subsp(:,:dim_subsp)), t%dr_subsp(:,:dim_subsp)) &
        + matmul(transpose(t%dr_subsp(:,:dim_subsp)), t%df_subsp(:,:dim_subsp)))
    call dsyev('v', 'l', dim_subsp, t%h_evec_subsp(:dim_subsp, :dim_subsp)&
      , dim_subsp, t%h_eval(:dim_subsp), t%work, t%lwork, info)
    if (info  /= 0 ) stop 'h_eval dsyev'

    ! compute eq. 15
    t%h_evec = 0.d0
    do i = 1, dim_subsp
      do k = 1, dim_subsp
        t%h_evec(:, i) = t%h_evec(:, i) + t%h_evec_subsp(k, i) * t%dr_subsp(:, k)
      end do
    end do

    ! compute eq. 20
    do j = 1, dim_subsp
      t%res_temp = - t%h_eval(j) * t%h_evec(:, j)
      do k = 1, dim_subsp
        t%res_temp = t%res_temp + t%h_evec_subsp(k, j) * t%df_subsp(:, k)
      end do
      t%res(j) = norm2(t%res_temp)
    end do

    ! modify eigenvalues (eq. 18)
    do i = 1, dim_subsp
      t%h_eval(i) = sqrt(t%h_eval(i)**2 + t%res(i)**2)
    end do

    !print*, 'meigenvaluse', t%h_eval(:dim_subsp)

    ! decompose gradient (eq. 16)
    t%dir_of_descent = df_dx
    do i = 1, dim_subsp
      t%dir_of_descent = t%dir_of_descent &
        - sum(t%h_evec(:, i)*df_dx) * t%h_evec(:, i)
    end do
    t%dir_of_descent = t%dir_of_descent * t%alpha

    ! apply preconditioning to remaining gradient (eq. 21)
    do i = 1, dim_subsp
      t%dir_of_descent = t%dir_of_descent & 
        + sum(df_dx * t%h_evec(:, i)) * t%h_evec(:, i) / t%h_eval(i)
    end do

    t%dir_of_descent = - t%dir_of_descent

  end if

  dir_of_descent = t%dir_of_descent
  t%expected_positions = x + t%dir_of_descent
  t%prev_f = f_of_x
  t%prev_df_dx = df_dx
  
end subroutine sqnm_step

function get_lower_bound(t) result(lower_bound)
  !! calculates an energy uncertainty (see eq. 20 of vc-sqnm paper)
  !! The estimate is only accurate when the optimization is converged.
  class(sqnm_optimizer) :: t
  real(c_double) :: lower_bound
  if ( t%nhist == 0 ) then
    lower_bound = 0.d0
    print*, 'no estimate of a lower bound can be given at this point.'
  else
    lower_bound = t%prev_f - .5d0 * norm2(t%prev_df_dx)**2 / t%h_eval(1)
  end if

  end function get_lower_bound

end module sqnm

module vcsqnm
  use sqnm
  use iso_c_binding
  implicit none

private :: invertalat_lattice_per_opt

  type optimizer_periodic
  !! Implementation of the vc-sqnm method. More informations about the algorithm can be found here: https://arxiv.org/abs/2206.07339
    real(c_double) :: initial_lattice(3, 3)
    real(c_double) :: initial_lattice_inv(3, 3)
    real(c_double) :: lattice_transformer(3, 3)
    real(c_double) :: lattice_transformer_inv(3, 3)
    type(sqnm_optimizer) :: sqnm_opt
    integer(c_int) :: nat
    integer(c_int) :: ndim
    real(c_double) :: initial_step_size
    real(c_double) :: w
    real(c_double) :: f_sdt_deviation
    contains
    procedure :: initialize_optimizer
    procedure :: optimizer_step
    procedure :: get_lower_energy_bound
    procedure :: close_optimizer
  end type optimizer_periodic

  type(optimizer_periodic) :: vcsqnm_opt
  !! this optimizer can be used if static access to an optimizer object is required

contains
  
  subroutine initialize_optimizer(t, nat, init_lat, initial_step_size, nhist_max &
      , lattice_weigth, alpha0, eps_subsp)
    !! This subroutine is used to set up the optimizer obect.
    class(optimizer_periodic) :: t
    integer(c_int), intent(in) :: nat
    !! Number of atoms
    real(c_double), intent(in) :: init_lat(3, 3)
    !! initial lattice vectors a = init_lat(:, 1)
    real(c_double), intent(in) :: initial_step_size
    !! initial step size. default is 1.0. For systems with hard bonds (e.g. C-C) use a value between and 1.0 and
    !! 2.5. If a system only contains weaker bonds a value up to 5.0 may speed up the convergence.
    integer(c_int), intent(in) :: nhist_max
    !! Maximal number of steps that will be stored in the history list. 
    !! Use a value between 3 and 20. Must be <= than 3*nat.
    real(c_double), intent(in) :: lattice_weigth
    !! weight or size of the supercell that is used to transform lattice derivatives. Use a value between 1 and 2. 
    !! Default is 2.
    real(c_double), intent(in) :: alpha0
    !! Lower limit on the step size. 1.e-2 is the default.
    real(c_double), intent(in) :: eps_subsp
    !! Lower limit on linear dependencies of basis vectors in history list. Default 1.e-4.
    !! Increase this parameter if energy or forces contain noise.

    integer(c_int) :: i

    t%nat = nat
    t%ndim = 3*nat + 9
    t%initial_lattice = init_lat
    t%initial_step_size = initial_step_size
    t%w = lattice_weigth
    t%f_sdt_deviation = 0.d0

    call invertalat_lattice_per_opt(t%initial_lattice, t%initial_lattice_inv)
    t%lattice_transformer = 0.d0
    do i = 1, 3
      t%lattice_transformer(i, i) = 1.d0 / norm2(t%initial_lattice(:, i))
    end do
    t%lattice_transformer = t%lattice_transformer * t%w * sqrt(dble(t%nat))
    call invertalat_lattice_per_opt(t%lattice_transformer, t%lattice_transformer_inv)

    call t%sqnm_opt%initialize_sqnm(t%ndim, nhist_max, t%initial_step_size, alpha0, eps_subsp)

  end subroutine initialize_optimizer

  subroutine close_optimizer(t)
    !! This subroutine is used to free the memory of the optimizer object.
    class(optimizer_periodic) :: t

    call t%sqnm_opt%close_sqnm()

  end subroutine close_optimizer

  subroutine optimizer_step(t, r, alat, epot, f, deralat)
    !! Calculates new atomic coordinates that are closer to the local minimum. Variable cell shape optimization.
    !! This function should be used the following way:
    !! 1. calculate energies, forces and stress tensor at positions r and lattice vectors a, b, c.
    !! 2. call the step function to update positions r and lattice vectors.
    !! 3. repeat.
    class(optimizer_periodic) :: t
    real(c_double), intent(inout) :: r(3, t%nat)
    !! Positions of the atoms
    real(c_double), intent(inout) :: alat(3, 3)
    !! lattice vectors a = alat(:, 1)
    real(c_double), intent(in) :: epot
    !! potential energy
    real(c_double), intent(in) :: f(3, t%nat)
    !! forces 
    real(c_double), intent(in) :: deralat(3, 3)
    !! negative derivative of the pot, energy w.r. to the lattice vectors

    real(c_double) :: q(3, t%nat)
    real(c_double) :: df_dq(3, t%nat)
    real(c_double) :: a_tilde(3, 3)
    real(c_double) :: df_da_tilde(3, 3)
    real(c_double) :: a_inv(3, 3)
    real(c_double) :: q_and_lat(3, t%nat + 3)
    real(c_double) :: dq_and_dlat(3, t%nat + 3)
    real(c_double) :: dd(3, t%nat + 3)
    real(c_double) :: fnoise

    fnoise = norm2(sum(f, dim=2)) / sqrt(3.d0 * t%nat)
    if ( t%f_sdt_deviation == 0.d0 ) then
      t%f_sdt_deviation = fnoise
    else
      t%f_sdt_deviation = .8d0 * t%f_sdt_deviation + .2d0 * fnoise
    end if
    if ( t%f_sdt_deviation > 0.2 * maxval(abs(f)) ) then
      print*, "Noise in force is larger than 0.2 times the larges force component. Convergence cannot be guaranteed."
    end if

    call invertalat_lattice_per_opt(alat, a_inv)

    ! transform atom coordinates and derivatives
    q = matmul(matmul(t%initial_lattice, a_inv), r)
    df_dq = - matmul(matmul(alat, t%initial_lattice_inv), f)

    ! transform lattice and derivatives
    a_tilde = matmul(alat, t%lattice_transformer)
    df_da_tilde = - matmul(deralat, t%lattice_transformer_inv)

    q_and_lat(:, :t%nat) = q
    q_and_lat(:, t%nat+1:t%nat+3) = a_tilde
    dq_and_dlat(:, :t%nat) = df_dq
    dq_and_dlat(:, t%nat+1:t%nat+3) = df_da_tilde

    call t%sqnm_opt%sqnm_step(q_and_lat, epot, dq_and_dlat, dd)
    !print*, 'dd', norm2(dd)

    q_and_lat = q_and_lat + dd

    q = q_and_lat(:, :t%nat)
    a_tilde = q_and_lat(:, t%nat+1:t%nat+3)

    alat = matmul(a_tilde, t%lattice_transformer_inv)
    r = matmul(matmul(alat, t%initial_lattice_inv), q)
  
  end subroutine optimizer_step

  subroutine invertalat_lattice_per_opt(alat, alatinv)
    !Invert alat matrix
    implicit real*8(a - h, o - z)
    dimension alat(3, 3), alatinv(3, 3)

    div = (alat(1, 1)*alat(2, 2)*alat(3, 3) - alat(1, 1)*alat(2, 3)*alat(3, 2) - &
           alat(1, 2)*alat(2, 1)*alat(3, 3) + alat(1, 2)*alat(2, 3)*alat(3, 1) + &
           alat(1, 3)*alat(2, 1)*alat(3, 2) - alat(1, 3)*alat(2, 2)*alat(3, 1))
    div = 1.d0/div
    alatinv(1, 1) = (alat(2, 2)*alat(3, 3) - alat(2, 3)*alat(3, 2))*div
    alatinv(1, 2) = -(alat(1, 2)*alat(3, 3) - alat(1, 3)*alat(3, 2))*div
    alatinv(1, 3) = (alat(1, 2)*alat(2, 3) - alat(1, 3)*alat(2, 2))*div
    alatinv(2, 1) = -(alat(2, 1)*alat(3, 3) - alat(2, 3)*alat(3, 1))*div
    alatinv(2, 2) = (alat(1, 1)*alat(3, 3) - alat(1, 3)*alat(3, 1))*div
    alatinv(2, 3) = -(alat(1, 1)*alat(2, 3) - alat(1, 3)*alat(2, 1))*div
    alatinv(3, 1) = (alat(2, 1)*alat(3, 2) - alat(2, 2)*alat(3, 1))*div
    alatinv(3, 2) = -(alat(1, 1)*alat(3, 2) - alat(1, 2)*alat(3, 1))*div
    alatinv(3, 3) = (alat(1, 1)*alat(2, 2) - alat(1, 2)*alat(2, 1))*div
  end subroutine invertalat_lattice_per_opt
  
  function get_lower_energy_bound(t) result(lower_bound)
    !! calculates an energy uncertainty (see eq. 20 of vc-sqnm paper)
    !! The estimate is only accurate when the optimization is converged.
    class(optimizer_periodic) :: t
    real(c_double) :: lower_bound
  
    lower_bound = t%sqnm_opt%get_lower_bound()
  
  end function get_lower_energy_bound

end module vcsqnm


