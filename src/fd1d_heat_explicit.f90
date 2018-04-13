!> Andrew Williams, University of Bristol
!> Solves the one dimensional heat diffusion equation
!> \( \(frac{\partial H}{\partial t}
!>        \kappa\frac{\partial^{2} H}{\partial z^{2}} = f(x) \)
program fd1d_heat_explicit_prb
  use :: types_mod, only: dp
  use RHS_mod
  use CFL_mod
  use IO_mod
  use Solver_mod

  implicit none

  integer, parameter :: T_NUM = 201
  integer, parameter :: X_NUM = 21

  real (kind=dp) :: cfl
  real (kind=dp) :: dt

! error code for failed allocations
  integer :: alloc_err
  integer :: dealloc_err

  real (kind=dp), allocatable :: h(:)
  real (kind=dp), allocatable :: h_new(:)
! the "matrix" stores all x-values for all t-values
! remember Fortran is column major, meaning that rows are contiguous
!  real (kind=dp) :: hmat(x_num, t_num)
  real (kind=dp), allocatable :: hmat(:, :)
  integer :: i
  integer :: j
  real (kind=dp) :: k

!  real (kind=dp) :: t(t_num)
  real (kind=dp), allocatable :: t(:)
  real (kind=dp) :: t_max
  real (kind=dp) :: t_min
!  real (kind=dp) :: x(x_num)
  real (kind=dp), allocatable :: x(:)
  real (kind=dp) :: x_max
  real (kind=dp) :: x_min

  allocate(h(x_num), stat=alloc_err)

  if (alloc_err /= 0) then
    print *, "Failed to allocate h"
    stop
  end if

  allocate(h_new(x_num), stat=alloc_err)

  if (alloc_err /= 0) then
    print *, "Failed to allocate h_new"
    stop
  end if

  allocate(hmat(x_num, t_num), stat=alloc_err)

  if (alloc_err /= 0) then
    print *, "Failed to allocate hmat"
    stop
  end if

  allocate(t(t_num), stat=alloc_err)

  if (alloc_err /= 0) then
    print *, "Failed to allocate t"
    stop
  end if

  allocate(x(x_num), stat=alloc_err)

  if (alloc_err /= 0) then
    print *, "Failed to allocate x"
    stop
  end if

  print *, "Hey, I allocated all the dynamic arrays!"

  write (*, '(a)') ' '
  write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
  write (*, '(a)') '  FORTRAN77 version.'
  write (*, '(a)') '  Test the FD1D_HEAT_EXPLICIT library.'

  write (*, '(a)') ' '
  write (*, '(a)') 'FD1D_HEAT_EXPLICIT_PRB:'
  write (*, '(a)') '  Normal end of execution.'
  write (*, '(a)') ' '

  write (*, '(a)') ' '
  write (*, '(a)') 'FD1D_HEAT_EXPLICIT_TEST01:'
  write (*, '(a)') '  Compute an approximate solution to the time-dependent'
  write (*, '(a)') '  one dimensional heat equation:'
  write (*, '(a)') ' '
  write (*, '(a)') '    dH/dt - K * d2H/dx2 = f(x,t)'
  write (*, '(a)') ' '
  write (*, '(a)') '  Run a simple test case.'

! heat coefficient
  k = 0.002e+00_dp

! the x-range values
  x_min = 0.0e+00_dp
  x_max = 1.0e+00_dp
! x_num is the number of intervals in the x-direction
  call r8vec_linspace(x_min, x_max, x)

! the t-range values. integrate from t_min to t_max
  t_min = 0.0e+00_dp
  t_max = 80.0e+00_dp

! t_num is the number of intervals in the t-direction
  dt = (t_max-t_min)/real(t_num-1, kind=dp)
  call r8vec_linspace(t_min, t_max, t)

! get the CFL coefficient
  call fd1d_heat_explicit_cfl(k, t_num, t_min, t_max, x_num, x_min, x_max, &
    cfl)

  if (0.5e+00_dp<=cfl) then
    write (*, '(a)') ' '
    write (*, '(a)') 'FD1D_HEAT_EXPLICIT_CFL - Fatal error!'
    write (*, '(a)') '  CFL condition failed.'
    write (*, '(a)') '  0.5 <= K * dT / dX / dX = CFL.'
    stop
  end if

! set the initial condition
  do j = 1, x_num
    h(j) = 50.0e+00_dp
  end do

! set the bounday condition
  h(1) = 90.0e+00_dp
  h(x_num) = 70.0e+00_dp

! initialise the matrix to the initial condition
  do i = 1, x_num
    hmat(i, 1) = h(i)
  end do

! the main time integration loop 
  do j = 2, t_num
    call fd1d_heat_explicit(x, t(j-1), dt, cfl, h, h_new)

    do i = 1, x_num
      hmat(i, j) = h_new(i)
      h(i) = h_new(i)
    end do
  end do

! write data to files
  call r8mat_write('h_test01.txt', hmat)
  call r8vec_write('t_test01.txt', t)
  call r8vec_write('x_test01.txt', x)

! Deallocate all my dynamic memory
  deallocate(h, stat=dealloc_err)

  if (dealloc_err /= 0) then
    print *, "Failed to deallocate h"
    stop
  end if

  deallocate(h_new, stat=dealloc_err)

  if (dealloc_err /= 0) then
    print *, "Failed to deallocate h_new"
    stop
  end if

  deallocate(hmat, stat=dealloc_err)

  if (dealloc_err /= 0) then
    print *, "Failed to deallocate hmat"
    stop
  end if

  deallocate(t, stat=dealloc_err)

  if (dealloc_err /= 0) then
    print *, "Failed to deallocate t"
    stop
  end if

  deallocate(x, stat=dealloc_err)

  if (dealloc_err /= 0) then
    print *, "Failed to deallocate x"
    stop
  end if

  print *, "Successfully deallocated dynamic memory"

end program
