module Solver_mod
  use :: types_mod, only: dp
  use RHS_mod
  implicit none

  public

contains

  subroutine fd1d_heat_explicit(x, t, dt, cfl, h, h_new)
    implicit none

    real (kind=dp), intent(in) :: cfl
    real (kind=dp), intent(in) :: dt
    real (kind=dp), intent(in) :: h(:)
    real (kind=dp), intent(out) :: h_new(:)
    integer :: j
    real (kind=dp), intent(in) :: t
    real (kind=dp), intent(in) :: x(:)
    real (kind=dp) :: f(size(x))

    do j = 1, size(x)
      f(j) = func(j, x)
    end do

    h_new(1) = 0.0e+00_dp

    do j = 2, size(x) - 1
      h_new(j) = h(j) + dt*f(j) + cfl*(h(j-1)-2.0e+00_dp*h(j)+h(j+1))
    end do

! set the boundary conditions again
    h_new(1) = 90.0e+00_dp
    h_new(size(x)) = 70.0e+00_dp
  end subroutine

end module Solver_mod
