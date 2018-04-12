module RHS_mod
  use :: types_mod, only: dp
  implicit none

  public

contains
    function func(j, x) result (d)
    implicit none

    integer, intent(in) :: j
    real (kind=dp) :: d
    real (kind=dp), dimension(:), intent(in) :: x

    d = 0.0e+00_dp
  end function func

end module RHS_mod
