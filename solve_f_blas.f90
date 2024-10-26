module solve_f
  implicit none

   interface
    function ddot(n, dx, incx, dy, incy) result(dotprod)
      integer, intent(in) :: n, incx, incy
      real(8), intent(in) :: dx(*), dy(*)
      real(8) :: dotprod
    end function ddot
      end interface

contains

  subroutine solveGd(nbIter, beta, HB, d, H_t, size1, size2)
    implicit none
    integer, intent(in) :: nbIter, size1, size2
    real(8), intent(in) :: beta
    real(8), intent(in) :: HB(size1, size2), d(size1), H_t(size2, size2)
    real(8) :: r(size2), h(size2), temp_r(size2), temp_h(size2), mu
    integer :: i

    r = 0.0d0
    h = 0.0d0

    call dgemv('N', size1, size2, 1.0d0, HB, size1, d, 1, 0.0d0, temp_r, 1)

    do i = 1, nbIter
        call dsymv('U', size2, 1.0d0, H_t, size2, h, 1, 0.0d0, temp_h, 1)

        r = -(1.0d0 - beta) * temp_r - temp_h

        call dsymv('U', size2, 1.0d0, H_t, size2, r, 1, 0.0d0, temp_h, 1)
        
        mu = ddot(size2, r, 1, r, 1) / ddot(size2, r, 1, temp_h, 1)

        h = h - mu * r
    end do
  end subroutine solveGd


 subroutine solvegdc(nbIter, beta, HB, d, H_t, size1, size2)
    implicit none
    integer, intent(in) :: nbIter
    real(8), intent(in) :: beta
    real(8), intent(in) :: HB(size1, size2), d(size1), H_t(size2, size2)
    integer, intent(in) :: size1, size2

    real(8) :: h(size2), r(size2), p(size2), r_old(size2), H_tp(size2), temp_r(size1)
    real(8) :: mu, bet
    integer :: iter

    h = 0.0d0

    call dgemv('N', size1, size2, 1.0d0, HB, size1, d, 1, 0.0d0, temp_r, 1)
    call dsymv('U', size2, 1.0d0, H_t, size2, h, 1, 0.0d0, r, 1)

    r = (1.0d0 - beta) * temp_r - r
    p = r

    do iter = 1, nbIter
        call dsymv('U', size2, 1.0d0, H_t, size2, p, 1, 0.0d0, H_tp, 1)

        mu = ddot(size2, r, 1, p, 1) / ddot(size2, p, 1, H_tp, 1)

        h = h + mu * p

        r_old = r

        r = r - mu * H_tp

        bet = ddot(size2, r, 1, r, 1) / ddot(size2, r_old, 1, r_old, 1)

        p = r + bet * p
    end do
  end subroutine solvegdc
end module solve_f
