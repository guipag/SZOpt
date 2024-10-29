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

  subroutine compute_Ht(H_t, HB, HD, beta, lambda, size1, size2)
    implicit none
    integer, intent(in) :: size1, size2
    real(8), intent(in) :: beta, lambda
    real(8), intent(in) :: HB(size1, size2), HD(size1, size2)
    real(8), intent(out) :: H_t(size2, size2)
    real(8) :: temp_HB(size2, size2), temp_HD(size2, size2)
    integer :: i

    ! Compute temp_HB = HB^T * HB
    call dgemm('T', 'N', size2, size2, size1, 1.0d0, HB, size1, HB, size1, 0.0d0, temp_HB, size2)

    ! Compute temp_HD = HD^T * HD
    call dgemm('T', 'N', size2, size2, size1, beta, HD, size1, HD, size1, 0.0d0, temp_HD, size2)

    ! Combine results: H_t = (1 - beta) * temp_HB + temp_HD
    H_t = (1.0d0 - beta) * temp_HB + temp_HD

    ! Add lambda * I to H_t
    do i = 1, size2
        H_t(i, i) = H_t(i, i) + lambda
    end do
  end subroutine compute_Ht

  subroutine solveGd(HB, HD, d, beta, lambda, nbIter)
    implicit none
    real(8), intent(in) :: HB(:, :), HD(:, :), d(:)
    real(8), intent(in) :: beta, lambda
    integer, intent(in), optional :: nbIter

    integer :: size1, size2, i, max_iter
    real(8) :: mu
    real(8), allocatable :: H_t(:, :), r(:), h(:), temp_r(:), temp_h(:)

    size1 = size(HB, dim=1)
    size2 = size(HB, dim=2)

    allocate(H_t(size2, size2), r(size2), h(size2), temp_r(size2), temp_h(size2))

    ! Set the maximum number of iterations
    max_iter = merge(size2, nbIter, present(nbIter))

    ! Compute H_t
    call compute_Ht(H_t, HB, HD, beta, lambda, size1, size2)

    ! Initialize vectors
    r = 0.0d0
    h = 0.0d0

    ! Compute initial residual r = HB * d
    call dgemv('N', size1, size2, 1.0d0, HB, size1, d, 1, 0.0d0, temp_r, 1)

    do i = 1, max_iter
        ! Matrix-vector product temp_h = H_t * h
        call dsymv('U', size2, 1.0d0, H_t, size2, h, 1, 0.0d0, temp_h, 1)

        ! Update residual r
        r = -(1.0d0 - beta) * temp_r - temp_h

        ! Matrix-vector product temp_h = H_t * r
        call dsymv('U', size2, 1.0d0, H_t, size2, r, 1, 0.0d0, temp_h, 1)

        ! Compute step size mu
        mu = ddot(size2, r, 1, r, 1) / ddot(size2, r, 1, temp_h, 1)

        ! Update h
        h = h - mu * r
    end do

    ! Clean up
    deallocate(H_t, r, h, temp_r, temp_h)
  end subroutine solveGd

  subroutine solvegdc(HB, HD, d, beta, lambda, nbIter)
    implicit none
    real(8), intent(in) :: HB(:, :), HD(:, :), d(:)
    real(8), intent(in) :: beta, lambda
    integer, intent(in), optional :: nbIter

    real(8), allocatable :: H_t(:, :), h(:), r(:), p(:), r_old(:), H_tp(:)
    real(8) :: mu, bet

    integer :: size1, size2, k, max_iter

    size1 = size(HB, dim=1)
    size2 = size(HB, dim=2)

    allocate(H_t(size2, size2), h(size2), r(size2), p(size2), r_old(size2), H_tp(size2))

    ! Set the maximum number of iterations
    max_iter = merge(size2, nbIter, present(nbIter))

    ! Compute H_t
    call compute_Ht(H_t, HB, HD, beta, lambda, size1, size2)

    ! Initialize vectors
    h = 0.0d0
    r = 0.0d0
    p = 0.0d0

    ! Compute initial residual temp_r = HB * d
    call dgemv('N', size1, size2, 1.0d0, HB, size1, d, 1, 0.0d0, r, 1)

    ! Compute initial r = (1 - beta) * temp_r - H_t * h
    call dsymv('U', size2, 1.0d0, H_t, size2, h, 1, 0.0d0, r_old, 1)
    r = (1.0d0 - beta) * r - r_old
    p = r

    do k = 1, max_iter
        ! Matrix-vector product H_tp = H_t * p
        call dsymv('U', size2, 1.0d0, H_t, size2, p, 1, 0.0d0, H_tp, 1)

        ! Compute step size mu
        mu = ddot(size2, r, 1, p, 1) / ddot(size2, p, 1, H_tp, 1)

        ! Update h
        h = h + mu * p

        ! Save current r for beta calculation
        r_old = r

        ! Update residual r
        r = r - mu * H_tp

        ! Compute beta
        bet = ddot(size2, r, 1, r, 1) / ddot(size2, r_old, 1, r_old, 1)

        ! Update p
        p = r + bet * p
    end do

    ! Clean up
    deallocate(H_t, h, r, p, r_old, H_tp)
  end subroutine solvegdc

end module solve_f
