program test_orthogonal
  use polynom
  implicit none
  real(dp) :: x(4), v(4,4), gram(4,4)
  type(polylist_t) :: basis
  integer :: j

  x = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
  basis = orthogonal_polynomials(x, degree=3, normalized=.true.)
  if (basis%size() /= 4) error stop 'wrong basis size'
  do j = 1, 4
    v(:,j) = basis%item(j)%evaluate(x)
  end do
  gram = matmul(transpose(v), v)
  call assert_identity(gram, 1.0e-10_dp)
  call assert_close(basis%item(1)%coef, [0.5_dp], 1.0e-12_dp)
  call assert_close(basis%item(2)%coef, [-0.5916079783099616_dp, 0.3380617018914066_dp], 1.0e-10_dp)

  basis = orthogonal_polynomials(x, degree=3, normalized=.false.)
  do j = 1, 4
    if (abs(basis%item(j)%coef(size(basis%item(j)%coef)) - 1.0_dp) > 1.0e-10_dp) then
      error stop 'unnormalized basis is not monic'
    end if
  end do

  print *, 'test_orthogonal: PASS'
contains
  subroutine assert_identity(a, tol)
    real(dp), intent(in) :: a(:,:), tol
    integer :: i, j
    do j = 1, size(a,2)
      do i = 1, size(a,1)
        if (i == j) then
          if (abs(a(i,j) - 1.0_dp) > tol) error stop 'not identity'
        else
          if (abs(a(i,j)) > tol) error stop 'not orthogonal'
        end if
      end do
    end do
  end subroutine assert_identity
  subroutine assert_close(a, b, tol)
    real(dp), intent(in) :: a(:), b(:), tol
    if (size(a) /= size(b)) error stop 'shape mismatch'
    if (maxval(abs(a-b)) > tol) error stop 'values differ'
  end subroutine assert_close
end program test_orthogonal
