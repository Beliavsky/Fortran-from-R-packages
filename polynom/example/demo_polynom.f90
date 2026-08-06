program demo_polynom
  use polynom
  implicit none
  type(polynomial_t) :: p, shifted
  type(polylist_t) :: basis
  complex(dp), allocatable :: roots(:)
  real(dp) :: x(4), v(4,4)
  integer :: i

  p = poly_from_roots([1.0_dp, 2.0_dp, 3.0_dp, 4.0_dp, 5.0_dp])
  shifted = change_origin(p, 3.0_dp)
  roots = polynomial_roots(p)

  print '(a)', 'p(x) = ' // p%to_string()
  print '(a)', 'p(x + 3) = ' // shifted%to_string()
  print '(a)', 'roots:'
  do i = 1, size(roots)
    print '(2f12.6)', real(roots(i), dp), aimag(roots(i))
  end do

  x = [0.0_dp, 1.0_dp, 2.0_dp, 4.0_dp]
  basis = orthogonal_polynomials(x, degree=3)
  do i = 1, 4
    v(:,i) = basis%item(i)%evaluate(x)
  end do
  print '(a,es12.4)', 'maximum orthonormality error: ', &
    maxval(abs(matmul(transpose(v), v) - identity4()))
contains
  pure function identity4() result(a)
    real(dp) :: a(4,4)
    integer :: j
    a = 0.0_dp
    do j = 1, 4
      a(j,j) = 1.0_dp
    end do
  end function identity4
end program demo_polynom
