program demo_orthopolynom
  use orthopolynom
  implicit none
  type(polylist_t) :: p
  type(monic_recurrence_t) :: mr
  type(real_vector_list_t) :: roots
  integer :: k

  p = legendre_polynomials(5)
  print '(a)', 'Legendre polynomials P_0 through P_5:'
  do k = 1, p%size()
    print '(a,i0,a,a)', 'P_', k - 1, '(x) = ', p%item(k)%to_string()
  end do

  mr = monic_polynomial_recurrences(chebyshev_t_recurrences(5))
  roots = polynomial_roots(mr)
  print '(a)', ''
  print '(a)', 'Roots of T_5:'
  print '(*(f12.8,1x))', roots%item(6)%value
end program demo_orthopolynom
