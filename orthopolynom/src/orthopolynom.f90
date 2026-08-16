module orthopolynom
  use polynom, only : dp, polynomial_t, polylist_t, polynomial, derivative, integral_polynomial, &
    definite_integral, operator(+), operator(-), operator(*), operator(/), operator(**)
  use orthopolynom_types
  use orthopolynom_core
  use orthopolynom_families
  implicit none
  public
end module orthopolynom
