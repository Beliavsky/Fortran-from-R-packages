! Public convenience module for the rmoo Fortran translation.
module rmoo
  use ga_kinds, only : dp
  use ga_random, only : ga_seed
  use rmoo_pareto
  use rmoo_reference
  use rmoo_operators
  use rmoo_survival
  use rmoo_core
  implicit none
  public
end module rmoo
