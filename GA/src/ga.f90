! Public convenience module for GA-fortran.
module ga
  use ga_kinds, only : dp
  use ga_random, only : ga_seed
  use ga_utils, only : decimal2binary, binary2decimal, binary2gray, gray2binary
  use ga_utils, only : ga_pmutation, garun, fitness_summary, optim_probsel
  use ga_utils, only : repair_solution, reflect_solution
  use ga_operators
  use ga_core
  use ga_islands
  implicit none
  public
  interface ga_optimize
    module procedure ga_real
    module procedure ga_binary
    module procedure ga_permutation
  end interface ga_optimize
end module ga
