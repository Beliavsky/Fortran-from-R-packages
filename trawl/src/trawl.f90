module trawl
  use trawl_kinds, only : dp
  use trawl_types
  use trawl_rng, only : set_trawl_seed
  use trawl_functions
  use trawl_distributions
  use trawl_statistics, only : empirical_acf,sample_mean,sample_variance
  use trawl_fit
  use trawl_intersection
  use trawl_simulation
  implicit none
  public
end module trawl
