program test_tuning_features
  use mlr_kinds, only : dp, i8
  use mlr_rng, only : rng_state, rng_seed
  use mlr_types, only : tune_result, feature_select_result
  use mlr_tuning, only : grid_search, random_search
  use mlr_feature_selection, only : feature_select_exhaustive, feature_select_forward
  implicit none
  real(dp)::grid(5,1),lower(1),upper(1)
  type(tune_result)::tr
  type(feature_select_result)::fr
  type(rng_state)::rng
  integer::i
  do i=1,5;grid(i,1)=real(i-3,dp);end do
  call grid_search(grid,obj,.true.,tr)
  if(abs(tr%par(1)-1.0_dp)>1.0e-12_dp)error stop 'grid'
  lower=-5.0_dp;upper=5.0_dp;call rng_seed(rng,1_i8);call random_search(lower,upper,1000,obj,.true.,rng,tr)
  if(abs(tr%par(1)-1.0_dp)>0.05_dp)error stop 'random'
  call feature_select_exhaustive(4,fsobj,.true.,fr)
  if(.not.fr%selected(2).or.count(fr%selected)/=1)error stop 'exhaustive fs'
  call feature_select_forward(4,fsobj,.true.,fr)
  if(.not.fr%selected(2).or.count(fr%selected)/=1)error stop 'forward fs'
  print *, 'test_tuning_features: PASS'
contains
  real(dp) function obj(par)
    real(dp),intent(in)::par(:);obj=(par(1)-1.0_dp)**2
  end function
  real(dp) function fsobj(sel)
    logical,intent(in)::sel(:)
    fsobj=real(count(sel),dp)
    if(sel(2))fsobj=fsobj-3.0_dp
  end function
end program test_tuning_features
