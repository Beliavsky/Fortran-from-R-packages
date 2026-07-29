! SPDX-License-Identifier: GPL-3.0-or-later
! Copyright (C) 2018 Marius Hofert, Erik Hintz and Christiane Lemieux
module nvmix_mixing
  use nvmix_kinds, only : dp
  use nvmix_types, only : mix_constant,mix_inverse_gamma,mix_pareto,mix_gamma
  use nvmix_special, only : inverse_gamma_quantile,gamma_quantile
  use nvmix_random, only : random_uniform,random_gamma
  implicit none
  private
  public :: mixing_quantile,mixing_random,mixing_mean,mixing_mean_sqrt
contains
  real(dp) function mixing_quantile(u,family,parameter) result(w)
    real(dp), intent(in) :: u,parameter
    integer, intent(in) :: family
    select case(family)
    case(mix_constant); w=1.0_dp
    case(mix_inverse_gamma); w=inverse_gamma_quantile(u,parameter)
    case(mix_pareto); w=(1.0_dp-u)**(-1.0_dp/parameter)
    case(mix_gamma); w=gamma_quantile(u,parameter,1.0_dp/parameter)
    case default; w=1.0_dp
    end select
    w=max(w,tiny(1.0_dp))
  end function
  real(dp) function mixing_random(family,parameter) result(w)
    real(dp), intent(in) :: parameter
    integer, intent(in) :: family
    select case(family)
    case(mix_constant); w=1.0_dp
    case(mix_inverse_gamma); w=1.0_dp/random_gamma(0.5_dp*parameter,2.0_dp/parameter)
    case(mix_pareto); w=(1.0_dp-random_uniform())**(-1.0_dp/parameter)
    case(mix_gamma); w=random_gamma(parameter,1.0_dp/parameter)
    case default; w=1.0_dp
    end select
    w=max(w,tiny(1.0_dp))
  end function
  pure real(dp) function mixing_mean(family,parameter) result(v)
    real(dp), intent(in) :: parameter
    integer, intent(in) :: family
    select case(family)
    case(mix_constant,mix_gamma); v=1.0_dp
    case(mix_inverse_gamma)
      if(parameter>2.0_dp)then; v=parameter/(parameter-2.0_dp); else; v=huge(1.0_dp); end if
    case(mix_pareto)
      if(parameter>1.0_dp)then; v=parameter/(parameter-1.0_dp); else; v=huge(1.0_dp); end if
    case default; v=1.0_dp
    end select
  end function
  pure real(dp) function mixing_mean_sqrt(family,parameter) result(v)
    real(dp), intent(in) :: parameter
    integer, intent(in) :: family
    select case(family)
    case(mix_constant); v=1.0_dp
    case(mix_inverse_gamma)
      if(parameter>1.0_dp)then
        v=sqrt(0.5_dp*parameter)*exp(log_gamma(0.5_dp*(parameter-1.0_dp))-log_gamma(0.5_dp*parameter))
      else; v=huge(1.0_dp); end if
    case(mix_pareto)
      if(parameter>0.5_dp)then; v=parameter/(parameter-0.5_dp); else; v=huge(1.0_dp); end if
    case(mix_gamma)
      v=exp(log_gamma(parameter+0.5_dp)-log_gamma(parameter))/sqrt(parameter)
    case default; v=1.0_dp
    end select
  end function
end module nvmix_mixing
