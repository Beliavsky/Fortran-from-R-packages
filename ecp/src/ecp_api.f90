module ecp_api
  use ecp_kinds, only: dp
  use ecp_types, only: cp_result, agglo_result, divisive_result
  use ecp_energy, only: get_within, get_between
  use ecp_cp3o, only: e_cp3o, e_cp3o_delta, ks_cp3o, ks_cp3o_delta
  use ecp_agglo, only: e_agglo
  use ecp_divisive, only: e_divisive
  use ecp_kernel, only: kcpa
  implicit none
  private

  public :: dp, cp_result, agglo_result, divisive_result
  public :: get_within, get_between, e_cp3o, e_cp3o_delta
  public :: ks_cp3o, ks_cp3o_delta, e_agglo, e_divisive, kcpa
end module ecp_api
