module nfcp
  use nfcp_types
  use nfcp_math, only : nfcp_a_t, nfcp_covariance, nfcp_seasonality, normal_cdf, normal_quantile
  use nfcp_parameters, only : pack_parameters, unpack_parameters, default_parameter_bounds
  use nfcp_kalman, only : nfcp_kalman_filter
  use nfcp_forecast, only : spot_price_forecast, futures_price_forecast
  use nfcp_simulation, only : spot_price_simulate, futures_price_simulate
  use nfcp_options, only : european_option_value, american_option_value
  use nfcp_analysis, only : tsfit_volatility
  use nfcp_stitch, only : stitch_contract_numbers, stitch_by_maturity
  use nfcp_mle, only : nfcp_fit_mle
  implicit none
  public
end module nfcp
