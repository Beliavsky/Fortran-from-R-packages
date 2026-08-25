module mnormt
  use mnormt_special, only: dp, normal_pdf, normal_cdf, normal_quantile, student_t_pdf, student_t_cdf
  use mnormt_linalg, only: pd_solve
  use mnormt_core
  use mnormt_truncated
  use mnormt_moments
  implicit none
  public
end module mnormt
