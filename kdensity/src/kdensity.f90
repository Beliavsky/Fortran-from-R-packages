module kdensity
  use kdensity_kinds
  use kdensity_types
  use kdensity_math, only : normal_pdf, normal_cdf, normal_quantile, adaptive_integral
  use kdensity_starts, only : get_start, supported_starts
  use kdensity_kernels, only : get_kernel, supported_kernels
  use kdensity_bandwidths
  use kdensity_core
  implicit none
  public
end module kdensity
