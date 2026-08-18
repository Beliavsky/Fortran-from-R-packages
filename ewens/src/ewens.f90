! SPDX-License-Identifier: MIT
module ewens
  use ewens_kinds, only : dp, i8
  use ewens_distribution, only : dewens, dewens_log, dewens_counts, dewens_counts_log, &
    dewens_k, dewens_k_log, ewens_k_exact
  use ewens_sampling, only : ewens_seed, rewens, gcrp, rgem
  use ewens_estimation, only : ewens_mle, ewens_mle_nk, ewens_score
  use ewens_math, only : number_of_classes, class_size_spectrum
  implicit none
  private
  public :: dp, i8
  public :: dewens, dewens_log, dewens_counts, dewens_counts_log
  public :: dewens_k, dewens_k_log, ewens_k_exact
  public :: ewens_seed, rewens, gcrp, rgem
  public :: ewens_mle, ewens_mle_nk, ewens_score
  public :: number_of_classes, class_size_spectrum
end module ewens
