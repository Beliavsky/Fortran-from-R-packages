! SPDX-License-Identifier: BSD-3-Clause
! Modern Fortran computational translation of waveslim.
module waveslim
  use waveslim_kinds
  use waveslim_status
  use waveslim_types
  use waveslim_filters
  use waveslim_math
  use waveslim_transform_1d
  use waveslim_packet
  use waveslim_transform_nd
  use waveslim_dualtree
  use waveslim_hilbert_stats
  use waveslim_extended
  use waveslim_statistics
  use waveslim_denoise
  use waveslim_long_memory
  implicit none
  public
end module waveslim
