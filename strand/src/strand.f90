! SPDX-License-Identifier: GPL-3.0-only
! Upstream authors: Jeff Enos, David Kane, and strand contributors.
! Numerical translation of strand 0.2.3.
module strand
  use strand_kinds, only : dp
  use strand_types
  use strand_stats
  use strand_data
  use strand_optimizer
  use strand_simulation
  implicit none
  public
end module strand
